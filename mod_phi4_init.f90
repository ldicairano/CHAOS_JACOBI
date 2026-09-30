module mod_phi4_init
  use mod_kinds,           only: dp
  use mod_types_phi4,      only: Phi4Params, Phi4State, RNGState
  use mod_random,          only: rng_uniform, rng_gauss
  use mod_phi4_model,      only: compute_potential, compute_kinetic
  use mod_bins_phi4,       only: init_bins, reset_bins
  use mod_metropolis_phi4, only: vv_step_md_tangent_jlc
  implicit none
  private

  public :: allocate_state
  public :: init_configuration
  public :: init_bins_autorange
contains

  subroutine allocate_state(par, st)
    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(inout) :: st

    allocate(st%phi(par%N))
    allocate(st%pi(par%N))
    allocate(st%force(par%N))

    allocate(st%dphi(par%N))
    allocate(st%dpi_t(par%N))

    allocate(st%jphi(par%N))
    allocate(st%jpi(par%N))

    allocate(st%gradV(par%N))
    allocate(st%hdiag(par%N))

    st%phi   = 0.0_dp
    st%pi    = 0.0_dp

    st%force = 0.0_dp
    st%gradV = 0.0_dp
    st%hdiag = 0.0_dp

    st%dphi  = 0.0_dp
    st%dpi_t = 0.0_dp

    st%jphi  = 0.0_dp
    st%jpi   = 0.0_dp
  end subroutine allocate_state



  subroutine init_configuration(par, st, rng)
    type(Phi4Params), intent(inout) :: par
    type(Phi4State),  intent(inout) :: st
    type(RNGState),   intent(inout) :: rng

    integer :: i, itry
    real(dp) :: amp, Ktarget

    ! We need V(phi) < Etot, otherwise no real momenta exist.
    amp = 1.0_dp

    do itry = 1, 20
      do i = 1, par%N
        st%phi(i) = amp * (2.0_dp*rng_uniform(rng) - 1.0_dp)
      end do

      call compute_potential(par, st)

      if (st%V < par%Etot) exit
      amp = 0.5_dp * amp
    end do

    if (st%V >= par%Etot) then
      ! Safe fallback: phi = 0
      st%phi = 0.0_dp
      call compute_potential(par, st)
    end if

    if (st%V >= par%Etot) then
      write(*,*) "ERROR: could not build a microcanonical initial condition with V < Etot."
      write(*,*) "       V =", st%V, " Etot =", par%Etot
      stop
    end if

    ! Random Gaussian momenta, then rescale to K = Etot - V
    do i = 1, par%N
      st%pi(i) = rng_gauss(rng)
    end do

    call compute_kinetic(par, st)

    Ktarget = par%Etot - st%V
    call rescale_momenta_to_target_kinetic(st, Ktarget)

    call compute_kinetic(par, st)
    call compute_potential(par, st)

    ! Consistency check
    if (abs((st%K + st%V) - par%Etot) > 1.0e-10_dp * max(1.0_dp, abs(par%Etot))) then
      write(*,*) "WARNING: micro init energy mismatch:"
      write(*,*) "  K+V =", st%K + st%V, " Etot =", par%Etot
    end if

    call init_variation_vectors(par, st, rng)
    call reset_geometric_diagnostics(st)
    call initialize_md_geometry(par, st)
  end subroutine init_configuration


  subroutine rescale_momenta_to_target_kinetic(st, Ktarget)
    type(Phi4State), intent(inout) :: st
    real(dp),        intent(in)    :: Ktarget

    real(dp) :: K0, fac

    if (Ktarget <= 0.0_dp) then
      write(*,*) "ERROR: target kinetic energy <= 0 in micro init."
      stop
    end if

    K0 = 0.5_dp * sum(st%pi * st%pi)
    if (K0 <= 0.0_dp) then
      write(*,*) "ERROR: zero initial kinetic energy before rescaling."
      stop
    end if

    fac   = sqrt(Ktarget / K0)
    st%pi = fac * st%pi
  end subroutine rescale_momenta_to_target_kinetic


  subroutine init_variation_vectors(par, st, rng)
    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(inout) :: st
    type(RNGState),   intent(inout) :: rng

    integer :: i
    real(dp) :: nd, nj

    do i = 1, par%N
      st%dphi(i)  = rng_gauss(rng)
      st%dpi_t(i) = rng_gauss(rng)

      st%jphi(i)  = rng_gauss(rng)
      st%jpi(i)   = rng_gauss(rng)
    end do

    nd = sqrt(sum(st%dphi*st%dphi) + sum(st%dpi_t*st%dpi_t))
    nj = sqrt(sum(st%jphi*st%jphi) + sum(st%jpi*st%jpi))

    st%dphi  = st%dphi  / nd
    st%dpi_t = st%dpi_t / nd

    st%jphi = st%jphi / nj
    st%jpi  = st%jpi  / nj
  end subroutine init_variation_vectors


  subroutine reset_geometric_diagnostics(st)
    type(Phi4State), intent(inout) :: st

    st%jacW    = 0.0_dp
    st%gradV2  = 0.0_dp
    st%lapV    = 0.0_dp
    st%jacR    = 0.0_dp
    st%jacRreg = 0.0_dp

    st%lya_logsum = 0.0_dp
    st%jlc_logsum = 0.0_dp
    st%lya_time   = 0.0_dp
    st%jlc_time   = 0.0_dp
  end subroutine reset_geometric_diagnostics



  subroutine initialize_md_geometry(par, st)
    use mod_kinds, only: dp
    use mod_types_phi4, only: Phi4Params, Phi4State
    implicit none

    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(inout) :: st

    integer  :: i, N
    real(dp) :: Ni, pref
    real(dp) :: qi, qi2, qi4
    real(dp) :: s1, s2, s4
    real(dp) :: dVi, hdiag_i, Vloc
    real(dp) :: grad2, lap, W, Rreg

    N    = par%N
    Ni   = real(N, dp)
    pref = 4.0_dp * par%coup / real(N - 1, dp)

    s1 = 0.0_dp
    s2 = 0.0_dp
    s4 = 0.0_dp

    do i = 1, N
      qi  = st%phi(i)
      qi2 = qi*qi
      qi4 = qi2*qi2
      s1  = s1 + qi
      s2  = s2 + qi2
      s4  = s4 + qi4
    end do

    st%S1   = s1
    st%S2   = s2
    st%S4   = s4
    st%M    = s1 / Ni
    st%phi2 = s2 / Ni
    st%phi4 = s4 / Ni

    Vloc = (par%lambda/24.0_dp) * s4 - 0.5_dp * par%mu * s2 &
         + 2.0_dp * par%coup / real(N - 1, dp) * (Ni*s2 - s1*s1)

    st%V = Vloc

    grad2 = 0.0_dp
    lap   = 0.0_dp

    do i = 1, N
      qi  = st%phi(i)
      qi2 = qi*qi

      dVi = (par%lambda/6.0_dp) * qi*qi2 - par%mu * qi + pref * (Ni*qi - s1)
      hdiag_i = 0.5_dp * par%lambda * qi2 - par%mu + pref * Ni

      st%gradV(i) = dVi
      st%force(i) = -dVi
      st%hdiag(i) = hdiag_i

      grad2 = grad2 + dVi*dVi
      lap   = lap   + (hdiag_i - pref)
    end do

    st%gradV2 = grad2
    st%lapV   = lap

    st%K = 0.5_dp * sum(st%pi * st%pi)

    W       = par%Etot - st%V
    st%jacW = W

    if (W > 0.0_dp) then
      Rreg       = 2.0_dp * W * lap - real(N - 6, dp) * grad2
      st%jacRreg = Rreg
      st%jacR    = real(N - 1, dp) * Rreg / (4.0_dp * W * W)
    else
      st%jacRreg = - real(N - 6, dp) * grad2
      st%jacR    = huge(1.0_dp)
    end if
  end subroutine initialize_md_geometry


  subroutine init_bins_autorange(par, st, bins, n_warm, frac, &
    nbins_chi, chi_log, &
    nbins_rj, nbins_rreg, nbins_gradV2, nbins_lapV)
use mod_kinds, only: dp
use mod_types_phi4, only: Phi4Params, Phi4State, Phi4Bins
implicit none

type(Phi4Params), intent(in)    :: par
type(Phi4State),  intent(inout) :: st
type(Phi4Bins),   intent(inout) :: bins

integer,  intent(in) :: n_warm
real(dp), intent(in) :: frac
integer,  intent(in) :: nbins_chi, nbins_rj, nbins_rreg, nbins_gradV2, nbins_lapV
logical,  intent(in) :: chi_log

integer  :: iw
real(dp) :: chi_min, chi_max
real(dp) :: rj_min,  rj_max
real(dp) :: rr_min,  rr_max
real(dp) :: g_min,   g_max
real(dp) :: l_min,   l_max
real(dp) :: dchi, drj, drr, dg, dl

! ---------- init ranges ----------
chi_min = huge(1.0_dp);  chi_max = -huge(1.0_dp)
rj_min  = huge(1.0_dp);  rj_max  = -huge(1.0_dp)
rr_min  = huge(1.0_dp);  rr_max  = -huge(1.0_dp)
g_min   = huge(1.0_dp);  g_max   = -huge(1.0_dp)
l_min   = huge(1.0_dp);  l_max   = -huge(1.0_dp)

! ---------- warm-up loop ----------
do iw = 1, n_warm
call vv_step_md_tangent_jlc(par, st)

! chi = W = E - V
if (st%jacW > -huge(1.0_dp) .and. st%jacW < huge(1.0_dp)) then
chi_min = min(chi_min, st%jacW)
chi_max = max(chi_max, st%jacW)
end if

! R_J
if (st%jacR > -huge(1.0_dp) .and. st%jacR < huge(1.0_dp)) then
rj_min = min(rj_min, st%jacR)
rj_max = max(rj_max, st%jacR)
end if

! K2
if (st%k2_jv > -huge(1.0_dp) .and. st%k2_jv < huge(1.0_dp)) then
rr_min = min(rr_min, st%k2_jv)
rr_max = max(rr_max, st%k2_jv)
end if

! gradV2
if (st%gradV2 > -huge(1.0_dp) .and. st%gradV2 < huge(1.0_dp)) then
g_min = min(g_min, st%gradV2)
g_max = max(g_max, st%gradV2)
end if

! lapV
if (st%lapV > -huge(1.0_dp) .and. st%lapV < huge(1.0_dp)) then
l_min = min(l_min, st%lapV)
l_max = max(l_max, st%lapV)
end if
end do

! ---------- safety fallback if something never updated ----------
if (chi_min > chi_max) then
chi_min = 1.0e-12_dp
chi_max = max(1.0_dp, abs(par%Etot))
end if
if (rj_min > rj_max) then
rj_min = -1.0_dp; rj_max = 1.0_dp
end if
if (rr_min > rr_max) then
rr_min = -1.0_dp; rr_max = 1.0_dp
end if
if (g_min > g_max) then
g_min = 0.0_dp; g_max = 1.0_dp
end if
if (l_min > l_max) then
l_min = -1.0_dp; l_max = 1.0_dp
end if

! ---------- enlarge by margin (energy-adapted via warmup) ----------

! chi = W >= 0 : use additive margin, clamp to small positive
dchi = max(chi_max - chi_min, 1.0e-12_dp)
chi_min = max(1.0e-14_dp, chi_min - frac*dchi)
chi_max = chi_max + frac*dchi

! R_J : symmetric range around 0 if it crosses 0 (more stable across energy)
if (rj_min < 0.0_dp .and. rj_max > 0.0_dp) then
  drj = max(abs(rj_min), abs(rj_max))
  drj = max(drj, 1.0e-12_dp)
  rj_min = - (1.0_dp + frac) * drj
  rj_max =   (1.0_dp + frac) * drj
else
  drj = max(rj_max - rj_min, 1.0e-12_dp)
  rj_min = rj_min - frac*drj
  rj_max = rj_max + frac*drj
end if

! RR slot (jacRreg or k2_jv): symmetric if crosses 0
if (rr_min < 0.0_dp .and. rr_max > 0.0_dp) then
  drr = max(abs(rr_min), abs(rr_max))
  drr = max(drr, 1.0e-12_dp)
  rr_min = - (1.0_dp + frac) * drr
  rr_max =   (1.0_dp + frac) * drr
else
  drr = max(rr_max - rr_min, 1.0e-12_dp)
  rr_min = rr_min - frac*drr
  rr_max = rr_max + frac*drr
end if

! gradV2 >= 0 : clamp to 0 and use additive margin
dg = max(g_max - g_min, 1.0e-12_dp)
g_min = max(0.0_dp, g_min - frac*dg)
g_max = g_max + frac*dg

! lapV : often changes sign -> symmetric if crosses 0
if (l_min < 0.0_dp .and. l_max > 0.0_dp) then
  dl = max(abs(l_min), abs(l_max))
  dl = max(dl, 1.0e-12_dp)
  l_min = - (1.0_dp + frac) * dl
  l_max =   (1.0_dp + frac) * dl
else
  dl = max(l_max - l_min, 1.0e-12_dp)
  l_min = l_min - frac*dl
  l_max = l_max + frac*dl
end if


! ---------- build bins using learned ranges ----------
call init_bins(bins, &
nbins_chi   = nbins_chi,  chi_min = chi_min, chi_max = chi_max, chi_log = chi_log, &
nbins_rj    = nbins_rj,   rj_min  = rj_min,  rj_max  = rj_max, &
nbins_rreg  = nbins_rreg, rreg_min= rr_min,  rreg_max= rr_max, &
nbins_gradV2= nbins_gradV2, gradV2_min = g_min, gradV2_max = g_max, &
nbins_lapV  = nbins_lapV, lapV_min = l_min, lapV_max = l_max )

! IMPORTANT: start production with empty histograms
call reset_bins(bins)

end subroutine init_bins_autorange

end module mod_phi4_init