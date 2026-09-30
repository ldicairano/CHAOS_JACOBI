module mod_metropolis_phi4

    ! JLC in natural coordinates, rewritten in physical time t.
    !
    ! Inputs:
    !   J(:)    = Jacobi field components
    !   Jdot(:) = dJ/dt
    !
    ! Output:
    !   accJ(:) = d^2 J / dt^2
    !
    ! Formula used:
    !
    ! accJ_i =
    !   [ v_i (gradV·Jdot) - gradV_i (v·Jdot) ] / W
    ! + (gradV·J)/W^2 * [ v_i (gradV·v) - 1/2 gradV_i v^2 ]
    ! + [ v_i (J·H v) - 1/2 v^2 (H J)_i ] / W
    !
    ! ==============================
    ! One full MD step:
    !   physical dynamics
    !   tangent dynamics
    !   JLC in natural coordinates (time-t form)
    !
    ! Scheme:
    !   - evaluate at phi_n
    !   - tangent half-kick
    !   - JLC half-kick
    !   - physical half-kick
    !   - drift phi, dphi, jphi
    !   - evaluate at phi_{n+1}
    !   - physical half-kick
    !   - tangent half-kick
    !   - JLC half-kick


  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4State, Phi4Obs, RNGState
  implicit none

  public :: vv_step_md_tangent_jlc

contains

  subroutine vv_step_md_tangent_jlc(par, st)
    use mod_kinds,      only: dp
    use mod_types_phi4, only: Phi4Params, Phi4State
    implicit none
  
    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(inout) :: st
  
    integer  :: i, N
    real(dp) :: dt, hdt, Ni, pref
  
    ! --- old-state contractions ---
    real(dp) :: sumd_old, sumJ_old, sumpi_old
    real(dp) :: v2_old
    real(dp) :: dot_g_J_old, dot_g_v_old, dot_g_Jdot_old, dot_v_Jdot_old, dot_J_Hv_old
    real(dp) :: Hj_i, accJ_i, diag_i
  
    ! --- new-state raw sums ---
    real(dp) :: s1_new, s2_new, s4_new
    real(dp) :: qi, qi2, qi4
  
    ! --- new-state full evaluation + second half-kick ---
    real(dp) :: V_new, grad2_new, lap_new, W_new, Rreg_new
    real(dp) :: sumd_new, sumJ_new, sumpi_new
    real(dp) :: v2_new
    real(dp) :: dot_g_J_new, dot_g_v_new, dot_g_Jdot_new, dot_v_Jdot_new, dot_J_Hv_new
    real(dp) :: dot_J_v_new, dot_v_Hv_new
    real(dp) :: alpha_jv, sumJperp_new, dot_g_Jperp_new, dot_Jperp_Hv_new
    real(dp) :: jperp_i, Hjperp_i, accJ_geom_i
    real(dp) :: JJperp, JaccGeom
    real(dp) :: dVi, hdiag_i
    real(dp) :: K_new
    real(dp) :: norm_jlc, norm_tan
    
    N    = par%N
    Ni   = real(N, dp)
    dt   = par%dt_md
    hdt  = 0.5_dp * dt
    pref = 4.0_dp * par%coup / real(N - 1, dp)
  
    !===========================================================
    ! LOOP 1: old-state contractions from cached arrays
    !===========================================================
    sumd_old       = 0.0_dp
    sumJ_old       = 0.0_dp
    sumpi_old      = 0.0_dp
    v2_old         = 0.0_dp
    dot_g_J_old    = 0.0_dp
    dot_g_v_old    = 0.0_dp
    dot_g_Jdot_old = 0.0_dp
    dot_v_Jdot_old = 0.0_dp
    dot_J_Hv_old   = 0.0_dp
  
    do i = 1, N
      sumd_old       = sumd_old       + st%dphi(i)
      sumJ_old       = sumJ_old       + st%jphi(i)
      sumpi_old      = sumpi_old      + st%pi(i)
      v2_old         = v2_old         + st%pi(i)    * st%pi(i)
      dot_g_J_old    = dot_g_J_old    + st%gradV(i) * st%jphi(i)
      dot_g_v_old    = dot_g_v_old    + st%gradV(i) * st%pi(i)
      dot_g_Jdot_old = dot_g_Jdot_old + st%gradV(i) * st%jpi(i)
      dot_v_Jdot_old = dot_v_Jdot_old + st%pi(i)    * st%jpi(i)
      dot_J_Hv_old   = dot_J_Hv_old   + st%hdiag(i) * st%jphi(i) * st%pi(i)
    end do
  
    dot_J_Hv_old = dot_J_Hv_old - pref * sumJ_old * sumpi_old
  
    !===========================================================
    ! LOOP 2: first half-kicks + drift
    !===========================================================
    do i = 1, N
  
      ! tangent half-kick
      diag_i = st%hdiag(i)
      st%dpi_t(i) = st%dpi_t(i) + hdt * ( -diag_i * st%dphi(i) + pref * sumd_old )
  
        Hj_i = st%hdiag(i) * st%jphi(i) - pref * sumJ_old
        accJ_i = ( st%pi(i) * dot_g_Jdot_old - st%gradV(i) * dot_v_Jdot_old ) / st%jacW &
               + ( dot_g_J_old / (st%jacW*st%jacW) ) * ( st%pi(i) * dot_g_v_old - 0.5_dp * st%gradV(i) * v2_old ) &
               + ( st%pi(i) * dot_J_Hv_old - 0.5_dp * v2_old * Hj_i ) / st%jacW

      st%jpi(i) = st%jpi(i) + hdt * accJ_i
  
      ! physical half-kick
      st%pi(i) = st%pi(i) + hdt * st%force(i)
  
      ! drift
      st%phi(i)  = st%phi(i)  + dt * st%pi(i)
      st%dphi(i) = st%dphi(i) + dt * st%dpi_t(i)
      st%jphi(i) = st%jphi(i) + dt * st%jpi(i)
  
    end do
  
    !===========================================================
    ! LOOP 3: new-state raw sums
    !===========================================================
    s1_new = 0.0_dp
    s2_new = 0.0_dp
    s4_new = 0.0_dp
  
    do i = 1, N
      qi  = st%phi(i)
      qi2 = qi*qi
      qi4 = qi2*qi2
  
      s1_new = s1_new + qi
      s2_new = s2_new + qi2
      s4_new = s4_new + qi4
    end do
  
    st%S1   = s1_new
    st%S2   = s2_new
    st%S4   = s4_new
    st%M    = s1_new / Ni
    st%phi2 = s2_new / Ni
    st%phi4 = s4_new / Ni
  
    V_new = (par%lambda/24.0_dp) * s4_new - 0.5_dp * par%mu * s2_new &
          + 2.0_dp * par%coup / real(N - 1, dp) * (Ni*s2_new - s1_new*s1_new)
  
    !===========================================================
    ! LOOP 4: new-state evaluation + second half-kicks + K
    !===========================================================
    grad2_new      = 0.0_dp
    lap_new        = 0.0_dp
    sumd_new       = 0.0_dp
    sumJ_new       = 0.0_dp
    sumpi_new      = 0.0_dp
    v2_new         = 0.0_dp
    dot_g_J_new    = 0.0_dp
    dot_g_v_new    = 0.0_dp
    dot_g_Jdot_new = 0.0_dp
    dot_v_Jdot_new = 0.0_dp
    dot_J_Hv_new   = 0.0_dp
    dot_J_v_new    = 0.0_dp
    dot_v_Hv_new   = 0.0_dp

    ! First part of loop 4: evaluate new state and all needed contractions
    do i = 1, N
      qi  = st%phi(i)
      qi2 = qi*qi
  
      dVi     = (par%lambda/6.0_dp) * qi*qi2 - par%mu * qi + pref * (Ni*qi - s1_new)
      hdiag_i = 0.5_dp * par%lambda * qi2 - par%mu + pref * Ni
  
      st%gradV(i) = dVi
      st%force(i) = -dVi
      st%hdiag(i) = hdiag_i
  
      grad2_new      = grad2_new      + dVi * dVi
      lap_new        = lap_new        + (hdiag_i - pref)
      sumd_new       = sumd_new       + st%dphi(i)
      sumJ_new       = sumJ_new       + st%jphi(i)
      sumpi_new      = sumpi_new      + st%pi(i)
      v2_new         = v2_new         + st%pi(i)    * st%pi(i)
      dot_g_J_new    = dot_g_J_new    + dVi         * st%jphi(i)
      dot_g_v_new    = dot_g_v_new    + dVi         * st%pi(i)
      dot_g_Jdot_new = dot_g_Jdot_new + dVi         * st%jpi(i)
      dot_v_Jdot_new = dot_v_Jdot_new + st%pi(i)    * st%jpi(i)
      dot_J_Hv_new   = dot_J_Hv_new   + hdiag_i     * st%jphi(i) * st%pi(i)
      dot_J_v_new    = dot_J_v_new    + st%jphi(i) * st%pi(i)
      dot_v_Hv_new   = dot_v_Hv_new   + hdiag_i    * st%pi(i) * st%pi(i)
    end do
  
    dot_J_Hv_new = dot_J_Hv_new - pref * sumJ_new * sumpi_new
    dot_v_Hv_new = dot_v_Hv_new - pref * sumpi_new * sumpi_new

    alpha_jv = dot_J_v_new / v2_new

    sumJperp_new     = sumJ_new     - alpha_jv * sumpi_new
    dot_g_Jperp_new  = dot_g_J_new  - alpha_jv * dot_g_v_new
    dot_Jperp_Hv_new = dot_J_Hv_new - alpha_jv * dot_v_Hv_new
  
    st%V      = V_new
    st%gradV2 = grad2_new
    st%lapV   = lap_new
  
    W_new   = par%Etot - V_new
    st%jacW = W_new

    st%s_jac = st%s_jac + 2.0_dp * W_new * par%dt_md
  
    st%d_bdry_jac = (2.0_dp*sqrt(2.0_dp)/3.0_dp) * st%jacW**1.5_dp / sqrt(st%gradV2)

      Rreg_new   = 2.0_dp * W_new * lap_new - real(N - 6, dp) * grad2_new
      st%jacRreg = Rreg_new
      st%jacR    = real(N - 1, dp) * Rreg_new / (4.0_dp * W_new * W_new)
  
      K_new    = 0.0_dp
      JJperp   = 0.0_dp
      JaccGeom = 0.0_dp
  
    ! Second part of loop 4: second half-kicks + kinetic energy
    do i = 1, N

      ! tangent second half-kick
      st%dpi_t(i) = st%dpi_t(i) + hdt * ( -st%hdiag(i) * st%dphi(i) + pref * sumd_new )
  
      ! JLC second half-kick
      ! IMPORTANT:
      ! here we still use the intermediate velocity pi^{n+1/2},
      ! consistent with v2_new, dot_g_v_new, dot_J_Hv_new, etc.
      Hj_i = st%hdiag(i) * st%jphi(i) - pref * sumJ_new
      accJ_i = ( st%pi(i) * dot_g_Jdot_new - st%gradV(i) * dot_v_Jdot_new ) / W_new &
             + ( dot_g_J_new / (W_new*W_new) ) * ( st%pi(i) * dot_g_v_new - 0.5_dp * st%gradV(i) * v2_new ) &
             + ( st%pi(i) * dot_J_Hv_new - 0.5_dp * v2_new * Hj_i ) / W_new

      ! --- sectional curvature K^(2)(J,v): use only the geometric J-part, with J projected orthogonal to v
      jperp_i  = st%jphi(i) - alpha_jv * st%pi(i)
      Hjperp_i = st%hdiag(i) * jperp_i - pref * sumJperp_new

      accJ_geom_i = ( dot_g_Jperp_new / (W_new*W_new) ) * &
                    ( st%pi(i) * dot_g_v_new - 0.5_dp * st%gradV(i) * v2_new ) &
                  + ( st%pi(i) * dot_Jperp_Hv_new - 0.5_dp * v2_new * Hjperp_i ) / W_new

      JJperp   = JJperp   + jperp_i * jperp_i
      JaccGeom = JaccGeom + jperp_i * accJ_geom_i

      st%jpi(i) = st%jpi(i) + hdt * accJ_i
  
      ! physical second half-kick
      st%pi(i) = st%pi(i) + hdt * st%force(i)
  
      ! kinetic energy after full second half-kick
      K_new = K_new + 0.5_dp * st%pi(i) * st%pi(i)
  
    end do
  
    st%K = K_new

    st%k2_jv = - JaccGeom / (4.0_dp * W_new * W_new * JJperp)
  
    norm_jlc = sqrt(sum(st%jphi*st%jphi) + sum(st%jpi*st%jpi))
    st%jlc_logsum = st%jlc_logsum + log(norm_jlc)
    st%jphi = st%jphi / norm_jlc
    st%jpi  = st%jpi  / norm_jlc
    st%jlc_time = st%jlc_time + par%dt_md

    norm_tan = sqrt(sum(st%dphi*st%dphi) + sum(st%dpi_t*st%dpi_t))
    st%lya_logsum = st%lya_logsum + log(norm_tan)
    st%dphi  = st%dphi  / norm_tan
    st%dpi_t = st%dpi_t / norm_tan
    st%lya_time = st%lya_time + par%dt_md

  end subroutine vv_step_md_tangent_jlc


end module mod_metropolis_phi4