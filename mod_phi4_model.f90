module mod_phi4_model
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4State
  implicit none
contains

  pure real(dp) function onsite(par, phi) result(v0)
    ! Onsite potential as in the paper:
    !   V_loc(phi) = (lambda/4!) * phi^4  -  (mu^2/2) * phi^2
    ! NOTE: par%mu is interpreted as mu^2 (paper notation).
    type(Phi4Params), intent(in) :: par
    real(dp), intent(in) :: phi
    real(dp) :: phi2
    phi2 = phi*phi
    v0 = (par%lambda/24.0_dp) * (phi2*phi2) - 0.5_dp*par%mu*phi2
  end function onsite


  pure real(dp) function v_mf(par, S1, S2) result(vm)
    ! Mean-field coupling term (paper mean-field case):
    !   V_MF = (2J/(N-1)) * ( N*S2 - S1^2 )
    type(Phi4Params), intent(in) :: par
    real(dp), intent(in) :: S1, S2
    real(dp) :: Ndp
    Ndp = real(par%N, dp)
    vm  = (2.0_dp*par%coup/real(par%N-1,dp)) * ( Ndp*S2 - S1*S1 )
  end function v_mf


  subroutine compute_potential(par, st)
    ! Full potential V = sum onsite + V_MF, with coefficients/signs as in the paper.
    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(inout) :: st
    integer :: i
    real(dp) :: Vloc, S1, S2, S4
    real(dp) :: p, p2

    Vloc = 0.0_dp
    S1   = 0.0_dp
    S2   = 0.0_dp
    S4   = 0.0_dp

    do i = 1, par%N
      p  = st%phi(i)
      p2 = p*p
      Vloc = Vloc + onsite(par, p)
      S1   = S1 + p
      S2   = S2 + p2
      S4   = S4 + p2*p2
    end do

    st%S1 = S1
    st%S2 = S2
    st%S4 = S4

    st%V  = Vloc + v_mf(par, S1, S2)

    st%M    = st%S1 / real(par%N, dp)
    st%phi2 = st%S2 / real(par%N, dp)
    st%phi4 = st%S4 / real(par%N, dp)
  end subroutine compute_potential





  subroutine compute_kinetic(par, st)
    ! K = 1/2 sum pi^2 (canonical only).
    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(inout) :: st
    integer :: i
    real(dp) :: Kloc

    if (.not. allocated(st%pi)) then
      st%K = 0.0_dp
      return
    end if

    Kloc = 0.0_dp
    do i = 1, par%N
      Kloc = Kloc + 0.5_dp * st%pi(i) * st%pi(i)
    end do
    st%K = Kloc
  end subroutine compute_kinetic

end module mod_phi4_model