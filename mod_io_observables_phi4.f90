module mod_io_observables_phi4
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4State, Phi4Obs, IOParams
  implicit none
  private

  public :: reset_obs, update_obs, reinflate

contains


  subroutine reset_obs(obs)
    type(Phi4Obs), intent(inout) :: obs
    obs = Phi4Obs()
  end subroutine reset_obs


  subroutine update_obs(par, st, obs, n_MC)
    type(Phi4Params), intent(in)    :: par
    type(Phi4State),  intent(in)    :: st
    type(Phi4Obs),    intent(inout) :: obs
    integer,          intent(in)    :: n_MC

    real(dp) :: Nf, kin_per_dof
    real(dp) :: m, e, a1, a2, a3
    real(dp) :: lam_tan, lam_jlc

    Nf = real(par%N, dp)

    m = st%M
    e = (st%K + st%V) / real(par%N, dp)
    kin_per_dof = st%K / Nf

    !--------------------------------------------------
    ! Standard thermo / order-parameter accumulators
    !--------------------------------------------------
    obs%av_mag  = obs%av_mag  + m
    obs%av_mag2 = obs%av_mag2 + m*m
    obs%av_mag4 = obs%av_mag4 + m*m*m*m

    obs%av_pot  = obs%av_pot + st%V 
    obs%av_kin  = obs%av_kin + st%K 

    obs%av_en   = obs%av_en  + e
    obs%av_en2  = obs%av_en2 + e*e

    obs%binder = 1.0_dp - obs%c_mag4 / (3.0_dp * obs%c_mag2 * obs%c_mag2)

    obs%c_mag  = obs%av_mag  / real(n_MC, dp)
    obs%c_mag2 = obs%av_mag2 / real(n_MC, dp)
    obs%c_mag4 = obs%av_mag4 / real(n_MC, dp)

    obs%c_pot = obs%av_pot / real(n_MC, dp)
    obs%c_kin = obs%av_kin / real(n_MC, dp)

    obs%c_en  = obs%av_en  / real(n_MC, dp)
    obs%c_en2 = obs%av_en2 / real(n_MC, dp)

    !--------------------------------------------------
    ! Microcanonical kinetic estimators
    !--------------------------------------------------
      obs%av_kin_1 = obs%av_kin_1 + 1.0_dp / kin_per_dof
      obs%av_kin_2 = obs%av_kin_2 + 1.0_dp / (kin_per_dof*kin_per_dof)
      obs%av_kin_3 = obs%av_kin_3 + 1.0_dp / (kin_per_dof*kin_per_dof*kin_per_dof)

      obs%c_kin_1 = obs%av_kin_1 / real(n_MC, dp)
      obs%c_kin_2 = obs%av_kin_2 / real(n_MC, dp)
      obs%c_kin_3 = obs%av_kin_3 / real(n_MC, dp)

    !--------------------------------------------------
    ! Microcanonical beta and entropy derivatives
    ! Now with N momentum dof, not N-1
    !--------------------------------------------------
      a1 = 0.5_dp - 1.0_dp/Nf
      a2 = 0.5_dp - 2.0_dp/Nf
      a3 = 0.5_dp - 3.0_dp/Nf
  
      obs%beta = a1 * obs%c_kin_1
  
      obs%ders_2 = Nf * ( a1*a2*obs%c_kin_2 - (a1*a1)*(obs%c_kin_1*obs%c_kin_1) )
  
      obs%ders_3 = Nf*Nf * ( &
           a1*a2*a3*obs%c_kin_3 &
         - 3.0_dp*(a1*a1)*a2*obs%c_kin_2*obs%c_kin_1 &
         + 2.0_dp*(a1*a1*a1)*(obs%c_kin_1*obs%c_kin_1*obs%c_kin_1) )

    !--------------------------------------------------
    ! Geometry: global averages
    !--------------------------------------------------
    obs%av_jacW    = obs%av_jacW    + st%jacW
    obs%av_jacW2   = obs%av_jacW2   + st%jacW * st%jacW

    obs%av_gradV2  = obs%av_gradV2  + st%gradV2
    obs%av_gradV22 = obs%av_gradV22 + st%gradV2 * st%gradV2

    obs%av_lapV    = obs%av_lapV    + st%lapV
    obs%av_lapV2   = obs%av_lapV2   + st%lapV * st%lapV

    obs%av_jacR    = obs%av_jacR    + st%jacR
    obs%av_jacR2   = obs%av_jacR2   + st%jacR * st%jacR

    obs%av_jacRreg  = obs%av_jacRreg  + st%jacRreg
    obs%av_jacRreg2 = obs%av_jacRreg2 + st%jacRreg * st%jacRreg

    !--------------------------------------------------
    ! Current Lyapunov estimates
    ! These are NOT averaged in time here: they are read
    ! from the running Benettin/JLC accumulators in st.
    !--------------------------------------------------

      lam_tan = st%lya_logsum / st%lya_time
      lam_jlc = st%jlc_logsum / st%jlc_time

    !--------------------------------------------------
    ! Running means
    !--------------------------------------------------




    obs%c_jacW    = obs%av_jacW    / real(n_MC, dp)
    obs%c_jacW2   = obs%av_jacW2   / real(n_MC, dp)

    obs%c_gradV2  = obs%av_gradV2  / real(n_MC, dp)
    obs%c_gradV22 = obs%av_gradV22 / real(n_MC, dp)

    obs%c_lapV    = obs%av_lapV    / real(n_MC, dp)
    obs%c_lapV2   = obs%av_lapV2   / real(n_MC, dp)

    obs%c_jacR    = obs%av_jacR    / real(n_MC, dp)
    obs%c_jacR2   = obs%av_jacR2   / real(n_MC, dp)

    obs%c_jacRreg  = obs%av_jacRreg  / real(n_MC, dp)
    obs%c_jacRreg2 = obs%av_jacRreg2 / real(n_MC, dp)

    obs%var_jacW    = obs%c_jacW2    - obs%c_jacW    * obs%c_jacW
    obs%var_gradV2  = obs%c_gradV22  - obs%c_gradV2  * obs%c_gradV2
    obs%var_lapV    = obs%c_lapV2    - obs%c_lapV    * obs%c_lapV
    obs%var_jacR    = obs%c_jacR2    - obs%c_jacR    * obs%c_jacR
    obs%var_jacRreg = obs%c_jacRreg2 - obs%c_jacRreg * obs%c_jacRreg

    obs%c_lya_tan = lam_tan
    obs%c_lya_jlc = lam_jlc

    obs%n_MC = n_MC
  end subroutine update_obs


  subroutine reinflate(obs)
    type(Phi4Obs), intent(inout) :: obs
    real(dp) :: n

    n = real(obs%n_MC, dp)

    obs%av_mag   = obs%c_mag  * n
    obs%av_mag2  = obs%c_mag2 * n
    obs%av_mag4  = obs%c_mag4 * n

    obs%av_pot   = obs%c_pot  * n
    obs%av_kin   = obs%c_kin  * n

    obs%av_en    = obs%c_en   * n
    obs%av_en2   = obs%c_en2  * n

    obs%av_kin_1 = obs%c_kin_1 * n
    obs%av_kin_2 = obs%c_kin_2 * n
    obs%av_kin_3 = obs%c_kin_3 * n

    obs%av_jacW     = obs%c_jacW     * n
    obs%av_jacW2    = obs%c_jacW2    * n
    obs%av_gradV2   = obs%c_gradV2   * n
    obs%av_gradV22  = obs%c_gradV22  * n
    obs%av_lapV     = obs%c_lapV     * n
    obs%av_lapV2    = obs%c_lapV2    * n
    obs%av_jacR     = obs%c_jacR     * n
    obs%av_jacR2    = obs%c_jacR2    * n
    obs%av_jacRreg  = obs%c_jacRreg  * n
    obs%av_jacRreg2 = obs%c_jacRreg2 * n

    obs%av_lya_tan  = obs%c_lya_tan * n
    obs%av_lya_jlc  = obs%c_lya_jlc * n
  end subroutine reinflate

end module mod_io_observables_phi4