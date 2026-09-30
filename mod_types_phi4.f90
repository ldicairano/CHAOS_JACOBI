module mod_types_phi4
  use mod_kinds, only: dp
  implicit none

  type :: Phi4Params
    integer :: L = 0
    integer :: N = 0

    ! --- legacy sampling id (for obs_*.dat naming) ---
    integer  :: n_samp = 0
    real(dp) :: en_in  = 0.0_dp
    real(dp) :: Etot   = 0.0_dp   ! total energy for microcanonical MD

    ! --- run control ---
    integer :: n_steps  = 0
    integer :: n_jump   = 1
    integer :: n_therm  = 0
    integer :: n_realiz = 1
    integer :: n_rest   = 1
    integer :: n_print  = 1000
    real(dp) :: cluster_time = 172800.0_dp

    logical :: time_dyn = .false.
    integer :: n_print_inst = 0

    ! --- model params ---
    real(dp) :: coup   = 1.0_dp
    real(dp) :: mu     = 1.0_dp
    real(dp) :: lambda = 6.0_dp
    logical  :: auto_coup = .true.

    ! --- MD control ---
    real(dp) :: dt_md = 1.0e-2_dp

    ! --- RNG ---
    integer :: seed = 12345
  end type Phi4Params


  type :: Phi4State
    ! Physical state
    real(dp), allocatable :: phi(:)
    real(dp), allocatable :: pi(:)
    real(dp), allocatable :: force(:)

    ! Tangent dynamics: (delta phi, delta pi)
    real(dp), allocatable :: dphi(:)
    real(dp), allocatable :: dpi_t(:)

    ! JLC field: (J_phi, J_pi)
    ! For now these are just allocated/initialized consistently.
    ! The actual JLC evolution will be implemented in the MD module.
    real(dp), allocatable :: jphi(:)
    real(dp), allocatable :: jpi(:)

    real(dp), allocatable :: gradV(:)
    real(dp), allocatable :: hdiag(:)

    ! Energies
    real(dp) :: V = 0.0_dp
    real(dp) :: K = 0.0_dp

    ! Raw sums
    real(dp) :: S1 = 0.0_dp
    real(dp) :: S2 = 0.0_dp
    real(dp) :: S4 = 0.0_dp

    ! Normalized macros
    real(dp) :: M    = 0.0_dp
    real(dp) :: phi2 = 0.0_dp
    real(dp) :: phi4 = 0.0_dp

    ! Jacobi / geometry diagnostics
    real(dp) :: jacW    = 0.0_dp   ! W = E - V  (micro) or generic effective value if needed
    real(dp) :: gradV2  = 0.0_dp
    real(dp) :: lapV    = 0.0_dp
    real(dp) :: jacR    = 0.0_dp
    real(dp) :: jacRreg = 0.0_dp

    ! Running logs for Lyapunov/JLC renormalization
    real(dp) :: lya_logsum = 0.0_dp
    real(dp) :: jlc_logsum = 0.0_dp
    integer  :: lya_count  = 0
    integer  :: jlc_count  = 0

    real(dp) :: lya_time   = 0.0_dp
    real(dp) :: jlc_time   = 0.0_dp

    real(dp) :: s_jac = 0.0_dp

    real(dp) :: d_bdry_jac = 0.0_dp

    real(dp) :: k2_jv = 0.0_dp

  end type Phi4State


  type :: Phi4Obs
    integer :: n_acc = 0
    integer :: n_att = 0
    real(dp) :: acc_rate = 0.0_dp

    integer :: n_MC = 0

    ! ---------------------------------------------------------
    ! Raw accumulated sums
    ! ---------------------------------------------------------
    real(dp) :: av_mag   = 0.0_dp
    real(dp) :: av_mag2  = 0.0_dp
    real(dp) :: av_mag4  = 0.0_dp

    real(dp) :: av_pot   = 0.0_dp
    real(dp) :: av_kin   = 0.0_dp
    real(dp) :: av_kin_1 = 0.0_dp
    real(dp) :: av_kin_2 = 0.0_dp
    real(dp) :: av_kin_3 = 0.0_dp

    real(dp) :: av_en    = 0.0_dp
    real(dp) :: av_en2   = 0.0_dp

    ! ---------------------------------------------------------
    ! Running means / cumulants
    ! ---------------------------------------------------------
    real(dp) :: c_mag   = 0.0_dp
    real(dp) :: c_mag2  = 0.0_dp
    real(dp) :: c_mag4  = 0.0_dp
    real(dp) :: binder  = 0.0_dp

    real(dp) :: c_pot   = 0.0_dp
    real(dp) :: c_kin   = 0.0_dp
    real(dp) :: c_kin_1 = 0.0_dp
    real(dp) :: c_kin_2 = 0.0_dp
    real(dp) :: c_kin_3 = 0.0_dp

    real(dp) :: c_en    = 0.0_dp
    real(dp) :: c_en2   = 0.0_dp
    real(dp) :: cv      = 0.0_dp

    ! ---------------------------------------------------------
    ! Microcanonical derivatives
    ! ---------------------------------------------------------
    real(dp) :: beta   = 0.0_dp
    real(dp) :: ders_2 = 0.0_dp
    real(dp) :: ders_3 = 0.0_dp

    ! ---------------------------------------------------------
    ! Extra magnetization observables (legacy / optional)
    ! ---------------------------------------------------------
    real(dp) :: av_m   = 0.0_dp
    real(dp) :: av_m2  = 0.0_dp
    real(dp) :: c_m    = 0.0_dp
    real(dp) :: c_m2   = 0.0_dp
    real(dp) :: chi    = 0.0_dp

    ! ---------------------------------------------------------
    ! Geometric observables: accumulated sums
    ! ---------------------------------------------------------
    real(dp) :: av_jacW     = 0.0_dp
    real(dp) :: av_jacW2    = 0.0_dp

    real(dp) :: av_gradV2   = 0.0_dp
    real(dp) :: av_gradV22  = 0.0_dp

    real(dp) :: av_lapV     = 0.0_dp
    real(dp) :: av_lapV2    = 0.0_dp

    real(dp) :: av_jacR     = 0.0_dp
    real(dp) :: av_jacR2    = 0.0_dp

    real(dp) :: av_jacRreg  = 0.0_dp
    real(dp) :: av_jacRreg2 = 0.0_dp

    ! ---------------------------------------------------------
    ! Geometric observables: running means
    ! ---------------------------------------------------------
    real(dp) :: c_jacW     = 0.0_dp
    real(dp) :: c_jacW2    = 0.0_dp

    real(dp) :: c_gradV2   = 0.0_dp
    real(dp) :: c_gradV22  = 0.0_dp

    real(dp) :: c_lapV     = 0.0_dp
    real(dp) :: c_lapV2    = 0.0_dp

    real(dp) :: c_jacR     = 0.0_dp
    real(dp) :: c_jacR2    = 0.0_dp

    real(dp) :: c_jacRreg  = 0.0_dp
    real(dp) :: c_jacRreg2 = 0.0_dp

    ! ---------------------------------------------------------
    ! Geometric variances
    ! ---------------------------------------------------------
    real(dp) :: var_jacW    = 0.0_dp
    real(dp) :: var_gradV2  = 0.0_dp
    real(dp) :: var_lapV    = 0.0_dp
    real(dp) :: var_jacR    = 0.0_dp
    real(dp) :: var_jacRreg = 0.0_dp

    ! ---------------------------------------------------------
    ! Lyapunov estimators
    ! ---------------------------------------------------------
    real(dp) :: av_lya_tan = 0.0_dp
    real(dp) :: av_lya_jlc = 0.0_dp

    real(dp) :: c_lya_tan  = 0.0_dp
    real(dp) :: c_lya_jlc  = 0.0_dp
  end type Phi4Obs


  type :: Phi4Bins
  ! =========================================================
  ! Histogram of chi = E - V
  ! =========================================================
  integer :: nbins_chi = 0
  real(dp) :: chi_min = 0.0_dp
  real(dp) :: chi_max = 0.0_dp
  logical  :: chi_log = .false.

  integer :: nsamples_chi = 0
  integer :: n_bad_chi    = 0
  integer,  allocatable :: count_chi(:)

  ! Conditional averages in bins of chi
  real(dp), allocatable :: sum_Rj_vs_chi(:)
  real(dp), allocatable :: sum_Rj2_vs_chi(:)

  real(dp), allocatable :: sum_Rreg_vs_chi(:)
  real(dp), allocatable :: sum_Rreg2_vs_chi(:)

  real(dp), allocatable :: sum_gradV2_vs_chi(:)
  real(dp), allocatable :: sum_gradV22_vs_chi(:)

  real(dp), allocatable :: sum_lapV_vs_chi(:)
  real(dp), allocatable :: sum_lapV2_vs_chi(:)

  ! =========================================================
  ! Histogram of R_J
  ! =========================================================
  integer :: nbins_rj = 0
  real(dp) :: rj_min = 0.0_dp
  real(dp) :: rj_max = 0.0_dp

  integer :: nsamples_rj = 0
  integer :: n_bad_rj    = 0
  integer,  allocatable :: count_rj(:)

  ! =========================================================
  ! Histogram of R_J^reg
  ! =========================================================
  integer :: nbins_rreg = 0
  real(dp) :: rreg_min = 0.0_dp
  real(dp) :: rreg_max = 0.0_dp

  integer :: nsamples_rreg = 0
  integer :: n_bad_rreg    = 0
  integer,  allocatable :: count_rreg(:)

  ! =========================================================
  ! Histogram of gradV2 = |grad V|^2
  ! =========================================================
  integer :: nbins_gradV2 = 0
  real(dp) :: gradV2_min = 0.0_dp
  real(dp) :: gradV2_max = 0.0_dp

  integer :: nsamples_gradV2 = 0
  integer :: n_bad_gradV2    = 0
  integer,  allocatable :: count_gradV2(:)

  ! =========================================================
  ! Histogram of lapV = Delta V
  ! =========================================================
  integer :: nbins_lapV = 0
  real(dp) :: lapV_min = 0.0_dp
  real(dp) :: lapV_max = 0.0_dp

  integer :: nsamples_lapV = 0
  integer :: n_bad_lapV    = 0
  integer,  allocatable :: count_lapV(:)

  real(dp) :: sum_k2_global  = 0.0_dp
  real(dp) :: sum_k22_global = 0.0_dp

end type Phi4Bins


  type :: IOParams
    character(len=256) :: out_dir = "out"
    character(len=256) :: restart_dir = "."
  end type IOParams


  type :: RNGState
    integer :: seed = 12345
  end type RNGState

end module mod_types_phi4