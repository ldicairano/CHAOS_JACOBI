program main_phi4
  use mod_kinds,              only: dp
  use mod_types_phi4,         only: Phi4Params, IOParams, RNGState, Phi4State, Phi4Obs, Phi4Bins
  use mod_cli_phi4,           only: parse_cli_phi4, cli_msg_restart_already_complete
  use mod_input_phi4,         only: read_input_file
  use mod_random,             only: rng_init
  use mod_phi4_init,          only: allocate_state, init_configuration, init_bins_autorange
  use mod_metropolis_phi4,    only: vv_step_md_tangent_jlc
  use mod_io_observables_phi4,only: reset_obs, update_obs, reinflate
  use mod_io_append_dispatch, only: print_observables, append_inst_micro, clear_output_files1, clear_output_files2
  use mod_io_restart_phi4,    only: read_restart, write_restart, read_last_nmc_from_dat
  use mod_bins_phi4,          only: init_bins, reset_bins, update_bins, write_bins
  implicit none

!               ██╗██████╗      ███╗   ███╗███████╗    ██████╗██╗  ██╗██╗██╗  ██╗
!             ████║██╔══███╗    ████╗ ████║██╔════╝    ██╔═██║██║  ██║██║██║  ██║
!             ╚═██║██║  ███║    ██╔████╔██║█████╗      ██████║███████║██║███████║
!               ██║██║  ███║    ██║╚██╔╝██║██╔══╝      ██╔═══╝██╔══██║██║╚════██║
!               ██║██████╔╝     ██║ ╚═╝ ██║██║         ██║    ██║  ██║██║     ██║
!               ╚═╝╚═════╝      ╚═╝     ╚═╝╚═╝         ╚═╝    ╚═╝  ╚═╝╚═╝     ╚═╝

!================================================================================
!
!   ██╗  ██╗ █████╗ ███╗   ███╗██╗██╗  ████████╗ ██████╗ ███╗   ██╗██╗ █████╗ ███╗  ██╗
!   ██║  ██║██╔══██╗████╗ ████║██║██║     ██╔══╝██╔═══██╗████╗  ██║██║██╔══██╗████╗ ██║
!   ███████║███████║██╔████╔██║██║██║     ██║   ██║   ██║██╔██╗ ██║██║███████║██╔██╗██║
!   ██╔══██║██╔══██║██║╚██╔╝██║██║██║     ██║   ██║   ██║██║╚██╗██║██║██╔══██║██║╚████║
!   ██║  ██║██║  ██║██║ ╚═╝ ██║██║███████╗██║   ╚██████╔╝██║ ╚████║██║██║  ██║██║ ╚███║
!   ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝╚══════╝╚═╝    ╚═════╝ ╚═╝  ╚═══╝╚═╝╚═╝  ╚═╝╚═╝  ╚══╝
!
!                        ██████╗██╗  ██╗ █████╗  ██████╗ ███████╗                             
!                       ██╔════╝██║  ██║██╔══██╗██╔═══██╗██╔════╝
!                       ██║     ███████║███████║██║   ██║███████╗
!                       ██║     ██╔══██║██╔══██║██║   ██║╚════██║
!                      ╚██████╗ ██║  ██║██║  ██║╚██████╔╝███████║
!                        ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝ ╚══════╝
!
!                                      AND
!
!                       ██╗ █████╗  ██████╗ ██████╗ ██████╗ ██╗
!                       ██║██╔══██╗██╔════╝██╔═══██╗██╔══██╗██║
!                       ██║███████║██║     ██║   ██║██████╔╝██║
!                  ██   ██║██╔══██║██║     ██║   ██║██╔══██╗██║
!                  ╚█████╔╝██║  ██║╚██████╗╚██████╔╝██████╔╝██║
!                   ╚════╝ ╚═╝  ╚═╝ ╚═════╝ ╚═════╝ ╚═════╝ ╚═╝
!
!                  ███╗   ███╗███████╗████████╗██████╗ ██╗ ██████╗
!                  ████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██║██╔════╝
!                  ██╔████╔██║█████╗     ██║   ██████╔╝██║██║     
!                  ██║╚██╔╝██║██╔══╝     ██║   ██╔══██╗██║██║     
!                  ██║ ╚═╝ ██║███████╗   ██║   ██║  ██║██║╚██████╗
!                  ╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝ ╚═════╝
!
!================================================================================

!                                  ϕ^4   (MEAN-FIELD, 1D)
!                                   Luxembourg, 3 Mar 2026
!  //============================================================================\\
!  |                        1D  MEAN-FIELD  ϕ^4  MODEL                            |
!  |                     (CHAIN  +  ALL-TO-ALL  TERM)                             |
!  |                    MOLECULAR DYNAMICS MD                                     |
!  |                    TANGENT DYNAMICS AND JACOBI-LEVI CIVITA EQUATION          |
!  \\===========================================================================//

! ------------------------------------------------
!  Hamiltonian (typical form):
!
! H = sum_i [ p_i^2/2 ] +
!     sum_i [ (lambda/4!) phi_i^4 - (mu^2/2) phi_i^2 ] +
!     (2J/(N-1)) * ( N*sum_i phi_i^2 - (sum_i phi_i)^2 )    (all-to-all mean-field)
!      NOTE: in code par%mu is mu^2 (paper notation).
!

  type(Phi4Params) :: par
  type(IOParams)   :: io
  type(RNGState)   :: rng
  type(Phi4State)  :: st
  type(Phi4State) :: st_saved
  type(Phi4Obs)    :: obs
  type(Phi4Bins)   :: bins

  character(len=256) :: input_file
  integer :: n_samp_cli, L_cli, n_steps_cli, n_jump_cli, n_realiz_cli, n_rest_cli
  integer :: step, step0, n_print, n_sweep
  integer :: n_MC
  integer :: nmc_dat
  logical :: ok_dat
  real(dp) :: ekin
  real(dp) :: t_in, t_out, time_tot, time_config_in, time_config_out

  time_tot = 0.0_dp
  CALL CPU_TIME(time_config_in)

  ! ---------------------------
  ! CLI parsing
  !   ./phi4_mmc [-i input.inp] n_samp L n_steps n_jump n_realiz n_rest
  ! ---------------clear_output_files------------
  call parse_cli_phi4(input_file, n_samp_cli, L_cli, n_steps_cli, n_jump_cli, n_realiz_cli, n_rest_cli)

  ! ---------------------------
  ! Read input file (namelist)
  ! ---------------------------
  call read_input_file(trim(input_file), par, io, rng)

  ! ---------------------------
  ! Override from CLI (legacy-compatible)
  ! ---------------------------
  par%L       = L_cli
  par%N       = par%L

  ! ---------------------------
  ! Auto-normalize mean-field coupling using runtime N
  ! Target: -(1/(4N)) (sum_i q_i)^2
  ! In this code the MF piece is:  -(2*coup/(N-1)) (sum_i q_i)^2
  ! => choose coup = (N-1)/(8N)
  ! ---------------------------
  if (par%auto_coup) then
    if (par%N <= 1) stop "ERROR: auto_coup requires N>1"
    par%coup = real(par%N - 1, dp) / (8.0_dp * real(par%N, dp))
  end if

  par%n_steps = n_steps_cli
  par%n_jump  = n_jump_cli
  n_print     = par%n_print
  ! External run identifier (trajectory id)
  par%n_realiz = n_realiz_cli
  
  ! Segment index for chained restarts
  par%n_rest   = n_rest_cli

  ! --- legacy sampling id + energy-per-site ---
  par%n_samp = n_samp_cli


  par%en_in = 0.0001_dp * real(par%n_samp, dp)
  par%Etot  = real(par%N, dp) * par%en_in

  ! ---------------------------
  ! Remove old output
  ! ---------------------------
  if (par%n_rest == 1) then
    call clear_output_files1(io, par, par%n_realiz)
  end if

  ! ---------------------------
  ! RNG + allocate
  ! ---------------------------
  call rng_init(rng, par%seed)

  call allocate_state(par, st)

  ! ---------------------------
  ! Realizations loop
  ! ---------------------------

  ! initialize the observables
    call reset_obs(obs)
    n_MC = 0
    step0 = 0

    ! Set initial conditions or restart the trajectory
    if (par%n_rest == 1) then
      call init_configuration(par, st, rng)
      n_MC = 0
      st_saved = st

      call init_bins_autorange(par, st, bins, &
      n_warm = 20000, frac = 1.0_dp, &
      nbins_chi = 500, chi_log = .true., &
      nbins_rj = 300, nbins_rreg = 100, nbins_gradV2 = 300, nbins_lapV = 300)
      
      st = st_saved
      call reset_bins(bins)
      call reset_obs(obs)
      n_MC = 0
    else
      call read_restart(io, par, st, obs, bins, rng, step0, par%n_realiz)   ! auto: reads only if n_rest>1
      call reinflate(obs)  !
      n_MC = obs%n_MC
    end if

      ! Reset dynamical diagnostics contaminated by the warmup used
      ! only to determine histogram ranges.
    st%lya_logsum = 0.0_dp
    st%jlc_logsum = 0.0_dp
    st%lya_time   = 0.0_dp
    st%jlc_time   = 0.0_dp
    st%lya_count  = 0
    st%jlc_count  = 0
    st%s_jac      = 0.0_dp

    if (par%n_rest == 1) then
      call clear_output_files2(io, par, par%n_realiz)
    end if


    ! Check if the total length has been already reached
if (par%n_rest > 1) then
  if (step0 >= par%n_steps) then

    call read_last_nmc_from_dat(io, par, par%n_realiz, par%n_rest-1, nmc_dat, ok_dat)

    if (ok_dat .and. nmc_dat < par%n_steps) then
!      write(*,*) "WARNING: restart step0 seems complete (", step0, ")."
!      write(*,*) "         Overriding from .dat last n_MC = ", nmc_dat

      step0    = nmc_dat
      n_MC     = nmc_dat
      obs%n_MC = nmc_dat
    else
      call cli_msg_restart_already_complete(par, step0)
      stop 0
    end if

  end if
end if

    call cpu_time(time_config_out)
    time_tot = time_config_out - time_config_in



    ! ---- INTEGRATION ----
    do step = step0 + 1, par%n_steps
      call cpu_time(t_in)

      ! ---- MONTE CARLO STEP ----
      DO n_sweep = 1, par%n_jump
         call vv_step_md_tangent_jlc(par, st)
      END DO
      ! ---------------------
      call cpu_time(t_out)
      time_tot = time_tot + abs(t_out - t_in)

      ! --- early stop close to walltime ---
      if (time_tot >= 0.95_dp * par%cluster_time) then
        obs%n_MC = n_MC
        call print_observables(io, par, obs, st, time_tot, n_MC)
        exit
      end if
      n_MC = n_MC + 1
      ekin = st%K

      !    --- measurement schedule ---
      call update_obs(par, st, obs, n_MC)
      call update_bins(st, bins)

      if (par%time_dyn) then
        if ( mod(n_MC, par%n_print_inst) == 0) then
          call append_inst_micro(io, par, st, obs, n_MC)
        end if
      end if

      !    --- printing schedule ---
      if (mod(n_MC, n_print) == 0) then
        call print_observables(io, par, obs, st, time_tot, n_MC)
        call write_bins(io, par, bins)
      end if

    end do

    obs%n_MC = n_MC
    call write_restart(io, par, st, obs, bins, rng, par%n_steps, par%n_realiz)
    call write_bins(io, par, bins)
end program main_phi4
