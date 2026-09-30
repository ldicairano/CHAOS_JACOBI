module mod_input_phi4
  use mod_kinds,              only: dp
  use mod_types_phi4,         only: Phi4Params, IOParams, RNGState
  implicit none
  private
  public :: read_input_file

contains

  subroutine read_input_file(fname, par, io, rng)
    character(len=*), intent(in) :: fname
    type(Phi4Params), intent(inout) :: par
    type(IOParams),   intent(inout) :: io
    type(RNGState),   intent(inout) :: rng

    integer :: u, ios

    ! --------- mirrors (defaults from structs) ----------
    real(dp) :: cluster_time

    real(dp) :: coup, mu, lambda
    logical  :: auto_coup, time_dyn

    integer :: n_steps, n_jump, n_therm, n_realiz, n_print

    character(len=256) :: out_dir, restart_dir

    integer :: seed, n_print_inst

    ! --------- namelists (match input.inp groups) ----------
    namelist /system/ cluster_time
    namelist /model/  coup, mu, lambda, auto_coup
    namelist /run/    n_steps, n_jump, n_therm, n_realiz, n_print, time_dyn, n_print_inst
    namelist /ioctl/  out_dir, restart_dir
    namelist /rngctl/ seed

    ! Defaults
    cluster_time = par%cluster_time

    coup = par%coup; mu = par%mu; lambda = par%lambda
    auto_coup = par%auto_coup

    n_steps  = par%n_steps
    n_jump   = par%n_jump
    n_therm  = par%n_therm
    n_realiz = par%n_realiz
    n_print  = par%n_print

    n_print_inst = par%n_print_inst
    time_dyn     = par%time_dyn

    out_dir     = io%out_dir
    restart_dir = io%restart_dir

    seed = par%seed

    open(newunit=u, file=trim(fname), status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR: cannot open input file: ", trim(fname)
      stop
    end if

    ! Read each group independently (order in file does not matter)
    rewind(u); read(u, nml=system, iostat=ios)
    rewind(u); read(u, nml=model,  iostat=ios)
    rewind(u); read(u, nml=run,    iostat=ios)
    rewind(u); read(u, nml=ioctl,  iostat=ios)
    rewind(u); read(u, nml=rngctl, iostat=ios)

    close(u)

    ! Normalize ensemble

    par%cluster_time = cluster_time

    par%coup = coup
    par%mu = mu
    par%lambda = lambda
    par%auto_coup = auto_coup

    par%n_steps  = n_steps
    par%n_jump   = n_jump
    par%n_therm  = n_therm
    par%n_realiz = n_realiz
    par%n_print  = n_print

    par%n_print_inst = n_print_inst
    par%time_dyn     = time_dyn

    par%seed = seed
    rng%seed = seed

    io%out_dir     = out_dir
    io%restart_dir = restart_dir

  end subroutine read_input_file



  
end module mod_input_phi4