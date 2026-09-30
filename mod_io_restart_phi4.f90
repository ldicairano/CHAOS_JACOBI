module mod_io_restart_phi4
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4State, Phi4Obs, Phi4Bins, IOParams, RNGState
  implicit none
  private

  public :: read_restart, write_restart
  public :: read_last_nmc_from_dat

  integer, parameter :: RESTART_VERSION = 20260304

contains

  subroutine ensure_dir(path)
    character(len=*), intent(in) :: path
    integer :: istat
    if (len_trim(path) == 0) return
    call execute_command_line("mkdir -p " // trim(path), exitstat=istat)
    if (istat /= 0) then
      write(*,*) "ERROR: cannot create directory: ", trim(path)
      stop
    end if
  end subroutine ensure_dir


  subroutine restart_fname(io, par, ir, nrest, fname)
    type(IOParams),   intent(in)  :: io
    type(Phi4Params), intent(in)  :: par
    integer,          intent(in)  :: ir, nrest
    character(len=*), intent(out) :: fname
    character(len=256) :: base

    write(base,'("restart_micro_",I0,"_",I0,"_",I0,"_",I0,".bin")') &
      par%n_samp, ir, par%N, nrest

    fname = trim(io%restart_dir)//"/"//trim(base)
  end subroutine restart_fname


  subroutine write_restart(io, par, st, obs, bins, rng, step, ir, nrest_out)
    integer, intent(in), optional :: nrest_out
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    type(Phi4State),  intent(in) :: st
    type(Phi4Obs),    intent(in) :: obs
    type(Phi4Bins),   intent(in) :: bins
    type(RNGState),   intent(in) :: rng
    integer,          intent(in) :: step
    integer,          intent(in) :: ir

    integer :: nrest_write, u
    character(len=512) :: fname

    nrest_write = par%n_rest
    if (present(nrest_out)) nrest_write = nrest_out

    call ensure_dir(io%restart_dir)
    call restart_fname(io, par, ir, nrest_write, fname)

    open(newunit=u, file=trim(fname), access="stream", form="unformatted", &
         status="replace", action="write")

    write(u) RESTART_VERSION

    ! --- minimal params ---
    write(u) par%L, par%N
    write(u) par%coup, par%mu, par%lambda
    write(u) par%Etot
    write(u) par%n_samp, par%en_in
    write(u) par%dt_md
    write(u) nrest_write

    ! --- run position ---
    write(u) step
    write(u) rng%seed

    ! --- observables ---
    write(u) obs%n_MC

    write(u) obs%av_mag,   obs%av_mag2,   obs%av_mag4
    write(u) obs%av_pot,   obs%av_kin
    write(u) obs%av_kin_1, obs%av_kin_2,  obs%av_kin_3
    write(u) obs%av_en,    obs%av_en2

    write(u) obs%av_jacW,    obs%av_jacW2
    write(u) obs%av_gradV2,  obs%av_gradV22
    write(u) obs%av_lapV,    obs%av_lapV2
    write(u) obs%av_jacR,    obs%av_jacR2
    write(u) obs%av_jacRreg, obs%av_jacRreg2

    write(u) obs%av_lya_tan, obs%av_lya_jlc

    write(u) obs%c_mag,   obs%c_mag2,   obs%c_mag4,   obs%binder
    write(u) obs%c_pot,   obs%c_kin
    write(u) obs%c_kin_1, obs%c_kin_2,  obs%c_kin_3
    write(u) obs%c_en,    obs%c_en2
    write(u) obs%beta,    obs%ders_2,   obs%ders_3

    write(u) obs%c_jacW,    obs%c_jacW2,    obs%var_jacW
    write(u) obs%c_gradV2,  obs%c_gradV22,  obs%var_gradV2
    write(u) obs%c_lapV,    obs%c_lapV2,    obs%var_lapV
    write(u) obs%c_jacR,    obs%c_jacR2,    obs%var_jacR
    write(u) obs%c_jacRreg, obs%c_jacRreg2, obs%var_jacRreg

    write(u) obs%c_lya_tan, obs%c_lya_jlc

    ! --- bins metadata ---
    write(u) bins%nbins_chi, bins%chi_min, bins%chi_max, bins%chi_log
    write(u) bins%nbins_rj,  bins%rj_min,  bins%rj_max
    write(u) bins%nbins_rreg, bins%rreg_min, bins%rreg_max

    ! --- bins diagnostics ---
    write(u) bins%n_bad_chi, bins%n_bad_rj, bins%n_bad_rreg

    ! --- bins arrays ---
    write(u) bins%count_chi
    write(u) bins%sum_Rj_vs_chi
    write(u) bins%sum_Rreg_vs_chi
    write(u) bins%sum_gradV2_vs_chi
    write(u) bins%sum_lapV_vs_chi

    write(u) bins%count_rj
    write(u) bins%count_rreg


    ! --- state scalars ---
    write(u) st%V, st%K
    write(u) st%S1, st%S2, st%S4
    write(u) st%M, st%phi2, st%phi4

    write(u) st%jacW, st%gradV2, st%lapV, st%jacR, st%jacRreg
    write(u) st%lya_logsum, st%lya_time
    write(u) st%jlc_logsum, st%jlc_time

    ! --- state arrays ---
    write(u) st%phi
    write(u) st%pi
    write(u) st%force
    write(u) st%gradV
    write(u) st%hdiag

    write(u) st%dphi
    write(u) st%dpi_t

    write(u) st%jphi
    write(u) st%jpi

    close(u)
  end subroutine write_restart


  subroutine read_restart(io, par, st, obs, bins, rng, step0, ir)
    type(IOParams),   intent(in)    :: io
    type(Phi4Params), intent(inout) :: par
    type(Phi4State),  intent(inout) :: st
    type(Phi4Obs),    intent(inout) :: obs
    type(Phi4Bins),   intent(inout) :: bins
    type(RNGState),   intent(inout) :: rng
    integer,          intent(out)   :: step0
    integer,          intent(in)    :: ir

    integer :: u, ios, ver
    integer :: Lr, Nr, nrestr
    real(dp) :: coupr, mur, lambdar
    real(dp) :: Etotr, dtmdr
    integer  :: nsampr
    real(dp) :: eninr
    character(len=512) :: fname

    step0 = 0

    if (par%n_rest <= 1) return

    call restart_fname(io, par, ir, par%n_rest-1, fname)

    open(newunit=u, file=trim(fname), access="stream", form="unformatted", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) then
      write(*,*) "ERROR: restart file not found for n_rest-1."
      write(*,*) "  expected: ", trim(fname)
      stop
    end if

    read(u) ver
    if (ver /= RESTART_VERSION) then
      write(*,*) "ERROR: restart version mismatch. file=", ver, " expected=", RESTART_VERSION
      stop
    end if

    read(u) Lr, Nr
    read(u) coupr, mur, lambdar
    read(u) Etotr
    read(u) nsampr, eninr
    read(u) dtmdr
    read(u) nrestr

    if (nrestr /= par%n_rest-1) then
      write(*,*) "ERROR: restart mismatch (n_rest). file:", nrestr, " expected:", par%n_rest-1
      stop
    end if
    if (Lr /= par%L .or. Nr /= par%N) then
      write(*,*) "ERROR: restart mismatch (L,N). file:", Lr, Nr, " input:", par%L, par%N
      stop
    end if
    if (abs(coupr-par%coup) > 1.0e-12_dp .or. abs(mur-par%mu) > 1.0e-12_dp .or. &
        abs(lambdar-par%lambda) > 1.0e-12_dp) then
      write(*,*) "ERROR: restart mismatch (model params)."
      stop
    end if

    par%Etot   = Etotr
    par%n_samp = nsampr
    par%en_in  = eninr
    par%dt_md  = dtmdr

    read(u) step0
    read(u) rng%seed

    ! --- observables ---
    read(u) obs%n_MC

    read(u) obs%av_mag,   obs%av_mag2,   obs%av_mag4
    read(u) obs%av_pot,   obs%av_kin
    read(u) obs%av_kin_1, obs%av_kin_2,  obs%av_kin_3
    read(u) obs%av_en,    obs%av_en2

    read(u) obs%av_jacW,    obs%av_jacW2
    read(u) obs%av_gradV2,  obs%av_gradV22
    read(u) obs%av_lapV,    obs%av_lapV2
    read(u) obs%av_jacR,    obs%av_jacR2
    read(u) obs%av_jacRreg, obs%av_jacRreg2

    read(u) obs%av_lya_tan, obs%av_lya_jlc

    read(u) obs%c_mag,   obs%c_mag2,   obs%c_mag4,   obs%binder
    read(u) obs%c_pot,   obs%c_kin
    read(u) obs%c_kin_1, obs%c_kin_2,  obs%c_kin_3
    read(u) obs%c_en,    obs%c_en2
    read(u) obs%beta,    obs%ders_2,   obs%ders_3

    read(u) obs%c_jacW,    obs%c_jacW2,    obs%var_jacW
    read(u) obs%c_gradV2,  obs%c_gradV22,  obs%var_gradV2
    read(u) obs%c_lapV,    obs%c_lapV2,    obs%var_lapV
    read(u) obs%c_jacR,    obs%c_jacR2,    obs%var_jacR
    read(u) obs%c_jacRreg, obs%c_jacRreg2, obs%var_jacRreg

    read(u) obs%c_lya_tan, obs%c_lya_jlc


    if (allocated(bins%count_chi)) then
      if (size(bins%count_chi) /= bins%nbins_chi) then
        deallocate(bins%count_chi, bins%sum_Rj_vs_chi, bins%sum_Rreg_vs_chi, &
                   bins%sum_gradV2_vs_chi, bins%sum_lapV_vs_chi)
      end if
    end if
    if (.not. allocated(bins%count_chi)) then
      allocate(bins%count_chi(bins%nbins_chi))
      allocate(bins%sum_Rj_vs_chi(bins%nbins_chi))
      allocate(bins%sum_Rreg_vs_chi(bins%nbins_chi))
      allocate(bins%sum_gradV2_vs_chi(bins%nbins_chi))
      allocate(bins%sum_lapV_vs_chi(bins%nbins_chi))
    end if

    if (allocated(bins%count_rj)) then
      if (size(bins%count_rj) /= bins%nbins_rj) deallocate(bins%count_rj)
    end if
    if (.not. allocated(bins%count_rj)) allocate(bins%count_rj(bins%nbins_rj))

    if (allocated(bins%count_rreg)) then
      if (size(bins%count_rreg) /= bins%nbins_rreg) deallocate(bins%count_rreg)
    end if
    if (.not. allocated(bins%count_rreg)) allocate(bins%count_rreg(bins%nbins_rreg))

    ! --- bins metadata ---
    read(u) bins%nbins_chi, bins%chi_min, bins%chi_max, bins%chi_log
    read(u) bins%nbins_rj,  bins%rj_min,  bins%rj_max
    read(u) bins%nbins_rreg, bins%rreg_min, bins%rreg_max

    ! --- bins diagnostics ---
    read(u) bins%n_bad_chi, bins%n_bad_rj, bins%n_bad_rreg



    ! --- state scalars ---
    read(u) st%V, st%K
    read(u) st%S1, st%S2, st%S4
    read(u) st%M, st%phi2, st%phi4

    read(u) st%jacW, st%gradV2, st%lapV, st%jacR, st%jacRreg
    read(u) st%lya_logsum, st%lya_time
    read(u) st%jlc_logsum, st%jlc_time

    ! --- state arrays ---
    read(u) st%phi
    read(u) st%pi
    read(u) st%force
    read(u) st%gradV
    read(u) st%hdiag

    read(u) st%dphi
    read(u) st%dpi_t

    read(u) st%jphi
    read(u) st%jpi

    close(u)

    ! restore derived means if needed
    st%M    = st%S1 / real(par%N, dp)
    st%phi2 = st%S2 / real(par%N, dp)
    st%phi4 = st%S4 / real(par%N, dp)
  end subroutine read_restart


  subroutine read_last_nmc_from_dat(io, par, ir, nrest, nmc_last, found)
    use mod_kinds,      only: dp
    use mod_types_phi4, only: IOParams, Phi4Params
    implicit none
    type(IOParams),   intent(in)  :: io
    type(Phi4Params), intent(in)  :: par
    integer,          intent(in)  :: ir, nrest
    integer,          intent(out) :: nmc_last
    logical,          intent(out) :: found

    character(len=512) :: fname, line
    integer :: u, ios, ios2
    real(dp) :: t_dummy
    integer  :: nmc_tmp

    found    = .false.
    nmc_last = 0

    write(fname,'(A,"/obs_micro_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, nrest, par%N

    open(newunit=u, file=trim(fname), status="old", action="read", iostat=ios)
    if (ios /= 0) return

    do
      read(u,'(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      read(line, *, iostat=ios2) t_dummy, nmc_tmp
      if (ios2 == 0) then
        nmc_last = nmc_tmp
        found = .true.
      end if
    end do

    close(u)
  end subroutine read_last_nmc_from_dat

end module mod_io_restart_phi4