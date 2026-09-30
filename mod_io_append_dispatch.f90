module mod_io_append_dispatch
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4Params, Phi4Obs, IOParams, Phi4State
  implicit none
  private
  public :: print_observables, clear_output_files1, clear_output_files2
  public :: append_inst_micro


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


  subroutine print_observables(io, par, obs, st, time_tot, n_MC)
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    type(Phi4State),  intent(in) :: st
    type(Phi4Obs),    intent(in) :: obs
    real(dp),         intent(in) :: time_tot
    integer,          intent(in) :: n_MC

    character(len=256) :: fname
    integer :: u

    ! Columns:
    !  1  time_tot
    !  2  n_MC
    !  3  en_in
    !  4  dt_md
    !  5  c_en
    !  6  c_en2
    !  7  c_kin
    !  8  c_pot
    !  9  c_kin_1
    ! 10  c_kin_2
    ! 11  c_kin_3
    ! 12  beta
    ! 13  ders_2
    ! 14  ders_3
    ! 15  c_mag
    ! 16  c_mag2
    ! 17  c_mag4
    ! 18  binder
    ! 19  c_lya_tan
    ! 20  c_lya_jlc
    ! 21  c_jacW
    ! 22  var_jacW
    ! 23  c_gradV2
    ! 24  var_gradV2
    ! 25  c_lapV
    ! 26  var_lapV
    ! 27  c_jacR
    ! 28  var_jacR
    ! 29  c_jacRreg
    ! 30  var_jacRreg
    ! 31  jacR
    ! 32  jacRreg
    ! 33  s_jac
    ! 34  jacW
    ! 35  d_bdry_jac
    ! 36  k2_jv

    call ensure_dir(io%out_dir)

    write(fname,'(A,"/obs_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N

    open(newunit=u, file=trim(fname), status="unknown", action="write", position="append")

    write(u,'(ES22.14,1X,I12    ,1X,ES22.14,1X,ES22.14,1X,' // &
             'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14,1X,' // &
             'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14,1X,' // &
             'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14,1X,' // &
             'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14,1X,' // &
             'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14,1X,' // &
             'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14,1X,' // &
             'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14,1X,' // &
             'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14)') &
      time_tot, n_MC, par%en_in, par%dt_md, &
      obs%c_en, obs%c_en2, obs%c_kin, obs%c_pot, &
      obs%c_kin_1, obs%c_kin_2, obs%c_kin_3, &
      obs%beta, obs%ders_2, obs%ders_3, &
      obs%c_mag, obs%c_mag2, obs%c_mag4, obs%binder, &
      obs%c_lya_tan, obs%c_lya_jlc, &
      obs%c_jacW, obs%var_jacW, &
      obs%c_gradV2, obs%var_gradV2, &
      obs%c_lapV, obs%var_lapV, &
      obs%c_jacR, obs%var_jacR, &
      obs%c_jacRreg, obs%var_jacRreg, &
      st%jacR, st%jacRreg, st%s_jac, &
      st%jacW, st%d_bdry_jac, st%k2_jv

    close(u)
  end subroutine print_observables


  subroutine clear_output_files1(io, par, ir)
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    integer,          intent(in) :: ir

    integer :: u
    character(len=512) :: fname

    call ensure_dir(io%out_dir)

    ! --------------------------------------------------
    ! Main observables file
    ! --------------------------------------------------
    write(fname,'(A,"/obs_micro_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, par%n_rest, par%N
    open(newunit=u, file=trim(fname), status="replace", action="write")
    close(u)

  end subroutine clear_output_files1



  subroutine clear_output_files2(io, par, ir)
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    integer,          intent(in) :: ir

    integer :: u
    character(len=512) :: fname

    call ensure_dir(io%out_dir)

    ! --------------------------------------------------
    ! Main observables file
    ! --------------------------------------------------
    write(fname,'(A,"/obs_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, par%n_rest, par%N
    open(newunit=u, file=trim(fname), status="replace", action="write")
    close(u)

    ! --------------------------------------------------
    ! If you also want to clear the bins / hist files
    ! uncomment or keep these if they are produced in the same run
    ! --------------------------------------------------

    write(fname,'(A,"/bins_chi_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, par%n_rest, par%N
    open(newunit=u, file=trim(fname), status="replace", action="write")
    close(u)

    write(fname,'(A,"/hist_rj_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, par%n_rest, par%N
    open(newunit=u, file=trim(fname), status="replace", action="write")
    close(u)

    write(fname,'(A,"/hist_k2_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, par%n_rest, par%N
    open(newunit=u, file=trim(fname), status="replace", action="write")
    close(u)

    write(fname,'(A,"/hist_gradV2_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, par%n_rest, par%N
    open(newunit=u, file=trim(fname), status="replace", action="write")
    close(u)

    write(fname,'(A,"/hist_lapV_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, par%n_rest, par%N
    open(newunit=u, file=trim(fname), status="replace", action="write")
    close(u)

    write(fname,'(A,"/bins_diag_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, ir, par%n_rest, par%N
    open(newunit=u, file=trim(fname), status="replace", action="write")
    close(u)

  end subroutine clear_output_files2


  subroutine append_inst_micro(io, par, st, obs, step)
    use mod_kinds,      only: dp
    use mod_types_phi4, only: IOParams, Phi4Params, Phi4State, Phi4Obs
    implicit none
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    type(Phi4State),  intent(in) :: st
    type(Phi4Obs),    intent(in) :: obs
    integer,          intent(in) :: step
  
    integer :: u
    character(len=512) :: fname
  
    call ensure_dir(io%out_dir)
  
    ! File name: instantaneous time series
    write(fname,'(A,"/obs_time_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N
  
    open(newunit=u, file=trim(fname), status="unknown", action="write", position="append")
  
    ! One line: step + instantaneous values (NO running averages here)
    write(u,*) par%dt_md*step, &
      obs%c_en, par%en_in, obs%beta, &
      merge(st%lya_logsum/st%lya_time, 0.0_dp, st%lya_time>0.0_dp), &
      merge(st%jlc_logsum/st%jlc_time, 0.0_dp, st%jlc_time>0.0_dp), &
      st%jacW, st%jacR, st%jacRreg, st%s_jac, &
      st%d_bdry_jac, st%k2_jv
  
    close(u)
  end subroutine append_inst_micro


end module mod_io_append_dispatch