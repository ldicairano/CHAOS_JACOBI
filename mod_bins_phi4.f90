module mod_bins_phi4
  use mod_kinds,      only: dp
  use mod_types_phi4, only: Phi4State, Phi4Bins, IOParams, Phi4Params
  implicit none
  private

  public :: init_bins
  public :: reset_bins
  public :: update_bins
  public :: write_bins

contains

  subroutine init_bins(bins, &
                       nbins_chi, chi_min, chi_max, chi_log, &
                       nbins_rj, rj_min, rj_max, &
                       nbins_rreg, rreg_min, rreg_max, &
                       nbins_gradV2, gradV2_min, gradV2_max, &
                       nbins_lapV, lapV_min, lapV_max)
    type(Phi4Bins), intent(inout) :: bins

    integer,  intent(in) :: nbins_chi, nbins_rj, nbins_rreg, nbins_gradV2, nbins_lapV
    real(dp), intent(in) :: chi_min, chi_max, rj_min, rj_max, rreg_min, rreg_max
    real(dp), intent(in) :: gradV2_min, gradV2_max, lapV_min, lapV_max
    logical,  intent(in) :: chi_log

    bins%nbins_chi = nbins_chi
    bins%chi_min   = chi_min
    bins%chi_max   = chi_max
    bins%chi_log   = chi_log

    bins%nbins_rj = nbins_rj
    bins%rj_min   = rj_min
    bins%rj_max   = rj_max

    bins%nbins_rreg = nbins_rreg
    bins%rreg_min   = rreg_min
    bins%rreg_max   = rreg_max

    bins%nbins_gradV2 = nbins_gradV2
    bins%gradV2_min   = gradV2_min
    bins%gradV2_max   = gradV2_max

    bins%nbins_lapV = nbins_lapV
    bins%lapV_min   = lapV_min
    bins%lapV_max   = lapV_max

    allocate(bins%count_chi(bins%nbins_chi))
    allocate(bins%sum_Rj_vs_chi(bins%nbins_chi))
    allocate(bins%sum_Rj2_vs_chi(bins%nbins_chi))
    allocate(bins%sum_Rreg_vs_chi(bins%nbins_chi))
    allocate(bins%sum_Rreg2_vs_chi(bins%nbins_chi))
    allocate(bins%sum_gradV2_vs_chi(bins%nbins_chi))
    allocate(bins%sum_gradV22_vs_chi(bins%nbins_chi))
    allocate(bins%sum_lapV_vs_chi(bins%nbins_chi))
    allocate(bins%sum_lapV2_vs_chi(bins%nbins_chi))

    allocate(bins%count_rj(bins%nbins_rj))
    allocate(bins%count_rreg(bins%nbins_rreg))
    allocate(bins%count_gradV2(bins%nbins_gradV2))
    allocate(bins%count_lapV(bins%nbins_lapV))

    call reset_bins(bins)
  end subroutine init_bins


  subroutine reset_bins(bins)
    type(Phi4Bins), intent(inout) :: bins

    bins%nsamples_chi    = 0
    bins%nsamples_rj     = 0
    bins%nsamples_rreg   = 0
    bins%nsamples_gradV2 = 0
    bins%nsamples_lapV   = 0

    bins%n_bad_chi    = 0
    bins%n_bad_rj     = 0
    bins%n_bad_rreg   = 0
    bins%n_bad_gradV2 = 0
    bins%n_bad_lapV   = 0

    bins%count_chi = 0
    bins%sum_Rj_vs_chi     = 0.0_dp
    bins%sum_Rj2_vs_chi    = 0.0_dp
    bins%sum_Rreg_vs_chi   = 0.0_dp
    bins%sum_Rreg2_vs_chi  = 0.0_dp
    bins%sum_gradV2_vs_chi = 0.0_dp
    bins%sum_gradV22_vs_chi= 0.0_dp
    bins%sum_lapV_vs_chi   = 0.0_dp
    bins%sum_lapV2_vs_chi  = 0.0_dp

    bins%count_rj     = 0
    bins%count_rreg   = 0
    bins%count_gradV2 = 0
    bins%count_lapV   = 0

    bins%sum_k2_global  = 0.0_dp
    bins%sum_k22_global = 0.0_dp

    bins%sum_k2_global  = 0.0_dp
    bins%sum_k22_global = 0.0_dp
  end subroutine reset_bins


  integer function bin_index_linear(x, xmin, xmax, nbins) result(ibin)
    real(dp), intent(in) :: x, xmin, xmax
    integer,  intent(in) :: nbins
    real(dp) :: t

    if (x < xmin .or. x >= xmax) then
      ibin = 0
      return
    end if

    t = (x - xmin) / (xmax - xmin)
    ibin = 1 + int(t * real(nbins, dp))

    if (ibin < 1) ibin = 1
    if (ibin > nbins) ibin = nbins
  end function bin_index_linear


  integer function bin_index_log(x, xmin, xmax, nbins) result(ibin)
    real(dp), intent(in) :: x, xmin, xmax
    integer,  intent(in) :: nbins
    real(dp) :: lx, lmin, lmax, t

    if (x <= 0.0_dp) then
      ibin = 0
      return
    end if
    if (x < xmin .or. x >= xmax) then
      ibin = 0
      return
    end if

    lx   = log10(x)
    lmin = log10(xmin)
    lmax = log10(xmax)

    t = (lx - lmin) / (lmax - lmin)
    ibin = 1 + int(t * real(nbins, dp))

    if (ibin < 1) ibin = 1
    if (ibin > nbins) ibin = nbins
  end function bin_index_log


  subroutine update_bins(st, bins)
    type(Phi4State), intent(in)    :: st
    type(Phi4Bins),  intent(inout) :: bins

    integer :: ibin
    real(dp) :: chi

    ! =====================================================
    ! Histogram of chi and conditional averages vs chi
    ! =====================================================
    chi = st%jacW

    if (chi > 0.0_dp) then
      if (bins%chi_log) then
        ibin = bin_index_log(chi, bins%chi_min, bins%chi_max, bins%nbins_chi)
      else
        ibin = bin_index_linear(chi, bins%chi_min, bins%chi_max, bins%nbins_chi)
      end if

      if (ibin > 0) then
        bins%nsamples_chi    = bins%nsamples_chi + 1
        bins%count_chi(ibin) = bins%count_chi(ibin) + 1

        bins%sum_Rj_vs_chi(ibin)      = bins%sum_Rj_vs_chi(ibin)      + st%jacR
        bins%sum_Rj2_vs_chi(ibin)     = bins%sum_Rj2_vs_chi(ibin)     + st%jacR*st%jacR

        bins%sum_Rreg_vs_chi(ibin)    = bins%sum_Rreg_vs_chi(ibin)    + st%k2_jv
        bins%sum_Rreg2_vs_chi(ibin)   = bins%sum_Rreg2_vs_chi(ibin)   + st%k2_jv*st%k2_jv

        bins%sum_k2_global  = bins%sum_k2_global  + st%k2_jv
        bins%sum_k22_global = bins%sum_k22_global + st%k2_jv*st%k2_jv

        bins%sum_gradV2_vs_chi(ibin)  = bins%sum_gradV2_vs_chi(ibin)  + st%gradV2
        bins%sum_gradV22_vs_chi(ibin) = bins%sum_gradV22_vs_chi(ibin) + st%gradV2*st%gradV2

        bins%sum_lapV_vs_chi(ibin)    = bins%sum_lapV_vs_chi(ibin)    + st%lapV
        bins%sum_lapV2_vs_chi(ibin)   = bins%sum_lapV2_vs_chi(ibin)   + st%lapV*st%lapV
      else
        bins%n_bad_chi = bins%n_bad_chi + 1
      end if
    else
      bins%n_bad_chi = bins%n_bad_chi + 1
    end if

    ! =====================================================
    ! Histogram of R_J
    ! =====================================================
    if (st%jacR > -huge(1.0_dp) .and. st%jacR < huge(1.0_dp)) then
      ibin = bin_index_linear(st%jacR, bins%rj_min, bins%rj_max, bins%nbins_rj)
      if (ibin > 0) then
        bins%nsamples_rj     = bins%nsamples_rj + 1
        bins%count_rj(ibin)  = bins%count_rj(ibin) + 1
      else
        bins%n_bad_rj = bins%n_bad_rj + 1
      end if
    else
      bins%n_bad_rj = bins%n_bad_rj + 1
    end if

    ! =====================================================
    ! Histogram of R_J^reg
    ! =====================================================
    if (st%k2_jv > -huge(1.0_dp) .and. st%k2_jv < huge(1.0_dp)) then
      ibin = bin_index_linear(st%k2_jv, bins%rreg_min, bins%rreg_max, bins%nbins_rreg)
      if (ibin > 0) then
        bins%nsamples_rreg      = bins%nsamples_rreg + 1
        bins%count_rreg(ibin)   = bins%count_rreg(ibin) + 1
      else
        bins%n_bad_rreg = bins%n_bad_rreg + 1
      end if
    else
      bins%n_bad_rreg = bins%n_bad_rreg + 1
    end if

    ! =====================================================
    ! Histogram of gradV2
    ! =====================================================
    if (st%gradV2 > -huge(1.0_dp) .and. st%gradV2 < huge(1.0_dp)) then
      ibin = bin_index_linear(st%gradV2, bins%gradV2_min, bins%gradV2_max, bins%nbins_gradV2)
      if (ibin > 0) then
        bins%nsamples_gradV2       = bins%nsamples_gradV2 + 1
        bins%count_gradV2(ibin)    = bins%count_gradV2(ibin) + 1
      else
        bins%n_bad_gradV2 = bins%n_bad_gradV2 + 1
      end if
    else
      bins%n_bad_gradV2 = bins%n_bad_gradV2 + 1
    end if

    ! =====================================================
    ! Histogram of lapV
    ! =====================================================
    if (st%lapV > -huge(1.0_dp) .and. st%lapV < huge(1.0_dp)) then
      ibin = bin_index_linear(st%lapV, bins%lapV_min, bins%lapV_max, bins%nbins_lapV)
      if (ibin > 0) then
        bins%nsamples_lapV      = bins%nsamples_lapV + 1
        bins%count_lapV(ibin)   = bins%count_lapV(ibin) + 1
      else
        bins%n_bad_lapV = bins%n_bad_lapV + 1
      end if
    else
      bins%n_bad_lapV = bins%n_bad_lapV + 1
    end if
  end subroutine update_bins


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


  subroutine write_bins(io, par, bins)
    type(IOParams),   intent(in) :: io
    type(Phi4Params), intent(in) :: par
    type(Phi4Bins),   intent(in) :: bins

    integer :: u, i
    character(len=256) :: fname
    real(dp) :: xL, xR, xC, prob
    real(dp) :: meanRj, meanRj2, varRj
    real(dp) :: meanK2, meanK22, varK2
    real(dp) :: meanGrad, meanGrad2, varGrad
    real(dp) :: meanLap, meanLap2, varLap

    call ensure_dir(io%out_dir)

    ! =====================================================
    ! Bins in chi: occupancy + conditional means/variances
    ! =====================================================
    write(fname,'(A,"/bins_chi_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N

    open(newunit=u, file=trim(fname), status="replace", action="write")

    do i = 1, bins%nbins_chi
      if (bins%chi_log) then
        xL = 10.0_dp**( log10(bins%chi_min) + (real(i-1,dp)/real(bins%nbins_chi,dp)) * &
             (log10(bins%chi_max)-log10(bins%chi_min)) )
        xR = 10.0_dp**( log10(bins%chi_min) + (real(i  ,dp)/real(bins%nbins_chi,dp)) * &
             (log10(bins%chi_max)-log10(bins%chi_min)) )
      else
        xL = bins%chi_min + (real(i-1,dp)/real(bins%nbins_chi,dp)) * (bins%chi_max - bins%chi_min)
        xR = bins%chi_min + (real(i  ,dp)/real(bins%nbins_chi,dp)) * (bins%chi_max - bins%chi_min)
      end if
      xC = 0.5_dp*(xL + xR)

      if (bins%nsamples_chi > 0) then
        prob = real(bins%count_chi(i),dp) / real(bins%nsamples_chi,dp)
      else
        prob = 0.0_dp
      end if

      if (bins%count_chi(i) > 0) then
        meanRj   = bins%sum_Rj_vs_chi(i)  / real(bins%count_chi(i),dp)
        meanRj2  = bins%sum_Rj2_vs_chi(i) / real(bins%count_chi(i),dp)
        varRj    = meanRj2 - meanRj*meanRj

        meanK2  = bins%sum_Rreg_vs_chi(i)  / real(bins%count_chi(i),dp)
        meanK22 = bins%sum_Rreg2_vs_chi(i) / real(bins%count_chi(i),dp)
        varK2   = meanK22 - meanK2*meanK2

        meanGrad  = bins%sum_gradV2_vs_chi(i)  / real(bins%count_chi(i),dp)
        meanGrad2 = bins%sum_gradV22_vs_chi(i) / real(bins%count_chi(i),dp)
        varGrad   = meanGrad2 - meanGrad*meanGrad

        meanLap  = bins%sum_lapV_vs_chi(i)  / real(bins%count_chi(i),dp)
        meanLap2 = bins%sum_lapV2_vs_chi(i) / real(bins%count_chi(i),dp)
        varLap   = meanLap2 - meanLap*meanLap
      else
        meanRj = 0.0_dp;    varRj = 0.0_dp
        meanK2 = 0.0_dp;    varK2 = 0.0_dp
        meanGrad = 0.0_dp;  varGrad = 0.0_dp
        meanLap = 0.0_dp;   varLap = 0.0_dp
      end if

      write(u,'(ES22.14,1X,ES22.14,1X,ES22.14,1X,I12,1X,ES22.14,1X,'// &
               'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14,1X,'// &
               'ES22.14,1X,ES22.14,1X,ES22.14,1X,ES22.14)') &
        xL, xR, xC, bins%count_chi(i), prob, &
        meanRj, varRj, &
        meanK2, varK2, &
        meanGrad, varGrad, &
        meanLap, varLap
    end do

    close(u)

    ! =====================================================
    ! Histogram of R_J
    ! =====================================================
    write(fname,'(A,"/hist_rj_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N

    open(newunit=u, file=trim(fname), status="replace", action="write")
    do i = 1, bins%nbins_rj
      xL = bins%rj_min + (real(i-1,dp)/real(bins%nbins_rj,dp)) * (bins%rj_max - bins%rj_min)
      xR = bins%rj_min + (real(i  ,dp)/real(bins%nbins_rj,dp)) * (bins%rj_max - bins%rj_min)
      xC = 0.5_dp*(xL + xR)

      if (bins%nsamples_rj > 0) then
        prob = real(bins%count_rj(i),dp) / real(bins%nsamples_rj,dp)
      else
        prob = 0.0_dp
      end if

      write(u,'(ES22.14,1X,ES22.14,1X,ES22.14,1X,I12,1X,ES22.14)') xL, xR, xC, bins%count_rj(i), prob
    end do
    close(u)

    ! =====================================================
    ! Histogram of R_J^reg
    ! =====================================================
    write(fname,'(A,"/hist_k2_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N

    open(newunit=u, file=trim(fname), status="replace", action="write")
    do i = 1, bins%nbins_rreg
      xL = bins%rreg_min + (real(i-1,dp)/real(bins%nbins_rreg,dp)) * (bins%rreg_max - bins%rreg_min)
      xR = bins%rreg_min + (real(i  ,dp)/real(bins%nbins_rreg,dp)) * (bins%rreg_max - bins%rreg_min)
      xC = 0.5_dp*(xL + xR)

      if (bins%nsamples_rreg > 0) then
        prob = real(bins%count_rreg(i),dp) / real(bins%nsamples_rreg,dp)
      else
        prob = 0.0_dp
      end if

      write(u,'(ES22.14,1X,ES22.14,1X,ES22.14,1X,I12,1X,ES22.14)') xL, xR, xC, bins%count_rreg(i), prob
    end do
    close(u)

    ! =====================================================
    ! Histogram of gradV2
    ! =====================================================
    write(fname,'(A,"/hist_gradV2_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N

    open(newunit=u, file=trim(fname), status="replace", action="write")
    do i = 1, bins%nbins_gradV2
      xL = bins%gradV2_min + (real(i-1,dp)/real(bins%nbins_gradV2,dp)) * (bins%gradV2_max - bins%gradV2_min)
      xR = bins%gradV2_min + (real(i  ,dp)/real(bins%nbins_gradV2,dp)) * (bins%gradV2_max - bins%gradV2_min)
      xC = 0.5_dp*(xL + xR)

      if (bins%nsamples_gradV2 > 0) then
        prob = real(bins%count_gradV2(i),dp) / real(bins%nsamples_gradV2,dp)
      else
        prob = 0.0_dp
      end if

      write(u,'(ES22.14,1X,ES22.14,1X,ES22.14,1X,I12,1X,ES22.14)') xL, xR, xC, bins%count_gradV2(i), prob
    end do
    close(u)

    ! =====================================================
    ! Histogram of lapV
    ! =====================================================
    write(fname,'(A,"/hist_lapV_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N

    open(newunit=u, file=trim(fname), status="replace", action="write")
    do i = 1, bins%nbins_lapV
      xL = bins%lapV_min + (real(i-1,dp)/real(bins%nbins_lapV,dp)) * (bins%lapV_max - bins%lapV_min)
      xR = bins%lapV_min + (real(i  ,dp)/real(bins%nbins_lapV,dp)) * (bins%lapV_max - bins%lapV_min)
      xC = 0.5_dp*(xL + xR)

      if (bins%nsamples_lapV > 0) then
        prob = real(bins%count_lapV(i),dp) / real(bins%nsamples_lapV,dp)
      else
        prob = 0.0_dp
      end if

      write(u,'(ES22.14,1X,ES22.14,1X,ES22.14,1X,I12,1X,ES22.14)') xL, xR, xC, bins%count_lapV(i), prob
    end do
    close(u)

    ! =====================================================
    ! Diagnostics
    ! =====================================================
    write(fname,'(A,"/bins_diag_",I0,"_",I0,"_",I0,"_",I0,".dat")') &
      trim(io%out_dir), par%n_samp, par%n_realiz, par%n_rest, par%N

    open(newunit=u, file=trim(fname), status="replace", action="write")
    write(u,'("nsamples_chi    ",I12)') bins%nsamples_chi
    write(u,'("nsamples_rj     ",I12)') bins%nsamples_rj
    write(u,'("nsamples_rreg   ",I12)') bins%nsamples_rreg
    write(u,'("nsamples_gradV2 ",I12)') bins%nsamples_gradV2
    write(u,'("nsamples_lapV   ",I12)') bins%nsamples_lapV
    write(u,'("n_bad_chi       ",I12)') bins%n_bad_chi
    write(u,'("n_bad_rj        ",I12)') bins%n_bad_rj
    write(u,'("n_bad_rreg      ",I12)') bins%n_bad_rreg
    write(u,'("n_bad_gradV2    ",I12)') bins%n_bad_gradV2
    write(u,'("n_bad_lapV      ",I12)') bins%n_bad_lapV

    if (bins%nsamples_rreg > 0) then
      write(u,'("mean_k2        ",ES22.14)') bins%sum_k2_global / real(bins%nsamples_rreg,dp)
      write(u,'("var_k2         ",ES22.14)') bins%sum_k22_global / real(bins%nsamples_rreg,dp) - &
                                             (bins%sum_k2_global / real(bins%nsamples_rreg,dp))**2
    else
      write(u,'("mean_k2        ",ES22.14)') 0.0_dp
      write(u,'("var_k2         ",ES22.14)') 0.0_dp
    end if

    close(u)
  end subroutine write_bins


  ! subroutine update_ranges(st, rmin, rmax, gmin, gmax, lmin, lmax, rrmin, rrmax)
  !   use mod_kinds, only: dp
  !   use mod_types_phi4, only: Phi4State
  !   implicit none
  !   type(Phi4State), intent(in) :: st
  !   real(dp), intent(inout) :: rmin, rmax, gmin, gmax, lmin, lmax, rrmin, rrmax
  
  !   ! R_J
  !   if (st%jacR > -huge(1.0_dp) .and. st%jacR < huge(1.0_dp)) then
  !     rmin = min(rmin, st%jacR)
  !     rmax = max(rmax, st%jacR)
  !   end if
  
  !   ! Replace Rreg with your chosen observable (e.g. k2_jv) if you did that swap
  !   if (st%k2_jv > -huge(1.0_dp) .and. st%k2_jv < huge(1.0_dp)) then
  !     rrmin = min(rrmin, st%k2_jv)
  !     rrmax = max(rrmax, st%k2_jv)
  !   end if
  
  !   ! gradV2
  !   if (st%gradV2 > -huge(1.0_dp) .and. st%gradV2 < huge(1.0_dp)) then
  !     gmin = min(gmin, st%gradV2)
  !     gmax = max(gmax, st%gradV2)
  !   end if
  
  !   ! lapV
  !   if (st%lapV > -huge(1.0_dp) .and. st%lapV < huge(1.0_dp)) then
  !     lmin = min(lmin, st%lapV)
  !     lmax = max(lmax, st%lapV)
  !   end if
  ! end subroutine update_ranges


end module mod_bins_phi4