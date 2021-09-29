module cls_detector
  use cls_c3_data, only: c3_data
  use cls_line_text, only: line_max
  use mod_mpi, only: get_mpi_task_id
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  type detector
     private
     
     integer :: n_stations
     character(line_max), allocatable :: station_names(:)

     !integer :: n_histo = 9000
     !double precision :: thred = 0.95d0

    type(C_PTR)                            :: plan_r2c, plan_c2r
    real(C_DOUBLE), allocatable            :: r_tmp(:), r2_tmp(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: c_tmp(:), c2_tmp(:)


   contains
     !procedure :: calc_stats => detector_calc_stats
     procedure :: calc_correlogram => detector_calc_correlogram
     procedure :: run_cross_corr => detector_run_cross_corr
  end type detector
  
  interface detector
     module procedure init_detector
  end interface detector
  
contains
  
  !-------------------------------------------------------------------------
  type(detector) function init_detector(station_names) result(self)
    character(line_max), intent(in) :: station_names(:)
    integer :: n_stations
    
    n_stations = size(station_names)
    self%n_stations = n_stations
    allocate(self%station_names(n_stations))
    self%station_names = station_names

    
    return 
  end function init_detector
  
  !-------------------------------------------------------------------------

  subroutine detector_calc_correlogram(self, env, t_win, t_step)
    class(detector), intent(inout) :: self
    type(c3_data), intent(in) ::  env(:)
    double precision, intent(in) :: t_win, t_step
    integer :: i_sta, j_sta, n_pair, n_sta, i, j
    integer :: id_start, id_end, id
    integer :: count

    n_sta = self%n_stations
    if (size(env) /= n_sta) then
       error stop "ERROR: invalid dimension of instance array"
    end if

    
    do i_sta = 1, n_sta
       if (env(i_sta)%get_n_cmps() /= 1) then
          print *, "n_cmps=", env(i_sta)%get_n_cmps()
          error stop "ERROR: multiple components is not supported yet" // &
               & " (detector_calc_correlogram)"
       end if
    end do

    if (t_win < 0.d0) then
       error stop "ERROR: t_win must be > 0 (detector_calc_correlogram)"
    end if
    
    n_pair = n_sta * (n_sta - 1) / 2
    
    call get_mpi_task_id(n_pair, id_start, id_end)
    
    do id = id_start, id_end
       ! Assign station pair
       count = 0
       i_sta = 0
       j_sta = 0
       pair_assignment: do i = 1, n_sta - 1
          do j = i + 1, n_sta
             count = count + 1
             if (count == id) then
                i_sta = i
                j_sta = j
                exit pair_assignment
             end if
          end do
       end do pair_assignment
       
       print *, id, trim(self%station_names(i_sta)), trim(self%station_names(j_sta))
       call self%run_cross_corr(env(i_sta), env(j_sta), t_win, t_step)
       
    end do
    

    
    return 
  end subroutine detector_calc_correlogram
    
  !-------------------------------------------------------------------------

  subroutine detector_run_cross_corr(self, env1, env2, t_win, t_step)
    class(detector), intent(inout) :: self
    type(c3_data), intent(in) :: env1, env2
    double precision, intent(in) :: t_win, t_step
    integer :: n_smp, n_smp2, n, n_step, i, n_win, i1, i2, n2
    double precision :: dt, dt2, l1, l2
    double precision, allocatable :: x1(:), x2(:), cc(:)
    complex(kind(0d0)), allocatable :: c1(:), c2(:)
    logical :: first_flag

    n_smp = env1%get_n_smp()
    n_smp2 = env2%get_n_smp()
    if (n_smp /= n_smp2) then
       error stop "ERROR: n_smp /= n_smp2 "
    end if
    
    dt = env1%get_dt()
    dt2 = env2%get_dt()
    if (dt /= dt2) then
       error stop "ERROR dt /= dt2"
    end if
    write(*,*)dt
    n = nint(t_win / dt)
    allocate(self%r_tmp(n))
    allocate(self%c_tmp(n))
    allocate(self%r2_tmp(n))
    allocate(self%c2_tmp(n))
    self%plan_r2c = fftw_plan_dft_r2c_1d(n, self%r_tmp, self%c_tmp, &
         & FFTW_ESTIMATE)
    self%plan_c2r = fftw_plan_dft_c2r_1d(n, self%c2_tmp, self%r2_tmp, &
         & FFTW_ESTIMATE)
    
    n_step = nint(t_step / dt) 
    n_win = int((n_smp - n) / n_step)
    n2 = n / 2
    if (n2 + n2 /= n) then
       error stop "n2 + n2 /= n"
    end if

    allocate(x1(n))
    allocate(x2(n))
    allocate(c1(n))
    allocate(c2(n))
    allocate(cc(n))

    ! Main
    i1 = 1
    i2 = n
    do i = 1, n_win
       write(*,*) i, "/", n_win
       
       x1(1:n) = env1%extract_data(i1, i2, 1)
       x1(1:n) = x1(1:n) - sum(x1(1:n)) / n
       l1 = sqrt(sum(x1**2))
       x1 = x1 / l1
       self%r_tmp = x1
       call fftw_execute_dft_r2c(self%plan_r2c, self%r_tmp, self%c_tmp)
       c1 = self%c_tmp
       
       x2(1:n) = env2%extract_data(i1, i2, 1)
       x2(1:n) = x2(1:n) - sum(x2(1:n)) / n
       l2 = sqrt(sum(x2**2))
       x2 = x2 / l2
       self%r_tmp = x2
       call fftw_execute_dft_r2c(self%plan_r2c, self%r_tmp, self%c_tmp)
       c2 = self%c_tmp
       
       self%c2_tmp = c1(1:n) * conjg(c2(1:n))
       call fftw_execute_dft_c2r(self%plan_c2r, self%c2_tmp, self%r2_tmp)
       cc(1:n/2) = self%r2_tmp(n/2+1:n) / n
       cc(n/2+1:n) = self%r2_tmp(1:n/2) / n
       
       block 
         integer :: j
         do j = 1, n
            write(111,*)i, (j - n/2 - 1) * dt, cc(j)
         end do
       end block

       i1 = i1 + n_step
       i2 = i2 + n_step

       
    end do
    deallocate(self%r_tmp)
    deallocate(self%c_tmp)

    return 
  end subroutine detector_run_cross_corr

  !-------------------------------------------------------------------------
  !
  !subroutine detector_calc_stats(self, i_sta, mean, stdv)
  !  class(detector), intent(inout) :: self
  !  integer, intent(in) :: i_sta
  !  double precision, allocatable, intent(out) :: mean(:), stdv(:)
  !  double precision, allocatable :: pcnt(:)
  !  integer :: n_cmps, n_smp, i_cmp, io, i, k, count
  !  integer, allocatable :: histo(:,:)
  !  character(line_max) :: out_file
  !
  !  if (i_sta < 0 .or. i_sta > self%n_stations) then
  !     error stop "ERROR: invalid station ID is given to detector_calc_stats"
  !  end if
  !  
  !  n_cmps = self%c3(i_sta)%get_n_cmps()
  !  n_smp  = self%c3(i_sta)%get_n_smp()
  !  allocate(histo(self%n_histo, n_cmps))
  !  allocate(mean(n_cmps))
  !  allocate(stdv(n_cmps))
  !  
  !  mean = self%c3(i_sta)%calc_mean()
  !
  !  ! histogram
  !  histo = self%c3(i_sta)%calc_histogram(self%n_histo)
  !  do i_cmp = 1, n_cmps
  !     print *, "mean :", i_cmp, mean(i_cmp)
  !     
  !  end do
  !
  !  out_file = trim(self%station_names(i_sta)) // "_1.histo" 
  !  open(newunit=io, file=out_file, status='unknown')
  !  do i = 1, self%n_histo
  !     write(io, *)i, histo(i, 1)
  !  end do
  !  close(io)
  !
  !
  !  ! Rough estimation of percentile
  !  allocate(pcnt(n_cmps))
  !  k = int((1.d0 - self%thred) * n_smp)
  !
  !  do i_cmp = 1, n_cmps
  !     count = 0
  !     do i = self%n_histo, 1, -1
  !        count = count + histo(i, i_cmp)
  !        if (count >= k) then
  !           pcnt(i_cmp) = i
  !           exit
  !        end if
  !     end do
  !     print *, "pcnt :",  trim(self%station_names(i_sta)), " ", &
  !          & i_cmp, pcnt(i_cmp)
  !  end do
  !
  !
  !  return 
  !end subroutine detector_calc_stats

  !-------------------------------------------------------------------------
  
end module cls_detector
