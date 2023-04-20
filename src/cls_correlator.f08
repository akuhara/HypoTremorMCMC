module cls_correlator
  use cls_c3_data, only: c3_data
  use cls_line_text, only: line_max
  use mod_mpi
  use mod_signal_process, only: apply_taper
  use mod_sort, only: quick_sort
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  type correlator
     private
     
     integer :: n_stations
     character(line_max), allocatable :: station_names(:)

     type(C_PTR)                            :: plan_r2c, plan_c2r
     real(C_DOUBLE), allocatable            :: r_tmp(:), r2_tmp(:)
     complex(C_DOUBLE_COMPLEX), allocatable :: c_tmp(:), c2_tmp(:)

     double precision :: t_win
     double precision :: t_step
     double precision :: dt
     
     integer :: n_smp      ! # of samples in whole envelope
     integer :: n_win      ! # of time windows
     integer :: n_step     ! # of samples in a window interval
     integer :: n          ! # of samples in a time window


     double precision :: alpha 
     integer :: n_detect = 0
     
     logical :: debug = .false.

   contains
     !procedure :: calc_stats => correlator_calc_stats
     procedure :: calc_correlogram => correlator_calc_correlogram
     procedure :: run_cross_corr   => correlator_run_cross_corr
     procedure :: measure_lag_time => correlator_measure_lag_time
     procedure :: eval_correlogram => correlator_eval_correlogram
     
  end type correlator
  
  interface correlator
     module procedure init_correlator
  end interface correlator
  
contains
  
  !-------------------------------------------------------------------------
  type(correlator) function init_correlator(station_names, t_win, t_step, &
       & dt, n_smp, alpha) result(self)
    character(line_max), intent(in) :: station_names(:)
    double precision, intent(in) :: t_win, t_step, dt, alpha
    integer, intent(in) :: n_smp
    integer :: n_stations
    
    n_stations     = size(station_names)
    self%n_stations = n_stations
    
    allocate(self%station_names(n_stations))
    self%station_names = station_names
    
    self%t_win  = t_win
    self%t_step = t_step
    self%dt     = dt
    self%n_smp  = n_smp
    self%alpha  = alpha

    if (t_win < 0.d0) then
       error stop "ERROR: t_win must be > 0 (init_correlator)"
    end if
    if (t_step < 0.d0) then
       error stop "ERROR: t_step must be > 0 (init_correlator)"
    end if

    self%n      = nint(self%t_win  / self%dt)
    self%n_step = nint(self%t_step / self%dt) 
    self%n_win  = int((self%n_smp - self%n)  / self%n_step)
    
    allocate(self%r_tmp(self%n))
    allocate(self%c_tmp(self%n))
    allocate(self%r2_tmp(self%n))
    allocate(self%c2_tmp(self%n))

    self%plan_r2c = fftw_plan_dft_r2c_1d(self%n, self%r_tmp, self%c_tmp, &
         & FFTW_ESTIMATE)
    self%plan_c2r = fftw_plan_dft_c2r_1d(self%n, self%c2_tmp, self%r2_tmp, &
         & FFTW_ESTIMATE)
    
    
    return 
  end function init_correlator

  !-------------------------------------------------------------------------

  subroutine correlator_calc_correlogram(self, env)
    class(correlator), intent(inout) :: self
    type(c3_data), intent(in) ::  env(:)
    integer :: i_sta, j_sta, n_pair, n_sta, i, j, io, ierr
    integer :: id_start, id_end, id, rank
    integer :: count

    n_sta = self%n_stations
    if (size(env) /= n_sta) then
       error stop "ERROR: invalid dimension of instance array"
    end if

    
    do i_sta = 1, n_sta
       if (env(i_sta)%get_n_cmps() /= 1) then
          print *, "n_cmps=", env(i_sta)%get_n_cmps()
          error stop "ERROR: multiple components is not supported yet" // &
               & " (correlator_calc_correlogram)"
       end if
    end do


    
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
       
       ! Check
       if (   self%n_smp /= env(i_sta)%get_n_smp()  .or. &
            & self%n_smp /= env(j_sta)%get_n_smp()) then
          error stop "ERROR: invalid n_smp"
       end if
       if (   self%dt /= env(i_sta)%get_dt()    .or. &
            & self%dt /= env(j_sta)%get_dt()) then
          error stop "ERROR: invalid dt"
       end if
       
       call self%run_cross_corr( &
            & env1   = env(i_sta), &
            & env2   = env(j_sta), &
            & sta1   = self%station_names(i_sta), &
            & sta2   = self%station_names(j_sta))
       
    end do
    
    
    return 
  end subroutine correlator_calc_correlogram
    
  !-------------------------------------------------------------------------

  subroutine correlator_run_cross_corr(self, env1, env2, sta1, sta2)
    class(correlator), intent(inout) :: self
    type(c3_data), intent(in) :: env1, env2
    character(line_max), intent(in) :: sta1, sta2
    integer :: n, i, i1, i2, n2, io, j, ierr, io2, i_corr_max
    character(line_max) :: out_file, out_file2
    double precision :: l1, l2, cc_max
    double precision, allocatable :: x1(:), x2(:), cc(:,:)
    complex(kind(0d0)), allocatable :: c1(:), c2(:)
    
    
    n = self%n
    n2 = n / 2
    if (n2 + n2 /= n) then
       error stop "n2 + n2 /= n"
    end if

    allocate(x1(n))
    allocate(x2(n))
    allocate(c1(n))
    allocate(c2(n))
    allocate(cc(n,self%n_win))
    
    ! Main
    i1 = 1
    i2 = n
    write(out_file, '(a)') trim(sta1) // "." // trim(sta2) // ".corr"
    write(out_file2, '(a)') trim(sta1) // "." // trim(sta2) // ".max_corr"
    print *, trim(out_file), "n_win=", self%n_win, "n=", n
    open(newunit=io, file=out_file, access="stream", form="unformatted", &
         & iostat=ierr, status="replace")
    open(newunit=io2, file=out_file2, access="stream", form="unformatted", &
         & iostat=ierr, status="replace")
    !open(newunit=io, file=out_file, form="formatted", &
    !     & iostat=ierr, status="replace")
    if (ierr /= 0) then
       error stop "ERROR: cannot create correlation output file"
    end if
    
    do i = 1, self%n_win
       
       x1(1:n) = env1%extract_data(i1, i2, 1)
       x1(1:n) = apply_taper(n, x1)
       x1 = x1 - sum(x1) / n
       l1 = sqrt(sum(x1**2))
       self%r_tmp = x1 / l1
       
       self%c_tmp = (0.d0, 0.d0) ! Not sure for reason, but seems necessary
       call fftw_execute_dft_r2c(self%plan_r2c, self%r_tmp(1:n), self%c_tmp(1:n))
       c1 = self%c_tmp
       
       x2(1:n) = env2%extract_data(i1, i2, 1)
       x2(1:n) = apply_taper(n, x2)
       x2 = x2 - sum(x2) / n
       l2 = sqrt(sum(x2**2))
       self%r_tmp = x2 / l2

       self%c_tmp = (0.d0, 0.d0) ! Not sure for reason, but seems necessary
       call fftw_execute_dft_r2c(self%plan_r2c, self%r_tmp(1:n), self%c_tmp(1:n))
       c2 = self%c_tmp
       
       
       self%c2_tmp = conjg(c1)*c2
       call fftw_execute_dft_c2r(self%plan_c2r, self%c2_tmp(1:n), self%r2_tmp(1:n))
       cc(1:n/2,i) = self%r2_tmp(n/2+1:n) / n
       cc(n/2+1:n,i) = self%r2_tmp(1:n/2) / n
       cc_max = maxval(cc(:,i))
       
       !if (cc_max(i) >= self%cc_thred) then
       !   self%n_detect = self%n_detect + 1
       !   i_corr_max = maxloc(cc,dim=1)
       !   i_corr_max = i_corr_max - (n/2 + 1)
       !   !Note : i_corr_max = 1 corresponds to zero-lag time
       !   call self%measure_lag_time(x1, x2, l1, l2, i_corr_max, sta1, sta2)
       !end if

       ! Output to file
       do j = 1, n
          write(io)(i-1)*self%t_step + 0.5d0 * self%t_win, &
               & (j - n/2 - 1) * self%dt, cc(j,i)
       end do
       
       write(io2)(i-1)*self%t_step + 0.5d0 * self%t_win, cc_max

       i1 = i1 + self%n_step
       i2 = i2 + self%n_step
    end do
    close(io)
    close(io2)

    !call self%eval_correlogram(cc, sta1, sta2)
    

    return 
  end subroutine correlator_run_cross_corr

  !-------------------------------------------------------------------------

  subroutine correlator_eval_correlogram(self, cc, sta1, sta2)
    class(correlator), intent(inout) :: self
    double precision, intent(in) :: cc(self%n, self%n_win)
    character(line_max), intent(in) :: sta1, sta2
    double precision :: tmp(self%n_win * self%n), cc_thred
    integer :: i, io, ierr
    character(line_max) :: out_file
    

    
    tmp = [cc]
    call quick_sort(tmp, 1, self%n*self%n_win)
    
    cc_thred = tmp(int(self%n*self%n_win*self%alpha))

    
    write(out_file,'(A)')trim(sta1) // "." //trim(sta2) // ".thred"
    open(newunit=io, file=out_file, status="replace", form='formatted', &
         & iostat=ierr)
    if (ierr /= 0) then
       error stop "ERROR: cannot creat threshold output file" 
    end if
    
    write(io,*)cc_thred
    close(io)


    return 
  end subroutine correlator_eval_correlogram
    
    
  !-------------------------------------------------------------------------
  
  subroutine correlator_measure_lag_time(self, x1, x2, l1, l2, i_corr_max, &
       & sta1, sta2)
    class(correlator), intent(inout) :: self
    double precision, intent(in) :: x1(self%n), x2(self%n)
    double precision, intent(in) :: l1, l2
    integer, intent(in) :: i_corr_max
    character(line_max), intent(in) :: sta1, sta2
    character(line_max) :: out_file
    double precision :: x_sum(self%n), sigma
    double precision :: err(self%n), acf(self%n), l(self%n, self%n)
    integer :: j, io, i
    double precision, parameter :: pi = acos(-1.d0)
    
    ! Note : i_corr_max = 1 corresponds to zero-lag time

    x_sum = 0.d0
    err = 0.d0
    !write(*,*)"I_corr_max=", i_corr_max
    do j = 1, self%n
       if (j+i_corr_max < 1 .or. j+i_corr_max > self%n) cycle
       x_sum(j) = (x1(j)/ l1 + x2(j+i_corr_max) / l2) * 0.5d0
       err(j) = x_sum(j) - x1(j) / l1
    end do
    err = err / sqrt(sum(err**2))
    sigma = sqrt(sum(err**2) / self%n)

    self%c_tmp = (0.d0, 0.d0) ! Not sure for reason, but seems necessary
    self%r_tmp = err
    call fftw_execute_dft_r2c(self%plan_r2c, self%r_tmp, self%c_tmp)
    
    self%c2_tmp = conjg(self%c_tmp)*self%c_tmp
    call fftw_execute_dft_c2r(self%plan_c2r, self%c2_tmp, self%r2_tmp)
    acf = 0.d0
    do j = 1, self%n
       if (self%r2_tmp(j) < 0.d0) exit
       acf(j) = self%r2_tmp(j) / self%n
    end do

    do concurrent (i=1:self%n)
       do j = 1, self%n
          l(j,i) = acf(abs(j-i) + 1) * sigma * sigma
       end do
    end do

    if (self%debug) then
       write(out_file,'(a,I6.6,a)') trim(sta1) // "." // trim(sta2) // ".", &
            self%n_detect, ".dat"
       open(newunit=io, file=out_file, status="replace")
       do j = 1, self%n
          write(io,*)j, x1(j) / l1
       end do
       write(io,*)
       write(io,*)
       do j = 1, self%n
          write(io,*)j-i_corr_max, x2(j) / l2
       end do
       write(io,*)
       write(io,*)
       do j = abs(i_corr_max) + 1, self%n
          write(io,*)j, x_sum(j)
       end do
       write(io,*)
       write(io,*)
       do j = 1, self%n
          write(io,*)j, acf(j) * sigma 
       end do
       
       close(io)
    end if


    return 
  end subroutine correlator_measure_lag_time

  !-------------------------------------------------------------------------
  
end module cls_correlator
