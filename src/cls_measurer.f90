module cls_measurer
  use mod_mpi
  use cls_line_text, only: line_max
  use mod_signal_process, only: apply_taper
  use, intrinsic :: iso_fortran_env, only: iostat_end
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'


  type measurer
     private

     integer :: n_sta
     integer :: n_pair
     integer :: n_win
     integer :: n_step
     integer :: n
     integer :: n_detected
     integer, allocatable :: win_id(:)
     double precision :: dt
     double precision :: t_win
     double precision :: t_step
     character(line_max), allocatable :: station_names(:)
     double precision, allocatable :: sta_y(:)
     double precision, allocatable :: sta_x(:)
     
     character(line_max), allocatable :: pair(:,:)
     logical :: verb

     double precision, allocatable :: cc_thred(:)
     integer :: n_pair_thred
     logical, allocatable :: detected(:,:)

     type(C_PTR)                            :: plan_r2c, plan_c2r
     real(C_DOUBLE), allocatable            :: r_tmp(:), r2_tmp(:)
     complex(C_DOUBLE_COMPLEX), allocatable :: c_tmp(:), c2_tmp(:)

     
   contains
     procedure :: check_files => measurer_check_files
     procedure :: scan_cc     => measurer_scan_cc
     procedure :: measure_lag_time => measurer_measure_lag_time
     procedure :: optimize_cc => measurer_optimize_cc
     procedure :: optimize_amp => measurer_optimize_amp
  end type measurer

  
  interface measurer
     module procedure init_measurer
  end interface measurer

contains

  !-------------------------------------------------------------------
  
  type(measurer) function init_measurer(station_names, sta_x, sta_y, &
       & t_win, t_step, n_pair_thred, verb) result(self)
    character(line_max), intent(in) :: station_names(:)
    double precision, intent(in) :: sta_x(:), sta_y(:)
    double precision, intent(in) :: t_win, t_step
    integer, intent(in) :: n_pair_thred
    logical, intent(in) :: verb
    integer :: i, j, k
    
    self%n_sta = size(station_names)
    allocate(self%station_names(self%n_sta))
    allocate(self%sta_x(self%n_sta))
    allocate(self%sta_y(self%n_sta))
    self%station_names = station_names
    self%sta_x = sta_x
    self%sta_y = sta_y
    self%n_pair = self%n_sta * (self%n_sta - 1) / 2
    self%verb = verb
    self%t_win = t_win
    self%t_step = t_step
    self%n_pair_thred = n_pair_thred
    
    allocate(self%pair(2, self%n_pair))
    k = 0
    do i = 1, self%n_sta - 1
       do j = i + 1, self%n_sta
          k = k + 1
          self%pair(1, k) = station_names(i)
          self%pair(2, k) = station_names(j)
       end do
    end do

    allocate(self%cc_thred(self%n_pair))

    call self%check_files()

    allocate(self%r_tmp(self%n))
    allocate(self%r2_tmp(self%n))
    allocate(self%c_tmp(self%n))
    allocate(self%c2_tmp(self%n))
    self%plan_r2c = fftw_plan_dft_r2c_1d(self%n, self%r_tmp, self%c_tmp, &
         & FFTW_ESTIMATE)
    self%plan_c2r = fftw_plan_dft_c2r_1d(self%n, self%c2_tmp, self%r2_tmp, &
         & FFTW_ESTIMATE)

    return 
  end function init_measurer
  
  !-------------------------------------------------------------------
  
  subroutine measurer_check_files(self)
    class(measurer), intent(inout) :: self
    integer :: i, ierr, io, n_smp, prev_n_smp
    character(line_max) :: target_file, sta1, sta2
    logical :: is_ok
    double precision :: t1, t2, prev_dt, cc


    do i = 1, self%n_sta
       write(target_file, '(a)') trim(self%station_names(i)) // ".merged.env"
       if (self%verb) then
          print *, trim(target_file)
       end if
       inquire(file=target_file, EXIST=is_ok)
       if (.not. is_ok) then
          write(0,*) "ERROR: " // trim(target_file) // "does not exsit"
          error stop 
       end if
       ! Get dt
       open(newunit=io, file=target_file, status="old", form="unformatted", &
            & access="direct", recl=8)
       read(io, rec=1) t1
       read(io, rec=3) t2
       close(io)
       self%dt = t2 - t1
       if (i > 1 .and. abs(self%dt - prev_dt) > 1.0e-8) then
          error stop "invalid delta in envelope file"
       end if
          
       prev_dt = self%dt
       
    end do
    self%n_step = nint(self%t_step / self%dt)
    self%n      = nint(self%t_win  / self%dt)

    do i = 1, self%n_pair
       sta1 = self%pair(1, i)
       sta2 = self%pair(2, i)

       ! Max corr file
       write(target_file, '(a)') &
            & trim(sta1) // "." // trim(sta2) // ".max_corr"
       if (self%verb) then
          print *, trim(target_file)
       end if
       inquire(file=target_file, EXIST=is_ok)
       if (.not. is_ok) then
          write(0,*) "ERROR: " // trim(target_file) // "does not exsit"
          error stop 
       end if
       
       open(newunit=io, file=target_file, status="old", form="unformatted", &
            & access="stream")
       self%n_win = 0
       do 
          read(io, iostat=ierr)t1, cc
          if (ierr == iostat_end) then
             exit
          end if
          self%n_win = self%n_win + 1
       end do
       close(io)

       if (i > 1 .and. prev_n_smp /= self%n_win) then
          error stop "ERROR: invalid n_win in corr_max file"
       end if
       prev_n_smp = self%n_win
       
       ! Threshold file
       write(target_file, '(a)') &
            & trim(sta1) // "." // trim(sta2) // ".thred"
       if (self%verb) then
          print *, trim(target_file)
       end if
       inquire(file=target_file, EXIST=is_ok)
       if (.not. is_ok) then
          write(0,*) "ERROR: " // trim(target_file) // "does not exsit"
          error stop 
       end if
    end do
    allocate(self%detected(self%n_win, self%n_pair))
    
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    return 
  end subroutine measurer_check_files
  
  !-------------------------------------------------------------------

  subroutine measurer_scan_cc(self)
    class(measurer), intent(inout) :: self
    integer :: id1, id2, i, io, i_win, ierr, rank
    character(line_max) :: sta1, sta2, thred_file, max_corr_file
    character(line_max) :: detected_win_file
    double precision :: cc, t
    logical, allocatable :: ltmp(:,:)

    self%detected = .false.
    call get_mpi_task_id(self%n_pair, id1, id2)
    do i = id1, id2
       sta1 = self%pair(1, i)
       sta2 = self%pair(2, i)
       thred_file = trim(sta1) // "." // trim(sta2) // ".thred"
       open(newunit=io, file=thred_file, form="formatted", status="old")
       read(io,*) self%cc_thred(i)
       !print *, trim(sta1), " -- ", trim(sta2), "  : CC threshold = ", &
       !     & self%cc_thred(i)
       close(io)

       max_corr_file = trim(sta1) // "." // trim(sta2) // ".max_corr"
       open(newunit=io, file=max_corr_file, form="unformatted", access="stream")
       do i_win = 1, self%n_win
          read(io)t, cc
          if (cc >= self%cc_thred(i)) then
             self%detected(i_win, i) = .true.
          end if
       end do
       close(io)
    end do
    
    
    allocate(ltmp(self%n_win, self%n_pair))
    ltmp = self%detected
    call mpi_allreduce(ltmp, self%detected, size(ltmp), MPI_LOGICAL, MPI_LOR, &
         & MPI_COMM_WORLD, ierr)
    call get_mpi_task_id(self%n_win, id1, id2)

    self%win_id = [integer :: ]
    do i_win = 1, self%n_win
       if (count(self%detected(i_win,:)) > self%n_pair_thred) then
          !print *, "detected! ", count(self%detected(i_win,:)), i_win
          self%win_id = [self%win_id, i_win]
       end if
    end do
    self%n_detected = size(self%win_id)

    if (self%verb) then
       print *, "# of detected events: ", self%n_detected, "out of ", self%n_win
    end if
    
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
    if (rank == 0) then
       detected_win_file = "detected_win.dat"
       open(newunit=io, file=detected_win_file, status="replace", iostat=ierr,&
       & form="formatted")
       do i = 1, self%n_detected
          write(io,*)self%win_id(i), (self%win_id(i)-1)*self%t_step + 0.5d0 * self%t_win
       end do
       close(io)
    end if

    return 
  end subroutine measurer_scan_cc

  !-------------------------------------------------------------------

  subroutine measurer_measure_lag_time(self)
    class(measurer), intent(inout) :: self
    integer :: i1, i2, id, i, io, ierr, j, ista, j1, j2, io2
    double precision, allocatable :: tmp(:, :)
    double precision :: fac, t(self%n_sta), t_stdv(self%n_sta)
    double precision :: amp(self%n_sta), amp_stdv(self%n_sta)
    character(line_max) :: trace_file, env_file, time_file
    
    !print *, self%win_id
    
    allocate(tmp(self%n, self%n_sta))
    call get_mpi_task_id(self%n_detected, i1, i2)

    do i = i1, i2
       id = self%win_id(i)
       j1 = (id - 1) * self%n_step + 1
       j2 = j1 + self%n - 1
       !print *, j1, j2
       write(trace_file,'(A,I5.5,A)') "trace.", id, ".dat"
       open(newunit=io, file=trace_file, form="formatted", status="replace", &
            & iostat=ierr)
       
       
       if (ierr /= 0) error stop "cannot create trace file"
       
       do ista = 1, self%n_sta
       
          env_file = trim(self%station_names(ista)) // ".merged.env"
          open(newunit=io2, file=env_file, status="old", form="unformatted", &
               & access="direct", recl=8)
          do j = j1, j2
             read(io2, rec=2*j)tmp(j-j1+1, ista)
          end do
          close(io2)
          tmp(:,ista) = tmp(:,ista) - sum(tmp(:,ista))/self%n
          fac = max(-minval(tmp(:,ista)),maxval(tmp(:,ista)))
          do j = 1, self%n
             write(io,*)(j - 1) * self%dt, tmp(j,ista) / fac + ista
          end do
          write(io,*)
       end do
       call self%optimize_cc(tmp, t, t_stdv)
       !call self%optimize_amp(tmp, t, amp, amp_stdv)

       write(io,*)
       
       do ista = 1, self%n_sta
          do j = 1, self%n
             fac = max(-minval(tmp(:,ista)),maxval(tmp(:,ista)))
             write(io,*)(j - 1) * self%dt - t(ista), tmp(j,ista) / fac + ista
          end do
          write(io,*)
       end do

       write(io,*)
       
       call self%optimize_amp(tmp, t, amp, amp_stdv)

       ! Unnormalized envelopes
       do ista = 1, self%n_sta
          do j = 1, self%n
             write(io,*)(j - 1) * self%dt - t(ista), tmp(j,ista)
          end do
          write(io,*)
       end do
       write(io,*)

       ! Normalized envelopes
       do ista = 1, self%n_sta
          do j = 1, self%n
             write(io,*)(j - 1) * self%dt - t(ista), &
                  & tmp(j,ista) / exp(amp(ista))
          end do
          write(io,*)
       end do
       write(io,*)
       
       close(io)

       write(time_file,'(A,I5.5,A)')"rel_time.", id, ".dat"
       open(newunit=io, file=time_file, form="formatted", status="replace", &
            & iostat=ierr)
       do ista = 1, self%n_sta
          write(io, *) self%sta_x(ista), self%sta_y(ista), t(ista), t_stdv(ista), &
               & amp(ista), amp_stdv(ista)
       end do
       close(io)
       
    end do
    
    return 
  end subroutine measurer_measure_lag_time


  !-------------------------------------------------------------------
  subroutine measurer_optimize_amp(self, x, t, amp, amp_stdv)
    class(measurer), intent(inout) :: self
    double precision, intent(in) :: x(self%n, self%n_sta)
    double precision, intent(in) :: t(self%n_sta)
    double precision, intent(out) :: amp(self%n_sta)
    double precision, intent(out) :: amp_stdv(self%n_sta)
    double precision :: x2(self%n, self%n_sta), sxy
    double precision :: rel_log_amp(self%n_sta, self%n_sta), sxx(self%n_sta)
    integer :: i, j, it

    
    x2 = 0.d0
    rel_log_amp = 0.d0
    do i = 1, self%n_sta
       it = nint(t(i) / self%dt)
       do j = 1, self%n
          if (j+it < 1 .or. j+it > self%n) cycle
          x2(j, i) = x(j+it, i)
       end do
       sxx(i) = sum(x2(:,i)**2)
    end do
    
    do i = 1, self%n_sta - 1
       do j = i + 1, self%n_sta
          sxy = sum(x2(:,i) * x2(:,j))
          if (sxy < 0.d0) then
             amp = 0.d0
             amp_stdv = 0.d0
             return
          end if
          rel_log_amp(i,j) = log( sxy / sxx(i))
          rel_log_amp(j,i) = - rel_log_amp(i, j)
       end do
    end do
    
    amp = 0.d0
    do i = 1, self%n_sta
       do j = 1, self%n_sta
          amp(i) = amp(i) - rel_log_amp(i,j)
       end do
       amp(i) = amp(i) / self%n_sta
    end do

    amp_stdv = 0.d0
    do i = 1, self%n_sta
       do j = 1, self%n_sta
          if (i==j) cycle 
          amp_stdv(i) = amp_stdv(i) + (amp(j) - amp(i) - rel_log_amp(i,j))**2
          
       end do
    end do
    amp_stdv = sqrt(amp_stdv /(self%n_sta - 2)) 
    
    return 
  end subroutine measurer_optimize_amp

  !-------------------------------------------------------------------
  
  subroutine measurer_optimize_cc(self, x, t, t_stdv)
    class(measurer), intent(inout) :: self
    double precision, intent(in) :: x(self%n, self%n_sta)
    double precision, intent(out) :: t(self%n_sta)
    double precision, intent(out) :: t_stdv(self%n_sta)
    double precision :: a(self%n_sta, self%n_sta)
    double precision :: lag_t(self%n_sta, self%n_sta)
    double precision :: tmp(self%n_sta - 1)
    double precision :: l
    complex(kind(0d0)) :: cx(self%n, self%n_sta)
    integer :: i, j, ilag_t


    
    do i = 1, self%n_sta
       l = sum(x(1:self%n,i)**2)
       self%r_tmp(1:self%n) = apply_taper(self%n, x(1:self%n, i))
       self%r_tmp(1:self%n) = self%r_tmp(1:self%n) / l
       self%c_tmp = (0.d0, 0.d0)
       call fftw_execute_dft_r2c(self%plan_r2c, self%r_tmp, self%c_tmp)
       cx(1:self%n, i) = self%c_tmp(1:self%n)
    end do

    lag_t = 0.d0
    do i = 1, self%n_sta - 1
       do j = i + 1, self%n_sta
          self%c2_tmp(1:self%n) = conjg(cx(1:self%n,i)) * cx(1:self%n,j)
          self%r2_tmp = 0.d0
          call fftw_execute_dft_c2r(self%plan_c2r, self%c2_tmp, self%r2_tmp)
          ilag_t = maxloc(self%r2_tmp, dim=1)
          if (ilag_t <= self%n/2) then
             lag_t(i, j) = (ilag_t - 1) * self%dt
          else
             lag_t(i, j) = (ilag_t - self%n - 1) * self%dt
          end if
          lag_t(j, i) = -lag_t(i, j)
          !print *, lag_t(i, j)
       end do
    end do

    t = 0.d0
    do i = 1, self%n_sta
       do j = 1, self%n_sta
          t(i) = t(i) - lag_t(i,j)
       end do
       t(i) = t(i) / self%n_sta
    end do
    

    ! standard deviation
    t_stdv = 0.d0
    do i = 1, self%n_sta
       do j = 1, self%n_sta
          if (i==j) cycle 
          t_stdv(i) = t_stdv(i) + (t(j) - t(i) - lag_t(i,j))**2
       end do
    end do
    t_stdv = sqrt(t_stdv /(self%n_sta - 2))
    
    return 
  end subroutine measurer_optimize_cc
  
end module cls_measurer
