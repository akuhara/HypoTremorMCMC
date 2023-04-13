module cls_convertor
  use cls_c3_data, only: c3_data
  use cls_line_text, only: line_max
  use mod_signal_process, only: apply_taper
  use, intrinsic :: iso_c_binding
  implicit none 
  include 'fftw3.f03'

  type convertor
     private
     type(c3_data)              :: c3
     type(c3_data)              :: c3_out
     double precision           :: t_win
     double precision           :: t_window_interval
     double precision           :: dt
     integer                    :: n
     integer                    :: n_cmps
     character(line_max), allocatable :: filenames(:,:)
     character(line_max)        :: station_name
     character(line_max), allocatable :: comps(:)
     double precision, allocatable :: sta_amp_fac(:)
     
     ! Band pass
     double precision :: f1 = 1.d0
     double precision :: f2 = 3.d0
     double precision :: f3 = 8.d0
     double precision :: f4 = 10.d0
     
     ! resampling
     integer :: n_sps = 1
     
     ! Work space for FFTW
     type(C_PTR)                            :: plan_r2c, plan_c2c, plan_c2r
     real(C_DOUBLE), allocatable            :: r_tmp(:), r2_tmp(:)
     complex(C_DOUBLE_COMPLEX), allocatable :: c_tmp(:), c2_tmp(:)
     complex(C_DOUBLE_COMPLEX), allocatable :: c_out_tmp(:)
     complex(C_DOUBLE_COMPLEX), allocatable :: c_in_tmp(:)
     
     !
     integer :: i_read_file

   contains
     procedure :: convert => convertor_convert
     procedure :: band_pass => convertor_band_pass
     procedure :: envelope => convertor_envelope
     procedure :: detrend => convertor_detrend
     !procedure :: taper   => convertor_taper
     procedure :: rectangle_smoothing => convertor_rectangle_smoothing
     procedure :: get_c3_out => convertor_get_c3_out
     procedure :: merge_components => convertor_merge_components
  end type convertor
  
  interface convertor
     module procedure init_convertor
  end interface convertor
     
contains
  
  !---------------------------------------------------------------------
  
  type(convertor) function init_convertor(t_win, &
       & filenames, comps, station_name, sta_amp_fac) result(self)
    double precision, intent(in) :: t_win
    character(line_max), intent(in) :: filenames(:,:)
    character(line_max), intent(in) :: station_name
    character(line_max), intent(in) :: comps(:)
    double precision, intent(in) :: sta_amp_fac(:)
    integer :: n_files
    

    self%t_win = t_win
    self%t_window_interval = t_win * 0.5d0
    n_files = size(filenames(:,1))
    self%n_cmps = size(comps)
    allocate(self%filenames(n_files, self%n_cmps))
    self%filenames = filenames
    self%station_name = station_name
    allocate(self%comps(self%n_cmps))
    self%comps = comps
    allocate(self%sta_amp_fac(self%n_cmps))
    self%sta_amp_fac = sta_amp_fac
    
    
    return 
  end function init_convertor
  
  !---------------------------------------------------------------------

  subroutine convertor_convert(self)
    class(convertor), intent(inout) :: self
    double precision, allocatable :: tmp(:, :), processed(:,:), merged(:, :)
    integer :: n2, n4, io, it, j, irow, i_cmp, n_len, n_fac, n_fac_mod
    integer :: i_start, i_end
    logical :: first_flag, last_flag
    logical, parameter :: debug = .true.
    

    
    ! Read the ealiest time series
    ! Time increment (delta) is read from SAC hedaer here and is stored in
    ! c3%dt 
    self%c3 = c3_data(files=self%filenames(1, 1:self%n_cmps)) 

    self%i_read_file = 1 ! Number of files already read
    
    ! Get time incremnt  
    self%dt = self%c3%get_dt() ! From SAC header
    self%n = nint(self%t_win / self%dt) ! Number of elements in a 
                                        !  single time window
    n_fac = nint(1.d0 / self%dt / self%n_sps) ! decimate
    n_fac_mod = 0
    self%c3_out = c3_data(dt=self%dt * n_fac, n_cmps=1)
    n2 = self%n / 2
    if (n2 + n2 /= self%n) then
       error stop "ERROR: n2 + n2 /= self%n"
    end if
    n4 = self%n / 4
    if (n4 * 4 /= self%n) then
       error stop "ERROR: n4 * 4 /= self%n"
    end if
    
    allocate(tmp(self%n, self%n_cmps))
    allocate(processed(self%n, self%n_cmps))
    allocate(merged(self%n, 1))
    allocate(self%r_tmp(self%n))
    allocate(self%c_tmp(self%n))
    allocate(self%c_in_tmp(self%n))
    allocate(self%c_out_tmp(self%n))
    allocate(self%r2_tmp(self%n))
    allocate(self%c2_tmp(self%n))
    self%plan_r2c = fftw_plan_dft_r2c_1d(self%n, self%r_tmp, self%c_tmp, &
         & FFTW_ESTIMATE)
    self%plan_c2r = fftw_plan_dft_c2r_1d(self%n, self%c2_tmp, self%r2_tmp, &
         & FFTW_ESTIMATE)
    self%plan_c2c = fftw_plan_dft_1d(self%n, self%c_in_tmp, self%c_out_tmp, &
         & FFTW_BACKWARD, FFTW_ESTIMATE)

    
    tmp(1:self%n, 1:self%n_cmps) = 0.d0
    it = 0
    irow = 0
    
    write(*,*)"n_stored = ", self%c3%get_n_smp()
    ! n_smp is a total number of elemets stored in c3_data


    
    write(*,*)self%filenames
    first_flag = .true.
    last_flag  = .false.
    if (debug) open(newunit=io, file="windows.txt", status="replace")
    
    processing: do 

       ! << Append data, if necessary >>
       if (self%c3%get_n_smp() < self%n) then
          read_file: do 
             ! Read new SAC files
             if (self%i_read_file + 1 > size(self%filenames(:,1))) exit read_file
             call self%c3%enqueue_from_files(&
                  & self%filenames(self%i_read_file + 1,:))
             self%i_read_file = self%i_read_file + 1
             if (self%c3%get_n_smp() >= self%n) exit read_file
          end do read_file
       end if

       ! << Check length >>
       if (self%c3%get_n_smp() >= n2) then
          n_len = self%n
       else
          last_flag = .true.
          n_len = n2 + self%c3%get_n_smp()
       end if
       
       
       ! << Processing >> 
       if (.not. last_flag) then

          ! Dequeue data
          if (.not. first_flag .and. .not. last_flag) then
             
             ! * Get half and slide
             tmp(1:n2,1:self%n_cmps) = tmp(n2+1:self%n,1:self%n_cmps)
             tmp(n2+1:self%n,1:self%n_cmps) = self%c3%dequeue_data(n2)
          else if (first_flag) then
             
             ! * Initial dequeue
             tmp(1:self%n,1:self%n_cmps) = self%c3%dequeue_data(self%n)

          else
             ! This block seems not necessary (2023/4/13)
             ! * Last dequeue
             tmp(1:n2,1:self%n_cmps) = tmp(n2+1:self%n,1:self%n_cmps)
             tmp(n2+1:n_len,1:self%n_cmps) = self%c3%dequeue_data(n_len-n2)
          end if
          
          
          ! Main signal processing
          processed(1:n_len, 1:self%n_cmps) = tmp(1:n_len, 1:self%n_cmps)
          do i_cmp = 1, self%n_cmps
             processed(1:n_len, i_cmp) = &
                  & self%detrend(n_len, processed(1:n_len,i_cmp))
             processed(1:n_len, i_cmp) = &
                  & apply_taper(n_len, processed(1:n_len,i_cmp))
             processed(1:n_len, i_cmp) = &
                  & self%envelope(n_len, processed(1:n_len, i_cmp))
             processed(1:n_len, i_cmp) = &
                  & self%rectangle_smoothing(n_len, processed(1:n_len, i_cmp), &
                  & int(1.5d0 / self%dt)) ! 2.5
             processed(1:n_len, i_cmp) = &
                  & self%rectangle_smoothing(n_len, processed(1:n_len, i_cmp), &
                  & int(1.5d0 / self%dt))
          end do
          merged(1:n_len, 1) = self%merge_components(n_len, self%n_cmps, &
               & processed(1:n_len, 1:self%n_cmps))
       end if
       
       ! << Make output >> ! This section has an issue of expanding memory size
       
       if (.not. first_flag .and. .not. last_flag) then
          i_start = n4 + 1
          i_end   = self%n - n4
       else if (first_flag) then
          i_start = 1
          i_end   = self%n - n4
       else
          i_start = n4 + 1
          i_end = n_len
       end if
       call self%c3_out%enqueue_data(merged(i_start+n_fac_mod:i_end:n_fac,:)) 
       n_fac_mod = n_fac - mod(i_end - i_start - n_fac_mod + 1, n_fac)
       if (n_fac_mod == n_fac) then
          n_fac_mod = 0
       end if

       ! << Debug output >>
       if (debug)  then
          block
            double precision :: amp_fac = 4000.d0
            do j = 1, self%n
               write(io,*)(it + j) * self%dt, tmp(j, 1) / amp_fac + irow
            end do
            write(io,*)
            !do j = 1, self%n
            !   write(io,*)(it + j) * self%dt, merged(j, 1) / amp_fac + irow
            !end do
            it = it + n2
            irow = irow + 1
            write(io,*)
          end block
       end if
       
       if (last_flag) then
          exit processing
       end if
    end do processing
    if (debug) close(io)
    
    ! Decimate
    block 
      integer :: n_out, i, n_fac, i_cmp
      double precision, allocatable :: x_out(:,:)
      double precision :: dt_out
      character(line_max) :: out_file

      n_out = self%c3_out%get_n_smp()
      allocate(x_out(1:n_out,1))
      x_out = self%c3_out%get_data()
      dt_out = self%c3_out%get_dt()
      write(out_file,'(a,a,a)')trim(self%station_name) // ".", &
           & "merged", '.env' 
      open(newunit=io, file=out_file, access='stream', &
           & form='unformatted', status='replace')
      do i = 1, n_out
         write(io) (i-1) * dt_out, x_out(i, 1)
      end do
      close(io)
    end block


    return 
  end subroutine convertor_convert

  !---------------------------------------------------------------------
  
  function convertor_rectangle_smoothing(self, n, x, n_half) result(x_out)
    class(convertor), intent(inout) :: self
    integer, intent(in) :: n 
    double precision, intent(in) :: x(n)
    integer, intent(in) :: n_half
    double precision :: x_out(n), s
    integer :: i

    s = sum(x(1:n_half))
    do i = 1, n_half
       s = s + x(n_half+i)
       x_out(i) = s / dble(n_half+i)
    end do
    do i = n_half+1, n-n_half
       s = s + x(n_half+i) - x(i - n_half)
       x_out(i) = s / dble(2*n_half + 1)
    end do
    do i = n-n_half+1, n
       s = s - x(i-n_half)
       x_out(i) = s / (n + n_half - i + 1)
    end do


    return 
  end function convertor_rectangle_smoothing

  !---------------------------------------------------------------------
  
  function convertor_band_pass(self, n, x) result(x_out)
    class(convertor), intent(inout) :: self
    integer, intent(in) :: n
    complex(kind(0d0)), intent(in) :: x(n)
    complex(kind(0d0)) :: x_out(n)
    double precision :: df, flt
    integer :: if1, if2, if3, if4, i
    double precision, parameter :: pi = acos(-1.d0)

    df = 1.d0 / (n * self%dt)
    if1 = nint(self%f1 / df) + 1
    if2 = nint(self%f2 / df) + 1
    if3 = nint(self%f3 / df) + 1
    if4 = nint(self%f4 / df) + 1
    
    do i = 1, n
       if (i < if1) then
          flt = 0.d0
       else if (i < if2) then
          flt = 0.5d0 * (1.d0 - cos((i-if1)*pi/(if2-if1)))
       else if (i < if3) then
          flt = 1.d0
       else if (i < if4) then
          flt = 0.5d0 * (1.d0 + cos((i-if3)*pi/(if4-if3)))
       else
          flt = 0.d0
       end if
       x_out(i) = x(i) * flt
    end do
    
    return 
  end function convertor_band_pass
  
  !---------------------------------------------------------------------

  function convertor_envelope(self, n, x) result(env)
    class(convertor), intent(inout) :: self
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision :: env(n)
    



    self%r_tmp(1:n) = x(1:n)
    call fftw_execute_dft_r2c(self%plan_r2c, self%r_tmp, self%c_tmp)
    self%c_tmp = self%band_pass(n, self%c_tmp)
    self%c_in_tmp(1) = (0.d0, 0.d0)
    self%c_in_tmp(2:n/2+1) = 2.d0 * self%c_tmp(2:n/2+1)
    self%c_in_tmp(n/2+2:n) = (0.d0, 0.d0)
    
    call fftw_execute_dft(self%plan_c2c, self%c_in_tmp, self%c_out_tmp)
    env(1:n) = abs(self%c_out_tmp(1:n) / n) 
    
    return 
  end function convertor_envelope

  !---------------------------------------------------------------------

  function convertor_detrend(self, n, x) result(x_out)
    class(convertor), intent(in) :: self
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision :: x_out(n)
    double precision :: sxy, sxx, x_mean, y_mean, a, b
    integer :: i
    
    x_mean = 0.5d0 * (1.d0 + dble(n))
    y_mean = sum(x) / dble(n)
    sxx = 0.d0
    sxy = 0.d0

    do i = 1, n
       sxx = sxx + (dble(i) - x_mean)**2
       sxy = sxy + (x(i) - y_mean) * (dble(i) - x_mean)
    end do

    a = sxy / sxx
    b = y_mean - a * x_mean
    
    do concurrent (i = 1:n)
       x_out(i) = x(i) - (a * dble(i) + b)
    end do
    
    return 
  end function convertor_detrend

  !!---------------------------------------------------------------------
  !
  !function convertor_taper(self, n, x) result(x_out)
  !  class(convertor), intent(inout) :: self
  !  integer, intent(in) :: n
  !  double precision, intent(in) :: x(1:n)
  !  double precision :: x_out(1:n)
  !  double precision, parameter :: pcnt = 0.05d0, pi = acos(-1.d0)
  !  integer :: nleng, i
  !  double precision :: fac
  !  
  !
  !  nleng = int(n * pcnt)
  !  x_out = x
  !  do concurrent (i = 1:nleng)
  !     fac = 0.5d0 * (1.d0 - cos((i - 1) * pi / nleng))
  !     x_out(i) = x(i) * fac
  !     x_out(n - i + 1) = x(n - i + 1) * fac
  !  end do
  ! 
  !  return 
  !end function convertor_taper
  !
  !!---------------------------------------------------------------------

  type(c3_data) function convertor_get_c3_out(self) result(c3_out)
    class(convertor), intent(in) :: self
    
    c3_out = self%c3_out
    
    return 
  end function convertor_get_c3_out
  
  !---------------------------------------------------------------------

  function convertor_merge_components(self, n, n_cmps, x) result(x_out)
    class(convertor), intent(inout) :: self
    integer, intent(in) :: n, n_cmps
    double precision, intent(in) :: x(1:n,1:n_cmps)
    double precision :: x_out(1:n)
    double precision :: s
    integer :: i, i_cmp
    
    do concurrent (i = 1:n)
       x_out(i) = sqrt(sum((x(i,1:n_cmps) * self%sta_amp_fac(1:n_cmps))**2))
    end do
          
    return 
  end function convertor_merge_components
  
  !---------------------------------------------------------------------
  
end module cls_convertor
