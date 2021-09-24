module cls_convertor
  use cls_c3_data
  use cls_line_text
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
     procedure :: taper   => convertor_taper
     procedure :: rectangle_smoothing => convertor_rectangle_smoothing
  end type convertor
  
  interface convertor
     module procedure init_convertor
  end interface convertor
     
contains
  
  !---------------------------------------------------------------------
  
  type(convertor) function init_convertor(n_cmps, t_win, &
       & filenames) result(self)
    integer, intent(in) :: n_cmps
    double precision, intent(in) :: t_win
    character(line_max), intent(in) :: filenames(:,:)
    integer :: n_files
    
    self%t_win = t_win
    self%t_window_interval = t_win * 0.5d0
    self%n_cmps = n_cmps
    n_files = size(filenames(:,1))
    allocate(self%filenames(n_files, self%n_cmps))
    self%filenames = filenames
    
    return 
  end function init_convertor
  
  !---------------------------------------------------------------------

  subroutine convertor_convert(self)
    class(convertor), intent(inout) :: self
    double precision, allocatable :: tmp(:, :), processed(:,:)
    integer :: n2, n4, io, it, j, irow
    logical :: first_window
    logical, parameter :: debug = .false.
    
    self%i_read_file = 1
    self%c3 = c3_data(files=self%filenames(1, 1:self%n_cmps))
    self%dt = self%c3%get_dt()
    self%n = nint(self%t_win / self%dt)
    self%c3_out = c3_data(dt=self%dt, n_cmps=self%n_cmps)
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
    write(*,*)self%filenames
    first_window = .true.
    if (debug) open(newunit=io, file="windows.txt", status="replace")
    do 
       if (self%c3%get_n_smp() < self%n) then
          self%i_read_file = self%i_read_file + 1
          write(*,*)"i_read_file=", self%i_read_file
          if (self%i_read_file > size(self%filenames(:,1))) exit

          write(*,*)"enqueued from <", &
               & trim(self%filenames(self%i_read_file,1)), " / ",&
               & trim(self%filenames(self%i_read_file,2)), " / ",&
               & trim(self%filenames(self%i_read_file,3)), "> "
          call self%c3%enqueue_from_files(&
               & self%filenames(self%i_read_file,:))
          
          write(*,*)"n_stored = ", self%c3%get_n_smp()
               
          
       end if
       if (.not. first_window) then
          tmp(1:n2,1:self%n_cmps) = tmp(n2+1:self%n,1:self%n_cmps)
          tmp(n2+1:self%n,1:self%n_cmps) = self%c3%dequeue_data(n2)
       else
          tmp(1:self%n,1:self%n_cmps) = self%c3%dequeue_data(self%n)
       end if

       write(*,*)"dequeued: n_sotred = ", self%c3%get_n_smp()

       
       ! Main 
       processed(1:self%n, 1:self%n_cmps) = tmp(1:self%n, 1:self%n_cmps)
       block 
         integer :: i_cmp
         do i_cmp = 1, self%n_cmps
            processed(1:self%n, i_cmp) = &
                 & self%detrend(processed(1:self%n,i_cmp))
            processed(1:self%n, i_cmp) = &
                 & self%taper(processed(1:self%n,i_cmp))
            processed(1:self%n, i_cmp) = &
                 & self%envelope(processed(1:self%n, i_cmp))
            processed(1:self%n, i_cmp) = &
                 & self%rectangle_smoothing(processed(1:self%n, i_cmp), &
                 & int(2.5d0 / self%dt))
            processed(1:self%n, i_cmp) = &
                 & self%rectangle_smoothing(processed(1:self%n, i_cmp), &
                 & int(2.5d0 / self%dt))
        end do
       end block
       
       ! output
       if (.not. first_window) then
          call self%c3_out%enqueue_data(&
               & processed(n4+1:self%n-n4,1:self%n_cmps))
       else
          call self%c3_out%enqueue_data(&
               & processed(1:self%n-n4,1:self%n_cmps))
       end if
       
       if (debug)  then
          block
            double precision :: amp_fac = 4000.d0
            do j = 1, self%n
               write(io,*)(it + j) * self%dt, tmp(j, 1) / amp_fac + irow
            end do
            write(io,*)
            do j = 1, self%n
               write(io,*)(it + j) * self%dt, processed(j, 1) / amp_fac + irow
            end do
            it = it + n2
            irow = irow + 1
            write(io,*)
          end block
       end if
       
       first_window = .false.
    end do
    if (debug) close(io)
    
    ! End 
    call self%c3_out%enqueue_data(&
         & processed(self%n-n4+1:self%n,1:self%n_cmps))


    block 
      integer :: n_out, i, n_fac, i_cmp
      double precision, allocatable :: x_out(:,:)
      character(line_max) :: out_file
      

      n_out = self%c3_out%get_n_smp()
      allocate(x_out(1:n_out,1:self%n_cmps))
      x_out = self%c3_out%get_data()
      n_fac = nint(1.d0 / self%dt / self%n_sps) ! decimate
      do i_cmp = 1, self%n_cmps
         write(out_file,'(a,I1,a)')'output_', i_cmp, '.dat' 
         open(newunit=io, file=out_file, access='stream', &
              & form='unformatted', status='replace')
         do i = 1, n_out
            if (mod(i,n_fac) == 1) then
               write(io) (i-1) * self%dt, x_out(i, i_cmp)
            end if
         end do
         close(io)
      end do

    end block


    return 
  end subroutine convertor_convert

  !---------------------------------------------------------------------
  
  function convertor_rectangle_smoothing(self, x, n_half) result(x_out)
    class(convertor), intent(inout) :: self
    double precision, intent(in) :: x(self%n)
    integer, intent(in) :: n_half
    double precision :: x_out(self%n), s
    integer :: i

    s = sum(x(1:n_half))
    do i = 1, n_half
       s = s + x(n_half+i)
       x_out(i) = s / dble(n_half+i)
    end do
    do i = n_half+1, self%n-n_half
       s = s + x(n_half+i) - x(i - n_half)
       x_out(i) = s / dble(2*n_half + 1)
    end do
    do i = self%n-n_half+1, self%n
       s = s - x(i-n_half)
       x_out(i) = s / (self%n + n_half - i + 1)
    end do


    return 
  end function convertor_rectangle_smoothing

  !---------------------------------------------------------------------
  
  function convertor_band_pass(self, x) result(x_out)
    class(convertor), intent(inout) :: self
    complex(kind(0d0)), intent(in) :: x(self%n)
    complex(kind(0d0)) :: x_out(self%n)
    double precision :: df, flt
    integer :: if1, if2, if3, if4, i
    double precision, parameter :: pi = acos(-1.d0)

    df = 1.d0 / (self%n * self%dt)
    if1 = nint(self%f1 / df) + 1
    if2 = nint(self%f2 / df) + 1
    if3 = nint(self%f3 / df) + 1
    if4 = nint(self%f4 / df) + 1
    
    do i = 1, self%n
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

  function convertor_envelope(self, x) result(env)
    class(convertor), intent(inout) :: self
    double precision, intent(in) :: x(self%n)
    double precision :: env(self%n)
    



    self%r_tmp(1:self%n) = x(1:self%n)
    call fftw_execute_dft_r2c(self%plan_r2c, self%r_tmp, self%c_tmp)
    self%c_tmp = self%band_pass(self%c_tmp)
    self%c_in_tmp(1) = (0.d0, 0.d0)
    self%c_in_tmp(2:self%n/2+1) = 2.d0 * self%c_tmp(2:self%n/2+1)
    self%c_in_tmp(self%n/2+2:self%n) = (0.d0, 0.d0)
    
    call fftw_execute_dft(self%plan_c2c, self%c_in_tmp, self%c_out_tmp)
    env(1:self%n) = abs(self%c_out_tmp(1:self%n) / self%n) 
    
    return 
  end function convertor_envelope

  !---------------------------------------------------------------------

  function convertor_detrend(self, x) result(x_out)
    class(convertor), intent(in) :: self
    double precision, intent(in) :: x(1:self%n)
    double precision :: x_out(1:self%n)
    double precision :: sxy, sxx, x_mean, y_mean, a, b
    integer :: i

    x_mean = 0.5d0 * (1.d0 + dble(self%n))
    y_mean = sum(x) / dble(self%n)
    sxx = 0.d0
    sxy = 0.d0

    do i = 1, self%n
       sxx = sxx + (dble(i) - x_mean)**2
       sxy = sxy + (x(i) - y_mean) * (dble(i) - x_mean)
    end do

    a = sxy / sxx
    b = y_mean - a * x_mean
    
    do concurrent (i = 1:self%n)
       x_out(i) = x(i) - (a * dble(i) + b)
    end do
    
    return 
  end function convertor_detrend

  !---------------------------------------------------------------------

  function convertor_taper(self, x) result(x_out)
    class(convertor), intent(inout) :: self
    double precision, intent(in) :: x(1:self%n)
    double precision :: x_out(1:self%n)
    double precision, parameter :: pcnt = 0.05d0, pi = acos(-1.d0)
    integer :: nleng, i
    double precision :: fac
    

    nleng = int(self%n * pcnt)
    x_out = x
    do concurrent (i = 1:nleng)
       fac = 0.5d0 * (1.d0 - cos((i - 1) * pi / nleng))
       x_out(i) = x(i) * fac
       x_out(self%n - i + 1) = x(self%n - i + 1) * fac
    end do
   
    return 
  end function convertor_taper
  
end module cls_convertor
