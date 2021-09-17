module cls_convertor
  use cls_c3_data
  use cls_line_text
  use, intrinsic :: iso_c_binding
  implicit none 
  include 'fftw3.f03'

  type convertor
     private
     type(c3_data), allocatable :: c3(:)
     type(c3_data), allocatable :: c3_out(:)
     double precision           :: t_win
     double precision           :: t_window_interval
     double precision           :: dt
     integer                    :: n
     integer                    :: n_sta
     integer                    :: n_cmps
     character(line_max), allocatable :: filenames(:,:,:)
     
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
     procedure :: envelope => convertor_envelope
     procedure :: detrend => convertor_detrend
     procedure :: taper   => convertor_taper
  end type convertor
  
  interface convertor
     module procedure init_convertor
  end interface convertor
     
contains
  
  !---------------------------------------------------------------------
  
  type(convertor) function init_convertor(n_sta, n_cmps, t_win, &
       & filenames) result(self)
    integer, intent(in) :: n_sta, n_cmps
    double precision, intent(in) :: t_win
    character(line_max), intent(in) :: filenames(:,:,:)
    integer :: n_files
    
    self%t_win = t_win
    self%t_window_interval = t_win * 0.5d0
    self%n_sta = n_sta
    self%n_cmps = n_cmps
    n_files = size(filenames(:,1,1))
    allocate(self%c3(self%n_sta))
    allocate(self%c3_out(self%n_sta))
    allocate(self%filenames(n_files, self%n_cmps, self%n_sta))
    self%filenames = filenames
    
    return 
  end function init_convertor
  
  !---------------------------------------------------------------------

  subroutine convertor_convert(self, i_sta)
    class(convertor), intent(inout) :: self
    integer, intent(in) :: i_sta
    double precision, allocatable :: tmp(:, :), processed(:,:)
    integer :: n2, n4, io, it, j, irow
    logical :: first_taper
    
    self%i_read_file = 1
    self%c3(i_sta) = c3_data(files=self%filenames(1, 1:3, i_sta))
    self%dt = self%c3(i_sta)%get_dt()
    self%n = nint(self%t_win / self%dt)
    self%c3_out(i_sta) = c3_data(dt=self%dt, n_cmps=self%n_cmps)
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

    
    tmp(1:self%n, 1:3) = 0.d0
    it = 0
    irow = 0
    
    write(*,*)"n_stored = ", self%c3(i_sta)%get_n_smp()
    write(*,*)self%filenames
    first_taper = .true.
    open(newunit=io, file="test.txt", status="replace")
    do 
       if (self%c3(i_sta)%get_n_smp() < self%n) then
          self%i_read_file = self%i_read_file + 1
          write(*,*)"i_read_file=", self%i_read_file
          if (self%i_read_file > size(self%filenames(:,1,i_sta))) exit
          write(*,*)"enqueued from <", &
               & trim(self%filenames(self%i_read_file,1,i_sta)), " / ",&
               & trim(self%filenames(self%i_read_file,2,i_sta)), " / ",&
               & trim(self%filenames(self%i_read_file,3,i_sta)), "> "
          call self%c3(i_sta)%enqueue_from_files(&
               & self%filenames(self%i_read_file,:,i_sta))
          
          write(*,*)"n_stored = ", self%c3(i_sta)%get_n_smp()
               
          
       end if
       tmp(1:n2,1:self%n_cmps) = tmp(n2+1:self%n,1:self%n_cmps)
       tmp(n2+1:self%n,1:self%n_cmps) = self%c3(i_sta)%dequeue_data(n2)

       write(*,*)"dequeued: n_sotred = ", self%c3(i_sta)%get_n_smp()

       
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
         end do
       end block
       
       ! output
       if (.not. first_taper) then
          call self%c3_out(i_sta)%enqueue_data(&
               & processed(n4+1:self%n-n4,1:self%n_cmps))
       else
          call self%c3_out(i_sta)%enqueue_data(&
               & processed(1:self%n-n4,1:self%n_cmps))
          first_taper = .false.
       end if
       
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
       
    end do
    close(io)
    
    return 
  end subroutine convertor_convert

  !---------------------------------------------------------------------

  function convertor_envelope(self, x) result(env)
    class(convertor), intent(inout) :: self
    double precision, intent(in) :: x(self%n)
    double precision :: env(self%n)



    self%r_tmp(1:self%n) = x(1:self%n)
    call fftw_execute_dft_r2c(self%plan_r2c, self%r_tmp, self%c_tmp)
    
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
