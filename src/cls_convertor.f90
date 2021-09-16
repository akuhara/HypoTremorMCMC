module cls_convertor
  use cls_c3_data
  use cls_line_text
  use, intrinsic :: iso_c_binding
  implicit none 
  include 'fftw3.f03'

  type convertor
     private
     type(c3_data), allocatable :: c3(:)
     double precision           :: t_win
     double precision           :: t_window_interval
     double precision           :: dt
     integer                    :: n
     integer                    :: n_sta
     integer                    :: n_cmps
     character(line_max), allocatable :: filenames(:,:,:)
     
     ! Work space for FFTW
     type(C_PTR)                            :: ifftw
     real(C_DOUBLE), allocatable            :: r_tmp(:)
     complex(C_DOUBLE_COMPLEX), allocatable :: c_tmp(:)
     
     !
     integer :: i_read_file

   contains
     procedure :: calc_envelope => convertor_calc_envelope

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
    allocate(self%filenames(n_files, self%n_cmps, self%n_sta))
    self%filenames = filenames
    
    return 
  end function init_convertor
  
  !---------------------------------------------------------------------

  subroutine convertor_calc_envelope(self, i_sta)
    class(convertor), intent(inout) :: self
    integer, intent(in) :: i_sta
    double precision :: tmp(self%n, 3)
    integer :: n2, io, it, j, irow
    
    self%i_read_file = 1
    self%c3(i_sta) = c3_data(self%filenames(1, 1:3, i_sta))
    self%dt = self%c3(i_sta)%get_dt()
    self%n = nint(self%t_win / self%dt)
    n2 = self%n / 2
    if (n2 + n2 /= self%n) then
       error stop "ERROR: n2 + n2 /= self%n"
    end if
    
    self%ifftw = fftw_plan_dft_r2c_1d(self%n, self%r_tmp, self%c_tmp, &
         & FFTW_ESTIMATE)
    
    tmp(1:self%n, 1:3) = 0.d0
    it = 0
    irow = 0
    
    write(*,*)"n_stored = ", self%c3(i_sta)%get_n_smp()
    write(*,*)self%filenames
    
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
          call self%c3(i_sta)%enqueue_data(&
               & self%filenames(self%i_read_file,:,i_sta))
          
          write(*,*)"n_stored = ", self%c3(i_sta)%get_n_smp()
               
          
       end if
       tmp(1:n2,1:self%n_cmps) = tmp(n2+1:self%n,1:self%n_cmps)
       tmp(n2+1:self%n,1:self%n_cmps) = self%c3(i_sta)%dequeue_data(n2)
       it = it + n2
       write(*,*)"dequeued: n_sotred = ", self%c3(i_sta)%get_n_smp()
       do j = 1, self%n
          write(io,*)(it + j) + self%dt, tmp(j,1) / 30000.d0 + irow
       end do
       irow = irow + 1
       write(io,*)
    end do
    close(io)
    
    return 
  end subroutine convertor_calc_envelope

  !---------------------------------------------------------------------

end module cls_convertor
