module cls_measurer
  use cls_line_text, only: line_max
  implicit none
  
  type measurer
     private

     integer :: n_sta
     integer :: n_pair
     character(line_max), allocatable :: pair(:,:)
     logical :: verb
     
   contains
     procedure :: check_files => measurer_check_files
     
  end type measurer

  
  interface measurer
     module procedure init_measurer
  end interface measurer

contains

  !-------------------------------------------------------------------
  
  type(measurer) function init_measurer(station_names, verb) result(self)
    character(line_max), intent(in) :: station_names(:)
    logical, intent(in) :: verb
    integer :: i, j, k
    self%n_sta = size(station_names)
    self%n_pair = self%n_sta * (self%n_sta - 1) / 2
    self%verb = verb
    
    allocate(self%pair(2, self%n_pair))
    k = 0
    do i = 1, self%n_sta - 1
       do j = i + 1, self%n_sta
          k = k + 1
          self%pair(1, k) = station_names(i)
          self%pair(2, k) = station_names(j)
       end do
    end do

    call self%check_files()
    
    return 
  end function init_measurer
  
  !-------------------------------------------------------------------
  
  subroutine measurer_check_files(self)
    class(measurer), intent(inout) :: self
    integer :: i, ierr
    character(line_max) :: target_file, sta1, sta2
    logical :: is_ok
    
    do i = 1, self%n_pair
       sta1 = self%pair(1, i)
       sta2 = self%pair(2, i)
       write(target_file, '(a)') &
            & trim(sta1) // "." // trim(sta2) // ".max_corr"
       if (self%verb) then
          print *, target_file
       end if
       inquire(file=target_file, EXIST=is_ok)
       if (.not. is_ok) then
          write(0,*) "ERROR: " // trim(target_file) // "does not exsit"
          error stop 
       end if
       
    end do
    
    
    
    return 
  end subroutine measurer_check_files
  
  !-------------------------------------------------------------------
  
end module cls_measurer
