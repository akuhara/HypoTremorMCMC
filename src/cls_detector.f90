module cls_detector
  use cls_c3_data
  implicit none
  
  type detector
     private
     
     integer :: n_stations
     type(c3_data), allocatable :: c3(:)

   contains
     procedure :: calc_stats => detector_calc_stats
  end type detector
  
  interface detector
     module procedure init_detector
  end interface detector
  
contains
  
  !-------------------------------------------------------------------------
  type(detector) function init_detector(c3) result(self)
    type(c3_data), intent(in) :: c3(:)
    integer :: n_stations
    
    n_stations = size(c3)
    self%n_stations = n_stations
    allocate(self%c3(n_stations))
    self%c3(:) = c3(:)
    
    return 
  end function init_detector
  
  !-------------------------------------------------------------------------
  
  subroutine detector_calc_stats(self, i_sta, mean, stdv)
    class(detector), intent(inout) :: self
    integer, intent(in) :: i_sta
    double precision, allocatable, intent(out) :: mean(:), stdv(:)
    integer :: n_cmps, n_smp, i_cmp

    if (i_sta < 0 .or. i_sta > self%n_stations) then
       error stop "ERROR: invalid station ID is given to detector_calc_stats"
    end if
    
    n_cmps = self%c3(i_sta)%get_n_cmps()
    n_smp  = self%c3(i_sta)%get_n_smp()

    allocate(mean(n_cmps))
    allocate(stdv(n_cmps))
    
    mean = self%c3(i_sta)%calc_mean()
    do i_cmp = 1, n_cmps
       print *, "mean :", i_cmp, mean(i_cmp)

    end do

    return 
  end subroutine detector_calc_stats

  !-------------------------------------------------------------------------
  
end module cls_detector
