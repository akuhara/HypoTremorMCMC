module cls_detector
  use cls_c3_data, only: c3_data
  use cls_line_text, only: line_max
  implicit none
  
  type detector
     private
     
     integer :: n_stations
     type(c3_data), allocatable :: c3(:)
     character(line_max), allocatable :: station_names(:)

     integer :: n_histo = 9000
     double precision :: thred = 0.95d0

   contains
     procedure :: calc_stats => detector_calc_stats
  end type detector
  
  interface detector
     module procedure init_detector
  end interface detector
  
contains
  
  !-------------------------------------------------------------------------
  type(detector) function init_detector(c3, station_names) result(self)
    type(c3_data), intent(in) :: c3(:)
    character(line_max), intent(in) :: station_names(:)
    integer :: n_stations
    
    n_stations = size(c3)
    self%n_stations = n_stations
    allocate(self%c3(n_stations))
    self%c3(:) = c3(:)
    allocate(self%station_names(n_stations))
    self%station_names = station_names

    return 
  end function init_detector
  
  !-------------------------------------------------------------------------
  
  subroutine detector_calc_stats(self, i_sta, mean, stdv)
    class(detector), intent(inout) :: self
    integer, intent(in) :: i_sta
    double precision, allocatable, intent(out) :: mean(:), stdv(:)
    double precision, allocatable :: pcnt(:)
    integer :: n_cmps, n_smp, i_cmp, io, i, k, count
    integer, allocatable :: histo(:,:)
    character(line_max) :: out_file

    if (i_sta < 0 .or. i_sta > self%n_stations) then
       error stop "ERROR: invalid station ID is given to detector_calc_stats"
    end if
    
    n_cmps = self%c3(i_sta)%get_n_cmps()
    n_smp  = self%c3(i_sta)%get_n_smp()
    allocate(histo(self%n_histo, n_cmps))
    allocate(mean(n_cmps))
    allocate(stdv(n_cmps))
    
    mean = self%c3(i_sta)%calc_mean()

    ! histogram
    histo = self%c3(i_sta)%calc_histogram(self%n_histo)
    do i_cmp = 1, n_cmps
       print *, "mean :", i_cmp, mean(i_cmp)
       
    end do

    out_file = trim(self%station_names(i_sta)) // "_1.histo" 
    open(newunit=io, file=out_file, status='unknown')
    do i = 1, self%n_histo
       write(io, *)i, histo(i, 1)
    end do
    close(io)


    ! Rough estimation of percentile
    allocate(pcnt(n_cmps))
    k = int((1.d0 - self%thred) * n_smp)

    do i_cmp = 1, n_cmps
       count = 0
       do i = self%n_histo, 1, -1
          count = count + histo(i, i_cmp)
          if (count >= k) then
             pcnt(i_cmp) = i
             exit
          end if
       end do
       print *, "pcnt :",  trim(self%station_names(i_sta)), " ", &
            & i_cmp, pcnt(i_cmp)
    end do


    return 
  end subroutine detector_calc_stats

  !-------------------------------------------------------------------------
  
end module cls_detector
