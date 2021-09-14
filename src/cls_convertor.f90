module cls_convertor
  use cls_c3_data
  use cls_param
  use cls_line_text
  implicit none 
  
  type convertor
     private
     type(param) :: para
     type(c3_data), allocatable :: c3(:)
     
   contains
     procedure :: feed_data => convertor_feed_data
  end type convertor
  
  interface convertor
     module procedure init_convertor
  end interface convertor
     
contains
  
  !---------------------------------------------------------------------
  
  type(convertor) function init_convertor(para) result(self)
    type(param), intent(in) :: para
    
    self%para = para
    call self%feed_data()


    block
      integer :: i_sta
      character(line_max), allocatable :: stations(:)
      character(line_max) :: file
      
      allocate(stations(self%para%get_n_stations()))
      stations = self%para%get_stations()
      
      do i_sta = 1, para%get_n_stations()
         file = trim(stations(i_sta)) // '.txt'
         call self%c3(i_sta)%output_data(file)
      end do
    end block

    return 
  end function init_convertor

  !---------------------------------------------------------------------
  
  subroutine convertor_feed_data(self)
    class(convertor), intent(inout) :: self
    
    allocate(self%c3(self%para%get_n_stations()))
    block
      integer :: i_sta, i_id
      do i_sta = 1, self%para%get_n_stations()
         self%c3(i_sta) = c3_data([self%para%get_filename(1, 1, i_sta), &
              &          self%para%get_filename(1, 2, i_sta), &
              &          self%para%get_filename(1, 3, i_sta)])
         do i_id = 2, self%para%get_n_data_id()
            call self%c3(i_sta)%append_data(&
                 &    [self%para%get_filename(i_id, 1, i_sta), &
                 &     self%para%get_filename(i_id, 2, i_sta), &
                 &     self%para%get_filename(i_id, 3, i_sta)])
         end do
      end do
    end block
    
    return 
  end subroutine convertor_feed_data

end module cls_convertor
