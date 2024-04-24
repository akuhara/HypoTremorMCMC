module cls_selector
  use cls_line_text, only: line_max
  use mod_regress, only: linear_regression, weighted_corr
  implicit none
  
  type selector
     private

     ! stations
     integer :: n_sta
     character(line_max), allocatable :: station_names(:)
     double precision, allocatable :: sta_x(:), sta_y(:), sta_z(:)
     
     ! Assumed focal depth
     double precision :: z_guess
     
     ! Optimized data
     double precision, allocatable :: t(:), t_err(:)
     double precision, allocatable :: a(:), a_err(:)

     ! Source-receiver distances
     double precision, allocatable :: d(:, :) ! d(receiver_index, source_index)
     
   contains
     procedure :: eval_wave_propagation => selector_eval_wave_propagation
     
  end type selector
  
  interface selector
     module procedure init_selector
  end interface selector
     
     
     
   contains
     
     !------------------------------------------------------------------------
     
     type(selector) function init_selector(station_names, &
          & sta_x, sta_y, sta_z, z_guess) result(self)
       character(line_max), intent(in) :: station_names(:)
       double precision, intent(in) :: sta_x(:), sta_y(:), sta_z(:)
       double precision, intent(in) :: z_guess
       integer :: i, j
       
       ! Set stations
       self%n_sta = size(station_names)
       allocate(self%station_names(self%n_sta))
       allocate(self%sta_x(self%n_sta))
       allocate(self%sta_y(self%n_sta))
       allocate(self%sta_z(self%n_sta))
       self%station_names = station_names
       self%sta_x = sta_x
       self%sta_y = sta_y
       self%sta_z = sta_z
       
       ! Event
       self%z_guess = z_guess

       ! Set soure-receiver distances
       allocate(self%d(self%n_sta, self%n_sta))
       do concurrent (i = 1:self%n_sta)
          do concurrent (j = 1:self%n_sta)
             self%d(j, i) = sqrt((sta_x(j) - sta_x(i))** 2 + &
                  &              (sta_y(j) - sta_y(i))** 2 + &
                  &              (sta_z(j) - z_guess)** 2 )
          end do
       end do
       
       return 
     end function init_selector
     
     !------------------------------------------------------------------------
     
     subroutine selector_eval_wave_propagation(self, id, vs, t0, b, a0, &
          & cc_t, cc_a)
       class(selector), intent(inout) :: self
       integer, intent(in) :: id
       double precision, intent(out) :: vs, b, t0, a0, cc_t, cc_a
       integer :: ios, io, i, io2
       integer :: i_sta_nearest(1)
       character(line_max) :: opt_file, out_file
       double precision :: t(self%n_sta), a(self%n_sta)
       double precision :: w_t(self%n_sta), w_a(self%n_sta)
       double precision :: t_err(self%n_sta), a_err(self%n_sta)
       double precision :: dummy1, dummy2, dummy3
       double precision :: slope, intercept, res
       
       write(opt_file, '(A,I6.6,A)')"opt_data.", id, ".dat"
       open(newunit=io, file=opt_file, iostat=ios, status='old')
       print *, "Now working on ", trim(opt_file)
       do i = 1, self%n_sta
          read(io, *) dummy1, dummy2, dummy3, t(i), t_err(i), a(i), a_err(i)
       end do
       close(io)

       ! Guess the nearest station from the amplitude
       i_sta_nearest = maxloc(a)
       
       ! geometrical spereading correction
       a(1:self%n_sta) = a(1:self%n_sta) + &
            & log(self%d(1:self%n_sta, i_sta_nearest(1)))
       

       write(out_file,'(A,i6.6,A)')"dist_plot.", id, ".dat"
       open(newunit=io2, status="replace", file=out_file, iostat=ios)
       do i = 1, self%n_sta
          write(io2, *)self%d(i,i_sta_nearest(1)), t(i), t_err(i), a(i), a_err(i)
       end do
       close(io2)

       ! Regress
       w_t(1:self%n_sta) = 1.d0 / t_err(1:self%n_sta)**2
       call linear_regression(self%n_sta, self%d(:, i_sta_nearest(1)), &
            & t,  w_t, slope, intercept, res)
       vs = 1.0 / slope
       t0 = intercept

       w_a(1:self%n_sta) = 1.d0 / a_err(1:self%n_sta)**2
       call linear_regression(self%n_sta, self%d(:, i_sta_nearest(1)), &
            & a,  w_a, slope, intercept, res)
       b = -1.0 * slope
       a0 = intercept

       ! Correlation coefficient
       call weighted_corr(self%n_sta, self%d(:, i_sta_nearest(1)), &
            & t, w_t, cc_t)
       call weighted_corr(self%n_sta, self%d(:, i_sta_nearest(1)), &
            & a, w_a, cc_a)

       return 
     end subroutine selector_eval_wave_propagation

     !------------------------------------------------------------------------

     

   end module cls_selector
