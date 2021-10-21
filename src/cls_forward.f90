module cls_forward
  use cls_obs_data, only: obs_data
  use cls_model, only: model
  implicit none
  double precision, parameter, private :: log_2pi_half = 0.5d0 * log(2.d0 * acos(-1.d0))
  type forward
     private
     integer :: n_events
     integer :: n_sta
     double precision, allocatable :: sta_x(:), sta_y(:), sta_z(:)
     double precision, allocatable :: dt_stdv(:,:,:), dt_obs(:,:,:)
     double precision, allocatable :: log_dt_stdv(:,:,:)
     logical, allocatable :: data_used(:,:,:)

     

   contains
     procedure :: calc_log_likelihood => forward_calc_log_likelihood
     procedure :: calc_travel_time => forward_calc_travel_time
     procedure :: calc_dt => forward_calc_dt
  end type forward

  interface forward
     module procedure init_forward
  end interface forward

contains

  !------------------------------------------------------------------------------

  type(forward) function init_forward(n_sta, n_events, sta_x, sta_y, sta_z, &
       & obs) result(self)
    integer, intent(in) :: n_sta, n_events
    double precision, intent(in) :: sta_x(:), sta_y(:), sta_z(:)
    type(obs_data), intent(in) :: obs
    integer :: i, j, k

    self%n_sta = n_sta
    allocate(self%sta_x(self%n_sta))
    allocate(self%sta_y(self%n_sta))
    allocate(self%sta_z(self%n_sta))
    self%sta_x = sta_x
    self%sta_y = sta_y
    self%sta_z = sta_z

    self%n_events = n_events
    !self%obs = obs

    allocate(self%dt_obs(self%n_sta, self%n_sta, self%n_events))
    allocate(self%dt_stdv(self%n_sta, self%n_sta, self%n_events))
    allocate(self%log_dt_stdv(self%n_sta, self%n_sta, self%n_events))
    allocate(self%data_used(self%n_sta, self%n_sta, self%n_events))
    
    self%dt_obs = obs%get_dt_obs()
    self%dt_stdv = obs%get_dt_stdv()

    !self%data_used = .false.
    do concurrent (i=1:self%n_events)
       do concurrent (j=1:self%n_sta-1)
          do concurrent (k=j+1:self%n_sta)
             if (self%dt_stdv(k,j,i) > 1.d-16) then
                self%log_dt_stdv(k, j, i) = log(self%dt_stdv(k, j, i))
                !self%data_used(k, j, i) = .true.
             else
                self%log_dt_stdv(k, j, i) =1.d0
                self%dt_stdv(k, j, i) = 1.d0
                !self%data_used(k, j, i) = .false.
             end if
    
          end do
       end do
    end do


    
    return 
  end function init_forward

  !------------------------------------------------------------------------------

  subroutine forward_calc_travel_time(self, hypo, t_corr, vs, t_syn)
    class(forward), intent(inout) :: self
    type(model), intent(in) :: hypo, t_corr, vs
    double precision, intent(out) :: t_syn(self%n_sta, self%n_events)
    integer :: i, j
    double precision :: x, y, z, beta, tc
    
    beta = vs%get_x(1)
    do i = 1, self%n_events
       x = hypo%get_x(3*i-2)
       y = hypo%get_x(3*i-1)
       z = hypo%get_x(3*i)
       
       do j = 1, self%n_sta
          tc = t_corr%get_x(j)
          t_syn(j,i) = sqrt((x - self%sta_x(j))**2 &
               &          + (y - self%sta_y(j))**2 &
               &          + (z - self%sta_z(j))**2) &
               & / beta - tc
       end do
    end do
    
    return 
  end subroutine forward_calc_travel_time
    
  !------------------------------------------------------------------------------
  
  subroutine forward_calc_log_likelihood(self, hypo, t_corr, vs, log_likelihood)
    class(forward), intent(inout) :: self
    type(model), intent(in) :: hypo, t_corr, vs
    double precision, intent(out) :: log_likelihood
    double precision :: t_syn(self%n_sta, self%n_events)
    double precision :: dt_syn(self%n_sta, self%n_sta, self%n_events)
    double precision :: tmp(self%n_sta, self%n_sta, self%n_events)
    integer :: i, j,  k
    call self%calc_travel_time(hypo, t_corr, vs, t_syn)
    call self%calc_dt(t_syn, dt_syn)

    log_likelihood = 0.d0
    !return 

    !tmp(:,:,:) = -(self%dt_obs(:,:,:) - dt_syn(:,:,:)) ** 2 / &
    !     & (2.d0 * self%dt_stdv(:,:,:)) &
    !     & - log_2pi_half &
    !     & - self%log_dt_stdv(:,:,:)
    !log_likelihood = sum(tmp, mask=self%data_used)

    
    do i = 1, self%n_events
       do j = 1, self%n_sta - 1
          do k = j + 1, self%n_sta
             if (self%dt_stdv(k,j,i) == 0.d0) cycle
             if (self%dt_stdv(k,j,i) > 3.d0) cycle
             log_likelihood = log_likelihood - &
                  & (self%dt_obs(k, j, i) - dt_syn(k, j, i))**2 / &
                  & (2.d0 * self%dt_stdv(k, j, i)**2) - log_2pi_half &
            & -self%log_dt_stdv(k,j,i)
            
         end do
       end do
    end do

    return 
  end subroutine forward_calc_log_likelihood

  !------------------------------------------------------------------------------
  
  subroutine forward_calc_dt(self, t_syn, dt_syn)
    class(forward), intent(inout) :: self
    double precision, intent(in) :: t_syn(self%n_sta, self%n_events)
    double precision, intent(out) :: dt_syn(self%n_sta, self%n_sta, self%n_events)
    integer :: i, j, k

    do i = 1, self%n_events
       do j = 1, self%n_sta - 1
          do k = j + 1, self%n_sta
             dt_syn(k, j, i) = t_syn(j,i) - t_syn(k,i)
          end do
       end do
    end do
    
    return 
  end subroutine forward_calc_dt
  
end module cls_forward
