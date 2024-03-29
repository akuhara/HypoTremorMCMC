module cls_forward
  use cls_obs_data, only: obs_data
  use cls_model, only: model
  implicit none
  double precision, parameter, private :: log_2pi_half = 0.5d0 * log(2.d0 * acos(-1.d0))
  type forward
     private
     integer :: n_events
     integer :: n_sta
     integer :: n_tri
     double precision, allocatable :: sta_x(:), sta_y(:), sta_z(:)
     double precision, allocatable :: t_stdv(:,:), t_obs(:,:), t_precision(:,:)
     double precision, allocatable :: a_stdv(:,:), a_obs(:,:), a_precision(:,:)
     double precision, allocatable :: log_dt_stdv(:,:,:)     
     double precision, allocatable :: log_t_stdv(:,:)
     double precision, allocatable :: log_a_stdv(:,:)
     
     logical :: use_amp
     logical :: use_time

     
   contains
     procedure :: calc_log_likelihood => forward_calc_log_likelihood
     procedure :: partially_update_log_likelihood => &
          & forward_partially_update_log_likelihood
     procedure :: calc_travel_time => forward_calc_travel_time
     procedure :: calc_amp => forward_calc_amp
     procedure :: calc_travel_time_single => forward_calc_travel_time_single
     procedure :: calc_amp_single => forward_calc_amp_single
  end type forward

  interface forward
     module procedure init_forward
  end interface forward

contains

  !------------------------------------------------------------------------------

  type(forward) function init_forward(n_sta, n_events, sta_x, sta_y, sta_z, &
       & obs, use_amp, use_time) result(self)
    integer, intent(in) :: n_sta, n_events
    double precision, intent(in) :: sta_x(:), sta_y(:), sta_z(:)
    type(obs_data), intent(in) :: obs
    logical, intent(in) :: use_amp, use_time
    
    integer :: i, j, k

    self%n_sta = n_sta
    allocate(self%sta_x(self%n_sta))
    allocate(self%sta_y(self%n_sta))
    allocate(self%sta_z(self%n_sta))
    self%sta_x = sta_x
    self%sta_y = sta_y
    self%sta_z = sta_z
    
    self%n_events = n_events
    
    self%use_amp = use_amp
    self%use_time = use_time
    
    allocate(self%t_obs(self%n_sta, self%n_events))
    allocate(self%t_stdv(self%n_sta, self%n_events))
    allocate(self%a_obs(self%n_sta, self%n_events))
    allocate(self%a_stdv(self%n_sta, self%n_events))
    allocate(self%log_t_stdv(self%n_sta, self%n_events))
    allocate(self%log_a_stdv(self%n_sta, self%n_events))
    allocate(self%t_precision(self%n_sta, self%n_events))
    allocate(self%a_precision(self%n_sta, self%n_events))
    
    self%t_obs = obs%get_t_obs()
    self%t_stdv = obs%get_t_stdv()
    self%a_obs  = obs%get_a_obs()
    self%a_stdv = obs%get_a_stdv()
    
    do concurrent (i=1:self%n_events)
       do concurrent (j=1:self%n_sta)
          if (self%t_stdv(j,i) > 1.d-16) then
             self%log_t_stdv(j, i) = log(self%t_stdv(j, i))
             self%t_precision(j, i) = 1.d0 / self%t_stdv(j, i)**2
             self%log_a_stdv(j, i) = log(self%a_stdv(j, i))
             self%a_precision(j, i) = 1.d0 / self%a_stdv(j, i)**2
          else
             self%log_t_stdv(j, i) =1.d0
             self%t_stdv(j, i) = 1.d0
             self%t_precision(j, i) = 1.d0
             self%log_a_stdv(j, i) =1.d0
             self%a_stdv(j, i) = 1.d0
             self%a_precision(j, i) = 1.d0
          end if
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
    double precision :: x, y, z, beta, tc, t_mean
    
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
    
    ! Demean
    do concurrent (i = 1:self%n_events)
       t_mean = &
            & sum(                                    &
            &      self%t_precision(1:self%n_sta,i) * &
            &      (t_syn(1:self%n_sta,i) -      &
            &       self%t_obs(1:self%n_sta,i))            &
            &     )                                   &
            & / sum(self%t_precision(1:self%n_sta,i))
       t_syn(1:self%n_sta,i) = t_syn(1:self%n_sta,i) - t_mean
    end do
    
    
    
    return 
  end subroutine forward_calc_travel_time

  !------------------------------------------------------------------------------

  subroutine forward_calc_travel_time_single(self, evt_id, hypo, t_corr, vs, t_syn)
    class(forward), intent(inout) :: self
    integer, intent(in) :: evt_id
    type(model), intent(in) :: hypo, t_corr, vs
    double precision, intent(out) :: t_syn(self%n_sta)
    integer :: j
    double precision :: x, y, z, beta, tc, t_mean
    
    beta = vs%get_x(1)
    x = hypo%get_x(3*evt_id-2)
    y = hypo%get_x(3*evt_id-1)
    z = hypo%get_x(3*evt_id)
    
    do j = 1, self%n_sta
       tc = t_corr%get_x(j)
       t_syn(j) = sqrt((x - self%sta_x(j))**2 &
            &          + (y - self%sta_y(j))**2 &
            &          + (z - self%sta_z(j))**2) &
            & / beta - tc
       
    end do
    
    
    ! Demean
    t_mean = &
         & sum(                                    &
         &      self%t_precision(1:self%n_sta, evt_id) * &
         &      (t_syn(1:self%n_sta) -      &
         &       self%t_obs(1:self%n_sta, evt_id))            &
         &     )                                   &
         & / sum(self%t_precision(1:self%n_sta, evt_id))
    t_syn(1:self%n_sta) = t_syn(1:self%n_sta) - t_mean
    
    
    
    
    return 
  end subroutine forward_calc_travel_time_single
    
  !------------------------------------------------------------------------------

  subroutine forward_calc_amp(self, hypo, a_corr, qs, vs, a_syn)
    class(forward), intent(inout) :: self
    type(model), intent(in) :: hypo, a_corr, qs, vs
    double precision, intent(out) :: a_syn(self%n_sta, self%n_events)
    integer :: i, j
    double precision :: x, y, z, q, ac, a_mean, beta, d
    double precision, parameter :: pi = acos(-1.d0)
    double precision, parameter :: freq = 5.d0
    
    q = qs%get_x(1)
    beta = vs%get_x(1)
    do i = 1, self%n_events
       x = hypo%get_x(3*i-2)
       y = hypo%get_x(3*i-1)
       z = hypo%get_x(3*i)
       
       do j = 1, self%n_sta
          ac = a_corr%get_x(j)
          d = sqrt((x - self%sta_x(j))**2 &
               & + (y - self%sta_y(j))**2 &
               & + (z - self%sta_z(j))**2)
          a_syn(j,i) = - d * pi * freq / ( q * beta) - log(d) - ac
       end do
    end do
    
    ! Demean
    do concurrent (i = 1:self%n_events)
       a_mean = &
            & sum(                                    &
               &      self%a_precision(1:self%n_sta,i) * &
               &      (a_syn(1:self%n_sta,i) -      &
               &       self%a_obs(1:self%n_sta,i))            &
               &     )                                   &
               & / sum(self%a_precision(1:self%n_sta,i))
       a_syn(1:self%n_sta,i) = a_syn(1:self%n_sta,i) - a_mean
    end do
    
    
    return 
  end subroutine forward_calc_amp
    
  !------------------------------------------------------------------------------

  subroutine forward_calc_amp_single(self, evt_id, hypo, a_corr, qs, vs, a_syn)
    class(forward), intent(inout) :: self
    integer, intent(in) :: evt_id
    type(model), intent(in) :: hypo, a_corr, qs, vs
    double precision, intent(out) :: a_syn(self%n_sta)
    integer :: j
    double precision :: x, y, z, q, ac, a_mean, beta, d
    double precision, parameter :: pi = acos(-1.d0)
    double precision, parameter :: freq = 5.d0
    
    q = qs%get_x(1)
    beta = vs%get_x(1)
    x = hypo%get_x(3*evt_id-2)
    y = hypo%get_x(3*evt_id-1)
    z = hypo%get_x(3*evt_id)
    
    do j = 1, self%n_sta
       ac = a_corr%get_x(j)
       d = sqrt((x - self%sta_x(j))**2 &
            & + (y - self%sta_y(j))**2 &
            & + (z - self%sta_z(j))**2)
       a_syn(j) = - d * pi * freq / ( q * beta) - log(d) - ac
    end do
    
    
    ! Demean
    a_mean = &
         & sum(                                    &
         &      self%a_precision(1:self%n_sta, evt_id) * &
         &      (a_syn(1:self%n_sta) -      &
         &       self%a_obs(1:self%n_sta, evt_id))            &
         &     )                                   &
         & / sum(self%a_precision(1:self%n_sta, evt_id))
    a_syn(1:self%n_sta) = a_syn(1:self%n_sta) - a_mean

    
    
    return 
  end subroutine forward_calc_amp_single
    
  !------------------------------------------------------------------------------
  
  subroutine forward_calc_log_likelihood(self, hypo, t_corr, vs, a_corr, qs, &
       & log_likelihood)
    class(forward), intent(inout) :: self
    type(model), intent(in) :: hypo, t_corr, vs, a_corr, qs
    double precision, intent(out) :: log_likelihood
    double precision :: t_syn(self%n_sta, self%n_events)
    double precision :: a_syn(self%n_sta, self%n_events)
    integer :: i, j
    
    log_likelihood = 0.d0

    if (self%use_time) then
       call self%calc_travel_time(hypo, t_corr, vs, t_syn)
       do i = 1, self%n_events
          do j = 1, self%n_sta 
             log_likelihood = log_likelihood - &
                  & (self%t_obs(j, i) - t_syn(j, i))**2 / &
                  & (2.d0 * self%t_stdv(j, i)**2) - log_2pi_half &
                  & -self%log_t_stdv(j,i)
          end do
       end do
    end if
    if (self%use_amp) then
       call self%calc_amp(hypo, a_corr, qs, vs, a_syn)
       do i = 1, self%n_events
          do j = 1, self%n_sta 
             log_likelihood = log_likelihood - &
                  & (self%a_obs(j, i) - a_syn(j, i))**2 / &
                  & (2.d0 * self%a_stdv(j, i)**2) - log_2pi_half &
                  & -self%log_a_stdv(j,i)
          end do
       end do
    end if
    
    return 
  end subroutine forward_calc_log_likelihood

  !------------------------------------------------------------------------------

  subroutine forward_partially_update_log_likelihood(self, evt_id, hypo_old, &
       & log_likelihood_old, hypo, t_corr, vs, a_corr, qs, &
       & log_likelihood)
    class(forward), intent(inout) :: self
    integer, intent(in) :: evt_id
    type(model), intent(in) :: hypo_old, hypo, t_corr, vs
    type(model), intent(in) :: a_corr, qs
    double precision, intent(in) :: log_likelihood_old
    double precision, intent(out) :: log_likelihood
    double precision :: t_syn(self%n_sta), a_syn(self%n_sta)
    integer :: j

    log_likelihood = log_likelihood_old
    if (self%use_time) then
       call self%calc_travel_time_single(evt_id, &
            & hypo_old, t_corr, vs, t_syn)
       do j = 1, self%n_sta 
          log_likelihood = log_likelihood + &
               & (self%t_obs(j, evt_id) - t_syn(j))**2 / &
               & (2.d0 * self%t_stdv(j, evt_id)**2) + log_2pi_half &
               &  + self%log_t_stdv(j,evt_id)
       end do
       
       call self%calc_travel_time_single(evt_id, &
            & hypo, t_corr, vs, t_syn)
       do j = 1, self%n_sta 
          log_likelihood = log_likelihood - &
               & (self%t_obs(j, evt_id) - t_syn(j))**2 / &
               & (2.d0 * self%t_stdv(j, evt_id)**2) - log_2pi_half &
               &  - self%log_t_stdv(j,evt_id)
       end do
    end if

    if (self%use_amp) then
       call self%calc_amp_single(evt_id, &
            & hypo_old, a_corr, qs, vs, a_syn)
       do j = 1, self%n_sta 
          log_likelihood = log_likelihood + &
               & (self%a_obs(j, evt_id) - a_syn(j))**2 / &
               & (2.d0 * self%a_stdv(j, evt_id)**2) + log_2pi_half &
               &  + self%log_a_stdv(j,evt_id)
       end do
       
       call self%calc_amp_single(evt_id, &
            & hypo, a_corr, qs, vs, a_syn)
       do j = 1, self%n_sta 
          log_likelihood = log_likelihood - &
               & (self%a_obs(j, evt_id) - a_syn(j))**2 / &
               & (2.d0 * self%a_stdv(j, evt_id)**2) - log_2pi_half &
               &  - self%log_a_stdv(j,evt_id)
       end do
    end if
    

    return 
  end subroutine forward_partially_update_log_likelihood
    
  
  !------------------------------------------------------------------------------
 
end module cls_forward
