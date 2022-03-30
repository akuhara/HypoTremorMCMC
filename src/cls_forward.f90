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
     double precision, allocatable :: t_stdv(:,:), t_obs(:,:)
     double precision, allocatable :: dt_stdv(:,:,:), dt_obs(:,:,:)
     double precision, allocatable :: log_dt_stdv(:,:,:)     
     double precision, allocatable :: log_t_stdv(:,:)
     double precision, allocatable :: cov(:,:)
     double precision, allocatable :: inv_cov(:,:,:)
     double precision, allocatable :: log_det_cov(:)
     
     logical :: use_laplace
     logical :: forward_diff
     logical :: use_covariance

   contains
     procedure :: calc_log_likelihood => forward_calc_log_likelihood
     procedure :: calc_travel_time => forward_calc_travel_time
     procedure :: calc_dt => forward_calc_dt
     procedure :: init_cov => forward_init_cov
  end type forward

  interface forward
     module procedure init_forward
  end interface forward

contains

  !------------------------------------------------------------------------------

  type(forward) function init_forward(n_sta, n_events, sta_x, sta_y, sta_z, &
       & obs, use_laplace, use_covariance, forward_diff) result(self)
    integer, intent(in) :: n_sta, n_events
    double precision, intent(in) :: sta_x(:), sta_y(:), sta_z(:)
    type(obs_data), intent(in) :: obs
    logical, intent(in) :: use_laplace, forward_diff, use_covariance
    integer :: i, j, k

    self%n_sta = n_sta
    allocate(self%sta_x(self%n_sta))
    allocate(self%sta_y(self%n_sta))
    allocate(self%sta_z(self%n_sta))
    self%sta_x = sta_x
    self%sta_y = sta_y
    self%sta_z = sta_z
    self%use_laplace = use_laplace
    self%forward_diff = forward_diff
    self%use_covariance = use_covariance
    
    self%n_events = n_events
    !self%obs = obs

    allocate(self%t_obs(self%n_sta, self%n_events))
    allocate(self%t_stdv(self%n_sta, self%n_events))
    allocate(self%log_t_stdv(self%n_sta, self%n_events))
    
    
    if (self%forward_diff) then
       allocate(self%dt_obs(self%n_sta, self%n_sta, self%n_events))
       allocate(self%dt_stdv(self%n_sta, self%n_sta, self%n_events))
       allocate(self%log_dt_stdv(self%n_sta, self%n_sta, self%n_events))
    end if
    
    if (self%use_covariance) then
       allocate(self%cov(self%n_sta, self%n_sta))
       allocate(self%inv_cov(self%n_sta, self%n_sta, self%n_events))
       allocate(self%log_det_cov(self%n_events))
       if (self%forward_diff) then
          error stop "use_covariance and forward_diff cannot coexit"
       end if
    end if
    


    if (self%forward_diff) then
       self%dt_obs = obs%get_dt_obs()
       self%dt_stdv = obs%get_dt_stdv()
    else
       self%t_obs = obs%get_t_obs()
       self%t_stdv = obs%get_t_stdv()
    end if
       
    if (self%use_covariance) then
       call self%init_cov()
    end if


    !self%data_used = .false.
    if (self%forward_diff) then
       do concurrent (i=1:self%n_events)
          do concurrent (j=1:self%n_sta-1)
             do concurrent (k=j+1:self%n_sta)
                if (self%dt_stdv(k,j,i) > 1.d-16) then
                   self%log_dt_stdv(k, j, i) = log(self%dt_stdv(k, j, i))
                else
                   self%log_dt_stdv(k, j, i) =1.d0
                   self%dt_stdv(k, j, i) = 1.d0
                end if
             end do
          end do
       end do
    else
       do concurrent (i=1:self%n_events)
          do concurrent (j=1:self%n_sta)
             if (self%t_stdv(j,i) > 1.d-16) then
                self%log_t_stdv(j, i) = log(self%t_stdv(j, i))
             else
                self%log_t_stdv(j, i) =1.d0
                self%t_stdv(j, i) = 1.d0
             end if
          end do
       end do
    end if

    
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
    if (.not. self%forward_diff) then
       do concurrent (i = 1:self%n_events)
          t_mean = sum(t_syn(1:self%n_sta,i)) / self%n_sta
          t_syn(1:self%n_sta,i) = t_syn(1:self%n_sta,i) - t_mean
       end do
    end if
    
    
    return 
  end subroutine forward_calc_travel_time
    
  !------------------------------------------------------------------------------
  
  subroutine forward_calc_log_likelihood(self, hypo, t_corr, vs, log_likelihood)
    class(forward), intent(inout) :: self
    type(model), intent(in) :: hypo, t_corr, vs
    double precision, intent(out) :: log_likelihood
    double precision :: t_syn(self%n_sta, self%n_events)
    double precision :: dt_syn(self%n_sta, self%n_sta, self%n_events)
    double precision :: tmp(self%n_sta, self%n_sta, self%n_events), b
    integer :: i, j,  k
    call self%calc_travel_time(hypo, t_corr, vs, t_syn)
    
    if (self%forward_diff) then
       call self%calc_dt(t_syn, dt_syn)
    end if

    log_likelihood = 0.d0
    !return 

    
    if (self%forward_diff) then
       if (.not. self%use_laplace) then
          do i = 1, self%n_events
             do j = 1, self%n_sta - 1
                do k = j + 1, self%n_sta
                   if (self%dt_stdv(k,j,i) == 0.d0) cycle
                   log_likelihood = log_likelihood - &
                        & (self%dt_obs(k, j, i) - dt_syn(k, j, i))**2 / &
                        & (2.d0 * self%dt_stdv(k, j, i)**2) - log_2pi_half &
                        & -self%log_dt_stdv(k,j,i)
                   
                end do
             end do
          end do
       else
          do i = 1, self%n_events
             do j = 1, self%n_sta - 1
                do k = j + 1, self%n_sta
                   if (self%dt_stdv(k,j,i) == 0.d0) cycle
                   log_likelihood = log_likelihood - &
                        & abs(self%dt_obs(k, j, i) - dt_syn(k, j, i)) / &
                        & (self%dt_stdv(k,j,i)/sqrt(2.d0)) &
                        & - 0.5d0 * log(2.d0) - self%log_dt_stdv(k,j,i) 
                end do
             end do
          end do
       end if
    else
       if (.not. self%use_laplace) then
          do i = 1, self%n_events
             do j = 1, self%n_sta 
                if (self%t_stdv(j,i) == 0.d0) cycle
                !if (self%t_stdv(j,i) > 1.d0) cycle
                log_likelihood = log_likelihood - &
                     & (self%t_obs(j, i) - t_syn(j, i))**2 / &
                     & (2.d0 * self%t_stdv(j, i)**2) - log_2pi_half &
                     & -self%log_t_stdv(j,i)
             end do
          end do
       else
          do i = 1, self%n_events
             do j = 1, self%n_sta 
                
                log_likelihood = log_likelihood - &
                     abs(self%t_obs(j, i) - t_syn(j, i)) / &
                     & (self%t_stdv(j,i)/sqrt(2.d0)) &
                     & - 0.5d0 * log(2.d0) - self%log_t_stdv(j,i) 
             end do
          end do
       end if
    end if
    
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
 
  !------------------------------------------------------------------------------
 
  subroutine forward_init_cov(self)
    class(forward), intent(inout) :: self
    integer :: i, j, k
    double precision :: r, r2, s2, si2, sj2
    double precision :: tmp(self%n_sta, self%n_sta)
    double precision, allocatable :: work(:)
    integer :: ipiv(self%n_sta)
    integer :: info
    integer :: lwork
    
    r = 1.d0 / dble(self%n_sta)
    r2 = 1.d0 - r

    do k = 1, self%n_events
       s2 = sum(self%t_stdv(1:self%n_sta, k)**2)
       
       
       do i = 1, self%n_sta
          si2 = self%t_stdv(i,k) * self%t_stdv(i,k)
          self%cov(i,i) = &
               & r2 * r2 * si2 + (s2 - si2) * r * r
       end do
       
       do i = 1, self%n_sta-1
          si2 = self%t_stdv(i,k) * self%t_stdv(i,k)
          do j = i+1, self%n_sta
             sj2 = self%t_stdv(j,k) * self%t_stdv(j,k)
             self%cov(j,i) = - r * r2 * si2 - r * r2 * sj2
             self%cov(i,j) = self%cov(j,i)
          end do
       end do
       
       tmp = self%cov
       
       ! Cholesky decomposition
       call dgetrf(self%n_sta, self%n_sta, tmp, self%n_sta, ipiv, info)
       if (info /= 0) then
          error stop "ERROR in Cholesky decomposition"
       end if
       
       ! Determinant
       self%log_det_cov(k) = 0.d0
       do i = 1, self%n_sta
          self%log_det_cov(k) = self%log_det_cov(k) + log(abs(tmp(i,i)))
       end do
       print *, "det = ", self%log_det_cov(k), exp(self%log_det_cov(k))
       
       allocate(work(1))
       call dgetri(self%n_sta, tmp, self%n_sta, ipiv, work, -1, info)
       if (info /= 0) then
          error stop "ERROR in dgetri (lwork query)"
       end if
       lwork = int(real(work(1)) + 0.5d0)
       deallocate(work)
       allocate(work(lwork))
       call dgetri(self%n_sta, tmp, self%n_sta, ipiv, work, lwork, info)
       if (info /= 0) then
          error stop "ERROR in dgetri (2nd)"
       end if
       
       self%inv_cov(:,:,k) = tmp
       tmp = matmul(tmp, self%cov)

       do i = 1, 3
          print *, i, i, tmp(i,i)
       end do


          
       
    end do
    

    stop
    return 
  end subroutine forward_init_cov
 
end module cls_forward
