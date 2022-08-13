module cls_mcmc
  use mod_random
  use cls_model
  
  implicit none 
  
  type mcmc 
     private
     integer :: n_events, n_sta
     type(model) :: hypo, t_corr, vs, a_corr, qs
     integer :: n_proposal_type = 5
     integer :: n_iter

     integer :: i_iter
     integer, allocatable :: n_propose(:)
     integer, allocatable :: n_accept(:)
     integer :: i_proposal_type
     integer :: target_id
     double precision :: log_likelihood
     double precision, allocatable :: likelihood_saved(:)
     double precision, allocatable :: temp_saved(:)
     double precision :: temp = 1.d0
     logical :: is_accepted

     double precision :: p_hypo, p_vs, p_t_corr, p_qs, p_a_corr
     
   contains
     procedure :: propose_model => mcmc_propose_model
     procedure :: judge_model   => mcmc_judge_model
     procedure :: one_step_summary => mcmc_one_step_summary
     procedure :: set_temp => mcmc_set_temp
     procedure :: get_temp => mcmc_get_temp
     procedure :: get_log_likelihood => mcmc_get_log_likelihood
     procedure :: get_likelihood_saved => mcmc_get_likelihood_saved
     procedure :: get_temp_saved => mcmc_get_temp_saved
     procedure :: get_n_propose => mcmc_get_n_propose
     procedure :: get_n_accept => mcmc_get_n_accept
     procedure :: get_n_iter => mcmc_get_n_iter
     procedure :: get_hypo => mcmc_get_hypo
     procedure :: get_vs => mcmc_get_vs
     procedure :: get_qs => mcmc_get_qs
     procedure :: get_t_corr => mcmc_get_t_corr
     procedure :: get_a_corr => mcmc_get_a_corr
     procedure :: write_out_hypo => mcmc_write_out_hypo
     procedure :: write_out_t_corr => mcmc_write_out_t_corr
     procedure :: write_out_vs => mcmc_write_out_vs
     procedure :: write_out_a_corr => mcmc_write_out_a_corr
     procedure :: write_out_qs => mcmc_write_out_qs
  end type mcmc
  
  interface mcmc
     module procedure init_mcmc
  end interface mcmc


contains

  !---------------------------------------------------------------------
  
  type(mcmc) function init_mcmc(hypo, t_corr, vs, a_corr, qs, n_iter, &
       & solve_t_corr, solve_vs, solve_a_corr, solve_qs) result(self)
    type(model), intent(in) :: hypo, t_corr, vs
    type(model), intent(in) :: qs, a_corr
    integer, intent(in) :: n_iter
    logical, intent(in) :: solve_t_corr, solve_vs
    logical, intent(in) :: solve_a_corr, solve_qs

    self%hypo     = hypo
    self%n_events = hypo%get_nx() / 3
    self%t_corr   = t_corr
    self%n_sta    = t_corr%get_nx()
    self%a_corr   = a_corr
    self%vs       = vs
    self%qs       = qs
    
    allocate(self%n_propose(self%n_proposal_type))
    allocate(self%n_accept(self%n_proposal_type))
    self%n_propose = 0
    self%n_accept  = 0
    self%n_iter    = n_iter
    self%i_iter    = 0
    self%log_likelihood = -9.d+300

    ! set proposal probability
    self%p_vs = 0.d0
    self%p_t_corr = 0.d0
    self%p_qs = 0.d0
    self%p_a_corr = 0.d0
    if (solve_vs) then
       self%p_vs = 0.025d0
    end if
    if (solve_t_corr) then
       self%p_t_corr = 0.025d0
    end if
    if (solve_qs) then
       self%p_qs = 0.025d0
    end if
    if (solve_a_corr) then
       self%p_a_corr = 0.025d0
    end if
    self%p_hypo = 1.d0 - self%p_vs - self%p_t_corr &
         & - self%p_qs - self%p_a_corr
    
    return 
  end function init_mcmc

  !---------------------------------------------------------------------
  
  subroutine mcmc_propose_model(self, hypo_proposed, t_corr_proposed, &
       & vs_proposed, a_corr_proposed, qs_proposed, &
       & log_prior_ratio, prior_ok, prop_type, target_id)
    class(mcmc), intent(inout) :: self
    type(model), intent(out) :: hypo_proposed, t_corr_proposed, vs_proposed
    type(model), intent(out) :: a_corr_proposed, qs_proposed
    double precision, intent(out) :: log_prior_ratio
    integer, intent(out) :: prop_type, target_id
    integer :: itype, id, icmp
    logical, intent(out) :: prior_ok
    double precision :: a_select

    hypo_proposed   = self%hypo
    t_corr_proposed = self%t_corr
    vs_proposed     = self%vs
    a_corr_proposed = self%a_corr
    qs_proposed     = self%qs
    
    a_select = rand_u()

    prior_ok = .true. 
    
    if (a_select < self%p_vs) then
       ! Perturb Vs
       call vs_proposed%perturb(1, log_prior_ratio, prior_ok)
       self%i_proposal_type = 1
       target_id = 0
    else if (a_select < self%p_vs + self%p_t_corr) then
       ! Perturb t_corr
       id   = int(rand_u() * self%n_sta) + 1
       call t_corr_proposed%perturb(id, log_prior_ratio, prior_ok)
       self%i_proposal_type = 2
       target_id = id
    else if (a_select < self%p_vs + self%p_t_corr + self%p_qs) then
       ! Perturb Qs
       call qs_proposed%perturb(1, log_prior_ratio, prior_ok)
       self%i_proposal_type = 3
       target_id = 0
    else if (a_select < self%p_vs + self%p_t_corr + self%p_qs &
         & + self%p_a_corr) then
       ! Perturb a_corr
       id   = int(rand_u() * self%n_sta) + 1
       call a_corr_proposed%perturb(id, log_prior_ratio, prior_ok)
       self%i_proposal_type = 4
       target_id = id
    else 
       ! Perturb hypocenter
       id   = int(rand_u() * self%n_events) + 1
       icmp = int(rand_u() * 3) 
       call hypo_proposed%perturb(3*id-icmp, log_prior_ratio, prior_ok)
       self%i_proposal_type = 5
       target_id = 3 * id - icmp
    end if
    prop_type = self%i_proposal_type
    
    
    ! Count proposal
    self%n_propose(self%i_proposal_type) = &
         & self%n_propose(self%i_proposal_type) + 1
    
    return
  end subroutine mcmc_propose_model

  !---------------------------------------------------------------------

  subroutine mcmc_judge_model(self, hypo, t_corr, vs, a_corr, qs, &
       & log_likelihood, log_prior_ratio, prior_ok)
    class(mcmc), intent(inout) :: self
    type(model), intent(in) :: hypo, vs, t_corr, qs, a_corr
    double precision, intent(in) :: log_likelihood, log_prior_ratio
    logical, intent(in) :: prior_ok
    double precision :: ratio
    double precision :: r
    double precision, parameter :: eps = epsilon(1.d0)
    

    self%is_accepted = .false.
    if (prior_ok) then
       ratio = (log_likelihood - self%log_likelihood) / self%temp
       ratio = ratio + log_prior_ratio 
       
       r = rand_u()
       if (r >= eps) then
          if (log(r) <= ratio) then
             self%is_accepted = .true.
          end if
       end if
    end if


    if (self%is_accepted) then
       ! Accept model
       self%hypo           = hypo
       self%t_corr         = t_corr
       self%vs             = vs
       self%a_corr         = a_corr
       self%qs             = qs
       self%log_likelihood = log_likelihood
       self%n_accept(self%i_proposal_type) = &
         & self%n_accept(self%i_proposal_type) + 1
    end if

    ! Adds iteration counter
    self%i_iter = self%i_iter + 1

    ! Save models etc.
    !if (self%diagnostic_mode) then
    !   self%likelihood_saved(self%i_iter) = self%log_likelihood
    !   self%k_saved(self%i_iter) = self%tm%get_k()
    !   self%temp_saved(self%i_iter) = self%temp
    !end if
    

    return 
  end subroutine mcmc_judge_model
  
  !---------------------------------------------------------------------

  subroutine mcmc_one_step_summary(self)
    class(mcmc), intent(in) :: self
    write(*,*)"Iteration    : ", self%i_iter, "/", self%n_iter
    write(*,*)"Proposal type: ", self%i_proposal_type
    write(*,*)"Accepted     : ", self%is_accepted
    write(*,*)"Likelihood   : ", self%log_likelihood
    write(*,*)
  end subroutine mcmc_one_step_summary

  !---------------------------------------------------------------------
  
  subroutine mcmc_set_temp(self, temp)
    class(mcmc), intent(inout) :: self
    double precision, intent(in) :: temp

    self%temp = temp

    return 
  end subroutine mcmc_set_temp

  !---------------------------------------------------------------------

  double precision function mcmc_get_temp(self) result(temp)
    class(mcmc), intent(in) :: self

    temp = self%temp
    
    return 
  end function mcmc_get_temp

  !---------------------------------------------------------------------

  double precision function mcmc_get_log_likelihood(self) result(l)
    class(mcmc), intent(in) :: self
    
    l = self%log_likelihood
    
    return 
  end function mcmc_get_log_likelihood
  
  !---------------------------------------------------------------------

  function mcmc_get_likelihood_saved(self) result(l)
    class(mcmc), intent(in) :: self
    double precision :: l(self%n_iter)
    
    l(:) = self%likelihood_saved(:)
    
    return 
  end function mcmc_get_likelihood_saved
    
  !---------------------------------------------------------------------
  
  function mcmc_get_temp_saved(self) result(t)
    class(mcmc), intent(in) :: self
    double precision :: t(self%n_iter)
    
    t(:) = self%temp_saved(:)
    
    return 
  end function mcmc_get_temp_saved
  
  !---------------------------------------------------------------------

  function mcmc_get_n_accept(self) result(n_accept)
    class(mcmc), intent(in) :: self
    integer :: n_accept(self%n_proposal_type)

    n_accept = self%n_accept

    return 
  end function mcmc_get_n_accept
  
  !---------------------------------------------------------------------

  function mcmc_get_n_propose(self) result(n_propose)
    class(mcmc), intent(in) :: self
    integer :: n_propose(self%n_proposal_type)

    n_propose = self%n_propose

    return 
  end function mcmc_get_n_propose
  
  !---------------------------------------------------------------------

  function mcmc_get_n_iter(self) result(n_iter)
    class(mcmc), intent(in) :: self
    integer :: n_iter

    n_iter = self%n_iter

    return 
  end function mcmc_get_n_iter
  
  !---------------------------------------------------------------------

  function mcmc_get_hypo(self) result(hypo)
    class(mcmc), intent(in) :: self
    type(model) :: hypo
    
    hypo = self%hypo
    
    return 
  end function mcmc_get_hypo

  !---------------------------------------------------------------------

  function mcmc_get_vs(self) result(vs)
    class(mcmc), intent(in) :: self
    type(model) :: vs
    
    vs = self%vs
    
    return 
  end function mcmc_get_vs

  !---------------------------------------------------------------------

  function mcmc_get_qs(self) result(qs)
    class(mcmc), intent(in) :: self
    type(model) :: qs
    
    qs = self%qs
    
    return 
  end function mcmc_get_qs

  !---------------------------------------------------------------------

  function mcmc_get_t_corr(self) result(t_corr)
    class(mcmc), intent(in) :: self
    type(model) :: t_corr
    
    t_corr = self%t_corr
    
    return 
  end function mcmc_get_t_corr

  !---------------------------------------------------------------------

  function mcmc_get_a_corr(self) result(a_corr)
    class(mcmc), intent(in) :: self
    type(model) :: a_corr
    
    a_corr = self%a_corr
    
    return 
  end function mcmc_get_a_corr
  
  !---------------------------------------------------------------------
  
  function mcmc_write_out_hypo(self) result(out)
    class(mcmc), intent(inout) :: self
    double precision :: out(3*self%n_events)
    
    out = self%hypo%get_all_x()
    
    return 
  end function mcmc_write_out_hypo

  !---------------------------------------------------------------------

  function mcmc_write_out_t_corr(self) result(out)
    class(mcmc), intent(inout) :: self
    double precision :: out(self%n_sta)
    
    out = self%t_corr%get_all_x()
    
    return 
  end function mcmc_write_out_t_corr

  !---------------------------------------------------------------------

  double precision function mcmc_write_out_vs(self) result(out)
    class(mcmc), intent(inout) :: self
    
    out = self%vs%get_x(1)
    
    return 
  end function mcmc_write_out_vs

  !---------------------------------------------------------------------

  function mcmc_write_out_a_corr(self) result(out)
    class(mcmc), intent(inout) :: self
    double precision :: out(self%n_sta)
    
    out = self%a_corr%get_all_x()
    
    return 
  end function mcmc_write_out_a_corr

  !---------------------------------------------------------------------

  double precision function mcmc_write_out_qs(self) result(out)
    class(mcmc), intent(inout) :: self
    
    out = self%qs%get_x(1)
    
    return 
  end function mcmc_write_out_qs

  !---------------------------------------------------------------------
  
end module cls_mcmc
