module cls_model
  use mod_random
  implicit none 
  
  type model
     private
     
     integer :: nx  = 0 ! Number of parameters
     
     double precision, allocatable :: x(:)
     double precision, allocatable :: mu(:), sigma(:)
     double precision, allocatable :: step_size(:) ! STDEV

     logical :: verb = .false.
     
   contains
     procedure :: set_prior => model_set_prior
     procedure :: set_perturb => model_set_perturb
     procedure :: get_nx => model_get_nx
     procedure :: get_x => model_get_x
     procedure :: get_all_x => model_get_all_x
     procedure :: set_x => model_set_x
     procedure :: generate_model => model_generate_model
     procedure :: perturb => model_perturb
     procedure :: display => model_display
     
  end type model
  
  interface model
     module procedure init_model
  end interface model
  
contains

  !---------------------------------------------------------------------
  
  type(model) function init_model(nx, verb) result(self)
    integer, intent(in) :: nx
    logical, intent(in), optional :: verb
    
    if (present(verb)) then
       self%verb = verb
    end if
    if (self%verb) then
       write(*,'(A)')"<< Initialize hyper model parameters >>"
    end if
    
    ! Get number of model parameters
    self%nx = nx
    allocate(self%x(nx))
    allocate(self%mu(nx))
    allocate(self%sigma(nx))
    allocate(self%step_size(nx))
    
    if (self%verb) write(*,'(A,I4)')" nx=", nx
    if (self%verb) write(*,*)
    
    return 
  end function init_model
  
  !---------------------------------------------------------------------
  
  subroutine model_set_prior(self, i, mu, sigma)
    class(model), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: mu, sigma
    
    self%mu(i) = mu
    self%sigma(i) = sigma
    
    return 
  end subroutine model_set_prior

  !---------------------------------------------------------------------
  
  subroutine model_set_perturb(self, i, step_size)
    class(model), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: step_size

    self%step_size(i) = step_size
    
    return 
  end subroutine model_set_perturb
  
  !---------------------------------------------------------------------

  integer function model_get_nx(self) result(nx)
    class(model), intent(in) :: self
    
    nx = self%nx

    return 
  end function model_get_nx

  !---------------------------------------------------------------------
  
  double precision function model_get_x(self, iparam) result(x)
    class(model), intent(in) :: self
    integer, intent(in) :: iparam
    
    x = self%x(iparam)
    
    return 
  end function model_get_x

  !---------------------------------------------------------------------

  function model_get_all_x(self) result(all_x)
    class(model), intent(in) :: self
    double precision :: all_x(self%nx)
    
    all_x = self%x
    
    return 
  end function model_get_all_x

  !---------------------------------------------------------------------

  subroutine model_set_x(self, iparam, x)
    class(model), intent(inout) :: self
    integer, intent(in) :: iparam
    double precision, intent(in) :: x

    self%x(iparam) = x
    
    return 
  end subroutine model_set_x

  !---------------------------------------------------------------------

  subroutine model_generate_model(self)
    class(model), intent(inout) :: self
    integer :: i

    do i = 1, self%nx
       self%x(i) = self%mu(i) + rand_g() * self%sigma(i)
    end do

    return 
  end subroutine model_generate_model
  
  !---------------------------------------------------------------------

  subroutine model_perturb(self, i, log_prior_ratio)
    class(model), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(out) :: log_prior_ratio
    double precision :: x_new, x_old


    x_old = self%x(i)
    x_new = x_old + rand_g() * self%step_size(i)
    self%x(i) = x_new

    log_prior_ratio = &
         & -((x_new- self%mu(i))**2 - (x_old-self%mu(i))**2) / &
         & (2.d0 * self%sigma(i) * self%sigma(i))

    return 
  end subroutine model_perturb

  !---------------------------------------------------------------------
  
  subroutine model_display(self)
    class(model), intent(inout) :: self
    integer :: i
    
    do i = 1, self%nx
       write(*,*)"------------------------------"
       write(*,*)"Parameter", i
       write(*,*)self%x(i)
    end do

    write(*,*)
    
    return 
  end subroutine model_display

  !---------------------------------------------------------------------


end module cls_model
  
