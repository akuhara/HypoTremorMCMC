module cls_obs_data
  use cls_line_text, only: line_max
  implicit none 
  type obs_data
     private
     integer :: n_events
     integer :: n_sta
     integer, allocatable :: win_id(:)
     logical :: verb
     
     ! relative arrival time
     double precision, allocatable :: t_obs(:,:)
     double precision, allocatable :: t_stdv(:,:)
     double precision, allocatable :: dt_obs(:,:,:)
     double precision, allocatable :: dt_stdv(:,:,:)
     double precision, allocatable :: sta_x(:), sta_y(:)
     
     
   contains
     procedure :: read_obs_files => obs_data_read_obs_files
     procedure :: make_initial_guess => obs_data_make_initial_guess
     procedure :: get_dt_obs => obs_data_get_dt_obs
     procedure :: get_dt_stdv => obs_data_get_dt_stdv
  end type obs_data
  
  interface obs_data
     module procedure init_obs_data
  end interface obs_data

contains
  
  !-------------------------------------------------------------------------
  
  type(obs_data) function init_obs_data(win_id, n_sta, &
       & sta_x, sta_y, verb) result(self)
    integer, intent(in) :: win_id(:), n_sta
    double precision, intent(in) :: sta_x(:), sta_y(:)
    logical, intent(in) :: verb

    self%n_events = size(win_id)
    allocate(self%win_id(self%n_events))
    self%win_id = win_id
    self%verb = verb


    self%n_sta = n_sta
    allocate(self%t_obs(self%n_sta, self%n_events))
    allocate(self%t_stdv(self%n_sta, self%n_events))
    allocate(self%dt_obs(self%n_sta, self%n_sta, self%n_events))
    allocate(self%dt_stdv(self%n_sta, self%n_sta, self%n_events))

    allocate(self%sta_x(self%n_sta))
    allocate(self%sta_y(self%n_sta))
    self%sta_x = sta_x
    self%sta_y = sta_y
    
    
    call self%read_obs_files()
    
    return 
  end function init_obs_data

  !-------------------------------------------------------------------------
  
  subroutine obs_data_read_obs_files(self )
    class(obs_data), intent(inout) :: self
    integer :: i, io, ierr, j, k
    character(line_max) :: obs_file
    double precision :: dummy1, dummy2, dummy3, dummy4

    if (self%verb) then
       print *, "<< Reading obs files>>"
    end if
    do i = 1, self%n_events
       write(obs_file,'(A,I5.5,A)')"rel_time.", self%win_id(i), ".dat"
       open(newunit=io, file=obs_file, status="old", iostat=ierr)
       if (ierr /= 0) then
          print *, trim(obs_file)
          error stop "ERROR: obs_file is not found"
       end if

       do j = 1, self%n_sta
          read(io,*) dummy1, dummy2, self%t_obs(j,i), self%t_stdv(j,i), &
               & dummy3, dummy4
          if (i == 1 .and. self%verb) then
             print *, "T=", self%t_obs(j, i), "T_stdv=", self%t_stdv(j,i)
          end if
       end do

       ! Make dt
       do j = 1, self%n_sta - 1
          do k = j + 1, self%n_sta
             self%dt_obs(k, j, i) = self%t_obs(j, i) - self%t_obs(k, i)
             self%dt_obs(j, k, i) = -self%dt_obs(k, j, i)
             self%dt_stdv(k, j, i) = sqrt(self%t_stdv(j, i)**2 + &
                  & self%t_stdv(k, i)**2)
             self%dt_stdv(j, k, i) = self%dt_stdv(k, j, i)
          end do
       end do
       close(io)
       
    end do
    

    return 
  end subroutine obs_data_read_obs_files
  
  !-------------------------------------------------------------------------
  
  subroutine obs_data_make_initial_guess(self, x_mu, y_mu)
    class(obs_data), intent(inout) :: self
    double precision, intent(out) :: x_mu(self%n_events), y_mu(self%n_events)
    integer :: i, ista
    
    
    do i = 1, self%n_events
       ista = minloc(self%t_obs(:,i), dim=1)
       x_mu(i) = self%sta_x(ista)
       y_mu(i) = self%sta_y(ista)
    end do
    
    
    return 
  end subroutine obs_data_make_initial_guess

  !-------------------------------------------------------------------------
  
  function obs_data_get_dt_obs(self) result(dt_obs)
    class(obs_data), intent(in) :: self
    double precision :: dt_obs(self%n_sta, self%n_sta, self%n_events)
    
    dt_obs  = self%dt_obs
    
    return 
  end function obs_data_get_dt_obs

  !-------------------------------------------------------------------------
  
  function obs_data_get_dt_stdv(self) result(dt_stdv)
    class(obs_data), intent(in) :: self
    double precision :: dt_stdv(self%n_sta, self%n_sta, self%n_events)
    
    dt_stdv  = self%dt_stdv
    
    return 
  end function obs_data_get_dt_stdv
    
  !-------------------------------------------------------------------------
    
    
  
end module cls_obs_data
