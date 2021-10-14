program main
  use mod_mpi
  use mod_random
  use cls_model, only: model
  use cls_obs_data, only: obs_data
  use cls_param, only:param
  use cls_line_text, only: line_max
  use, intrinsic :: iso_fortran_env, only: iostat_end
  
  implicit none
  integer :: n_args, ierr, rank, n_procs
  type(model) :: hypo
  type(model) :: t_corr
  type(param) :: para
  type(obs_data) :: obs
  
  character(line_max) :: param_file, win_id_file
  logical :: verb
  integer, allocatable :: win_id(:)
  integer :: i, io, id, n_events
  double precision :: dummy
  
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, n_procs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)


  if (rank == 0) then
     verb = .true.
  else
     verb = .false.
  end if

  ! Get argument
  n_args = command_argument_count()
  if (n_args /= 1) then
     error stop "USAGE: hypo_tremor_mcmc [parameter file]"
  end if
  call get_command_argument(1, param_file)

  
  ! Read parameter file
  para = param(param_file, verb=verb)
  call mpi_barrier(MPI_COMM_WORLD, ierr)


  !
  call init_random(1111, 222222, 313131, 55555555, rank)
  
  ! Events
  win_id = [ integer :: ]
  win_id_file = "detected_win.dat"
  open(newunit=io, file=win_id_file, status="old", iostat=ierr)
  if (ierr /= 0) then
     error stop "cannot open detected_win.dat"
  end if
  do 
     read(io, *, iostat=ierr) id, dummy
     if (ierr == iostat_end) exit 
     win_id = [win_id, id]
  end do
  n_events = size(win_id)
  close(io)
  
  ! Observed data
  obs = obs_data( &
       & win_id = win_id, &
       & n_sta  = para%get_n_stations(), &
       & verb   = verb)



  ! MCMC model parameters
  ! ** Station correction
  t_corr = model(&
       & nx   = para%get_n_stations(), &
       & verb = verb)
  do i = 1, para%get_n_stations()
     call t_corr%set_prior(i=i, mu=0.d0, sigma=1.d0) 
  end do
  call t_corr%generate_model()
  call t_corr%display()
  
  ! ** Hypocenters
  hypo = model(&
       & nx   = 3 * n_events, &
       & verb = verb)
    
  
  
  call mpi_finalize(ierr)

  stop
end program main
