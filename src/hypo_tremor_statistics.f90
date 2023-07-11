program main
  use mod_mpi
  use cls_param, only: param
  use cls_line_text, only: line_max
  use cls_statistics, only: statistics
  implicit none
  
  type(param) :: para
  type(statistics) :: st
  character(line_max) :: param_file
  integer :: n_procs, rank, ierr, n_args
  logical :: verb

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
  para = param(param_file, verb=verb, from_where="select")
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Check if MPI process number is as intended
  if (verb .and. n_procs /= para%get_n_procs()) then
     write(*,*) "ERROR: n_procs in parameter file must be " // &
          & "equal to that is given in the command line"
     call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
  end if
  
  write(*,*)verb
  st = statistics( &
       & n_procs       = n_procs, &
       & n_iter        = para%get_n_iter(), &
       & n_burn        = para%get_n_burn(), &
       & n_interval    = para%get_n_interval(), &
       & n_cool        = para%get_n_cool(), &
       & station_names = para%get_stations(), &
       & sta_x         = para%get_sta_x(), &
       & sta_y         = para%get_sta_y(), &
       & sta_z         = para%get_sta_z(), &
       & verb          = verb &
       & )

  call st%estimate_vs_qs()
  call st%estimate_corr_factors()
  call st%estimate_hypo()

  call mpi_finalize(ierr)
  stop
end program main
