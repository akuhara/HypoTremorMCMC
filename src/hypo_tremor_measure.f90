program main
  use mod_mpi
  use cls_param, only: param
  use cls_line_text, only: line_max
  use cls_measurer, only: measurer
  implicit none 
  integer :: n_args, ierr, rank, n_procs, n_pair
  type(param) :: para
  type(measurer) :: msr
  character(line_max) :: param_file
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
  para = param(param_file, verb=verb)
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  
  msr = measurer(&
       & station_names = para%get_stations(), & 
       & verb          = verb)

  
  call mpi_finalize(ierr)
  


  stop
end program main
