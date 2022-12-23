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
       & sta_x         = para%get_sta_x(), &
       & sta_y         = para%get_sta_y(), &
       & sta_z         = para%get_sta_z(), &
       & t_win         = para%get_t_win_corr(), &
       & t_step        = para%get_t_step_corr(), &
       & n_pair_thred  = para%get_n_pair_thred(), &
       & verb          = verb)
  call msr%scan_cc()
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  call msr%measure_lag_time()
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  call mpi_finalize(ierr)
  
  stop
end program main
