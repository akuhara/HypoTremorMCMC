program hypo_tremor_mcmc
  use mod_mpi
  use cls_param
  use cls_line_text, only: line_max
  use cls_convertor
  implicit none 
  integer :: n_args, ierr, rank, n_procs, id_start, id_end
  character(line_max) :: param_file
  type(param) :: para
  type(convertor) :: conv
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

  ! Get task ID
  call get_mpi_task_id(para%get_n_stations(), id_start, id_end, debug=.true.)
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  
  ! Convert raw data to CF
  
  do i_sta = id_start, id_end
     conv = convertor(&
          & n_cmps = para%get_n_cmps(), &
          & t_win  = para%get_t_win() , &
          & filenames = para%get_filenames(1))
     call conv%convert(i_sta)
  end do
  
  call mpi_finalize(ierr)
  
  stop
end program hypo_tremor_mcmc

