program hypo_tremor_convert
  use mod_mpi
  use cls_param, only: param
  use cls_line_text, only: line_max
  use cls_convertor, only: convertor
  implicit none 
  integer :: n_args, ierr, rank, n_procs, id_start, id_end, i_sta, n_pair
  integer, allocatable :: rank_in_charge(:)
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
  para = param(param_file, verb=verb, from_where="convert")
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  
  !Check if MPI process number is as intended
  if (verb .and. n_procs /= para%get_n_procs()) then
     write(*,*) "ERROR: n_procs in parameter file must be " // &
          & "equal to that is given in the command line"
     call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
  end if


  ! Get task ID
  allocate(rank_in_charge(para%get_n_stations()))
  call get_mpi_task_id(para%get_n_stations(), id_start, id_end, rank_in_charge, .false.)
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Convert raw data to smoothed envelope
  do i_sta = id_start, id_end
     conv = convertor(&
          & t_win        = para%get_t_win_conv() ,        &
          & filenames    = para%get_filenames(i_sta),     &
          & station_name = para%get_station(i_sta),       &
          & sta_amp_fac  = para%get_sta_amp_fac(i_sta)    &
          & )

     call conv%convert()
  end do
  

  call mpi_finalize(ierr)
  stop
end program hypo_tremor_convert

