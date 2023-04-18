program hypo_tremor_mcmc
  use mod_mpi
  use cls_param, only: param
  use cls_line_text, only: line_max
  use cls_convertor, only: convertor
  use cls_c3_data, only: c3_data
  use cls_detector, only: detector
  use cls_measurer, only: measurer
  implicit none 
  integer :: n_args, ierr, rank, n_procs, id_start, id_end, i_sta, n_pair
  integer, allocatable :: rank_in_charge(:)
  character(line_max) :: param_file
  type(param) :: para
  type(convertor) :: conv
  type(c3_data), allocatable :: env(:)
  type(detector) :: dtct
  type(measurer) :: msr
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
  para = param(param_file, verb=verb, from_where="optimize")
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Get task ID
  allocate(rank_in_charge(para%get_n_stations()))
  call get_mpi_task_id(para%get_n_stations(), id_start, id_end, &
       & rank_in_charge, debug=.true.)
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Convert raw data to smoothed envelope
  allocate(env(para%get_n_stations()))
  do i_sta = id_start, id_end
     conv = convertor(&
          & t_win        = para%get_t_win_conv() ,        &
          & filenames    = para%get_filenames(i_sta),     &
          & comps        = para%get_comps(),              &
          & station_name = para%get_station(i_sta),       &
          & sta_amp_fac  = para%get_sta_amp_fac(i_sta)    &
          & )

     call conv%convert()
     env(i_sta) = conv%get_c3_out()
  end do
  
  ! Gather results of all processes
  block 
    integer :: i_sta, ierr, n_smp, n_cmps
    double precision :: dt
    double precision, allocatable :: data(:,:)

    do i_sta = 1, para%get_n_stations()
       if (rank == rank_in_charge(i_sta)) then
          n_smp  = env(i_sta)%get_n_smp()
          n_cmps = env(i_sta)%get_n_cmps()
          dt     = env(i_sta)%get_dt()
          data   = env(i_sta)%get_data()
       end if
       call mpi_bcast(n_smp, 1, MPI_INTEGER4, rank_in_charge(i_sta), &
            & MPI_COMM_WORLD, ierr)
       call mpi_bcast(n_cmps, 1, MPI_INTEGER4, rank_in_charge(i_sta), &
            & MPI_COMM_WORLD, ierr)
       call mpi_bcast(dt, 1, MPI_DOUBLE_PRECISION, rank_in_charge(i_sta), &
            & MPI_COMM_WORLD, ierr)
       if (rank /= rank_in_charge(i_sta)) then
          allocate(data(n_smp, n_cmps))
          env(i_sta) = c3_data(dt=dt, n_cmps=n_cmps)
       end if
       call mpi_bcast(data, n_smp * n_cmps, MPI_DOUBLE_PRECISION, &
            & rank_in_charge(i_sta), MPI_COMM_WORLD, ierr)
       if (rank /= rank_in_charge(i_sta)) then
          call env(i_sta)%enqueue_data(data)
       end if
       deallocate(data)
    end do
    
  end block
  

  ! Detect

  dtct = detector(&
       & station_names = para%get_stations(),    &
       & t_win         = para%get_t_win_corr(),  &
       & t_step        = para%get_t_step_corr(), &
       & dt            = env(1)%get_dt(),        &
       & n_smp         = env(1)%get_n_smp(),     &
       & alpha         = para%get_alpha()        &
       & )
  call dtct%calc_correlogram(env=env)
  
  call mpi_barrier(MPI_COMM_WORLD, ierr)


  ! Time- and amplitude-difference optimization
  !  Parallelized inside
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
end program hypo_tremor_mcmc

