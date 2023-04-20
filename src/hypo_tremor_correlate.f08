program main
  use mod_mpi
  use cls_line_text, only: line_max
  use cls_correlator, only: correlator
  use cls_param, only: param
  use cls_c3_data, only: c3_data
  use, intrinsic :: iso_fortran_env, only: iostat_end
  implicit none 
  integer :: n_args, ierr, rank, n_procs, id_start, id_end, i_sta
  integer :: ios, io_env, n_smp, i
  integer, allocatable :: rank_in_charge(:)
  character(line_max) :: param_file, env_file
  character(line_max), allocatable :: stnm(:)
  type(correlator) :: corr
  type(param) :: para
  type(c3_data), allocatable :: env(:)
  double precision, allocatable :: data(:,:)
  double precision :: t1, t2, amp, dt
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
  para = param(param_file, verb=verb, from_where="correlate")
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  
  ! Get task ID
  allocate(rank_in_charge(para%get_n_stations()))
  call get_mpi_task_id(para%get_n_stations(), id_start, id_end, &
       & rank_in_charge, debug=.false.)
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! Read envelope files
  allocate(stnm(para%get_n_stations()))
  allocate(env(para%get_n_stations()))
  stnm = para%get_stations()
  do i_sta = id_start, id_end

     write(env_file,'(a,a)')trim(stnm(i_sta)), ".merged.env"
     open(newunit=io_env, file=env_file, status='old', iostat=ios, &
          & access='stream', form='unformatted')
     if (ios /= 0) then
        write(*,*)"ERROR: cannot open ", trim(env_file)
        call mpi_finalize(ierr)
        stop
     end if
     print *, "Now reading ", trim(env_file)
     n_smp = 0
     do
        read(io_env, iostat=ios) t1 ! Time
        read(io_env, iostat=ios) amp ! Amp
        if (ios == iostat_end) exit
        n_smp = n_smp + 1
     end do
     allocate(data(n_smp, 1))
     rewind(io_env)
     t1 = 0.d0
     do i = 1, n_smp
        read(io_env) t2
        read(io_env) data(i,1)
        dt = t2 - t1
        t1 = t2
     end do
     close(io_env)

     ! Sotre envelope data into object
     env(i_sta) = c3_data(n_cmps=1, dt=dt)
     call env(i_sta)%enqueue_data(data)
     
     deallocate(data)
  end do
  
  
  ! Broadcast n_smp and dt from rank 0 process
  !  Just in case there are processes not reading envelope files at all
  call mpi_bcast(n_smp, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  
  do i_sta = 1, para%get_n_stations()
     allocate(data(n_smp, 1))
     
     ! Store envelope data to be broadcasted
     if (rank == rank_in_charge(i_sta)) then
        data   = env(i_sta)%get_data()
     end if
     
     ! Broadcast envelope data from from a process in charge to all processes
     call mpi_bcast(data, n_smp, MPI_DOUBLE_PRECISION, &
          & rank_in_charge(i_sta), MPI_COMM_WORLD, ierr)

     ! Store data into envelope object
     if (rank /= rank_in_charge(i_sta)) then
        env(i_sta) = c3_data(dt=dt, n_cmps=1)
        call env(i_sta)%enqueue_data(data)
     end if
     
     deallocate(data)
  end do
  
  ! Detect
  corr = correlator(&
       & station_names = para%get_stations(),    &
       & t_win         = para%get_t_win_corr(),  &
       & t_step        = para%get_t_step_corr(), &
       & dt            = env(1)%get_dt(),        &
       & n_smp         = env(1)%get_n_smp(),     &
       & alpha         = para%get_alpha()        &
       & )
  call corr%calc_correlogram(env=env)
  
  

  call mpi_finalize(ierr)
  
  stop
end program main
