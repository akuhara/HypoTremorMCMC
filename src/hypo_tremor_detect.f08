program main
  use mod_mpi
  use cls_line_text, only: line_max
  use cls_detector, only: detector
  use cls_param, only: param
  use cls_c3_data, only: c3_data
  use, intrinsic :: iso_fortran_env, only: iostat_end
  implicit none 
  integer :: n_args, ierr, rank, n_procs, id_start, id_end, i_sta
  integer :: ios, io_env
  integer, allocatable :: rank_in_charge(:)
  character(line_max) :: param_file, env_file
  character(line_max), allocatable :: stnm(:)
  type(detector) :: dtct
  type(param) :: para
  type(c3_data), allocatable :: env(:)
  double precision, allocatable :: data(:,:)
  double precision :: tmp
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

  ! Read envelope files
  allocate(stnm(para%get_n_stations()))
  stnm = para%get_stations()
  do i_sta = id_start, id_end
     allocate(data(1,1))
     write(env_file,'(a,a)')trim(stnm(i_sta)), ".merged.env"
     open(newunit=io_env, file=env_file, status='old', iostat=ios, &
          & access='stream', form='unformatted')
     if (ios /= 0) then
        write(*,*)"ERROR: cannot open ", trim(env_file)
        call mpi_finalize(ierr)
        stop
     end if
     do
        read(io_env, iostat=ios) tmp ! Time
        read(io_env, iostat=ios) tmp ! Amp
        data(:,1) = [data(:,1), tmp]
     end do
     close(io_env)
     deallocate(data)
  end do
  

  call mpi_finalize(ierr)
  
  stop
end program main
