program main
  use mod_mpi
  use cls_param, only: param
  use cls_line_text, only: line_max
  use cls_selector, only: selector
  implicit none   
  integer :: i, i1, i2
  integer :: n_args, ierr, n_procs, rank, io, ios, n_win
  integer, allocatable :: win_id(:)
  double precision :: dummy
  character(line_max) :: param_file
  type(param) :: para
  type(selector) :: slct
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

  ! Read detected-window file
  open(newunit=io, file='detected_win.dat', status='old', iostat=ios)
  if (ios /= 0) then
     call mpi_finalize(ierr)
     error stop "ERROR: detected_win.dat is not found"
  end if
  n_win = 0
  do 
     read(io, *, iostat=ios)
     if (ios /= 0) exit
     n_win = n_win + 1
  end do
  
  allocate(win_id(n_win))

  rewind(io)
  do i = 1, n_win
     read(io, *)win_id(i), dummy
  end do
  close(io)

  call get_mpi_task_id(n_win, i1, i2)

  slct = selector(&
       & win_id =win_id(i1:i2), &
       & station_names = para%get_stations(), &
       & sta_x = para%get_sta_x(), &
       & sta_y = para%get_sta_y(), &
       & sta_z = para%get_sta_z(), &
       & vs_min = para%get_vs_min(), &
       & vs_max = para%get_vs_max(), &
       & b_min = para%get_b_min(), &
       & b_max = para%get_b_max() &
       & )

  
  
end program main
