program main
  use mod_mpi
  use param: only: param
  use cls_line_text, only: line_max
  implicit none   
  integer :: n_args, ierr, n_procs, rank, io, ios, n_win
  type(param) :: para
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
     if (iostat /= 0) exit
     n_win = n_win + 1
  end do
  rewind(io)
  close(io)
  
end program main
