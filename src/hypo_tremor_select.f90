program main
  use mod_mpi
  use cls_param, only: param
  use cls_line_text, only: line_max
  use cls_selector, only: selector
  implicit none   
  integer :: i, i1, i2
  integer :: n_args, ierr, n_procs, rank, io, ios, n_win, io2, io3
  integer, allocatable :: win_id(:), rank_in_charge(:)
  double precision, allocatable :: dummy(:)
  double precision, allocatable :: vs(:), b(:), t0(:), a0(:), cc_t(:), cc_a(:)
  double precision, allocatable :: vs_all(:), b_all(:), t0_all(:), a0_all(:)
  double precision, allocatable :: cc_t_all(:), cc_a_all(:)
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

  ! Check if MPI process number is as intended
  if (verb .and. n_procs /= para%get_n_procs()) then
     write(*,*) "ERROR: n_procs in parameter file must be " // &
          & "equal to that is given in the command line"
     call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
  end if

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
  
  allocate(win_id(n_win), dummy(n_win))

  rewind(io)
  do i = 1, n_win
     read(io, *)win_id(i), dummy(i)
  end do
  close(io)
  

  ! Allocate arrays for results
  allocate(vs(n_win), b(n_win), t0(n_win), a0(n_win))
  allocate(cc_t(n_win), cc_a(n_win))
  allocate(vs_all(n_win), b_all(n_win))
  allocate(t0_all(n_win), a0_all(n_win))
  allocate(cc_t_all(n_win), cc_a_all(n_win))

  vs(:)  = 0.d0
  b(:)   = 0.d0
  t0(:) = 0.d0
  a0(:)  = 0.d0

  allocate(rank_in_charge(n_win))
  call get_mpi_task_id(n_win, i1, i2, rank_in_charge, .false.)
  
  slct = selector(&
       & station_names = para%get_stations(), &
       & sta_x   = para%get_sta_x(), &
       & sta_y   = para%get_sta_y(), &
       & sta_z   = para%get_sta_z(), &
       & z_guess = para%get_z_guess() &
       & )
  
  do i = i1, i2
     call slct%eval_wave_propagation(win_id(i), vs(i), t0(i), b(i), a0(i), &
          & cc_t(i), cc_a(i))
  end do
  
  call mpi_reduce(vs(1:n_win), vs_all(1:n_win), n_win, MPI_DOUBLE_PRECISION, &
       & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call mpi_reduce(t0(1:n_win), t0_all(1:n_win), n_win, MPI_DOUBLE_PRECISION, &
       & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call mpi_reduce(b(1:n_win), b_all(1:n_win), n_win, MPI_DOUBLE_PRECISION, &
       & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call mpi_reduce(a0(1:n_win), a0_all(1:n_win), n_win, MPI_DOUBLE_PRECISION, &
       & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call mpi_reduce(cc_t(1:n_win), cc_t_all(1:n_win), n_win, MPI_DOUBLE_PRECISION, &
       & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  call mpi_reduce(cc_a(1:n_win), cc_a_all(1:n_win), n_win, MPI_DOUBLE_PRECISION, &
       & MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (rank == 0) then
     open(newunit=io2, status="replace", file="regress.dat", iostat=ios)
     if (ios /= 0) then
        call mpi_finalize(ierr)
        error stop "ERROR: cannot create regress.dat"
     end if

     open(newunit=io3, status="replace", file="selected_win.dat", iostat=ios)
     if (ios /= 0) then
        call mpi_finalize(ierr)
        error stop "ERROR: cannot create selected_win.dat"
     end if
     
     do i = 1, n_win
        write(io2,*)win_id(i), vs_all(i), b_all(i), t0_all(i), a0_all(i), &
             & cc_t_all(i), cc_a_all(i)
        if (   vs_all(i) >= para%get_vs_min() .and. &
             & vs_all(i) <= para%get_vs_max() .and. &
             & b_all(i)  >= para%get_b_min() .and. &
             & b_all(i)  <= para%get_b_max()) then
           write(io3, *)win_id(i), dummy(i)
        end if
     end do

     close(io2)
  end if
  call mpi_finalize(ierr)
  
  stop
end program main
