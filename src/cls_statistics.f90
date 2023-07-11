module cls_statistics
  use mod_sort
  use cls_line_text, only: line_max
  use mod_mpi
  implicit none
  type statistics
     private

     ! MPI processes
     integer :: n_procs, rank
     
     ! # of model sampled by MCMC
     integer :: n_mod
     
     ! Stations
     integer :: n_sta
     character(line_max), allocatable :: station_names(:)
     double precision, allocatable :: sta_x(:), sta_y(:), sta_z(:)
     
     integer:: n_evt
     integer, allocatable :: win_id(:)
     
     ! sampled values
     double precision, allocatable :: t_corr(:,:), a_corr(:,:)
     double precision, allocatable :: vs(:), qs(:)
     double precision, allocatable :: hypo_x(:,:), hypo_y(:,:), hypo_z(:,:)
     
     ! verb
     logical :: verb

   contains
     procedure :: estimate_hypo => statistics_estimate_hypo
     procedure :: estimate_corr_factors => statistics_estimate_corr_factors
     procedure :: estimate_vs_qs => statistics_estimate_vs_qs
     procedure :: read_hypo_file => statistics_read_hypo_file
     procedure :: read_corr_file => statistics_read_corr_file
     procedure :: read_v_file => statistics_read_v_file
     
     
  end type statistics
  
  interface statistics
     module procedure init_statistics
  end interface statistics
  
contains
  
  !-----------------------------------------------------------------------
  
  type(statistics) function init_statistics(n_procs, &
       & n_iter, n_burn, n_interval, n_cool, &
       & station_names, sta_x, sta_y, sta_z, verb) &
       & result(self)
    integer, intent(in) :: n_procs
    integer, intent(in) :: n_iter, n_burn, n_interval, n_cool
    character(line_max), intent(in) :: station_names(:)
    double precision, intent(in) :: sta_x(:), sta_y(:), sta_z(:)
    logical, intent(in), optional :: verb
    integer :: ierr, io, ios, i
    double precision :: dummy
    
    self%n_procs = n_procs
    call mpi_comm_rank(MPI_COMM_WORLD, self%rank, ierr)
    
    self%n_mod = (n_iter - n_burn)  * n_procs * n_cool / n_interval
    
    self%n_sta = size(station_names(:))
    allocate(self%station_names(self%n_sta))
    allocate(self%sta_x(self%n_sta))
    allocate(self%sta_y(self%n_sta))
    allocate(self%sta_z(self%n_sta))
    self%station_names = station_names
    self%sta_x = sta_x
    self%sta_y = sta_y
    self%sta_z = sta_z

    allocate(self%t_corr(self%n_mod, self%n_sta))
    allocate(self%a_corr(self%n_mod, self%n_sta))
    allocate(self%vs(self%n_mod))
    allocate(self%qs(self%n_mod))

    open(newunit=io, file="selected_win.dat", status="old", &
         & form="formatted", iostat=ios)
    if (self%verb .and. ios /= 0) then
       write(0,*)"ERROR: selected_win.dat is not found."
       call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
    end if
    self%n_evt = 0
    do 
       read(io,*,iostat=ios)
       if (ios /= 0) exit
       self%n_evt = self%n_evt + 1
    end do



    allocate(self%hypo_x(self%n_mod, self%n_evt))
    allocate(self%hypo_y(self%n_mod, self%n_evt))
    allocate(self%hypo_z(self%n_mod, self%n_evt))
    allocate(self%win_id(self%n_evt))

    rewind(io)
    do i = 1, self%n_evt
       read(io, *)self%win_id(i), dummy
    end do
    close(io)

    self%verb = .false.
    if (present(verb)) then
       self%verb = verb
    end if

    return 
  end function init_statistics


  !-----------------------------------------------------------------------
  
  subroutine statistics_estimate_hypo(self)
    class(statistics), intent(inout) :: self
    integer :: ierr, i_evt, n_mod, io, ios
    integer :: il, im, iu
    character(line_max) :: hypo_file
    character(line_max) :: out_file

    out_file = "hypo.stat"

    write(hypo_file, '(a,I2.2,a)')"hypo.", self%rank, ".out"
    
    call self%read_hypo_file(hypo_file, self%hypo_x, self%hypo_y, &
         & self%hypo_z)

    il = 0.025 * self%n_mod
    im = 0.5 * self%n_mod
    iu = 0.975 * self%n_mod
    if (self%rank == 0) then

       do i_evt = 1, self%n_evt
          write(*,*)"sort", i_evt
          call quick_sort(self%hypo_x(:,i_evt), 1, self%n_mod)
          call quick_sort(self%hypo_y(:,i_evt), 1, self%n_mod)
          call quick_sort(self%hypo_z(:,i_evt), 1, self%n_mod)
       end do

       open(newunit=io, file=out_file, status="replace", form="formatted", &
            & iostat=ios)
       write(io, '(A)') "# window ID, x (50%), x (2.5%) " // &
            & "x (97.5%), y (50%), y (2.5%), y (97.5%)" // &
            & "z (50 %), z (2.5%), z (97.5%)"
       do i_evt = 1, self%n_evt
          write(io,'(I9,9F11.6)')self%win_id(i_evt), &
               & self%hypo_x(im, i_evt), &
               & self%hypo_x(il, i_evt), self%hypo_x(iu, i_evt), &
               & self%hypo_y(im, i_evt), &
               & self%hypo_y(il, i_evt), self%hypo_y(iu, i_evt), &
               & self%hypo_z(im, i_evt), &
               & self%hypo_z(il, i_evt), self%hypo_z(iu, i_evt) 
       end do
       close(io)
       
    end if


    
    return 
  end subroutine statistics_estimate_hypo
  
  !-----------------------------------------------------------------------
  
  subroutine statistics_read_hypo_file(self, hypo_file, &
       & hypo_x_all, hypo_y_all, hypo_z_all)
    class(statistics), intent(in) :: self
    character(line_max), intent(in) :: hypo_file
    double precision, intent(out) :: hypo_x_all(:,:), hypo_y_all(:,:), &
         & hypo_z_all(:,:)
    double precision :: hypo_x(self%n_mod, self%n_evt)
    double precision :: hypo_y(self%n_mod, self%n_evt)
    double precision :: hypo_z(self%n_mod, self%n_evt)
    integer :: io, ios, ierr, i_mod, i_evt, i_dummy, j_mod, n_recv
    integer :: i_rank
    integer, allocatable :: ista(:)
    
    if (self%verb) then
       write(*,'(a)')"<< Now reading " // trim(hypo_file) // " >>"
    end if
    open(newunit=io, file=hypo_file, iostat=ios, status='old', &
         & access="stream", form="unformatted")
    if (ios /= 0) then
       write(0,*)"ERROR: cannot open ", trim(hypo_file)
       call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
    end if
    
    i_mod = 0
    do 
       read(io, iostat=ios)i_dummy, (hypo_x(i_mod+1, i_evt), &
            & hypo_y(i_mod+1, i_evt), hypo_z(i_mod+1, i_evt), &
            & i_evt=1,self%n_evt)
       if (ios /= 0) exit 
       i_mod = i_mod + 1
    end do
    close(io)
    
    ! Gather
    allocate(ista(MPI_STATUS_SIZE))
    do i_evt = 1, self%n_evt
       if (self%rank == 0) then
          hypo_x_all(1:i_mod, i_evt) = hypo_x(1:i_mod, i_evt)
          hypo_y_all(1:i_mod, i_evt) = hypo_y(1:i_mod, i_evt)
          hypo_z_all(1:i_mod, i_evt) = hypo_z(1:i_mod, i_evt)
          j_mod = i_mod
          do i_rank = 1, self%n_procs - 1
             call mpi_recv(n_recv, 1, MPI_INTEGER, i_rank, 0, &
                  & MPI_COMM_WORLD, ista, ierr)
             call mpi_recv(hypo_x_all(j_mod+1, i_evt), n_recv, &
                  & MPI_DOUBLE_PRECISION, &
                  & i_rank, 1, MPI_COMM_WORLD, ista, ierr)
             call mpi_recv(hypo_y_all(j_mod+1, i_evt), n_recv, &
                  & MPI_DOUBLE_PRECISION, &
                  & i_rank, 2, MPI_COMM_WORLD, ista, ierr)
             call mpi_recv(hypo_z_all(j_mod+1, i_evt), n_recv, &
                  & MPI_DOUBLE_PRECISION, &
                  & i_rank, 3, MPI_COMM_WORLD, ista, ierr)
             j_mod = j_mod + n_recv
          end do
       else
          call mpi_send(i_mod, 1, MPI_INTEGER, 0, 0, &
               & MPI_COMM_WORLD, ista, ierr)
          call mpi_send(hypo_x(1, i_evt), i_mod, MPI_DOUBLE_PRECISION, &
               & 0, 1, MPI_COMM_WORLD, ista, ierr)
          call mpi_send(hypo_y(1, i_evt), i_mod, MPI_DOUBLE_PRECISION, &
               & 0, 2, MPI_COMM_WORLD, ista, ierr)
          call mpi_send(hypo_z(1, i_evt), i_mod, MPI_DOUBLE_PRECISION, &
               & 0, 3, MPI_COMM_WORLD, ista, ierr)
       end if
    end do
    


    return 
  end subroutine statistics_read_hypo_file
  
  !-----------------------------------------------------------------------
  
  subroutine statistics_estimate_corr_factors(self)
    class(statistics), intent(inout) :: self
    integer :: ierr, i_sta, n_mod, io, ios
    integer :: il, im, iu
    character(line_max) :: t_corr_file, a_corr_file
    character(line_max) :: out_file
    
    out_file = "station_corrections.stat"

    write(t_corr_file, '(a,I2.2,a)')"t_corr.", self%rank, ".out"
    write(a_corr_file, '(a,I2.2,a)')"a_corr.", self%rank, ".out"
    
    call self%read_corr_file(t_corr_file, self%t_corr)
    call self%read_corr_file(a_corr_file, self%a_corr)
    
    il = 0.025 * self%n_mod
    im = 0.5 * self%n_mod
    iu = 0.975 * self%n_mod
    if (self%rank == 0) then
       do i_sta = 1, self%n_sta
          write(*,*)"sort", i_sta
          call quick_sort(self%t_corr(:,i_sta), 1, self%n_mod)
          call quick_sort(self%a_corr(:,i_sta), 1, self%n_mod)

       end do

       open(newunit=io, file=out_file, status="replace", form="formatted", &
            & iostat=ios)
       write(io, '(A)') "# station name, t_corr (50%), t_corr (2.5%) " // &
            & "t_corr (97.5%), a_corr (50%), a_corr (2.5%), a_corr (97.5%)" 
       do i_sta = 1, self%n_sta
          write(io,'(A12,6F11.6)')trim(self%station_names(i_sta)), &
               & self%t_corr(im, i_sta), &
               & self%t_corr(il, i_sta), self%t_corr(iu, i_sta), &
               & self%a_corr(im, i_sta), &
               & self%a_corr(il, i_sta), self%a_corr(iu, i_sta) 
       end do
       close(io)
       
    end if
    

    
    
    return 
  end subroutine statistics_estimate_corr_factors

  !-----------------------------------------------------------------------
  
  subroutine statistics_estimate_vs_qs(self)
    class(statistics), intent(inout) :: self
    integer :: ierr, io, ios, il, iu, im
    character(line_max) :: vs_file, qs_file, out_file

    out_file = "uniform_structure.stat"
    
    write(vs_file, '(a,I2.2,a)')"vs.", self%rank, ".out"

    write(qs_file, '(a,I2.2,a)')"qs.", self%rank, ".out"
    

    call self%read_v_file(vs_file, self%vs)
    call self%read_v_file(qs_file, self%qs)
    
    il = 0.025 * self%n_mod
    im = 0.5 * self%n_mod
    iu = 0.975 * self%n_mod
    
    if (self%rank == 0) then
       
       call quick_sort(self%vs(:), 1, self%n_mod)
       call quick_sort(self%qs(:), 1, self%n_mod)
       open(newunit=io, file=out_file, status="replace", form="formatted", &
            & iostat=ios)
       write(io, '(A)') "# Vs (50%), Vs (2.5%) " // &
            & "Vs (97.5%), Qs (50%), Qs (2.5%), Qs (97.5%)" 
       write(io,'(6F11.6)') &
            & self%vs(im), &
            & self%vs(il), self%vs(iu), &
            & self%qs(im), &
            & self%qs(il), self%qs(iu) 
       close(io)
    end if


    return 
  end subroutine statistics_estimate_vs_qs

  !-----------------------------------------------------------------------
  

  subroutine statistics_read_v_file(self, v_file, v_all)
    class(statistics), intent(in) :: self
    character(line_max), intent(in) :: v_file
    double precision, intent(out) :: v_all(:)
    double precision :: v(self%n_mod)
    integer :: io, ios, ierr, i_mod, i_dummy, j_mod, i
    integer :: n_procs, i_rank, n_recv
    integer, allocatable :: ista(:)
    
    ! Read MCMC samples from file
    if (self%verb) then
       write(*,'(a)')"<< Now reading " // trim(v_file) // " >>"
    end if
    open(newunit=io, file=v_file, iostat=ios, status='old', &
         & access="stream", form="unformatted")
    if (ios /= 0) then
       write(0,*)"ERROR: cannot open ", trim(v_file)
       call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
    end if
    i_mod = 0
    do 
       read(io, iostat=ios)i_dummy, v(i_mod+1)
       if (ios /= 0) exit 
       i_mod = i_mod + 1
    end do
    close(io)
    
    
    ! Gather
    allocate(ista(MPI_STATUS_SIZE))
    if (self%rank == 0) then
       v_all(1:i_mod) = v(1:i_mod)
       j_mod = i_mod
       do i_rank = 1, self%n_procs - 1
          call mpi_recv(n_recv, 1, MPI_INTEGER, i_rank, 0, &
               & MPI_COMM_WORLD, ista, ierr)
          call mpi_recv(v_all(j_mod+1), n_recv, MPI_DOUBLE_PRECISION, &
               & i_rank, 1, MPI_COMM_WORLD, ista, ierr)
          j_mod = j_mod + n_recv
       end do
    else
       call mpi_send(i_mod, 1, MPI_INTEGER, 0, 0, &
            & MPI_COMM_WORLD, ista, ierr)
       call mpi_send(v(1), i_mod, MPI_DOUBLE_PRECISION, &
            & 0, 1, MPI_COMM_WORLD, ista, ierr)
    end if
    
    return 
  end subroutine statistics_read_v_file
    
    


  !-----------------------------------------------------------------------
  
  subroutine statistics_read_corr_file(self, corr_file, corr_all)
    class(statistics), intent(in) :: self
    character(line_max), intent(in) :: corr_file
    double precision, intent(out) :: corr_all(:,:)
    double precision :: corr(self%n_mod, self%n_sta)
    integer :: io, ios, ierr, i_mod, i_sta, i_dummy, j_mod, n_recv
    integer :: i_rank
    integer, allocatable :: ista(:)
    
    if (self%verb) then
       write(*,'(a)')"<< Now reading " // trim(corr_file) // " >>"
    end if
    open(newunit=io, file=corr_file, iostat=ios, status='old', &
         & access="stream", form="unformatted")
    if (ios /= 0) then
       write(0,*)"ERROR: cannot open ", trim(corr_file)
       call mpi_abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
    end if
    
    i_mod = 0
    do 
       read(io, iostat=ios)i_dummy, (corr(i_mod+1, i_sta), i_sta=1,self%n_sta)
       if (ios /= 0) exit 
       i_mod = i_mod + 1
    end do
    close(io)

    ! Gather
    allocate(ista(MPI_STATUS_SIZE))
    do i_sta = 1, self%n_sta
       if (self%rank == 0) then
          corr_all(1:i_mod, i_sta) = corr(1:i_mod, i_sta)
          j_mod = i_mod
          do i_rank = 1, self%n_procs - 1
             call mpi_recv(n_recv, 1, MPI_INTEGER, i_rank, 0, &
                  & MPI_COMM_WORLD, ista, ierr)
             call mpi_recv(corr_all(j_mod+1, i_sta), n_recv, &
                  & MPI_DOUBLE_PRECISION, &
                  & i_rank, 1, MPI_COMM_WORLD, ista, ierr)
             j_mod = j_mod + n_recv
          end do
       else
          call mpi_send(i_mod, 1, MPI_INTEGER, 0, 0, &
               & MPI_COMM_WORLD, ista, ierr)
          call mpi_send(corr(1, i_sta), i_mod, MPI_DOUBLE_PRECISION, &
               & 0, 1, MPI_COMM_WORLD, ista, ierr)
       end if
    end do
    


    return 
  end subroutine statistics_read_corr_file

  !-----------------------------------------------------------------------

end module cls_statistics
