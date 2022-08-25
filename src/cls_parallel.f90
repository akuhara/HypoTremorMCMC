module cls_parallel
  use cls_mcmc
  use mod_mpi
  use mod_random
  implicit none 

  type parallel
     private
     integer :: n_proc 
     integer :: rank   
     integer :: n_chain
     type(mcmc), allocatable :: mc(:) 
     logical :: verb = .false.
   contains
     procedure :: set_mc => parallel_set_mc
     procedure :: get_mc => parallel_get_mc
     procedure :: get_rank => parallel_get_rank
     procedure :: swap_temperature => parallel_swap_temperature
     procedure :: select_pair => parallel_select_pair
     procedure :: output_proposal => parallel_output_proposal

  end type parallel

  interface parallel
     module procedure init_parallel
  end interface parallel

  private pack_pair_info, unpack_pair_info, judge_swap
  private pack_mc_info, unpack_mc_info
  
contains
  
  !---------------------------------------------------------------------
  
  type(parallel) function init_parallel(n_proc, rank, n_chain, verb) &
       & result(self)
    integer, intent(in) :: n_proc, rank, n_chain
    logical, intent(in), optional :: verb

    if (present(verb)) then
       self%verb = verb
    end if

    if (self%verb) then
       write(*,'(A)')"<< Initialize prallel MCMC >>"
    end if

    self%n_proc = n_proc
    self%rank = rank
    self%n_chain = n_chain
    
    allocate(self%mc(n_chain))

    if (self%verb) then
       write(*,'(A,I5,A,I5,A,I5)')"# of total MCMC chains = ", &
            & n_proc, " proc. * ", n_chain, " chains/proc. = ", &
            & n_chain * n_proc
       write(*,*)
    end if
    
    return 
  end function init_parallel
  
 
  !---------------------------------------------------------------------

  subroutine parallel_set_mc(self, i, mc)
    class(parallel), intent(inout) :: self
    integer, intent(in) :: i
    type(mcmc), intent(in) :: mc
    
    self%mc(i) = mc
    
    return 
  end subroutine parallel_set_mc

  !---------------------------------------------------------------------

  type(mcmc) function parallel_get_mc(self, i) result(mc)
    class(parallel), intent(inout) :: self
    integer, intent(in) :: i
    
    mc = self%mc(i)

    return 
  end function parallel_get_mc
    
  !---------------------------------------------------------------------

  integer function parallel_get_rank(self) result(rank)
    class(parallel), intent(in) :: self

    rank = self%rank
    
    return 
  end function parallel_get_rank

  !---------------------------------------------------------------------

  subroutine parallel_swap_temperature(self, verb)
    class(parallel), intent(inout) :: self
    logical, intent(in), optional :: verb
    integer :: ipack(4), ierr, status(MPI_STATUS_SIZE)
    integer :: rank1, rank2, chain1, chain2
    double precision :: temp1, temp2, l1, l2, rpack(2)
    type(mcmc) :: mc1, mc2
    logical :: is_accepted

    if (self%rank == 0) then
       ipack = self%select_pair()
    end if
    call mpi_bcast(ipack, 4, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*)"ERROR: while MPI_BCAST"
       call mpi_finalize(ierr)
       stop
    end if
    call unpack_pair_info(ipack=ipack, rank1=rank1, rank2=rank2, &
         & chain1=chain1, chain2=chain2)
    
    if (self%rank == rank1 .and. self%rank == rank2) then
       mc1   = self%get_mc(chain1)
       temp1 = mc1%get_temp()
       l1    = mc1%get_log_likelihood()
       mc2   = self%get_mc(chain2)
       temp2 = mc2%get_temp()
       l2    = mc2%get_log_likelihood()
       
       call judge_swap(temp1, temp2, l1, l2, is_accepted)
       
       if (is_accepted) then
          call mc2%set_temp(temp1)
          call mc1%set_temp(temp2)
          call self%set_mc(chain1, mc1)
          call self%set_mc(chain2, mc2)
       end if
       
       !if (present(verb) .and. verb) then
       !   write(*,*)"----"
       !   write(*,*)"Rank1    :", rank1
       !   write(*,*)"Chain1   :", chain1
       !   write(*,*)"Temp1    : ", temp1
       !   write(*,*)"Rank2    :", rank2 
       !   write(*,*)"Chain2   :", chain2
       !   write(*,*)"Temp2    : ", temp2
       !   write(*,*)"Accepted :", is_accepted
       !   write(*,*)"----"
       !end if
       
    else if (self%rank == rank1) then
       mc1   = self%get_mc(chain1)
       temp1 = mc1%get_temp()
       l1    = mc1%get_log_likelihood()
       call mpi_recv(rpack, 2, MPI_DOUBLE_PRECISION, rank2, 1111, &
            & MPI_COMM_WORLD, status, ierr)
       if (ierr /= MPI_SUCCESS) then
          write(0,*)"ERROR: while MPI_RECV rpack 1"
          call mpi_finalize(ierr)
          stop
       end if
       call unpack_mc_info(rpack, temp=temp2, likelihood=l2)

       call judge_swap(temp1, temp2, l1, l2, is_accepted)
       
       if (is_accepted) then
          rpack = pack_mc_info(temp=temp1, likelihood=-999.d0)
          call mc1%set_temp(temp2)
          call self%set_mc(chain1, mc1)
       else
          rpack = pack_mc_info(temp=temp2, likelihood=-999.d0)
       end if
       call mpi_send(rpack, 2, MPI_DOUBLE_PRECISION, rank2, 2222, &
            & MPI_COMM_WORLD, status, ierr)
       if (ierr /= MPI_SUCCESS) then
          write(0,*)"ERROR: while MPI_SEND rpack 1"
          call mpi_finalize(ierr)
          stop
       end if
       !if (present(verb) .and. verb) then
       !   write(*,*)"----"
       !   write(*,*)"Rank1    :", rank1
       !   write(*,*)"Chain1   :", chain1
       !   write(*,*)"Temp1    : ", temp1
       !   write(*,*)"Rank2    :", rank2 
       !   write(*,*)"Chain2   :", chain2
       !   write(*,*)"Temp2    : ", temp2
       !   write(*,*)"Accepted :", is_accepted
       !   write(*,*)"----"
       !end if
    else if (self%rank == rank2) then
       ! Sender
       mc2   = self%get_mc(chain2)
       temp2 = mc2%get_temp()
       l2    = mc2%get_log_likelihood()
       rpack = pack_mc_info(temp=temp2, likelihood=l2)
       call mpi_send(rpack, 2, MPI_DOUBLE_PRECISION, rank1, 1111, &
            & MPI_COMM_WORLD, status, ierr)
       if (ierr /= MPI_SUCCESS) then
          write(0,*)"ERROR: while MPI_SEND rpack 2"
          call mpi_finalize(ierr)
          stop
       end if
       call mpi_recv(rpack, 2, MPI_DOUBLE_PRECISION, rank1, 2222, &
            & MPI_COMM_WORLD, status, ierr)
       if (ierr /= MPI_SUCCESS) then
          write(0,*)"ERROR: while MPI_RECV rpack 2"
          call mpi_finalize(ierr)
          stop
       end if
       call unpack_mc_info(rpack, temp=temp2, likelihood=l2)
       call mc2%set_temp(temp2)
       call self%set_mc(chain2, mc2)
    end if
    
    return 
  end subroutine parallel_swap_temperature
    
  !---------------------------------------------------------------------

  function parallel_select_pair(self) result(ipack)
    class(parallel), intent(in) :: self
    integer :: ipack(4)
    integer :: rank1, rank2, chain1, chain2
    integer :: i1, i2
    
    i1 = int(rand_u() * self%n_proc * self%n_chain)
    do 
       i2 = int(rand_u() * self%n_proc * self%n_chain)
       if (i1 /= i2) exit
    end do
    rank1 = int(i1 / self%n_chain)
    rank2 = int(i2 / self%n_chain)
    chain1 = mod(i1, self%n_chain) + 1
    chain2 = mod(i2, self%n_chain) + 1
    
    ipack = pack_pair_info(rank1=rank1, rank2=rank2, &
         & chain1=chain1, chain2=chain2)
    
    return 
  end function parallel_select_pair

  !---------------------------------------------------------------------

  subroutine parallel_output_proposal(self, filename, label)
    class(parallel), intent(inout) :: self
    character(*), intent(in) :: filename, label(:)
    type(mcmc) :: mc
    integer, allocatable :: n_accept_all(:), n_propose_all(:)
    integer, allocatable :: n_accept_sum(:), n_propose_sum(:)
    integer :: io, i, n, ierr
   
    n = size(label)
    allocate(n_accept_all(n), n_propose_all(n))
    allocate(n_accept_sum(n), n_propose_sum(n))
    n_accept_all(:) = 0
    n_propose_all(:) = 0

    ! Gather within the same node
    do i = 1, self%n_chain
       mc = self%get_mc(i)
       n_accept_all = n_accept_all + mc%get_n_accept()
       n_propose_all = n_propose_all + mc%get_n_propose()
    end do
    
    call mpi_reduce(n_accept_all, n_accept_sum, n, MPI_INTEGER4, MPI_SUM, &
         & 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(n_propose_all, n_propose_sum, n, MPI_INTEGER4, MPI_SUM, &
         & 0, MPI_COMM_WORLD, ierr)
    ! Output
    if (self%rank == 0) then
       open(newunit = io, file = filename, status = "unknown", &
            & iostat = ierr)
       do i = 1, n
          write(io, '(A,2I10)')'"' // label(i) // '"', &
               & n_propose_sum(i), n_accept_sum(i)
       end do
       close(io)
    end if
    
    return 
  end subroutine parallel_output_proposal
  
  !---------------------------------------------------------------------

  subroutine judge_swap(temp1, temp2, l1, l2, is_accepted)
    double precision, intent(in) :: temp1, temp2, l1, l2
    logical, intent(out) :: is_accepted
    double precision :: del_s
    double precision :: r
    double precision, parameter :: eps = epsilon(1.d0)
    
    del_s = (l2 - l1) * (1.d0 / temp1 - 1.d0 / temp2)
    is_accepted = .false.
    r = rand_u()
    if (r >= eps) then
       if(log(r) <= del_s) then
          is_accepted = .true.
       end if
    end if
    
    return 
  end subroutine judge_swap
  
  !---------------------------------------------------------------------

  function pack_pair_info(rank1, rank2, chain1, chain2) result(ipack)
    integer, intent(in) :: rank1, rank2, chain1, chain2
    integer :: ipack(4)
    
    ipack(1) = rank1
    ipack(2) = rank2
    ipack(3) = chain1
    ipack(4) = chain2
    
    return 
  end function pack_pair_info

  !---------------------------------------------------------------------

  subroutine unpack_pair_info(ipack, rank1, rank2, chain1, chain2)
    integer, intent(in) ::ipack(4)
    integer, intent(out) :: rank1, rank2, chain1, chain2
    
    rank1 = ipack(1)
    rank2 = ipack(2)
    chain1 = ipack(3)
    chain2 = ipack(4)
    
    return 
  end subroutine unpack_pair_info
  
  !---------------------------------------------------------------------

  function pack_mc_info(temp, likelihood) result(rpack)
    double precision, intent(in) :: temp, likelihood
    double precision :: rpack(2)
    
    rpack(1) = temp
    rpack(2) = likelihood

    return 
  end function pack_mc_info

  !---------------------------------------------------------------------

  subroutine unpack_mc_info(rpack, temp, likelihood)
    double precision, intent(in) :: rpack(2)
    double precision, intent(out) :: temp, likelihood
    
    temp = rpack(1)
    likelihood = rpack(2)

    return 
  end subroutine unpack_mc_info

  !---------------------------------------------------------------------


end module cls_parallel
