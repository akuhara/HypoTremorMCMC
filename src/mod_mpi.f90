module mod_mpi
  implicit none
  include 'mpif.h'

  public get_mpi_task_id

contains
  
  subroutine get_mpi_task_id(n_all_task, id_start, id_end, debug)
    integer, intent(in) :: n_all_task
    integer, intent(out) :: id_start, id_end
    logical, optional, intent(in) :: debug
    integer :: ierr, n_procs, rank, n_task, i
    
    call mpi_comm_size(MPI_COMM_WORLD, n_procs, ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank,    ierr)
    
    n_task = n_all_task / n_procs
    
    id_start = n_task * rank + 1
    id_end   = id_start + n_task - 1
    id_start = id_start + min(rank,   mod(n_all_task, n_procs))
    id_end   = id_end   + min(rank+1, mod(n_all_task, n_procs))

    if (present(debug) .and. debug) then
       if (rank == 0) then
          print *, "#      Rank  :     start ID  --       end ID"
       end if
       call mpi_barrier(MPI_COMM_WORLD, ierr)

       do i = 0, n_procs - 1
          if (i == rank) then
             print *, rank, " : ", id_start, " -- ", id_end
          end if
          call mpi_barrier(MPI_COMM_WORLD, ierr)
       end do
    end if
    
    return 
  end subroutine get_mpi_task_id
  
end module mod_mpi
  
