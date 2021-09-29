!=======================================================================
!   Hypo_tremor_MCMC
!   Locating tremors via MCMC
!   Copyright (C) 2021 Takeshi Akuhara
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
!   Contact information
!
!   Email  : akuhara @ eri. u-tokyo. ac. jp 
!   Address: Earthquake Research Institute, The Univesity of Tokyo
!           1-1-1, Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan
!
!=======================================================================
!
! module cls_c3_data
!   This module defines a class 'c3_data' that handles three-component
!   seismogram data. Different instance should be used for each station.
!
!   USAGE:
!     
!     type(c3_data) :: c3
!     c3 = c3_data(['hoge.BHE', 'hoge.BHN', 'hoge.BHZ']) 
!     call c3%enqueue_from_files(['hoge2.BHE', 'hoge2.BHN', 'hoge3.BHZ'])
!     call c3%output_data('hoge.txt')
!=======================================================================
module cls_c3_data
  implicit none 
  
  type c3_data
     private
     integer                       :: n_smp 
     integer                       :: n_cmps
     double precision              :: dt
     double precision, allocatable :: data(:,:) ! size(n_smp, n_cmps)

     double precision              :: fac_memory_alloc = 4.d0
    
   contains
     procedure :: read_sac     => c3_data_read_sac
     procedure :: enqueue_from_files => c3_data_enqueue_from_files
     procedure :: enqueue_data => c3_data_enqueue_data
     procedure :: dequeue_data => c3_data_dequeue_data
     procedure :: output_data  => c3_data_output_data
     procedure :: decimate_data => c3_data_decimate_data
     procedure :: get_n_smp    => c3_data_get_n_smp
     procedure :: get_n_cmps   => c3_data_get_n_cmps
     procedure :: get_dt       => c3_data_get_dt
     procedure :: get_data     => c3_data_get_data
     procedure :: extract_data => c3_data_extract_data
     procedure :: calc_mean     => c3_data_calc_mean
     procedure :: calc_histogram => c3_data_calc_histogram
  end type c3_data

  interface c3_data
     module procedure init_c3_data
  end interface c3_data

contains
  
  !---------------------------------------------------------------------
  
  type(c3_data) function init_c3_data(files, n_cmps, dt) &
       & result(self)
    character(*), intent(in), optional :: files(:)
    integer, intent(in), optional :: n_cmps
    double precision, intent(in), optional :: dt
    
    if (present(n_cmps)) then
       self%n_cmps = n_cmps
    end if
    if (present(dt)) then
       self%dt = dt
    end if
    if (present(files)) then
       self%n_cmps = size(files)
       call self%read_sac(files)
       return 
    end if

    if (.not. present(files) .and. &
         & (.not. present(n_cmps) .or. .not. present(dt))) then
       error stop "ERROR: cannot initialize c3_data"
    end if

    

    return 
  end function init_c3_data
  
  !---------------------------------------------------------------------

  subroutine c3_data_read_sac(self, files)
    class(c3_data), intent(inout) :: self
    character(*), intent(in) :: files(:)
    integer :: io(3), ierr, npts, i, prev_npts, j
    real :: delta_4, tmp_4
    double precision, allocatable :: tmp(:, :)
    integer :: n

    n = self%n_cmps
    
    ! Open file
    do i = 1, n
       write(*,'(A,1X,A,1X,A)') "<< Reading ", trim(files(i)), ">>"
       open(newunit=io(i), file=files(i), status="old", iostat=ierr, &
            & access="direct", recl=4)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot open: ", trim(files(i))
          stop
       end if
    end do
    
    ! Read headers
    do i = 1, n
       ! Get dt
       read(io(i), rec=1) delta_4
       if (.not. allocated(self%data)) then
          self%dt = dble(delta_4)
       else 
          if (abs(self%dt -dble(delta_4)) > 1.e-6) then
             write(0,*)"ERROR: error in SAC header delta"
             stop
          end if
       end if
       
       ! Get n_smp
       read(io(i), rec=80) npts
       if (.not. allocated(self%data)) then
          ! Initial
          allocate(self%data(npts, n))
          prev_npts = npts
          self%n_smp = 0
       else if (i == 1) then
          ! Enqueue
          allocate(tmp(self%n_smp + npts, n))
          tmp(1:self%n_smp, 1:n) = self%data(1:self%n_smp, 1:n)
          call move_alloc(from=tmp, to=self%data)
          prev_npts = npts
       else
          if (prev_npts /= npts) then
             write(0,*)"ERROR: invalid npts"
             stop
          end if
       end if
    end do

    ! Read data
    do i = 1, n
       do j = 1, npts 
          read(io(i), rec=158+j) tmp_4
          self%data(self%n_smp+j, i) = dble(tmp_4)
       end do
    end do
    self%n_smp = self%n_smp + npts
    
    do i = 1, n
       close(io(i))
    end do
    
    
    return 
  end subroutine c3_data_read_sac

  !---------------------------------------------------------------------

  subroutine c3_data_enqueue_from_files(self, files)
    class(c3_data), intent(inout) :: self
    character(*), intent(in) :: files(:)

    call self%read_sac(files)

    return 
  end subroutine c3_data_enqueue_from_files
    
  !---------------------------------------------------------------------
  
  subroutine c3_data_enqueue_data(self, x)
    class(c3_data), intent(inout) :: self
    double precision, intent(in) :: x(:,:)
    double precision, allocatable :: tmp(:,:)    
    integer :: n, n_space, n_current, n_new

    !write(*,*)"size = ", size(x(1,:)), "n_cmps=", self%n_cmps
    if (size(x(1,:)) /= self%n_cmps) then
       error stop "ERROR: illegal component number in queued data"
    end if
    
    n = size(x(:,1))
    n_current = self%n_smp
    if (allocated(self%data)) then
       n_space = size(self%data(:,1))
    else
       n_space = 0
    end if
    
    if (n + n_current <= n_space) then
       
       ! Append data without memory allocation
       self%data(self%n_smp+1:self%n_smp+n,1:self%n_cmps) = &
            & x(1:n,1:self%n_cmps)

       self%n_smp = self%n_smp + n
    else if (n_space > 0) then

       ! Append data with memory allocation
       
       ! * Increase space
       n_new = max(self%n_smp + n, int(self%fac_memory_alloc * self%n_smp))
       allocate(tmp(1:n_new,1:self%n_cmps))
       ! * Copy stored data
       tmp(1:self%n_smp,1:self%n_cmps) = &
            & self%data(1:self%n_smp,1:self%n_cmps)
       ! * Append new data
       tmp(self%n_smp+1:self%n_smp+n, 1:self%n_cmps) =&
            & x(1:n, 1:self%n_cmps)
       call move_alloc(from=tmp, to=self%data)
       self%n_smp = self%n_smp + n
    else
       
       ! ---Initial memory allocation ---
       allocate(self%data(1:n,1:self%n_cmps))
       self%data(1:n, 1:self%n_cmps) = x(1:n, 1:self%n_cmps)
       self%n_smp = n
    end if

    return 
  end subroutine c3_data_enqueue_data

  !---------------------------------------------------------------------
  
  function c3_data_dequeue_data(self, n) result(data_out)
    class(c3_data), intent(inout) :: self
    integer, intent(in) :: n
    double precision :: data_out(n, self%n_cmps)
    double precision :: tmp(self%n_smp - n,self%n_cmps)

    if (n > self%n_smp) then
       error stop "ERROR: data length is not enough in queue"
    end if

    data_out(1:n, 1:self%n_cmps) = self%data(1:n, 1:self%n_cmps)
    tmp(1:self%n_smp - n, 1:self%n_cmps) = &
         & self%data(n + 1:self%n_smp, 1:self%n_cmps)
    self%data  = tmp
    self%n_smp = self%n_smp - n
    
    return 
  end function c3_data_dequeue_data
  
  !---------------------------------------------------------------------
  
  subroutine c3_data_output_data(self, out_file)
    class(c3_data), intent(inout) :: self
    character(*), intent(in) :: out_file
    integer :: io, ierr, i, j

    open(newunit=io, file=out_file, status="replace", iostat=ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot create ", trim(out_file)
       stop
    end if

    do i = 1, self%n_cmps
       do j = 1, self%n_smp
          write(io,*)(j-1) * self%dt, self%data(j, i)
       end do
       write(io,*)
       write(io,*)
    end do

    close(io)
    return 
  end subroutine c3_data_output_data

  !---------------------------------------------------------------------

  subroutine c3_data_decimate_data(self, n_fac) 
    class(c3_data), intent(inout) :: self
    integer, intent(in) :: n_fac
    double precision :: tmp(self%n_smp / n_fac, self%n_cmps)
    integer :: i, j
    j = 1
    do i = 1, self%n_smp
       if (mod(i, n_fac) == 1) then
          tmp(j, 1:self%n_cmps) = self%data(i, 1:self%n_cmps)
          j = j +1
       end if
    end do

    deallocate(self%data)
    allocate(self%data(self%n_smp / n_fac, self%n_cmps))
    self%data = tmp
    self%n_smp = self%n_smp / n_fac
    self%dt = self%dt * n_fac
    
    return 
  end subroutine c3_data_decimate_data
  
  !---------------------------------------------------------------------
  
  integer function c3_data_get_n_smp(self) result(n_smp)
    class(c3_data), intent(in) :: self
    
    n_smp = self%n_smp
    
    return 
  end function c3_data_get_n_smp

  !---------------------------------------------------------------------

  integer function c3_data_get_n_cmps(self) result(n_cmps)
    class(c3_data), intent(in) :: self
    
    n_cmps = self%n_cmps
    
    return 
  end function c3_data_get_n_cmps

  !---------------------------------------------------------------------

  double precision function c3_data_get_dt(self) result(dt)
    class(c3_data), intent(in) :: self
    
    dt = self%dt
    
    return 
  end function c3_data_get_dt

  !---------------------------------------------------------------------

  function c3_data_get_data(self) result(data)
    class(c3_data), intent(in) :: self
    double precision :: data(self%n_smp, self%n_cmps)
    
    data = self%data(1:self%n_smp, 1:self%n_cmps)
    
    return 
  end function c3_data_get_data

  !---------------------------------------------------------------------

  function c3_data_extract_data(self, i1, i2, k) result(data)
    class(c3_data), intent(in) :: self
    integer, intent(in) :: i1, i2, k
    double precision :: data(i2 - i1 + 1)

    data = self%data(i1:i2, k)
    
    return 
  end function c3_data_extract_data

  !---------------------------------------------------------------------

  function c3_data_calc_mean(self) result(mean)
    class(c3_data), intent(inout) :: self
    double precision :: mean(self%n_cmps)
    integer :: i_cmp

    do concurrent (i_cmp = 1:self%n_cmps)
       mean(i_cmp) = sum(self%data(1:self%n_smp,i_cmp))
       mean(i_cmp) = mean(i_cmp) / self%n_smp
    end do
    
    return 
  end function c3_data_calc_mean

  !---------------------------------------------------------------------

  function c3_data_calc_histogram(self, n) result(histo)
    class(c3_data), intent(inout) :: self
    integer, intent(in) :: n
    integer :: histo(n, self%n_cmps)
    double precision :: v_min, v_max, w
    integer :: i, j, idx
    
    v_min = minval(self%data(:,:))
    v_max = maxval(self%data(:,:)) * 1.001d0
    w = (v_max - v_min) / n
    histo(1:n, 1:self%n_cmps) = 0
    
    do concurrent(i = 1:self%n_cmps)
       do concurrent (j = 1:self%n_smp)
          idx = int((self%data(j,i) - v_min) / w) + 1
          histo(idx, i) = histo(idx, i) + 1
       end do
    end do
    

    return 
  end function c3_data_calc_histogram
  
  !---------------------------------------------------------------------

end module cls_c3_data
    
