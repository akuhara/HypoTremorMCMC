!=======================================================================
!   Hypo_tremor_MCMC
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
module cls_param
  use cls_line_text
  implicit none 
  
  type param
     private
     character(line_max) :: param_file
     
     ! Components
     integer :: n_cmps
     
     ! Stations
     character(line_max) :: station_file
     character(line_max), allocatable :: stations(:)
     integer :: n_stations
     double precision, allocatable :: x(:), y(:), z(:)

     ! Data ID
     character(line_max) :: data_id_file
     integer :: n_data_id
     character(line_max), allocatable :: data_id(:)
     
     ! Components
     character(line_max) :: comp_file
     character(line_max), allocatable :: comps(:)

     ! File format
     character(line_max) :: data_dir, filename_format
     character(line_max), allocatable :: filenames(:,:,:)

     ! Time window for converting envelope
     double precision :: t_win_conv

     ! Time window for correlation detection
     double precision :: t_win_corr
     double precision :: t_step_corr


     ! Detection criteria
     integer :: n_pair_thred

     ! Verbose
     logical :: verb = .false.

     

   contains
     procedure :: read_param_file => param_read_param_file
     procedure :: read_station_file => param_read_station_file
     procedure :: read_data_id_file => param_read_data_id_file
     procedure :: read_comp_file => param_read_comp_file
     procedure :: set_value  => param_set_value
     procedure :: get_n_cmps => param_get_n_cmps
     procedure :: get_n_stations => param_get_n_stations
     procedure :: get_stations => param_get_stations
     procedure :: get_station => param_get_station
     procedure :: get_x => param_get_x
     procedure :: get_y => param_get_y
     procedure :: get_n_data_id => param_get_n_data_id
     procedure :: get_data_id => param_get_data_id
     procedure :: get_comps => param_get_comps
     procedure :: get_data_dir => param_get_data_dir
     procedure :: get_filename_format => param_get_filename_format
     procedure :: make_filenames => param_make_filenames
     procedure :: get_filenames => param_get_filenames
     procedure :: get_t_win_conv => param_get_t_win_conv
     procedure :: get_t_win_corr => param_get_t_win_corr
     procedure :: get_t_step_corr => param_get_t_step_corr
     procedure :: get_n_pair_thred => param_get_n_pair_thred
  end type param
  
  interface param
     module procedure init_param
  end interface param

contains
  
  !---------------------------------------------------------------------
  
  type(param) function init_param(param_file, verb) result(self)
    character(len=*), intent(in) :: param_file
    logical, intent(in), optional :: verb

    if (present(verb)) then
       self%verb = verb
    end if
    
    ! Read parmeter file
    if (self%verb) then
       write(*,'(3A)')"<< Reading parameters from ", &
            & trim(param_file), " >>"
    end if    
    self%param_file = param_file
    call self%read_param_file()
    if (self%verb) then
       write(*,*)
    end if

    ! Read station file
    if (self%verb) &
         & write(*,'(3A)')"<< Reading ", trim(self%station_file), " >>"
    call self%read_station_file()
    if (self%verb) write(*,*)

    ! Read data ID file
    if (self%verb) &
         & write(*,'(3A)')"<< Reading ", trim(self%data_id_file), " >>"
    call self%read_data_id_file()
    if (self%verb) write(*,*)

    ! Read component file
    if (self%verb) &
         & write(*,'(3A)')"<< Reading ", trim(self%comp_file), " >>"
    call self%read_comp_file()
    if (self%verb) write(*,*)

    ! Make filenames
    call self%make_filenames()
    !write(*,*)self%filenames

    return 
  end function init_param

  !---------------------------------------------------------------------
    
  subroutine param_read_param_file(self)
    class(param), intent(inout) :: self
    character(len=line_max) :: line, name, val
    integer :: ierr, io
    type(line_text) :: lt
    logical :: is_ok

    open(newunit = io, file = self%param_file, &
         & status = 'old', iostat = ierr)
    if (ierr /= 0) then
       if (self%verb) then
          write(0,*) "ERROR: cannot open ", trim(self%param_file)
       end if
       stop
    end if
    
    do 
       read(io, '(a)', iostat=ierr) line
       if (ierr /= 0) then
          exit
       end if
       lt = init_line_text(line)
       call lt%read_value(name, val, is_ok)
       if (is_ok) then
          call self%set_value(name, val)
       end if
    end do
    close(io)
    
    return 
  end subroutine param_read_param_file

  !---------------------------------------------------------------------

  subroutine param_read_station_file(self)
    class(param), intent(inout) :: self
    integer :: io, ierr, n_sta, i
    
    open(newunit=io, file=self%station_file, status="old", iostat=ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open", trim(self%station_file)
       stop
    end if
    
    ! Count lines
    n_sta = 0
    do
       read(io, *, iostat=ierr)
       if (ierr /= 0) exit
       n_sta = n_sta + 1
    end do

    ! Allocate
    self%n_stations = n_sta
    allocate(self%stations(n_sta))
    allocate(self%x(n_sta))
    allocate(self%y(n_sta)) 
    allocate(self%z(n_sta))
    ! Read lines
    rewind(io)
    do i = 1, n_sta
       read(io, *) self%stations(i), self%x(i), self%y(i), self%z(i)
       if (self%verb) then
          write(*,'(1x,a,3F9.3)')trim(self%stations(i)), &
               & self%x(i), self%y(i), self%z(i)
       end if
    end do
    close(io)

    return 
  end subroutine param_read_station_file

  !---------------------------------------------------------------------

  subroutine param_read_data_id_file(self)
    class(param), intent(inout) :: self
    integer :: io, ierr, n_id, i
    
    open(newunit=io, file=self%data_id_file, status="old", iostat=ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open", trim(self%data_id_file)
       stop
    end if
    
    ! Count lines
    n_id = 0
    do
       read(io, *, iostat=ierr)
       if (ierr /= 0) exit
       n_id = n_id + 1
    end do
    
    ! Allocate
    self%n_data_id = n_id
    allocate(self%data_id(n_id))

    ! Read data ID
    rewind(io)
    do i = 1, n_id
       read(io, '(a)')self%data_id(i)
       if (self%verb) write(*,'(1x,a)') trim(self%data_id(i))
    end do
    close(io)
    
    return 
  end subroutine param_read_data_id_file

  !---------------------------------------------------------------------  
  
  subroutine param_read_comp_file(self)
    class(param), intent(inout) :: self
    integer :: io, ierr, i
    

    open(newunit=io, file=self%comp_file, status="old", iostat=ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open", trim(self%comp_file)
       stop
    end if
    
    self%n_cmps = 0
    do 
       read(io, *, iostat=ierr)
       if (ierr /= 0) exit
       self%n_cmps = self%n_cmps + 1
    end do
    
    rewind(io)
    allocate(self%comps(self%n_cmps))
    
    do i = 1, self%n_cmps
       read(io, '(a)') self%comps(i)
       if (self%verb) write(*,'(1x,a)') trim(self%comps(i))
    end do
    close(io)
    
    return 
  end subroutine param_read_comp_file

  !---------------------------------------------------------------------

  subroutine param_set_value(self, name, val)
    class(param), intent(inout) :: self
    character(len=*), intent(in) :: name, val
    
    if (self%verb) then
       write(*,*)trim(name), " <- ", trim(val)
    end if

    if (name == "station_file") then
       self%station_file = val
    else if (name == "data_id_file") then
       self%data_id_file = val
    else if (name == "comp_file") then
       self%comp_file = val
    else if (name == "data_dir") then
       self%data_dir = val
    else if (name == "filename_format") then
       self%filename_format = val
    else if (name == "t_win_conv") then
       read(val,*) self%t_win_conv
    else if (name == "t_win_corr") then
       read(val,*) self%t_win_corr
    else if (name == "t_step_corr") then
       read(val,*) self%t_step_corr
    else if (name == "n_pair_thred") then
       read(val,*) self%n_pair_thred
    else
       if (self%verb) then
          write(0,*)"ERROR: Invalid parameter name"
          write(0,*)"        : ", name, "  (?)"
       end if
       stop
    end if
    return 
  end subroutine param_set_value

  !---------------------------------------------------------------------

  integer function param_get_n_cmps(self) result(n_cmps)
    class(param), intent(in) :: self
    
    n_cmps = self%n_cmps

    return 
  end function param_get_n_cmps

  !---------------------------------------------------------------------

  integer function param_get_n_stations(self) result(n_stations)
    class(param), intent(in) :: self
    
    n_stations = self%n_stations

    return 
  end function param_get_n_stations

  !---------------------------------------------------------------------

  function param_get_stations(self) result(stations)
    class(param), intent(in) :: self
    character(line_max) :: stations(self%n_stations)

    stations = self%stations
    
    return 
  end function param_get_stations

  !---------------------------------------------------------------------
  
  character(line_max) function param_get_station(self, id) result(station)
    class(param), intent(in) :: self
    integer, intent(in) :: id
    
    if (id < 1 .or. id > self%n_stations) then
       error stop "ERROR: invalid station ID is given to param_get_station"
    end if

    station = self%stations(id)
    
    return 
  end function param_get_station

  !---------------------------------------------------------------------

  function param_get_x(self) result(x)
    class(param), intent(in) :: self
    double precision :: x(self%n_stations)

    x = self%x
    
    return 
  end function param_get_x

  !---------------------------------------------------------------------

  function param_get_y(self) result(y)
    class(param), intent(in) :: self
    double precision :: y(self%n_stations)

    y = self%y
    
    return 
  end function param_get_y

  !---------------------------------------------------------------------
  
  integer function param_get_n_data_id(self) result(n_data_id)
    class(param), intent(in) :: self
    
    n_data_id = self%n_data_id

    return 
  end function param_get_n_data_id
  
  !---------------------------------------------------------------------

  function param_get_data_id(self) result(data_id)
    class(param), intent(in) :: self
    character(line_max) :: data_id(self%n_data_id)
    
    data_id = self%data_id
    
    return 
  end function param_get_data_id

  !---------------------------------------------------------------------

  function param_get_comps(self) result(comps)
    class(param), intent(in) :: self
    character(line_max) :: comps(self%n_cmps)

    comps = self%comps

    return 
  end function param_get_comps

  !---------------------------------------------------------------------
  
  character(line_max) function param_get_data_dir(self) result(data_dir)
    class(param), intent(in) :: self
    
    data_dir = self%data_dir
    
    return 
  end function param_get_data_dir
  
  !---------------------------------------------------------------------

  character(line_max) function param_get_filename_format(self) &
       & result(filename_format)
    class(param), intent(in) :: self
    
    filename_format = self%filename_format
    
    return 
  end function param_get_filename_format
  
  !---------------------------------------------------------------------
  
  subroutine param_make_filenames(self)
    class(param), intent(inout) :: self
    type(line_text) :: lt
    character(line_max), allocatable :: list(:)
    
    allocate(self%filenames(self%n_data_id, self%n_cmps, self%n_stations))
    lt = line_text(self%filename_format, separator="+")
    call lt%make_list(list)

    block
      integer :: i_sta, i_cmp, i_id, i_list
      do i_sta = 1, self%n_stations
         do i_cmp = 1, self%n_cmps
            do i_id = 1, self%n_data_id
               self%filenames(i_id, i_cmp, i_sta) = &
                    & trim(self%data_dir) // "/"
               do i_list = 1, size(list)
                  if (trim(list(i_list)) == "$ID") then
                     self%filenames(i_id, i_cmp, i_sta) = &
                          trim(self%filenames(i_id, i_cmp, i_sta)) // &
                          & trim(self%data_id(i_id))
                  else if (trim(list(i_list)) == "$STA") then
                     self%filenames(i_id, i_cmp, i_sta) = &
                          trim(self%filenames(i_id, i_cmp, i_sta)) // &
                          & trim(self%stations(i_sta))
                  else if (trim(list(i_list)) == "$COMP") then
                     self%filenames(i_id, i_cmp, i_sta) = &
                          trim(self%filenames(i_id, i_cmp, i_sta)) // &
                          & trim(self%comps(i_cmp))
                  else 
                     self%filenames(i_id, i_cmp, i_sta) = &
                          trim(self%filenames(i_id, i_cmp, i_sta)) // &
                          trim(list(i_list))
                  end if
               end do
            end do
         end do
      end do
    end block

    return 
  end subroutine param_make_filenames
  
  !---------------------------------------------------------------------
  
  function param_get_filenames(self, i_sta) result(filenames)
    class(param), intent(inout) :: self
    integer, intent(in) :: i_sta
    character(line_max) :: filenames(self%n_data_id, self%n_cmps)
    
    filenames = self%filenames(:,:,i_sta)
    
    return 
  end function param_get_filenames

  !---------------------------------------------------------------------
  
  double precision function param_get_t_win_conv(self) result(t_win_conv)
    class(param), intent(in) :: self
    
    t_win_conv = self%t_win_conv

    return 
  end function param_get_t_win_conv

  !---------------------------------------------------------------------

  double precision function param_get_t_win_corr(self) result(t_win_corr)
    class(param), intent(in) :: self
    
    t_win_corr = self%t_win_corr

    return 
  end function param_get_t_win_corr

  !---------------------------------------------------------------------
  
  double precision function param_get_t_step_corr(self) result(t_step_corr)
    class(param), intent(in) :: self
    
    t_step_corr = self%t_step_corr

    return 
  end function param_get_t_step_corr

  !---------------------------------------------------------------------

  integer function param_get_n_pair_thred(self) result(n_pair_thred)
    class(param), intent(in) :: self
    
    n_pair_thred = self%n_pair_thred

    return 
  end function param_get_n_pair_thred

  !---------------------------------------------------------------------
  

end module cls_param
