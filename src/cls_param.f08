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
     integer :: n_cmps = 2
     
     ! Stations
     character(line_max) :: station_file
     character(line_max), allocatable :: stations(:)
     integer :: n_stations
     double precision, allocatable :: sta_x(:), sta_y(:), sta_z(:)
     double precision, allocatable :: sta_amp_fac(:,:)

     ! Data ID
     character(line_max) :: data_id_file
     integer :: n_data_id
     character(line_max), allocatable :: data_id(:)
     
     ! Components
     character(line_max) :: cmps(2)

     ! File format
     character(line_max) :: data_dir, filename_format
     character(line_max), allocatable :: filenames(:,:,:)

     ! Time window for converting envelope
     double precision :: t_win_conv

     ! Time window for correlation detection
     double precision :: t_win_corr
     double precision :: t_step_corr
     

     ! Detection criteria
     double precision :: alpha
     integer :: n_pair_thred

     ! Selection criteria 
     double precision :: vs_min, vs_max
     double precision :: b_min, b_max

     ! Event assumed during quality control
     double precision :: z_guess

     ! MCMC
     integer :: n_iter, n_burn, n_interval
     integer :: n_chains
     integer :: n_cool
     double precision :: temp_high
     double precision :: prior_width_xy, prior_width_z
     double precision :: prior_xy, prior_z
     double precision :: prior_width_vs, prior_vs
     double precision :: prior_width_qs, prior_qs
     double precision :: prior_t_corr = 0.d0
     double precision :: prior_width_t_corr
     double precision :: prior_a_corr = 0.d0
     double precision :: prior_width_a_corr
     double precision :: step_size_xy, step_size_z
     double precision :: step_size_vs, step_size_t_corr
     double precision :: step_size_qs, step_size_a_corr
     
     ! Inversion setting
     logical :: solve_vs
     logical :: solve_t_corr
     logical :: solve_qs
     logical :: solve_a_corr

     logical :: use_amp
     logical :: use_time

     ! Verbose
     logical :: verb = .false.
     
     ! For checking parameters 
     integer :: n_given
     character(line_max) :: given(100)
     character(line_max) :: from_where
     character(line_max), dimension(7) :: param_convert = &
          & [ character(line_max) :: &
          & "station_file", "data_dir", "data_id_file", &
          & "cmp1", "cmp2", "filename_format", "t_win_conv" ]

   contains
     procedure :: read_param_file => param_read_param_file
     procedure :: check_params => param_check_params
     procedure :: read_station_file => param_read_station_file
     procedure :: read_data_id_file => param_read_data_id_file
     procedure :: set_value  => param_set_value
     procedure :: get_n_cmps => param_get_n_cmps
     procedure :: get_n_stations => param_get_n_stations
     procedure :: get_stations => param_get_stations
     procedure :: get_station => param_get_station
     procedure :: get_sta_y => param_get_sta_y
     procedure :: get_sta_x => param_get_sta_x
     procedure :: get_sta_z => param_get_sta_z
     procedure :: get_sta_amp_fac => param_get_sta_amp_fac
     procedure :: get_n_data_id => param_get_n_data_id
     procedure :: get_data_id => param_get_data_id
     procedure :: get_cmps => param_get_cmps
     procedure :: get_data_dir => param_get_data_dir
     procedure :: get_filename_format => param_get_filename_format
     procedure :: make_filenames => param_make_filenames
     procedure :: get_filenames => param_get_filenames
     procedure :: get_t_win_conv => param_get_t_win_conv
     procedure :: get_t_win_corr => param_get_t_win_corr
     procedure :: get_t_step_corr => param_get_t_step_corr
     procedure :: get_n_pair_thred => param_get_n_pair_thred
     procedure :: get_alpha => param_get_alpha
     procedure :: get_vs_min => param_get_vs_min
     procedure :: get_vs_max => param_get_vs_max
     procedure :: get_b_min => param_get_b_min
     procedure :: get_b_max => param_get_b_max
     procedure :: get_z_guess => param_get_z_guess
     procedure :: get_n_iter => param_get_n_iter
     procedure :: get_n_burn => param_get_n_burn
     procedure :: get_n_interval => param_get_n_interval
     procedure :: get_n_chains => param_get_n_chains
     procedure :: get_n_cool => param_get_n_cool
     procedure :: get_temp_high => param_get_temp_high
     procedure :: get_prior_width_xy => param_get_prior_width_xy
     procedure :: get_prior_width_z => param_get_prior_width_z
     procedure :: get_prior_z => param_get_prior_z
     procedure :: get_prior_vs => param_get_prior_vs
     procedure :: get_prior_width_vs => param_get_prior_width_vs
     procedure :: get_prior_qs => param_get_prior_qs
     procedure :: get_prior_width_qs => param_get_prior_width_qs
     procedure :: get_prior_t_corr => param_get_prior_t_corr
     procedure :: get_prior_width_t_corr => param_get_prior_width_t_corr
     procedure :: get_prior_a_corr => param_get_prior_a_corr
     procedure :: get_prior_width_a_corr => param_get_prior_width_a_corr
     procedure :: get_step_size_z => param_get_step_size_z
     procedure :: get_step_size_xy => param_get_step_size_xy
     procedure :: get_step_size_vs => param_get_step_size_vs
     procedure :: get_step_size_qs => param_get_step_size_qs
     procedure :: get_step_size_t_corr => param_get_step_size_t_corr
     procedure :: get_step_size_a_corr => param_get_step_size_a_corr
     procedure :: get_solve_vs => param_get_solve_vs
     procedure :: get_solve_t_corr => param_get_solve_t_corr
     procedure :: get_solve_qs => param_get_solve_qs
     procedure :: get_solve_a_corr => param_get_solve_a_corr
     procedure :: get_use_amp => param_get_use_amp
     procedure :: get_use_time => param_get_use_time
     
  end type param
  
  interface param
     module procedure init_param
  end interface param

contains
  
  !---------------------------------------------------------------------
  
  type(param) function init_param(param_file, verb, from_where) result(self)
    character(len=*), intent(in) :: param_file
    logical, intent(in)      :: verb
    character(*), intent(in) :: from_where
    
    self%verb = verb
    
    self%from_where = from_where
    
    ! Read parmeter file
    if (self%verb) then
       write(*,'(3A)')"<< Reading parameters from ", &
            & trim(param_file), " >>"
    end if    
    self%param_file = param_file
    self%n_given = 0
    call self%read_param_file()
    if (self%verb) then
       write(*,*)
    end if
    
    ! Check parameters
    call self%check_params()

    ! Read station file
    if (self%verb) &
         & write(*,'(3A)')"<< Reading ", trim(self%station_file), " >>"
    call self%read_station_file()
    if (self%verb) write(*,*)
    
    if (from_where == "convert") then
       ! Read data ID file
       if (self%verb) &
            & write(*,'(3A)')"<< Reading ", trim(self%data_id_file), " >>"
       call self%read_data_id_file()
       if (self%verb) write(*,*)
    end if

    if (from_where == "convert") then
       call self%make_filenames()
    end if
    

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

  subroutine param_check_params(self)
    class(param), intent(inout) :: self
    integer :: i
    
    if (self%from_where == "convert") then
       
       do i = 1, size(self%param_convert)
          if (.not. any(self%given == self%param_convert(i))) then
             
             if (self%verb) write(*,*)"ERROR: ", trim(self%param_convert(i)), &
                  & " is not given."
             stop
          end if
       end do
    end if

    
    
    return 
  end subroutine param_check_params

  !---------------------------------------------------------------------

  subroutine param_read_station_file(self)
    class(param), intent(inout) :: self
    integer :: io, ierr, n_sta, i, j
    
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
    allocate(self%sta_y(n_sta))
    allocate(self%sta_x(n_sta)) 
    allocate(self%sta_z(n_sta))
    allocate(self%sta_amp_fac(self%n_cmps, n_sta))
    
    ! Read lines
    rewind(io)
    do i = 1, n_sta
       read(io, *) self%stations(i), &
            & self%sta_x(i), self%sta_y(i), self%sta_z(i), &
            & (self%sta_amp_fac(j,i), j= 1, self%n_cmps)
       if (self%verb) then
          write(*,'(1x,a,3F9.3)')trim(self%stations(i)), &
               & self%sta_x(i), self%sta_y(i), self%sta_z(i)
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
    else if (name == "cmp1") then
       self%cmps(1) = val
    else if (name == "cmp2") then
       self%cmps(2) = val
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
    else if (name == "alpha") then
       read(val,*) self%alpha
    else if (name == "vs_min") then
       read(val,*) self%vs_min
    else if (name == "vs_max") then
       read(val,*) self%vs_max
    else if (name == "b_min") then
       read(val,*) self%b_min
    else if (name == "b_max") then
       read(val,*) self%b_max
    else if (name == "z_guess") then
       read(val,*) self%z_guess
    else if (name == "n_iter") then
       read(val,*) self%n_iter
    else if (name == "n_burn") then
       read(val,*) self%n_burn
    else if (name == "n_interval") then
       read(val,*) self%n_interval
    else if (name == "n_chains") then
       read(val,*) self%n_chains
    else if (name == "n_cool") then
       read(val,*) self%n_cool
    else if (name == "temp_high") then
       read(val,*) self%temp_high
    else if (name == "prior_width_xy") then
       read(val,*) self%prior_width_xy
    else if (name == "prior_width_z") then
       read(val,*) self%prior_width_z
    else if (name == "prior_z") then
       read(val,*) self%prior_z
    else if (name == "prior_vs") then
       read(val,*) self%prior_vs
    else if (name == "prior_width_vs") then
       read(val,*) self%prior_width_vs
    else if (name == "prior_qs") then
       read(val,*) self%prior_qs
    else if (name == "prior_width_qs") then
       read(val,*) self%prior_width_qs
    else if (name == "prior_t_corr") then
       read(val,*) self%prior_t_corr
    else if (name == "prior_width_t_corr") then
       read(val,*) self%prior_width_t_corr
    else if (name == "prior_a_corr") then
       read(val,*) self%prior_a_corr
    else if (name == "prior_width_a_corr") then
       read(val,*) self%prior_width_a_corr
    else if (name == "step_size_xy") then
       read(val,*) self%step_size_xy
    else if (name == "step_size_z") then
       read(val,*) self%step_size_z
    else if (name == "step_size_vs") then
       read(val,*) self%step_size_vs
    else if (name == "step_size_t_corr") then
       read(val,*) self%step_size_t_corr
    else if (name == "step_size_qs") then
       read(val,*) self%step_size_qs
    else if (name == "step_size_a_corr") then
       read(val,*) self%step_size_a_corr
    else if (name == "solve_vs") then
       read(val,*) self%solve_vs
    else if (name == "solve_qs") then
       read(val,*) self%solve_qs
    else if (name == "solve_t_corr") then
       read(val,*) self%solve_t_corr
    else if (name == "solve_a_corr") then
       read(val,*) self%solve_a_corr
    else if (name == "use_amp") then
       read(val,*) self%use_amp
    else if (name == "use_time") then
       read(val,*) self%use_time
    else
       if (self%verb) then
          write(0,*)"ERROR: Invalid parameter name"
          write(0,*)"        : ", name, "  (?)"
       end if
       stop
    end if
    
    self%n_given = self%n_given + 1
    self%given(self%n_given) = name
    
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

  function param_get_sta_y(self) result(sta_y)
    class(param), intent(in) :: self
    double precision :: sta_y(self%n_stations)

    sta_y = self%sta_y
    
    return 
  end function param_get_sta_y

  !---------------------------------------------------------------------

  function param_get_sta_x(self) result(sta_x)
    class(param), intent(in) :: self
    double precision :: sta_x(self%n_stations)

    sta_x = self%sta_x
    
    return 
  end function param_get_sta_x

  !---------------------------------------------------------------------

  function param_get_sta_z(self) result(sta_z)
    class(param), intent(in) :: self
    double precision :: sta_z(self%n_stations)

    sta_z = self%sta_z
    
    return 
  end function param_get_sta_z

  !---------------------------------------------------------------------

  function param_get_sta_amp_fac(self, id) result(sta_amp_fac)
    class(param), intent(in) :: self
    integer, intent(in) :: id
    double precision :: sta_amp_fac(self%n_cmps)

    sta_amp_fac = self%sta_amp_fac(:, id)
    
    return 
  end function param_get_sta_amp_fac

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

  function param_get_cmps(self) result(cmps)
    class(param), intent(in) :: self
    character(line_max) :: cmps(self%n_cmps)

    cmps = self%cmps

    return 
  end function param_get_cmps

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
                  else if (trim(list(i_list)) == "$CMP") then
                     self%filenames(i_id, i_cmp, i_sta) = &
                          trim(self%filenames(i_id, i_cmp, i_sta)) // &
                          & trim(self%cmps(i_cmp))
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

  double precision function param_get_alpha(self) result(alpha)
    class(param), intent(in) :: self
    
    alpha = self%alpha

    return 
  end function param_get_alpha

  !---------------------------------------------------------------------
  
  double precision function param_get_vs_min(self) result(vs_min)
    class(param), intent(in) :: self
    
    vs_min = self%vs_min

    return 
  end function param_get_vs_min

  !---------------------------------------------------------------------

  double precision function param_get_vs_max(self) result(vs_max)
    class(param), intent(in) :: self
    
    vs_max = self%vs_max

    return 
  end function param_get_vs_max

  !---------------------------------------------------------------------

  double precision function param_get_b_min(self) result(b_min)
    class(param), intent(in) :: self
    
    b_min = self%b_min

    return 
  end function param_get_b_min

  !---------------------------------------------------------------------

  double precision function param_get_b_max(self) result(b_max)
    class(param), intent(in) :: self
    
    b_max = self%b_max

    return 
  end function param_get_b_max

  !---------------------------------------------------------------------

  double precision function param_get_z_guess(self) result(z_guess)
    class(param), intent(in) :: self
    
    z_guess = self%z_guess

    return 
  end function param_get_z_guess

  !---------------------------------------------------------------------

  integer function param_get_n_iter(self) result(n_iter)
    class(param), intent(in) :: self
    
    n_iter = self%n_iter

    return 
  end function param_get_n_iter

  !---------------------------------------------------------------------

  integer function param_get_n_interval(self) result(n_interval)
    class(param), intent(in) :: self
    
    n_interval = self%n_interval

    return 
  end function param_get_n_interval

  !---------------------------------------------------------------------

  integer function param_get_n_burn(self) result(n_burn)
    class(param), intent(in) :: self
    
    n_burn = self%n_burn

    return 
  end function param_get_n_burn

  !---------------------------------------------------------------------

  integer function param_get_n_chains(self) result(n_chains)
    class(param), intent(in) :: self
    
    n_chains = self%n_chains

    return 
  end function param_get_n_chains

  !---------------------------------------------------------------------
  
  integer function param_get_n_cool(self) result(n_cool)
    class(param), intent(in) :: self
    
    n_cool = self%n_cool

    return 
  end function param_get_n_cool

  !---------------------------------------------------------------------

  double precision function param_get_temp_high(self) result(temp_high)
    class(param), intent(in) :: self
    
    temp_high = self%temp_high

    return 
  end function param_get_temp_high
  
  !---------------------------------------------------------------------

  double precision function param_get_prior_width_xy(self) &
       & result(prior_width_xy)
    class(param), intent(in) :: self
    
    prior_width_xy = self%prior_width_xy

    return 
  end function param_get_prior_width_xy
  
  !---------------------------------------------------------------------
  
  double precision function param_get_prior_width_z(self) &
       & result(prior_width_z)
    class(param), intent(in) :: self
    
    prior_width_z = self%prior_width_z

    return 
  end function param_get_prior_width_z
  
  !---------------------------------------------------------------------
  
  double precision function param_get_prior_z(self) &
       & result(prior_z)
    class(param), intent(in) :: self
    
    prior_z = self%prior_z

    return 
  end function param_get_prior_z


  !---------------------------------------------------------------------

  double precision function param_get_prior_width_vs(self) &
       & result(prior_width_vs)
    class(param), intent(in) :: self
    
    prior_width_vs = self%prior_width_vs

    return 
  end function param_get_prior_width_vs
  
  !---------------------------------------------------------------------
  
  double precision function param_get_prior_vs(self) &
       & result(prior_vs)
    class(param), intent(in) :: self
    
    prior_vs = self%prior_vs

    return 
  end function param_get_prior_vs

  !---------------------------------------------------------------------

  double precision function param_get_prior_width_qs(self) &
       & result(prior_width_qs)
    class(param), intent(in) :: self
    
    prior_width_qs = self%prior_width_qs

    return 
  end function param_get_prior_width_qs
  
  !---------------------------------------------------------------------
  
  double precision function param_get_prior_qs(self) &
       & result(prior_qs)
    class(param), intent(in) :: self
    
    prior_qs = self%prior_qs

    return 
  end function param_get_prior_qs

  !---------------------------------------------------------------------

  double precision function param_get_prior_width_t_corr(self) &
       & result(prior_width_t_corr)
    class(param), intent(in) :: self
    
    prior_width_t_corr = self%prior_width_t_corr

    return 
  end function param_get_prior_width_t_corr
  
  !---------------------------------------------------------------------
  
  double precision function param_get_prior_t_corr(self) &
       & result(prior_t_corr)
    class(param), intent(in) :: self
    
    prior_t_corr = self%prior_t_corr

    return 
  end function param_get_prior_t_corr

  !---------------------------------------------------------------------

  double precision function param_get_prior_width_a_corr(self) &
       & result(prior_width_a_corr)
    class(param), intent(in) :: self
    
    prior_width_a_corr = self%prior_width_a_corr

    return 
  end function param_get_prior_width_a_corr
  
  !---------------------------------------------------------------------
  
  double precision function param_get_prior_a_corr(self) &
       & result(prior_a_corr)
    class(param), intent(in) :: self
    
    prior_a_corr = self%prior_a_corr

    return 
  end function param_get_prior_a_corr

  !---------------------------------------------------------------------
  
  double precision function param_get_step_size_xy(self) &
       & result(step_size_xy)
    class(param), intent(in) :: self
    
    step_size_xy = self%step_size_xy

    return 
  end function param_get_step_size_xy
  
  !---------------------------------------------------------------------

  double precision function param_get_step_size_z(self) &
       & result(step_size_z)
    class(param), intent(in) :: self
    
    step_size_z = self%step_size_z

    return 
  end function param_get_step_size_z


  !---------------------------------------------------------------------

    double precision function param_get_step_size_t_corr(self) &
       & result(step_size_t_corr)
    class(param), intent(in) :: self
    
    step_size_t_corr = self%step_size_t_corr

    return 
  end function param_get_step_size_t_corr
  
  !---------------------------------------------------------------------

    double precision function param_get_step_size_a_corr(self) &
       & result(step_size_a_corr)
    class(param), intent(in) :: self
    
    step_size_a_corr = self%step_size_a_corr

    return 
  end function param_get_step_size_a_corr
  
  !---------------------------------------------------------------------
  
  double precision function param_get_step_size_vs(self) &
       & result(step_size_vs)
    class(param), intent(in) :: self
    
    step_size_vs = self%step_size_vs

    return 
  end function param_get_step_size_vs

  !---------------------------------------------------------------------
  
  double precision function param_get_step_size_qs(self) &
       & result(step_size_qs)
    class(param), intent(in) :: self
    
    step_size_qs = self%step_size_qs

    return 
  end function param_get_step_size_qs
  
  !---------------------------------------------------------------------

  logical function param_get_solve_vs(self) &
       & result(solve_vs)
    class(param), intent(in) :: self
    
    solve_vs = self%solve_vs

    return 
  end function param_get_solve_vs
  
  !---------------------------------------------------------------------

  logical function param_get_solve_t_corr(self) &
       & result(solve_t_corr)
    class(param), intent(in) :: self
    
    solve_t_corr = self%solve_t_corr

    return 
  end function param_get_solve_t_corr
  
  !---------------------------------------------------------------------

  logical function param_get_solve_qs(self) &
       & result(solve_qs)
    class(param), intent(in) :: self
    
    solve_qs = self%solve_qs

    return 
  end function param_get_solve_qs
  
  !---------------------------------------------------------------------

  logical function param_get_solve_a_corr(self) &
       & result(solve_a_corr)
    class(param), intent(in) :: self
    
    solve_a_corr = self%solve_a_corr

    return 
  end function param_get_solve_a_corr
  
  !---------------------------------------------------------------------
  
  logical function param_get_use_amp(self) &
       & result(use_amp)
    class(param), intent(in) :: self
    
    use_amp = self%use_amp

    return 
  end function param_get_use_amp
  
  !---------------------------------------------------------------------

  logical function param_get_use_time(self) &
       & result(use_time)
    class(param), intent(in) :: self
    
    use_time = self%use_time

    return 
  end function param_get_use_time
  
  !---------------------------------------------------------------------

end module cls_param
