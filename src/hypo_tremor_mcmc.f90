program main
  use mod_mpi
  use mod_random
  use cls_model, only: model
  use cls_obs_data, only: obs_data
  use cls_param, only:param
  use cls_line_text, only: line_max
  use cls_mcmc, only: mcmc
  use cls_parallel, only: parallel
  use cls_forward, only: forward
  use, intrinsic :: iso_fortran_env, only: iostat_end
  
  implicit none
  integer :: n_args, ierr, rank, n_procs
  type(model) :: t_corr, hypo, vs
  type(model) :: t_corr_tmp, hypo_tmp, vs_tmp
  type(forward) :: fwd
  type(param) :: para
  type(obs_data) :: obs
  type(mcmc) :: mc
  type(parallel) :: pt
  
  character(line_max) :: param_file, win_id_file
  character(line_max) :: hypo_file, t_corr_file, vs_file
  logical :: verb, prior_ok
  integer, allocatable :: win_id(:)
  integer :: i, j, io, id, n_events, io_hypo, io_t_corr, io_vs
  double precision, allocatable :: x_mu(:), y_mu(:)
  double precision :: dummy
  double precision :: temp, log_prior_ratio, log_likelihood
  double precision, parameter :: eps = 1.0d-8


  !------------------------------------------------------------------------
  ! 1. Initialize
  !------------------------------------------------------------------------
  
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
  para = param(param_file, verb=verb)
  call mpi_barrier(MPI_COMM_WORLD, ierr)


  !
  call init_random(5551111, 453222, 4444431,6765, rank)
  
  ! Events
  win_id = [ integer :: ]
  win_id_file = "detected_win.dat"
  open(newunit=io, file=win_id_file, status="old", iostat=ierr)
  if (ierr /= 0) then
     error stop "cannot open detected_win.dat"
  end if
  do 
     read(io, *, iostat=ierr) id, dummy
     if (ierr == iostat_end) exit 
     win_id = [win_id, id]
  end do
  n_events = size(win_id)
  close(io)
  
  ! Observed data
  obs = obs_data( &
       & win_id = win_id, &
       & n_sta  = para%get_n_stations(), &
       & sta_x  = para%get_sta_x(), &
       & sta_y  = para%get_sta_y(), &
       & verb   = verb)
  allocate(x_mu(n_events))
  allocate(y_mu(n_events))
  call obs%make_initial_guess(x_mu, y_mu)
  
  ! forward
  fwd = forward(&
       & n_sta          = para%get_n_stations(), &
       & n_events       = n_events, &
       & sta_x          = para%get_sta_x(), &
       & sta_y          = para%get_sta_y(), &
       & sta_z          = para%get_sta_z(), &
       & obs            = obs,               &
       & use_laplace    = para%get_use_laplace(), &
       & use_covariance = para%get_use_covariance(), &
       & forward_diff   = para%get_forward_diff() &
       & )
       

  ! Initialize parallel MCMC chains
  pt = parallel( &
       & n_proc  = n_procs, &
       & rank    = rank,    &
       & n_chain = para%get_n_chains(), &
       & verb    = verb)
  
  do j = 1, para%get_n_chains()
     ! MCMC model parameters
     ! ** Station correction
     t_corr = model(&
          & nx   = para%get_n_stations(), &
          & verb = verb)
     if (para%get_solve_t_corr()) then
        do i = 1, para%get_n_stations()
           call t_corr%set_prior(i=i, mu=para%get_prior_t_corr(), &
                & sigma=para%get_prior_width_t_corr()) 
           call t_corr%set_perturb(i=i, step_size=para%get_step_size_t_corr())
        end do
        call t_corr%generate_model()
     else
        do i = 1, para%get_n_stations()
           call t_corr%set_x(i, para%get_prior_t_corr())
        end do
     end if
     if (verb) call t_corr%display()
     
     ! ** Hypocenters
     hypo = model(&
          & nx   = 3 * n_events, &
          & verb = verb)
     do i = 1, n_events
        call hypo%set_prior(i=3*i-2, mu=x_mu(i), sigma=para%get_prior_width_xy())
        call hypo%set_prior(i=3*i-1, mu=y_mu(i), sigma=para%get_prior_width_xy())
        call hypo%set_prior(i=3*i  , mu=0.d0,  sigma=para%get_prior_z(), &
             & prior_type=1)
        call hypo%set_perturb(i=3*i-2, step_size=para%get_step_size_xy())
        call hypo%set_perturb(i=3*i-1, step_size=para%get_step_size_xy())
        call hypo%set_perturb(i=3*i  , step_size=para%get_step_size_z())
     end do
     call hypo%generate_model()
     if (verb) call hypo%display()
     
     ! ** Velocity 
     vs = model(nx = 1, verb = verb)
     call vs%set_prior(i=1, mu=para%get_prior_vs(), sigma=para%get_prior_width_vs())
     call vs%set_perturb(i=1, step_size=para%get_step_size_vs())
     call vs%set_x(1, para%get_prior_vs())
     if (verb) call vs%display()
     
     mc = mcmc(&
          & hypo         = hypo,   &
          & t_corr       = t_corr, &
          & vs           = vs,     &
          & n_iter       = para%get_n_iter(), &
          & solve_t_corr = para%get_solve_t_corr(), &
          & solve_vs     = para%get_solve_vs() &
          & )
     
     ! Set temperatures
     if (j <= para%get_n_cool()) then
        call mc%set_temp(1.d0)
     else
        temp = exp((rand_u() *(1.d0 - eps) + eps) &
             & * log(para%get_temp_high()))
        call mc%set_temp(temp)
     end if
     call pt%set_mc(j, mc)

  end do
  
  !------------------------------------------------------------------------
  ! 2. MCMC
  !------------------------------------------------------------------------
  write(hypo_file, '(A,I2.2,A)')"hypo.", rank, ".out"
  write(t_corr_file, '(A,I2.2,A)')"t_corr.", rank, ".out"
  write(vs_file, '(A,I2.2,A)')"vs.", rank, ".out"
  open(newunit=io_hypo, file=hypo_file, status="replace", &
       & access="stream", form="unformatted")
  open(newunit=io_t_corr, file=t_corr_file, status="replace", &
       & access="stream", form="unformatted")
  open(newunit=io_vs, file=vs_file, status="replace", &
       & access="stream", form="unformatted")
!  open(newunit=io_vs, file=vs_file, status="replace", form="formatted")
  print *, "start MCMC"
  do i = 1, para%get_n_iter()
     do j = 1, para%get_n_chains()
        mc = pt%get_mc(j)

        ! Proposal
        call mc%propose_model(hypo_tmp, t_corr_tmp, vs_tmp, &
             & log_prior_ratio, prior_ok)

        ! Forward 
        if (prior_ok) then
           call fwd%calc_log_likelihood(hypo_tmp, t_corr_tmp, vs_tmp, &
                & log_likelihood)
        end if

        ! Judge
        call mc%judge_model(hypo_tmp, t_corr_tmp, vs_tmp, &
             & log_likelihood, log_prior_ratio, prior_ok)
        call pt%set_mc(j, mc)

        if (verb .and. j== 1) call mc%one_step_summary()
        
        ! Recording
        if (mc%get_temp() < 1.d0 + eps .and. i > para%get_n_burn() &
             & .and. mod(i, para%get_n_interval()) == 1) then
           write(io_vs)i,  mc%write_out_vs()
           write(io_hypo)i, mc%write_out_hypo()
           write(io_t_corr)i, mc%write_out_t_corr()
        end if
     end do
     ! Swap Temp.
     call pt%swap_temperature(verb=verb)
  end do

  close(io_hypo)
  close(io_t_corr)
  close(io_vs)

  call mpi_barrier(MPI_COMM_WORLD, ierr)
  call mpi_finalize(ierr)

  stop
end program main
