program hypo_tremor_mcmc
  use cls_param
  use cls_line_text, only: line_max
  use cls_convertor
  implicit none 
  integer :: n_args
  character(line_max) :: param_file
  type(param) :: para
  type(convertor) :: conv
  logical :: verb
  
  verb = .true.
  
  ! Get argument
  n_args = command_argument_count()
  if (n_args /= 1) then
     error stop "USAGE: hypo_tremor_mcmc [parameter file]"
  end if
  call get_command_argument(1, param_file)

  ! Read parameter file
  para = param(param_file, verb=verb)

  ! Convert raw data to CF
  conv = convertor(para)
  
  stop
end program hypo_tremor_mcmc

