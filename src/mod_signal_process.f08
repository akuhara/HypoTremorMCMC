module mod_signal_process
  implicit none
  private
  double precision, parameter :: pi = acos(-1.d0)
  
  public apply_taper
contains
  !---------------------------------------------------------------------

  function apply_taper(n, x) result(x_out)
    integer, intent(in) :: n
    double precision, intent(in) :: x(1:n)
    double precision :: x_out(1:n)
    double precision, parameter :: pcnt = 0.05d0
    integer :: nleng, i
    double precision :: fac


    nleng = int(n * pcnt)
    x_out = x
    do concurrent (i = 1:nleng)
       fac = 0.5d0 * (1.d0 - cos((i - 1) * pi / nleng))
       x_out(i) = x(i) * fac
       x_out(n - i + 1) = x(n - i + 1) * fac
    end do

    return
  end function apply_taper

  !---------------------------------------------------------------------  
  
end module mod_signal_process
