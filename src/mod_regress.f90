module mod_regress
  implicit none
  
contains
  subroutine linear_regression(n, x, y, w, a, b, err)
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), y(n), w(n)
    double precision, intent(out) :: a, b, err
    integer :: i
    double precision :: sumx, sumy, sumw, sumxy, sumx2, d
    
    sumx = 0.d0
    sumy = 0.d0
    sumw = 0.d0
    sumxy = 0.d0
    sumx2 = 0.d0
    
    do i = 1, n
       sumx = sumx + x(i) * w(i)
       sumy = sumy + y(i) * w(i)
       sumw = sumw + w(i)
       sumxy = sumxy + x(i) * y(i) * w(i)
       sumx2 = sumx2 + x(i) * x(i) * w(i)
    end do
    
    d = sumw * sumx2 - sumx * sumx
    
    a = (sumw * sumxy - sumx * sumy) / d
    b = (sumx2 * sumy - sumx * sumxy) / d
    
    err = 0
    do i = 1, n
       err = err + w(i) * (y(i) - a * x(i) - b) * (y(i) - a * x(i) - b)
    end do
    err = sqrt(err / (n - 2))
    
    return 
  end subroutine linear_regression
  
  subroutine weighted_corr(n, x, y, w, cc)
    integer, intent(in) :: n
    double precision, intent(in) :: x(n), y(n), w(n)
    double precision, intent(out) :: cc
    double precision :: mean_x, mean_y, sum_w, s_xx, s_yy, s_xy
    
    
    sum_w = sum(w(1:n))
    mean_x = sum(x(1:n) * w(1:n)) / sum_w
    mean_y = sum(y(1:n) * w(1:n)) / sum_w
    
    s_xx = sum((x(1:n)-mean_x)**2)
    s_yy = sum((y(1:n)-mean_y)**2)
    s_xy = sum((x(1:n)-mean_x) * (y(1:n)-mean_y))

    cc = s_xy / sqrt(s_xx * s_yy)
    
    return 
  end subroutine weighted_corr

end module mod_regress
