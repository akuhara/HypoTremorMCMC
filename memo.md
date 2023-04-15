```
SUBROUTINE linear_regression(N, x, y, w, a, b, err)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN) :: x(N), y(N), w(N)
  REAL*8, INTENT(OUT) :: a, b, err
  REAL*8 :: sumx, sumy, sumw, sumxy, sumx2, d

  sumx = 0
  sumy = 0
  sumw = 0
  sumxy = 0
  sumx2 = 0

  DO i=1,N
     sumx = sumx + x(i) * w(i)
     sumy = sumy + y(i) * w(i)
     sumw = sumw + w(i)
     sumxy = sumxy + x(i) * y(i) * w(i)
     sumx2 = sumx2 + x(i) * x(i) * w(i)
  END DO

  d = sumw * sumx2 - sumx * sumx

  a = (sumw * sumxy - sumx * sumy) / d
  b = (sumx2 * sumy - sumx * sumxy) / d

  err = 0
  DO i=1,N
     err = err + w(i) * (y(i) - a * x(i) - b) * (y(i) - a * x(i) - b)
  END DO
  err = sqrt(err / (N - 2))

END SUBROUTINE linear_regression
```

This subroutine takes the following inputs:

N: the length of the data vectors x, y, and w
x: the vector of x-values
y: the vector of y-values
w: the vector of measurement errors
a: a variable to store the slope of the linear fit
b: a variable to store the intercept of the linear fit
err: a variable to store the standard error of the fit
The subroutine uses the weighted least squares method to fit a line to the data. The weights are given by the vector w. The slope and intercept of the line are stored in the variables a and b, respectively. The standard error of the fit is stored in the variable err.

To use this subroutine, you can call it from your main program like this:

```
INTEGER :: N
REAL*8 :: x(N), y(N), w(N), a, b, err

! read in data from file or generate it

CALL linear_regression(N, x, y, w, a, b, err)

! do something with the results
```

Note that this code assumes that the measurement errors are normally distributed and that the errors in the y-values are independent. If these assumptions are not met, a different method may be more appropriate.
