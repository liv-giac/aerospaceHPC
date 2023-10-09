PROGRAM heateq
IMPLICIT NONE

INTEGER :: i, n, nt
REAL(KIND=8) :: time, dx, dt, a, b, error
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: T, To, x
READ(*,*) n
ALLOCATE (T(0:n+1),To(0:n+1),x(0:n+1))
nt = 4 * n * n
dx = 1.d0 / (n+1)
time = 0.d0
dt = 1.d0 /nt
a= dt/(dx*dx)
b = -2.d0 * a
error = 0.d0

To = 0.d0

DO i=0,n+1
    x(i)=i*dx
END DO
DO WHILE (time < 1.d0)
    time = time + dt
    T(0) = 0.d0
    T(1:n)=To(1:n)*(1.d0 + b) + a * (To(0:n-1)+To(2:n+1))+dt*sin(x(1:n))*(sin(time)+cos(time))

    T(n+1) = sin(1.d0)*sin(time)
    To=T

END DO
error=NORM2(T-sin(x)*sin(time))/n
WRITE(*,*) time, error
END PROGRAM heateq
