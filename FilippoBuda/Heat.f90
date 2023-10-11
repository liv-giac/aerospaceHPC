PROGRAM timedep
    IMPLICIT NONE

    INTEGER, PARAMETER :: nx = 10, nt = 400
    REAL(KIND = 8) :: len, dx, time, dt, timestep
    
    REAL(KIND = 8) :: sincos, dt_dx2

    REAL(KIND = 8) :: sum, error

    INTEGER :: i

    REAL(KIND = 8), DIMENSION(nx) :: told, tnew


    len = 1
    time = 1

    dx = len / (nx-1)
    dt = 1.d0 / nt

    WRITE(*,*) 'dx = ' ,dx
    WRITE(*,*) 'dt = ' ,dt 
    
    dt_dx2 = dt / (dx**2)

    told = 0
    tnew = 0
    timestep = 0

    DO WHILE (timestep .LE. time)
        timestep = timestep + dt

        tnew(1) = 0

        sincos = sin(timestep) + cos(timestep)

        DO i=2, nx-1
            tnew(i) = told(i) + dt_dx2 * (told(i-1) - 2 * told(i) + told(i+1)) + dt * sin( (i-1) * dx) * sincos
        END DO

        tnew(nx) = sin(1.d0) * sin(timestep)

        told = tnew
    END DO

    sum = 0

    DO i = 1, SIZE(told)
        sum = sum + ( (told(i) - sin((i-1)*dx)*sin(1.d0)) )**2
    END DO

    error = sqrt(sum) / nx
    
    WRITE(*,*)'error =', error
    ! WRITE(*,*) timestep, NORM2(Tnew-sin(len)*sin(timestep - dt))/nx

END PROGRAM



