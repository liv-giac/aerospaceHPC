PROGRAM TimeDependentSolver
    IMPLICIT NONE

    INTEGER :: n, i
    REAL(KIND=8) :: a, b, dt, dx, T0, T1, nt, error, t_dep, tt,start_time, end_time
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: T, sol
    REAL(KIND=8), PARAMETER :: Time = 1.d0
    REAL(KIND=8), PARAMETER :: X = 1.0

    WRITE(*, *) "Insert n:"
    ! Read the value of n from user input
    READ(*, *) n
        nt = 4.0 * n * n
        WRITE(*,*) "nt: ", nt

        error = 0.0

        dx = X / REAL(n + 1)
        dt = Time / REAL(nt)

        a = dt / (dx * dx)
        b = -2.d0 * a

        ! Allocate arrays
        ALLOCATE(T(1:n), sol(1:n))

        ! Initialize T array
        T = 0.0
        tt=0.0
        
        CALL CPU_TIME(start_time)
        DO WHILE (tt<Time)
            
            t_dep = SIN(REAL(tt)) + COS(REAL(tt))
            T0 = 0.d0
            T1 = sin(X) * sin(REAL(tt))
            

            DO i = 1, n
                sol(i) = T(i) + b * T(i) + dt * SIN(REAL(i) * dx) * t_dep

                IF (i == 1) THEN
                    sol(i) = sol(i) + a * (T0 + T(i + 1))
                END IF

                IF (i == n) THEN
                    sol(i) = sol(i) + a * (T(i - 1) + T1)
                END IF

                IF (i > 1 .AND. i < n) THEN
                    sol(i) = sol(i) + a * (T(i - 1) + T(i + 1))
                END IF
            END DO

            ! Update T array
            T = sol
            tt=tt+dt
        END DO
        CALL CPU_TIME(end_time)
        ! Exact Sol
        DO i = 1, n
            sol(i) = SIN(REAL(i * dx)) * SIN(Time)
        END DO
        ! Calculate error
        error = 0.0
        DO i = 1, n
            error = error + (sol(i) - T(i))**2
        END DO
        error = SQRT(error) / REAL(n)

        WRITE(*,*) "Error: ", error, ". Error / n:",error/REAL(n),"."
        WRITE(*,*)
        WRITE(*,*) "Time elapsed: ", (end_time - start_time),&
        ". Time elapsed scaled by n: ", (end_time - start_time)/ REAL(n),"."
        ! Deallocate arrays
        DEALLOCATE(T, sol)


END PROGRAM TimeDependentSolver
