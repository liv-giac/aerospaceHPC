PROGRAM heateq
    IMPLICIT NONE
    
    INTEGER :: i, n
    REAL(KIND=8) :: dx, a, b, error,w,start_time, end_time
    REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: T, x, rhs, diag
    
    WRITE(*, *) "Insert n:"
    ! Read the value of n from user input
    READ(*, *) n
    
    ! Allocate arrays with indices from 1 to n
    ALLOCATE(T(1:n), x(1:n), rhs(1:n), diag(1:n))
    
    dx = 1.d0 / (n+1)
    a = -1.d0 / (dx*dx)
    b = 2.d0 * a
    error = 0.d0
    
    ! Initialize arrays
    DO i = 1, n
        x(i) = dx *i
        diag(i) = -b
    END DO
    
    ! Initialize rhs
    rhs = SIN(x)

    rhs(n) = rhs(n) - SIN(1.d0) * a
    CALL CPU_TIME(start_time)

    ! Solve tridiagonal system
    DO i = 2, n
        w = a/diag(i-1)
        diag(i) = diag(i) - a * w 
        rhs(i) = rhs(i) - w * rhs(i-1)
    END DO
    
    ! Back substitution
    T(n) = rhs(n) / diag(n)
    DO i = n-1, 1, -1
        T(i) = (rhs(i) - a * T(i+1)) / diag(i)
    END DO
    
    CALL CPU_TIME(end_time)
    ! Calculate error
    DO i = 1, n
        error = error + (T(i) - SIN(x(i)))**2
    END DO
    error = SQRT(error)/ REAL(n)
    
    ! Print the error
    WRITE(*, *) "Error: ", error, ". Error /n: ",error/REAL(n),". Time elapsed: ", end_time - start_time,". Time elapsed scaled by n: ",(end_time-start_time)/REAL(n)
    DEALLOCATE(T, diag, rhs,x)
END PROGRAM heateq
