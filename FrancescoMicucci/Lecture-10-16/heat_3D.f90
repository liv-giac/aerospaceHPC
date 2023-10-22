PROGRAM Heat3D

    IMPLICIT NONE

    INTEGER, PARAMETER :: nx = 100, nt = 100, N = nx * nx * nx
    INTEGER :: i, j, k, index
    REAL(KIND=8) :: time = 0, dx, dt, diag, upperDiag, underDiag, coeff, w, b
    REAL(KIND=8) :: start, end, error
    REAL(KIND=8), DIMENSION(1:N) :: rhs, sol, previous_sol, exact_sol
    REAL(KIND=8), DIMENSION(nx) :: xsin
    
    dx = 1.d0 / (nx + 1)
    dt = 1.d0 / nt
    coeff = dt / (dx * dx)
    diag = 1.d0 + coeff
    upperDiag = -0.5d0 * coeff
    underDiag = -0.5d0 * coeff
    sol = 0.d0
    w = underDiag/diag
    b = diag - underDiag*upperDiag/diag

    DO i = 1, nx
        xsin(i) = sin(i * dx)
    END DO

    ! Set exact solution
    DO k = 1, nx
        DO j = 1, nx
            DO i = 1, nx
                index = i + (j - 1) * nx + (k - 1) * nx * nx
                exact_sol(index) = xsin(i)*xsin(j)*xsin(k)*sin(1.d0)
            END DO
        END DO
    END DO

    CALL cpu_time(start)

    DO WHILE (time < 1.d0)
        time = time + dt
        previous_sol = sol

        ! Set rhs
        DO k = 1, nx 
            DO j = 1, nx 
                DO i = 1, nx 
                    index = i + (j - 1) * nx + (k - 1) * nx * nx
                    rhs(index) = dt * (xsin(i)*xsin(j)*xsin(k)*(3.d0*sin(time+dt/2.d0)+cos(time+dt/2.d0)) &
                                - 6.d0/(dx*dx)*sol(index))
                    IF (i > 1) THEN
                        rhs(index) = rhs(index) + coeff * sol(i-1+(j-1)*nx+(k-1)*nx*nx)
                    END IF
                    IF (j > 1) THEN
                        rhs(index) = rhs(index) + coeff * sol(i+(j-2)*nx+(k-1)*nx*nx)
                    END IF
                    IF (k > 1) THEN
                        rhs(index) = rhs(index) + coeff * sol(i+(j-1)*nx+(k-2)*nx*nx)
                    END IF
                    IF (i < nx) THEN
                        rhs(index) = rhs(index) + coeff * sol(i+1+(j-1)*nx+(k-1)*nx*nx)
                    ELSE
                        rhs(index) = rhs(index) + coeff * (sin(1.d0)*xsin(j)*xsin(k)*sin(time))
                    END IF
                    IF (j < nx) THEN
                        rhs(index) = rhs(index) + coeff * sol(i+j*nx+(k-1)*nx*nx)
                    ELSE
                        rhs(index) = rhs(index) + coeff * (sin(1.d0)*xsin(i)*xsin(k)*sin(time))
                    END IF
                    IF (k < nx) THEN
                        rhs(index) = rhs(index) + coeff * sol(i+(j-1)*nx+k*nx*nx)
                    ELSE
                        rhs(index) = rhs(index) + coeff * (sin(1.d0)*xsin(i)*xsin(j)*sin(time)) 
                    END IF
                END DO      
            END DO
        END DO

        ! Thomas algorithm along x
        DO k = 1, nx 
            DO j = 1, nx 
                DO i = 2, nx 
                    index = i + (j - 1) * nx + (k - 1) * nx * nx
                    rhs(index) = rhs(index) - w * rhs(index - 1) 
                END DO
                index = nx + (j - 1) * nx + (k - 1) * nx * nx
                sol(index) = rhs(index) / b
                DO i = nx - 1, 1, -1
                    index = i + (j - 1) * nx + (k - 1) * nx * nx
                    sol(index) = (rhs(index) - upperDiag * sol(index + 1))/b    
                END DO
            END DO
        END DO

        ! Thomas algorithm along y
        rhs = sol
        DO k = 1, nx 
            DO i = 1, nx 
                DO j = 2, nx 
                    index = i + (j - 1) * nx + (k - 1) * nx * nx
                    if(j == nx) then
                        rhs(index) = rhs(index) - underDiag*sin(1.d0)*xsin(i)*xsin(k)*(sin(time+dt)-sin(time))
                    end if
                    rhs(index) = rhs(index) - w * rhs(index - nx) 
                END DO
                index = i + (nx - 1) * nx + (k - 1) * nx * nx
                sol(index) = rhs(index) / b
                DO j = nx - 1, 1, -1
                    index = i + (j - 1) * nx + (k - 1) * nx * nx
                    sol(index) = (rhs(index) - upperDiag * sol(index + nx))/b    
                END DO
            END DO
        END DO
        
        ! Thomas algorithm along z
        rhs = sol
        DO i = 1, nx 
            DO j = 1, nx 
                DO k = 2, nx 
                    index = i + (j - 1) * nx + (k - 1) * nx * nx
                    if(k == nx) then
                        rhs(index) = rhs(index) - underDiag*sin(1.d0)*xsin(j)*xsin(i)*(sin(time+dt)-sin(time))
                    end if
                    rhs(index) = rhs(index) - w * rhs(index - nx*nx) 
                END DO
                index = i + (j - 1) * nx + (nx - 1) * nx * nx
                sol(index) = rhs(index) / b
                DO k = nx - 1, 1, -1
                    index = i + (j - 1) * nx + (k - 1) * nx * nx
                    sol(index) = (rhs(index) - upperDiag * sol(index + nx * nx))/b    
                END DO
            END DO
        END DO
        
        sol = sol + previous_sol
        time = time + dt
    END DO
    CALL cpu_time(end)
    error = NORM2(sol - exact_sol)/N 
    
    WRITE(*,*) "Time:  ", end-start  
    WRITE(*,*) "Error: ", error

END PROGRAM Heat3D