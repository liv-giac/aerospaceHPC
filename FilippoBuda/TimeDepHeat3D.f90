PROGRAM jimdouglas
    IMPLICIT NONE
    INTEGER, PARAMETER :: nx = 30
    INTEGER, PARAMETER :: nt = 100
    REAL(KIND = 8), PARAMETER :: T = 1.0 ! time domain
    REAL(KIND = 8), PARAMETER :: X = 1.0 ! space domain
    REAL(KIND = 8) :: dx , dt, timestep = 0 ! intervals and current timestep
    REAL(KIND = 8) :: mat_coeff, dt_dx2 ! matrix A coefficient
    REAL(KIND = 8) :: error,start, finish
    REAL(KIND = 8) :: w = 0 ! real coefficient for thomas algorithm
    !	 a - sub-diagonal (means it is the diagonal below the main diagonal)
    !	 b - the main diagonal
    !	 c - sup-diagonal (means it is the diagonal above the main diagonal)
    !	 d - right part -> RHS
    !	 x - the answer
    !	 n - number of equations
    REAL(KIND = 8),DIMENSION(:), ALLOCATABLE :: a,b,c,d
    REAL(KIND = 8),DIMENSION(:), ALLOCATABLE :: f
    REAL(KIND = 8), DIMENSION(:), ALLOCATABLE :: told,tnew,real_solution

    INTEGER :: i,j,k,dim

    dx = X / ( 1 + nx)
    dt = T / nt
    dt_dx2 = dt / (dx**2)
    mat_coeff =  dt_dx2 / 2

    dim = nx**3

    allocate(a(1:dim)); allocate(b(1:dim)); allocate(c(1:dim)); allocate(d(1:dim));
    allocate(f(1:dim)); allocate(told(1:dim)); allocate(tnew(1:dim));
    allocate(real_solution(1:dim));
    
    a = - mat_coeff ! sub-diagonal
    b = 1 + 2 * mat_coeff ! diagonal
    c = - mat_coeff ! upper-diagonal

    ! w0_ijk = T (x_i, y_j, z_k, t= 0) = 0 for each i,j,k

    told = 0
    
    CALL CPU_TIME(start)

    timestep = timestep + dt
    DO WHILE(timestep .LE. T)
        ! initializing rhs 
        DO k=1, nx
            DO j=1, nx
                DO i=1, nx
                    d(index(i,j,k)) = dt * forcing_term(i,j,k,timestep + dt/2)
                    d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(-6 * told(index(i,j,k)))
                    ! i-1
                    IF(i .NE. 1) THEN
                        d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(told(index(i-1,j,k)))
                    END IF
                    ! i+1
                    IF(i .NE. nx) THEN
                        d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(told(index(i+1,j,k)))
                    ELSE
                        d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(exact_sol(i+1,j,k,timestep)) ! boundary condition
                    END IF
                    ! j-1
                    IF (j .NE. 1) THEN
                        d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(told(index(i,j-1,k)))
                    END IF
                    ! j+1
                    IF (j .NE. nx) THEN
                        d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(told(index(i,j+1,k)))
                    ELSE
                        d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(exact_sol(i,j+1,k,timestep)) ! boundary condition
                    END IF
                    ! k-1
                    IF (k .NE. 1) THEN
                        d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(told(index(i,j,k-1)))
                    END IF
                    ! k+1
                    IF (k .NE. nx) THEN
                        d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(told(index(i,j,k+1)))
                    ELSE 
                        d(index(i,j,k)) = d(index(i,j,k)) + dt_dx2*(exact_sol(i,j,k+1,timestep)) ! boundary condition
                    END IF
                END DO
            d(index(nx,j,k)) = d(index(nx,j,k)) + mat_coeff *(exact_sol(nx+1,j,k,timestep+dt) - exact_sol(nx+1,j,k,timestep)) ! boundary condition
            ! now you solve a nx x nx linear system inside the big linear system in order to consider only an x row at a time
            ! we do it inside the i cycle in order to have the j and k fixed
            DO i = 2, nx
                w = a(index(i,j,k)) / b(index(i-1,j,k))
                b(index(i,j,k)) = b(index(i,j,k)) - w * c(index(i-1,j,k))
                d(index(i,j,k)) = d(index(i,j,k)) - w * d(index(i-1,j,k))  
            END DO
            tnew(index(nx,j,k)) = d(index(nx,j,k)) / b(index(nx,j,k))
            DO i = nx-1, 1, -1
                tnew(index(i,j,k)) = (d(index(i,j,k)) - c(index(i,j,k)) * tnew(index(i+1,j,k))) / b(index(i,j,k))
            END DO
            ! end of thomas algorithm    
            END DO
        END DO

        ! THOMAS ON Y
        d = tnew
        ! instead of reassigning every value of told i can just slide trough it fixing different indexes
        b = 1 + 2 * mat_coeff ! reassigning b
        DO k=1,nx
            DO j=2,nx ! now i solve on j
                DO i=1,nx
                    IF (j .EQ. nx) THEN
                        d(index(i,j,k)) = d(index(i,j,k)) + mat_coeff * (exact_sol(i,nx+1,k,timestep+dt) - exact_sol(i,nx+1,k,timestep))
                    END IF
                    w = a(index(i,j,k)) / b(index(i,j-1,k))
                    b(index(i,j,k)) = b(index(i,j,k)) - w * c(index(i,j-1,k))
                    d(index(i,j,k)) = d(index(i,j,k)) - w * d(index(i,j-1,k))
                END DO
            END DO
        END DO
        DO k=1,nx
            DO j=nx-1,1,-1
                DO i=1,nx
                    tnew(index(i,nx,k)) = d(index(i,nx,k)) / b(index(i,nx,k))                    
                    tnew(index(i,j,k)) = (d(index(i,j,k)) - c(index(i,j,k)) * tnew(index(i,j+1,k))) / b(index(i,j,k))
                END DO
            END DO
        END DO
        
        ! THOMAS ON Z
        d = tnew
        b = 1 + 2 * mat_coeff
        DO k=2,nx
            DO j=1,nx
                DO i=1, nx
                    IF(k .EQ. nx) THEN
                        d(index(i,j,k)) = d(index(i,j,k)) + mat_coeff * (exact_sol(i,j,nx+1,timestep+dt) - exact_sol(i,j,nx+1,timestep))
                    END IF
                    w = a(index(i,j,k)) / b(index(i,j,k-1))
                    b(index(i,j,k)) = b(index(i,j,k)) - w * c(index(i,j,k-1))
                    d(index(i,j,k)) = d(index(i,j,k)) - w * d(index(i,j,k-1))
                END DO
            END DO
        END DO
        DO k=nx-1,1,-1
            DO j=1,nx
                DO i=1,nx
                    tnew(index(i,j,nx)) = d(index(i,j,nx)) / b(index(i,j,nx))
                    tnew(index(i,j,k)) = (d(index(i,j,k)) - c(index(i,j,k)) * tnew(index(i,j,k+1))) / b(index(i,j,k))
                END DO
            END DO
        END DO
        b = 1 + 2 * mat_coeff

        told = tnew + told
        timestep = timestep + dt
    END DO

    CALL CPU_TIME(finish)

    PRINT '("Time = ",f15.10," seconds.")',(finish-start)

    do k = 1, Nx
        do j = 1, Nx
            do i = 1, Nx
                real_solution(index(i, j, k)) = exact_sol(i,j,k,1.d0)
            end do
        end do
    end do

    error = NORM2(told - real_solution)/dim
    
    WRITE(*,*) "Error= ",error
    


    CONTAINS
        INTEGER FUNCTION index(i,j,k)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i,j,k

            index = i + (j-1)*nx + (k-1)*nx*nx

        END FUNCTION index

        REAL FUNCTION exact_sol(i,j,k,t)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i,j,k
            REAL(KIND = 8), INTENT(IN) :: t

            exact_sol = sin(i*dx) * sin(j*dx) * sin(k*dx) * sin(t)

        END FUNCTION exact_sol

        REAL FUNCTION forcing_term(i,j,k,t)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i,j,k
            REAL(KIND = 8), INTENT(IN) :: t

            forcing_term = (cos(t) + 3*sin(t))*sin(i*dx)*sin(j*dx)*sin(k*dx)

        END FUNCTION forcing_term
END PROGRAM
