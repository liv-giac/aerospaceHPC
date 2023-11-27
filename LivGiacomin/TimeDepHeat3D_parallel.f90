PROGRAM jimdouglas
    USE decomp_2d
    USE mpi
    IMPLICIT NONE
    INTEGER, PARAMETER :: nx = 30
    INTEGER, PARAMETER :: nt = 1
    INTEGER,DIMENSION(2) :: coords_proc
    REAL(mytype), PARAMETER :: T = 1.0 ! time domain
    REAL(mytype), PARAMETER :: X = 1.0 ! space domain
    REAL(mytype) :: d_x , dt, timestep = 0 ! intervals and current timestep
    REAL(mytype) :: mat_coeff, dt_dx2 ! matrix A coefficient
    REAL(mytype) :: error,start, finish
    REAL(mytype) :: w = 0 ! real coefficient for thomas algorithm
    !	 a - sub-diagonal (means it is the diagonal below the main diagonal)
    !	 b - the main diagonal
    !	 c - sup-diagonal (means it is the diagonal above the main diagonal)
    !	 d - right part -> RHS
    !	 x - the answer
    !	 n - number of equations
    REAL(mytype),DIMENSION(:), ALLOCATABLE :: a,b,c
    REAL(mytype),DIMENSION(:,:,:), ALLOCATABLE :: dx,dy,dz,told,real_solution
    REAL(mytype),DIMENSION(:,:), ALLOCATABLE :: toldn, tolde, toldw, tolds, btoldew,btoldns
    !REAL(mytype),DIMENSION(:), ALLOCATABLE :: f

    INTEGER process_Rank, size_Of_Cluster, ierror, tag
    INTEGER :: i,j,k,dim,P_row,P_col,dest_proc,source_proc

    CALL MPI_INIT(ierror)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, process_Rank, ierror)

    
    P_col=FLOOR(SQRT(dble(size_Of_Cluster)))
    DO WHILE (MOD(size_Of_Cluster,P_col) .NE. 0)
        P_col=P_col-1
    END DO
    P_row=size_Of_Cluster/P_col

    CALL decomp_2d_init(nx,nx,nx,P_row,P_col)
    IF (process_Rank .EQ. 0) WRITE (*,*) P_row, P_col

    CALL MPI_CART_COORDS(DECOMP_2D_COMM_CART_Z,process_Rank,2,coords_proc,ierror)

    CALL MPI_Cart_shift(DECOMP_2D_COMM_CART_Z, 0, 1, source_proc,dest_proc,ierror)


    WRITE (*,*) process_Rank, coords_proc, source_proc, dest_proc



    d_x = X / ( 1 + nx)
    dt = T / nt
    dt_dx2 = dt / (d_x**2)
    mat_coeff =  dt_dx2 / 2

    dim = nx**3

    allocate(a(1:nx)); allocate(b(1:nx)); allocate(c(1:nx)); 
    
    allocate(told(zsize(1),zsize(2),zsize(3)));

    allocate(dx(xsize(1),xsize(2),xsize(3)));
    allocate(dy(ysize(1),ysize(2),ysize(3)));
    allocate(dz(zsize(1),zsize(2),zsize(3)));
    allocate(tolde(zsize(2),zsize(3)))
    allocate(toldw(zsize(2),zsize(3)))
    allocate(btoldew(zsize(2),zsize(3)))
    allocate(toldn(zsize(1),zsize(3)))
    allocate(tolds(zsize(1),zsize(3)))
    allocate(btoldns(zsize(1),zsize(3)))
    
    


    allocate(real_solution(zsize(1),zsize(2),zsize(3)));
    
    a = - mat_coeff ! sub-diagonal
    b = 1 + 2 * mat_coeff ! diagonal
    c = - mat_coeff ! upper-diagonal

    ! w0_ijk = T (x_i, y_j, z_k, t= 0) = 0 for each i,j,k

    told = 0
    tolde = 0
    tolds = 0
    toldn = 0
    toldw = 0
    btoldew = 0
    btoldns = 0
    CALL CPU_TIME(start)

    timestep = timestep + dt
    DO WHILE(timestep .LE. T)
        ! initializing rhs 
        DO k=zstart(3), zend(3)
            DO j=zstart(2),zend(2)
                DO i=zstart(1),zend(1)
                    dz(i,j,k) = dt * forcing_term(i,j,k,timestep + dt/2)
                    dz(i,j,k) = dz(i,j,k) + dt_dx2*(-6 * told(i,j,k))
                    ! i-1
                    IF(i .NE. 1) THEN
                        dz(i,j,k) = dz(i,j,k) + dt_dx2*(toldw(j,k)) !nb
                    END IF
                    ! i+1
                    IF(i .NE. nx) THEN
                        dz(i,j,k) = dz(i,j,k) + dt_dx2*(tolde(j,k)) !nb
                    ELSE
                        dz(i,j,k) = dz(i,j,k) + dt_dx2*(exact_sol(i+1,j,k,timestep)) ! boundary condition
                    END IF
                    ! j-1
                    IF (j .NE. 1) THEN
                        dz(i,j,k) = dz(i,j,k) + dt_dx2*(toldn(i,k)) !nb
                    END IF
                    ! j+1
                    IF (j .NE. nx) THEN
                        dz(i,j,k) = dz(i,j,k) + dt_dx2*(tolds(i,k)) !nb
                    ELSE
                        dz(i,j,k) = dz(i,j,k) + dt_dx2*(exact_sol(i,j+1,k,timestep)) ! boundary condition
                    END IF
                    ! k-1
                    IF (k .NE. 1) THEN
                        dz(i,j,k) = dz(i,j,k) + dt_dx2*(told(i,j,k-1))
                    END IF
                    ! k+1
                    IF (k .NE. nx) THEN
                        dz(i,j,k) = dz(i,j,k) + dt_dx2*(told(i,j,k+1))
                    ELSE 
                        dz(i,j,k) = dz(i,j,k) + dt_dx2*(exact_sol(i,j,k+1,timestep)) ! boundary condition
                    END IF
            
                END DO
            dz(nx,j,k) = dz(nx,j,k) + mat_coeff *(exact_sol(nx+1,j,k,timestep+dt) - exact_sol(nx+1,j,k,timestep)) ! boundary condition

            END DO
        END DO

        !THOMAS ON Z
        CALL threed_thomas(a,b,c,dz)


        !!rotate
        ! THOMAS ON Y

        CALL transpose_z_to_y(dz, dy)
        ! instead of reassigning every value of told i can just slide trough it fixing different indexes
        b = 1 + 2 * mat_coeff ! reassigning b
        CALL threed_thomas(a,b,c,dy)

        ! THOMAS ON Z
        CALL transpose_y_to_x(dy, dx)
        b = 1 + 2 * mat_coeff
        CALL threed_thomas(a,b,c,dx)
        b = 1 + 2 * mat_coeff
        CALL transpose_x_to_y(dx, dy) !NOT OPTIMAL
        CALL transpose_y_to_z(dy, dz)
        told = dz + told
        btoldew = told(zend(1),zstart(2):zend(2),zstart(3):zend(3))
        CALL MPI_SEND(btoldew, size(btoldew,1)*size(btoldew,2), mytype, dest_proc, dest_proc, DECOMP_2D_COMM_CART_Z,ierror)
        CALL MPI_RECV(btoldew, size(btoldew,1)*size(btoldew,2), mytype, source_proc, source_proc, DECOMP_2D_COMM_CART_Z,ierror)
        timestep = timestep + dt
        CALL 
    END DO

    CALL CPU_TIME(finish)

    PRINT '("Time = ",f15.10," seconds.")',(finish-start)

    do k = zstart(3),zend(3)
        do j = zstart(2),zend(2)
            do i = zstart(1), zend(1)
                real_solution(i, j, k) = exact_sol(i,j,k,timestep)
            end do
        end do
    end do

    error = NORM2(told - real_solution)/dim
    
    WRITE(*,*) "Error= ",error
    
    CALL decomp_2d_finalize()
    CALL MPI_FINALIZE(ierror)

    CONTAINS

        SUBROUTINE threed_thomas(a,b,c,d)
            IMPLICIT NONE
            REAL(mytype),DIMENSION(:) :: a,b,c
            REAL(mytype),DIMENSION(:,:,:) :: d

            DO k=2,nx
                        !IF(k .EQ. nx) THEN
                        !    d(:,:,k) = d(:,:,k) + mat_coeff * (exact_sol(i,j,nx+1,timestep+dt) - exact_sol(i,j,nx+1,timestep))
                        !END IF
                        w = a(k) / b(k-1)
                        b(k) = b(k) - w * c(k-1)
                        d(:,:,k) = d(:,:,k) - w * d(:,:,k-1)
             
            END DO
    
            d(:,:,nx) = d(:,:,nx) / b(nx)
            DO k=nx-1,1,-1      
                        d(:,:,k) = (d(:,:,k) - c(k) * d(:,:,k+1)) / b(k)
            END DO
            

        END SUBROUTINE threed_thomas



         FUNCTION exact_sol(i,j,k,t) RESULT (exact_s)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i,j,k
            REAL(mytype), INTENT(IN) :: t
            REAL(mytype) :: exact_s
            exact_s = sin(i*d_x) * sin(j*d_x) * sin(k*d_x) * sin(t)

        END FUNCTION exact_sol

        REAL FUNCTION forcing_term(i,j,k,t)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i,j,k
            REAL(mytype), INTENT(IN) :: t

            forcing_term = (cos(t) + 3*sin(t))*sin(i*d_x)*sin(j*d_x)*sin(k*d_x)

        END FUNCTION forcing_term
END PROGRAM
