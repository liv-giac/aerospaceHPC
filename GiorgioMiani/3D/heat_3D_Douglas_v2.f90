
PROGRAM Heat_Eq_TD_1D
    IMPLICIT NONE

    INTEGER ,PARAMETER :: Nx = 100, Nt = 100
    INTEGER, PARAMETER ::dim = Nx*Nx*Nx
    INTEGER, EXTERNAL :: findex
    REAL(KIND = 8) ,PARAMETER :: deltax = 1.d0/ (1 + Nx), deltat = 1.d0 / Nt, deltax2 = deltaX*deltaX
    REAL(KIND = 8) ,PARAMETER :: diag = 1.d0 +  deltaT/(deltaX*deltaX)
    REAL(KIND = 8), PARAMETER :: extra_diag = -deltaT/(2*deltaX*deltaX)
    REAL(KIND = 8) :: t
    REAL(KIND = 8), DIMENSION(1:Nx,1:Nx,1:Nx) :: rhs, old_solution, solution, b, real_solution, w
    REAL(KIND = 8), DIMENSION(0:Nx + 1, 0:Nx + 1, 0:Nx + 1) :: xsin
    real(KIND = 8) :: a, sin_1
    REAL(KIND = 8) :: error = 0.d0
    REAL :: start, finish, startx, finishx, totalx, starty, finishy, totaly, startz, finishz, totalz
    INTEGER ::  i
    INTEGER :: j
    INTEGER :: k

    do k = 0, Nx + 1
        do j = 0, Nx + 1
            do i = 0, Nx + 1
                xsin(i, j, k) = sin(i * deltaX)*sin(j * deltaX)*sin(k * deltaX)
            end do
        end do
    end do
    sin_1 = sin(1.d0)

    old_solution(1:Nx, 1:Nx, 1:Nx) = 0.d0
    a = extra_diag
    t = 0.d0 + deltaT

    call cpu_time(start)

    DO WHILE ( t <= 1.d0)

        call cpu_time(startx)
        
        !RHS
        !THOMAS ON X
        b(1:Nx, 1:Nx, 1:Nx) = diag

 

        rhs(1:Nx, 1:Nx, 1:Nx) = deltaT*xsin(1:Nx,1:Nx,1:Nx)*(3*sin(t + deltaT/2) + cos(t + deltaT/2)) - (1/(deltax2))*(6 * old_solution(1:Nx, 1:Nx, 1:Nx))

    

        rhs(Nx, 1:Nx, 1:Nx) = rhs(Nx, 1:Nx, 1:Nx) + deltaT*(1/(deltax2))*old_solution(Nx - 1, 1:Nx, 1:Nx) + deltaT*(1/(deltax2))*(xsin(Nx,1:Nx,1:Nx)*sin(t))
        rhs(1, 1:Nx, 1:Nx) = rhs(1, 1:Nx, 1:Nx) + deltaT*(1/(deltax2))*old_solution(2, 1:Nx, 1:Nx)
        rhs(2:Nx - 1, 1:Nx, 1:Nx) = rhs(2:Nx - 1, 1:Nx, 1:Nx) + deltaT*(1/(deltax2))*old_solution(1:Nx - 2, 1:Nx, 1:Nx) + deltaT*(1/(deltax2))*old_solution(3:Nx, 1:Nx, 1:Nx)

        rhs(1:Nx, Nx, 1:Nx) = rhs(1:Nx, Nx, 1:Nx) + deltaT*(1/(deltax2))*old_solution(1:Nx, Nx - 1, 1:Nx) + deltaT*(1/(deltax2))*(xsin(1:Nx,Nx,1:Nx)*sin(t))
        rhs(1:Nx, 1, 1:Nx) = rhs(1:Nx, 1, 1:Nx) + deltaT*(1/(deltax2))*old_solution(1:Nx, 2, 1:Nx)
        rhs(1:Nx, 2:Nx - 1, 1:Nx) = rhs(1:Nx, 2:Nx - 1, 1:Nx) + deltaT*(1/(deltax2))*old_solution(1:Nx, 1:Nx - 2, 1:Nx) + deltaT*(1/(deltax2))*old_solution(1:Nx, 3:Nx, 1:Nx)

        rhs(1:Nx, 1:Nx, Nx) = rhs(1:Nx, 1:Nx, Nx) + deltaT*(1/(deltax2))*old_solution(1:Nx, 1:Nx, Nx - 1) + deltaT*(1/(deltax2))*(xsin(1:Nx,1:Nx,Nx)*sin(t))
        rhs(1:Nx, 1:Nx, 1) = rhs(1:Nx, 1:Nx, 1) + deltaT*(1/(deltax2))*old_solution(1:Nx, 1:Nx, 2)
        rhs(1:Nx, 1:Nx, 2:Nx - 1) = rhs(1:Nx, 1:Nx, 2:Nx - 1) + deltaT*(1/(deltax2))*old_solution(1:Nx, 1:Nx, 1:Nx - 2) + deltaT*(1/(deltax2))*old_solution(1:Nx, 1:Nx, 3:Nx)



        rhs(Nx, 1:Nx, 1:Nx) = rhs(Nx, 1:Nx, 1:Nx) - extra_diag*xsin(Nx,1:Nx,1:Nx)*(sin(t + deltaT) - sin(t))
        w(1:Nx - 1, 1:Nx, 1:Nx) = a/b(1:Nx - 1, 1:Nx, 1:Nx)
        rhs(2:Nx, 1:Nx, 1:Nx) = rhs(2:Nx, 1:Nx, 1:Nx) - w(1:Nx - 1, 1:Nx, 1:Nx)*rhs(1:Nx - 1, 1:Nx, 1:Nx)
        b(2:Nx, 1:Nx, 1:Nx) = b(2:Nx, 1:Nx, 1:Nx) - w(1:Nx - 1, 1:Nx, 1:Nx)*a

        solution(Nx, 1:Nx, 1:Nx) = rhs(Nx, 1:Nx, 1:Nx)/b(Nx, 1:Nx, 1:Nx)
        do i = Nx - 1, 1, -1
            solution(i, 1:Nx, 1:Nx) = (rhs(i, 1:Nx, 1:Nx) - a * solution(i + 1, 1:Nx, 1:Nx))/b(i, 1:Nx, 1:Nx)
        end do

        call cpu_time(finishx)

        totalx = totalx + finishx - startx

        call cpu_time(starty)

        !THOMAS ON Y
        rhs(1:Nx, 1:Nx, 1:Nx) = solution(1:Nx, 1:Nx, 1:Nx)
        b(1:Nx, 1:Nx, 1:Nx) = diag
        
        rhs(1:Nx, Nx, 1:Nx) = rhs(1:Nx, Nx, 1:Nx) + deltaT*(1/(deltax2))*old_solution(1:Nx, Nx - 1, 1:Nx) + deltaT*(1/(deltax2))*(xsin(1:Nx,Nx,1:Nx)*sin(t))
        w(1:Nx, 1:Nx - 1, 1:Nx) = a/b(1:Nx, 1:Nx - 1, 1:Nx)
        rhs(1:Nx, 2:Nx, 1:Nx) = rhs(1:Nx, 2:Nx, 1:Nx) - w(1:Nx, 1:Nx - 1, 1:Nx)*rhs(1:Nx, 1:Nx - 1, 1:Nx)
        b(1:Nx, 2:Nx, 1:Nx) = b(1:Nx, 2:Nx, 1:Nx) - w(1:Nx, 1:Nx - 1, 1:Nx)*a

        solution(1:Nx, Nx, 1:Nx) = rhs(1:Nx, Nx, 1:Nx)/b(1:Nx, Nx, 1:Nx)
        do j = Nx - 1, 1, -1
            solution(1:Nx, j, 1:Nx) = (rhs(1:Nx, j, 1:Nx) - a * solution(1:Nx, j + 1, 1:Nx))/b(1:Nx, j, 1:Nx)
        end do

        call cpu_time(finishy)

        totaly = totaly + finishy - starty

        call cpu_time(startz)
        !THOMAS ON Z
        rhs(1:Nx, 1:Nx, 1:Nx) = solution(1:Nx, 1:Nx, 1:Nx)
        b(1:Nx, 1:Nx, 1:Nx) = diag

        rhs(1:Nx, 1:Nx, Nx) = rhs(1:Nx, 1:Nx, Nx) + deltaT*(1/(deltax2))*old_solution(1:Nx, 1:Nx, Nx - 1) + deltaT*(1/(deltax2))*(xsin(1:Nx,1:Nx,Nx)*sin(t))
        w(1:Nx, 1:Nx, 1:Nx - 1) = a/b(1:Nx, 1:Nx, 1:Nx - 1)
        rhs(1:Nx, 1:Nx, 2:Nx) = rhs(1:Nx, 1:Nx, 2:Nx) - w(1:Nx, 1:Nx, 1:Nx - 1)*rhs(1:Nx, 1:Nx, 1:Nx - 1)
        b(1:Nx, 1:Nx, 2:Nx) = b(1:Nx, 1:Nx, 2:Nx) - w(1:Nx, 1:Nx, 1:Nx - 1)*a

        solution(1:Nx, 1:Nx, Nx) = rhs(1:Nx, 1:Nx, Nx)/b(1:Nx, 1:Nx, Nx)
        do k = Nx - 1, 1, -1
            solution(1:Nx, 1:Nx, k) = (rhs(1:Nx, 1:Nx, k) - a * solution(1:Nx, 1:Nx, k + 1))/b(1:Nx, 1:Nx, k)
        end do 
 

        old_solution = solution + old_solution

        t = t + deltaT

        call cpu_time(finishz)

        totalz = totalz + finishz - startz

    END DO


    call cpu_time(finish)

    Write (*,*) "Time x axis: "
    Write (*,*) totalx

    Write (*,*) "Time y axis: "
    Write (*,*) totaly

    Write (*,*) "Time z axis: "
    Write (*,*) totalz

    Write (*,*) "Total time: "
    Write (*,*) finish - start
    
    real_solution(1:Nx, 1:Nx, 1:Nx) = xsin(1:Nx,1:Nx,1:Nx)*sin(1.d0)

    error = NORM2(old_solution - real_solution)/dim
    
    WRITE(*,*) "Error of "
    WRITE(*,*) error

END

integer function findex(o,p,q, N)
    integer, intent (in) :: o, p, q , N
    findex = o + (p - 1) * N + (q - 1) * N * N
end function findex
