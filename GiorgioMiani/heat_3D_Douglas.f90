
PROGRAM Heat_Eq_TD_1D
    
    IMPLICIT NONE

    INTEGER ,PARAMETER :: Nx = 100, Nt = 100
    INTEGER, PARAMETER ::dim = Nx*Nx*Nx
    INTEGER, EXTERNAL :: findex
    REAL(KIND = 8) ,PARAMETER :: deltax = 1.d0/ (1 + Nx), deltat = 1.d0 / Nt, deltax2 = deltaX*deltaX
    REAL(KIND = 8) ,PARAMETER :: diag = 1.d0 +  deltaT/(deltaX*deltaX)
    REAL(KIND = 8), PARAMETER :: extra_diag = -deltaT/(2*deltaX*deltaX)
    REAL(KIND = 8) :: t
    REAL(KIND = 8) :: w
    REAL(KIND = 8), DIMENSION(1:dim) :: rhs, old_solution, solution, b, real_solution
    real(kind = 8) :: a
    REAL(KIND = 8) :: error = 0.d0
    real(kind = 8), DIMENSION(Nx) :: xsin
    REAL :: start1, finish1, total1, start2, finish2, total2, start3, finish3, total3
    INTEGER ::  i
    INTEGER :: j
    INTEGER :: k

    DO i = 1, Nx
        xsin(i) = sin(i * deltaX)
    END DO

    old_solution = 0.d0
    b = diag
    a = extra_diag
    t = 0.d0 + deltaT

    !call cpu_time(start)

    DO WHILE ( t <= 1.d0)

        call cpu_time(start1)
        
        !RHS
        !THOMAS ON X
        b = diag
        do k = 1, Nx
            do j = 1, Nx
                do i = 1, Nx
                    rhs(findex(i, j, k, Nx)) = deltaT*((xsin(i)*xsin(j)*xsin(k))*(3*sin(t + deltaT/2) + cos(t + deltaT/2)) - (1/(deltax2)*(6 * old_solution(findex(i, j, k, Nx)))))
                    if ( i > 1 ) then
                        rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i - 1, j, k, Nx))
                    end if
                    if ( i < Nx) then
                        rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i + 1, j, k, Nx))
                    else
                        rhs(findex(Nx, j, k, Nx)) = rhs(findex(Nx, j, k, Nx)) + deltaT*(1/(deltax2))*(sin(1.d0)*xsin(j)*xsin(k)*sin(t))
                    end if
                    if ( j > 1 ) then
                        rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, j - 1, k, Nx))
                    end if
                    if ( j < Nx ) then
                        rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, j + 1, k, Nx))
                    else
                        rhs(findex(i, Nx, k, Nx)) = rhs(findex(i, Nx, k, Nx)) + deltaT*(1/(deltax2))*(sin(1.d0)*xsin(i)*xsin(k)*sin(t))                
                    end if
                    if ( k > 1 ) then
                        rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i,j, k - 1, Nx))
                    end if
                    if ( k < Nx ) then
                        rhs(findex(i, j, k, Nx)) = rhs(findex(i, j, k, Nx)) + deltaT*(1/(deltax2))*old_solution(findex(i, j, k + 1, Nx))
                    else
                        rhs(findex(i, j, Nx, Nx)) = rhs(findex(i, j, Nx, Nx)) + deltaT*(1/(deltax2))*(sin(1.d0)*xsin(j)*xsin(i)*sin(t))                
                    end if
                end do 
                !BCS
                rhs(findex(Nx, j, k, Nx)) = rhs(findex(Nx, j, k, Nx)) - extra_diag*sin(1.d0)*xsin(j)*xsin(k)*(sin(t + deltaT) - sin(t))                
                do i = 2, Nx
                    w = a/b(findex(i - 1, j, k, Nx))  
                    rhs(findex(i, j, k ,Nx)) = rhs(findex(i, j, k ,Nx)) - w*rhs(findex(i - 1, j, k, Nx))                      
                    b(findex(i, j, k, Nx)) = b(findex(i, j, k, Nx)) - w*a
                end do                
                solution(findex(Nx, j, k, Nx)) = rhs(findex(Nx, j, k, Nx))/b(findex(Nx, j, k, Nx)) 
                do i = Nx - 1,1, -1 
                    solution(findex(i, j, k, Nx)) = (rhs(findex(i, j, k, Nx)) - a * solution(findex(i + 1, j, k, Nx)))/b(findex(i, j, k, Nx))                  
                end do               
            end do
        end do 

        call cpu_time(finish1)

        total1 = total1 + finish1 - start1

        call cpu_time(start2)

        !THOMAS ON Y
        rhs = solution
        b = diag
        do k = 1, Nx
            do j = 2, Nx
                do i = 1, Nx
                    if(j == Nx) then
                        !BCS
                        rhs(findex(i, Nx, k, Nx)) = rhs(findex(i, Nx, k, Nx)) - extra_diag*sin(1.d0)*xsin(i)*xsin(k)*(sin(t + deltaT) - sin(t))
                    end if
                    w = a/b(findex(i, j - 1, k, Nx))  
                    rhs(findex(i, j, k ,Nx)) = rhs(findex(i, j, k ,Nx)) - w*rhs(findex(i, j - 1, k, Nx))                      
                    b(findex(i, j, k, Nx)) = b(findex(i, j, k, Nx)) - w*a
                end do   
            end do
        end do
        
        do k = 1, Nx
            do j = Nx - 1, 1, -1
                do i = 1, Nx
                    solution(findex(i, Nx, k, Nx)) = rhs(findex(i, Nx, k, Nx))/b(findex(i, Nx, k, Nx))  
                    solution(findex(i, j, k, Nx)) = (rhs(findex(i, j, k, Nx)) - a * solution(findex(i, j + 1, k, Nx)))/b(findex(i, j, k, Nx))                  
                end do
            end do
        end do

        call cpu_time(finish2)

        total2 = total2 + finish2 - start2

        call cpu_time(start3)
        !THOMAS ON Z
        rhs = solution
        b = diag
        do k = 2, Nx
            do j = 1, Nx
                do i = 1, Nx
                    if ( k == Nx ) then
                        !BCS
                        rhs(findex(i, j,Nx, Nx)) = rhs(findex(i, j, Nx, Nx)) - extra_diag*sin(1.d0)*xsin(j)*xsin(i)*(sin(t + deltaT) - sin(t)) 
                    end if
                    w = a/b(findex(i, j,k - 1, Nx))  
                    rhs(findex(i, j, k ,Nx)) = rhs(findex(i, j, k ,Nx)) - w*rhs(findex(i, j, k - 1, Nx))                      
                    b(findex(i, j, k, Nx)) = b(findex(i, j, k, Nx)) - w*a
                end do         
            end do
        end do


        do k = Nx -1, 1, -1
            do j = 1, Nx
                do i = 1, Nx
                    solution(findex(i, j, Nx, Nx)) = rhs(findex(i, j, Nx, Nx))/b(findex(i, j, Nx, Nx))   
                    solution(findex(i, j, k, Nx)) = (rhs(findex(i, j, k, Nx)) - a * solution(findex(i, j, k + 1, Nx)))/b(findex(i, j, k, Nx))                  
                end do
            end do
        end do    
 

        old_solution = solution + old_solution

        t = t + deltaT

        call cpu_time(finish3)

        total3 = total3 + finish3 - start3

    END DO


    !call cpu_time(finish)

    Write (*,*) "Time x axis: "
    Write (*,*) total1

    Write (*,*) "Time y axis: "
    Write (*,*) total2

    Write (*,*) "Time z axis: "
    Write (*,*) total3

    Write (*,*) "Total time: "
    Write (*,*) total1 + total2 + total3
    
    do k = 1, Nx
        do j = 1, Nx
            do i = 1, Nx
                real_solution(findex(i, j, k, Nx)) = xsin(i)*xsin(j)*xsin(k)*sin(1.d0)
            end do
        end do
    end do

    error = NORM2(old_solution - real_solution)/dim
    
    WRITE(*,*) "Error of "
    WRITE(*,*) error

END

integer function findex(o,p,q, N)
    integer, intent (in) :: o, p, q , N
    findex = o + (p - 1) * N + (q - 1) * N * N
end function findex
