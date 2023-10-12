PROGRAM heat34
    IMPLICIT NONE
    INTEGER, PARAMETER :: n = 3
    INTEGER :: matn
    REAL(KIND = 8),  DIMENSION(:,:), ALLOCATABLE :: mat
    REAL(KIND = 8) :: dx, dt, len , time
    INTEGER :: i, j
    REAL(KIND = 8) :: one_dt, three_dx2, one_2dx2
    INTEGER :: mod3,mod4

    len = 1.d0
    time = 1.d0
    dx = 0.1d0 !len / (n-1)
    dt = 0.1d0 !time / 10 !chemmerda

    one_dt = 1.d0 / dt
    three_dx2 = 3.d0 / (dx**2)
    one_2dx2 = - 1.d0 / ( 2 * (dx**2))
 
    
    matn = n**3
    allocate(mat(matn,matn))
    mat = 0

    DO i = 1, matn
        ! diagonal T_i,j,k
        mat(i,i) = one_dt + three_dx2
        
        mod3 = MOD(i,n)
        mod4 = MOD(i,n+1)

        ! upper second diagonal T_i+1,j,k
        IF( (i .LE. matn-1) .AND. (mod3 .EQ. 0) ) THEN   
            mat(i,i+1) = 0.d0
        ELSE IF ( (i .LE. matn-1) .AND. (mod3 .NE.  0) ) THEN
            mat(i,i+1) = one_2dx2 
        END IF

        ! lower second diagonal T_i-1,j,k
        IF((i .GE. 2) .AND. (mod4 .EQ. 0) ) THEN
            mat(i, i-1) = 0.d0
        ELSE IF ((i .GE. 2) .AND. (mod4 .NE. 0) ) THEN
            mat(i, i-1) = one_2dx2
        END IF

        ! right diagonals T_i,j+1,k
        IF ( i .LE. (matn - 3)) THEN
            mat (i, i + 3) = one_2dx2
        END IF

        ! rightmost diagonal T_i,j,k+1
        IF ( i .LE. (matn - n**2)) THEN
            mat (i, i + 9) = one_2dx2
        END IF

        ! leftmost diagonal T_i,j,k-1
        IF (i .GE. 4) THEN
            mat (i, i - 3) = one_2dx2
        END IF
        
        ! left diagonal T_i,j-1,k
        IF (i .GE. 10) THEN
            mat (i, i - 9) = one_2dx2
        END IF

    END DO

    DO i=1, matn
        DO j=1, matn
            WRITE(*, '(I4,X)', advance='no') NINT((mat(i,j)))
        END DO
        WRITE(*,*) ''
    END DO

END PROGRAM