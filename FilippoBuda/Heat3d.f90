PROGRAM heat34
    IMPLICIT NONE
    INTEGER, PARAMETER :: n = 5
    INTEGER :: matn
    REAL(KIND = 8),  DIMENSION(:,:), ALLOCATABLE :: mat
    REAL(KIND = 8) :: dx, dt, len , time
    INTEGER :: i, j
    REAL(KIND = 8) :: one_dt, three_dx2, one_2dx2
    INTEGER :: modn, modn_n

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
        
        modn = MOD(i,n)
        modn_n = MOD(i, (n*n) + 1) / n+1

        ! upper second diagonal T_i+1,j,k
        IF( (i .LE. matn-1) .AND. (modn .EQ. 0) ) THEN   
            mat(i,i+1) = 0.d0
        ELSE IF ( (i .LE. matn-1) .AND. (modn .NE.  0) ) THEN
            mat(i,i+1) = one_2dx2 
        END IF

        ! lower second diagonal T_i-1,j,k
        IF((i .GE. 2) .AND.  ( modn .EQ. 1 ) ) THEN
            mat(i, i-1) = 0.d0
        ELSE IF ((i .GE. 2) .AND. ( modn .NE. 1 )) THEN
            mat(i, i-1) = one_2dx2
        END IF

        ! right diagonals T_i,j+1,k !!!!!!!!!!!!!!!!!!!!!!
        IF ( (i .LE. (matn - n) ) .AND. (modn_n .EQ. n) ) THEN
            mat (i, i + n) = 0
        ELSE IF ( (i .LE. (matn - n )) .AND. (modn_n .NE. n) ) THEN
            mat(i, i + n) = one_2dx2
        END IF

        ! left diagonal T_i,j-1,k
        IF ((i .GE. 4) .AND. (modn_n .EQ. 1)  ) THEN
            mat (i, i - n) = 0
        ELSE IF ((i .GE. 4) .AND. (modn_n .NE. 1)  ) THEN
            mat (i, i - n) = one_2dx2
        END IF
        
        ! rightmost diagonal T_i,j,k+1
        IF ( i .LE. (matn - n**2)) THEN
            mat (i, i + n**2) = one_2dx2
        END IF

        ! leftmost diagonal T_i,j,k-1
        IF (i .GE. (n**2 + 1)) THEN
            mat (i, i - n**2) = one_2dx2
        END IF

    END DO

    DO i=1, matn
        DO j=1, matn
            WRITE(*, '(I3,X)', advance='no') NINT((mat(i,j)))
        END DO
        WRITE(*,*) ''
    END DO

END PROGRAM
