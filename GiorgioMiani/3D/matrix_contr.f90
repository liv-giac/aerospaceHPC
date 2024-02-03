program MATRIX

    ! use csv_file

    INTEGER, PARAMETER :: n = 3
    INTEGER, PARAMETER :: dim = n*n*n
    REAL(KIND = 4), PARAMETER :: deltax = 1.0/(n + 1)
    REAL(KIND = 4), DIMENSION(0:dim - 1, 0:dim - 1) :: A

    INTEGER :: i
    INTEGER :: j
    INTEGER :: k

    REAL(KIND = 4), PARAMETER :: diag_coef = 6.0/(deltaX*deltaX)
    REAL(KIND = 4), PARAMETER :: coef = -1.0/(deltaX*deltaX)

    A = 0.0

    DO i = 0,n-1
        DO j = 0,n-1
            DO k = 0,n-1
                A(i + j*n + k*n*n, i + j*n + k*n*n) = diag_coef
                if ( i > 0 ) then
                    A(i - 1 + j*n + k*n*n, i + j*n + k*n*n) = coef
                end if
                if ( i < n - 1 ) then
                    A(i + 1 + j*n + k*n*n, i + j*n + k*n*n) = coef
                end if
                if ( j > 0 ) then
                    A(i + (j - 1)*n + k*n*n, i + j*n + k*n*n) = coef
                end if
                if ( j < n - 1 ) then
                    A(i + (j + 1)*n + k*n*n, i + j*n + k*n*n) = coef
                end if
                if ( k > 0 ) then
                    A(i + j*n + (k - 1)*n*n, i + j*n + k*n*n) = coef
                end if
                if ( k < n - 1 ) then
                    A(i + j*n + (k + 1)*n*n, i + j*n + k*n*n) = coef
                end if
            END DO
        END DO
    END DO

    ! Write matrix A to a CSV file
    OPEN(UNIT=10, FILE='matrix_sparsity.csv', STATUS='REPLACE', ACTION='WRITE')
    ! DO i = 0, dim - 1
    !     WRITE(10, '(27F4.0)') A(i, :)
    ! END DO
    
    DO i = 0, dim - 1
        DO j = 0, dim - 1
            IF (A(i,j) /= 0.0) THEN
                WRITE(10 ,*) i,',', j
            END IF
        END DO
    END DO
    CLOSE(UNIT=10)


    WRITE(*,*) "Matrix A has been written to matrix.csv"
END

    
