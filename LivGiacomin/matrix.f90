PROGRAM matrix
    IMPLICIT NONE
    
    INTEGER :: n, NT, m, i, j, k
    REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: mat
    REAL(KIND=8) :: dx, dt, a, b
    READ(*,*) n
    NT = n * n * n
    ALLOCATE (mat(NT, NT))
 
    mat = 0.0d0
    dx = 0.1
    dt = 0.17
    a = 1.0 / dt 
    b = 1.0 / (dx *dx)

    DO k = 1, n
        DO j = 1, n
          DO i = 1, n
            m = i + n * (j - 1) + n * n * (k - 1)
            mat(m, m) = a + 3.0 * b

            IF (i > 1) mat(m, m-1) = b / 2.0 ! Left (x)
            IF (i < n) mat(m, m+1) = b / 2.0  ! Right (x)
            IF (j > 1) mat(m, m-n) = b / 2.0  ! Below (y)
            IF (j < n) mat(m, m+n) = b / 2.0  ! Above (y)
            IF (k > 1) mat(m, m-n*n) = b / 2.0  ! In Front (z)
            IF (k < n) mat(m, m+n*n) = b / 2.0  ! Behind (z)
          END DO
        END DO
      END DO

      DO i = 1, NT
        DO j = 1, NT
          WRITE(*, '(F6.2, 2X)', advance='no') mat(i, j)
        END DO
        PRINT *
      END DO
    
    END PROGRAM matrix
    