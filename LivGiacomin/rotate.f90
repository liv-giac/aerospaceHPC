PROGRAM rotate
    IMPLICIT NONE
    INTEGER, PARAMETER :: Nx = 3
    INTEGER, PARAMETER :: Ny = 3
    INTEGER, PARAMETER :: Nz = 3
    INTEGER, PARAMETER :: N(3) = [Nx, Ny, Nz]
    
    CHARACTER(len=20) :: source(Nx*Ny*Nz)
    CHARACTER(len=20) :: destination(Nx*Ny*Nz)

    INTEGER :: i, j, k, h
    CHARACTER :: value

    IF (COMMAND_ARGUMENT_COUNT() == 0) then
        ! No command-line argument provided, prompt the user
        WRITE(*, *) "Enter the value of h (0, 1, or 2): "
        READ(value, *) h
      ELSE
        ! Get h from command-line argument
        CALL get_command_argument(1, value)
        READ(value, *) h
        ! Check IF the argument is a valid integer
        IF (h /= 0 .and. h /= 1 .and. h /= 2) then
          WRITE(*, *) "Invalid value for h. Please enter 0, 1, or 2."
          STOP
        END IF
      END IF

    ! Initialize source array
    DO i = 1, Nz
        DO j = 1, Ny
            DO k = 1, Nx
                source(k + (j - 1) * Nx + (i - 1) * Nx * Ny) = 'T_' // ACHAR(k+IACHAR('0')) // &
                '_' // ACHAR(j+IACHAR('0')) // '_' // ACHAR(i+IACHAR('0'))
                WRITE(*,*) source(k + (j - 1) * Nx + (i - 1) * Nx * Ny)
            END DO
        END DO
    END DO

    !Rotation
    WRITE (*,*)
    IF (h == 0) WRITE (*,*) "Now rotating around axis x."
    IF (h == 1) WRITE (*,*) "Now rotating around axis y."
    IF (h == 2) WRITE (*,*) "Now rotating around axis z."
    WRITE (*,*)

    DO i = 1, N(mod(h+2,3)+1)
        DO j = 1, N(mod(h+1,3)+1)
            DO k = 1, N(mod(h,3)+1)
                destination(j + (i-1) * N(mod(h+1,3)+1) + (k-1) * N(mod(h+1,3)+1) * N(mod(h+2,3)+1)) = &
                    source(k + (j - 1) * N(mod(h,3)+1) + (i - 1) * N(mod(h,3)+1) * N(mod(h+1,3)+1))
                WRITE(*,*) destination(j + (i - 1) * N(mod(h+1,3)+1) + (k-1) * N(mod(h+1,3)+1) * N(mod(h+2,3)+1))
            END DO
        END DO
    END DO


END PROGRAM rotate
