! This module contains the subroutine write_data which writes our data onto
! a file in binary format. This subroutine is based on the write_netcdf_array.f90
! subroutine from Workshop 7.

MODULE write_netcdf

  USE ISO_FORTRAN_ENV
  USE netcdf

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE write_data(grid_final, grid_init, mag, xaxis, yaxis, taxis, init, &
    beta, J, filename, ierr)

    ! Declaring the axes/variables that will written to the file
    INTEGER(INT32), INTENT(IN), DIMENSION(:,:) :: grid_final, grid_init
    INTEGER(INT32), INTENT(IN), DIMENSION(:) :: xaxis, yaxis, taxis
    REAL(REAL64), INTENT(IN), DIMENSION(:) :: mag

    ! Declaring the variables that will be  as global attributes
    REAL(REAL64), INTENT(IN) :: beta, J
    CHARACTER(LEN=5), INTENT(IN) :: init

    ! Parameters representing the rank our variables, useful
    INTEGER, PARAMETER :: ndims = 2, tdims = 1

    ! Declaring the variables sizes, sizet that store the numerical lengths of
    ! our variables and the variable dim_ids which store the dimension ids for
    ! our grid data
    INTEGER, DIMENSION(ndims) :: sizes, dim_ids
    INTEGER, DIMENSION(tdims) :: sizet

    ! Declaring string names for the dimensions of our variables
    CHARACTER(LEN=1), DIMENSION(ndims) :: dims=(/"x", "y"/)
    CHARACTER(LEN=1) :: dimt = "t"

    ! Declaring the variable ids for each of our variables
    INTEGER(INT32) :: x_id, y_id, t_id, file_id, dim_idt
    INTEGER(INT32) :: grid_final_id, grid_init_id, mag_id

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(INT32) :: ierr, i

    ! Obtain the dimensions of our variables
    sizes = SHAPE(grid_final)
    sizet = SHAPE(taxis)

    ! Create the netCDF file, overwriting a file if it has the same name
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    ! This indicates a writing step error, prompting a return from subroutine
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! We define the dimensions for anything related to the grid here
    ! These dimensions are used for the rank 2 arrays and their axes
    DO i = 1, ndims
      ierr = nf90_def_dim(file_id, dims(i), sizes(i), dim_ids(i))
      IF (ierr /= nf90_noerr) THEN
        PRINT*, TRIM(nf90_strerror(ierr))
        RETURN
      END IF
    END DO

    ! We define the dimension of the magnetization and its axis here
    ierr = nf90_def_dim(file_id, dimt, sizet(1), dim_idt)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Define the variable types for our axes for x, y and t, respectively
    ierr = nf90_def_var(file_id, "x", NF90_INT, dim_ids(1), x_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id, "y", NF90_INT, dim_ids(2), y_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id, "t", NF90_INT, dim_idt, t_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Define the variable types for our three Ising model variables
    ierr = nf90_def_var(file_id, "final_grid", NF90_INT, dim_ids, grid_final_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id, "initial_grid", NF90_INT, dim_ids, grid_init_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_def_var(file_id, "net_magnetization", NF90_REAL, dim_idt, mag_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Adding global attribute describing data
    ierr = nf90_put_att(file_id, NF90_GLOBAL, "Description", "Initial grid, final &
    &grid after T_timesteps and net magnetization history for Ising model simulation, &
    & program file names and command line input given below.")
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Adding global attribute listing the code files used for simulation
    ierr = nf90_put_att(file_id, NF90_GLOBAL, "CodeFiles", "ising.f90, sweetener.f90, &
    & write_netcdf.f90, command_line.f90, random_mod.f90, ising.sh")
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Add global attributes consisting of our five run data items
    ! This includes the total iterations
    ierr = nf90_put_att(file_id, NF90_GLOBAL, "N_gridsize", sizes(1))
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "T_timesteps", sizet(1)-1)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "init_config", init)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "beta_temp", beta)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_put_att(file_id, NF90_GLOBAL, "J_interstr", J)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Finish metadata definitions
    ierr = nf90_enddef(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Write the three axes to file
    ierr = nf90_put_var(file_id, x_id, xaxis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_put_var(file_id, y_id, yaxis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ierr = nf90_put_var(file_id, t_id, taxis)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Write the final grid state
    ierr = nf90_put_var(file_id, grid_final_id, grid_final)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Write the final grid state
    ierr = nf90_put_var(file_id, grid_init_id, grid_init)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Write magnetization
    ierr = nf90_put_var(file_id, mag_id, mag)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

    ! Close the file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) THEN
      PRINT*, TRIM(nf90_strerror(ierr))
      RETURN
    END IF

  END SUBROUTINE write_data

END MODULE write_netcdf
