! This module contains subroutines and functions that will check the command
! line input for reasonability, implement the Ising Spin Model and output the
! results. It has been modified since Assignment 2 through incorporating feedback
! and provided Ising model solutions.

MODULE sweetener

  ! We import the modules that are utilized by these subroutines/functions
  USE command_line
  USE random_mod
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  CONTAINS

  ! This subroutine collects/parses through each of the command line arguments
  ! to ensure that they meet the minimal criteria for the program to run and are
  ! sensible in terms of the Ising Spin Model. Default parameters are utilized if
  ! command line input is unreasonable.
  SUBROUTINE parse(N, T, init, beta, J)
    INTEGER(INT32) :: N, T
    CHARACTER(LEN=5) :: init
    REAL(REAL64) :: beta, J
    LOGICAL :: success, exists

    ! We utilize the parse_args subroutine to collect and store all the
    ! command line arguments
    CALL parse_args

    ! get_arg is used on each argument individually to access them. If input is
    ! unreasonable, a default value is set for each of the parameters.
    success = get_arg("N", N, exists=exists)
    IF (.NOT. success) THEN
      N = 50
      PRINT*,"N=50 default was used for grid size as there was no/incorrect &
      &command line input for N."
      PRINT*,""
    END IF

    success = get_arg("T", T, exists=exists)
    IF (.NOT. success) THEN
      T = 1000000
      PRINT*,"T=1000000 default was used for the number of timesteps as there &
      &was no/incorrect command line input for T."
      PRINT*,""
    END IF

    ! The following checks to make sure that init is set to equal to one of
    ! {a, A, f, F, r, R} on the command line. If init is equal to anything else,
    ! a random initial grid is utilized.
    success = get_arg("init", init)
    IF(success)THEN
      IF (init == "A" .OR. init == "a" .OR. init == "F" .OR. init == "f" .OR. &
      init == "R" .OR. init == "r") THEN
        CONTINUE
      END IF
    ELSE
      init = "R"
      PRINT*,"init=R default was used for the initial state as there was &
      &no/incorrect command line input for init."
      PRINT*,""
    END IF

    success = get_arg("beta", beta)
    IF(.NOT. success) THEN
      beta = 2.0_REAL64
      PRINT*,"beta=2.0 default was used for the temperature as there was&
      &no/incorrect command line input for beta."
      PRINT*,""
    END IF

    success = get_arg("J", J)
    IF(.NOT. success) THEN
      J = 1.0_REAL64
      PRINT*,"J=1.0 default was used for the interaction strength as there was &
      &no/incorrect command line input for J."
      PRINT*,""
    END IF
  END SUBROUTINE parse

  ! This initializes the NxN grid with the selected state using one of the
  ! three subroutines that follow
  SUBROUTINE initialize(grid_in, grid_copy, N, init)
    INTEGER(INT32), INTENT(IN) :: N
    INTEGER(INT32), DIMENSION(:, :), ALLOCATABLE :: grid_in, grid_copy
    CHARACTER(LEN=5), INTENT(IN) :: init

    ALLOCATE(grid_in(1:N, 1:N))
    ALLOCATE(grid_copy(1:N, 1:N))

    IF(init == "F" .OR. init == "f") THEN
      CALL init_ferromagnet(grid_in, N)
    ELSE IF(init == "R" .OR. init == "r") THEN
      CALL init_random(grid_in, N)
    ELSE IF(init == "A" .OR. init == "a") THEN
      CALL init_alternating(grid_in, N)
    END IF

    ! This makes a copy of the initial configuration for later
    grid_copy = grid_in
  END SUBROUTINE initialize

  ! This subroutine intializes the grid with all states as +1
  SUBROUTINE init_ferromagnet(grid_in, N)
    INTEGER(INT32), INTENT(IN) :: N
    INTEGER(INT32), DIMENSION(1:N, 1:N), INTENT(OUT) :: grid_in
    grid_in = 1
  END SUBROUTINE init_ferromagnet

  ! This subroutine initializes the grid cells randomly with +/- 1's
  ! except along the boundary. The DO statement here is incorporated from the
  ! model solution.
  SUBROUTINE init_random(grid_in, N)
    INTEGER(INT32), INTENT(IN) :: N
    INTEGER(INT32) :: i, j
    INTEGER(INT32), DIMENSION(1:N, 1:N), INTENT(OUT) :: grid_in

    ! Setting the grid equal to 1 everywhere, then will randomize the interior
    grid_in = 1

    ! For each non-boundary cell, a random number is drawn and cell is set to -1
    ! if it is < 0.5
    DO i = 2, N-1
      DO j = 2, N-1
        IF(random() < 0.5) grid_in(i, j) = -1
        ! Realized this else statement was unnecessary
        ! ELSE
        !   grid_in(i, j) = 1
        ! DO
      END DO
    END DO
  END SUBROUTINE init_random

  ! This subroutine initializes the grid cells with alternating +/- 1's
  ! similar to a chessboard layout except along the boundary. The double DO
  ! statement below is incorporated from the model solution.
  SUBROUTINE init_alternating(grid_in, N)
    INTEGER(INT32), INTENT(IN) :: N
    INTEGER(INT32) :: i, j
    INTEGER(INT32), DIMENSION(1:N, 1:N), INTENT(OUT) :: grid_in

    ! Setting the grid equal to 1 everywhere, then will alternate the interior
    grid_in = 1

    ! Alternating +/- 1's are intialized using the sum of the cell indices as
    ! the exponent for -1
    DO i = 2, N-1
      DO j = 2, N-1
        ! grid_in(i, j) = (-1)**(MOD(i+j, 2))
        IF(MOD(i+j, 2)==0) grid_in(i,j) = -1
      END DO
    END DO
  END SUBROUTINE init_alternating

  ! This function calculates the magnetization as the average value of the grid
  ! cells
  REAL(REAL64) FUNCTION magnetization(grid, N)
    INTEGER(INT32), INTENT(IN) :: N
    INTEGER(INT32), DIMENSION(1:N, 1:N), INTENT(IN) :: grid
    magnetization = sum(REAL(grid, KIND=REAL64))/(REAL(size(grid), KIND=REAL64))
  END FUNCTION magnetization

  ! This function calculates the change in energy if the spin of the randomly
  ! selected cell were flipped by looking at the spin of neighboring four cells.
  REAL(REAL64) FUNCTION delta_E(grid, N, i, k, J)
    INTEGER(INT32), INTENT(IN) :: N, i, k
    INTEGER(INT32), DIMENSION(1:N, 1:N), INTENT(IN) :: grid
    REAL(REAL64) :: J
    delta_E = J * REAL(grid(i, k) * (grid(i-1, k) + grid(i+1, k) + grid(i, k-1) &
    + grid(i, k+1)), KIND=REAL64)
  END FUNCTION delta_E

  ! This simulates the Ising model for one timestep. Two random numbers between
  ! 2 and N-1 are randomly chosen as the indices for the selected grid cell for
  ! each timestep.
  SUBROUTINE simulate(grid, N, beta, J)
    INTEGER(INT32) :: N
    INTEGER(INT32), DIMENSION(1:N, 1:N) :: grid
    INTEGER(INT32) :: rand_cell_i, rand_cell_j
    REAL(REAL64) :: beta, J, deltaE, P

    ! Add 2 to each index as we don't want to select cells from the boundary
    rand_cell_i = INT(random()*(N-2), KIND=INT32) + 2
    rand_cell_j = INT(random()*(N-2), KIND=INT32) + 2

    ! We obtain the change in energy of the selected cell's spin were flipped
    ! and flip the spin as given by the Ising Spin Model algorithm
    deltaE = delta_E(grid, N, rand_cell_i, rand_cell_j, J)
    IF (deltaE < 0) THEN
      grid(rand_cell_i, rand_cell_j) = -grid(rand_cell_i, rand_cell_j)
    ELSE IF (deltaE > 0) THEN
      P = exp(-beta*deltaE)
      IF(P > random()) THEN
        grid(rand_cell_i, rand_cell_j) = -grid(rand_cell_i, rand_cell_j)
      END IF
    END IF
  END SUBROUTINE simulate

  ! This subroutine creates an axis array when given the dimension and starting value
  ! The values are incremented in steps of 1, similar to numpy.arange(val0,dim+1)
  SUBROUTINE define_axis(axis, dim, val0)
    INTEGER(INT32) :: i, dim, val0
    INTEGER(INT32), DIMENSION(:), ALLOCATABLE :: axis
    ALLOCATE(axis(dim))
    axis = (/(i, i = val0, dim, 1)/)
  END SUBROUTINE define_axis

END MODULE sweetener
