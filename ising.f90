! Written by Steven Tseng
! This program simulates the Ising Spin Model and writes the results to an output
! file 'ising.nc' in the netcdf format.

! An accompanying bash script 'ising.sh' links the necessary libraries, compiles
! the program files and runs the simulation at the values given below in the
! example terminal command. For reproducibility, the initial grid along with the
! command line inputs are written to the output file; these may be utilized as
! initial conditions.

! Example terminal command with arguments for running the program:
! ./ising N=50 T=1000000 init=R beta=2 J=1
! In the above example, the Ferromagnet starts from an ordered state
! and transitions into an unordered state. As shown, the user will enter
! their desired square grid size (N), the number of timesteps to
! iterate (T), the initial state (init), the temperature (beta) and
! the spin-spin interaction strength (J).

! Citations: For the two subroutines that initialize for alternating and random grids,
! I incorporated codes from the model solutions in the DO loops that determine
! if a cell should be set to -1. For writing the data, I used the write_netcdf_array.f90
! module and added additional codes to write all the axes/data to file. Also, as before,
! I am using the modules command_line and random_mod to access command line arguments
! and get random numbers.

PROGRAM main

  ! We import the modules that are utilized by this program.
  USE sweetener
  USE write_netcdf
  USE ISO_FORTRAN_ENV

  IMPLICIT NONE

  ! We declare the variables used in the program.
  ! The first five variables are input through command line arguments.
  ! N is the size of the NxN grid
  ! T is the number of timesteps
  ! init is one of {a or A: Alternating; f or F: Ferromagnet; r or R: Random}
  ! beta is the temperature
  ! J is the spin-spin interaction strength
  INTEGER(INT32) :: N, T
  CHARACTER(LEN=5) :: init
  REAL(REAL64) :: beta, J
  ! Here we declare the variables that will be stored on the datafile:
  ! two allocatable rank two arrays for the initial grid and final grid, and
  ! one allocatable rank one array for the net magnetization at each timestep
  INTEGER(INT32), DIMENSION(:, :), ALLOCATABLE :: grid, grid_init
  REAL(REAL64), DIMENSION(:), ALLOCATABLE :: mag
  ! Three rank one arrays are declared to be used as the axes for x, y, t
  INTEGER(INT32), DIMENSION(:), ALLOCATABLE :: xaxis, yaxis, taxis
  ! ierr variable here is passed to write_data subroutine to store error codes
  INTEGER(INT32) :: ierr, count = 0
  ! The data file that will be output at the end is named here
  CHARACTER(LEN=30) :: filename = "ising.nc"

  ! We use the subroutine 'parse' to check our command line arguments to make
  ! sure they are suitable for this simulation
  CALL parse(N, T, init, beta, J)

  ! Initialize the grid with the selected/default state whilst storing a copy of
  ! the initial grid for comparison purposes
  CALL initialize(grid, grid_init, N, init)

  ! We allocate the net magnetization array with dimension one greater than the
  ! total number of timesteps and then set the first value to the net magnetization
  ! of the initial grid
  ALLOCATE(mag(T+1))
  mag(1) = magnetization(grid_init, N)

  ! Iterate the Ising Spin Model for T timesteps
  DO WHILE(count < T)
    CALL simulate(grid, N, beta, J)
    mag(count+2) = magnetization(grid, N)
    count = count + 1
  END DO

  ! The axes for x, y, t are created here
  CALL define_axis(xaxis, N, 1)
  CALL define_axis(yaxis, N, 1)
  CALL define_axis(taxis, T, 0)

  ! The data/axes are written to file here
  CALL write_data(grid, grid_init, mag, xaxis, yaxis, taxis, init, beta, J, &
  filename, ierr)

  ! We deallocate allocatble array variables
  DEALLOCATE(grid, grid_init)
  DEALLOCATE(mag, xaxis, yaxis, taxis)

END PROGRAM main
