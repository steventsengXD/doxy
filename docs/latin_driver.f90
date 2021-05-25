!> @brief Main driver for simulation
!!
!! More description to be added
!!
!> @param n_electrons_in Number of electrons
program main_driver
  use shared_constants
  use read_data
  use basis_functions
  use param_search
  use mcmc
  !> use omp_lib
  implicit none

  !> User Inputs required - the user must define all of these
  !> Will obtain these from either txt file or a gui that generates the txt file
  CHARACTER(LEN=20) :: filename = "init_params.txt"
  integer :: n_electrons_in !> @param n_electrons_in Number of electrons
  integer :: n_atoms_in !> Number of atoms
  real(dp) :: bond_length !> Bond length for 2 atoms
  real(dp), allocatable, dimension(:,:) :: atom_coords_in !> Coordinates of the atoms. Shape (3,n_atoms)
  integer :: n_trials !> number of trials in the hypercube search
  integer :: n_MCMC_steps !> Total number of steps for each MCMC run

  !> Optional user Inputs - these may be defined by the user, but have good default values below
  integer :: n_basis_functions_per_atom_in !> Number of linear terms in the single-electron wavefunction per atom
  integer :: search_seed !> Seed for the latin hypercube
  integer :: n_Jastrow_in !> Number of dofs in the Jastrow (interaction) term
  real(dp) :: fd_length_in !> Lengthscale of the finite difference code
  real(dp) :: Jastrow_b_length_in !> Inverse lengthscale of nuclear-electron interaction
  real(dp) :: Jastrow_d_length_in !> Inverse lengthscale of electron-electron interaction


  !>MCMC variables
  real(dp), allocatable, dimension(:) :: x_0 !> MCMC initial coordinates
  integer ::  n_burned, thinning_interval
  integer ::  e_code
  real(dp) :: s, E, norm_coeff
  real(dp), allocatable, dimension(:,:) :: mcmc_run


  !> Outputs

  real(dp), dimension(:), allocatable :: trial_energies

  !> Internal variables
  real(dp), dimension(:, :), allocatable :: trials !> array of trials
  integer :: i, loop !> Loop variable
  real(dp) :: cpu_start_time, cpu_finish_time !> CPU timing
  integer :: real_start_time,real_finish_time, rate !> Real timing

  !> Start timers
  call cpu_time(cpu_start_time)
  call system_clock(real_start_time,rate)

  call readtxt(filename, n_electrons_in, n_atoms_in, bond_length, \
    n_trials, n_MCMC_steps, n_basis_functions_per_atom_in, search_seed, \
    n_Jastrow_in, fd_length_in, Jastrow_b_length_in, Jastrow_d_length_in)

  allocate( atom_coords_in(3,n_atoms_in) )
  if (n_atoms_in == 1) then
    atom_coords_in(:,1)=[0.0_dp,0.0_dp,0.0_dp] ! 1 atom at origin
  else if (n_atoms_in == 2) then
    atom_coords_in(:,1)=[bond_length/2.0_dp,0.0_dp,0.0_dp] !2 atoms evenly placed about origin
    atom_coords_in(:,2)=[-bond_length/2.0_dp,0.0_dp,0.0_dp]
  end if

! Print to see if data was properly read into fortran
  PRINT *, n_electrons_in
  PRINT *, n_atoms_in
  PRINT *, bond_length
  PRINT *, n_trials
  PRINT *, n_MCMC_steps
  PRINT *, n_basis_functions_per_atom_in
  PRINT *, search_seed
  PRINT *, n_Jastrow_in
  PRINT *, fd_length_in
  PRINT *, Jastrow_b_length_in
  PRINT *, Jastrow_d_length_in
  PRINT *, atom_coords_in(:,1)
  PRINT *, atom_coords_in(:,2)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! REPLACE FOLLOWING WITH INPUT CODE !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! This code is commented out due to enabling txt file input

    ! Commented out and replaced with gui/txt input which follows
    ! ! Assign Required Inputs
    ! n_electrons_in = 2
    ! n_atoms_in = 2
    !
    ! ! Can either take in bond length or the actual coordinates.
    !
    ! bond_length = 2.0_dp
    !
    ! allocate( atom_coords_in(3,n_atoms_in) )
    ! if (n_atoms_in == 1) then
    !   atom_coords_in(:,1)=[0.0_dp,0.0_dp,0.0_dp] ! 1 atom at origin
    ! else if (n_atoms_in == 2) then
    !   atom_coords_in(:,1)=[bond_length/2.0_dp,0.0_dp,0.0_dp] !2 atoms evenly placed about origin
    !   atom_coords_in(:,2)=[-bond_length/2.0_dp,0.0_dp,0.0_dp]
    ! end if
    !
    ! n_trials = 10 !> This is just a test value - should be much larger, and scale with number of dofs
    ! n_MCMC_steps = 10000000 !> This is just a test value, could be larger, try 10^7

    !> Assign optional arguments - these are sensible default values
    ! search_seed = 132
    ! n_basis_functions_per_atom_in = 1 !> For slater 1 is sensible. Increase for Gaussians if we implement them
    ! n_Jastrow_in = 3 !> This is low. Would like to run with 7, but that might be too high for testing
    ! fd_length_in = 0.1_dp !> Not sure on this
    ! Jastrow_b_length_in = 1.0_dp !> This is fine, but could experiment with this
    ! Jastrow_d_length_in = 1.0_dp !> This is fine, but could experiment with this


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! END OF INPUT CODE!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> Initialisation of basis and trials

  ! call initialise_basis(n_electrons_in, n_basis_functions_per_atom_in, n_atoms_in, atom_coords_in,&
  !                       n_Jastrow_in, fd_length_in, Jastrow_b_length_in, Jastrow_d_length_in)
  !
  !
  ! allocate(trials(n_trials, number_dofs))
  ! call  latin_hypercube(trials,dof_bounds, n_trials, search_seed)
  !
  ! allocate(trial_energies(size(trials,1)))
  !
  ! !> MCMC setup
  ! n_burned = n_MCMC_steps/10 !> Integer division
  ! thinning_interval = n_MCMC_steps/10000
  ! allocate(x_0(n_space_dims))
  ! x_0(1:3) = atom_coords_in(:,1) + 0.2_dp !> Displace starting electrons from atoms
  ! if (n_electrons_in == 2) then
  !   x_0(4:6) = atom_coords_in(:,2) + 0.2_dp
  ! end if
  ! allocate(mcmc_run((n_MCMC_steps-n_burned)/thinning_interval+1,n_space_dims))
  ! norm_coeff = real(thinning_interval,dp)/real(n_MCMC_steps-n_burned,dp)
  !
  !
  ! !> Main loop over trials in parameter space
  ! !$OMP parallel do default(shared) private(i,s,e_code,mcmc_run, E,loop)
  ! do i=1, size(trials,1)
  !   !> Run MCMC to generate coordinates in space
  !   call mcmc_adapt(s, log_density, x_0, 10000, 0.1_dp, e_code, 0.4_dp, 0.03_dp, 500.0_dp, 100, trials(i,:),n_space_dims)
  !   call mcmc_sample(mcmc_run, log_density, x_0, n_MCMC_steps, n_burned, thinning_interval, s,e_code,&
  !    trials(i,:),n_space_dims)
  !
  !
  !   !> Compute energy by summing over positions
  !   E=0.0_dp
  !   do loop=1,(n_MCMC_steps-n_burned)/thinning_interval
  !     E = E + reduced_hamiltonian(mcmc_run(loop,:),trials(i,:))
  !   end do
  !
  !   trial_energies(i) = norm_coeff*E
  !   print*, "completed trial" ,i," out of ", size(trials,1)
  !   print*, "energy = ",trial_energies(i)
  ! end do
  !
  ! call find_best_params(trial_energies,trials)
  ! call param_wall_check(dof_bounds)
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!! OUTPUT CODE GOES HERE!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print*,"best trial=", best_trial
  ! print*,"min energy=", trial_energies(best_trial),"Hartree Energy"
  ! print*,"min energy=", electronvolt*trial_energies(best_trial),"eV"
  ! print*,"best parameters=", best_params
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!! OUTPUT CODE FINISH   !!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! !> deallocation
  ! deallocate(trials,trial_energies,atom_coords_in,x_0,mcmc_run)
  ! call deinitialise_basis
  !
  ! !> Report Timings
  ! call cpu_time(cpu_finish_time)
  ! call system_clock(real_finish_time,rate)
  ! print '("Total cpu runtime = ",f7.3," seconds.")',cpu_finish_time-cpu_start_time
  ! print'("Total real runtime = ",f7.3," seconds.")',real(real_finish_time-real_start_time,kind=dp)/real(rate,kind=dp)

end program main_driver
