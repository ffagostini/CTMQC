program main
  use wigner_distribution
  use tools
  use time_evolution
  use read_BOfiles
  use variables
  implicit none

  call read_system_variables()

  call initialize_system_vars
  call initialize_dynamics_vars
  call read_dynamics_variables()

  call initialize_trajectory_vars

  call read_external_files

  if(el_basis=="diabatic") call change_basis

  call readBOfiles

  call initial_conditions

  call evolution

  call finalize

  contains

  subroutine read_system_variables(unit)

    integer,intent(in),optional :: unit
    integer :: ioerr
    integer :: unit_loc=5

    if(present(unit)) unit_loc=unit

    ioerr=0
    read(unit_loc,system,iostat=ioerr) !< System variables are read 
    !! in the first block of the input file.
    if(ioerr/=0) print*,'error reading system variables'

    !> System variables are used to define the dimensions of the potential energy surfaces
    !! and nonadiabatic coupling arrays, of the xyz grids, and, if necessary, of the
    !! diabatic Hamiltonian.

    npairs=(nstates*(nstates-1))/2

  end subroutine read_system_variables


  subroutine read_dynamics_variables(unit)

    integer,intent(in),optional :: unit
    integer :: ioerr,i_dof
    integer :: unit_loc=5

    if(present(unit)) unit_loc=unit

    ioerr=0
    read(unit_loc,dynamics,iostat=ioerr)
    if(ioerr/=0) print*,'error reading dynamics variables'

    nsteps=anint(final_time/dt)

    do i_dof=1,n_dof
      r0(i_dof)=r_init(i_dof)
      k0(i_dof)=k_init(i_dof)
      sigma(i_dof)=sigma_init(i_dof)
      if (sigma(i_dof)==0) sigma(i_dof)=(20.0_dp/k0(i_dof))
      mass(i_dof)=mass_input(i_dof)
    end do

    !> Dynamics variables are read in the second block of the input file.

  end subroutine read_dynamics_variables


  subroutine read_external_files(unit)

    integer,intent(in),optional :: unit
    integer :: ioerr
    integer :: unit_loc=5

    !> Read external files:

    if(present(unit)) unit_loc=unit

    ioerr=0
    read(unit_loc,external_files,iostat=ioerr) !< The locations of external files
    !! are read in the third block of the input file.
    if(ioerr/=0) print*,'error reading paths to external files'

  end subroutine read_external_files


end program main
