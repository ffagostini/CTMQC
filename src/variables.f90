module variables
  use kinds
  implicit none

  !>@param hbar
  !!@param Im_unit

  real(kind=dp),parameter :: kB=0.0000031577464_dp
  real(kind=dp),parameter :: hbar=1.0_dp !< Atomic units are used throughout the program
  real(kind=dp),parameter :: zero=0.0000000010_dp
  real(kind=dp),parameter :: M_parameter=6.0_dp !<
  !< The M_parameter is used to determine "how far" each trajectory has
  !! to search to find its neighbours, and determine at each time-step
  !! the variance of the associated Gaussian. It is set to 6.0 by default.
  complex(kind=dp) :: Im_unit !< Imaginary unit
  real(kind=dp) :: PI

  character(len=100) :: model_system="unknown"
  !< The value is set to "unknown" by default but it can be set to tully_1, 
  !! tully_2, or tully_3 in the input file.
  !! It does not affect the calculations, but it appears in the output as summary of 
  !! of the computation.
  integer :: n_dof=1 !< The number of degrees of freedom, from 1 to 3 available.
  integer :: x_points=4001 !< The number of grid points along the x-direction, usually 
  !! determined by the potential files, but it is necessary to indicate it in the input file.
  integer :: y_points=1 !< Same as x_points, but for the y-direction.
  integer :: z_points=1 !< Same as x_points, but for the z-direction.
  integer :: nstates=2 !< Number of electronic states included.
  character(len=100) :: el_basis="adiabatic" !< Both adiabatic and diabatic basis can be
  !! specified. The program works with adiabatic electronic states, therefore it should
  !! be indicated in the input file if the diabatic Hamiltonian is given, such that
  !! the diagonalization routine is invoked.
  character(len=1) :: dia_to_ad="n" !< Request if the electronic populations should
  !! written in output also in the diabatic basis.
  real(kind=dp),allocatable :: transformation_matrix(:,:,:,:,:)!< Transformation matrix from
  !! diabatic to adiabatic basis.
  integer :: npairs=1 !< The number of pairs of electronic states, without double counting.
  namelist /system/ model_system,n_dof,x_points, &
    y_points,z_points,nstates,el_basis,dia_to_ad

  real(kind=dp),allocatable :: BOpes(:,:,:,:)
  real(kind=dp),allocatable :: na_coup(:,:,:,:,:,:)
  real(kind=dp),allocatable :: x_grid(:)
  real(kind=dp),allocatable :: y_grid(:)
  real(kind=dp),allocatable :: z_grid(:)
  real(kind=dp),allocatable :: Hel(:,:,:,:,:)

  character(len=100) :: algorithm='CTMQC'
  real(kind=dp) :: dt=0.1_dp
  real(kind=edp) :: final_time=0._edp
  real(kind=dp) :: r_init(100)
  real(kind=dp) :: k_init(100)
  real(kind=dp) :: mass_input(100)
  real(kind=dp) :: sigma_init(100)
  integer :: ntraj=100
  integer :: nsteps=100
  integer(kind=di) :: dump=1
  integer :: initial_BOstate=0
  integer :: initial_DIAstate=1
  real(kind=dp) :: viscosity
  real(kind=dp) :: temperature
  namelist /dynamics/ algorithm,final_time,dt,dump,initial_BOstate, &
    initial_DIAstate,ntraj,r_init,k_init,sigma_init,mass_input, &
    viscosity,temperature

  real(kind=dp),allocatable :: r0(:)
  real(kind=dp),allocatable :: r02(:)
  real(kind=dp),allocatable :: k0(:)
  real(kind=dp),allocatable :: mass(:)
  real(kind=dp),allocatable :: sigma(:)
  real(kind=dp),allocatable :: BOforce(:,:,:)
  real(kind=dp),allocatable :: coup(:,:,:,:)
  real(kind=dp),allocatable :: BOenergy(:,:)
  real(kind=dp),allocatable :: BO_pop(:)
  real(kind=dp),allocatable :: BO_coh(:)
  real(kind=dp),allocatable :: initial_positions(:,:)
  real(kind=dp),allocatable :: initial_momenta(:,:)
  real(kind=dp),allocatable :: weight(:)
  real(kind=dp),allocatable :: store_gamma(:,:)
  real(kind=dp),allocatable :: tdpes(:)
  real(kind=dp),allocatable :: density(:)
  real(kind=dp),allocatable :: vv_param(:,:)

  character(len=100) :: path_to_potentials="./tests/"
  character(len=100) :: positions_file=""
  character(len=100) :: momenta_file=""
  namelist /external_files/ path_to_potentials, &
    positions_file,momenta_file


end module variables

