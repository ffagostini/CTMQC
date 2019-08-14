module tools
  use variables
  use kinds
  implicit none

  contains

  subroutine initialize_system_vars

    integer :: check

    PI=acos(-1.0_dp)
    Im_unit=(0.0_dp,1.0_dp)

    allocate(x_grid(x_points),stat=check)
    if(check/=0) print*,'error allocation x_grid'
    allocate(y_grid(y_points),stat=check)
    if(check/=0) print*,'error allocation y_grid'
    allocate(z_grid(z_points),stat=check)
    if(check/=0) print*,'error allocation z_grid'

    allocate(BOpes(x_points,y_points,z_points,nstates),stat=check)
    if(check/=0) print*,'error allocation BOpes'
    allocate(na_coup(n_dof,x_points,y_points,z_points,nstates,nstates),stat=check)
    if(check/=0) print*,'error allocation na_coup'

    allocate(Hel(nstates,nstates,x_points,y_points,z_points),stat=check)
    if(check/=0) print*,'error allocation Hel'
    allocate(transformation_matrix(nstates,nstates,x_points,y_points,z_points),stat=check)
    if(check/=0) print*,'error allocation transformation_matrix'

    allocate(density(x_points),stat=check)
    if(check/=0) print*,'error density'

  end subroutine initialize_system_vars


  subroutine initialize_dynamics_vars

    integer :: check

    allocate(r0(n_dof),stat=check)
    if(check/=0) print*,'error allocatin r0'
    allocate(r02(n_dof),stat=check)
    if(check/=0) print*,'error allocatin r02'
    allocate(k0(n_dof),stat=check)
    if(check/=0) print*,'error allocation k0'
    allocate(mass(n_dof),stat=check)
    if(check/=0) print*,'error allocation mass'
    allocate(sigma(n_dof),stat=check)
    if(check/=0) print*,'error allocation sigma'
    allocate(vv_param(n_dof,3),stat=check)
    if(check/=0) print*,'error allocation vv_param'

  end subroutine initialize_dynamics_vars


  subroutine initialize_trajectory_vars
  
    integer :: check

    allocate(BOenergy(ntraj,nstates),stat=check)
    if(check/=0) print*,'error BOenergy'
    allocate(BOforce(ntraj,n_dof,nstates),stat=check)
    if(check/=0) print*,'error BOforce'
    allocate(coup(ntraj,n_dof,nstates,nstates),stat=check)
    if(check/=0) print*,'error coup'

    allocate(BO_pop(nstates),stat=check)
    if(check/=0) print*,'error BO_pop'
    allocate(BO_coh(npairs),stat=check)
    if(check/=0) print*,'error BO_coh'

    allocate(initial_positions(ntraj,n_dof),stat=check)
    if(check/=0) print*,'error initial_positions'
    allocate(initial_momenta(ntraj,n_dof),stat=check)
    if(check/=0) print*,'error initial_momenta'

    allocate(store_gamma(n_dof,ntraj),stat=check)
    if(check/=0) print*,'error store_gamma'

    allocate(tdpes(ntraj),stat=check)
    if(check/=0) print*,'error tdpes'

  end subroutine initialize_trajectory_vars


  subroutine finalize

    integer :: check

    deallocate(coup,stat=check)
    if(check/=0) print*,'error coup'
    deallocate(BOenergy,stat=check)
    if(check/=0) print*,'error BOenergy'
    deallocate(BOforce,stat=check)
    if(check/=0) print*,'error BOforce'

    deallocate(initial_positions,stat=check)
    if(check/=0) print*,'error initial_positions'
    deallocate(initial_momenta,stat=check)
    if(check/=0) print*,'error initial_momenta'
    deallocate(sigma,stat=check)
    if(check/=0) print*,'error sigma'
    deallocate(store_gamma,stat=check)
    if(check/=0) print*,'error store_gamma'

    deallocate(tdpes,stat=check)
    if(check/=0) print*,'error tdpes'
    deallocate(density,stat=check)
    if(check/=0) print*,'error density'

    deallocate(x_grid,stat=check)
    if(check/=0) print*,'error x_grid'
    deallocate(y_grid,stat=check)
    if(check/=0) print*,'error y_grid'
    deallocate(z_grid,stat=check)
    if(check/=0) print*,'error z_grid'

    deallocate(BOpes,stat=check)
    if(check/=0) print*,'error BOpes'
    deallocate(na_coup,stat=check)
    if(check/=0) print*,'error na_coup'

    deallocate(Hel,stat=check)
    if(check/=0) print*,'error Hel'
    deallocate(transformation_matrix)
    if(check/=0) print*,'error transformation_matrix'

    deallocate(BO_pop,stat=check)
    if(check/=0) print*,'error BO_pop'
    deallocate(BO_coh,stat=check)
    if(check/=0) print*,'error BO_coh'

    deallocate(r0,stat=check)
    if(check/=0) print*,'error r0'
    deallocate(r02,stat=check)
    if(check/=0) print*,'error r02'
    deallocate(k0,stat=check)
    if(check/=0) print*,'error k0'
    deallocate(mass,stat=check)
    if(check/=0) print*,'error mass'
    deallocate(vv_param,stat=check)
    if(check/=0) print*,'error allocation vv_param'

  end subroutine finalize


  subroutine generate_random_seed

    integer,allocatable :: seed(:)
    integer :: n,i,clock

    call random_seed
    call random_seed(size=n)

    allocate(seed(n))
    !call system_clock(count=clock)
    !seed=clock+37*(/(i-1,i=1,n)/)
    !seed=1
    do i = 1,n
      seed(i) = 12345 + i
    end do
    call random_seed(PUT=seed)
          
    deallocate(seed)

  end subroutine generate_random_seed


  subroutine integrator_parameters

    integer :: i_dof

    do i_dof=1,n_dof
      !vv_param(i_dof,1)=exp(-viscosity*(dt/2.0_dp)/(mass(i_dof)))
      !vv_param(i_dof,2)=(mass(i_dof)/viscosity)* &
      !  (1._dp-vv_param(i_dof,1)) / (dt/2.0_dp)
      !vv_param(i_dof,3)=sqrt(mass(i_dof)*kB*temperature)* &
      !  sqrt(1._dp-(vv_param(i_dof,1))**2) / (dt/2.0_dp)
      vv_param(i_dof,1)=exp(-viscosity*(dt/1.0_dp)/(mass(i_dof)))
      vv_param(i_dof,2)=(mass(i_dof)/viscosity)* &
        (1._dp-vv_param(i_dof,1)) / (dt/1.0_dp)
      vv_param(i_dof,3)=sqrt(mass(i_dof)*kB*temperature)* &
        sqrt(1._dp-(vv_param(i_dof,1))**2) / (dt/1.0_dp)
    end do

  end subroutine integrator_parameters

  subroutine change_basis

    use linear_algebra
    use read_diabatic_properties

    call read_diabatic_hamiltonian

    if(n_dof==1) call diabatic_to_adiabatic_1d(Hel(:,:,:,1,1))
    if(n_dof==2) call diabatic_to_adiabatic_2d(Hel(:,:,:,:,1))
    if(n_dof==3) call diabatic_to_adiabatic_3d(Hel(:,:,:,:,:))

  end subroutine change_basis


  subroutine smoothing(funct,x,dim,var)

    integer,intent(in) :: dim
    real(kind=dp),intent(in) :: x(dim)
    real(kind=dp),intent(inout) :: funct(dim)
    real(kind=dp),intent(in),optional :: var
    real(kind=dp) :: delta,my_var
    real(kind=dp),allocatable :: f_in(:),f_out(:),g(:)
    integer :: xi,xj,dim_3sigma

    delta=x(2)-x(1)
    if(.not.present(var)) then
      my_var=0.10_dp
    else
      my_var=var
    end if

    dim_3sigma=1+anint(2.0_dp*3.0_dp*my_var/delta)

    allocate(f_in(dim),f_out(dim),g(-(dim_3sigma-1)/2:(dim_3sigma-1)/2))

    f_in=funct
    f_out=0.0_dp

    do xi=1,dim
      do xj=-(dim_3sigma-1)/2,(dim_3sigma-1)/2
        if(xi+xj>=1 .and. xi+xj<=dim) then
          g(xj)=sqrt(1.0_dp/2.0_dp/PI/my_var**2)* &
            exp(-(x(xi+xj)-x(xi))**2/2.0_dp/my_var**2)
        else
          g(xj)=0.0_dp
        end if
      end do
      do xj=-(dim_3sigma-1)/2,(dim_3sigma-1)/2
        if(xi-xj>=1 .and. xi-xj<=dim) &
          f_out(xi)=f_out(xi)+delta*f_in(xi-xj)*g(xj)
      end do
    end do

    funct=f_out

    deallocate(f_in,f_out,g)

  end subroutine smoothing

end module tools
