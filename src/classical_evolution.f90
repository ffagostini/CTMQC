module classical_evolution
  use variables
  use kinds
  use wigner_distribution
  implicit none

  real(kind=dp) :: mydgaugeterm

  contains

  subroutine velocity_verlet(x,v,force)

    real(kind=dp),intent(inout) :: x(n_dof),v(n_dof)
    real(kind=dp),intent(in) :: force(n_dof)
    real(kind=dp) :: old_x(n_dof),old_v(n_dof),old_old_v(n_dof)

    old_x=x
    old_old_v=v

    old_v=VV_velocity(old_old_v,force)
    x=VV_position(old_x,old_v)
    v=VV_velocity(old_v,force)

  end subroutine velocity_verlet


  function VV_position(x,v) result(my_x)

    real(kind=dp),intent(in) :: x(n_dof),v(n_dof)
    real(kind=dp) :: my_x(n_dof)

    my_x=x+dt*v

  end function VV_position


  function VV_velocity(v,force) result(my_v)

    real(kind=dp),intent(in) :: v(n_dof)
    real(kind=dp),intent(in)    :: force(n_dof)
    real(kind=dp) :: my_v(n_dof)

    if(model_system=="marcus") then
      !my_v=vv_param(:,1)*v+0.50_dp*dt*force/mass
      my_v=vv_param(:,1)*v+dt*force/mass
    else
      !my_v=v+0.50_dp*dt*force/mass
      my_v=v+dt*force/mass
    end if

  end function VV_velocity


  subroutine non_adiabatic_force(coeff,force,acc_force,k_li,trajlabel,velocity)

    integer,intent(in) :: trajlabel
    complex(kind=dp),intent(in) :: coeff(nstates)
    real(kind=dp),intent(in) :: acc_force(n_dof,nstates),&
      k_li(nstates,nstates),velocity(n_dof)
    real(kind=dp),intent(out) :: force(n_dof)
    integer :: i,j,i_dof,check,nrand
    complex(kind=dp),allocatable :: my_rho(:,:)
    real(kind=dp) :: noise,tmp_vector(2),variance
    real(kind=dp),allocatable :: xi(:)

    allocate(my_rho(nstates,nstates),stat=check)
    if(check/=0) print*,'error 1 my_rho in non_adiabatic_force'

    force=0.0_dp
    do i=1,nstates
      do j=1,nstates
        my_rho(i,j)=conjg(coeff(i))*coeff(j)
      end do
    end do

    !Ehrenfest like force terms
    do i_dof=1,n_dof
      force(i_dof)=0.0_dp
      do i=1,nstates
        force(i_dof)=force(i_dof)+ &
          BOforce(trajlabel,i_dof,i)*real(my_rho(i,i),kind=dp)
      end do
      do i=1,nstates
        do j=i+1,nstates
          force(i_dof)=force(i_dof)- &
            2.0_dp*real(my_rho(i,j),kind=dp)* &
            (BOenergy(trajlabel,j)-BOenergy(trajlabel,i))* &
            coup(trajlabel,i_dof,i,j)
        end do
      end do
    end do

    !Force with the quantum momentum
    if(algorithm=="CTMQC") then
      do i_dof=1,n_dof
        do i=1,nstates
          do j=1,nstates
            force(i_dof)=force(i_dof)+ &
              0.5_dp*k_li(i,j)*(acc_force(i_dof,j)-acc_force(i_dof,i))* &
              (real(my_rho(i,i),kind=dp)*real(my_rho(j,j),kind=dp))
          end do
        end do
      end do
    end if

    !Random and viscous force
    if(model_system=="marcus") then
      force=force*vv_param(:,2)
      nrand=2
      allocate(xi(nrand),stat=check)
      if(check/=0) print*,'error xi'
      call random_number(xi)
      do i_dof=1,n_dof
        variance=1.0_dp!sqrt(2.0_dp*viscosity*mass(i_dof)*kB*temperature/dt)
        tmp_vector=gaussian_distribution(xi,nrand,variance,0.0_dp,nrand)
        noise=tmp_vector(1)
        force(i_dof)=force(i_dof)- &
          !viscosity*mass(i_dof)*velocity(i_dof)+
          vv_param(i_dof,3)*noise
      end do
      deallocate(xi)
    end if

    deallocate(my_rho,stat=check)
    if(check/=0) print*,'error 2 my_rho in non_adiabatic_force'

  end subroutine non_adiabatic_force

end module classical_evolution





