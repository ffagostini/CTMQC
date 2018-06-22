module coefficients_evolution
  use variables
  use kinds
  implicit none

  contains

  subroutine RK4_coeff(v,coeff,k_li,trajlabel)

    integer,intent(in) :: trajlabel
    real(kind=dp),intent(in) :: v(n_dof),k_li(nstates,nstates)
    complex(kind=dp),intent(inout) :: coeff(nstates)
    integer :: i
    complex(kind=dp) :: k1,k2,k3,k4,&
      variation(nstates),kfunction,my_coeff(nstates)
    real(kind=dp) :: normalization

    variation=cmplx(0.0_dp,0.0_dp)
    my_coeff=coeff

    statesloop: do i=1,nstates
      kfunction=cmplx(0.0_dp,0.0_dp)
      k1=dt*cdot(i,kfunction,v,my_coeff,k_li,trajlabel)
      kfunction=0.50_dp*k1
      k2=dt*cdot(i,kfunction,v,my_coeff,k_li,trajlabel)
      kfunction=0.50_dp*(-1.0_dp+sqrt(2.0_dp))*k1+&
        (1.0_dp-0.5_dp*sqrt(2.0_dp))*k2
      k3=dt*cdot(i,kfunction,v,my_coeff,k_li,trajlabel)
      kfunction=-0.50_dp*sqrt(2.0_dp)*k2+&
        (1.0_dp+0.5_dp*sqrt(2.0_dp))*k3
      k4=dt*cdot(i,kfunction,v,my_coeff,k_li,trajlabel)
      variation(i)=(k1+(2.0_dp-sqrt(2.0_dp))*k2+&
        (2.0_dp+sqrt(2.0_dp))*k3+k4)/6.0_dp
    end do statesloop

    coeff=my_coeff+variation

    normalization=0.0_dp
    do i=1,nstates
      normalization=normalization+(abs(coeff(i)))**2
    end do
    !write(6,*) normalization
    !!if(abs(normalization-1.0_dp)>zero) &
    coeff=coeff/sqrt(normalization)

  end subroutine RK4_coeff


  function cdot(state,kfunction,v,coeff,k_li,trajlabel)

    integer,intent(in) :: state,trajlabel
    real(kind=dp),intent(in) :: v(n_dof),k_li(nstates,nstates)
    complex(kind=dp),intent(in) :: coeff(nstates),kfunction
    complex(kind=dp) :: cdot,nonadiabatic_sum
    complex(kind=dp),allocatable :: my_coeff(:)
    integer :: i

    allocate(my_coeff(nstates))

    do i=1,nstates
      if(i==state) then
        my_coeff(i)=coeff(i)+kfunction
      else
        my_coeff(i)=coeff(i)
      end if
    end do

    nonadiabatic_sum=cmplx(0.0_dp,0.0_dp)
    do i=1,nstates
      if(i/=state) nonadiabatic_sum=nonadiabatic_sum+ &
          my_coeff(i)*dot_product(coup(trajlabel,:,state,i),v(:))
    end do

    cdot=-Im_unit*my_coeff(state)*BOenergy(trajlabel,state)-nonadiabatic_sum

    do i=1,nstates
      if(i/=state) then
        cdot=cdot+0.250_dp*(k_li(i,state)-k_li(state,i))* &
          (abs(my_coeff(i)))**2*my_coeff(state)
      end if
    end do

    deallocate(my_coeff)

  end function cdot


end module coefficients_evolution
