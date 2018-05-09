module electronic_problem
  use variables
  use kinds
  implicit none

  !---------------------------------------------------------------------------
  !> @author
  !> Federica Agostini, University Paris-Sud.
  !
  ! DESCRIPTION:
  !> Electronic structure calculation at the position of the trajectory
  ! REVISION HISTORY:
  ! Version of March 22nd 2018
  !
  !> @param[in] x Position at which adiabatic potential energies, forces and
  !! nonadiabatic couplings are required.
  !> @param[in] trajlabel Labels the trajectory at the position x.
  !> @param istar 3-dimensional vector of integers indicating the position of x in the
  !! energy and nonadiabatic couplings arrays.
  !> @param i_dof Labels the degree of freedom.
  !> @return The arrays BOenergy, BOforce, coup will contain the values of energies,
  !! forces, and nonadiabatic couplings at the position x of the trajectory trajlabel.
  !---------------------------------------------------------------------------

  contains

  subroutine BOproblem(x,trajlabel)

    integer,intent(in) :: trajlabel
    real(kind=dp),intent(in) :: x(n_dof)
    integer :: istar(3)
    integer :: i_dof

    call locate_x(x,istar)

    BOenergy(trajlabel,:)=BOpes(istar(1),istar(2),istar(3),:)
    do i_dof=1,n_dof
      coup(trajlabel,i_dof,:,:)=na_coup(i_dof,istar(1),istar(2),istar(3),:,:)
    end do
    BOforce(trajlabel,:,:)=hf_forces(istar)

  end subroutine BOproblem

  !< @param[in] x Position at which adiabatic potential energies, forces and
  !! nonadiabatic couplings are required.
  !! @param[in,out] istar 3-dimensional vector of integers indicating the position of x in the
  !! energy and nonadiabatic couplings arrays.
  !! @param i_dof Labels the degree of freedom.

  subroutine locate_x(x,istar)

    real(kind=dp),intent(in) :: x(n_dof)
    integer,intent(out) :: istar(3)
    integer :: i_dof

    istar=1

    do i_dof=1,n_dof
      if(i_dof==1) &
          istar(i_dof)=1+anint((x(i_dof)-x_grid(1))/(x_grid(2)-x_grid(1)))
      if(i_dof==2) &
        istar(i_dof)=1+anint((x(i_dof)-y_grid(1))/(y_grid(2)-y_grid(1)))
      if(i_dof==3) &
        istar(i_dof)=1+anint((x(i_dof)-z_grid(1))/(z_grid(2)-z_grid(1)))
    end do

  end subroutine locate_x


  function hf_forces(istar)

    integer,intent(in) :: istar(3)
    real(kind=dp) :: hf_forces(n_dof,nstates),c1,c2,dx,dy,dz
    integer :: x,y,z,n

    c1=1.0_dp/12.0_dp
    c2=8.0_dp/12.0_dp
    hf_forces=0.0_dp

    x=istar(1)
    y=istar(2)
    z=istar(3)

    if(n_dof==1) then
      dx=x_grid(2)-x_grid(1)
      if(x>=3 .and. x<=x_points-2) then
        do n=1,nstates
          hf_forces(1,n)=-(c1*(BOpes(x-2,y,z,n)-BOpes(x+2,y,z,n))+ &
            c2*(BOpes(x+1,y,z,n)-BOpes(x-1,y,z,n)))/dx
        end do
      end if
    end if

    if(n_dof==2) then
      dx=x_grid(2)-x_grid(1)
      dy=y_grid(2)-y_grid(1)
      if(x>=3 .and. x<=x_points-2 .and. y>=3 .and. y<=y_points-2) then
        hf_forces(1,:)=-(c1*(BOpes(x-2,y,z,:)-BOpes(x+2,y,z,:))+ &
          c2*(BOpes(x+1,y,z,:)-BOpes(x-1,y,z,:)))/dx
        hf_forces(2,:)=-(c1*(BOpes(x,y-2,z,:)-BOpes(x,y+2,z,:))+ &
          c2*(BOpes(x,y+1,z,:)-BOpes(x,y-1,z,:)))/dy
      end if
    end if

    if(n_dof==3) then
      dx=x_grid(2)-x_grid(1)
      dy=y_grid(2)-y_grid(1)
      dz=z_grid(2)-z_grid(1)
      if(x>=3 .and. x<=x_points-2 .and. y>=3 .and. y<=y_points-2) then
        hf_forces(1,:)=-(c1*(BOpes(x-2,y,z,:)-BOpes(x+2,y,z,:))+ &
          c2*(BOpes(x+1,y,z,:)-BOpes(x-1,y,z,:)))/dx
        hf_forces(2,:)=-(c1*(BOpes(x,y-2,z,:)-BOpes(x,y+2,z,:))+ &
          c2*(BOpes(x,y+1,z,:)-BOpes(x,y-1,z,:)))/dy
        hf_forces(3,:)=-(c1*(BOpes(x,y,z-2,:)-BOpes(x,y,z+2,:))+ &
          c2*(BOpes(x,y,z+1,:)-BOpes(x,y,z-1,:)))/dz
      end if
    end if

  end function hf_forces


end module electronic_problem
