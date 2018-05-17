module wigner_distribution
  use variables
  use tools
  implicit none

  contains

  subroutine initial_conditions

    real(kind=dp),allocatable :: xi(:)
    integer :: check,nrand,i,j,ios

    call generate_random_seed

    if(trim(positions_file)=="") then
      nrand=ntraj
      if(mod(ntraj,2)/=0) nrand=ntraj+1
      allocate(xi(nrand),stat=check)
      if(check/=0) print*,'error xi'
      do i=1,n_dof
        call random_number(xi)
        initial_positions(:,i)= &
          gaussian_distribution(xi,nrand,sigma(i)/sqrt(2.0_dp),r0(i),ntraj)
      end do
      deallocate(xi)
    else
      open(26,file=trim(positions_file),status="old", &
        form="formatted",action="read",iostat=ios)
      if(ios/=0) print*,'error opening file of positions'
      do i=1,ntraj
        read(26,*) initial_positions(i,:)
      end do
      close(26)
    end if

    if(trim(momenta_file)=="") then
      nrand=ntraj
      if(mod(ntraj,2)/=0) nrand=ntraj+1
      allocate(xi(nrand),stat=check)
      if(check/=0) print*,'error xi'
      do i=1,n_dof
        call random_number(xi)
        initial_momenta(:,i)= &
          gaussian_distribution(xi,nrand,hbar/sigma(i)/sqrt(2.0_dp),k0(i),ntraj)
      end do
      deallocate(xi)
    else
      open(26,file=trim(momenta_file),status="old", &
        form="formatted",action="read",iostat=ios)
      if(ios/=0) print*,'error opening file of momenta'
      do i=1,ntraj
        read(26,*) initial_momenta(i,:)
        initial_momenta(i,:)=initial_momenta(i,:)*mass
      end do
      close(26)
    end if

    !output initial conditions
    open(23,file='./output/initial_conditions.dat')
    do i=1,ntraj
      write(23,'(2000f14.4)') initial_positions(i,:),initial_momenta(i,:)
    end do
    close(23)

  end subroutine initial_conditions
  

  function gaussian_distribution(xi,nrand,var,x0,my_nrand) result(y)

    integer,intent(in) :: nrand,my_nrand
    real(kind=dp),intent(in) :: xi(nrand),var,x0
    real(kind=dp) :: y_tmp(nrand),y(my_nrand)
    integer :: i

    do i=1,nrand,2
      y_tmp(i)=sqrt(-2.0_dp*log(xi(i+1)))*cos(2.0_dp*PI*xi(i))
      y_tmp(i+1)=sqrt(-2.0_dp*log(xi(i+1)))*sin(2.0_dp*PI*xi(i))
    end do

    y_tmp=y_tmp*var+x0

    if (nrand/=my_nrand) then
      do i=1,my_nrand
        y(i)=y_tmp(i)
      end do
    else
      y=y_tmp
    end if

  end function gaussian_distribution

end module wigner_distribution
