module output
  use variables
  use kinds
  implicit none

  contains

  subroutine plot(BOsigma,Rcl,Vcl,time)

    integer,intent(in) :: time
    real(kind=dp),intent(inout) :: Rcl(ntraj,n_dof),Vcl(ntraj,n_dof)
    complex(kind=dp),intent(in) :: BOsigma(ntraj,nstates,nstates)
    integer :: i,j,itraj,index_ij

    if(n_dof==1) then
      call plot_coefficients(BOsigma,Rcl(:,1),time)
      call plot_histogram(Rcl(:,1),time)
      call plot_density(Rcl(:,1),time)
    end if

    do itraj=1,ntraj
      call compute_energy(BOsigma(itraj,:,:),BOenergy(itraj,:),itraj)
    end do
    call plot_R_P_E(Rcl,Vcl,time)

    if(time==0) call initialize_output

    BO_pop=0.0_dp
    do itraj=1,ntraj
      do i=1,nstates
        BO_pop(i)=BO_pop(i)+real(BOsigma(itraj,i,i),kind=dp)
      end do
    end do
    BO_pop=BO_pop/dble(ntraj)

    BO_coh=0.0_dp
    index_ij=0
    do i=1,nstates
      do j=i+1,nstates
        index_ij=index_ij+1
        do itraj=1,ntraj
          BO_coh(index_ij)=BO_coh(index_ij)+ &
            (real(BOsigma(itraj,i,i),kind=dp))* &
            (real(BOsigma(itraj,j,j),kind=dp))
        end do
      end do
    end do
    BO_coh=BO_coh/dble(ntraj)

    write(88,'(f14.4,100f14.8)') dble(time)*dt,BO_coh
    write(89,'(f14.4,100f14.8)') dble(time)*dt,BO_pop

    if(time==nsteps) call finalize_output

  end subroutine plot


  subroutine plot_coefficients(BOsigma,Rcl,time)

    integer,intent(in) :: time
    real(kind=dp),intent(in) :: Rcl(ntraj)
    complex(kind=dp),intent(in) :: BOsigma(ntraj,nstates,nstates)
    character(len=3) :: idx
    character(len=100) :: filename
    integer :: itraj,ios

    write(idx,'(i3.3)') time/dump
    filename="./output/coeff/coeff."//trim(idx)//".dat"

    open(128,file=trim(filename),status="replace", &
      form="formatted",action="write",iostat=ios)
    if(ios/=0) print*,'error opening coefficients file'
    write(128,*) "#Postion, Coefficients: Real part and Imaginary part"

    do itraj=1,ntraj
      write(128,'(100f14.8)') Rcl(itraj), &
        (real(BOsigma(itraj,:,:))),(aimag(BOsigma(itraj,:,:)))
    end do

    close(128)

  end subroutine plot_coefficients


  subroutine plot_R_P_E(Rcl,Vcl,time)

    integer,intent(in) :: time
    real(kind=dp),intent(in) :: Rcl(ntraj,n_dof), &
      Vcl(ntraj,n_dof)
    character(len=3) :: idx
    character(len=100) :: filename
    integer :: itraj,ios

    write(idx,'(i3.3)') time/dump
    filename="./output/trajectories/RPE."//trim(idx)//".dat"

    open(128,file=trim(filename),status="replace", &
      form="formatted",action="write",iostat=ios)
    if(ios/=0) print*,'error opening RPE file'
    write(128,*) "#Postions, Momenta, TDPES (GI part)"

    do itraj=1,ntraj
      write(128,'(10f14.8)') Rcl(itraj,:), &
        mass*Vcl(itraj,:),tdpes(itraj)
    end do

    close(128)

  end subroutine plot_R_P_E


  subroutine plot_histogram(Rcl,time)

    integer,intent(in) :: time
    real(kind=dp),intent(in) :: Rcl(ntraj)
    character(len=3) :: idx
    character(len=100) :: filename
    integer :: itraj,nbeads,ibead,ios
    real(kind=dp) :: dpos
    real(kind=dp),allocatable :: histo(:)

    nbeads=200
    allocate(histo(nbeads))

    dpos=(x_grid(x_points)-x_grid(1))/dble(nbeads)

    histo=0
    do itraj=1,ntraj
      ibead=1+anint((Rcl(itraj)-x_grid(1))/dpos)
      histo(ibead)=histo(ibead)+1.0_dp
    end do

    write(idx,'(i3.3)') time/dump
    filename="./output/histo/histo."//trim(idx)//".dat"

    open(128,file=trim(filename),status="replace", &
      form="formatted",action="write",iostat=ios)
    if(ios/=0) print*,'error opening histogram file'
    write(128,*) "#Postion, histogram"

    do ibead=1,nbeads
      write(128,'(9f14.8)') x_grid(1)+dble(ibead-1)*dpos, &
        dble(histo(ibead))/(dble(sum(histo))*dpos)
    end do

    close(128)

    deallocate(histo)

  end subroutine plot_histogram


  subroutine plot_density(Rcl,time)
    use tools

    integer,intent(in) :: time
    real(kind=dp),intent(in) :: Rcl(ntraj)
    character(len=3) :: idx
    character(len=100) :: filename
    integer :: itraj,x,ios
    real(kind=dp),allocatable :: my_gamma(:)

    allocate(my_gamma(ntraj))
    density=0.0_dp

    if(time==0) then
      my_gamma=minval(sigma)
    else
      my_gamma=store_gamma(1,:)
    end if

    do x=1,x_points
      do itraj=1,ntraj
        density(x)=density(x)+ &
          1.0_dp/sqrt(PI*(my_gamma(itraj))**2)* &
          exp(-(x_grid(x)-Rcl(itraj))**2/(my_gamma(itraj))**2)/dble(ntraj)
      end do
    end do

    write(idx,'(i3.3)') time/dump
    filename="./output/density/density."//trim(idx)//".dat"
    open(128,file=trim(filename),status="replace", &
    form="formatted",action="write",iostat=ios)
    if(ios/=0) print*,'error opening density file'
    write(128,*) "#Postion, density"
    do x=1,x_points
      write(128,*) x_grid(x),density(x)
    end do
    close(128)

    call smoothing(density,x_grid,x_points,var=0.5D0)

    write(idx,'(i3.3)') time/dump
    filename="./output/density/smooth_density."//trim(idx)//".dat"
    open(128,file=trim(filename),status="replace", &
      form="formatted",action="write",iostat=ios)
    if(ios/=0) print*,'error opening density smooth file'
    write(128,*) "#Postion, smooth density"
    do x=1,x_points
      write(128,*) x_grid(x),density(x)
    end do
    close(128)

    deallocate(my_gamma)

  end subroutine plot_density

  subroutine diabatic_output(Rcl,BOcoeff,time)
    use electronic_problem

    real(kind=dp),intent(in) :: Rcl(ntraj,n_dof)
    complex(kind=dp),intent(in) :: BOcoeff(ntraj,nstates)
    integer,intent(in) :: time
    integer :: check,i,itraj,istar(3)
    real(kind=dp),allocatable :: DIA_pop(:)
    complex(kind=dp),allocatable :: tmp_vector(:)

    if(time==0) call initialize_diabatic_output

    allocate(DIA_pop(nstates),stat=check)
    if(check/=0) print*,'error allocation DIA_pop'
    allocate(tmp_vector(nstates),stat=check)
    if(check/=0) print*,'error allocation tmp_vector'

    DIA_pop=0.0_dp

    do itraj=1,ntraj
      call locate_x(Rcl(itraj,:),istar)
      tmp_vector=matmul(transformation_matrix(:,:, &
        istar(1),istar(2),istar(3)),BOcoeff(itraj,:))
      do i=1,nstates
        DIA_pop(i)=DIA_pop(i)+real(conjg(tmp_vector(i))*tmp_vector(i))
      end do
    end do
    DIA_pop=DIA_pop/dble(ntraj)

    write(87,'(f14.4,100f14.8)') dble(time)*dt,DIA_pop

    deallocate(DIA_pop,stat=check)
    if(check/=0) print*,'error DIA_pop'
    deallocate(tmp_vector,stat=check)
    if(check/=0) print*,'error tmp_vector'

    if(time==nsteps) call finalize_diabatic_output

  end subroutine diabatic_output


  subroutine compute_energy(my_rho,e_BO,trajlabel)

    integer,intent(in) :: trajlabel
    complex(kind=dp),intent(in) :: my_rho(nstates,nstates)
    real(kind=dp),intent(in) :: e_BO(nstates)
    integer :: i

    tdpes(trajlabel)=0.0_dp

    do i=1,nstates
      tdpes(trajlabel)=tdpes(trajlabel)+ &
        real(my_rho(i,i),kind=dp)*e_BO(i)
    end do

  end subroutine compute_energy


  subroutine initialize_output

    integer :: ios

    open(88,file="./output/BO_coherences.dat",status="replace", &
      form="formatted",action="write",iostat=ios)
    if(ios/=0) print*,'error opening BO_coherences.dat'
    write(88,*) "#Time, Coherences"

    open(89,file="./output/BO_population.dat",status="replace", &
      form="formatted",action="write",iostat=ios)
    if(ios/=0) print*,'error opening BO_population.dat'
    write(89,*) "#Time, Populations"

  end subroutine initialize_output


  subroutine finalize_output

    close(89)
    close(88)

  end subroutine finalize_output


  subroutine initialize_diabatic_output

    integer :: ios

    open(87,file="./output/DIA_population.dat",status="replace", &
      form="formatted",action="write",iostat=ios)
    if(ios/=0) print*,'error opening DIA_population.dat'
    write(87,*) "#Time, Diabatic population"

  end subroutine initialize_diabatic_output


  subroutine finalize_diabatic_output

    close(87)

  end subroutine finalize_diabatic_output


end module output
