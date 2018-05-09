module linear_algebra
  use variables
  use kinds
  implicit none

  contains

  subroutine diabatic_to_adiabatic_1d(Hel)

    real(kind=8),intent(in) :: Hel(nstates,nstates,x_points)
    real(kind=8),allocatable :: Ebo(:,:),U(:,:,:)
    integer :: x,i,ioerr,lwork,dim_work,liwork,dim_iwork
    integer,allocatable :: iwork(:)
    real(kind=8),allocatable :: work(:)
    character(len=500) :: idx,filename

    allocate(Ebo(nstates,x_points))
    allocate(U(nstates,nstates,x_points))

    U(:,:,1)=Hel(:,:,1)
    dim_work=1
    dim_iwork=1
    allocate(work(dim_work),iwork(dim_iwork))
    lwork=-1
    liwork=-1
    call dsyevd('V','U',nstates,U(:,:,1),nstates,Ebo(:,1), &
    work,lwork,iwork,liwork,ioerr)
    if(ioerr/=0) print*,'error diagonalizing 1'
    dim_work=int(work(1))
    dim_iwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(dim_work),iwork(dim_iwork))
    do x=1,x_points
      U(:,:,x)=Hel(:,:,x)
      lwork=1+6*nstates+2*nstates**2
      liwork=3+5*nstates
      call dsyevd('V','U',nstates,U(:,:,x),nstates,Ebo(:,x), &
        work,lwork,iwork,liwork,ioerr)
      if(ioerr/=0) print*,'error diagonalizing 2'
      if(x>1) call check_overlap(U(:,:,x),U(:,:,x-1))
    end do
    deallocate(work,iwork)

    call non_adiabatic_couplings_1d(U)

    do i=1,nstates
      write(idx,"(i1.1)") i
      filename=trim(idx)//"_bopes.dat"
      open(300+i,file=trim(filename),status='replace',form='formatted', &
        action='write',iostat=ioerr)
      if(ioerr/=0) print*, 'error writing ',trim(filename)
      do x=1,x_points
        write(300+i,'(8f14.4)') Ebo(i,x),x_grid(x)
      end do
      close(300+i)
    end do

    deallocate(Ebo)
    deallocate(U)

  end subroutine diabatic_to_adiabatic_1d


  subroutine diabatic_to_adiabatic_2d(Hel)

    real(kind=8),intent(in) :: Hel(nstates,nstates,x_points,y_points)
    real(kind=8),allocatable :: Ebo(:,:,:),U(:,:,:,:)
    integer :: x,y,i,ioerr,lwork,dim_work,liwork,dim_iwork
    integer,allocatable :: iwork(:)
    real(kind=8),allocatable :: work(:)
    character(len=500) :: idx,filename

    allocate(Ebo(nstates,x_points,y_points))
    allocate(U(nstates,nstates,x_points,y_points))

    U(:,:,1,1)=Hel(:,:,1,1)
    dim_work=1
    dim_iwork=1
    allocate(work(dim_work),iwork(dim_iwork))
    lwork=-1
    liwork=-1
    call dsyevd('V','U',nstates,U(:,:,1,1),nstates,Ebo(:,1,1), &
    work,lwork,iwork,liwork,ioerr)
    if(ioerr/=0) print*,'error diagonalizing 1'
    dim_work=int(work(1))
    dim_iwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(dim_work),iwork(dim_iwork))

    do x=1,x_points
      do y=1,y_points
        U(:,:,x,y)=Hel(:,:,x,y)
        lwork=1+6*nstates+2*nstates**2
        liwork=3+5*nstates
        call dsyevd('V','U',nstates,U(:,:,x,y),nstates,Ebo(:,x,y), &
          work,lwork,iwork,liwork,ioerr)
        if(ioerr/=0) print*,'error diagonalizing 2'
        if(x>1) call check_overlap(U(:,:,x,y),U(:,:,x-1,y))
        if(y>1) call check_overlap(U(:,:,x,y),U(:,:,x,y-1))
      end do
    end do
    deallocate(work,iwork)

    call non_adiabatic_couplings_2d(U)

    do i=1,nstates
      write(idx,"(i1.1)") i
      filename=trim(idx)//"_bopes.dat"
      open(300+i,file=trim(filename),status='replace',form='formatted', &
        action='write',iostat=ioerr)
      if(ioerr/=0) print*, 'error writing ',trim(filename)
      do x=1,x_points
        do y=1,y_points
          write(300+i,'(8f14.4)') Ebo(i,x,y),x_grid(x),y_grid(y)
        end do
        write(300+i,*) ""
      end do
      close(300+i)
    end do

    deallocate(Ebo)
    deallocate(U)

  end subroutine diabatic_to_adiabatic_2d


  subroutine diabatic_to_adiabatic_3d(Hel)

    real(kind=8),intent(in),optional :: Hel(nstates,nstates,x_points,y_points,z_points)
    real(kind=8),allocatable :: Ebo(:,:,:,:),U(:,:,:,:,:)
    integer :: x,y,z,i,ioerr,lwork,dim_work,liwork,dim_iwork
    integer,allocatable :: iwork(:)
    real(kind=8),allocatable :: work(:)
    character(len=500) :: idx,filename

    allocate(Ebo(nstates,x_points,y_points,z_points))
    allocate(U(nstates,nstates,x_points,y_points,z_points))

    U(:,:,1,1,1)=Hel(:,:,1,1,1)
    dim_work=1
    dim_iwork=1
    allocate(work(dim_work),iwork(dim_iwork))
    lwork=-1
    liwork=-1
    call dsyevd('V','U',nstates,U(:,:,1,1,1),nstates,Ebo(:,1,1,1), &
    work,lwork,iwork,liwork,ioerr)
    if(ioerr/=0) print*,'error diagonalizing 1'
    dim_work=int(work(1))
    dim_iwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(dim_work),iwork(dim_iwork))

    do x=1,x_points
      do y=1,y_points
        do z=1,z_points
          U(:,:,x,y,z)=Hel(:,:,x,y,z)
          lwork=1+6*nstates+2*nstates**2
          liwork=3+5*nstates
          call dsyevd('V','U',nstates,U(:,:,x,y,z),nstates,Ebo(:,x,y,z), &
            work,lwork,iwork,liwork,ioerr)
          if(ioerr/=0) print*,'error diagonalizing 2'
          if(x>1) call check_overlap(U(:,:,x,y,z),U(:,:,x-1,y,z))
          if(y>1) call check_overlap(U(:,:,x,y,z),U(:,:,x,y-1,z))
          if(z>1) call check_overlap(U(:,:,x,y,z),U(:,:,x,y,z-1))
        end do
      end do
    end do
    deallocate(work,iwork)

    call non_adiabatic_couplings_3d(U)

    do i=1,nstates
      write(idx,"(i1.1)") i
      filename=trim(idx)//"_bopes.dat"
      open(300+i,file=trim(filename),status='replace',form='formatted', &
        action='write',iostat=ioerr)
      if(ioerr/=0) print*, 'error writing ',trim(filename)
      do x=1,x_points
        do y=1,y_points
          do z=1,z_points
            write(300+i,'(8f14.4)') Ebo(i,x,y,z),x_grid(x),y_grid(y),z_grid(z)
          end do
          write(300+i,*) ""
        end do
        write(300+i,*) ""
      end do
      close(300+i)
    end do

    deallocate(Ebo)
    deallocate(U)

  end subroutine diabatic_to_adiabatic_3d


  subroutine check_overlap(eigenv_out,eigenv_in)

    real(kind=8),intent(out) :: eigenv_out(nstates,nstates)
    real(kind=8),intent(in) :: eigenv_in(nstates,nstates)
    integer :: i

    do i=1,nstates
      if(eigenv_out(1,i)*eigenv_in(1,i)+eigenv_out(2,i)*eigenv_in(2,i)<0.0d0) &
        eigenv_out(:,i)=-eigenv_out(:,i)
    end do

  end subroutine check_overlap


  subroutine non_adiabatic_couplings_1d(U)

    real(kind=8),intent(in) :: U(nstates,nstates,x_points)
    real(kind=8),allocatable :: eigenv1(:,:),eigenv2(:,:), &
      deigenv2_dR(:,:),d_coup(:,:,:)
    real(kind=8) :: dx,c1,c2
    integer :: x,i,j,ioerr
    character(len=500) :: idx1,idx2,filename

    allocate(eigenv1(nstates,x_points),eigenv2(nstates,x_points), &
      deigenv2_dR(nstates,x_points),d_coup(nstates,nstates,x_points))

    d_coup=0.0d0
    c1=1.0_dp/12.0_dp
    c2=8.0_dp/12.0_dp

    dx=x_grid(2)-x_grid(1)

    do i=1,nstates
      do j=i+1,nstates
        do x=1,x_points
          eigenv1(:,x)=U(:,i,x)
          eigenv2(:,x)=U(:,j,x)
        end do
        deigenv2_dR=0.0_dp
        do x=3,x_points-2
          deigenv2_dR(:,x)=(c1*(eigenv2(:,x-2)-eigenv2(:,x+2))+ &
            c2*(eigenv2(:,x+1)-eigenv2(:,x-1)))/dx
        end do
        do x=1,x_points
          d_coup(i,j,x)=dot_product(eigenv1(:,x),deigenv2_dR(:,x))
          d_coup(i,j,x)=-d_coup(j,i,x)
        end do
      end do
    end do

    do i=1,nstates
      do j=i+1,nstates
        write(idx1,"(i1.1)") i
        write(idx2,"(i1.1)") j
        filename="nac1-"//trim(idx1)//trim(idx2)//"_x"
        open(134,file=trim(filename),status='replace',form='formatted', &
          action='write',iostat=ioerr)
        if(ioerr/=0) print*, 'error writing ',trim(filename)
        do x=1,x_points
          write(134,'(5f14.6)') d_coup(i,j,x),x_grid(x)
        end do
        close(134)
      end do
    end do

    deallocate(eigenv1,eigenv2,deigenv2_dR)

    deallocate(d_coup)

  end subroutine non_adiabatic_couplings_1d


  subroutine non_adiabatic_couplings_2d(U)

    real(kind=8),intent(in) :: U(nstates,nstates,x_points,y_points)
    real(kind=8),allocatable :: eigenv1(:,:,:),eigenv2(:,:,:), &
      deigenv2_dR(:,:,:),d_coup(:,:,:,:)
    real(kind=8) :: dx,dy,c1,c2
    integer :: x,y,i,j,ioerr
    character(len=500) :: idx1,idx2,filename

    allocate(eigenv1(nstates,x_points,y_points), &
      eigenv2(nstates,x_points,y_points), &
      deigenv2_dR(nstates,x_points,y_points),      &
      d_coup(nstates,nstates,x_points,y_points))

    c1=1.0_dp/12.0_dp
    c2=8.0_dp/12.0_dp

    d_coup=0.0d0
    dx=x_grid(2)-x_grid(1)

    do i=1,nstates
      do j=i+1,nstates
        do x=1,x_points
          do y=1,y_points
            eigenv1(:,x,y)=U(:,i,x,y)
            eigenv2(:,x,y)=U(:,j,x,y)
          end do
        end do
        deigenv2_dR=0.0_dp
        do x=3,x_points-2
          do y=1,y_points
            deigenv2_dR(:,x,y)=(c1*(eigenv2(:,x-2,y)-eigenv2(:,x+2,y))+ &
              c2*(eigenv2(:,x+1,y)-eigenv2(:,x-1,y)))/dx
          end do
        end do
        do x=1,x_points
          do y=1,y_points
            d_coup(i,j,x,y)=dot_product(eigenv1(:,x,y),deigenv2_dR(:,x,y))
            d_coup(i,j,x,y)=-d_coup(j,i,x,y)
          end do
        end do
      end do
    end do

    do i=1,nstates
      do j=i+1,nstates
        write(idx1,"(i1.1)") i
        write(idx2,"(i1.1)") j
        filename="nac1-"//trim(idx1)//trim(idx2)//"_x"
        open(134,file=trim(filename),status='replace',form='formatted', &
          action='write',iostat=ioerr)
        if(ioerr/=0) print*, 'error writing ',trim(filename)
        do x=1,x_points
          do y=1,y_points
            write(134,'(5f14.6)') d_coup(i,j,x,y),x_grid(x),y_grid(y)
          end do
          write(134,*) ""
        end do
        close(134)
      end do
    end do

    d_coup=0.0d0
    dy=y_grid(2)-y_grid(1)

    do i=1,nstates
      do j=i+1,nstates
        do x=1,x_points
          do y=1,y_points
            eigenv1(:,x,y)=U(:,i,x,y)
            eigenv2(:,x,y)=U(:,j,x,y)
          end do
        end do
        deigenv2_dR=0.0_dp
        do x=1,x_points
          do y=3,y_points-2
            deigenv2_dR(:,x,y)=(c1*(eigenv2(:,x,y-2)-eigenv2(:,x,y+2))+ &
              c2*(eigenv2(:,x,y+1)-eigenv2(:,x,y-1)))/dy
          end do
        end do
        do x=1,x_points
          do y=1,y_points
            d_coup(i,j,x,y)=dot_product(eigenv1(:,x,y),deigenv2_dR(:,x,y))
            d_coup(i,j,x,y)=-d_coup(j,i,x,y)
          end do
        end do
      end do
    end do

    do i=1,nstates
      do j=i+1,nstates
        write(idx1,"(i1.1)") i
        write(idx2,"(i1.1)") j
        filename="nac1-"//trim(idx1)//trim(idx2)//"_y"
        open(134,file=trim(filename),status='replace',form='formatted', &
          action='write',iostat=ioerr)
        if(ioerr/=0) print*, 'error writing ',trim(filename)
        do x=1,x_points
          do y=1,y_points
            write(134,'(5f14.6)') d_coup(i,j,x,y),x_grid(x),y_grid(y)
          end do
          write(134,*) ""
        end do
        close(134)
      end do
    end do

    deallocate(eigenv1,eigenv2,deigenv2_dR)

    deallocate(d_coup)

  end subroutine non_adiabatic_couplings_2d


  subroutine non_adiabatic_couplings_3d(U)

    real(kind=8),intent(in) :: U(nstates,nstates,x_points,y_points,x_points)
    real(kind=8),allocatable :: eigenv1(:,:,:,:),eigenv2(:,:,:,:), &
      deigenv2_dR(:,:,:,:),d_coup(:,:,:,:,:)
    real(kind=8) :: dx,dy,dz,c1,c2
    integer :: x,y,z,i,j,ioerr
    character(len=500) :: idx1,idx2,filename

    allocate(eigenv1(nstates,x_points,y_points,z_points), &
      eigenv2(nstates,x_points,y_points,z_points), &
      deigenv2_dR(nstates,x_points,y_points,z_points),      &
      d_coup(nstates,nstates,x_points,y_points,z_points))

    c1=1.0_dp/12.0_dp
    c2=8.0_dp/12.0_dp

    d_coup=0.0d0
    dx=x_grid(2)-x_grid(1)

    do i=1,nstates
      do j=i+1,nstates
        do x=1,x_points
          do y=1,y_points
            do z=1,z_points
              eigenv1(:,x,y,z)=U(:,i,x,y,z)
              eigenv2(:,x,y,z)=U(:,j,x,y,z)
            end do
          end do
        end do
        deigenv2_dR=0.0_dp
        do x=3,x_points-2
          do y=1,y_points
            do z=1,z_points
              deigenv2_dR(:,x,y,z)=(c1*(eigenv2(:,x-2,y,z)-eigenv2(:,x+2,y,z))+ &
                c2*(eigenv2(:,x+1,y,z)-eigenv2(:,x-1,y,z)))/dx
            end do
          end do
        end do
        do x=1,x_points
          do y=1,y_points
            do z=1,z_points
              d_coup(i,j,x,y,z)=dot_product(eigenv1(:,x,y,z),deigenv2_dR(:,x,y,z))
              d_coup(i,j,x,y,z)=-d_coup(j,i,x,y,z)
            end do
          end do
        end do
      end do
    end do

    do i=1,nstates
      do j=i+1,nstates
        write(idx1,"(i1.1)") i
        write(idx2,"(i1.1)") j
        filename="nac1-"//trim(idx1)//trim(idx2)//"_x"
        open(134,file=trim(filename),status='replace',form='formatted', &
          action='write',iostat=ioerr)
        if(ioerr/=0) print*, 'error writing ',trim(filename)
        do x=1,x_points
          do y=1,y_points
            do z=1,z_points
              write(134,'(5f14.6)') d_coup(i,j,x,y,z),x_grid(x),y_grid(y),z_grid(x)
            end do
            write(134,*) ""
          end do
          write(134,*) ""
        end do
        close(134)
      end do
    end do

    d_coup=0.0d0
    dy=y_grid(2)-y_grid(1)

    do i=1,nstates
      do j=i+1,nstates
        do x=1,x_points
          do y=1,y_points
            do z=1,z_points
              eigenv1(:,x,y,z)=U(:,i,x,y,z)
              eigenv2(:,x,y,z)=U(:,j,x,y,z)
            end do
          end do
        end do
        deigenv2_dR=0.0_dp
        do x=1,x_points
          do y=3,y_points-2
            do z=1,z_points
              deigenv2_dR(:,x,y,z)=(c1*(eigenv2(:,x,y-2,z)-eigenv2(:,x,y+2,z))+ &
                c2*(eigenv2(:,x,y+1,z)-eigenv2(:,x,y-1,z)))/dy
            end do
          end do
        end do
        do x=1,x_points
          do y=1,y_points
            do z=1,z_points
              d_coup(i,j,x,y,z)=dot_product(eigenv1(:,x,y,z),deigenv2_dR(:,x,y,z))
              d_coup(i,j,x,y,z)=-d_coup(j,i,x,y,z)
            end do
          end do
        end do
      end do
    end do

    do i=1,nstates
      do j=i+1,nstates
        write(idx1,"(i1.1)") i
        write(idx2,"(i1.1)") j
        filename="nac1-"//trim(idx1)//trim(idx2)//"_y"
        open(134,file=trim(filename),status='replace',form='formatted', &
          action='write',iostat=ioerr)
        if(ioerr/=0) print*, 'error writing ',trim(filename)
        do x=1,x_points
          do y=1,y_points
            do z=1,z_points
              write(134,'(5f14.6)') d_coup(i,j,x,y,z),x_grid(x),y_grid(y),z_grid(z)
            end do
            write(134,*) ""
          end do
          write(134,*) ""
        end do
        close(134)
      end do
    end do

    d_coup=0.0d0
    dz=z_grid(2)-z_grid(1)

    do i=1,nstates
      do j=i+1,nstates
        do x=1,x_points
          do y=1,y_points
            do z=1,z_points
              eigenv1(:,x,y,z)=U(:,i,x,y,z)
              eigenv2(:,x,y,z)=U(:,j,x,y,z)
            end do
          end do
        end do
        deigenv2_dR=0.0_dp
        do x=1,x_points
          do y=1,y_points
            do z=3,z_points-2
              deigenv2_dR(:,x,y,z)=(c1*(eigenv2(:,x,y,z-2)-eigenv2(:,x,y,z+2))+ &
                c2*(eigenv2(:,x,y,z+1)-eigenv2(:,x,y,z-1)))/dz
            end do
          end do
        end do
        do x=1,x_points
          do y=1,y_points
            do z=1,z_points
              d_coup(i,j,x,y,z)=dot_product(eigenv1(:,x,y,z),deigenv2_dR(:,x,y,z))
              d_coup(i,j,x,y,z)=-d_coup(j,i,x,y,z)
            end do
          end do
        end do
      end do
    end do

    do i=1,nstates
      do j=i+1,nstates
        write(idx1,"(i1.1)") i
        write(idx2,"(i1.1)") j
        filename="nac1-"//trim(idx1)//trim(idx2)//"_z"
        open(134,file=trim(filename),status='replace',form='formatted', &
          action='write',iostat=ioerr)
        if(ioerr/=0) print*, 'error writing ',trim(filename)
        do x=1,x_points
          do y=1,y_points
            do z=1,z_points
              write(134,'(5f14.6)') d_coup(i,j,x,y,z),x_grid(x),y_grid(y),z_grid(z)
            end do
            write(134,*) ""
          end do
          write(134,*) ""
        end do
        close(134)
      end do
    end do

    deallocate(eigenv1,eigenv2,deigenv2_dR)

    deallocate(d_coup)

  end subroutine non_adiabatic_couplings_3d


end module linear_algebra
