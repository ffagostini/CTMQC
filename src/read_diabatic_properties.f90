module read_diabatic_properties
  use variables
  use kinds
  implicit none

  character(len=500) :: path

  contains

  subroutine read_diabatic_hamiltonian

   path=trim(path_to_potentials)

   call read_gridpoints
   call read_hamiltonian

  end subroutine read_diabatic_hamiltonian


  subroutine read_gridpoints

    character(len=200) :: filename,basename
    real(kind=dp),allocatable :: tmp_energy(:)
    integer :: ioerr,x,y,z

    allocate(tmp_energy(nstates*nstates))
    
    basename=trim('electronic_hamiltonian.dat')
    filename=trim(path)//trim(basename)

    open(301,file=filename,status='old',form='formatted',&
      action='read',iostat=ioerr)
    if(ioerr/=0) print*,'problem reading ',filename

    xloop: do x=1,x_points
      if(n_dof==1) then
        read(301,*) tmp_energy,x_grid(x)
      end if
      yloop: do y=1,y_points
        if(n_dof==2) then
          read(301,*) tmp_energy,x_grid(x),y_grid(y)
        end if
        zloop: do z=1,z_points
          if(n_dof==3) then
            read(301,*) tmp_energy,x_grid(x),y_grid(y),z_grid(z)
          end if
        end do zloop
        if(n_dof==3) read(301,*)
      end do yloop
      if(n_dof==2) read(301,*)
    end do xloop

    close(301)

    deallocate(tmp_energy)

  end subroutine read_gridpoints


  subroutine read_hamiltonian

    character(len=200) :: filename,basename
    integer :: i,j,k,ioerr,x,y,z
    real(kind=dp) :: tmp(n_dof)
    real(kind=dp),allocatable :: tmp_energy(:)

    allocate(tmp_energy(nstates*nstates))

    basename=trim('electronic_hamiltonian.dat')
    filename=trim(path)//trim(basename)

    open(330,file=filename,status='unknown',form='formatted',&
      action='read',iostat=ioerr)
    if(ioerr/=0) print*,'error reading ',filename

    xloop: do x=1,x_points
      yloop: do y=1,y_points
        zloop: do z=1,z_points
          read(330,*) tmp_energy,tmp
          k=1
          do i=1,nstates
            do j=1,nstates
              Hel(i,j,x,y,z)=tmp_energy(k)
              k=k+1
            end do
          end do
        end do zloop
        if(n_dof==3) read(330,*)
      end do yloop
      if(n_dof==2) read(330,*)
    end do xloop

    close(330)

  end subroutine read_hamiltonian


  subroutine read_transformation_matrix_1d()

    integer :: x,i,j,k,ioerr
    character(len=100) :: basename,filename
    real(kind=dp),allocatable :: tmp_vector(:)
    real(kind=dp) :: tmp_scalar

    allocate(tmp_vector(nstates*nstates))

    path=trim(path_to_potentials)
    basename="transformation_matrix"
    filename=trim(path)//trim(basename)

    open(401,file=filename,status='unknown',form='formatted',&
      action='read',iostat=ioerr)
    if(ioerr/=0) print*,'error reading ',filename

    xloop: do x=1,x_points
      read(401,*) tmp_vector,tmp_scalar
      k=0
      do i=1,nstates
        do j=1,nstates
          k=k+1
          transformation_matrix(i,j,x,1,1)=tmp_vector(k) ! from adiabatic to diabatic basis
        end do
      end do
    end do xloop

    close(401)

    deallocate(tmp_vector)

  end subroutine read_transformation_matrix_1d


end module read_diabatic_properties
