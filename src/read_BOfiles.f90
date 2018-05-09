module read_BOfiles
  use variables
  use kinds
  implicit none

  character(len=500) :: path

  contains
  
  subroutine readBOfiles

    path=trim(path_to_potentials)

    call read_gridpoints
    call read_BOsurface
    call read_nacoup
  
  end subroutine readBOfiles


  subroutine read_gridpoints

    character(len=200) :: filename,basename
    real(kind=dp) :: tmp_energy
    integer :: ioerr,x,y,z
    
    basename=trim('1_bopes.dat')
    filename=trim(path)//trim(basename)

    open(301,file=filename,status='old',form='formatted',&
      action='read',iostat=ioerr)
    if(ioerr/=0) print*,'problem reading',filename

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

  end subroutine read_gridpoints


  subroutine read_BOsurface

    character(len=200) :: filename,basename
    character :: pes_idx
    integer :: i,ioerr,unit,x,y,z
    real(kind=dp) :: tmp(n_dof)

    basename='_bopes.dat'
    basename=trim(basename)

    statesloop: do i=1,nstates
      write(pes_idx,'(I1)') i
      filename=trim(path)//pes_idx//trim(basename)
      unit=310+i
      open(unit,file=filename,status='unknown',form='formatted',&
        action='read',iostat=ioerr)
      if(ioerr/=0) print*,'error reading',filename
      xloop: do x=1,x_points
        yloop: do y=1,y_points
          zloop: do z=1,z_points
            read(unit,*) BOpes(x,y,z,i),tmp
          end do zloop
          if(n_dof==3) read(unit,*)
        end do yloop
        if(n_dof==2) read(unit,*)
      end do xloop

      close(unit)

    end do statesloop

  end subroutine read_BOsurface


  subroutine read_nacoup
  
    character(len=200) :: filename,basename
    character :: na1,na2
    integer :: unit,ioerr,i,j,x,y,z,i_dof
    real(kind=dp) :: tmp(n_dof)
    character(len=2) :: idx
    
    basename='nac1-'
    basename=trim(basename)

    do i=1,nstates
      do j=i+1,nstates

        write(na1,'(I1)') i
        write(na2,'(I1)') j

        do i_dof=1,n_dof

          if(i_dof==1) idx="_x"
          if(i_dof==2) idx="_y"
          if(i_dof==3) idx="_z"

          filename=trim(path)//trim(basename)//na1//na2//trim(idx)
          unit=320+i+j
          open(unit,file=filename,status='unknown',form='formatted',&
            action='read',iostat=ioerr)
          if(ioerr/=0) print*,'error reading',filename
          xloop: do x=1,x_points
            yloop: do y=1,y_points
              zloop: do z=1,z_points
                read(unit,*) na_coup(i_dof,x,y,z,i,j),tmp
                na_coup(i_dof,x,y,z,j,i)=-na_coup(i_dof,x,y,z,i,j)
              end do zloop
              if(n_dof==3) read(unit,*)
            end do yloop
            if(n_dof==2) read(unit,*)
          end do xloop
          close(unit)

        end do

      end do
    end do

  end subroutine read_nacoup


end module read_BOfiles
