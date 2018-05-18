module time_evolution
  use omp_lib
  use variables
  use kinds
  use electronic_problem
  use coefficients_evolution
  use classical_evolution
  use coherence_corrections
  use output
  use tools
  implicit none

  real(kind=dp),allocatable :: Rcl(:,:),Vcl(:,:),classical_force(:)
  real(kind=dp),allocatable :: my_force(:,:,:),k_li(:,:,:)
  complex(kind=dp),allocatable :: BOsigma(:,:,:)
  complex(kind=dp),allocatable :: BOcoeff(:,:),DIAcoeff(:,:)

  contains

  subroutine evolution

    integer(kind=4) :: time
    integer :: itraj,i,j

    call input_summary

    call initialize_local_vars
    if(dia_to_ad=="y") call diabatic_output(Rcl,BOcoeff,time=0)
    call plot(BOsigma,Rcl,Vcl,time=0)

    if(model_system=="marcus") call integrator_parameters

    timeloop: do time=1,nsteps

      if(mod(time,dump)==0) write(*,'(a,1x,f14.2)') 'time=',dble(time)*dt

      !!!!$omp parallel do private(itraj,classical_force) &
      !!!!$omp shared(Rcl,Vcl,ntraj) &
      !!!!$omp shared(BOcoeff,my_force,k_li) &
      !!!!$omp default(none)
      trajsloop: do itraj=1,ntraj
        call BOproblem(Rcl(itraj,:),itraj)
        if(algorithm=="CTQMC" .or. algorithm=="CTeMQC") &
          call accumulated_BOforce(BOcoeff(itraj,:),my_force(itraj,:,:),itraj)
        call non_adiabatic_force(BOcoeff(itraj,:),classical_force, &
          my_force(itraj,:,:),k_li(itraj,:,:),itraj,Vcl(itraj,:))
        call RK4_coeff(Vcl(itraj,:),BOcoeff(itraj,:),k_li(itraj,:,:),itraj)
        call velocity_verlet(Rcl(itraj,:),Vcl(itraj,:),classical_force)
      end do trajsloop
      !!!!$omp end parallel do

      do i=1,nstates
        do j=1,nstates
          BOsigma(:,i,j)=conjg(BOcoeff(:,i))*BOcoeff(:,j)
        end do
      end do

      if(mod(time,dump)==0) then
        if(dia_to_ad=="y") call diabatic_output(Rcl,BOcoeff,time)
        call plot(BOsigma,Rcl,Vcl,time)
      end if

      if(trim(algorithm)=="CTMQC" .or. trim(algorithm)=="CTeMQC") &
        call quantum_momentum(Rcl,my_force,BOsigma,k_li,time)

    end do timeloop

    call finalize_local_vars

  end subroutine evolution
  

  subroutine input_summary

    write(6,"(a,1x,a)") "Model system:",trim(model_system)
    if(initial_BOstate/=0) write(6,"(a,i5)") "Initial BO state",initial_BOstate
    if(initial_DIAstate/=0) write(6,"(a,i5)") "Initial DIA state",initial_DIAstate
    write(6,"(a,i5,1x,a)") "Running",ntraj,"trajectories"
    write(6,"(a,f14.4,1x,a)") "centered at",r0,"a.u."
    write(6,"(a,f14.4,1x,a)") "with initial momentum",k0,"a.u."
    write(6,"(a,f14.4,1x,a)") "variance",sqrt(sigma),"a.u."
    write(6,"(a,es11.3,1x,a)") "for",final_time,"a.u."
    write(6,"(a,f14.4,1x,a)") "with time-step",dt,"a.u."
    write(6,"(a,1x,i8,1x,a)") "dumping every",dump,"steps"
    write(6,"(a,1x,i7)") "Total number of printed snapshots", &
      int((final_time/dt)/dble(dump))

    write(6,"(a,1x,a,1x,a)") "Starting",trim(algorithm),"dynamics..."

  end subroutine input_summary


  subroutine initialize_local_vars
    use electronic_problem
    use read_diabatic_properties

    integer :: itraj,i,j,istar(3)

    allocate(Rcl(ntraj,n_dof),Vcl(ntraj,n_dof), &
      classical_force(n_dof),BOsigma(ntraj,nstates,nstates), &
      BOcoeff(ntraj,nstates),my_force(ntraj,n_dof,nstates),  &
      k_li(ntraj,nstates,nstates),DIAcoeff(ntraj,nstates))

    my_force=0.0_dp

    do itraj=1,ntraj
      Rcl(itraj,:)=initial_positions(itraj,:)
      Vcl(itraj,:)=initial_momenta(itraj,:)/mass
    end do

    if(initial_BOstate/=0) then
      BOcoeff=cmplx(0.0_dp,0.0_dp)
      BOcoeff(:,initial_BOstate)=cmplx(1.0_dp,0.0_dp)
    end if
    if(initial_DIAstate/=0) then
      if(dia_to_ad=="y") then
        DIAcoeff=cmplx(0.0_dp,0.0_dp)
        DIAcoeff(:,initial_DIAstate)=cmplx(1.0_dp,0.0_dp)
        call read_transformation_matrix_1d() ! I get the transformation matrix
        !! from the adiabatic to the diabatic basis for different values of x,y,z
        do itraj=1,ntraj
          call locate_x(Rcl(itraj,:),istar)
          BOcoeff(itraj,:)=matmul(transpose(transformation_matrix(:,:, &
            istar(1),istar(2),istar(3))),DIAcoeff(itraj,:))
        end do
        call compute_diabatic_surfaces
      else
        write(6,*) "I need the transformation matrix from adiabatic to diabatic basis"
      end if
    end if

    do itraj=1,ntraj
      do i=1,nstates
        do j=1,nstates
          BOsigma(itraj,i,j)=conjg(BOcoeff(itraj,i))*BOcoeff(itraj,j)
        end do
      end do
      call BOproblem(Rcl(itraj,:),itraj)
    end do

  end subroutine initialize_local_vars


  subroutine finalize_local_vars

    integer :: check

    deallocate(Rcl,stat=check)
    if(check/=0) print*,'error Rcl'
    deallocate(Vcl,stat=check)
    if(check/=0) print*,'error Vcl'

    deallocate(classical_force,stat=check)
    if(check/=0) print*,'error classical_force'

    deallocate(BOsigma,stat=check)
    if(check/=0) print*,'error BOsigma'
    deallocate(BOcoeff,stat=check)
    if(check/=0) print*,'error BOcoeff'
    deallocate(DIAcoeff,stat=check)
    if(check/=0) print*,'error DIAcoeff'

    deallocate(my_force,stat=check)
    if(check/=0) print*,'error my_force'
    deallocate(k_li,stat=check)
    if(check/=0) print*,'error k_li'

  end subroutine finalize_local_vars

end module time_evolution
