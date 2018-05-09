module time_evolution
  use omp_lib
  use variables
  use kinds
  use electronic_problem
  use coefficients_evolution
  use classical_evolution
  use coherence_corrections
  use output
  implicit none

  real(kind=dp),allocatable :: Rcl(:,:),Vcl(:,:),classical_force(:)
  real(kind=dp),allocatable :: my_force(:,:,:),k_li(:,:,:)
  complex(kind=dp),allocatable :: BOsigma(:,:,:)
  complex(kind=dp),allocatable :: BOcoeff(:,:)

  contains

  subroutine evolution

    integer :: time,itraj,i,j

    call input_summary

    call initialize_local_vars
    call plot(BOsigma,Rcl,Vcl,time=0)

    timeloop: do time=1,nsteps

      if(mod(time-1,dump)==0) write(*,'(a,1x,f14.2)') 'time=',dble(time-1)*dt

      !!!!$omp parallel do private(itraj,classical_force) &
      !!!!$omp shared(Rcl,Vcl,ntraj) &
      !!!!$omp shared(BOcoeff,my_force,k_li) &
      !!!!$omp default(none)
      trajsloop: do itraj=1,ntraj
        call BOproblem(Rcl(itraj,:),itraj)
        call accumulated_BOforce(BOcoeff(itraj,:),my_force(itraj,:,:),itraj)
        call non_adiabatic_force(BOcoeff(itraj,:),classical_force, &
          my_force(itraj,:,:),k_li(itraj,:,:),itraj)
        call RK4_coeff(Vcl(itraj,:),BOcoeff(itraj,:),k_li(itraj,:,:),itraj)
        call velocity_verlet(Rcl(itraj,:),Vcl(itraj,:),classical_force)
      end do trajsloop
      !!!!$omp end parallel do

      do i=1,nstates
        do j=1,nstates
          BOsigma(:,i,j)=conjg(BOcoeff(:,i))*BOcoeff(:,j)
        end do
      end do

      if(mod(time,dump)==0) call plot(BOsigma,Rcl,Vcl,time)

      call quantum_momentum(Rcl,my_force,BOsigma,k_li,time)

    end do timeloop

    call finalize_local_vars

  end subroutine evolution
  

  subroutine input_summary

    write(6,"(a,1x,a)") "Model system ",trim(model_potential)
    write(6,"(a,i5)") "Initial BO state",initial_BOstate
    write(6,"(a,i5,1x,a)") "Running",ntraj,"trajectories"
    write(6,"(a,f14.2,1x,a)") "centered at",r0,"a.u."
    write(6,"(a,f14.2,1x,a)") "with initial momentum",k0,"a.u."
    write(6,"(a,f14.2,1x,a)") "for",final_time,"a.u."
    write(6,"(a,f14.4,1x,a)") "with time-step",dt,"a.u."
    write(6,"(a,i5,1x,a)") "dumping every",dump,"steps"

    write(6,"(a)") "Starting dynamics..."

  end subroutine input_summary


  subroutine initialize_local_vars

    integer :: itraj,i,j

    allocate(Rcl(ntraj,n_dof),Vcl(ntraj,n_dof), &
      classical_force(n_dof),BOsigma(ntraj,nstates,nstates), &
      BOcoeff(ntraj,nstates),my_force(ntraj,n_dof,nstates),  &
      k_li(ntraj,nstates,nstates))

    my_force=0.0_dp

    BOcoeff=cmplx(0.0_dp,0.0_dp)
    BOcoeff(:,initial_BOstate)=cmplx(1.0_dp,0.0_dp)

    do itraj=1,ntraj
      do i=1,nstates
        do j=1,nstates
          BOsigma(itraj,i,j)=conjg(BOcoeff(itraj,i))*BOcoeff(itraj,j)
        end do
      end do
      Rcl(itraj,:)=initial_positions(itraj,:)
      Vcl(itraj,:)=initial_momenta(itraj,:)/mass
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

    deallocate(my_force,stat=check)
    if(check/=0) print*,'error my_force'
    deallocate(k_li,stat=check)
    if(check/=0) print*,'error k_li'

  end subroutine finalize_local_vars

end module time_evolution
