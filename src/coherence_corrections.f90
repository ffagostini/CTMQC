module coherence_corrections
  use variables
  use kinds
  implicit none

  contains

  subroutine accumulated_BOforce(coeff,force,langevin_force,trajlabel)

    integer,intent(in) :: trajlabel
    complex(kind=dp),intent(in) :: coeff(nstates)
    real(kind=dp),intent(in) :: langevin_force(n_dof)
    real(kind=dp),intent(out) :: force(n_dof,nstates)
    integer :: i,j,i_dof,check
    real(kind=dp),parameter :: threshold=0.005_dp
    complex(kind=dp),allocatable :: rho(:,:)

    allocate(rho(nstates,nstates),stat=check)
    if(check/=0) print*,'error 1 rho in accumulated_BOforce'
    do i=1,nstates
      do j=1,nstates
        rho(i,j)=conjg(coeff(i))*coeff(j)
      end do
    end do

    do i=1,nstates
      if(abs(rho(i,i))>threshold .and. abs(rho(i,i))<1.0_dp-threshold) then
        do i_dof=1,n_dof
          force(i_dof,i)=force(i_dof,i)+ &
            dt*(nabla_dot_phase(rho,i,i_dof,trajlabel)+ &
            langevin_force(i_dof))
        end do
      else
        force(:,i)=0.0_dp
      end if
    end do

    deallocate(rho,stat=check)
    if(check/=0) print*,'error 2 rho in accumulated_BOforce'

  end subroutine accumulated_BOforce


  function nabla_dot_phase(rho,state,i_dof,trajlabel)

    integer,intent(in) :: state,i_dof,trajlabel
    complex(kind=8),intent(in) :: rho(nstates,nstates)
    real(kind=dp) :: nabla_dot_phase,mean_force
    integer :: i

    mean_force=0.0_dp
    do i=1,nstates
      mean_force=mean_force+real(rho(i,i),kind=dp)* &
        BOforce(trajlabel,i_dof,i)
     end do

     nabla_dot_phase=BOforce(trajlabel,i_dof,state)!-mean_force

  end function nabla_dot_phase


  subroutine quantum_momentum(Rcl,acc_force,BOsigma,k_li,time)
    use output

    integer,intent(in) :: time
    real(kind=dp),intent(in) :: Rcl(ntraj,n_dof), &
      acc_force(ntraj,n_dof,nstates)
    complex(kind=dp),intent(in) :: BOsigma(ntraj,nstates,nstates)
    real(kind=dp),intent(inout) :: k_li(ntraj,nstates,nstates)
    integer :: itraj,jtraj,i_dof,index_ij,istate,jstate
    real(kind=dp),allocatable :: avR(:),avR2(:),gamma(:,:), &
      g_i(:),prod_g_i(:,:),w_ij(:,:,:),slope_i(:,:),ratio(:,:,:), &
      num_old(:,:,:),num_new(:,:,:),num(:,:,:),denom(:,:),qmom(:,:,:)
    integer,allocatable :: ntr_j(:)
    real(kind=dp) :: dist2,dist_cutoff,threshold,a

    allocate(avR(n_dof),avR2(n_dof),gamma(n_dof,ntraj),ntr_j(n_dof), &
      g_i(ntraj),prod_g_i(ntraj,ntraj),w_ij(n_dof,ntraj,ntraj), &
      slope_i(n_dof,ntraj),ratio(n_dof,npairs,ntraj),num_old(n_dof,ntraj,npairs), &
      num_new(n_dof,ntraj,npairs),num(n_dof,ntraj,npairs),denom(n_dof,npairs), &
      qmom(n_dof,ntraj,npairs))
    gamma=0.0_dp

    if(time==1) then
      dist_cutoff=M_parameter*minval(sigma)/dble(ntraj)
      store_gamma=sqrt(minval(sigma)/dble(ntraj))!dist_cutoff
    end if

    threshold=M_parameter*minval(sigma)/dble(ntraj)

    itrajloop: do itraj=1,ntraj
      avR=0._dp
      avR2=0._dp
      ntr_j=0
      !In this loop I check the distance among the trajectories from the
      !selected trajectry (itraj) and I estimate the variance associated to the
      !gaussian on the trajectory (itraj)
      !Noticde: I have 3N_n variances for each trajectories
      dist2=0._dp
      i_dofloop: do i_dof=1,n_dof
        if(time/=1) dist_cutoff=M_parameter*store_gamma(i_dof,itraj)
        jtrajloop: do jtraj=1,ntraj
          dist2=(Rcl(itraj,i_dof)-Rcl(jtraj,i_dof))**2
          !if(sqrt(dist2) .le. min(M_parameter*store_gamma(i_dof,itraj),dist_cutoff)) then
          if(sqrt(dist2) .le. min(36._dp*store_gamma(i_dof,itraj),dist_cutoff)) then
          !if(sqrt(dist2) .le. min(store_gamma(i_dof,itraj),dist_cutoff)) then
            avR(i_dof)=avR(i_dof)+Rcl(jtraj,i_dof)
            !Average position of trajectories within a sphere
            avR2(i_dof)=avR2(i_dof)+Rcl(jtraj,i_dof)**2
            !Average squared position of trajectories within a sphere
            ntr_j(i_dof)=ntr_j(i_dof)+1
            !Number of trajectories within the sphere
          end if
        end do jtrajloop
        avR(i_dof)=avR(i_dof)/dble(ntr_j(i_dof))
        avR2(i_dof)=avR2(i_dof)/dble(ntr_j(i_dof))
        gamma(i_dof,itraj)=sqrt((avR2(i_dof)-avR(i_dof)**2))/dble(ntr_j(i_dof))
        !write(*,*) gamma(i_dof,itraj)
        if((gamma(i_dof,itraj))**2 .lt. threshold .or. ntr_j(i_dof)==1) then
          !if(gamma(i_dof,itraj)<0.00000001_dp) print*,'ciao'
          gamma(i_dof,itraj)=sqrt(threshold)!dist_cutoff!
          !if(gamma(i_dof,itraj)<0.00000001_dp) then
          !  print*,'ciao'
          !else
          !  write(6,*) gamma(i_dof,itraj),ntr_j(i_dof)
          !end if
          !Variance associated to each trajectory depending on the spreading
          !of the trajectories close to the trajectory (itraj)
          !if the variance is too small, I set it to a given value dist_cutoff
        end if
      end do i_dofloop
    end do itrajloop

    !gamma=1._dp
    store_gamma=gamma

    do itraj=1,ntraj
      g_i(itraj)=0.0_dp
      !This part computes the nuclear density as the sum of normalized
      !Gaussians centered at the positions of the trajectories (labeled jtraj)
      !with variances as computed above
      !!!!! I need to know the value of the nuclear density at the position
      !!!!! of the trajectory (itraj)
      do jtraj=1,ntraj
        prod_g_i(itraj,jtraj)=1.0_dp
        do i_dof=1,n_dof
          prod_g_i(itraj,jtraj)=prod_g_i(itraj,jtraj)*              &
           (dexp(-(Rcl(itraj,i_dof)-Rcl(jtraj,i_dof))**2/           &
           (gamma(i_dof,jtraj))**2/2._dp))*                              &
           (1.0_dp/sqrt(2.0_dp*PI*(gamma(i_dof,jtraj))**2))
        end do
        g_i(itraj)=g_i(itraj)+prod_g_i(itraj,jtraj)
      end do
    end do

    !W_ij => see SI of paper Min, Agostini, Tavernelli, Gross for its definition
    !This part computes W_ij, whose sum over j (jtraj in the loop below) is the
    !slope of the quantum momentum when a sum of Gaussians is used
    !to approximate the nuclear density
    do itraj=1,ntraj
      do jtraj=1,ntraj
        do i_dof=1,n_dof
          w_ij(i_dof,itraj,jtraj)=prod_g_i(itraj,jtraj)/ &
            2.0_dp/(gamma(i_dof,jtraj))**2/g_i(itraj)
        end do
      end do
    end do
    !The slope is calculated here as a sum over j of W_ij
    slope_i=0.0_dp
    do itraj=1,ntraj
      do i_dof=1,n_dof
        do jtraj=1,ntraj
          slope_i(i_dof,itraj)=slope_i(i_dof,itraj)- &
            w_ij(i_dof,itraj,jtraj)
        end do
      end do
    end do

    !Here I compute the center of the quantum momentum
    !See SI of paper Min, Agostini, Tavernelli, Gross for the expression Eq.(28)
    do i_dof=1,n_dof
      index_ij=0
      do istate=1,nstates
        do jstate=istate+1,nstates
          index_ij=index_ij+1
          ! denominator
          denom(i_dof,index_ij)=0._dp
          do jtraj=1,ntraj
            denom(i_dof,index_ij)=denom(i_dof,index_ij)+  &
              real(BOsigma(jtraj,istate,istate),kind=dp)* &
              real(BOsigma(jtraj,jstate,jstate),kind=dp)* &
              (acc_force(jtraj,i_dof,istate)-                   &
              acc_force(jtraj,i_dof,jstate))*slope_i(i_dof,jtraj)
          end do
        end do
      end do
    end do

    ratio=0.0_dp
    do itraj=1,ntraj
      do i_dof=1,n_dof
        index_ij=0
        do istate=1,nstates
          do jstate=istate+1,nstates
            index_ij=index_ij+1
            ! numerator
            num_old(i_dof,itraj,index_ij)=  &
              real(BOsigma(itraj,istate,istate),kind=dp)* &
              real(BOsigma(itraj,jstate,jstate),kind=dp)*     &
              (acc_force(itraj,i_dof,istate)-acc_force(itraj,i_dof,jstate))* &
              Rcl(itraj,i_dof)*slope_i(i_dof,itraj)
            ! ratio
            if(abs(denom(i_dof,index_ij)) .lt. 0.00000001_dp) then
              ratio(i_dof,index_ij,itraj)=0._dp
            else
              ratio(i_dof,index_ij,itraj)=num_old(i_dof,itraj,index_ij)/&
                denom(i_dof,index_ij)
            end if
          end do
        end do
      end do
    end do

    !Here I acutally compute the sum over I of Eq.(28) of SI of
    !paper Min, Agostini, Tavernelli, Gross
    do itraj=1,ntraj
      do i_dof=1,n_dof
        index_ij=0
        do istate=1,nstates
          do jstate=istate+1,nstates
            index_ij=index_ij+1
            num_old(i_dof,itraj,index_ij)=0._dp
            do jtraj=1,ntraj
              num_old(i_dof,itraj,index_ij)=num_old(i_dof,itraj,index_ij)+   &
                ratio(i_dof,index_ij,jtraj)
            end do
            if(abs(slope_i(i_dof,itraj)) .lt. 0.0000001_dp .or. &
              num_old(i_dof,itraj,index_ij) .eq. 0.0_dp) then
                num_old(i_dof,itraj,index_ij)=Rcl(itraj,i_dof)
            end if
          end do
        end do
      end do
    end do

    !Here the centers of the quantum momentum are cmoputed from the expression (21)
    !of SI of paper Min, Agostini, Tavernelli, Gross. This is just to be sure that
    !if the previous
    do itraj=1,ntraj
      do i_dof=1,n_dof
        index_ij=0
        do istate=1,nstates
          do jstate=istate+1,nstates
            index_ij=index_ij+1
            num_new(i_dof,itraj,index_ij)=0.0_dp
            if(abs(slope_i(i_dof,itraj))<0.00000001_dp) then
              num_new(i_dof,itraj,index_ij)=Rcl(itraj,i_dof)
            else
              do jtraj=1,ntraj
                num_new(i_dof,itraj,index_ij)=num_new(i_dof,itraj,index_ij)+  &
                  Rcl(jtraj,i_dof)*prod_g_i(itraj,jtraj)/2.0_dp/(gamma(i_dof,jtraj))**2/ &
                  g_i(itraj)/(-slope_i(i_dof,itraj))
              end do
            end if
          end do
        end do
      end do
    end do

    do itraj=1,ntraj
      do i_dof=1,n_dof
        index_ij=0
        do istate=1,nstates
          do jstate=istate+1,nstates
            index_ij=index_ij+1
            if(abs(num_old(i_dof,itraj,index_ij)-Rcl(itraj,i_dof))>M_parameter*sigma(i_dof)) then
              if(abs(num_new(i_dof,itraj,index_ij)-Rcl(itraj,i_dof))>M_parameter*sigma(i_dof)) then
            !if(abs(num_old(i_dof,itraj,index_ij)-Rcl(itraj,i_dof))>sigma(i_dof)) then
            !  if(abs(num_new(i_dof,itraj,index_ij)-Rcl(itraj,i_dof))>sigma(i_dof)) then
            !if(abs(num_old(i_dof,itraj,index_ij)-Rcl(itraj,i_dof))>1.0_dp) then
            !  if(abs(num_new(i_dof,itraj,index_ij)-Rcl(itraj,i_dof))>1.0_dp) then
                num(i_dof,itraj,index_ij)=Rcl(itraj,i_dof)
                !write(*,*) "Num 1"
              else
                num(i_dof,itraj,index_ij)=num_new(i_dof,itraj,index_ij)
                !write(*,*) "Num 2"
              end if
            else
              num(i_dof,itraj,index_ij)=num_old(i_dof,itraj,index_ij)
              !write(*,*) "Num 3"
            end if
          end do
        end do
      end do
    end do

    ! quantum momentum
    do itraj=1,ntraj
      do i_dof=1,n_dof
        index_ij=0
        do istate=1,nstates
          do jstate=istate+1,nstates
            index_ij=index_ij+1
            qmom(i_dof,itraj,index_ij)=slope_i(i_dof,itraj)*(Rcl(itraj,i_dof)- &
              num(i_dof,itraj,index_ij))
          end do
        end do
      end do
    end do

    index_ij=0
    do istate=1,nstates
      do jstate=istate+1,nstates
        index_ij=index_ij+1
        do itraj=1,ntraj
          k_li(itraj,istate,jstate)=0.0_dp
          k_li(itraj,jstate,istate)=0.0_dp
          do i_dof=1,n_dof
            !write(6,*) i_dof,itraj,index_ij,istate,jstate
            k_li(itraj,istate,jstate)=k_li(itraj,istate,jstate)+(2.0_dp/mass(i_dof))*  &
              qmom(i_dof,itraj,index_ij)*acc_force(itraj,i_dof,istate)
            k_li(itraj,jstate,istate)=k_li(itraj,jstate,istate)+(2.0_dp/mass(i_dof))*  &
              qmom(i_dof,itraj,index_ij)*acc_force(itraj,i_dof,jstate)
          end do
        end do
      end do
    end do

    deallocate(avR,avR2,gamma,ntr_j,g_i,prod_g_i,w_ij,slope_i, &
      ratio,num_old,num_new,num,denom,qmom)

  end subroutine quantum_momentum


end module coherence_corrections

