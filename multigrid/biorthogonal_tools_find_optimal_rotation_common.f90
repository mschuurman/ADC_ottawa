!   subroutine find_optimal_rotation_rk(tmat,umat,vcur)
!     complex(rk), intent(in)  :: tmat(:,:)  ! Subset overlap of the current left MOs with reference right MOs.
!     complex(rk), intent(in)  :: umat(:,:)  ! Subset overlap of the reference left MOs with current right MOs.
!     complex(rk), intent(out) :: vcur(:,:)  ! Unitary transformation of the right MOs.
!                                            ! Left MOs are to be transformed by conjg(vcur)
      !
      !  Keep the code type-generic beyond this point!
      ! 
      integer(ik), parameter        :: max_nr   = 200
      real(kind(tmat)), parameter   :: max_step = 0.50
      integer(ik)                   :: alloc, iter
      integer(ik)                   :: bs         ! Number of eigenvectors in the current degenerate block
      integer(ik)                   :: lg         ! Number of independent rotation parameters
      real(kind(tmat))              :: rcur       ! Length of the NR step
      real(kind(tmat))              :: gcur       ! Current gradient norm
      real(kind(tmat))              :: rgcur      ! Dot product of current gradient and NR step
      complex(kind(tmat))           :: sim        ! Current value of the similarity function
      real(kind(tmat)), allocatable :: grad(:,:)  ! Gradients of |sim|^2. Second dimension has to be 1
      real(kind(tmat)), allocatable :: hess(:,:)  ! Hessian of |sim|^2
      real(kind(tmat)), allocatable :: disp(:,:)  ! Optimization step
      !
      bs = size(tmat,dim=1)
      lg = bs*(bs-1) ! Number elements of rotation matrix above the diagonal, times 2
      allocate (grad(lg,1),hess(lg,lg),disp(lg,1),stat=alloc)
      if (alloc/=0) then
        call stop('biorthogonal_tools%find_optimal_rotation - no memory')
      end if
      !
      !  The initial guess is important: we have to start with something close to the solution,
      !  or iterations won't converge. A permutation matrix, bringing largest elements of the
      !  occupation matrix to the diagonal, seems to be a good choice. 
      !
      call guess_initial_vcur(vcur)
      !
      !  If our eigenvector was non-degenerate, there is nothing to optimize.
      !  If it is degenerate, we'll also have to find the optimal rotation.
      !
      if (lg>0) then
        !
        !  Before starting optimization, perturb the initial guess by a small
        !  random rotation. This should reduce the chance of hitting a non-
        !  optimal stationary point by accident.
        !
        call randomize_initial_guess(vcur)
        !
        !  Solve optimization problem using damped NR iteratoins.
        !
        nr_steps: do iter=1,max_nr
          call gradient_and_hessian(vcur,sim,grad(:,1),hess)
          gcur = sqrt(sum(abs(grad)**2))
          if (verbose>=3) then
            write (out,"('On iteration ',i0,' similarity = ',2f14.6)") iter, sim
            write (out,"('                        gradient norm is ',g14.6)") gcur
          end if
          disp = -grad
          call lapack_gelss(hess,disp)
          rcur = sqrt(sum(abs(disp)**2))
          if (rcur>max_step) then
            disp = (max_step/rcur)*disp
          end if
          rgcur = sum(disp*grad)
          if (verbose>=3) then
            write (out,"('                      NR step has length ',g14.6)") rcur
            write (out,"('                          grad . disp is ',g14.6)") rgcur
          endif
          if (rgcur<0) then 
            !
            !  We are going towards decreasing goal function; reverse
            !
            disp = -disp 
          end if
          call apply_rotation(disp(:,1),vcur)
          !
          !  Stop if further rotations could not possible change the goal function
          !
          if (gcur*rcur<spacing(abs(sim))) exit nr_steps
        end do nr_steps
        if (iter>max_nr) then
          write (out,"(/'WARNING: NR iterations did not converge after ',i0,' steps.')") iter
          write (out,"( 'WARNING:            Current rotation gradient ',g14.6)") gcur
          write (out,"( 'WARNING:                Current rotation step ',g14.6)") rcur
          write (out,"( 'WARNING: Continuing with unconverged transformation; expect trouble later'/)")
        end if
      end if
      !
      !  We are almost done: the per-orbital density matrices are now in the
      !  best possible alignment. However, this still leaves us with two
      !  possible choices of the orbital phases. Fortunately, this one is easy ...
      !
      call choose_orbital_phases(vcur)
      !
      deallocate (grad,hess,disp)
      !
      contains
      !
      subroutine guess_initial_vcur(vc)
        complex(kind(tmat)), intent(out) :: vc(:,:) ! Initial guess for the rotation matrix
        !
        integer(ik)         :: i, j
        complex(kind(tmat)) :: wt(bs,bs)  ! Occupation matrix
        logical             :: taken(bs)
        !
        !  Begin by constructing the initial occupations matrix
        !
        wt = umat * transpose(tmat)
        !
        !  For each orbital in vcur choose the best matching reference orbital
        !  Once 
        !
        vc    = 0
        taken = .false.
        select_curr: do i=1,bs
          j = maxloc(abs(wt(:,i)),dim=1,mask=.not.taken)
          taken(j) = .true.
          vc(i,j)  = 1
        end do select_curr
      end subroutine guess_initial_vcur
      !
      subroutine randomize_initial_guess(vcur)
        complex(kind(tmat)), intent(inout) :: vcur(:,:)
        !
        real(kind(tmat)), parameter :: magnitude = 1e-2
        integer(ik)         :: i, j
        real(kind(tmat))    :: rn_re, rn_im
        complex(kind(tmat)) :: wt(bs,bs)  ! anti-Hermitian rotation exponent
        complex(kind(tmat)) :: vt(bs,bs)  ! finite rotation matrix
        !
        wt = 0
        fill_rotation_i: do i=1,bs
          fill_rotation_j: do j=i+1,bs
            call random_number(rn_re)
            call random_number(rn_im)
            wt(i,j) = cmplx(rn_re,rn_im,kind=kind(tmat))
            wt(j,i) = -conjg(wt(i,j))
          end do fill_rotation_j
        end do fill_rotation_i
        wt = magnitude * wt
        call matrix_exponential(wt,vt)
        vcur = matmul(vcur,vt)
      end subroutine randomize_initial_guess
      !
      subroutine choose_orbital_phases(vcur)
        complex(kind(tmat)), intent(inout) :: vcur(:,:)
        !
        complex(kind(tmat)) :: ut(bs,bs)  ! Transformed ovelap matrix between left reference and right orbitals
        integer(ik)         :: j
        complex(kind(tmat)) :: sg
        !
        ut = matmul(umat,vcur)
        !
        !  It is OK to change the phase of the final vector any which way; choose it to
        !  match the phase of the reference vector.
        !
        scan_signs: do j=1,bs
          sg = abs(ut(j,j))/ut(j,j)
          vcur(:,j) = sg * vcur(:,j)
        end do scan_signs
      end subroutine choose_orbital_phases
      !
      !  Calculate gradient and hessian of the objective function at vcur
      !
      subroutine gradient_and_hessian(vc,sim,grad,hess)
        complex(kind(tmat)), intent(in)  :: vc  (:,:) ! Transformation matrix where function and its derivatives are needed
        complex(kind(tmat)), intent(out) :: sim       ! Similarity parameter
        real(kind(tmat)), intent(out)    :: grad(:)   ! Gradients of |sim|^2 wrt exponential rotation parameters. Elements are:
                                                      ! 1 = Re(w(1,2)); 2 = Im(w(1,2)); 3 = Re(w(1,3)), etc
        real(kind(tmat)), intent(out)    :: hess(:,:) ! Hessian of |sim|^2. Order of the parameters as in grad(:)
        !
        integer(ik)         :: i, j, ij  ! Indices and compound indices
        real(kind(tmat))    :: dw        ! Step size for numerical diffentiation
        real(kind(tmat))    :: tg(lg,4)  ! Temporary gradient vectors
        real(kind(tmat))    :: testg(lg) ! Gradient calculated using numerical differentiation; sanity checking
        complex(kind(tmat)) :: ts(4)     ! Just a dummy ...
        integer(ik)         :: itab(lg)  ! I indices
        integer(ik)         :: jtab(lg)  ! J indices
        !
        ij = 1
        fill_displacements_i: do i=1,bs
          fill_displacements_j: do j=i+1,bs
            itab(ij) = i
            jtab(ij) = j
            ij = ij + 2
          end do fill_displacements_j
        end do fill_displacements_i
        !
        !  Evaluate hessian using symmetric displacements
        !
        dw = 1e-4
        !$omp parallel do default(none) &
        !$omp& shared(bs,lg,vc,dw,itab,jtab,testg,hess,tmat) &
        !$omp& private(tg,ts,ij,i,j)
        scan_displacements: do ij=1,lg,2
          i = itab(ij)
          j = jtab(ij)
          !
          !  Calculate gradients at small real and imaginary displacements along w(i,j)
          !
          call displaced_gradient(vc,i,j,cmplx( dw,  0,kind=kind(tmat)),ts(1),tg(:,1))
          call displaced_gradient(vc,i,j,cmplx(-dw,  0,kind=kind(tmat)),ts(2),tg(:,2))
          call displaced_gradient(vc,i,j,cmplx(  0, dw,kind=kind(tmat)),ts(3),tg(:,3))
          call displaced_gradient(vc,i,j,cmplx(  0,-dw,kind=kind(tmat)),ts(4),tg(:,4))
          !
          !  Evaluate "Hessian" and numerical gradient
          !  Keep in mind that matrix multiplication does not commute, so that
          !  order of displacements matters and our "Hessian" is -not- a symmetric matrix
          !  This is not an issue for Newton-Rafson steps, however.
          !
          hess (:,ij+0) = (tg(:,1)-tg(:,2))/(2*dw)
          hess (:,ij+1) = (tg(:,3)-tg(:,4))/(2*dw)
          ts = abs(ts)**2
          testg(  ij+0) = real(ts(1)-ts(2),kind=kind(tmat))/(2*dw)
          testg(  ij+1) = real(ts(3)-ts(4),kind=kind(tmat))/(2*dw)
        end do scan_displacements
        !$omp end parallel do
        !
        !  Gradient at central point
        !
        call gradient(vc,sim,grad)
        !
        if (verbose>=4) then
          write (out,"('Numerical vs analytical gradients (non-zero)')")
          do i=1,lg
            if (maxval(abs((/testg(i),grad(i)/)))>1e-5) then
              write (out,"(5x,i5,2x,f20.12,2x,f20.12,2x,f20.12)") i, testg(i), grad(i), testg(i)-grad(i)
            end if
          end do
        end if
      end subroutine gradient_and_hessian
      !
      subroutine matrix_exponential(w,expw)
        complex(kind(tmat)), intent(in)  :: w   (:,:) ! Matrix to exponentiate; assumed to be small, so that
                                                      ! series expansion converges rapidly
        complex(kind(tmat)), intent(out) :: expw(:,:) ! Matrix exponential
        !
        integer(ik)         :: i
        complex(kind(tmat)) :: poww(bs,bs)
        real(kind(tmat))    :: eps, delta
        !
        poww = 0
        fill_diagonal: do i=1,bs
          poww(i,i) = 1
        end do fill_diagonal
        !
        expw = poww
        exponent_series: do i=1,200
          poww  = matmul(poww,w)/i
          expw  = expw + poww
          eps   = spacing(maxval(abs(expw)))
          delta = maxval(abs(poww))
          if (delta<eps) exit exponent_series
        end do exponent_series
        if (i>100) then
          call stop('biorthogonal_tools%find_optimal_rotation - matrix_exponential did not converge!')
        end if
        !
        !  Make sure the result is unitary
        !
        poww = matmul(expw,conjg(transpose(expw)))
        sub_diagonal: do i=1,bs
          poww(i,i) = poww(i,i) - 1
        end do sub_diagonal
        if (maxval(abs(poww))>spacing(real(1e6,kind(poww)))) then
          write (out,"('exp(w)*conjg(transpose(w)) - diag(1) = ')")
          write (out,*) poww
          call stop('biorthogonal_tools%find_optimal_rotation - matrix_exponential is non-unitary')
        end if
      end subroutine matrix_exponential
      ! 
      !  Apply exponential rotation step to the current transformation matrix
      !
      subroutine apply_rotation(step,vc)
        real(kind(tmat)), intent(in)       :: step(:) ! Finite rotation matrix
                                                      ! 1 = Re(w(1,2)); 2 = Im(w(1,2)); 3 = Re(w(1,3)), etc
        complex(kind(tmat)), intent(inout) :: vc(:,:) ! In: Initial rotation matrix
                                                      ! Out: Updated rotation matrix
        !
        integer(ik)         :: i, j, ij
        complex(kind(tmat)) :: wt(bs,bs)  ! anti-Hermitian rotation exponent
        complex(kind(tmat)) :: vt(bs,bs)  ! finite rotation matrix
        !
        wt = 0
        ij = 1
        fill_rotation_i: do i=1,bs
          fill_rotation_j: do j=i+1,bs
            wt(i,j) = cmplx(step(ij+0),step(ij+1),kind=rk)
            wt(j,i) = -conjg(wt(i,j))
            ij = ij + 2
          end do fill_rotation_j
        end do fill_rotation_i
        call matrix_exponential(wt,vt)
        vc = matmul(vc,vt)
      end subroutine apply_rotation
      !
      !  Calculate similarity function and gradient after perturbing current rotation
      !  matrix by a given w(i,j)
      !
      subroutine displaced_gradient(vt0,i,j,w,sim,grad)
        complex(kind(tmat)), intent(in)  :: vt0(:,:) ! Original rotation matrix
        integer(ik), intent(in)          :: i, j     ! Elements of the exponential rotation to perturb
        complex(kind(tmat)), intent(in)  :: w        ! Value of the perturbation
        complex(kind(tmat)), intent(out) :: sim      ! Similarity at displaced geometry
        real(kind(tmat)), intent(out)    :: grad(:)  ! Gradients of |sim|^2 at displaced geometry
        !
        complex(kind(tmat)) :: wt(bs,bs) ! Displacement matrix
        complex(kind(tmat)) :: vt(bs,bs) ! Rotation matrix
        !
        wt = 0 ; wt(i,j) = w ; wt(j,i) = -conjg(w)
        call matrix_exponential(wt,vt)
        vt = matmul(vt0,vt)
        call gradient(vt,sim,grad)
      end subroutine displaced_gradient
      !
      !  Calculate trace of the occupation matrix and gradient of its square modulus 
      !  with respect to upper diagonal of anti-Hermitian matrix w parameterizing the 
      !  rotation [ie vt = vt0*exp(w)] at w = 0
      !
      !  Note that since w is taken to be anti-Hermitian, trace of the occupation
      !  matrix is -not- an analytical function of the complex elements on the upper
      !  triangular part of matrix w. 
      !
      subroutine gradient(vt,sim,grad)
        complex(kind(tmat)), intent(in)  :: vt(:,:) ! Transformation matrix
        complex(kind(tmat)), intent(out) :: sim     ! Trace of the transformed simularity matrix
        real(kind(tmat)), intent(out)    :: grad(:) ! Gradients of abs(sim)**2. Elements of grad correspond to:
                                            ! 1 = Re(w(1,2)); 2 = Im(w(1,2)); 3 = Re(w(1,3)), etc
        !
        complex(kind(tmat)) :: ut(bs,bs)  ! umat(:,:) transformed according to vt()
        complex(kind(tmat)) :: tt(bs,bs)  ! tmat(:,:) transformed according to tt()
        complex(kind(tmat)) :: gt(bs,bs)  ! Gradients of sim with respect to -general- matrix in the exponential
        !
        integer(ik)         :: i, j, ij
        complex(kind(tmat)) :: df_dr, df_di   ! Gradient of sim wrt Re and Im parts of anti-hermitian matrix element
        complex(kind(tmat)) :: df2_dr, df2_di ! Gradient of |sim|^2; must in fact be real
        real(kind(tmat))    :: eps
        !
        ut = matmul(umat,vt)
        tt = matmul(conjg(transpose(vt)),tmat)
        !
        !  Calculate sim and its unrestricted gradients
        !
        sim = 0
        cmplx_gradient_j: do j=1,bs
          cmplx_gradient_i: do i=1,bs
            gt(i,j) = tt(j,j)*ut(j,i) - ut(i,i)*tt(j,i)
          end do cmplx_gradient_i
          sim = sim + tt(j,j)*ut(j,j)
        end do cmplx_gradient_j
        !
        ij = 1
        real_gradient_i: do i=1,bs
          real_gradient_j: do j=i+1,bs
            df_dr  =  gt(i,j) - gt(j,i)
            df_di  = (gt(i,j) + gt(j,i))*(0,1)
            df2_dr = conjg(sim)*df_dr + sim*conjg(df_dr)
            df2_di = conjg(sim)*df_di + sim*conjg(df_di)
            eps    = 100*spacing(maxval(abs((/sim,df_dr,df_di/))))
            if (any(abs(aimag((/df2_dr,df2_di/)))>eps)) then
              call stop('biorthogonal_tools%find_optimal_rotation - real reduction failed')
            end if
            grad(ij+0) = real(df2_dr,kind=rk)
            grad(ij+1) = real(df2_di,kind=rk)
            ij = ij + 2
          end do real_gradient_j
        end do real_gradient_i
      end subroutine gradient
!   end subroutine find_optimal_rotation_rk
