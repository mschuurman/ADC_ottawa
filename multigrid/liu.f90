 module liu
!
!  Iterative diagonalization, using Liu's variant of Davidson method, with
!  subspace compaction.
!
!  References:
!
!    [DAV] E.R. Davidson, J. Comput. Phys. 17, 87 (1975).
!    [LIU] B. Liu, pp. 49-52, In: C. Moler and I. Shavitt (eds) "Numerical 
!          Algorithms in Chemistry: Algebraic Methods", Laurence Berkeley 
!          Laboratory, Berkeley, California, 1978.
!    [VLP] J.H. van Lenthe, P. Pulay, J. Comput. Chem. 11, 1164 (1990).
!
!  Applied to solving for eigenstates on real-space grids, there is a
!  special quirk: Acting on a vector with the Hamiltonian is always
!  faster than transferring a vector from secondary storage. Therefore,
!  the code below will recompute H |psi> products, wherever practical.
!
!  We'll also try to make use of as many scratch fields as possible,
!  to reduce I/O. So, if there is enough memory for data fields in
!  the calling routine - give that to us!
!
!  We are looking for a n_roots of lowest eigenvalues and eigenvectors, 
!  given the Hamiltonian (defined by function func_pot) and the significance 
!  function (func_wgt) for the initial guess. 
!
!  1. We calculate M random starting guesses, using the weighting function
!     to modulate the amplitudes of the guess vectors. We can't use the 
!     weighting function directly as a guess - this will restrict our 
!     solutions to the same symmetry. 
!
!  2. Construct the optimal orthonormal basis from the trial vectors, using
!     singular value decomposition. Discard all vectors in excess of the
!     limit.
!
!  3. For all new trial vectors, multiply the trial vectors by the Hamiltonian,
!     and remember the results.
!
!  4. Form the subspace Hamiltonian, and diagonalize it
!
!  5. Calculate the residuals, and the correction vectors.
!
!  6. Append the correction vectors to the basis, and proceed to #2.
!
!  The disk space requirements (in vectors) are:
!
!   A. n_basis     for the reduced-space basis 
!   B. n_basis     for the Hamiltonian-basis products
!   C. n_roots     for the correction vectors
!   D. n_basis     for the reorthogonalized basis
!
   use accuracy
   use lapack
   use multigrid
   use qmech
   use fields
   use symmetry
   use timer
   implicit none
   private
   public LUeigenvalues, LUeigenvector
!
!  Important variables:
!
   real(rk), parameter      :: epsilon_root   = 1e-5 ! Convergence criteria for the roots
   real(rk), parameter      :: epsilon_stencil= 1e-5 ! Don't use stencils, if change in eigenvalue is
                                                     ! less than this.
   real(rk), parameter      :: epsilon_ortho  = 1e-1 ! Repeat orthogonalization if overlap eigenvalues
                                                     ! becomes smaller than this
   real(rk), parameter      :: epsilon_over   = 1e-2 ! Rebuild basis set if basis functions overlap
                                                     ! gets above this.
   integer(ik), parameter   :: verbose = 1           ! Verbosity level
   real(rk)                 :: mass                  ! Mass of the particle - for the preconditioner
   integer(ik)              :: n_roots               ! Number of roots we calculate residuals for
   integer(ik)              :: target_roots          ! Number of roots we need to converge
   integer(ik)              :: n_converge            ! Number of roots we are actually converging, as a
                                                     ! precaution against false convergence
   integer(ik)              :: n_basis               ! Current number of basis vectors in the reduced space
   integer(ik)              :: n_corr                ! Number of correction vectors we built
   integer(ik)              :: n_new                 ! Number of eigenvectors in ind_new
   integer(ik)              :: max_basis             ! Maximum number of basis vectors we'll use
!
   integer(ik)              :: n_scr                 ! Number of scratch muntigrid slots available
   integer(ik), pointer     :: scr(:)                ! Multigrid slots, avalable for our use
   integer(ik)              :: f_pot                 ! Multigrid slot, containing the potential
   integer(ik)              :: f_pre                 ! Multigrid slot, containing the preconditioner
!
   integer(ik), allocatable :: ind_bas  (:)          ! Positions of basis vectors in the scratch file
   integer(ik), allocatable :: ind_corr (:)          ! ... of the correction vectors
   integer(ik), allocatable :: ind_new  (:)          ! ... of the new basis vectors
   logical, allocatable     :: disk_used(:)          ! Flags for occupied disk slots.
!
   complex(rk), allocatable :: ovl     (:,:)         ! Overlap integrals
   complex(rk), allocatable :: h       (:,:)         ! Hamiltonian matrix in the reduced space
   complex(rk), allocatable :: vec     (:,:)         ! Subspace eigenvectors
   real(rk), allocatable    :: e       (:)           ! Subspace eigenvalues
   real(rk), allocatable    :: old_e   (:)           ! Subspace eigenvalues
   real(rk), allocatable    :: diff_e  (:)           ! Change in subspace eigenvalues
   logical, allocatable     :: use_stencil(:)        ! Use error diffusion stencils for this root
   real(rk), allocatable    :: res_norm(:)           ! Norms of the residuals
!
   contains
!
! Driver routines, visible from outside
!
   subroutine LUeigenvalues(mass_,scr_,eval,func_pot,func_wgt,func_pre)
     real(rk), intent(in)            :: mass_   ! Particle mass
     integer(ik), intent(in), target :: scr_(:) ! Scratch multigrid slots. The minimum of three slots is
                                                ! required, but more is better.
     real(rk), intent(out)   :: eval(:)  ! Space for the eigenvalues
     complex(rk), external   :: func_pot ! Function, generating multiplicative potential for the 
                                         ! Hamiltonian
     complex(rk), external   :: func_wgt ! Weighting function, large in the parts of space where 
                                         ! the solutions are expected to be large
     complex(rk), external   :: func_pre ! Function, generating the potential part of the preconditioner 
                                         ! for Davidson iteration. Most likely, this will be the same 
                                         ! as the multiplicative potential, but it does not need to. 
                                         ! The code in the solver will take care of adding the kinetic 
                                         ! energy density.
     !
     integer(ik) :: iter, root
     logical     :: depend, have_overlaps
     !
     call TimerStart('Block-Davidson solver')
     !
     !  Set global stuff for the solver
     !
     n_scr = size(scr_)                 ! How much scratch do we have?
     if (n_scr<5) then
       write (out,"(/'Davidson/Liu solver requires minimum of 5 scratch vectors. "// &
                  "Only ',i4,' were provided')") n_scr
       stop 'LUeigenvalues - no memory'
     end if
     !
     scr   => scr_                      ! Make scratch indices global
     mass  = mass_                      ! Particle mass for the kinetic part of precoditioner
     !
     target_roots = size(eval)          ! Number of roots desired
     n_converge   = target_roots + 1    ! Converging one extra root will help to
                                        ! avoid false convergence
     max_basis    = 3*n_converge + 4    ! Number of vectors in the subspace basis
     n_roots      = max_basis / 2       ! Number of solutions to calculate on each step
     if (verbose>=1) then
       write (out,"(/t5,'Davidson-Liu solver with error diffusion preconditioner'/)") 
       write (out,"(' Number of roots requested = ',i4)") target_roots
       write (out,"(' Number of roots converged = ',i4)") n_converge
       write (out,"(' Number of roots tracked   = ',i4)") n_roots
       write (out,"(' Max. subspace size        = ',i4)") max_basis
       write (out,"(' Required convegence (E)   = ',e10.3)") epsilon_root
       write (out,"()")
     end if
     !
     !  Bootstrap. Generate the potential, preconditioner, and initial guess.
     !             The potential and preconditioner will reduce the number
     !             of scratch slots by 2.
     !
     f_pot = scr(n_scr) ; n_scr = n_scr - 1
     f_pre = scr(n_scr) ; n_scr = n_scr - 1
     !
     ! Prepare the potential and preconditioner. There is a trade-off in 
     ! building them here, instead of every time they are needed. The on-
     ! demand construction can exploit the knowledge of the wavefuntion 
     ! extent, and avoid building the potential in areas where wavefunction 
     ! is zero. On the other hand, if the wavefunction extends over much 
     ! of the box, and evaluating the potential is expensive, we can save 
     ! a great deal of work by building the potential upfront.
     !
     call FieldInit(f_pot,func_pot)
     call FieldInit(f_pre,func_pre)
     !
     call allocate_memory
     call initial_guess(func_wgt)
     !
     iter        = 0
     e           = 0
     use_stencil = .true.
     depend      = .true.
     iteration: do
       iter = iter + 1
       !
       hmat_ortho: do
         !
         !  Basis overlap integrals are also calculated by subspace_hamiltonian.
         !  If we end up here due to an abort in subspace_hamiltonian, we can
         !  save a bit of work in orthonormalize_basis
         !
         have_overlaps = .false.
         !
         !  As a last resort, we'll do singular value decomposition of
         !  the basis set, and extract orthonormal vectors. This destroys
         !  the eigenstates, so it should only be done if everything
         !  else fails.
         !
         full_ortho: do while(depend)
           call orthonormalize_basis(depend,have_overlaps)
           have_overlaps = .false.
           if (.not.depend) exit full_ortho
           if (verbose>=1) then
             write (out,"(/'Subspace basis approaching linear dependence, re-orthogonalize')")
           end if
         end do full_ortho
         !
         !  Construct subspace Hamiltonian and subspace overlap matrix.
         !  If non-diagonal subspace overlaps become large, construction
         !  of the Hamiltonian will abort, and the basis will be reconstructed
         !  in orthonormalize_basis
         !
         call subspace_hamiltonian(depend)
         have_overlaps = .true.
         if (.not.depend) exit hmat_ortho
       end do hmat_ortho
       !
       old_e  = e
       call subspace_diagonalize
       diff_e = e - old_e
       !
       if (verbose>=1) then
         write (out,"(/'Iteration ',i4,' e (dec): ',(t25,3(f15.7,1x,a1,'(',e10.3,') ')))") &
                iter, (e(root), merge('S',' ',use_stencil(root)), &
                      -diff_e(root), root=1,n_roots)
       end if
       !
       !  Build eigenvectors for the n_roots lowest states. The new 
       !  eigenvectors are stored in ind_new(1:n_new)
       !
       call build_eigenvectors
       !
       !  Convergence test. We have to make sure enough iterations were
       !  made before making the test. Also, error diffusion limits
       !  variational freedom, so that if we used error diffusion
       !  preconditioner for any of the roots, we can't declare it
       !  converged.
       !
       if (iter>=3 .and. all(abs(diff_e(1:n_converge))<epsilon_root) .and. &
                  .not. any(use_stencil(1:n_converge)) ) then
           if (verbose>=1) then
             write(out,"(/' Davidson/Liu iterations converged after ',i4,' cycles'/)") &
                   iter
           end if
           exit iteration
       end if
       !
       !  Update stencilling flags
       !
       use_stencil = (iter<=1) .or. .not.(-diff_e<epsilon_stencil)
       !
       !  Build residuals - our basis set expansion vectors
       !
       call davidson_residuals
       !
       !  New correction vectors need to be orthogonalized with respect to
       !  the approximate solution vectors we calculated above. Because we 
       !  want to to keep the approximate eigenfunctions unperturbed, we'll
       !  ll use Gramm-Shmidt procedure here. Because linear dependencies 
       !  may be severe, this may need to be done more than once (see [LIU])
       !
       depend = .true.
       res_ortho: do while(depend)
         call orthonormalize_residuals(depend)
         if (.not.depend) exit res_ortho
       ! if (verbose>=1) then
       !   write (out,"(t4,'Residuals are linearly dependent, re-orthogonalize')") 
       ! end if
       end do res_ortho
       if (n_corr==0) then
         !
         !  Correction vectors collapced, assume convergence
         !
         if (verbose>=1) then
           write (out,"(/'All correction vectors vanish, assuming convergence')")
         end if
         exit iteration
       end if
       !
       !  The new basis set consists of the few lowest eigenstates of the
       !  reduced Hamiltonian, plus the residuals
       !
       n_basis = n_new + n_corr
       ind_bas(1:n_basis) = (/ ind_new(1:n_new), ind_corr(1:n_corr) /)
       call reset_disk(ind_bas(1:n_basis))
     end do iteration
     !
     !  Return the eigenvalues. The eigenvectors are stored in 
     !  ind_new(1:n_new), and can be retrived with LUeigenvectors
     !
     eval(:) = e(1:target_roots)
     !
     call TimerStop('Block-Davidson solver')
     call TimerReport
   end subroutine LUeigenvalues
!
   subroutine LUeigenvector(root,psi)
     integer(ik), intent(in) :: root  ! Desired root. Values less than 1
                                      ! will close all files and discard
                                      ! the roots.
     integer(ik), intent(in) :: psi   ! Multigrid slot for the root
     !
     if (root>=1) then
       if (root>n_new) then
         write (out,"('LUeigenvectors called for root ',i4,', but only '"// &
                    ",i4,' roots are available')") root, n_new
         stop 'LUeigenvectors - bad root'
       end if
       call FieldIO('READ',psi,slot=ind_new(root))
     else
       call release_memory
       n_new = 0
     end if
     !
   end subroutine LUeigenvector
!
! Internal service routines
!

   !
   ! Initial guess
   !
   subroutine initial_guess(f_wgt)
     complex(rk), external :: f_wgt
     !
     integer(ik)           :: n_guess
     integer(ik)           :: bi, slot
     real(rk)              :: guess_norm
     integer(ik)           :: symbas, symbas_count
     character(len=20)     :: symbas_name
     logical               :: randomize
     integer(ik)           :: psi1, psi2 ! Multigrid slots
     !
     call TimerStart('Davidson initial guess')
     symbas_count = SMgetTransformCount()
     n_guess = max(n_roots,min(symbas_count,max_basis))
     if (verbose>=1) then
       write (out,"(' Number of initial guess vectors   : ',i4)") n_guess
       write (out,"(' Number of symmetry representations: ',i4)") symbas_count
       write (out,"(' Number of randomized guesses      : ',i4)") max(0,n_guess-symbas_count)
       write (out,"()")
     end if
     if (n_guess<symbas_count) then
       write (out,"(/'WARNING - Convergence may be slow due to missing symmetries!'/)")
     end if
     !
     ! We have to cold-initialize the field at most once - this kicks in the 
     ! the cut-offs for the grid dimensions, and should make everything else
     ! faster.
     !
     psi1 = scr(1) ; psi2 = scr(2)
     call FieldInit(psi1,f_wgt)
     symbas    = 1
     randomize = .false.
     guess: do bi=1,n_guess
       call FieldCopy   (psi1,psi2)
       !
       !  Apply symmetry transformation, and possibly randomize
       !
       symbas_name = SMsetTransform(symbas)
       call FieldProcess(psi2,SMtransformField)
       if (randomize) then
         call FieldProcess(psi2,FLrandom)
       end if
       !
       !  Renormalize
       !
       call QMNormalize(psi2,1.0_rk,guess_norm)
       if (verbose>=1) then
         write (out,"(t5,'Guess vector ',i3,1x,a10,': Initial norm = ',f14.7)") &
                bi, '('//trim(symbas_name)//')', guess_norm
       end if
       slot = get_disk_slot()
       call FieldIO('WRITE',psi2,slot=slot)
       ind_bas(bi) = slot
       !
       !  Advance to next symop
       !
       symbas      = symbas + 1
       if (symbas>symbas_count) then
         randomize = .true.
         symbas    = 1
       end if
     end do guess
     n_basis = n_guess
     call TimerStop('Davidson initial guess')
   end subroutine initial_guess
   !
   ! Construct the optimal orthonormal basis
   !
   subroutine orthonormalize_basis(depend,have_overlap)
     logical, intent(out) :: depend       ! True if basis was linearly dependent
     logical, intent(in)  :: have_overlap ! True if overlap is already calculated
     !
     integer(ik) :: bi, bj
     integer(ik) :: new_bas, bas, slot
     real(rk)    :: e(n_basis), norm
     integer(ik) :: lblock                ! Number of vectors, processed in one pass
     integer(ik) :: psi1, psi2
     integer(ik) :: block_1, block_n      ! First and last vectors in a block
     !
     call TimerStart('Davidson basis SVD')
     !
     ! At this point, only the subspace basis is needed, discard
     ! all other vectors
     !
     call reset_disk(ind_bas(1:n_basis))
     !
     ! Assign scratch vectors
     !
     psi2   = scr(n_scr)           ! Use last scratch vector for kets
     lblock = min(n_scr-1,n_basis) ! Use the rest for bras
     !
     ! Calculate overlap integrals of the basis functions
     !
     if (.not.have_overlap) then
       call subspace_s_matrix(ovl,depend)
     end if
     !
     if (verbose>=3) then
       write (out,"(/t5,'Overlap integrals of the subspace basis functions'/)") 
       prt_i: do bi=1,n_basis
         write (out,"(5x,i4,2x,(t13,3(2x,2f9.5)))") bi, ovl(1:n_basis,bi)
         write (out,"()")
       end do prt_i
     end if
     !
     ! Choose the maximum spanning subspace. Replace the overlaps with
     ! the eigenvectors
     !
     call lapack_heev(ovl(1:n_basis,1:n_basis),e(1:n_basis))
     !
     ! *HEEV eigenvalues are in the ascending order, while we prefer
     ! the descending order. So, swap both the eigenvalues and eigenvectors
     !
     e  (  1:n_basis) = e  (  n_basis:1:-1)
     ovl(:,1:n_basis) = ovl(:,n_basis:1:-1)
     !
     ! We can use up to max_basis basis functions. Choose functions
     ! which span as much of the original space as possible. If we have
     ! less than that to begin with, we will simply reorthogonalize our
     ! basis set - this preserves numerical stability.
     !
     depend = .false.
     choose_newbas: do bas=1,n_basis
       !
       if (verbose>=3) then
         write (out,"(/t4,'Root ',i4,' eps = ',f12.4,' : ',(t34,4(f11.5)))") &
                bas, e(bas), ovl(1:n_basis,bas)
       end if
       !
       ! No matter what, we must choose at least n_roots functions. As long
       ! as no severe linear dependencies are present, we'll continue up to 
       ! n_basis
       !
       if (bas>n_roots .and. e(bas)<epsilon_ortho*e(1)) &
         exit choose_newbas
       !
       ! If we must continue for linearly-dependent functions, request 
       ! a re-run.
       !
       if (e(bas)<epsilon_ortho*e(1)) depend = .true.
     end do choose_newbas
     new_bas = bas - 1
     !
     if (verbose>=2) then
       write (out,"(/t3,'Selected ',i4,' basis vectors out of ',i4)") new_bas, n_basis
       write (out,"(/t3,'Overlap eigenvalues: ',(t25,5f12.4))") e(1:n_basis)
     end if
     !
     ! Rebuild basis functions we selected in the loop above
     !
     build_newbas: do block_1=1,new_bas,lblock
       block_n = min(block_1+lblock-1,new_bas)
       !
       build_zero_block: do bas=block_1,block_n
         psi1 = scr(bas+1-block_1)
         call FieldZero(psi1)
       end do build_zero_block
       !
       reortho: do bj=1,n_basis
         call FieldIO('READ',psi2,slot=ind_bas(bj))
         !
         reortho_block: do bas=block_1,block_n
           psi1 = scr(bas+1-block_1)
           call FieldAXPY(ovl(bj,bas),src=psi2,dst=psi1)
         end do reortho_block
       end do reortho
       !
       save_block: do bas=block_1,block_n
         psi1 = scr(bas+1-block_1)
         !
         ! Sanity check - the norm must be close to one.
         !
         call QMNormalize(psi1,1.0_rk,norm)
         if (verbose>=3) then
           write (out,"(t5,'Initial norm = ',f14.7)") norm
         end if
         slot = get_disk_slot()
         call FieldIO('WRITE',psi1,slot=slot)
         ind_new(bas) = slot
       end do save_block
       !
     end do build_newbas
     !
     ! Replace the basis with the new vectors, and we can go to the next cycle
     !
     n_basis = new_bas
     ind_bas(1:n_basis) = ind_new(1:n_basis)
     call reset_disk(ind_bas(1:n_basis))
     !
     call TimerStop('Davidson basis SVD')
   end subroutine orthonormalize_basis
   !
   !  Gramm-Shmidt orthonormalize residuals to all basis vectors
   !  We know what our eigenstates span the subspace of the basis.
   !  To improve convergence, we'll orthonormalize to the basis,
   !  -not- to the subspace. If this introduces significant numerical 
   !  noise, subspace_hamiltonian() will catch it, and rebuild the
   !  basis.
   !
   subroutine orthonormalize_residuals(depend)
     logical, intent(out) :: depend ! True if significant linear dependence 
                                    ! was found
     !
     integer(ik) :: corr, ci, bj, cj
     complex(rk) :: ovl
     real(rk)    :: norm
     integer(ik) :: lblock, block_1, block_n
     integer(ik) :: psi1, psi2
     !
     call TimerStart('Davidson residuals orthogonalization')
     ! 
     ! Partition scratch area
     !
     lblock = n_scr - 1
     !
     ! It is possible what some of the residuals will vanish
     ! completely after orthogonalization. Therefore, we'll
     ! keep a separate counter for the output residuals
     !
     depend = .false.
     corr   = 0
     cor_i: do block_1=1,n_corr,lblock
       block_n = min(block_1+lblock-1,n_corr)
       !
       if (verbose>=3) write (out,"()")
       !
       cor_read_block: do ci=block_1,block_n
         psi1 = scr(ci+1-block_1)
         call FieldIO('READ',psi1,slot=ind_corr(ci))
       end do cor_read_block
       !
       !  Orthogonalize with respect to the basis
       !
       psi2 = scr(n_scr)
       bas_j: do bj=1,n_basis
         call FieldIO('READ',psi2,slot=ind_bas(bj))
         !
         ortho_bas_block: do ci=block_1,block_n
           psi1 = scr(ci+1-block_1)
           ovl  = FieldConjgIntegrate(psi2,psi1)
           call FieldAXPY(-ovl,src=psi2,dst=psi1)
           !
           if (verbose>=3) then
             write (out,"(t4,'<cor ',i4,'| bas ',i4,'> = ',2f12.6)") &
                    ci, bj, ovl
           end if
         end do ortho_bas_block
       end do bas_j
       !
       !  Orthogonalize with respect to other corection vectors
       !  This can get a little tricky, as some of the vectors may 
       !  be in the memory buffer, or may vanish altogether. 
       !
       !  First, we'll do vectors which have already been spilled out
       !
       spilled_cor_j: do cj=1,corr
         call FieldIO('READ',psi2,slot=ind_corr(cj))
         !
         spilled_cor_block: do ci=block_1,block_n
           psi1 = scr(ci+1-block_1)
           ovl  = FieldConjgIntegrate(psi2,psi1)
           call FieldAXPY(-ovl,src=psi2,dst=psi1)
           !
           if (verbose>=3) then
             write (out,"(t4,'<cor ',i4,'| cor ',i4,'> = ',2f12.6)") &
                    cj, ci, ovl
           end if
         end do spilled_cor_block 
       end do spilled_cor_j
       !
       !  Now, we can do vectors still in memory, retiring them in the process
       !
       buf_cor_block: do ci=block_1,block_n
         psi1 = scr(ci+1-block_1)
         !
         !  Renormalize
         !
         call QMNormalize(psi1,1.0_rk,norm)
         !
         if (verbose>=2) then
           write (out,"(/4x,'Linearly independent norm of correction ',i4,': ',f12.6)") &
                  ci, norm
         end if
         !
         if (norm<=0) cycle buf_cor_block ! Correction vector vanished
         !
         if (abs(norm-1.0_rk)>=epsilon_ortho) depend = .true.
         !
         corr = corr + 1
         call FieldIO('WRITE',psi1,slot=ind_corr(corr))
         !
         !  Forward orthogonalize to remaining vectors in memory
         !
         buf_corr_forward: do cj=ci+1,block_n
           psi2 = scr(cj+1-block_1)
           ovl  = FieldConjgIntegrate(psi1,psi2)
           call FieldAXPY(-ovl,src=psi1,dst=psi2)
           !
           if (verbose>=3) then
             write (out,"(t4,'<cor ',i4,'| cor ',i4,'> = ',2f12.6)") &
                    ci, cj, ovl
           end if
         end do buf_corr_forward
       end do buf_cor_block
     end do cor_i
     !
     n_corr = corr
     !
     call TimerStop('Davidson residuals orthogonalization')
   end subroutine orthonormalize_residuals
   !
   ! Calculate subspace Hamiltonian matrix. For practical reasons, we'll
   ! split the calculation into five passes over the data: overlap matrix.
   ! potential matrix elements, and three passes for the Laplacian matrix
   ! elements (evaluated by splitting the Laplacian into three momentum
   ! contributions). Hopefully, this will make the code easier to read,
   ! at the expense of some additional I/O.
   !
   subroutine subspace_hamiltonian(depend)
     logical, intent(out)  :: depend    ! Construction of Hamiltonian aborted due 
                                        ! to linear dependencies in the basis
     !
     complex(rk) :: h_v(n_basis,n_basis)   ! Potential contribution to the Hamiltonian
     complex(rk) :: h_p(n_basis,n_basis,3) ! p_x,p_y,p_z contribution to the Hamiltonian
     integer(ik) :: ic, nc
     !
     call TimerStart('Davidson subspace Hamiltonian')
     !
     !  Calculate the ingredients
     !
     call subspace_s_matrix(ovl,depend)
     if (depend) then
       call TimerStop('Davidson subspace Hamiltonian')
       return
     end if
     call subspace_v_matrix(h_v)
     nc = FieldComponentCount()
     p_x: do ic=1,nc
       call subspace_pp_matrix(ic,h_p(:,:,ic))
     end do p_x
     !
     !  Combine the ingredients to form the complete Hamiltonian matrix
     !  Note that the sign in front of the h_p term is (+), not (-) as
     !  one might expect. The sign change is due to the fact that we
     !  have used \hat{p} on the bra vectors, instead of the \hat{p}^\dagger
     !  as we should have.
     !
     h = 0
     h(1:n_basis,1:n_basis) = (0.5_rk/mass)*sum(h_p(:,:,1:nc),dim=3) + h_v
     !
     call TimerStop('Davidson subspace Hamiltonian')
   end subroutine subspace_hamiltonian
   !
   !  Calculate overlap matrix within the subspace
   !
   subroutine subspace_s_matrix(s,depend)
     complex(rk), intent(out) :: s(:,:) ! Calculated overlap matrix
     logical, intent(out)     :: depend ! Does the overlap matrix contain 
                                        ! linear dependencies?
     !
     integer(ik) :: ket      ! Ket vector index
     integer(ik) :: f_ket    !  ... the field used to store it
     integer(ik) :: slot_ket !  ... and the corresponding disk slot
     integer(ik) :: bra      ! Bra vector index
     integer(ik) :: f_bra    !  ... the field used to store it
     integer(ik) :: slot_bra !  ... and the corresponding disk slot
     !
     integer(ik) :: lblock   ! Size of the ket vectors block used for I/O buffer
     integer(ik) :: ket_1    ! First ket vector in the current block
     integer(ik) :: ket_n    ! Last ket vector in the current block
     !
     call TimerStart('Davidson Overlap matrix')
     !
     ! Partition the scratch area
     !
     f_bra  = scr(n_scr)   ! Bra vector is always at the same place
     lblock = n_scr-1      ! Rest of the buffer is for the kets
     !
     s = 0
     ket_loop: do ket_1=1,n_basis,lblock
       ket_n = min(ket_1+lblock-1,n_basis)
       !
       ! Load as many ket vectors as memory will hold
       !
       ket_load_block: do ket=ket_1,ket_n
         f_ket    = scr(ket+1-ket_1)
         slot_ket = ind_bas(ket)
         call FieldIO('READ',f_ket,slot=slot_ket)
       end do ket_load_block
       !
       ! With the ket buffer full, go over all bras
       !
       bra_loop: do bra=1,n_basis
         slot_bra = ind_bas(bra)
         call FieldIO('READ',f_bra,slot=slot_bra)
         !
         braket_integrate: do ket=ket_1,ket_n
           f_ket = scr(ket+1-ket_1)
           s(bra,ket) = FieldConjgIntegrate(f_bra,f_ket)
           !
         end do braket_integrate
       end do bra_loop
       !
     end do ket_loop
     !
     !  Check overlap matrix for linear dependencies and Hermiticity
     !
     depend = .false.
     ket_check: do ket=1,n_basis
       if (abs(s(ket,ket)-1.0_rk)>=epsilon_over) then
         if (verbose>=2) then
           write (out,"(t4,'Vector ',i4,' unnormalized, rebuild')") bra
         end if
         depend = .true.
       end if
       bra_check: do bra=ket,n_basis
         if (bra/=ket .and. abs(s(bra,ket))>=epsilon_over) then
           if (verbose>=2) then
             write (out,"(t4,'Vectors ',i4,' and ',i4,' non-orthogonal, rebuild')") &
                    bra, ket
           end if
           depend = .true.
         end if
         if (abs(s(bra,ket)-conjg(s(ket,bra)))>100*spacing(1.0_rk)) then
           write (out,"(3x,'Bad S(',i3,',',i3,'): ',2(f14.7,1x,f14.7,3x),g14.7,1x,g14.7)") &
                  bra, ket, s(bra,ket), s(ket,bra), s(bra,ket)-conjg(s(ket,bra))
         end if
       end do bra_check
     end do ket_check
     !
     s = 0.5_rk * (s + transpose(conjg(s)))
     !
     call TimerStop('Davidson Overlap matrix')
   end subroutine subspace_s_matrix
   !
   !  Calculate matrix elements of the potential within the subspace
   !
   subroutine subspace_v_matrix(v)
     complex(rk), intent(out) :: v(:,:)   ! Calculated potental matrix elements
     !
     integer(ik) :: ket      ! Ket vector index
     integer(ik) :: f_ket    !  ... the field used to store it
     integer(ik) :: slot_ket !  ... and the corresponding disk slot
     integer(ik) :: bra      ! Bra vector index
     integer(ik) :: f_bra    !  ... the field used to store it
     integer(ik) :: slot_bra !  ... and the corresponding disk slot
     !
     integer(ik) :: f_scr    ! Scratch field
     !
     integer(ik) :: lblock   ! Size of the ket vectors block used for I/O buffer
     integer(ik) :: ket_1    ! First ket vector in the current block
     integer(ik) :: ket_n    ! Last ket vector in the current block
     !
     call TimerStart('Davidson Potential matrix')
     !
     ! Partition the scratch area
     !
     f_bra  = scr(n_scr)     ! Field for the bra vector. We also use this
     f_scr  = f_bra          ! field for scratch while loading kets
     lblock = n_scr-1        ! Rest of the buffer is for the v|ket> products
     !
     v = 0
     ket_loop: do ket_1=1,n_basis,lblock
       ket_n = min(ket_1+lblock-1,n_basis)
       !
       ! Load as many ket vectors as memory will hold, multiplying
       ! them with the potential along the way.
       !
       ket_load_block: do ket=ket_1,ket_n
         f_ket    = scr(ket+1-ket_1)
         slot_ket = ind_bas(ket)
         call FieldIO('READ',f_scr,slot=slot_ket)
         call FieldZero(f_ket)
         call FieldMulAdd(f_pot,f_scr,f_ket)
       end do ket_load_block
       !
       ! With the ket buffer full, go over all bras
       !
       bra_loop: do bra=1,n_basis
         slot_bra = ind_bas(bra)
         call FieldIO('READ',f_bra,slot=slot_bra)
         !
         braket_integrate: do ket=ket_1,ket_n
           f_ket = scr(ket+1-ket_1)
           v(bra,ket) = FieldConjgIntegrate(f_bra,f_ket)
           !
         end do braket_integrate
       end do bra_loop
       !
     end do ket_loop
     !
     ! Check matrix elements for Hermiticity
     !
     ket_check: do ket=1,n_basis
       bra_check: do bra=ket,n_basis
         if (abs(v(bra,ket)-conjg(v(ket,bra)))>100*spacing(1.0_rk)) then
           write (out,"(3x,'Bad V(',i3,',',i3,'): ',2(f14.7,1x,f14.7,3x),g14.7,1x,g14.7)") &
                  bra, ket, v(bra,ket), v(ket,bra), v(bra,ket)-conjg(v(ket,bra))
         end if
       end do bra_check
     end do ket_check
     !
     v = 0.5_rk * (v + transpose(conjg(v)))
     !
     call TimerStop('Davidson Potential matrix')
   end subroutine subspace_v_matrix
   !
   !  Calculate matrix elements of a given Laplacian component
   !
   subroutine subspace_pp_matrix(ic,pp)
     integer(ik), intent(in)  :: ic      ! Cartesian direction
     complex(rk), intent(out) :: pp(:,:) ! Calculated matrix elements
     !
     integer(ik) :: ket      ! Ket vector index
     integer(ik) :: f_ket    !  ... the field used to store it
     integer(ik) :: slot_ket !  ... and the corresponding disk slot
     integer(ik) :: bra      ! Bra vector index
     integer(ik) :: f_bra    !  ... the field used to store it
     integer(ik) :: slot_bra !  ... and the corresponding disk slot
     !
     integer(ik) :: f_scr    ! Scratch field
     !
     integer(ik) :: lblock   ! Size of the ket vectors block used for I/O buffer
     integer(ik) :: ket_1    ! First ket vector in the current block
     integer(ik) :: ket_n    ! Last ket vector in the current block
     !
     call TimerStart('Davidson Laplacian component matrix')
     !
     ! Partition the scratch area
     !
     f_bra  = scr(n_scr)   ! Field for the bra vector. 
     f_scr  = scr(n_scr-1) ! Scratch field, used to construct (d/dx)|ket>
                           ! and (d/dx)|bra> vectors
     lblock = n_scr-2      ! Rest of the buffer is for the (d/dx)|ket> 
     !
     ! write (out,"('f_bra = ',i4,' f_scr = ',i4)") f_bra, f_scr
     ! write (out,"('scr = ',20i4)") scr(1:lblock)
     !
     pp = (99.0_rk, 1.0_rk)
     ket_loop: do ket_1=1,n_basis,lblock
       ket_n = min(ket_1+lblock-1,n_basis)
       !
       ! Load as many ket vectors as memory will hold, applying the
       ! (d/dx) operator to them along the way.
       !
       ket_load_block: do ket=ket_1,ket_n
         f_ket    = scr(ket+1-ket_1)
         slot_ket = ind_bas(ket)
         call FieldIO('READ',f_scr,slot=slot_ket)
         call FieldGradientComponentRight(ic,f_scr,f_ket)
         ! write (out,"('ket ',i4,' is at ',i4)") ket, f_ket
       end do ket_load_block
       !
       ! With the ket buffer full, go over all bras. We need
       ! to apply the same operator to the bras.
       !
       bra_loop: do bra=1,n_basis
         slot_bra = ind_bas(bra)
         call FieldIO('READ',f_scr,slot=slot_bra)
         call FieldGradientComponentRight(ic,f_scr,f_bra)
         ! write (out,"('bra ',i4,' is at ',i4)") bra, f_bra
         !
         braket_integrate: do ket=ket_1,ket_n
           f_ket = scr(ket+1-ket_1)
           pp(bra,ket) = FieldConjgIntegrate(f_bra,f_ket)
           ! write (out,"('ket ',i4,' is at ',i4)") ket, f_ket
           ! write (out,"(3x,i3,1x,i3,1x,i1,2x,2f16.12)") bra, ket, ic, pp(bra,ket)
           ! write (out,"(                 t15,2f16.12)") FieldConjgIntegrate(f_ket,f_bra)
           !
         end do braket_integrate
       end do bra_loop
       !
     end do ket_loop
     !
     ! Check matrix elements for Hermiticity
     !
     ket_check: do ket=1,n_basis
       bra_check: do bra=ket,n_basis
         if (abs(pp(bra,ket)-conjg(pp(ket,bra)))>100*spacing(1.0_rk)) then
           write (out,"(3x,'Bad PP(',i3,',',i3,',',i1,'): ',2(f14.7,1x,f14.7,3x),g14.7,1x,g14.7)") &
                  bra, ket, ic, pp(bra,ket), pp(ket,bra), pp(bra,ket)-conjg(pp(ket,bra))
         end if
       end do bra_check
     end do ket_check
     !
     pp = 0.5_rk * (pp + transpose(conjg(pp)))
     !
     call TimerStop('Davidson Laplacian component matrix')
   end subroutine subspace_pp_matrix
   !
   ! Diagonalize subspace Hamiltonian matrix. This will solve generalized
   ! eigenvalue problem, to tolerate small losses of orthonormality in the
   ! basis set. The (classic) algorithm is:
   !
   !  -    -1/2    -1/2
   !  H = S     H S
   !
   !  - -   -
   !  H C = C E
   !
   !       -1/2 -
   !  C = S     C
   !
   subroutine subspace_diagonalize
     integer(ik) :: root
     complex(rk) :: smhalf(1:n_basis,1:n_basis)
     !
     !  Calculate S^(-1/2)
     !
     smhalf = ovl(1:n_basis,1:n_basis)
     call lapack_ginverse(smhalf,power_=0.25_rk)
     !
     if (verbose>=3) then
       write (out,"(/t15,'S^(-1/2) is:'/)") 
       p_shalf: do root=1,n_basis
         write (out,"(1x,i4,(t6,6(2f8.4,2x)))") root, smhalf(root,:)
         write (out,"()")
       end do p_shalf
     end if
     !
     vec(1:n_basis,1:n_basis) = matmul(matmul(smhalf,h(1:n_basis,1:n_basis)),smhalf)
     call lapack_heev(vec(1:n_basis,1:n_basis),e(1:n_basis))
     !
     if (verbose>=2) then
       write (out,"(/t3,'Subspace eigenvalues: ',(t26,5f14.7))") e(1:n_basis)
     end if
     if (verbose>=4) then
       write (out,"(/t4,'Orthogonal subspace eigenvectors:'/)")
       p_roots1: do root=1,n_basis
         write (out,"(2x,i4,1x,f12.5,2x,(t23,3(2f10.6,2x)))") &
                root, e(root), vec(1:n_basis,root)
         write (out,"()")
       end do p_roots1
     end if
     !
     !  Convert back to non-orthogonal basis
     !
     vec(1:n_basis,1:n_basis) = matmul(smhalf,vec(1:n_basis,1:n_basis))
     !
     if (verbose>=2) then
       write (out,"(/t4,'Non-orthogonal subspace eigenvectors:'/)")
       p_roots2: do root=1,n_basis
         write (out,"(2x,i4,1x,f12.5,2x,(t23,3(2f10.6,2x)))") &
                root, e(root), vec(1:n_basis,root)
         write (out,"()")
       end do p_roots2
     end if
   end subroutine subspace_diagonalize
   !
   ! Construct spatial eigenvectors for the low-energy eigenstates of the
   ! subspace Hamiltonian. These vectors will form the core of the subspace
   ! basis set for the next iteration
   !
   subroutine build_eigenvectors
     integer(ik) :: bj
     integer(ik) :: root, slot_root
     real(rk)    :: norm
     integer(ik) :: psi1, psi2
     integer(ik) :: lblock, block_1, block_n
     !
     call TimerStart('Davidson eigenvectors')
     !
     !  Partition scratch area
     !
     psi2   = scr(n_scr)
     lblock = n_scr - 1
     !
     if (verbose>=2) then
       write (out,"(/t4,'Building spatial vectors for ',i4,' reduced eigenstates')") &
              n_roots
     end if
     !
     if (n_roots>n_basis) then
       stop 'liu%build_eigenvectors - not enough roots'
     end if
     !
     roots: do block_1=1,n_roots,lblock
       block_n = min(block_1+lblock-1,n_roots)
       !
       roots_block_zero: do root=block_1,block_n
         psi1 = scr(root+1-block_1)
         call FieldZero(psi1)
       end do roots_block_zero
       !
       bas_j: do bj=1,n_basis
         call FieldIO('READ',psi2,slot=ind_bas(bj))
         !
         roots_block_add: do root=block_1,block_n
           psi1 = scr(root+1-block_1)
           call FieldAXPY(vec(bj,root),src=psi2,dst=psi1)
         end do roots_block_add
       end do bas_j
       !
       roots_block_save: do root=block_1,block_n
         psi1 = scr(root+1-block_1)
         call QMNormalize(psi1,1.0_rk,norm)
         if (verbose>=2) then
           write (out,"(t8,'State ',i4,' norm = ',f12.6)") root, norm
         end if
         slot_root = get_disk_slot()
         ind_new(root) = slot_root
         call FieldIO('WRITE',psi1,slot=slot_root)
       end do roots_block_save
     end do roots
     n_new = n_roots
     call TimerStop('Davidson eigenvectors')
   end subroutine build_eigenvectors
   !
   !  Calculate the residuals and correction vectors
   !
   subroutine davidson_residuals
     !
     integer(ik) :: root
     integer(ik) :: corr, slot_corr
     real(rk)    :: corr_norm
     integer(ik) :: psi1, psi2         ! Slots for wavefunctions fields
     integer(ik) :: tmp                ! Slot for scratch
     complex(rk) :: hexp
     !
     call TimerStart('Davidson residuals')
     !
     psi1 = scr(1)
     psi2 = scr(2)
     tmp  = scr(3)
     !
     ! Calculate residuals for the roots we are working with.
     ! Because we already have the eigenvectors, we'll simply 
     ! apply our Hamiltonian to them.
     !
     corr = 0
     roots: do root=1,n_roots
       !
       ! Get the basis vector
       !
       call FieldIO('READ',psi1,slot=ind_new(root))
       !
       ! Calculate H |psi> for this vector
       !
       call QMHpsi(mass,f_pot,psi=psi1,Hpsi=psi2)
       !
       ! Print H expectation value, compared to the subspace result
       !
       if (verbose>=1) then
         hexp = FieldConjgIntegrate(psi1,psi2)
         write (out,"(t10,'<',i4,'|H|',i4,'> = ',2f14.8,' (subspace = ',f14.8,')')") &
                root, root, hexp, e(root)
       end if
       !
       ! ... and now the residual
       !
       call FieldAXPY(cmplx(-e(root),0.0_rk,kind=rk),src=psi1,dst=psi2)
       !
       ! We have the residual, remember it's norm
       !
       res_norm(root) = FieldNorm(psi2)
       if (verbose>=2) then
         write (out,"(t5,'Residual norm ',i4,' = ',f14.8)") &
                root, res_norm(root)
       end if
       if (res_norm(root)<=0) then
         !
         ! Exact solution, nothing to be done
         !
         if (verbose>=1) then
           write (out,"('Got an exact eigenvalue for root ',i4)") root
         end if
         cycle roots 
       end if
       corr = corr + 1
       !
       ! Calculate the correction vector, by approximately inverting the
       ! residual
       !
       call invert_residual(res=psi2,scr=tmp,inv=psi1,en=e(root), &
                            stencil=use_stencil(root))
       !
       ! Correction vector is now in psi1, make sure it's normalized
       !
       call QMNormalize(psi1,1.0_rk,corr_norm)
       if (verbose>=3) then
         write (out,"(t8,'correction norm = ',f14.8)") corr_norm
       end if
       !
       ! Save the correction vector
       !
       slot_corr = get_disk_slot()
       ind_corr(corr) = slot_corr
       call FieldIO('WRITE',psi1,slot=slot_corr)
     end do roots
     n_corr = corr
     call TimerStop('Davidson residuals')
   end subroutine davidson_residuals
   !
   ! Calculate laplacian error diffusion preconditioner for the davidson-liu step. 
   !
   subroutine invert_residual(res,scr,inv,en,stencil)
     integer(ik), intent(in) :: res     ! Input: Residual to be inverted
                                        ! Output: destroyed
     integer(ik), intent(in) :: scr     ! Scratch field
     integer(ik), intent(in) :: inv     ! Output: Approximate inverse of the residual
     real(rk), intent(in)    :: en      ! Estimated eigenvalue
     logical, intent(in)     :: stencil ! Use error diffusion stencil
     !
     call TimerStart('Davidson preconditioner')
     !
     call FieldCopy(f_pre,scr)
     !
     ! Calculate approximate local inverse of the Hamiltonian
     !
     call FieldInvertHamiltonian(scr,res,en,mass)
     !
     ! Apply Laplacian error-diffusion stencil. Because stencilling 
     ! restricts variational freedom, final iterations must be done
     ! with stencilling turned off.
     !
     if (stencil) then
       call FieldDiffuseLaplacian(res,inv)
     else
       call FieldCopy(res,inv)
     end if
     !
     ! That's it
     !
     call TimerStop('Davidson preconditioner')
   end subroutine invert_residual
   !
   ! Allocate some memory
   !
   subroutine allocate_memory
     integer(ik) :: alloc
     !
     allocate (ind_bas(max_basis+n_roots),ind_corr(max_basis), &
               ind_new(max_basis),disk_used(3*max_basis+n_roots),ovl(max_basis+n_roots, &
               max_basis+n_roots),h(max_basis,max_basis),vec(max_basis,max_basis), &
               e(max_basis),old_e(max_basis),diff_e(max_basis),use_stencil(max_basis), &
               res_norm(max_basis),stat=alloc)
     if (alloc/=0) then
       write (out,"('liu%release_memory - error ',i8)") alloc
       stop 'liu%allocate_memory'
     end if
     disk_used = .false.
     if (verbose>=1) then
       write (out,"(' Number of disk slots      = ',i4/)") 3*max_basis+n_roots
     end if
     call FieldIO('OPEN')
   end subroutine allocate_memory
   !
   ! Release memory we used
   !
   subroutine release_memory
     integer(ik) :: alloc
     !
     deallocate (ind_bas,ind_corr,ind_new,disk_used,ovl,h,vec,e,res_norm, &
                 old_e,diff_e,use_stencil,stat=alloc)
     if (alloc/=0) then
       write (out,"('liu%release_memory - error ',i8)") alloc
       stop 'liu%release_memory'
     end if
     call FieldIO('CLOSE')
   end subroutine release_memory
   !
   ! Reset disk slots status to a known set of values
   !
   subroutine reset_disk(indices)
     integer(ik), intent(in) :: indices(:)
     !
     disk_used          = .false.
     disk_used(indices) = .true.
   end subroutine reset_disk
   !
   ! Release a disk slot
   !
   subroutine free_disk_slots(indices)
     integer(ik), intent(in) :: indices(:)
     !
     disk_used(indices) = .false.
   end subroutine free_disk_slots
   !
   ! Allocate an empty scratch space slot
   !
   function get_disk_slot() result(p)
     integer(ik) :: p
    
     scan: do p=1,size(disk_used)
       if (disk_used(p)) cycle scan
       disk_used(p) = .true.
       return
     end do scan
     write (out,"('Internal error - Run out of disk slots')") 
     stop 'liu%get_disk_slot - no slots'
   end function get_disk_slot
   !
 end module liu
