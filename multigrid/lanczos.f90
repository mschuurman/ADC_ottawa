 module lanczos
!
!  Iterative diagonalization, using an extremely naive implementation of
!  the Lanczos algorithm, with no subspace compaction and explicit orthogonalization.
!
!  This implementation is intended to be robust in the presence of numerical
!  noise, but is slow and I/O intensive. It will take an almost infinite
!  amount of time to find excited-state solutions; to get those, one would
!  have to get the ground state first, then shift it up in the eigenvalue 
!  spectrum.
!
   use accuracy
   use lapack
   use multigrid
   use qmech
   implicit none
   private
   public LCeigenvalues, LCeigenvector
!
! Interesting data entries
!
   integer(ik), parameter :: verbose = 1       ! Verbosity level
   integer(ik), parameter :: maxVectors = 24   ! Maximum number of basis vectors to use
   integer(ik), parameter :: maxCores   = 20   ! Maximum number of core states
   integer(ik)            :: nv                ! Current number of vectors
   integer(ik)            :: ncores            ! Actual number of core states
   real(rk)               :: mass              ! Mass of the QM particle
!
! Reduced eigenproblem - H matrix, eigenvalues, and eigenvectors
!
   complex(rk), save      :: lan_h   (maxVectors,maxVectors)
   real(rk), save         :: lan_eval(maxVectors) 
   complex(rk), save      :: lan_evec(maxVectors,maxVectors)
   integer(ik), save      :: f_scr   (2*maxVectors)
   integer(ik), save      :: f_cores (maxCores)
   logical, save          :: incore
!
! Fields used in this routine:
!
   integer(ik)            :: Hpot              ! One-electron potential
   integer(ik)            :: vec_cur           ! current Lanczos vector (v_{i})
   integer(ik)            :: vec_scr           ! space for various scratch vectors
   integer(ik)            :: Hvec              ! space for H v_{i}
!
   contains
!
! Driver routines, visible from outside
!
! LCeigenvalues calculates eigenvalues of a specified Hamiltonian. It leaves
! enough data lying around, to allow reconstruction of a desired eigenvector,
! with a later call to LCeigenvector
!
   subroutine LCeigenvalues(mass_,pot,guess,scr1,scr2,eval,cnv_e_,conv,f_scr_,f_cores_)
    real(rk), intent(in)    :: mass_                 ! Mass of the QM particle
    integer(ik), intent(in) :: pot                   ! Field containing multiplicative potential
                                                     ! The field is unchanged upon return
    integer(ik), intent(in) :: guess                 ! Field, containing the initial guess. 
                                                     ! It is destroyed upon return
    integer(ik), intent(in) :: scr1, scr2            ! Two scratch fields, for use inside LCeigenvalues
    real(rk), intent(out)   :: eval(:)               ! Eigenvalues. At least size(eval) eigenvalues will
                                                     ! be calculated and returned.
    real(rk), intent(in), optional :: cnv_e_         ! Desired convergence
    logical, intent(out), optional :: conv           ! Convergence flag
    integer(ik), intent(in), optional :: f_scr_(:)   ! In-core scratch fields
    integer(ik), intent(in), optional :: f_cores_(:) ! Core states - keep orthogonal to those

    integer(ik) :: neval
    real(rk)    :: change, eps, cnv_e
    logical     :: converged

    !
    !  Initialization
    !
    mass     = mass_
    neval    = size(eval)
    Hpot     = pot
    vec_cur  = guess
    vec_scr  = scr1
    Hvec     = scr2
    lan_h    = 0
    nv       = 1
    cnv_e    = 1e-7
    if (present(cnv_e_)) cnv_e = cnv_e_
    ncores   = 0
    if (present(f_cores_)) then
      ncores = size(f_cores_)
      f_cores(:ncores) = f_cores_
    end if
    !
    !  Decide whether we'll use in-core vector storage or I/O
    !
    incore = .false.
    if (present(f_scr_)) then
      if (size(f_scr_)>=size(f_scr)) then
        incore = .true.
        f_scr  = f_scr_(1:size(f_scr))
        write (out,"('lanczos: Using in-core vector storage')")
      end if
    end if
    !
    if (.not.incore) call FieldIO('OPEN')
    !
    !  Lanczos/Davidson iteration
    !
    converged = .false.
    do while (.not.converged .and. (nv<maxVectors))
      call diag_iteration
      call diag_eigensystem
      if (nv>max(neval,2)+1) then
        change = maxval(abs(lan_eval(1:neval)-eval(:)))
        eps    = max( spacing(maxval(lan_eval(1:nv-1))), cnv_e )
        if (verbose>=1) then
          write (out,"(' On iteration ',i4,' min E = ',f20.12,' (',f20.12,') max. chg = ',f20.12,' eps = ',f20.12)") &
            nv, lan_eval(1), lan_eval(2), change, eps
        end if
        if (change<eps) converged = .true.
      else
        if (verbose>=1) then
          write (out,"(' On iteration ',i4,' min E = ',f20.12)") nv, lan_eval(1)
        end if
      end if
      eval(:) = lan_eval(1:neval)
    end do
    if (.not.converged) write (out,"(8x,'*WARNING: Lanczos iterations not converged')")
    if (present(conv)) conv = converged
  end subroutine LCeigenvalues
!
! LCeigenvector can be used to obtain the final solution eigenvector, for
! one of the roots found in LCeigenvalues. It will also clean up scratch
! files, left around by LCeigenvalues.
!
  subroutine LCeigenvector(psi,scr,a_root,a_preserve)
    integer(ik), intent(in)           :: psi        ! Out: Optimized eigenfunction. Specifying
                                                    ! zero for psi will suppress calculation of
                                                    ! the wavefunction, and simply clean up the
                                                    ! scratch.
    integer(ik), intent(in)           :: scr        ! Scatch area, destroyed upon exit
    integer(ik), intent(in), optional :: a_root     ! Root to get wavefunction for; if not
                                                    ! specified, uses root 1 (ground state).
    logical, intent(in), optional     :: a_preserve ! Keep scratch files - another wavefunction
                                                    ! will be requested later.

    integer(ik) :: root, iv
    logical     :: preserve
    !
    !  Deal with the optional arguments
    !
    root = 1 
    if (present(a_root)) root = a_root
    preserve = .false.
    if (present(a_preserve)) preserve = a_preserve
    !
    !  Calculate the wavefunction. The coefficients are already in lan_evec; 
    !  all we need is to go over all trial wavefunctions, and add them up
    !
    if (psi>0) then
      if (verbose>=1) then
        write (out,"(' Calculating wavefunction for root ',i3,' (E=',f12.7,')')") &
               root, lan_eval(root)
        write (out,"(' Coefficients of the trial functons : ')")
        write (out,"((2x,4(1x,'(',f10.5,',',f10.5,')')))") lan_evec(1:nv-1,root)
      end if
      call FieldZero(psi)
      do iv=1,nv-1
        if (incore) then
          call FieldCopy(dst=scr,src=f_scr(2*iv-1))
        else
          call FieldIO('READ',scr,2*iv-1)
        end if
        call FieldAXPY(lan_evec(iv,root),scr,psi)
      end do
      call QMProject(f_cores(:ncores),psi)  ! Apply the projector
    end if
    !
    !  Clean up and exit
    !
    if (.not.preserve) then
      if (.not.incore) call FieldIO('CLOSE')
    end if
  end subroutine LCeigenvector
!
! Internal stuff
!
! diag_eigensystem - calculates eigenvalues and eigenvectors of the 
!      reduced space, in terms of the basis functions.
!
  subroutine diag_eigensystem
    integer(ik) :: nbas

    nbas = nv - 1
    !
    !  We have both the upper and lower diagonal elements of H; average them to
    !  reduce numerical noise.
    !
    lan_evec(1:nbas,1:nbas) = 0.5_rk * ( lan_h(1:nbas,1:nbas) + transpose(conjg(lan_h(1:nbas,1:nbas))) )
    call lapack_heev(lan_evec(1:nbas,1:nbas),lan_eval(1:nbas))

    if (verbose>=2) then
      write (out,"(' On iteration ',i4,' smallest reduced eigenvalue is ',f18.10)") nv, lan_eval(1)
    end if
    if (verbose>=2 .and. nbas>1) then
      write (out,"(' Remaining eigenvalues are: ')")
      write (out,"((5x,6(1x,f15.7)))") lan_eval(2:nbas)
    end if
  end subroutine diag_eigensystem
!
! diag_iteration - performs basic step of Lanczos/Davidson
!
  subroutine diag_iteration
    integer(ik) :: iv
    real(rk)    :: norm
    !
    !  0. Project and normalize current search vector
    !
    call QMProject(f_cores(:ncores),vec_cur)
    norm = FieldNorm(vec_cur)
    call FieldScale(vec_cur,cmplx(1.0_rk/norm,0.0_rk,kind=rk))
    if (verbose>=2) then
      write (out,"(' On iteration ',i4,' original guess norm = ',f15.11)") nv, norm
    end if
    !
    !  1. Calculate sub-diagonal entries of the reduced-space Hamiltonian.
    !     To do this, read all the previously calculated residuals (even
    !     records of the external data file), and integrate them with the 
    !     current basis vector. 
    !
    do iv=1,nv-1
      if (incore) then
        call FieldCopy(dst=vec_scr,src=f_scr(2*iv))
      else
        call FieldIO('READ',vec_scr,2*iv)
      end if
      lan_h(nv,iv) = FieldConjgIntegrate(vec_scr,vec_cur)
    end do
    !
    !  2. Apply the Hamiltonian to the current search vector, and calculate
    !     the next diagonal entry of the reduced-space Hamiltonian.
    !
    call QMHpsi(mass,Hpot,vec_cur,Hvec)
    call QMProject(f_cores(:ncores),Hvec)
    lan_h(nv,nv) = QMExpectation(vec_cur,Hvec)
    if (verbose>=2) then
      write (out,"(' Energy expectation value of the guess ',2f15.7)") lan_h(nv,nv)
    end if
    !
    !  3. Calculate the residual for the current search vector, and store
    !     basis vector and the residual for later use.
    !
    call FieldAXPY(-lan_h(nv,nv),vec_cur,Hvec)
    if (incore) then
      call FieldCopy(src=vec_cur,dst=f_scr(2*nv-1))
      call FieldCopy(src=Hvec,   dst=f_scr(2*nv)  )
    else
      call FieldIO('WRITE',vec_cur,2*nv-1)
      call FieldIO('WRITE',Hvec,   2*nv  )
    end if
    !
    !  4. Calculate superdiagonal elements of the reduced Hamiltonian, and
    !     use them to orthogonalize the new search vector wrt all previous
    !     search directions (odd records).
    !
    call FieldCopy(Hvec,vec_cur)
    do iv=1,nv-1
      if (incore) then
        call FieldCopy(dst=vec_scr,src=f_scr(2*iv-1))
      else
        call FieldIO('READ',vec_scr,2*iv-1)
      end if
      lan_h(iv,nv) = FieldConjgIntegrate(vec_scr,Hvec)
      call FieldAXPY(-lan_h(iv,nv),vec_scr,vec_cur)
    end do
    if (verbose>=3) then
      do iv=1,nv-1
        write (out,"(' H(',i3,',',i3,') = ',2f12.7,' H(',i3,',',i3,') = ',2f12.7,' delta = ',2f12.7)") &
               iv, nv, lan_h(iv,nv), nv, iv, lan_h(nv,iv), lan_h(iv,nv) - conjg(lan_h(nv,iv))
      end do
    end if
    !
    !  5. Advance to the next iteration.
    !
    nv = nv + 1
  end subroutine diag_iteration

 end module lanczos
