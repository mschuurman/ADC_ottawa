 module lanczos_complex
!
!  This module is derived from lanczos.f90
!
!  Iterative diagonalization, using an extremely naive implementation of
!  the Lanczos algorithm, with no subspace compaction and explicit 
!  orthogonalization. Hamiltonian matrix is allowed to become non-Hermitian.
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
   public LCeigenvaluesComplex, LCeigenvectorComplex
   !
   !  Interesting data entries
   !
   integer(ik), parameter :: verbose    = 10   ! Verbosity level
   integer(ik), parameter :: maxVectors = 8    ! Maximum number of basis vectors to use
   integer(ik), parameter :: maxFields  = 100  ! Maximum number of scratch fields to use
   integer(ik), save      :: nv                ! Current number of vectors
   !
   !  Reduced eigenproblem - H matrix, eigenvalues, and eigenvectors
   !
   complex(rk), save      :: lan_h   (maxVectors,maxVectors)
   complex(rk), save      :: lan_eval(maxVectors) 
   complex(rk), save      :: lan_evec(maxVectors,maxVectors)
   integer(ik), save      :: root
   !
   !  Fields used inside LCeigenvaluesComplex
   !
   integer(ik), save      :: fields(maxFields)  ! List of data fields available
   integer(ik), save      :: n_fields           ! Number of data fields available
   !
   contains
   !
   !  Driver routines, visible from outside
   !
   !  LCeigenvalues calculates eigenvalues of a specified Hamiltonian. It leaves
   !  enough data lying around, to allow reconstruction of a desired eigenvector,
   !  with a later call to LCeigenvector
   !
   !  Subroutine ham must take three arguments, defined as:
   !
   !    integer(ik) :: input     = input wavefunction
   !    integer(ik) :: hermit    = Terms which must be Hermitian resulting H|psi> product
   !    integer(ik) :: nothermit = Possibly non-Hermitian terms
   !
   subroutine LCeigenvaluesComplex(ham,eval,flds,cnv_e_,conv)
    external                       :: ham     ! Routine evaluating the Hamiltonian
                                              ! for a given right-hand side.
    complex(rk), intent(out)       :: eval    ! Eigenvalue
    integer(ik), intent(in)        :: flds(:) ! List of vectors, which can be used for
                                              ! scratch. The first vector should containt
                                              ! the initial wavefunction guess.
    real(rk), intent(in), optional :: cnv_e_  ! Desired convergence
    logical, intent(out), optional :: conv    ! Convergence flag
    !
    real(rk)    :: change, eps, cnv_e
    logical     :: converged
    !
    !  Initialization
    !
    n_fields           = min(maxFields,size(flds))
    fields(1:n_fields) = flds(1:n_fields)
    lan_h              = 0
    nv                 = 1
    cnv_e              = 1e-7
    if (present(cnv_e_)) cnv_e = cnv_e_
    call FieldIO('OPEN')
    !
    !  Lanczos/Davidson iteration
    !
    converged = .false.
    do while (.not.converged .and. (nv<maxVectors))
      call diag_iteration(ham)
      call diag_eigensystem
      if (nv>3) then
        change = abs(lan_eval(root)-eval)
        eps    = max( spacing(maxval(abs(lan_eval(1:nv-1)))), cnv_e )
        if (verbose>=1) then
          write (out,"(' On iteration ',i4,' min E = ',2g13.6,' chg = ',g13.6,' eps = ',g13.6)") &
            nv, lan_eval(root), change, eps
        end if
        if (change<eps) converged = .true.
      else
        if (verbose>=1) then
          write (out,"(' On iteration ',i4,' min E = ',2g13.6)") nv, lan_eval(root)
        end if
      end if
      eval = lan_eval(root)
    end do
    if (.not.converged) write (out,"(8x,'*WARNING: Lanczos iterations not converged')")
    if (present(conv)) conv = converged
  end subroutine LCeigenvaluesComplex
  !
  !  LCeigenvector can be used to obtain the final right solution eigenvector, 
  !  for the lowest root found in LCeigenvalues. It will also clean up scratch
  !  files, left around by LCeigenvalues.
  !
  subroutine LCeigenvectorComplex(psi,scr)
    integer(ik), intent(in)           :: psi        ! Out: Optimized eigenfunction. Specifying
                                                    ! zero for psi will suppress calculation of
                                                    ! the wavefunction, and simply clean up the
                                                    ! scratch.
    integer(ik), intent(in)           :: scr        ! Scratch area, destroyed upon exit
    !
    integer(ik) :: iv
    !
    !  Calculate the wavefunction. The coefficients are already in lan_evec; 
    !  all we need is to go over all trial wavefunctions, and add them up
    !
    if (psi>0) then
      if (verbose>=1) then
        write (out,"(' Coefficients of the trial functons : ')")
        write (out,"((2x,4(1x,'(',f10.5,',',f10.5,')')))") lan_evec(1:nv-1,root)
      end if
      call FieldZero(psi)
      do iv=1,nv-1
        call FieldIO('READ',scr,3*iv-2)
        call FieldAXPY(lan_evec(iv,root),scr,psi)
      end do
    end if
    !
    !  Clean up and exit
    !
    call FieldIO('CLOSE')
  end subroutine LCeigenvectorComplex
  !
  ! Internal stuff
  !
  ! diag_eigensystem - calculates eigenvalues and eigenvectors of the 
  !                    reduced space, in terms of the basis functions.
  !
  subroutine diag_eigensystem
    integer(ik) :: nbas, ir
    !
    nbas = nv - 1
    !
    if (verbose>=5) then
      write (out,"(/t3,'Reduced Hamiltonian is:'/)")
      ham_rows: do ir=1,nbas
        write (out,"(t3,8(2g11.4,2x))") lan_h(ir,1:nbas)
      end do ham_rows
      write (out,"()")
    end if
    !
    lan_evec(1:nbas,1:nbas) = lan_h(1:nbas,1:nbas)
    call lapack_geev(lan_evec(1:nbas,1:nbas),lan_eval(1:nbas))
    !
    !  GEEV routines return roots in no particular order; find the
    !  root with smallest real part
    !
    root = minloc(real(lan_eval(1:nbas),kind=rk),dim=1)
    !
    if (verbose>=5) then
      write (out,"(/t3,'Lowest solution is ',i3,' at ',2g14.7)") root, lan_eval(root)
      write (out,"(/t3,'Reduced Hamiltonian solutions are:'/)")
      write (out,"(t3,8(2g11.4,2x))") lan_eval(1:nbas)
      write (out,"()")
      vec_rows: do ir=1,nbas
        write (out,"(t3,8(2g11.4,2x))") lan_evec(ir,1:nbas)
      end do vec_rows
      write (out,"()")
    end if
    !
  end subroutine diag_eigensystem
  !
  !  diag_iteration - performs basic step of Lanczos/Davidson
  !
  subroutine diag_iteration(ham)
    external    :: ham
    integer(ik) :: iv
    real(rk)    :: norm
    integer(ik) :: vec_cur, vec_herm, vec_noth, vec_scr
    complex(rk) :: lan_herm(nv,nv) ! Hermitian part of the Hamiltonian - to be symmetrized
    complex(rk) :: ave, ort
    !
    vec_cur  = fields(1)
    vec_herm = fields(2)
    vec_noth = fields(3)
    vec_scr  = fields(4)
    if (n_fields<4) stop 'lanczos_complex - not enough scratch fields!'
    !
    !  0. Normalize current search vector
    !
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
    lan_herm = 0
    do iv=1,nv-1
      call FieldIO('READ',vec_herm,3*iv-1)
      call FieldIO('READ',vec_noth,3*iv-0)
      lan_herm(nv,iv) = FieldConjgIntegrate(vec_herm,vec_cur)
      lan_h   (nv,iv) = FieldConjgIntegrate(vec_noth,vec_cur)
    end do
    !
    !  2. Apply the Hamiltonian to the current search vector, and calculate
    !     the next diagonal entry of the reduced-space Hamiltonian.
    !
    call ham(vec_cur,vec_herm,vec_noth)
    lan_h(nv,nv) = real(FieldConjgIntegrate(vec_herm,vec_cur),rk) &
                      + FieldConjgIntegrate(vec_noth,vec_cur)
    if (verbose>=2) then
      write (out,"(' Energy expectation value of the guess ',2g15.7)") lan_h(nv,nv)
    end if
    !
    !  3. Calculate the residual for the current search vector, and store
    !     basis vector and the residual for later use.
    !
    call FieldAXPY(alpha=-lan_h(nv,nv),src=vec_cur,dst=vec_noth)
    call FieldIO('WRITE',vec_cur, 3*nv-2)
    call FieldIO('WRITE',vec_herm,3*nv-1)
    call FieldIO('WRITE',vec_noth,3*nv-0)
    !
    !  4. Calculate superdiagonal elements of the reduced Hamiltonian, and
    !     use them to orthogonalize the new search vector wrt all previous
    !     search directions (aka guess wavefunctions).
    !
    call FieldCopy(vec_herm,vec_cur)
    call FieldAXPY((1.0_rk,0.0_rk),vec_noth,vec_cur)
    do iv=1,nv-1
      call FieldIO('READ',vec_scr,3*iv-2)
      lan_herm(iv,nv) = FieldConjgIntegrate(vec_scr,vec_herm)
      lan_h   (iv,nv) = FieldConjgIntegrate(vec_scr,vec_noth)
      !
      ort = lan_herm(iv,nv) + lan_h   (iv,nv)
      ave = 0.5_rk*(lan_herm(iv,nv)+conjg(lan_herm(nv,iv)))
      lan_h(iv,nv) = lan_h(iv,nv) + ave
      lan_h(nv,iv) = lan_h(nv,iv) + conjg(ave)
      !
      call FieldAXPY(-ort,vec_scr,vec_cur)
    end do
    !
    !  5. Advance to the next iteration.
    !
    nv = nv + 1
  end subroutine diag_iteration

 end module lanczos_complex
