    !
    !  We now need a generalized inverse of the overlap matrix. 
    !  We could have used lapack_ginverse, but we'd like a bit
    !  more control - so we do our own here.
    !
!   subroutine invert_smat_r(nmo_null,smat,smhalf,sphalf,use_block_diag,eps_smat)
!     integer(ik), intent(out)       :: nmo_null        ! Number of null-space vectors
!     real(rk), intent(inout)        :: smat(:,:)       ! In: Raw overlap matrix
!                                                       ! Out: Overlap matrix with null-space vectors removed
!     real(rk), intent(out)          :: smhalf(:,:)     ! S**(-0.5), with null-space removed
!     real(rk), intent(out)          :: sphalf(:,:)     ! S**(+0.5), with null-space removed
!     logical, intent(in), optional  :: use_block_diag  ! Default is .true.
!     real(rk), intent(in), optional :: eps_smat        ! Desired null-space cutoff. Default is 0 (near-machice accuracy)
      !
      !  Keep the routine type-generic below this line
      !
      integer(ik)                   :: nbas
      logical                       :: block
      real(kind(smat))              :: eps_
      integer(ik)                   :: alloc, iev, nneg
      real(kind(smat)), allocatable :: eval(:), evec(:,:), evs(:,:), ecnd(:), einv(:), ehalf(:)
      real(kind(smat))              :: eps
      !
      call TimerStart('Invert AO overlap')
      !
      !  Deal with optional parameters
      !
      nbas  = size(smat,dim=1)
      block = .true.
      if (present(use_block_diag)) block = use_block_diag
      eps_  = 0
      if (present(eps_smat)) eps_ = eps_smat
      !
      allocate (eval(nbas),ecnd(nbas),einv(nbas),ehalf(nbas), &
                evec(nbas,nbas),evs(nbas,nbas),stat=alloc)
      if (alloc/=0) then
        write (out,"('Allocation failed in scf_tools%invert_smat. Error = ',i0)") alloc
        call flush(out)
        stop 'scf_tools%invert_smat - no memory'
      end if
      evec = smat
      if (block) then
        call block_syev(evec,eval)
      else
        call lapack_syev(evec,eval)
      end if
      !
      if (verbose>0) then
        write (out,"(/t5,'Overlap eigenvalues before removal of linear dependencies:'/)")
        write (out,"((10(1x,g14.7)))") eval
        write (out,"()")
      end if
      !
      eps = real(1e3_xrk,kind(smat))*spacing(maxval(eval))
      if (eps_>0) then
        if (eps_<eps) then
          write (out,"('WARNING: eps_smat specified on input (',g13.5,') is very tight')") eps_
          write (out,"('WARNING: smallest value sensible for this basis is about ',g13.5)") eps
          write (out,"('WARNING: expect numerical problems')")
        end if
        eps = eps_
      end if
      !
      nmo_null = count(eval(:)<=eps)
      nneg     = count(eval(:)<-eps)
      !
      scan_eigenvalues: do iev=1,nbas
        if (eval(iev)>eps) then
          einv (iev) = eval(iev)**real(-0.25_xrk,kind(smat))  ! We need S**-0.5
          ehalf(iev) = eval(iev)**real(+0.25_xrk,kind(smat))  ! ... and S**+0.5
          ecnd (iev) = sqrt(eval(iev))                        ! We'll project out nullspace explicitly
        else
          einv (iev) = 0
          ehalf(iev) = 0
          ecnd (iev) = 0
        end if
      end do scan_eigenvalues
      !
      write (out,"(/'Overlap matrix null-space threshold is ',g13.5)") eps
      write (out,"(1x,i0,' eigenvectors of overlap matrix are in the nullspace')") nmo_null
      write (out,"('This number includes ',i0,' negative overlap eigenvalues')") nneg
      write (out,"('The smallest AO overlap eigenvalue was ',g13.5)") minval(eval)
      write (out,"(' The largest AO overlap eigenvalue was ',g13.5)") maxval(eval)
      write (out,"()")
      !
      !  Condition the overlap matrix
      !
      evs    = evec * spread(ecnd,dim=1,ncopies=size(evec,dim=1))
      smhalf = mt_matmul(evs,transpose(evs))
      write (out,"('Maximum change in AO overlap upon null-space projection was ',g13.5)") &
             maxval(abs(smhalf-smat))
      smat   = smhalf
      ! 
      !  Construct S**(-1/2)
      !
      evs    = evec * spread(einv,dim=1,ncopies=size(evec,dim=1))
      smhalf = mt_matmul(evs,transpose(evs))
      ! 
      !  Construct S**(+1/2)
      !
      evs    = evec * spread(ehalf,dim=1,ncopies=size(evec,dim=1))
      sphalf = mt_matmul(evs,transpose(evs))
      !
      deallocate (eval,ecnd,einv,ehalf,evec,evs)
      call TimerStop('Invert AO overlap')
!   end subroutine invert_smat_r
