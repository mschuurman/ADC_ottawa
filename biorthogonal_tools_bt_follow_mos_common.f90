!   subroutine bt_follow_mos_rk(nmo_act,mo_occ,eps_eval,sphalf,mosg,mo_energy,mos)
!     integer(ik), intent(in)    :: nmo_act       ! Number of orbitals to include in the matching process
!     real(rk), intent(in)       :: eps_eval      ! Degeneracy threshold for eigenvalues. Zero or negative means machine accuracy
!     real(rk), intent(in)       :: sphalf(:,:)   ! Overlap matrix to the power +1/2
!     complex(rk), intent(in)    :: mosg(:,:,:)   ! Reference MOs; the last index is 1 of the left vectors; 2 for the right vectors
!     real(rk), intent(in)       :: mo_occ(:)     ! Orbital occupations, only used for printing
!     complex(rk), intent(inout) :: mo_energy(:)  ! In/Out: Orbital energies, same order as mos()
!     complex(rk), intent(inout) :: mos(:,:,:)    ! In: MOs to be aligned
!                                                 ! Out: MOs where degenerate subblocks are in maximal conincidence with the guess
      !
      !  Below this line we should be type-generic
      !
      integer(ik)                     :: alloc
      integer(ik)                     :: is, ie, iev3
      integer(ik)                     :: bs            ! Block size
      integer(ik)                     :: nmo           ! Total number of MOs
      integer(ik)                     :: nmo_null      ! Number of null-space MOs
      real(kind(mos))                 :: eps
      complex(kind(mos)), allocatable :: tmat(:,:)     ! Overlap of the current left MOs with reference right MOs.
      complex(kind(mos)), allocatable :: umat(:,:)     ! Overlap of the reference left MOs with current right MOs.
      logical, allocatable            :: ref_taken(:)  ! True if a particular reference MO has been claimed already
      integer(ik), allocatable        :: cur2ref(:)    ! For each current MO, gives the index of the reference MO
      complex(kind(mos)), allocatable :: ref_grid (:,:)! Projections from the reference space to orbital space
      complex(kind(mos)), allocatable :: ref_score(:)  ! Population score for each reference MO within the current block
      integer(ik), allocatable        :: ref_order(:)  ! Score-based reference MO ordering
      complex(kind(mos)), allocatable :: ref_score2(:) ! Additional scoring and ordering array when best-projection fails
      integer(ik), allocatable        :: ref_order2(:) ! 
      real(kind(mos)), allocatable    :: abs_score (:) ! -abs(ref_score) or -abs(ref_score2), for sorting
      complex(kind(mos)), allocatable :: vcur(:,:)     ! (Unitary) rotation matrix within the degenerate sub-block
      !
      call TimerStart('Follow reference MOs')
      !
      nmo      = size(mos,dim=2)
      nmo_null = nmo - nmo_act
      !
      !  Null-space MOs are supposed to be at the end of the list of the MOs at this
      !  point. We have a choice of either using all MOs, including the null-space
      !  in the process of similarity matching, or not.
      !
      allocate (tmat(nmo_act,nmo_act),umat(nmo_act,nmo_act), ref_grid(nmo_act,nmo_act), &
                ref_taken(nmo_act), cur2ref(nmo_act), ref_score(nmo_act), ref_order(nmo_act), &
                ref_score2(nmo_act), ref_order2(nmo_act), abs_score(nmo_act), stat=alloc)
      if (alloc/=0) then
        write (out,"('biorthogonal_tools%bt_follow_mos: Error ',i0,' allocating memory')") alloc
        call stop('biorthogonal_tools%bt_follow_mos - no memory (1)')
      end if
      !
      !  Calculate overlaps. An obvious way to perform the transformation
      !  is to do:
      !
      !    tmat = matmul(transpose(mos (:,:,1)),matmul(smat,mosg(:,:,2)))
      !    umat = matmul(transpose(mosg(:,:,1)),matmul(smat,mos (:,:,2)))
      !
      !  However, this transformation is prone to losing accuracy. Our
      !  way below is slower, but more accurate - especially for ill-conditioned
      !  overlap matrices (aren't they all ill-conditioned?)
      !
      tmat     = mt_matmul(mt_matmul(transpose(mos (:,:nmo_act,1)),sphalf),mt_matmul(sphalf,mosg(:,:nmo_act,2)))
      umat     = mt_matmul(mt_matmul(transpose(mosg(:,:nmo_act,1)),sphalf),mt_matmul(sphalf,mos (:,:nmo_act,2)))
      ref_grid = umat(:,:)*transpose(tmat(:,:))
      !
      !  Scan for degeneracies is helfully lifted from biorthogonalize_eigenvectors
      !
      ref_taken = .false.
      eps       = eigenvalue_degeneracy_epsilon(eps_eval,mo_energy,'follow_mos')
      is        = 1
      scan_for_degeneracies: do while (is<=nmo_act)
        iev3 = is
        find_identical_eigenvalues: do while(iev3<=nmo_act)
          if (abs(mo_energy(iev3)-mo_energy(is))>eps) exit find_identical_eigenvalues
          ie   = iev3
          iev3 = iev3 + 1
        end do find_identical_eigenvalues
        !
        !  We found a degenerate block (possibly of size 1); now do something about it
        !
        bs = ie - is + 1
        if (pick_best_references()) then
          !
          !  When pick_best_references succeeds, it places the list of reference MOs
          !  in ref_order(:bs). If it fails, we most likely have linearly-dependent
          !  MOs, for which tracking is meaningless. 
          !
          allocate (vcur(bs,bs),stat=alloc)
          if (alloc/=0) then
            call stop('biorthogonal_tools%bt_follow_mos - no memory (2)')
          end if
          call find_optimal_rotation(tmat(is:ie,ref_order(:bs)),umat(ref_order(:bs),is:ie),vcur)
          mos(:,is:ie,1) = mt_matmul(mos(:,is:ie,1),conjg(vcur))
          mos(:,is:ie,2) = mt_matmul(mos(:,is:ie,2),      vcur )
          deallocate (vcur)
        end if
        ! 
        is = ie + 1
      end do scan_for_degeneracies
      !
      !  Rearrange eigenvectors and eigenvalues to match order of the reference orbitals
      !  We are leaving out the null-space MOs
      !
      mos      (:,cur2ref(:),:) = mos(:,:nmo_act,:)
      mo_energy(  cur2ref(:)  ) = mo_energy(:nmo_act)
      !
      if (verbose>0) then
        !
        !  Report on the final orbital overlaps
        !
        tmat = mt_matmul(mt_matmul(transpose(mos (:,:nmo_act,1)),sphalf),mt_matmul(sphalf,mosg(:,:nmo_act,2)))
        umat = mt_matmul(mt_matmul(transpose(mosg(:,:nmo_act,1)),sphalf),mt_matmul(sphalf,mos (:,:nmo_act,2)))
        write (out,"((2x,a5,2x,a8,2x,a16,1x,a16,2(2x,a13,1x,a13)))") &
        ' MO ', ' Occ ', ' Re(eps)    ', ' Im(eps)    ', ' Re(<RefL|R>) ', ' Im(<RefL|R>) ', ' Re(<L|RefR>) ', ' Im(<L|RefR>) ', &
        '----', '-----', '---------   ', '---------   ', '--------------', '--------------', '--------------', '--------------'
        report_states: do is=1,nmo_act
          write (out,"(2x,i5,2x,f8.6,2x,g16.8,1x,g16.8,2(2x,f13.8,1x,f13.8))") &
                 is, mo_occ(is), mo_energy(is), umat(is,is), tmat(is,is)
        end do report_states
        write (out,"()")
      end if
      !
      deallocate (tmat,umat,ref_grid,ref_taken,cur2ref,ref_score,ref_order,ref_score2,ref_order2,abs_score)
      call bt_sanity_check(mos(:,:,:),sphalf=sphalf,msg='follow_mos',null_count=nmo_null)
      call TimerStop('Follow reference MOs')
      !
      contains
      !
      function pick_best_references() result(not_fail)
        complex(kind(mos)) :: cov_test(bs)
        logical            :: zeros   (bs), fail, not_fail
        integer(ik)        :: set_zero(bs)
        integer(ik)        :: set_not (bs)
        integer(ik)        :: cnt_zero, cnt_not, io
        !
        ref_score = sum(ref_grid(:,is:ie),dim=2)
        abs_score = -abs(ref_score)
        where (ref_taken) abs_score = 1 ! Place all orbitals already taken at the very back
        call order_keys(abs_score,ref_order)
        !
        !  There is a possible problem with this choice - we could possibly pick a set 
        !  of references which has essentially no overlap with some of the target 
        !  orbitals. If this happens, the only practical recourse is to separate out
        !  the orbitals which had no coverage, and treat them as a separate degenerate
        !  sub-block for selecting the references.
        !
        cov_test = sum(ref_grid(ref_order(:bs),is:ie),dim=1)
        zeros    = abs(cov_test)<0.1
        !
        if ((any(zeros) .and. verbose>=2) .or. verbose>=3) then
          write (out,"('Degenerate MOs ',i0,'-',i0)") is, ie
          write (out,"('     Energy range: ',4f16.8)") mo_energy((/is,ie/))
          write (out,"('Best overlap with: ',20i20)") ref_order(:bs)
          write (out,"('       Overlap is: ',20(2f10.6))") ref_score(ref_order(:bs))
          write (out,"('      Coverage is: ',20(2f10.6))") cov_test
        end if
        !
        fail = .false.
        if (any(zeros)) then
          !
          !  Build lists of zero and non-zero coverage, and look them up separately
          !
          cnt_zero = 0
          cnt_not  = 0
          build_sets: do io=is,ie
            if (zeros(io-is+1)) then
              cnt_zero = cnt_zero + 1 ; set_zero(cnt_zero) = io
            else
              cnt_not  = cnt_not  + 1 ; set_not (cnt_not ) = io
            end if
          end do build_sets
          !
          if (verbose>=2) then
            write (out,"('Set of non-zeros: ',20i5)") set_not (:cnt_not)
            write (out,"('    Set of zeros: ',20i5)") set_zero(:cnt_zero)
          endif
          !
          !  Look up best matches for orbitals which had coverage before
          !
          ref_score = sum(ref_grid(:,set_not(:cnt_not)),dim=2)
          abs_score = -abs(ref_score)
          where (ref_taken) abs_score = 1
          call order_keys(abs_score,ref_order)
          if (verbose>=2) then
            write (out,"('Best overlap of non-zeros with: ',20i20)") ref_order(:cnt_not)
            write (out,"('                    Overlap is: ',20(2f10.6))") ref_score(ref_order(:cnt_not))
          end if
          !
          !  Look up orbitals which had no coverage; make sure to exclude 
          !  the references picked up by the "covered" list.
          !
          ref_score2 = sum(ref_grid(:,set_zero(:cnt_zero)),dim=2)
          abs_score  = -abs(ref_score2)
          where (ref_taken) abs_score = 1
          abs_score(ref_order(:cnt_not)) = 1
          call order_keys(abs_score,ref_order2)
          if (verbose>=2) then
            write (out,"('Best overlap of zeros with: ',20i20)") ref_order2(:cnt_not)
            write (out,"('                Overlap is: ',20(2f10.6))") ref_score2(ref_order2(:cnt_not))
          end if
          !
          !  Merge the two lists, and repeat the coverage test. If we still
          !  have orbitals with zero coverage, we'll have to give up.
          !
          ref_order(cnt_not+1:cnt_not+cnt_zero) = ref_order2(:cnt_zero)
          !
          cov_test = sum(ref_grid(ref_order(:bs),is:ie),dim=1)
          fail     = any(abs(cov_test)<0.1)
          if (verbose>=2) then
            write (out,"('Revised references are: ',20i10)") ref_order(:bs)
            write (out,"('   Revised coverage is: ',20(2f10.6))") cov_test
            if (.not.fail) then
              write (out,"('Tracking succeded with difficulties.')")
              write (out,"('If this is a finite-field calculation, consider increasing the field in small(er) increments.'/)")
            end if
          end if
        end if
        !
        if (fail .and.verbose>=2) then
          write (out,"('Tracking coverage test failed. Expect meaningless results.')")
          write (out,"('If this is a finite-field calculation, try increasing the field in small(er) increments.'/)")
        end if
        if (any(ref_taken(ref_order(:bs)))) then
          write (out,"('biorthogonal_tools%bt_follow_mos - a reference orbital was picked twice.')") 
          call stop('biorthogonal_tools%bt_follow_mos - logic error in pick_best_references')
        end if
        ref_taken(ref_order(:bs)) = .true.
        cur2ref(is:ie) = ref_order(:bs)
        not_fail = .not.fail
      end function pick_best_references
!   end subroutine bt_follow_mos_rk
