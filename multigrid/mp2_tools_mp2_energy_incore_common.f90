! subroutine mp2_energy_incore_real(occ_det,occ_ref,eref,eps,int2e,mos,isa_de,eps_big,storage_mode)
!   integer(sik), intent(in)               :: occ_det(:)    ! Orbital occupation numbers for the MP2 reference determinant
!   real(rk), intent(in)                   :: occ_ref(:)    ! Orbital occupation numbers for the FON-SCF reference determinant
!   complex(rk), intent(in)                :: eref          ! Energy of the FON-SCF reference determinant
!   complex(rk), intent(in)                :: eps(:)        ! Orbital eigenvalues of the active orbitals
!   type(int2e_cache), intent(inout)       :: int2e         ! 2e AO integrals
!   complex(rk), intent(in)                :: mos(:,:,:)    ! Bi-orthogonal spin-MOs.
!   real(rk), intent(in)                   :: isa_de(:)     ! Intruder-state avoidance parameter(s) list
!   real(rk), intent(in), optional         :: eps_big       ! Threshold for picking out important orbitals; default is 1e-2
!   character(len=*), intent(in), optional :: storage_mode  ! Accuracy used for storing transformed integrals; can be
!                                                           ! 'real', 'quad', or 'as-is'
    !
    !  No explicit kind references below this point!
    !
    type(moint2e_cache)        :: m2e_vo             ! 2-electron integrals of the V,O,V,O type (LRLR)
    type(moint2e_cache)        :: m2e_ov             ! 2-electron integrals of the O,V,O,V type (LRLR)
                                                     ! We'll use a dirty trick for m2e_ov, and swap the L and R indices 
                                                     ! to make transformation more efficient.
    type(moint2e_cache)        :: m2e_dd             ! 2-electron integrals of the X,X,X,D type
    !                                                
    integer(ik)                :: nisa, isa          ! Number of ISA energy shifts and corresponding counter
    complex(kind(mos))         :: e_mp1              ! MP1 energy correction; this must be zero if there are no fractional occupations
    integer(ik)                :: norbs              ! Number of orbitals 
    integer(ik)                :: nocc               ! Number of occupied orbitals (aka electrons) in the MP2 reference determinant
    integer(ik)                :: nvir               ! Number of virtual orbitals in the MP2 referecne determinant
    integer(ik)                :: ndlt               ! Number of orbitals where occupation numbers differ in FON-SCF and MP2 references
    real(kind(mos)), parameter :: eps_occ = 1e-12_rk ! Treat occupation number as as identical if difference is smaller than this.
    integer(ik)                :: i, j               ! Orbital indices, actual MOs
    integer(ik)                :: pi, pj, pk, pl     ! Orbital indices within the occupied/virtual parts of the MP2 reference
    real(kind(mos))            :: wgt                
    integer(ik)                :: alloc              
    integer(ik), allocatable   :: occ(:)             ! List of MOs occupied in the MP2 reference determinant
    integer(ik), allocatable   :: vir(:)             ! List of MOs vacant in the MP2 reference determinant
    integer(ik), allocatable   :: dlt(:)             ! List of MOs where occupation numbers differ in FON-SCF and MP2 references
                                                     ! to the imaginary part of the energy. The choice is made -before- ISA is applied.
    real(kind(mos))            :: eps_bigplayer      ! Actual numerical threshold for eps_bigplayer
    character(len=20)          :: store              ! Integral storage accuracy
    !
    !  Items below are per-energy shift
    !
    complex(kind(mos)), allocatable   :: e_mp2s(:)       ! MP2(singles) energy correction; this must be zero if there are no fractional occupations
    complex(kind(mos)), allocatable   :: e_mp2d(:)       ! MP2(doubles) energy correction
    type(isa_statistics), allocatable :: isa_stat(:)     ! Summary of intruder-state avoidance activities
    logical, allocatable              :: bigplayer(:,:)  ! This MO contributed to a determinant which contributes more than eps_big*Im(eref)
    !
    call TimerStart('MP2 energy')
    norbs = size(occ_det)
    nocc  = sum (occ_det)
    nvir  = norbs - nocc
    ndlt  = count(abs(occ_det-occ_ref)>eps_occ)
    nisa  = size(isa_de)
    !
    write (out,"('MP2: ',i0,' orbitals, ',i0,' occupieds, ',i0,' virtuals, ',i0,' deltas, ',i0,' ISA shifts')") &
           norbs, nocc, nvir, ndlt, nisa
    !
    if (size(occ_ref)/=norbs .or. size(eps)/=norbs) stop 'mp2_tools%mp2_energy_incore - orbital count mismatch'
    if (any(occ_det<0) .or. any(occ_det>1)) stop 'mp2_tools%mp2_energy_incore - invalid occupations of MP2 reference'
    if (nisa<=0) stop 'mp2_tools%mp2_energy_incore - number of ISA shifts must be positive'
    !
    !  To enumerate the determinants, I need lists of occupied and virtual MOs in the MP2 reference
    !
    allocate (occ(nocc),vir(nvir),dlt(ndlt), &
             e_mp2s(nisa),e_mp2d(nisa),isa_stat(nisa),bigplayer(norbs,nisa),stat=alloc)
    if (alloc/=0) stop 'mp2_tools%mp2_energy_incore - allocation failed'
    pj = 0 ; pk = 0 ; pl = 0
    fill_occvir: do i=1,norbs
      if (occ_det(i)==1)                      then ; pj = pj + 1 ; occ(pj) = i ; end if
      if (occ_det(i)==0)                      then ; pk = pk + 1 ; vir(pk) = i ; end if
      if (abs(occ_det(i)-occ_ref(i))>eps_occ) then ; pl = pl + 1 ; dlt(pl) = i ; end if
    end do fill_occvir
    if (pj/=nocc) stop 'mp2_tools%mp2_energy_incore - occ count error'
    if (pk/=nvir) stop 'mp2_tools%mp2_energy_incore - vir count error'
    if (pl/=ndlt) stop 'mp2_tools%mp2_energy_incore - dlt count error'
    !
    bigplayer     = .false.
    eps_bigplayer = real(0.01_xrk,kind(mos)) * abs(aimag(eref))
    if (present(eps_big)) eps_bigplayer = eps_big * abs(aimag(eref))
    !
    !  2e integral transformation; we may need to do up to three separate blocks
    !
    store = 'as-is'
    if (present(storage_mode)) store = storage_mode
    call TimerStart('MP2 integral transformation')
    write (out,"(/t5,'Transforming VOVO block')")
    call transform_moint2e(int2e,'incore',mos(:,vir,1),mos(:,occ,2),mos(:,vir,1),mos(:,occ,2),m2e_vo,storage_mode=store)
    write (out,"(/t5,'Transforming OVOV block')")
    call transform_moint2e(int2e,'incore',mos(:,vir,2),mos(:,occ,1),mos(:,vir,2),mos(:,occ,1),m2e_ov,storage_mode=store)
    if (ndlt>0) then
      write (out,"(/t5,'Transforming XXXD block')")
      call transform_moint2e(int2e,'incore',mos(:,:,1),mos(:,:,2),mos(:,:,1),mos(:,dlt,2),m2e_dd,storage_mode=store)
    end if
    call TimerStop('MP2 integral transformation')
    !
    !  Figure out the first-order correction to the HF energy.
    !  It will be zero of the SCF part was done for the MP2 reference determinant.
    !  The code is copied from ci_tools%field_free_h_diagonal
    !  This part of MP2 energy does not depend on the ISA energy shift; do it once
    !
    bigplayer(dlt,:) = .true.           ! All orbitals in this block are by definition essential
    e_mp1 = sum(eps*(occ_det-occ_ref))  ! One-particle part of the energy
    hf_2e_i: do pi=1,ndlt
      i = dlt(pi)
      hf_2e_j: do pj=1,ndlt
        j = dlt(pj)
        if (i==j) cycle hf_2e_i
        wgt = real(0.5_xrk,kind(mos)) * (occ_det(i)*(occ_det(j)-occ_ref(j)) - (occ_det(i)-occ_ref(i))*occ_ref(j))
        select case (m2e_dd%ints_math)
          case default ; stop 'mp2_tools%mp2_energy_incore - bad ints_math'
          case ('real') ; e_mp1 = e_mp1 + wgt*(cmplx(m2e_dd%buffer_real(j,j,i,pi)-m2e_dd%buffer_real(i,j,j,pi),kind=kind(mos)))
          case ('quad') ; e_mp1 = e_mp1 + wgt*(cmplx(m2e_dd%buffer_quad(j,j,i,pi)-m2e_dd%buffer_quad(i,j,j,pi),kind=kind(mos)))
        end select
        isa_stat(:)%big_terms = isa_stat(:)%big_terms + 1
      end do hf_2e_j
    end do hf_2e_i
    !
    !$omp parallel do default(none) &
    !$omp& shared(nisa,nocc,nvir,ndlt,occ,vir,dlt,occ_det,occ_ref,eps,isa_de) &
    !$omp& shared(m2e_dd,m2e_vo,m2e_ov,eps_bigplayer,e_mp2s,e_mp2d,isa_stat,bigplayer) &
    !$omp& private(isa)
    scan_energy_shifts: do isa=1,nisa
      call incore_mp2_single_shift(nocc,nvir,ndlt,occ,vir,dlt,occ_det,occ_ref,eps,isa_de(isa), &
                                   m2e_dd,m2e_vo,m2e_ov,eps_bigplayer, &
                                   e_mp2s(isa),e_mp2d(isa),isa_stat(isa),bigplayer(:,isa))
    end do scan_energy_shifts
    !$omp end parallel do
    !
    write (out,"()")
      write (out,"('FON-SCF reference energy = ',g24.14,1x,g24.14)") eref
      write (out,"('   MP1 energy correction = ',g24.14,1x,g24.14)") e_mp1
    print_mp2_results: do isa=1,nisa
      write (out,"('        ISA energy shift = ',g24.14)") isa_de(isa)
      write (out,"('MP2(S) energy correction = ',g24.14,1x,g24.14)") e_mp2s(isa)
      write (out,"('MP2(D) energy correction = ',g24.14,1x,g24.14)") e_mp2d(isa)
      write (out,"('        Total MP2 energy = ',g24.14,1x,g24.14)") eref + e_mp1 + e_mp2s(isa) + e_mp2d(isa)
      call isa_report(real(isa_de(isa),kind=rk),isa_stat(isa)) ! It's OK to drop accuracy here - it's only for printing
    end do print_mp2_results
    write (out,"(/1x,i0,' occupieds (out of ',i0,') contributed more than ',g12.5,' to the lifetime')") &
           count(bigplayer(occ,1)), nocc, eps_bigplayer
    write (out,"( 1x,i0,' virtuals (out of ',i0,') contributed more than ',g12.5,' to the lifetime')") &
           count(bigplayer(vir,1)), nvir, eps_bigplayer
    write (out,"( 1x,i0,' determinants (out of ',i0,') contributed more than ',g12.5,' to the lifetime'/)") &
           isa_stat(1)%big_terms, isa_stat(1)%terms, eps_bigplayer
    call flush(out)
    !
    call destroy_moint2e(m2e_vo)
    call destroy_moint2e(m2e_ov)
    if (ndlt/=0) call destroy_moint2e(m2e_dd)
    deallocate (occ,vir,dlt,bigplayer,e_mp2s,e_mp2d,isa_stat)
    call TimerStop('MP2 energy')
! end subroutine mp2_energy_incore_real
