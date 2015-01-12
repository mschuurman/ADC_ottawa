! subroutine incore_mp2_single_shift(nocc,nvir,ndlt,occ,vir,dlt,occ_det,occ_ref,eps,isa_de,m2e_dd,m2e_vo,m2e_ov, &
!                                    eps_bigplayer,e_mp2s,e_mp2d,isa_stat,bigplayer)
    integer(ik), intent(in)             :: nocc           ! Number of occupied orbitals (aka electrons) in the MP2 reference determinant
    integer(ik), intent(in)             :: nvir           ! Number of virtual orbitals in the MP2 referecne determinant
    integer(ik), intent(in)             :: ndlt           ! Number of orbitals where occupation numbers differ in FON-SCF and MP2 references
    integer(ik), intent(in)             :: occ(:)         ! List of MOs occupied in the MP2 reference determinant
    integer(ik), intent(in)             :: vir(:)         ! List of MOs vacant in the MP2 reference determinant
    integer(ik), intent(in)             :: dlt(:)         ! List of MOs where occupation numbers differ in FON-SCF and MP2 references
    integer(sik), intent(in)            :: occ_det(:)     ! Orbital occupation numbers for the MP2 reference determinant
!   real(rk), intent(in)                :: occ_ref(:)     ! Orbital occupation numbers for the FON-SCF reference determinant
!   complex(rk), intent(in)             :: eps(:)         ! Orbital eigenvalues of the active orbitals
!   real(rk), intent(in)                :: isa_de         ! Value of the intruder-state avoidance parameter
    type(moint2e_cache), intent(in)     :: m2e_dd         ! 2-electron integrals of the X,X,X,D type
    type(moint2e_cache), intent(in)     :: m2e_vo         ! 2-electron integrals of the V,O,V,O type (LRLR)
    type(moint2e_cache), intent(in)     :: m2e_ov         ! 2-electron integrals of the O,V,O,V type (LRLR), O and V are swapped!
!   real(rk), intent(in)                :: eps_bigplayer  ! Actual numerical threshold for eps_bigplayer
!   complex(rk), intent(out)            :: e_mp2s         ! MP2(singles) energy correction; this must be zero if there are no fractional occupations
!   complex(rk), intent(out)            :: e_mp2d         ! MP2(doubles) energy correction
    type(isa_statistics), intent(inout) :: isa_stat       ! Summary of intruder-state avoidance activities
    logical, intent(inout)              :: bigplayer(:)   ! This MO contributed to a determinant which contributes more than eps_big*Im(eref)
                                                          ! to the imaginary part of the energy. The choice is made -before- ISA is applied.
    !
    !  Kind-agnostic below this line!
    !
    integer(ik)                         :: i, j, k, l     ! Orbital indices, actual MOs
    integer(ik)                         :: pi, pj, pk, pl ! Orbital indices within the occupied/virtual parts of the MP2 reference
    real(kind(eps))                     :: wgt            
    complex(kind(eps))                  :: vdr, vrd       ! Matrix elements for determinant<-reference and reference<-dterminant excitations
    complex(kind(eps))                  :: de0, de        ! Energy denominator before and after ISA
    !
    !  Single-excitation part of the MP2 energy. 
    !  This contribution must be zero if FON-SCF was for the MP2 reference determinant
    !
    e_mp2s    = 0
    mp2s_virtual: do pk=1,nvir
      k = vir(pk)
      mp2s_occupied: do pi=1,nocc
        i   = occ(pi)
        vdr = 0 ; vrd = 0
        mp2s_other: do pj=1,ndlt
          j = dlt(pj)
          if (j==i .or. j==k) cycle mp2s_other
          wgt = occ_det(j) - occ_ref(j)
          vdr = vdr + wgt*(get_int(m2e_dd,k,i,j,pj)-get_int(m2e_dd,j,i,k,pj))
          vrd = vrd + wgt*(get_int(m2e_dd,i,k,j,pj)-get_int(m2e_dd,j,k,i,pj))
        end do mp2s_other
        de0 = eps(i) - eps(k)
        de  = avoid_intruders(isa_de,de0,isa_stat)
        e_mp2s = e_mp2s + vrd * vdr / de
        !
        if (abs(de0)<=spacing(1._rk)) then ! Hard zero; automatically deemed important
          bigplayer((/i,k/)) = .true.      ! Note that our indices are guaranteed to be distinct
          isa_stat%big_terms = isa_stat%big_terms + 1
        else if (abs(aimag(vrd*vdr/de0))>=eps_bigplayer) then
          bigplayer((/i,k/)) = .true. 
          isa_stat%big_terms = isa_stat%big_terms + 1
        end if
      end do mp2s_occupied
    end do mp2s_virtual
    !
    !  Double-excitation part of the MP2 energy
    !
    e_mp2d = 0
    mp2d_vir_k: do pk=1,nvir
      k = vir(pk)
      mp2d_vir_l: do pl=pk+1,nvir
        l = vir(pl)
        mp2d_occ_i: do pi=1,nocc
          i = occ(pi)
          mp2d_occ_j: do pj=pi+1,nocc
            j   = occ(pj)
            vdr = get_int(m2e_vo,pk,pi,pl,pj) - get_int(m2e_vo,pl,pi,pk,pj)
            vrd = get_int(m2e_ov,pk,pi,pl,pj) - get_int(m2e_ov,pl,pi,pk,pj)
            de0 = eps(i) + eps(j) - eps(k) - eps(l)
            de  = avoid_intruders(isa_de,de0,isa_stat)
            ! write (out,"(4i4,' term = ',2f14.7,' vrd = ',2f14.7,' vdr = ',2f14.7,' de = ',2f14.7)") &
            !        i, j, k, l, vrd*vdr/de, vrd, vdr, de
            e_mp2d = e_mp2d + vrd*vdr/de
            !
            if (abs(de0)<=spacing(1._rk)) then ! Hard zero; automatically deemed important
              bigplayer((/i,j,k,l/)) = .true.  ! Note that our indices are guaranteed to be distinct
              isa_stat%big_terms = isa_stat%big_terms + 1
            else if (abs(aimag(vrd*vdr/de0))>=eps_bigplayer) then
              bigplayer((/i,j,k,l/)) = .true. 
              isa_stat%big_terms = isa_stat%big_terms + 1
            end if
          end do mp2d_occ_j
        end do mp2d_occ_i
      end do mp2d_vir_l
    end do mp2d_vir_k
    !
    contains
    !
    complex(kind(eps)) function get_int(i2e,i,j,k,l)
      type(moint2e_cache), intent(in) :: i2e        ! 2-e integrals
      integer(ik), intent(in)         :: i, j, k, l ! Indices 
      !
      select case (i2e%ints_math) 
        case default ; stop 'mp2_tools%incore_mp2_single_shift%get_int - bad ints_math'
        case ('real') ; get_int = cmplx(i2e%buffer_real(i,j,k,l),kind=kind(get_int))
        case ('quad') ; get_int = cmplx(i2e%buffer_quad(i,j,k,l),kind=kind(get_int))
      end select
    end function get_int
    !
    function avoid_intruders(isa_de,de,isa_stat) result (ade)
      real(kind(eps)), intent(in)         :: isa_de   ! Energy scale for intruder-state avoidance
      complex(kind(eps)), intent(in)      :: de       ! Energy denominator
      type(isa_statistics), intent(inout) :: isa_stat ! Counters for keeping track of what's going on
      complex(kind(eps))                  :: ade      ! Energy denominator, with the hard zero avoided
      !
      real(kind(eps)) :: x, rde
      !
      if (isa_de<=0) then
        ade = de
        return
      end if
      !
      x   = abs(de)
      rde = real(de,kind=kind(eps))
      !
      isa_stat%terms = isa_stat%terms + 1
           if (x<= real(0.1_xrk,kind(eps))*isa_de) then
        if (rde>=0) isa_stat%count_p001 = isa_stat%count_p001 + 1
        if (rde< 0) isa_stat%count_m001 = isa_stat%count_m001 + 1
      else if (x<= real(1.0_xrk,kind(eps))*isa_de) then
        if (rde>=0) isa_stat%count_p010 = isa_stat%count_p010 + 1
        if (rde< 0) isa_stat%count_m010 = isa_stat%count_m010 + 1
      else if (x<=real(10.0_xrk,kind(eps))*isa_de) then
        if (rde>=0) isa_stat%count_p100 = isa_stat%count_p100 + 1
        if (rde< 0) isa_stat%count_m100 = isa_stat%count_m100 + 1
      endif
      !
      if (x<=spacing(isa_de)) then
        !
        !  Hard zero, effectively drop the term
        !
        ade = real(1e20_xrk,kind(eps))
        isa_stat%hard_zeros = isa_stat%hard_zeros + 1
      else
        ade = de * (1 + (isa_de/x)**2)
      end if
    end function avoid_intruders
! end subroutine incore_mp2_single_shift
