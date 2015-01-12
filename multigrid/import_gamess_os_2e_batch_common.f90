!   subroutine os_2e_batch_real(what,v,xyz,sh_l,sh_ns,sh_np,p_zet,p_c,p_c_max,p_master,p_z_same,cutoff_2e)
!     type(gam_operator_data), intent(in) :: what              ! Operator parameters
!     real(rk), intent(out)               :: v(:,:,:,:)        ! Buffer for the results
!     real(ark), intent(in)               :: xyz(:,:)          ! Coordinates of the basis functions within each batch
!     integer(ik), intent(in)             :: sh_l (:)          ! Angular momentum of each shell within a batch is the same
!     integer(ik), intent(in)             :: sh_ns(:)          ! Number of shells in each batch
!     integer(ik), intent(in)             :: sh_np(:,:)        ! Number of primitives in each shell
!     real(ark), intent(in)               :: p_zet   (:,:,:)   ! Primitive exponents in each shell
!     real(ark), intent(in)               :: p_c     (:,:,:)   ! Primitive contraction coefficients in each shell
!     real(ark), intent(in)               :: p_c_max (:,:,:)   ! Absolute maximum of a contraction coefficient for all instances
!                                                              ! of a primitive. This value is only defined for the "master" entry
!     logical, intent(in)                 :: p_master(:,:,:)   ! "Master" copy of the exponent
!     integer(ik), intent(in)             :: p_z_same(:,:,:,:) ! Next contraction with an identical exponent
!                                                              ! First index: 1 = shell index; 2 = contraction within shell
!                                                              ! 2nd, 3rd, and 4th indices are as the 1st etc in p_zet/p_c/p_master
!     real(rk), intent(in)                :: cutoff_2e         ! Absolute accuracy required in the final integrals;
!                                                              ! Negative means retain full accuracy
!     !
      integer(ik)                 :: a, b, c, d      ! shell index within a batch
      integer(ik)                 :: i, j, k, l      ! primitive index within a shell
      integer(ik)                 :: p1(4)           ! Starting position of the desired block of primitive integrals
      integer(ik)                 :: p2(4)           ! Ending position of the desired block of primitive integrals
      integer(ik)                 :: sz(4)           ! Size of the desired block
      integer(ik)                 :: nrec            ! Maximum recursion order
      integer(ik)                 :: alloc           ! Allocation status, what else?
      real(kind(v)), allocatable  :: vb_r(:,:,:,:,:) ! Buffer for recursion
      real(kind(v))               :: wgt
      real(kind(v))               :: zet(4)          ! Primitive exponents
      real(kind(v))               :: w_cut           ! Absolute cutoff 
      logical                     :: zero_primitives ! True if all primitive integrals are negligible
      real(kind(v))               :: ang_c_kind(lbound(ang_c,dim=1):ubound(ang_c,dim=1))
      !
      !  Prepare angular factors of the right kind
      !
      ang_c_kind = real(ang_c,kind=kind(ang_c_kind))
      !
      !  All shells within the batch have the same angular momentum, so that the recursion buffer
      !  does not have to be resized.
      !
      p1   = ang_loc(sh_l)
      p2   = ang_loc(sh_l+1)-1
      sz   = p2-p1+1
      nrec = sum(sh_l)
      allocate (vb_r(0:p2(1),0:p2(2),0:p2(3),0:p2(4),0:nrec),stat=alloc)
      if (alloc/=0) then
        write (out,"('import_gamess%os_2e_batch: allocation failed with code ',i10,'. p2, nrec = ',5i8)") alloc, p2, nrec
        stop 'import_gamess%os_2e_batch: Buffer deallocation failed'
      end if
      !
      !  Figure out the cut-off point for the primitives; we need to include the angular
      !  scaling factors since these can vary by nearly an order of magnitude.
      !  The last angular scaling coefficient in a batch is always the largest one.
      !
      w_cut = cutoff_2e / product(ang_c_kind(p2))
      !
      v = 0
      !
      !  We do not want to process any of the primitives more than once; hence this
      !  very complicated loop. It gets even worse during accumulation phase ...
      !
      contraction_1s: do a=1,sh_ns(1)
        contraction_1p: do i=1,sh_np(a,1)
          if (.not.p_master(i,a,1)) cycle contraction_1p ! Will be handled elsewhere
          zet(1) = p_zet(i,a,1)
          contraction_2s: do b=1,sh_ns(2)
            contraction_2p: do j=1,sh_np(b,2)
              if (.not.p_master(j,b,2)) cycle contraction_2p ! Will be handled elsewhere
              zet(2) = p_zet(j,b,2)
              contraction_3s: do c=1,sh_ns(3)
                contraction_3p: do k=1,sh_np(c,3)
                  if (.not.p_master(k,c,3)) cycle contraction_3p ! Will be handled elsewhere
                  zet(3) = p_zet(k,c,3)
                  contraction_4s: do d=1,sh_ns(4)
                    contraction_4p: do l=1,sh_np(d,4)
                      if (.not.p_master(l,d,4)) cycle contraction_4p ! Will be handled elsewhere
                      zet(4) = p_zet(l,d,4)
                      wgt    = p_c_max(i,a,1)*p_c_max(j,b,2)*p_c_max(k,c,3)*p_c_max(l,d,4)
                      call os_2e_primitives(what,xyz,zet,p2,nrec,vb_r,w_cut/wgt,zero_primitives)
                      if (zero_primitives) cycle contraction_4p
                      call accumulate_all
                    end do contraction_4p
                  end do contraction_4s
                end do contraction_3p
              end do contraction_3s
            end do contraction_2p
          end do contraction_2s
        end do contraction_1p
      end do contraction_1s
      !
      deallocate (vb_r,stat=alloc)
      if (alloc/=0) then
        stop 'import_gamess%os_2e_batch: Buffer deallocation failed'
      end if
      !
      contains
      !
      !  Accumulate contibution of these primitives to all contractions containing them
      !
      subroutine accumulate_all
        integer(ik)   :: ea, eb, ec, ed ! "Equivalent" shell indices
        integer(ik)   :: ei, ej, ek, el ! "Equivalent" contraction indices
        integer(ik)   :: pi, pj, pk, pl ! Position of the contacted shell within the output buffer
        real(kind(v)) :: ca, cb, cc, cd ! Contraction coefficients
        real(kind(v)) :: wgt
        !
        !  We have to work through all instances of this primitive. We know there is
        !  at least one instance, but there could be more ...
        !
        ea = a ; ei = i
        primitive_1: do while (ea/=0)
          ca = p_c(ei,ea,1)
          pi = (ea-1)*sz(1) + 1
          eb = b ; ej = j
          primitive_2: do while (eb/=0)
            cb = p_c(ej,eb,2)
            pj = (eb-1)*sz(2) + 1
            ec = c ; ek = k
            primitive_3: do while (ec/=0)
              cc = p_c(ek,ec,3)
              pk = (ec-1)*sz(3) + 1
              ed = d ; el = l
              primitive_4: do while (ed/=0)
                cd = p_c(el,ed,4)
                pl = (ed-1)*sz(4) + 1
                wgt = ca*cb*cc*cd
                call accumulate(pi,pj,pk,pl,wgt)
                call next_instance(ed,el,4_ik)
              end do primitive_4
              call next_instance(ec,ek,3_ik)
            end do primitive_3
            call next_instance(eb,ej,2_ik)
          end do primitive_2
          call next_instance(ea,ei,1_ik)
        end do primitive_1
      end subroutine accumulate_all
      !
      subroutine next_instance(shell,primitive,ind)
        integer(ik), intent(inout) :: shell     ! Shell index within a batch
        integer(ik), intent(inout) :: primitive ! Primitive index within a shell
        integer(ik), intent(in)    :: ind       ! Centre index
        !
        integer(ik) :: a, i
        !
        a = shell ; i = primitive
        if (a<=0 .or. a>sh_ns  (ind)) stop 'import_gamess%os_2e_batch%next_instance - bad shell index'
        if (i<=0 .or. i>sh_np(a,ind)) stop 'import_gamess%os_2e_batch%next_instance - bad primitive index'
        !
        shell     = p_z_same(1,i,a,ind)
        primitive = p_z_same(2,i,a,ind)
      end subroutine next_instance
      !
      !  Accumulate contibution of these primitives to a single contaction
      !
      subroutine accumulate(pi,pj,pk,pl,wgt)
        integer(ik), intent(in)   :: pi, pj, pk, pl ! Starting positions within output buffer
        real(kind(v)), intent(in) :: wgt            ! Weight of this primitive
        !
        integer(ik) :: ip, jp, kp, lp  ! Indices within the recursion buffer
        integer(ik) :: iq, jq, kq, lq  ! Indices within the output buffer
        !
        add_1: do ip=p1(1),p2(1)
          iq = pi + ip - p1(1)
          add_2: do jp=p1(2),p2(2)
            jq = pj + jp - p1(2)
            add_3: do kp=p1(3),p2(3)
              kq = pk + kp - p1(3)
              add_4: do lp=p1(4),p2(4)
                lq = pl + lp - p1(4)
                v(iq,jq,kq,lq) = v(iq,jq,kq,lq) &
                               + wgt*ang_c_kind(ip)*ang_c_kind(jp)*ang_c_kind(kp)*ang_c_kind(lp) * vb_r(ip,jp,kp,lp,0)
              end do add_4
            end do add_3
          end do add_2
        end do add_1
      end subroutine accumulate
!   end subroutine os_2e_batch_real
