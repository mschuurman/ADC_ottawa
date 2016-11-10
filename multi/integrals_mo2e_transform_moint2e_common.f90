!   subroutine transform_moint2e_real(int2e,mode,mo_i,mo_j,mo_k,mo_l,moint2e,io_unit,l_block,storage_mode)
!     type(int2e_cache), intent(inout)       :: int2e        ! 2E integrals over spin-less AOs
!     character(len=*), intent(in)           :: mode         ! Mode of operation; either 'incore' or 'disk'
!     complex(rk), intent(in)                :: mo_i(:,:)    ! Spin-MOs to be used for the first integral index
!                                                            ! Index 1: all spin-alpha AOs, followed by all spin-beta AOs
!                                                            ! Index 2: MO index
!     complex(rk), intent(in)                :: mo_j(:,:)    ! Spin-MOs to be used for the second integral index
!     complex(rk), intent(in)                :: mo_k(:,:)    ! Spin-MOs to be used for the third integral index
!     complex(rk), intent(in)                :: mo_l(:,:)    ! Spin-MOs to be used for the fourth integral index
!                                                            ! Note that all mo_* arrays must have the same first dimension;
!                                                            ! the second dimension however does not need to be identical
!     type(moint2e_cache), intent(inout)     :: moint2e      ! MO integral descriptor
!     integer(ik), intent(in), optional      :: io_unit      ! I/O unit for storing MO integrals; not needed for mode='incore'
!     integer(ik), intent(in), optional      :: l_block      ! Number of the first-index MOs to transform at the same time 
!                                                            ! (this is an I/O optimization). The default is 5.
!     character(len=*), intent(in), optional :: storage_mode ! Override precision of the final transformed integrals;
!                                                            ! 'real', 'quad', or 'as-is'. The default is the same kind
!                                                            ! as used for the MOs (ie storage_mode='as-is')
      !
      !  Below this line we should not assume a specific kind of the real/complex types
      !
      integer(ik)               :: ios, alloc, kind_bytes, store_bytes
      integer(ik)               :: nao_spin, nao
      integer(ik)               :: l_count, l, l_pos, i, j, k
      integer(ik)               :: l_in_pass, lb1, lbn, l_pos1, l_posn
      integer(hik)              :: rec_len
      !
      complex(rk), pointer              :: buf_r  (:,:,:,:)   ! Real buffer for fully-transformed integrals; stored in moint2e
      complex(xrk), pointer             :: buf_q  (:,:,:,:)   ! Quad buffer for fully-transformed integrals; stored in moint2e
      complex(kind(mo_i)), allocatable  :: buf_l  (:,:,:,:,:) ! Integrals transformed over (a single value of) the last index 
                                                              ! The first three indices are spin-less AOs; the fourth index
                                                              ! is 1 (spin-alpha) or 2 (spin-beta); the fifth index is L
      complex(kind(mo_i)), allocatable  :: buf_kl (:,:,:)     ! Integrals transformed over two last indices
      complex(kind(mo_i)), allocatable  :: buf_jkl(:,:,:,:)   ! Integrals transformed over three last indices
                                                              ! Fourth index is spin
      !
      call TimerStart('AO>MO 2E transformation')
      !
      nao_spin = size(mo_i,dim=1) ! Must be a multiple of 2!
      nao      = nao_spin / 2
      if (nao_spin/=2*nao) stop 'integrals_mo2e%transform_moint2e - number of spin-AO is not even?!'
      if (size(mo_j,dim=1)/=nao_spin .or. &
          size(mo_k,dim=1)/=nao_spin .or. &
          size(mo_l,dim=1)/=nao_spin) stop 'integrals_mo2e%transform_moint2e - number of AOs is not the same for I/J/K/L'
      !
      moint2e%nmo = (/ size(mo_i,dim=2), size(mo_j,dim=2), size(mo_k,dim=2), size(mo_l,dim=2) /)
      !
      l_in_pass = 5
      if (present(l_block)) l_in_pass = l_block
      l_in_pass = min(max(1,l_in_pass),moint2e%nmo(4))
      !
      moint2e%mode = mode
      select case (moint2e%mode)
        case default
          write (out,"('MO integrals handling mode ',a,' is not recognized')") trim(moint2e%mode)
          stop 'integrals_mo2e%transform_moint2e - bad mode'
        case ('incore')
          l_count = moint2e%nmo(4)
        case ('disk')
          l_count = l_in_pass
      end select
      ! 
      !  It is tempting to use a switch here, but since rk and xrk are not necessarily
      !  distinct, this would be a mistake. This should be the last place we examine 
      !  kinds explicitly!
      !
      if (kind(mo_i)==rk) then
        moint2e%ints_math = 'real'
        kind_bytes        = rk_bytes
      else if (kind(mo_i)==xrk) then
        moint2e%ints_math = 'quad'
        kind_bytes        = xrk_bytes
      else
        stop 'integrals_mo2e%transform_moint2e - unsupported complex kind'
      end if
      nullify (buf_r,buf_q,moint2e%buffer_real,moint2e%buffer_quad)
      !
      !  Storage accuracy override, if present
      !
      store_bytes = kind_bytes
      if (present(storage_mode)) then
        select case (storage_mode)
          case default
            write (out,"('2E MO integrals storage mode ',a,' is not recognized')") trim(storage_mode)
            stop 'integrals_mo2e%transform_moint2e - bad storage mode'
          case ('as-is')
          case ('real')
            moint2e%ints_math = 'real'
            store_bytes       = rk_bytes
          case ('quad')
            moint2e%ints_math = 'quad'
            store_bytes       = xrk_bytes
        end select
      end if
      !
      write (out,"(/'Integral transformation will use ',i0,'-byte reals, and store the result in ',i0,'-byte reals.')") &
             kind_bytes, store_bytes
      if (store_bytes>kind_bytes) then
        write (out,"('WARNING: MO integrals transformed with less than stored accuracy.')")
      end if
      write (out,"( '2E MO integrals require ',f0.3,'-Gbyte buffer')") &
             2*store_bytes*l_count*product(real(moint2e%nmo(1:3),kind=rk))/(1024._rk**3)
      write (out,"( 'Integral transformation additionally requires ',f0.3,' Gbytes of memory')") &
            2*kind_bytes*(  l_in_pass*2*real(nao,kind=rk)**3 &
                          + moint2e%nmo(3)*real(nao,kind=rk)**2 &
                          + 2*nao*product(real(moint2e%nmo(2:3),kind=rk)))/(1024._rk**3)
      write (out,"( 'Transformation will include ',i0,' orbitals per pass')") l_in_pass
      call flush(out)
      !
      !  Integral buffer; may be kept if in-core
      !
      select case (moint2e%ints_math)
        case default ; stop 'integrals_mo2e%transform_moint2e - bad ints_math (1)'
        case ('real') ; allocate (buf_r(moint2e%nmo(1),moint2e%nmo(2),moint2e%nmo(3),l_count),stat=alloc)
        case ('quad') ; allocate (buf_q(moint2e%nmo(1),moint2e%nmo(2),moint2e%nmo(3),l_count),stat=alloc)
      end select
      if (alloc/=0) then 
        write (out,"('Error ',i0,' allocating memory')") alloc
        stop 'integrals_mo2e%transform_moint2e - memory allocation failed (1)'
      end if
      !
      !  Intermediates; will be released later
      !
      allocate (buf_l(nao,nao,nao,2,l_in_pass), buf_kl(nao,nao,moint2e%nmo(3)), &
                buf_jkl(nao,moint2e%nmo(2),moint2e%nmo(3),2), stat=alloc)
      if (alloc/=0) then 
        write (out,"('Error ',i0,' allocating memory')") alloc
        stop 'integrals_mo2e%transform_moint2e - memory allocation failed (2)'
      end if
      !
      !  Since all MO integral blocks are the same size, there is no harm in storing them in
      !  a direct-access file; this way, we can also fetch any integral block we need.
      !
      if (mode=='disk') then
        if (.not.present(io_unit)) stop 'integrals_mo2e%transform_moint2e - io_unit argument missing'
        select case (moint2e%ints_math)
          case default ; stop 'integrals_mo2e%transform_moint2e - bad ints_math (2)'
          case ('real') ; inquire (iolength=rec_len) buf_r(:,:,:,1)
          case ('quad') ; inquire (iolength=rec_len) buf_q(:,:,:,1)
        end select
        moint2e%io_unit = io_unit
        open (moint2e%io_unit,form='unformatted',action='readwrite',status='scratch',access='direct',recl=rec_len,iostat=ios)
        if (ios/=0) then
          write (out,"('Error ',i0,' creating scratch file for 2E MO integrals')") ios
          stop 'integrals_mo2e%transform_moint2e - error creating MO 2E integral file'
        end if
      end if
      !
      batches_l: do lb1=1,moint2e%nmo(4),l_in_pass
        lbn = min(moint2e%nmo(4),lb1+l_in_pass-1)
        write (out,"('Transforming integrals over MOs ',i0,'-',i0)") lb1, lbn
        call flush (out)
        l_pos1 = 1
        if (mode=='incore') l_pos1 = lb1
        l_posn = l_pos1 + (lbn-lb1)
        select case (trim(int2e%ints_math)//' '//moint2e%ints_math)
          case default; stop 'integrals_mo2e%transform_moint2e - unsupported kind of AO and/or MO integrals'
          !
          !  transform_one_orbital is a bit of a misnomer - we may be doing a block of orbitals on one pass ...
          ! 
          case ('real real') 
            call transform_one_orbital(int2e,nao,moint2e%nmo,mo_i,mo_j,mo_k, &
                                       mo_l(:,lb1:lbn),buf_r(:,:,:,l_pos1:l_posn),buf_jkl,buf_kl,buf_l,1._rk)
          case ('real quad') 
            call transform_one_orbital(int2e,nao,moint2e%nmo,mo_i,mo_j,mo_k, &
                                       mo_l(:,lb1:lbn),buf_q(:,:,:,l_pos1:l_posn),buf_jkl,buf_kl,buf_l,1._rk)
          case ('quad real') 
            call transform_one_orbital(int2e,nao,moint2e%nmo,mo_i,mo_j,mo_k, &
                                       mo_l(:,lb1:lbn),buf_r(:,:,:,l_pos1:l_posn),buf_jkl,buf_kl,buf_l,1._xrk)
          case ('quad quad') 
            call transform_one_orbital(int2e,nao,moint2e%nmo,mo_i,mo_j,mo_k, &
                                       mo_l(:,lb1:lbn),buf_q(:,:,:,l_pos1:l_posn),buf_jkl,buf_kl,buf_l,1._xrk)
        end select
        if (mode=='disk') then
          dump_orbitals_one_by_one: do l_pos=l_pos1,l_posn
            l = lb1 + (l_pos-l_pos1)
            select case (moint2e%ints_math)
              case default ; stop 'integrals_mo2e%transform_moint2e - bad ints_math (3)'
              case ('real') ; write (moint2e%io_unit,rec=l,iostat=ios) buf_r(:,:,:,l_pos)
              case ('quad') ; write (moint2e%io_unit,rec=l,iostat=ios) buf_q(:,:,:,l_pos)
            end select
            if (ios/=0) then
              write (out,"('Error ',i0,' writing transformed integrals')") ios
              stop 'integrals_mo2e%transform_moint2e - I/O error writing transformed integrals'
            end if
          end do dump_orbitals_one_by_one
        end if
        if (verbose>=2) then
          write (out,"(/t5,'2e integrals over active MOs'/)")
          write (out,"(4(1x,a4),2x,a20,1x,a20)") ' I ', ' J ', ' K ', ' L ', ' Re (IJ,KL) ', ' Im (I,J,K,L) ', &
                                                 '---', '---', '---', '---', '------------', '--------------'
          print_l: do l=lb1,lbn
            l_pos = l_pos1 + (l-lb1)
            print_k: do k=1,moint2e%nmo(3)
              print_j: do j=1,moint2e%nmo(2)
                print_i: do i=1,moint2e%nmo(1)
                  select case (moint2e%ints_math)
                    case default ; stop 'integrals_mo2e%transform_moint2e - bad ints_math (4)'
                    case ('real') ; write (out,"(4(1x,i4),2x,g20.12,1x,g20.12)") i, j, k, l, buf_r(i,j,k,l_pos)
                    case ('quad') ; write (out,"(4(1x,i4),2x,g20.12,1x,g20.12)") i, j, k, l, buf_q(i,j,k,l_pos)
                  end select
                end do print_i
              end do print_j
            end do print_k
          end do print_l
          write (out,"()")
        end if
      end do batches_l
      !
      deallocate (buf_l,buf_kl,buf_jkl)  ! buf(:,:,:,:) goes back to the user, and does not have to be destroyed here
      !
      if (moint2e%mode=='disk') then
        if (associated(buf_r)) deallocate (buf_r) ! We may have had multiple orbitals!
        if (associated(buf_q)) deallocate (buf_q) ! We may have had multiple orbitals!
        select case (moint2e%ints_math)
          case default ; stop 'integrals_mo2e%transform_moint2e - bad ints_math (5)'
          case ('real') ; allocate (buf_r(moint2e%nmo(1),moint2e%nmo(2),moint2e%nmo(3),1),stat=alloc)
          case ('quad') ; allocate (buf_q(moint2e%nmo(1),moint2e%nmo(2),moint2e%nmo(3),1),stat=alloc)
        end select
        if (alloc/=0) then 
          write (out,"('Error ',i0,' allocating memory')") alloc
          stop 'integrals_mo2e%transform_moint2e - memory allocation failed (3)'
        end if
        moint2e%mo_l = -1
      end if
      !
      !  One of the two pointers will have the actual buffer; the other will be NULL()
      !
      moint2e%buffer_real => buf_r
      moint2e%buffer_quad => buf_q
      !
      call TimerStop('AO>MO 2E transformation')
      ! call TimerReport
!   end subroutine transform_moint2e_real
