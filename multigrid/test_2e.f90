  subroutine driver
    use accuracy
    use timer
    use math
    use import_gamess
    !
    integer(ik)           :: info_2e(3)
    integer(ik)           :: nbatch, maxbatch, nbas
    integer(ik)           :: bi, bj, bk, bl
    integer(ik)           :: si, sj, sk, sl
    integer(ik)           :: zi, zj, zk, zl
    integer(ik)           ::  i,  j,  k,  l
    integer(ik)           :: ios
    real(ark)             :: compares, equals
    real(rk), allocatable :: buf (:,:,:,:)
    real(rk), allocatable :: ints(:,:,:,:)
  ! real(rk), allocatable :: tkin(:,:)
    real(rk)              :: ref
    !
    call accuracyInitialize
    call TimerStart('Test 2e')
    ref = MathFactorial(80)
    ref = MathDoubleFactorial(80)
    ref = MathLogFactorial(80)
    call gamess_load_orbitals(file='test_2e.dat')
    call gamess_2e_info('GET INFO',op_iparam=info_2e)
    nbatch   = info_2e(1)
    maxbatch = info_2e(2)
    nbas     = info_2e(3)
    write (out,"('Batches = ',i8,' max. batch = ',i8,' total functions = ',i8)") nbatch, maxbatch, nbas
    !
    !  Tag kinetic integrals along - it's too much work to do this properly ;-)
    !
    ! allocate (tkin(nbas,nbas))
    ! call gamess_1e_integrals('AO KINETIC',tkin)
    ! call gamess_print_1e_integrals(tkin)
    ! deallocate (tkin)
    !
    allocate (ints(nbas,nbas,nbas,nbas))
    !
    !$omp parallel default(none) private(bl,sl,zl,bk,sk,zk,bj,sj,zj,bi,si,zi,buf) shared(nbatch,ints,maxbatch)
    allocate (buf(maxbatch,maxbatch,maxbatch,maxbatch))
    !$omp do schedule(guided)
    batch_l: do bl=1,nbatch
      write (out,"('L-batch: ',i8)") bl
      call batch_size(bl,sl,zl)
      batch_k: do bk=1,nbatch
        call batch_size(bk,sk,zk)
        batch_j: do bj=1,nbatch
          call batch_size(bj,sj,zj)
          batch_i: do bi=1,nbatch
            call batch_size(bi,si,zi)
            call gamess_2e_integrals('AO 4C 1/R',buf,(/bi,bj,bk,bl/))
            ints(zi:zi+si-1,zj:zj+sj-1,zk:zk+sk-1,zl:zl+sl-1) = buf(:si,:sj,:sk,:sl)
            !
          end do batch_i
        end do batch_j
      end do batch_k
    end do batch_l
    !$omp end do
    deallocate (buf)
    !$omp end parallel
    !
    !  Symmetry check. Swapping i and j; k and l; and ij and kl should not change
    !  the integral. However, the actual recursion used was likely not the same,
    !  so there may be a bit of numerical noise
    !
    symm_l: do l=1,nbas
      symm_k: do k=1,nbas
        symm_j: do j=1,nbas
          symm_i: do i=1,nbas
            call sym_check(i,j,k,l, j,i,k,l)
            call sym_check(i,j,k,l, i,j,l,k)
            call sym_check(i,j,k,l, j,i,l,k)
            call sym_check(i,j,k,l, k,l,i,j)
            call sym_check(i,j,k,l, l,k,i,j)
            call sym_check(i,j,k,l, k,l,j,i)
            call sym_check(i,j,k,l, l,k,j,i)
          end do symm_i
        end do symm_j
      end do symm_k
    end do symm_l
    !
  ! print_l: do l=1,nbas
  !   print_k: do k=1,nbas
  !     print_j: do j=1,nbas
  !       print_i: do i=1,nbas
  !         write (out,"(4i8,2x,f25.15)") i, j, k, l, ints(i,j,k,l)
  !       end do print_i
  !     end do print_j
  !   end do print_k
  ! end do print_l
    write (out,"('Done with the symmetry test')")
    compares = 0 ; equals = 0 ;
    test_loop: do 
      read(input,*,iostat=ios) i, j, k, l, ref
      if (ios/=0) exit test_loop
      compares = compares + 1
!       write (out,"(4i8,2x,3f25.15)") i, j, k, l, ints(i,j,k,l), ref, ints(i,j,k,l)-ref
      if (abs(ints(i,j,k,l)-ref)>=1e-7_rk) then
        write (out,"(4i8,2x,3f25.15)") i, j, k, l, ints(i,j,k,l), ref, ints(i,j,k,l)-ref
      else
        equals = equals + 1
      end if
    end do test_loop
    write (out,"('Total comparisons = ',f14.0,' equals = ',f14.0)") compares, equals
    !
    deallocate (ints)
    call gamess_destroy
    call TimerStop('Test 2e')
    call TimerReport
    !
    contains
    !
    subroutine batch_size(ib,is,ip)
      integer(ik), intent(in)  :: ib
      integer(ik), intent(out) :: is, ip
      integer(ik)              :: info_2e(3)
      !
      info_2e(1) = ib
      call gamess_2e_info('GET BATCH INFO',op_iparam=info_2e)
      is = info_2e(2)  ! Batch size
      ip = info_2e(3)  ! Batch position
    end subroutine batch_size
    !
    subroutine sym_check(ai,aj,ak,al,  bi,bj,bk,bl)
      integer(ik), intent(in) :: ai,aj,ak,al
      integer(ik), intent(in) :: bi,bj,bk,bl
      !
      real(rk) :: v1, v2, eps
      !
      v1  = ints(ai,aj,ak,al)
      v2  = ints(bi,bj,bk,bl)
      eps = spacing(100._rk*max(1.0_rk,abs(v1),abs(v2)))
      if (abs(v1-v2)<=eps) return
      write (out,"(2(4i8,2x),3(1x,f25.15))") ai,aj,ak,al,  bi,bj,bk,bl,  v1,v2, v2-v1
    end subroutine sym_check
  end subroutine driver
