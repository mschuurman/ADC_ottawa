!   function mt_matmul_crrr(a,b,serial) result(c)
!     complex(rk), intent(in)       :: a(:,:) !
!     real   (rk), intent(in)       :: b(:,:) !
      logical, intent(in), optional :: serial ! Set to true to tile matrix multiplication without running in parallel
!     complex(rk)                   :: c(size(a,dim=1),size(b,dim=2))
!     !
      integer(ik)              :: is, ie, js, je, ks, ke ! Block indices - initial and starting
      integer(ik)              :: n, l, m
      integer(ik)              :: blocks, ib
      integer(ik)              :: alloc
      integer(ik), allocatable :: ijs(:,:)
      logical                  :: para
      !
      if (size(a,dim=2)/=size(b,dim=1) .or. size(c,dim=1)/=size(a,dim=1) .or. size(c,dim=2)/=size(b,dim=2)) then
        stop 'matrix_tools%mt_matmul - matrices not conformable'
      end if
      !
      para = .true.
      if (present(serial)) para = .not.serial
      n = size(a,dim=1) ; l = size(a,dim=2) ; m = size(b,dim=2)
      !
      !  Is it worth running in parallel?
      !
      blocks = ((n+stripe-1)/stripe) * ((m+stripe-1)/stripe)
      if (blocks<min_blocks .or. .not.para) then
        c = matmul(a,b)
        return
      end if
      !
      call TimerStart('mt_matmul')
      !
      !  Prepare worklist to maximize the number of work units available
      !
      allocate (ijs(2,blocks),stat=alloc)
      if (alloc/=0) stop 'matrix_tools%mt_matmul - allocation failed'
      ib = 0
      row_c: do is=1,n,stripe
        column_c: do js=1,m,stripe
          ib = ib + 1
          ijs(:,ib) = (/ is, js /)
        end do column_c
      end do row_c
      if (ib/=blocks) stop 'matrix_tools%mt_matmul - counting error'
      !
      c = 0
      column_a: do ks=1,l,stripe
      ke = min(ks+stripe-1,l)
        !
        !  Each block of c() inside this loop is updated independently, so that
        !  we can do a lock-free SMP parallelization now
        !
        !$omp parallel do default(none) shared(a,b,c,ks,ke,ijs,m,n,blocks) private(ib,is,ie,js,je) if(para)
        workblocks: do ib=1,blocks
          is = ijs(1,ib) ; ie = min(is+stripe-1,n)
          js = ijs(2,ib) ; je = min(js+stripe-1,m)
          c(is:ie,js:je) = c(is:ie,js:je) + matmul(a(is:ie,ks:ke),b(ks:ke,js:je))
        end do workblocks
        !$omp end parallel do
      end do column_a
      deallocate (ijs)
      call TimerStop('mt_matmul')
!   end function mt_matmul_crrr
