  program matmul_test
    use accuracy
    use timer
    !
    integer(ik) :: sz, niter, i
    !
    real(xrk), allocatable    :: tmp_r(:,:), tmp_i(:,:)
    complex(xrk), allocatable :: a(:,:), b(:,:), c(:,:)
    !
    read (input,*) sz, niter
    allocate (a(sz,sz), b(sz,sz), c(sz,sz), tmp_r(sz,sz), tmp_i(sz,sz))
    !
    call random_number(tmp_r) ; call random_number(tmp_i)
    a = cmplx(tmp_r,tmp_i,kind=xrk)
    call random_number(tmp_r) ; call random_number(tmp_i)
    b = cmplx(tmp_r,tmp_i,kind=xrk)
    !
    c = matmul(a,b) ! warmup
    call TimerStart('matmul')
    do i=1,niter
      c = matmul(a,b)
    end do
    call TimerStop('matmul')
    !
    call xgemm('N','N',sz,sz,sz,1.0_xrk,a,sz,b,sz,0.0_xrk,c,sz) ! warmup
    call TimerStart('xgemm')
    do i=1,niter
      call xgemm('N','N',sz,sz,sz,1.0_xrk,a,sz,b,sz,0.0_xrk,c,sz)
    end do
    call TimerStop('xgemm')
    !
!   call quad_zgemm('N','N',sz,sz,sz,1.0_xrk,a,sz,b,sz,0.0_xrk,c,sz) ! warmup
!   call TimerStart('quad_zgemm')
!   do i=1,niter
!     call quad_zgemm('N','N',sz,sz,sz,1.0_xrk,a,sz,b,sz,0.0_xrk,c,sz)
!   end do
!   call TimerStop('quad_zgemm')
    !
    call TimerReport
  end program matmul_test
