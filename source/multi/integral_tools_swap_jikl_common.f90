!   function swap_jikl_real(bi,p0,sz,p2,s2,a2e,a22) result(go)
!     real(rk), intent(in)     :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices (i<j) < (k<l)
!     real(rk), intent(out)    :: a22(:,:,:,:)   ! Index-permuted integrals
      integer(ik), intent(in)  :: bi(:)          ! Integral block index
      integer(ik), intent(in)  :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)  :: sz(:)          ! Size of the integral block
      integer(ik), intent(out) :: p2(4)          ! Same as p0, after index swapping 
      integer(ik), intent(out) :: s2(4)          ! Same as sz, after index swapping
      logical                  :: go             ! Set to true if swap is allowed
      !
      integer(ik)             :: ci, cj, ck, cl  ! Syntax sugar for bi(1), bi(2) etc
      integer(ik), parameter  :: code(4) = (/2,1,3,4/)
      integer(ik)             :: i, j, k, l
      !
      ci = bi(1) ; cj = bi(2) ; ck = bi(3) ; cl = bi(4)
      go = (ci/=cj)
      if (.not.go) return
      !
      p2 = p0(code) ; s2 = sz(code)
      do l=1,sz(4) ; do k=1,sz(3) ; do j=1,sz(2) ; do i=1,sz(1)
        a22(j,i,k,l) = a2e(i,j,k,l)
      end do ; end do ; end do ; end do
!   end function swap_jikl_real
