!   subroutine g_matrix_block_real(nao,bi,p0,sz,mxsz,a2e,rho,glocal)
      integer(ik), intent(in)    :: nao            ! Number of spin-less AOs
      integer(ik), intent(in)    :: bi(:)          ! Integral block index
      integer(ik), intent(in)    :: p0(:)          ! Initial position of the integral block
      integer(ik), intent(in)    :: sz(:)          ! Size of the integral block
      integer(ik), intent(in)    :: mxsz           ! Largest size of an integral block
!     real(rk), intent(in)       :: a2e(:,:,:,:)   ! Integrals block, for the canonical order of block indices
!                                                  ! (i<j) < (k<l)
!     complex(rk), intent(in)    :: rho   (:,:)    ! Shared density matrix
!     complex(rk), intent(inout) :: glocal(:,:)    ! Pre-thread copy of gmat; safe to update 
      !
      integer(ik)     :: p2(4)                     ! Same as p0, after index swapping 
      integer(ik)     :: s2(4)                     ! Same as sz, after index swapping
      real(kind(a2e)) :: a22(mxsz,mxsz,mxsz,mxsz)  ! Buffer for index-permuted integrals
      !
      !  This integral block may appear several times in G matrix construction
      !
                                             call accumulate_g(nao,p0,sz,mxsz,a2e,rho,glocal)
      if (swap_jikl(bi,p0,sz,p2,s2,a2e,a22)) call accumulate_g(nao,p2,s2,mxsz,a22,rho,glocal) 
      if (swap_ijlk(bi,p0,sz,p2,s2,a2e,a22)) call accumulate_g(nao,p2,s2,mxsz,a22,rho,glocal) 
      if (swap_jilk(bi,p0,sz,p2,s2,a2e,a22)) call accumulate_g(nao,p2,s2,mxsz,a22,rho,glocal) 
      if (swap_lkij(bi,p0,sz,p2,s2,a2e,a22)) call accumulate_g(nao,p2,s2,mxsz,a22,rho,glocal) 
      if (swap_klij(bi,p0,sz,p2,s2,a2e,a22)) call accumulate_g(nao,p2,s2,mxsz,a22,rho,glocal) 
      if (swap_klji(bi,p0,sz,p2,s2,a2e,a22)) call accumulate_g(nao,p2,s2,mxsz,a22,rho,glocal) 
      if (swap_lkji(bi,p0,sz,p2,s2,a2e,a22)) call accumulate_g(nao,p2,s2,mxsz,a22,rho,glocal) 
!   end subroutine g_matrix_block_real
