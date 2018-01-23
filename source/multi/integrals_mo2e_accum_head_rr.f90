      real(rk), intent(in)       :: a2e(:,:,:,:)  ! spinless 2-e integrals over AOs
      real(rk), intent(in)    :: mol(:)        ! spinless-AO coefficients for the fourth-index orbital
      real(rk), intent(inout) :: buf_l(:,:,:)  ! Spin-block of partially-transformed integrals
