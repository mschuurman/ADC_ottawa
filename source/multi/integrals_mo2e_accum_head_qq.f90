      real(xrk), intent(in)       :: a2e(:,:,:,:)  ! spinless 2-e integrals over AOs
      real(xrk), intent(in)    :: mol(:)        ! spinless-AO coefficients for the fourth-index orbital
      real(xrk), intent(inout) :: buf_l(:,:,:)  ! Spin-block of partially-transformed integrals
