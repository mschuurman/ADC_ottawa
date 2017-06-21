      real(rk), intent(in)       :: a2e(:,:,:,:)  ! spinless 2-e integrals over AOs
      complex(rk), intent(in)    :: mol(:)        ! spinless-AO coefficients for the fourth-index orbital
      complex(rk), intent(inout) :: buf_l(:,:,:)  ! Spin-block of partially-transformed integrals
