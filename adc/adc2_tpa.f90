module adc2tpamod

  use channels

contains

!#######################################################################

  subroutine adc2_tpa(gam)

    use constants
    use parameters
    use adc2common
    use fspace
    use misc
    use guessvecs
    use mp2
    use targetmatching
    
    implicit none

    integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
    integer                              :: i,ndim,ndims,ndimsf,&
                                            nout,ndimf,ndimd,noutf,&
                                            itmp
    integer*8                            :: noffd,noffdf
    real(d)                              :: e_init,e0,time
    real(d), dimension(:), allocatable   :: ener,mtm,mtmf,tmvec,osc_str,&
                                            vec_init,travec
    real(d), dimension(:,:), allocatable :: rvec,travec2
    type(gam_structure)                  :: gam


    print*,"FINISH WRITING THE TPA MODULE!"
    STOP
    
    return
    
  end subroutine adc2_tpa
    
!#######################################################################
  
end module adc2tpamod
