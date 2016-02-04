  module relaxmod

    use constants
    
    implicit none

    save

!-----------------------------------------------------------------------
! System definition
!-----------------------------------------------------------------------
    integer                            :: natm,ncoo
    real(d), dimension(:), allocatable :: xcoo

!-----------------------------------------------------------------------
! SCF parameters
!-----------------------------------------------------------------------
    integer :: scfiter

!-----------------------------------------------------------------------
! GAMESS Parameters
!-----------------------------------------------------------------------
    integer, parameter                          :: maxao=100
    integer, parameter                          :: maxprim=100
    integer                                     :: ncoo,difftype,contcent
    integer, dimension(5)                       :: ndiff
    integer, dimension(:), allocatable          :: naogms
    integer, dimension(:,:), allocatable        :: nprim,aotype,ilquant
    real(d), dimension(:,:,:), allocatable      :: aocoeff,aoexp
    real(d), dimension(:,:), allocatable        :: diffexp
    character(len=60)                           :: basname
    character(len=2), dimension(:), allocatable :: atlbl
    logical                                     :: lrungamess,ldiffcom

  end module relaxmod
