!
!  Cray SCIPORT/SCILIB interface for 2D and 3D complex to complex transforms
!
!  On input, we take uniformly spaced grids. On output, we
!  spit out frequency components (as you would expect from FFT ;-).
!  For each index, the frequencies increase monotonically from
!  -Pi to Pi. The zero-frequency component is found at element
!  1+(N-1)/2, where N is the number of points along a given dimension
!
 module fftw
   use accuracy
   implicit none

   private
!  public fftw_3d, fftw_2d
   public fftw_3d

   interface fftw_3d
     module procedure cfftw_3d
     module procedure zfftw_3d
   end interface ! fftw_3d

!  interface fftw_2d
!    module procedure cfftw_2d
!    module procedure zfftw_2d
!  end interface ! fftw_2d

   integer(ik), parameter :: verbose = 3

   contains

   subroutine cfftw_3d(n1_,n2_,n3_,a,b,invert)
     integer(ik), intent(in) :: n1_, n2_, n3_
     complex*8, intent(in)   :: a(n1_,n2_,n3_)
     complex*8, intent(out)  :: b(n1_,n2_,n3_)
     logical, intent(in)     :: invert
     !
     integer*4           :: isign, n1, n2, n3, i1x, i2x, i3x, &
                            i1y, i2y, i3y, ntable, nwork
     real*4              :: scale
     real*4, allocatable :: table(:), work(:)
     !
     n1 = n1_ ; n2 = n2_ ; n3 = n3_
     scale = 1
     i1x = 1 ; i2x = n1 ; i3x = n1*n2
     i1y = 1 ; i2y = n1 ; i3y = n1*n2
     ntable = 195 + 2*(n1+n2+n3)
     nwork  = 6*n1*n2*n3
     !
     allocate (table(ntable), work(nwork))
     !
     if (verbose>=1) then
       write (out,"('Activated cfft3d. n1= ',i8,' n2= ',i8,' n3= ',i8,' invert= ',l)") &
              n1, n2, n3, invert
     end if
     !
     !  Build up the tables
     !
     isign = 0
     call  cfft3d (isign,n1,n2,n3,scale,a,i1x,i2x,i3x,b,i1y,i2y,i3y,table,ntable,work,nwork)
     !
     isign = -1 ! Forward
     if (invert) isign = +1 ! Backward
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       b = cshift(b,+(n1-1)/2,dim=1)
       b = cshift(b,+(n2-1)/2,dim=2)
       b = cshift(b,+(n3-1)/2,dim=3)
     end if
     !
     call  cfft3d (isign,n1,n2,n3,scale,a,i1x,i2x,i3x,b,i1y,i2y,i3y,table,ntable,work,nwork)
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       b = cshift(b,-(n1-1)/2,dim=1)
       b = cshift(b,-(n2-1)/2,dim=2)
       b = cshift(b,-(n3-1)/2,dim=3)
     end if
     !
     deallocate (table,work)
   end subroutine cfftw_3d

   subroutine zfftw_3d(n1_,n2_,n3_,a,b,invert)
     integer(ik), intent(in) :: n1_, n2_, n3_
     complex*16, intent(in)  :: a(n1_,n2_,n3_)
     complex*16, intent(out) :: b(n1_,n2_,n3_)
     logical, intent(in)     :: invert
     !
     integer*4           :: isign, n1, n2, n3, i1x, i2x, i3x, &
                            i1y, i2y, i3y, ntable, nwork
     real*8              :: scale
     real*8, allocatable :: table(:), work(:)
     !
     n1 = n1_ ; n2 = n2_ ; n3 = n3_
     scale = 1
     i1x = 1 ; i2x = n1 ; i3x = n1*n2
     i1y = 1 ; i2y = n1 ; i3y = n1*n2
     ntable = 195 + 2*(n1+n2+n3)
     nwork  = 6*n1*n2*n3
     !
     allocate (table(ntable), work(nwork))
     !
     if (verbose>=1) then
       write (out,"('Activated zfft3d. n1= ',i8,' n2= ',i8,' n3= ',i8,' invert= ',l)") &
              n1, n2, n3, invert
     end if
     !
     isign = 0
     call  zfft3d (isign,n1,n2,n3,scale,a,i1x,i2x,i3x,b,i1y,i2y,i3y,table,ntable,work,nwork)
     !
     isign = -1 ! Forward
     if (invert) isign = +1 ! Backward
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       b = cshift(b,+(n1-1)/2,dim=1)
       b = cshift(b,+(n2-1)/2,dim=2)
       b = cshift(b,+(n3-1)/2,dim=3)
     end if
     !
     call  zfft3d (isign,n1,n2,n3,scale,a,i1x,i2x,i3x,b,i1y,i2y,i3y,table,ntable,work,nwork)
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       b = cshift(b,-(n1-1)/2,dim=1)
       b = cshift(b,-(n2-1)/2,dim=2)
       b = cshift(b,-(n3-1)/2,dim=3)
     end if
     !
     deallocate (table,work)
   end subroutine zfftw_3d

!  subroutine cfftw_2d(n1,n2,a,b,invert)
!    integer(ik), intent(in) :: n1, n2
!    complex, intent(in)     :: a(n1,n2) ! Due to the idiotic way FFTW 3 defines plans,
!    complex, intent(out)    :: b(n1,n2) ! these -can't- use Fortran-90 assumed shape
!    logical, intent(in)     :: invert
!    !
!    integer(hik) :: plan        ! Large enough to store a pointer
!    integer(ik)  :: direction
!    !
!    direction = FFTW_FORWARD
!    if (invert) direction = FFTW_BACKWARD
!    !
!    !  For the inverse transform, shift data into the order expected by FFTW
!    !
!    if (invert) then
!      b = cshift(b,+(n2-1)/2,dim=2)
!      b = cshift(b,+(n1-1)/2,dim=1)
!    end if
!    !
!*ps call sfftw_plan_dft_2d(plan, n1, n2, a, b, direction, FFTW_ESTIMATE) 
!*ps call sfftw_execute(plan) 
!*ps call sfftw_destroy_plan(plan)
!    !
!    if (.not.invert) then
!      !
!      !  FFTW returns data in the [0:2pi] range. We actually want to have
!      !  [-pi:pi] (it looks nicer!) so we shift it around.
!      !
!      b = cshift(b,-(n1-1)/2,dim=1)
!      b = cshift(b,-(n2-1)/2,dim=2)
!    end if
!  end subroutine cfftw_2d

!  subroutine zfftw_2d(n1,n2,a,b,invert)
!    integer(ik), intent(in)     :: n1, n2
!    double complex, intent(in)  :: a(n1,n2) ! Due to the idiotic way FFTW 3 defines plans,
!    double complex, intent(out) :: b(n1,n2) ! these -can't- use Fortran-90 assumed shape
!    logical, intent(in)         :: invert
!    !
!    integer(hik) :: plan        ! Large enough to store a pointer
!    integer(ik)  :: direction
!    !
!    direction = FFTW_FORWARD
!    if (invert) direction = FFTW_BACKWARD
!    !
!    !  For the inverse transform, shift data into the order expected by FFTW
!    !
!    if (invert) then
!      b = cshift(b,+(n1-1)/2,dim=1)
!      b = cshift(b,+(n2-1)/2,dim=2)
!    end if
!    !
!*ps call dfftw_plan_dft_2d(plan, n1, n2, a, b, direction, FFTW_ESTIMATE) 
!*ps call dfftw_execute(plan) 
!*ps call dfftw_destroy_plan(plan)
!    !
!    if (.not.invert) then
!      !
!      !  FFTW returns data in the [0:2pi] range. We actually want to have
!      !  [-pi:pi] (it looks nicer!) so we shift it around.
!      !
!      b = cshift(b,-(n1-1)/2,dim=1)
!      b = cshift(b,-(n2-1)/2,dim=2)
!    end if
!  end subroutine zfftw_2d

 end module fftw
