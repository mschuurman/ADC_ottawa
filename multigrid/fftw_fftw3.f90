!
!  WARNING: PRIOR TO MAY 12, 2008 ALL OUT-OF-PLACE INVERSE TRANSFORMS
!           WERE BROKEN. FORWARD AND IN-PLACE TRANSFORMS WERE OK.
!
!  FFTW 3 interface for 1D, 2D and 3D complex to complex transforms
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
   public fftw_3d, fftw_2d, fftw_1d, fftw_activate_threads

   interface fftw_3d
     module procedure cfftw_3d
     module procedure zfftw_3d
     module procedure cfftw_3d_inplace
     module procedure zfftw_3d_inplace
   end interface ! fftw_3d

   interface fftw_2d
     module procedure cfftw_2d
     module procedure zfftw_2d
     module procedure cfftw_2d_inplace
     module procedure zfftw_2d_inplace
   end interface ! fftw_2d

   interface fftw_1d
     module procedure cfftw_1d_inplace
     module procedure zfftw_1d_inplace
   end interface ! fftw_1d

   integer(ik), parameter :: FFTW_FORWARD        = -1
   integer(ik), parameter :: FFTW_BACKWARD       = +1
   integer(ik), parameter :: FFTW_MEASURE        =  0
   integer(ik), parameter :: FFTW_ESTIMATE       = 64
   integer(ik), parameter :: FFTW_DESTROY_INPUT  =  1
   integer(ik), parameter :: FFTW_PRESERVE_INPUT = 16

   integer(hik), save  :: threshold = -1_hik ! Number of elements needed to make threading worthwhile
   integer, save       :: nthreads  = -1     ! Number of threads available; must be of default kind

   contains
   !
   !  This subroutine is a doubly-edged sword: for small FFT sizes, it may
   !  actually decrease performance. Unfortunately, FFTW provides only an
   !  all-or-nothing control.
   !
   subroutine fftw_activate_threads(thread_threshold)
     !$ use OMP_LIB
     integer(hik), intent(in) :: thread_threshold ! Number of elements in array for threads
                                                  ! to be useful; negative values will disable
                                                  ! threading altogether
     !$ integer(ik) :: iret
     !$ !
     !$ !
     !$ threshold = thread_threshold
     !$ if (threshold<0) then
     !$   write (out,"('Will use a single thread for FFT')")
     !$ else
     !$   call sfftw_init_threads(iret)
     !$   if (iret==0) then
     !$     write (out,"('FFTW3 threads initialization (single) faile with error code ',i8)") iret
     !$     stop 'fftw_fftw3: FFTW threads initialization failed (single) !'
     !$   end if
     !$   !
     !$   call dfftw_init_threads(iret)
     !$   if (iret==0) then
     !$     write (out,"('FFTW3 threads initialization (double) faile with error code ',i8)") iret
     !$     stop 'fftw_fftw3: FFTW threads initialization failed (double) !'
     !$   end if
     !$   nthreads = omp_get_max_threads()
     !$   write (out,"('Will use up to ',i0,' threads for FFT')") nthreads
     !$   write (out,"('Parallel threshold is ',i0,' elements')") threshold
     !$ end if
   end subroutine fftw_activate_threads

   subroutine activate_threads(double,elements)
     logical, intent(in)      :: double   ! .True. if calling routine is in double precision
     integer(hik), intent(in) :: elements ! Size of the FFT
     !
     !$ if (threshold<0 .or. nthreads<=0) return ! Threading inactive; do nothing
     !$ if (elements<threshold) then
     !$   if (     double) call dfftw_plan_with_nthreads(1)
     !$   if (.not.double) call sfftw_plan_with_nthreads(1)
     !$ else
     !$   if (     double) call dfftw_plan_with_nthreads(nthreads)
     !$   if (.not.double) call sfftw_plan_with_nthreads(nthreads)
     !$ end if
   end subroutine activate_threads

   subroutine cfftw_3d(n1,n2,n3,a,b,invert)
     integer(ik), intent(in)     :: n1, n2, n3
     complex(srk), intent(inout) :: a(n1,n2,n3) ! Due to the idiotic way FFTW 3 defines plans,
     complex(srk), intent(out)   :: b(n1,n2,n3) ! these -can't- use Fortran-90 assumed shape
                                            ! a is left unchanged, but has to be made inout
     logical, intent(in)         :: invert
     !
     integer(hik) :: plan        ! Large enough to store a pointer
     integer(ik)  :: direction
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       a = cshift(a,+(n3-1)/2,dim=3)
       a = cshift(a,+(n2-1)/2,dim=2)
       a = cshift(a,+(n1-1)/2,dim=1)
     end if
     !
     call activate_threads(double=.false.,elements=product(int((/n1,n2,n3/),kind=hik)))
     call sfftw_plan_dft_3d(plan, n1, n2, n3, a, b, direction, FFTW_ESTIMATE+FFTW_PRESERVE_INPUT) 
     call sfftw_execute(plan) 
     call sfftw_destroy_plan(plan)
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       b = cshift(b,-(n1-1)/2,dim=1)
       b = cshift(b,-(n2-1)/2,dim=2)
       b = cshift(b,-(n3-1)/2,dim=3)
     else
       !
       !  Restore the input array to its original state
       !
       a = cshift(a,-(n1-1)/2,dim=1)
       a = cshift(a,-(n2-1)/2,dim=2)
       a = cshift(a,-(n3-1)/2,dim=3)
     end if
   end subroutine cfftw_3d

   subroutine cfftw_3d_inplace(n1,n2,n3,a,invert)
     integer(ik), intent(in)     :: n1, n2, n3
     complex(srk), intent(inout) :: a(n1,n2,n3) ! Due to the idiotic way FFTW 3 defines plans,
                                                ! a -can't- use Fortran-90 assumed shape
     logical(srk), intent(in)    :: invert
     !
     integer(hik) :: plan        ! Large enough to store a pointer
     integer(ik)  :: direction
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       a = cshift(a,+(n3-1)/2,dim=3)
       a = cshift(a,+(n2-1)/2,dim=2)
       a = cshift(a,+(n1-1)/2,dim=1)
     end if
     !
     call activate_threads(double=.false.,elements=product(int((/n1,n2,n3/),kind=hik)))
     call sfftw_plan_dft_3d(plan, n1, n2, n3, a, a, direction, FFTW_ESTIMATE) 
     call sfftw_execute(plan) 
     call sfftw_destroy_plan(plan)
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       a = cshift(a,-(n1-1)/2,dim=1)
       a = cshift(a,-(n2-1)/2,dim=2)
       a = cshift(a,-(n3-1)/2,dim=3)
     end if
   end subroutine cfftw_3d_inplace

   subroutine zfftw_3d(n1,n2,n3,a,b,invert)
     integer(ik), intent(in)     :: n1, n2, n3
     complex(drk), intent(inout) :: a(n1,n2,n3) ! Due to the idiotic way FFTW 3 defines plans,
     complex(drk), intent(out)   :: b(n1,n2,n3) ! these -can't- use Fortran-90 assumed shape
                                                     ! a is left unchanged, but has to be made inout
     logical, intent(in)              :: invert
     !
     integer(hik) :: plan        ! Large enough to store a pointer
     integer(ik)  :: direction
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       a = cshift(a,+(n1-1)/2,dim=1)
       a = cshift(a,+(n2-1)/2,dim=2)
       a = cshift(a,+(n3-1)/2,dim=3)
     end if
     !
     call activate_threads(double=.true.,elements=product(int((/n1,n2,n3/),kind=hik)))
     call dfftw_plan_dft_3d(plan, n1, n2, n3, a, b, direction, FFTW_ESTIMATE+FFTW_PRESERVE_INPUT) 
     call dfftw_execute(plan) 
     call dfftw_destroy_plan(plan)
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       b = cshift(b,-(n1-1)/2,dim=1)
       b = cshift(b,-(n2-1)/2,dim=2)
       b = cshift(b,-(n3-1)/2,dim=3)
     else
       !
       !  Restore the input array
       !
       a = cshift(a,-(n3-1)/2,dim=3)
       a = cshift(a,-(n2-1)/2,dim=2)
       a = cshift(a,-(n1-1)/2,dim=1)
     end if
   end subroutine zfftw_3d

   subroutine zfftw_3d_inplace(n1,n2,n3,a,invert)
     integer(ik), intent(in)     :: n1, n2, n3
     complex(drk), intent(inout) :: a(n1,n2,n3) ! Due to the idiotic way FFTW 3 defines plans,
                                                     ! these -can't- use Fortran-90 assumed shape
     logical, intent(in)         :: invert
     !
     integer(hik) :: plan        ! Large enough to store a pointer
     integer(ik)  :: direction
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       a = cshift(a,+(n1-1)/2,dim=1)
       a = cshift(a,+(n2-1)/2,dim=2)
       a = cshift(a,+(n3-1)/2,dim=3)
     end if
     !
     call activate_threads(double=.true.,elements=product(int((/n1,n2,n3/),kind=hik)))
     call dfftw_plan_dft_3d(plan, n1, n2, n3, a, a, direction, FFTW_ESTIMATE) 
     call dfftw_execute(plan) 
     call dfftw_destroy_plan(plan)
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       a = cshift(a,-(n1-1)/2,dim=1)
       a = cshift(a,-(n2-1)/2,dim=2)
       a = cshift(a,-(n3-1)/2,dim=3)
     end if
   end subroutine zfftw_3d_inplace

   subroutine cfftw_2d(n1,n2,a,b,invert,good_plan)
     integer(ik), intent(in)       :: n1, n2
     complex(srk), intent(inout)   :: a(n1,n2) ! Due to the idiotic way FFTW 3 defines plans,
     complex(srk), intent(out)     :: b(n1,n2) ! these -can't- use Fortran-90 assumed shape
                                               ! a is left unchanged, but has to be made inout
     logical, intent(in)           :: invert
     logical, intent(in), optional :: good_plan ! .True. requests a better-quality FFT plan
     !
     integer(hik) :: plan        ! Large enough to store a pointer
     integer(ik)  :: direction
     integer(ik)  :: flags
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       a = cshift(a,+(n2-1)/2,dim=2)
       a = cshift(a,+(n1-1)/2,dim=1)
     end if
     !
     flags = FFTW_PRESERVE_INPUT
     if (present(good_plan)) then
       if (good_plan) then
         flags = flags + FFTW_MEASURE
       else
         flags = flags + FFTW_ESTIMATE
       end if
     else
       flags = flags + FFTW_ESTIMATE
     end if
     !
     call activate_threads(double=.false.,elements=product(int((/n1,n2/),kind=hik)))
     call sfftw_plan_dft_2d(plan, n1, n2, a, b, direction, flags)
     call sfftw_execute(plan) 
     call sfftw_destroy_plan(plan)
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       b = cshift(b,-(n1-1)/2,dim=1)
       b = cshift(b,-(n2-1)/2,dim=2)
     else
       !
       !  Restore the input array
       !
       a = cshift(a,-(n1-1)/2,dim=1)
       a = cshift(a,-(n2-1)/2,dim=2)
     end if
   end subroutine cfftw_2d

   subroutine cfftw_2d_inplace(n1,n2,a,invert,good_plan)
     integer(ik), intent(in)       :: n1, n2
     complex(srk), intent(inout)   :: a(n1,n2)  ! Due to the idiotic way FFTW 3 defines plans ...
     logical, intent(in)           :: invert
     logical, intent(in), optional :: good_plan ! .True. requests a better-quality FFT plan
     !
     integer(hik) :: plan        ! Large enough to store a pointer
     integer(ik)  :: direction
     integer(ik)  :: flags
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       a = cshift(a,+(n2-1)/2,dim=2)
       a = cshift(a,+(n1-1)/2,dim=1)
     end if
     !
     flags = FFTW_PRESERVE_INPUT
     if (present(good_plan)) then
       if (good_plan) then
         flags = flags + FFTW_MEASURE
       else
         flags = flags + FFTW_ESTIMATE
       end if
     else
       flags = flags + FFTW_ESTIMATE
     end if
     !
     call activate_threads(double=.false.,elements=product(int((/n1,n2/),kind=hik)))
     call sfftw_plan_dft_2d(plan, n1, n2, a, a, direction, flags)
     call sfftw_execute(plan) 
     call sfftw_destroy_plan(plan)
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       a = cshift(a,-(n1-1)/2,dim=1)
       a = cshift(a,-(n2-1)/2,dim=2)
     end if
   end subroutine cfftw_2d_inplace

   subroutine zfftw_2d(n1,n2,a,b,invert,good_plan)
     integer(ik), intent(in)       :: n1, n2
     complex(drk), intent(inout)   :: a(n1,n2) ! Due to the idiotic way FFTW 3 defines plans,
     complex(drk), intent(out)     :: b(n1,n2) ! these -can't- use Fortran-90 assumed shape
                                                  ! a is left unchanged, but has to be made inout
     logical, intent(in)           :: invert
     logical, intent(in), optional :: good_plan ! .True. requests a better-quality FFT plan
     !
     integer(hik) :: plan        ! Large enough to store a pointer
     integer(ik)  :: direction
     integer(ik)  :: flags
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       a = cshift(a,+(n1-1)/2,dim=1)
       a = cshift(a,+(n2-1)/2,dim=2)
     end if
     !
     flags = FFTW_PRESERVE_INPUT
     if (present(good_plan)) then
       if (good_plan) then
         flags = flags + FFTW_MEASURE
       else
         flags = flags + FFTW_ESTIMATE
       end if
     else
       flags = flags + FFTW_ESTIMATE
     end if
     !
     call activate_threads(double=.true.,elements=product(int((/n1,n2/),kind=hik)))
     call dfftw_plan_dft_2d(plan, n1, n2, a, b, direction, flags)
     call dfftw_execute(plan) 
     call dfftw_destroy_plan(plan)
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       b = cshift(b,-(n1-1)/2,dim=1)
       b = cshift(b,-(n2-1)/2,dim=2)
     else
       !
       !  Restore the input array
       !
       a = cshift(a,-(n2-1)/2,dim=2)
       a = cshift(a,-(n1-1)/2,dim=1)
     end if
   end subroutine zfftw_2d

   subroutine zfftw_2d_inplace(n1,n2,a,invert,good_plan)
     integer(ik), intent(in)       :: n1, n2
     complex(drk), intent(inout)   :: a(n1,n2) 
     logical, intent(in)           :: invert
     logical, intent(in), optional :: good_plan ! .True. requests a better-quality FFT plan
     !
     integer(hik) :: plan        ! Large enough to store a pointer
     integer(ik)  :: direction
     integer(ik)  :: flags
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     ! write (out,"('zfftw_2d_inplace: ',i5,'x',i5,' shifts = ',i5,1x,i5)") &
     !        n1, n2, (n1-1)/2, (n2-1)/2
     ! write (out,"('zfftw_2d_inplace: invert = ',l3)") invert
     !
     !  For the inverse transform, shift data into the order expected by FFTW
     !
     if (invert) then
       a = cshift(a,+(n1-1)/2,dim=1)
       a = cshift(a,+(n2-1)/2,dim=2)
     end if
     !
     ! write (out,"('FFTW input array:')")
     ! do irow=1,n1
     !  write (out,"(1x,i5,(t8,5(f12.5,1x,f12.5,2x)))") irow, a(irow,1:n2)
     ! end do
     !
     flags = FFTW_PRESERVE_INPUT
     if (present(good_plan)) then
       if (good_plan) then
         flags = flags + FFTW_MEASURE
       else
         flags = flags + FFTW_ESTIMATE
       end if
     else
       flags = flags + FFTW_ESTIMATE
     end if
     !
     call activate_threads(double=.true.,elements=product(int((/n1,n2/),kind=hik)))
     call dfftw_plan_dft_2d(plan, n1, n2, a, a, direction, flags)
     call dfftw_execute(plan) 
     call dfftw_destroy_plan(plan)
     !
     ! write (out,"('FFTW return array:')")
     ! do irow=1,n1
     !  write (out,"(1x,i5,(t8,5(f12.5,1x,f12.5,2x)))") irow, a(irow,1:n2)
     ! end do
     !
     if (.not.invert) then
       !
       !  FFTW returns data in the [0:2pi] range. We actually want to have
       !  [-pi:pi] (it looks nicer!) so we shift it around.
       !
       a = cshift(a,-(n1-1)/2,dim=1)
       a = cshift(a,-(n2-1)/2,dim=2)
       !
       ! write (out,"('Rearranged array:')")
       ! do irow=1,n1
       !  write (out,"(1x,i5,(t8,5(f12.5,1x,f12.5,2x)))") irow, a(irow,1:n2)
       ! end do
     end if
   end subroutine zfftw_2d_inplace

   subroutine cfftw_1d_inplace(n1,n2,a,invert)
     integer(ik), intent(in)      :: n1, n2
     complex(srk), intent(inout)  :: a(n1,n2) ! Due to the idiotic way FFTW 3 defines plans ...
     logical, intent(in)          :: invert
     !
     integer(hik) :: plan         ! Large enough to store a pointer
     integer(ik)  :: direction, i !
     complex(srk) :: buf(n1)      ! Local buffer for 1D transforms
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     !  FFTW_MEASURE could take a long time; to be of any benefit, the array must be
     !  large enough to justify the cost
     !
     call activate_threads(double=.false.,elements=int(n1,kind=hik))
     if (n2>100) then
       call sfftw_plan_dft_1d(plan, n1, buf, buf, direction, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
     else
       call sfftw_plan_dft_1d(plan, n1, buf, buf, direction, FFTW_ESTIMATE+FFTW_DESTROY_INPUT) 
     end if
     column_loop: do i=1,n2
       if (invert) then
         buf = cshift(a(:,i),+(n1-1)/2)
       else
         buf = a(:,i)
       end if
       !
       call sfftw_execute(plan) 
       !
       if (invert) then
         a(:,i) = buf
       else
         a(:,i) = cshift(buf,-(n1-1)/2)
       end if
     end do column_loop
     call sfftw_destroy_plan(plan)
     !
   end subroutine cfftw_1d_inplace

   subroutine zfftw_1d_inplace(n1,n2,a,invert)
     integer(ik), intent(in)     :: n1, n2
     complex(drk), intent(inout) :: a(n1,n2) ! Due to the idiotic way FFTW 3 defines plans ...
     logical, intent(in)         :: invert
     !
     integer(hik) :: plan         ! Large enough to store a pointer
     integer(ik)  :: direction, i !
     complex(drk) :: buf(n1)      ! Local buffer for 1D transforms
     !
     direction = FFTW_FORWARD
     if (invert) direction = FFTW_BACKWARD
     !
     call activate_threads(double=.true.,elements=int(n1,kind=hik))
     if (n2>100) then
       call dfftw_plan_dft_1d(plan, n1, buf, buf, direction, FFTW_MEASURE+FFTW_DESTROY_INPUT) 
     else
       call dfftw_plan_dft_1d(plan, n1, buf, buf, direction, FFTW_ESTIMATE+FFTW_DESTROY_INPUT) 
     end if
     column_loop: do i=1,n2
       if (invert) then
         buf = cshift(a(:,i),+(n1-1)/2)  ! Circular shift to the left
       else
         buf = a(:,i)
       end if
       !
       call dfftw_execute(plan) 
       !
       if (invert) then
         a(:,i) = buf
       else
         a(:,i) = cshift(buf,-(n1-1)/2)  ! Circular shift to the right
       end if
     end do column_loop
     call dfftw_destroy_plan(plan)
     !
   end subroutine zfftw_1d_inplace

 end module fftw
