!
!  Evaluation of two-particle interaction potentials (aka 
!  general convolutions) on a grid.
!
!  The implementation relies on the following assumptions:
!
!  1. The interaction potential is non-singular and square-
!     integratable (or FFT won't work)
!  2. The interaction potential decays monotonically and faster
!     than 1/r once it extends outside of the 0,0,0 unit cell
!  3. The grid is rectangular and uniform, with periodic
!     boundary conditions imposed
!
!  Bugs: This module uses complex-to-complex FFT, even though
!        all quantities involved are real. This is done because
!        we already have C2C FFT, and the speed of the routines 
!        here is not critical.
!
module convolution
  use accuracy
  use timer
  use fftw
  implicit none
  private
  public ConvolutionT
  public convolution_initialize, convolution_destroy
  public convolution_rescale, convolution_uniform_shift, convolution_evaluate
  !
  type ConvolutionT
    private
    integer(ik)          :: verbose     ! Verbosity level for the convolution
    real(rk)             :: scale       ! Overall scaling factor for the potential
    integer(ik)          :: n(3)        ! Grid dimensions
    real(rk)             :: step(3)     ! Grid step
    complex(rk), pointer :: fft(:,:,:)  ! Fourier components of the convolution
    complex(rk), pointer :: scr(:,:,:)  ! Scratch area used for mean-field evaluation
  end type ! ConvolutionT
  !
  contains
  !
  subroutine convolution_initialize(cnv,v12,eps,step,npoints,scale,verbose)
    type(ConvolutionT), intent(out)   :: cnv        ! Convolution descriptor, to be filled out here
    real(rk), external                :: v12        ! Interaction potential, takes 1 real argument 
                                                    ! and returns a real
    real(rk), intent(in)              :: eps        ! Absolute cut-off point for v12
    real(rk), intent(in)              :: step(3)    ! Grid step for each direction
    integer(ik), intent(in)           :: npoints(3) ! Number of points along each dimension
    real(rk), intent(in), optional    :: scale      ! Scale of the interaction potential
    integer(ik), intent(in), optional :: verbose    ! Verbosity level
    !
    real(rk)                 :: r_cut               ! Distance at which v12 vanishes
    real(rk)                 :: v_cut               ! Value of v12 at the vanishing point
    integer(ik)              :: n12(3)              ! Dimensions of the extended grid
    integer(ik)              :: alloc
    complex(rk), allocatable :: v12_grid(:,:,:)     ! Extended grid, used to include periodic summation
                                                    ! with v12
    !
    call TimerStart('Convolution initialization')
    !
    !  Step 0: Fill out the descriptor and alocate space used in convolution
    !
    cnv%verbose = 0
    cnv%scale   = 1.0_rk
    cnv%n       = npoints
    cnv%step    = step
    if (present(verbose)) cnv%verbose = verbose
    if (present(scale  )) cnv%scale   = scale
    allocate (cnv%fft(npoints(1),npoints(2),npoints(3)), &
              cnv%scr(npoints(1),npoints(2),npoints(3)),stat=alloc)
    if (alloc/=0) then
      write (out,"('convolution_initialize: Error ',i8,' allocating FFT image. Bounds = ',3i8)") &
             alloc, npoints
      stop 'convolution%convolution_initialize - no memory for FFT image & scratch'
    end if
    !
    !  Step 1: figure out how far we need to go before the potential vanishes.
    !
    r_cut = v12_vanishing_point(v12,cnv%scale,eps,min(10._rk,0.5_rk*maxval(step*npoints)))
    if (cnv%verbose>=0) then
      v_cut = abs(v12(r_cut)) ! Try to avoid a possibility of recursive I/O
      write (out,"(/' Convolution function at ',g14.6,' bohr is below ',g12.4,' hartree')") &
             r_cut, v_cut
    end if
    !
    !  Step 2: create temporary grid for v12 sampling in real space.
    !          Grid spacing must match the simulation grid, so that
    !          we are only allowed to grow the grid by an integer factor.
    !
    n12 = int(r_cut / step,kind=ik)
    n12 = npoints * ((n12 + npoints - 1) / npoints)
    if (cnv%verbose>=0) then
      write (out,"(' Upper bounds for the auxiliary grid = ',3i8)") n12
      write (out,"(' Upper bounds for the physical grid  = ',3i8)") npoints
    end if
    allocate (v12_grid(n12(1),n12(2),n12(3)),stat=alloc)
    if (alloc/=0) then
      write (out,"('convolution_initialize: Error ',i8,' allocating auxiliary grid. Bounds = ',3i8)") &
             alloc, n12
      stop 'convolution%convolution_initialize - no memory for aux. grid'
    end if
    !
    !  Step 3: Populate and Fourier-transform the auxiliary grid
    !
    call populate_v12_grid(v12,cnv%scale,n12,step,v12_grid)
    call fftw_3d(n12(1),n12(2),n12(3),v12_grid,invert=.false.)
    !
    !  Step 4: Extract components of the transform which are commensurate
    !          with the physical grid.
    !
    call reap_v12_fft(n12,v12_grid,npoints,cnv%fft)
    !
    !  Step 5: Make sure that the resulting FFT image of the potential makes
    !          sense.
    !
    if (cnv%verbose>=0) then
      call check_extraction_sanity(cnv,v12)
    end if
    !
    !  Release the memory used by the extended grid, and we are done
    !
    deallocate (v12_grid,stat=alloc)
    if (alloc/=0) then
      write (out,"('convolution_initialize: Error ',i8,' releasing auxiliary grid.')") &
             alloc, n12
      stop 'convolution%convolution_initialize - deallocation error'
    end if
    !
    call TimerStop('Convolution initialization')
  end subroutine convolution_initialize
  !
  !  Change the overall scale of the existing convolution
  !
  subroutine convolution_rescale(cnv,scale)
    type(ConvolutionT), intent(inout) :: cnv    ! Convolution, must be already initialized
    real(rk), intent(in)              :: scale  ! Factor to apply to the convolution
    !
    cnv%fft   = cnv%fft * scale
    cnv%scale = cnv%scale * scale
  end subroutine convolution_rescale
  !
  !  Destroy allocatable fields in a convolution descriptor.
  !  
  subroutine convolution_destroy(cnv)
    type(ConvolutionT), intent(inout) :: cnv
    !
    integer(ik) :: alloc
    !
    deallocate (cnv%fft,cnv%scr,stat=alloc)
    if (alloc/=0) then
      write (out,"('convolution_destroy: error ',i8,' releasing data fields')") alloc
      stop 'convolution%convolution_destroy - deallocation error'
    end if
  end subroutine convolution_destroy
  !
  !  Return the integral of v12 over the entire space (uniform fluid),
  !  including the volume element.
  !
  function convolution_uniform_shift(cnv) result(v)
    type(ConvolutionT), intent(in) :: cnv ! Convolution descriptor
    real(rk)                       :: v
    !
    integer(ik) :: n0(3)
    !
    !  Return zero-frequency Fourier component of the transform.
    !
    n0 = 1+(cnv%n-1)/2
    v = product(cnv%step)*real(cnv%fft(n0(1),n0(2),n0(3)),kind=rk)
    !
  end function convolution_uniform_shift
  !
  !  Evaluate mean-field potential for the given density, including
  !  the grid volume element.
  !
  subroutine convolution_evaluate(cnv,dens,v)
    type(ConvolutionT), intent(inout) :: cnv  ! Convolution descriptor
    real(rk), intent(in)              :: dens(cnv%n(1),cnv%n(2),cnv%n(3)) ! Particle density
    real(rk), intent(out)             :: v   (cnv%n(1),cnv%n(2),cnv%n(3)) ! Mean-field potential
    !
    integer(ik) :: i1, i2, i3
    real(rk)    :: max_re, max_im
    !
    call TimerStart('Convolution evaluation')
    !
    !  Calculate FFT of the input density
    !
    cnv%scr = dens
    call fftw_3d(cnv%n(1),cnv%n(2),cnv%n(3),cnv%scr,invert=.false.)
    !
    !  Multiply by the image of the 2-particle interaction potential
    !
    cnv%scr = cnv%scr * cnv%fft
    !
    !  Transform back and normalize to the number of elements
    !
    call fftw_3d(cnv%n(1),cnv%n(2),cnv%n(3),cnv%scr,invert=.true.)
    cnv%scr = cnv%scr * (product(cnv%step)/product(cnv%n))
    !
    !  Copy the mean-field potential to the output, making sure it's real
    !  along the way
    !
    if (cnv%verbose>=3) then
      write (out,"(/3(1x,a4),1x,2(1x,a14))") &
             ' i1 ', ' i2 ', ' i3 ', ' Re(v) ', ' Im(v) '
    end if
    !
    max_re = 0
    max_im = 0
    !$omp parallel do private(i1,i2,i3) reduction(max:max_re,max_im)
    r3_loop: do i3=1,cnv%n(3)
      r2_loop: do i2=1,cnv%n(2)
        r1_loop: do i1=1,cnv%n(1)
          if (cnv%verbose>=3) then
            write (out,"(3(1x,i4),1x,2(1x,g14.7))") i1, i2, i3, cnv%scr(i1,i2,i3)
          end if
          max_re  = max(max_re,abs(real (cnv%scr(i1,i2,i3),kind=rk)))
          max_im  = max(max_im,abs(aimag(cnv%scr(i1,i2,i3))))
          v(i1,i2,i3) = real(cnv%scr(i1,i2,i3),kind=rk)
        end do r1_loop
      end do r2_loop
    end do r3_loop
    !
    if (max_im>=spacing(100*max_re)) then
      write (out,"('Mean-field potential is not real: Max(re)= ',g12.5,' Max(im)= ',g12.5)") &
             max_re, max_im
      stop 'convolution%convolution_evaluate - potential is not real enough'
    end if
    !
    call TimerStop('Convolution evaluation')
    !
  end subroutine convolution_evaluate
  !
  ! === Everything below is local, and should not be exported outside ===
  !
  function v12_vanishing_point(v12,scale,eps,r0) result (r_cut)
    real(rk), external   :: v12     ! Interaction potential
    real(rk), intent(in) :: scale   ! Overall scale
    real(rk), intent(in) :: eps     ! Cut-off point
    real(rk), intent(in) :: r0      ! Starting point for the search
    real(rk)             :: r_cut
    !
    real(rk), parameter :: step = 1.02_rk ! Increment in cut-off point
    real(rk)            :: v
    !
    r_cut = r0 / step
    grow_cutoff: do
      r_cut = r_cut * step
      v     = scale*v12(r_cut)
      if (abs(v)<eps) exit grow_cutoff
    end do grow_cutoff
  end function v12_vanishing_point
  !
  subroutine populate_v12_grid(v12,scale,n,step,v12_grid)
    real(rk), external       :: v12     ! Interaction potential
    real(rk), intent(in)     :: scale   ! Overall scale of the convolution
    integer(ik), intent(in)  :: n(3)    ! Aux. grid extent
    real(rk), intent(in)     :: step(3) ! Aux. grid spacing
    complex(rk), intent(out) :: v12_grid(n(1),n(2),n(3))
    !
    integer(ik) :: i1, i2, i3
    real(rk)    :: r1, r2, r3
    !
    !  The origin of the operator is at the index (1,1,1) -
    !  this matches our choice of the Fourier transform factors
    !
    !$omp parallel do private(i1,i2,i3,r1,r2,r3)
    g3_loop: do i3=1,n(3)
      r3 = min(i3-1,n(3)+1-i3) * step(3)
      g2_loop: do i2=1,n(2)
        r2 = min(i2-1,n(2)+1-i2) * step(2)
        g1_loop: do i1=1,n(1)
          r1 = min(i1-1,n(1)+1-i1) * step(1)
          v12_grid(i1,i2,i3) = scale*v12(sqrt(r1**2+r2**2+r3**2))
        end do g1_loop
      end do g2_loop
    end do g3_loop
    !$omp end parallel do
  end subroutine populate_v12_grid
  !
  !  Extracts components of V12 Fourier image, which are needed for
  !  evaluation of the potential on the physically periodic grid.
  !
  subroutine reap_v12_fft(n_aux,v_aux,n_phys,v_phys)
    integer(ik), intent(in)  :: n_aux (3)     ! Size of the auxiliary grid
    complex(rk), intent(in)  :: v_aux(:,:,:)  ! FFT result on commensurate containing grid
    integer(ik), intent(in)  :: n_phys(3)     ! Size of the physical grid
    complex(rk), intent(out) :: v_phys(:,:,:) ! Coefficients matching allowed physical oscillations
    !
    integer(ik) :: z_aux (3) ! Position of the zero momentum on the aux grid
    integer(ik) :: z_phys(3) ! ditto, physical grid
    integer(ik) :: scale (3) ! Ratio of aux/physical grid sizes
    integer(ik) :: a1, a2, a3, f1, f2, f3
    !
    z_aux  = 1 + (n_aux -1)/2
    z_phys = 1 + (n_phys-1)/2
    scale  = n_aux/n_phys
    !
    if (any(n_aux/=scale*n_phys)) then
      write (out,"('Auxiliary grid (',3i5,') is not commensurate to the physical grid (',3i5,')')") &
             n_aux, n_phys
      stop 'convolution%reap_v12_fft - incommensurate grids'
    end if
    !
    !$omp parallel do private(f1,f2,f3,a1,a2,a3)
    reap_3: do f3=1,n_phys(3)
      a3 = z_aux(3) + scale(3)*(f3-z_phys(3))
      reap_2: do f2=1,n_phys(2)
        a2 = z_aux(2) + scale(2)*(f2-z_phys(2))
        reap_1: do f1=1,n_phys(1)
          a1 = z_aux(1) + scale(1)*(f1-z_phys(1))
          v_phys(f1,f2,f3) = v_aux(a1,a2,a3)
         !write (out,"('[',3i5,'] = (',3i5,') = ',2g14.7)") &
         !       f1, f2, f3, a1, a2, a3, v_aux(a1,a2,a3)
        end do reap_1
      end do reap_2
    end do reap_3
    !$omp end parallel do
  end subroutine reap_v12_fft
  !
  !  Verify that the Fourier image of the interaction potential is
  !  at least somewhat sensible.
  !
  subroutine check_extraction_sanity(cnv,v12)
    type(ConvolutionT), intent(inout) :: cnv ! Convolution descriptor
    real(rk), external                :: v12 ! Interaction potential
    !
    integer(ik) :: i1, i2, i3
    real(rk)    :: r1, r2, r3, r12, v12_ref, vdiff
    real(rk)    :: max_re, max_im, max_vd
    !
    !  Begin by evaluating the inverse transform on the physical
    !  grid. The result must be: a) Purely real; and b) similar
    !  to the real-space interaction potential.
    !
    cnv%scr = cnv%fft
    call fftw_3d(cnv%n(1),cnv%n(2),cnv%n(3),cnv%scr,invert=.true.)
    cnv%scr = cnv%scr * (1._rk/product(cnv%n))
    !
    !  Accumulate some data on the inverse transform: the largest
    !  real part; the largest imaginary part; the max difference
    !  between the two potentials
    !
    if (cnv%verbose>=3) then
      write (out,"(/3(1x,a4),1x,3(1x,a10),1x,2(1x,a14))") &
             ' i1 ', ' i2 ', ' i3 ', '  x  ', '  y  ', '  z  ', &
             ' v_ref ', ' Re(v) ', ' Im(v) '
    end if
    !
    max_re = 0
    max_im = 0
    max_vd = 0
    !$omp parallel do private(i1,i2,i3,r1,r2,r3,r12,v12_ref,vdiff) &
    !$omp& reduction(max:max_re,max_im,max_vd)
    r3_loop: do i3=1,cnv%n(3)
      r3 = min(i3-1,cnv%n(3)+1-i3) * cnv%step(3)
      r2_loop: do i2=1,cnv%n(2)
        r2 = min(i2-1,cnv%n(2)+1-i2) * cnv%step(2)
        r1_loop: do i1=1,cnv%n(1)
          r1 = min(i1-1,cnv%n(1)+1-i1) * cnv%step(1)
          r12 = sqrt(r1**2+r2**2+r3**2)
          v12_ref = cnv%scale*v12(r12)
          if (cnv%verbose>=3) then
            write (out,"(3(1x,i4),1x,3(1x,f10.4),1x,g14.7,1x,2(1x,g14.7))") &
                   i1, i2, i3, r1, r2, r3, v12_ref, cnv%scr(i1,i2,i3)
          end if
          vdiff   = abs(v12_ref-cnv%scr(i1,i2,i3))
          max_re  = max(max_re,abs(real (cnv%scr(i1,i2,i3),kind=rk)))
          max_im  = max(max_im,abs(aimag(cnv%scr(i1,i2,i3))))
          max_vd  = max(max_vd,abs(v12_ref-real(cnv%scr(i1,i2,i3),kind=rk)))
        end do r1_loop
      end do r2_loop
    end do r3_loop
    !
    if (cnv%verbose>=2) then
      write (out,"(/' Back-transformation of the extended potential gives:')")
      write (out,"( '   Max(Re(V))       = ',g14.6)") max_re
      write (out,"( '   Max(Im(V))       = ',g14.6)") max_im
      write (out,"( '   Max(V_ext-V_min) = ',g14.6)") max_vd
      write (out,"()")
    end if
    !
    if (max_im>=spacing(100*max_re)) then
      write (out,"('Extended potential is not real: Max(re)= ',g12.5,' Max(im)= ',g12.5)") &
             max_re, max_im
      stop 'convolution%check_extraction_sanity - potential is not real enough'
    end if
    !
  end subroutine check_extraction_sanity
  !
end module convolution
