!
!  Subroutines used for integration of phase factors in eikonal
!  continuum wavefunctions (adiabatic approximation). For 
!  derivation of formulae see Mathematica notebook - 
!  "phase-interpolation-revised.nb" 
!
!  Calculation of eikonal phase corrections requires evaluation of
!  line integrals of the potential along the asymptotic direction
!  of the scattered wave travel (the k vector). For Cartesian product
!  grids, the problem would have been trivial if k were along one 
!  of the Cartesian directions. Since it is not, the easiest approach
!  seems to be the following:
!
!  0. Set up asymptotic boundary conditions at the edges of the 
!     integration volume
!  1. Interpolate the value of the integral in the base plane
!     (this is done using "Wu" matrix below).
!  2. Interpolate the integrand in the cube around the target
!     point, and calculate line integral from the base point
!     to the target point (this is done using "Wg" matrix below)
!
!  In practice, six special cases arise depending on the Cartesian
!  axis having the largest overlap with the k vector. These are
!  all handled in multigrid.f90 - the coefficients of the stencils
!  are defined identically for all six cases.
!
!  We choose to interpolate using products of Chebychev polynomials,
!  which is expected to give uniformly distributed errors. In practice,
!  the interpolation breaks down in the immediate vicinity of the
!  singularities, leading to characteristic artifacts in the phases
!  (the "bullet holes"). Since these artifacts appear in a small 
!  fraction of the grid points, and have wewll-behaved amplitudes,
!  the overall transition matrix elements are not affected.
!
module phase_interpolation
  use accuracy
  implicit none
  private
  public GetMatrixWg, GetMatrixWu
  
  contains
  !
  !  Functional interpolation using 2-nd order Chebyshev products - (3,3) grid
  !
  subroutine GetMatrixWg(rr,Wg)
    real(rk), intent(in)   :: rr(:)         ! Reduced displacement in the base plane - R(1:3)/gr%step(1:3)
                                            ! For some magical reason, declaring rr as rr(3) causes a
                                            ! segmentation violation in ifort 9.1 debugging compile. Go figure.
    real(rk), intent(out)  :: Wg(-1:1,-1:1) ! Interpolation matrix for the gradient value in the base plane
                                            ! Contraction with the values on grid should give the interpolated
                                            ! value at the intersection point.
    real(rk)    :: rx, ry ! Displacements
    !
    rx  = rr(1) 
    ry  = rr(2)
    !
    Wg(-1,-1) = (rx*ry)/4.0_rk - (rx**2*ry)/4.0_rk - (rx*ry**2)/4.0_rk + (rx**2*ry**2)/4.0_rk
    Wg(-1, 0) = -rx/2.0_rk + rx**2/2.0_rk + (rx*ry**2)/2.0_rk - (rx**2*ry**2)/2.0_rk
    Wg(-1, 1) = -(rx*ry)/4.0_rk + (rx**2*ry)/4.0_rk - (rx*ry**2)/4.0_rk + (rx**2*ry**2)/4.0_rk
    Wg( 0,-1) = -ry/2.0_rk + (rx**2*ry)/2.0_rk + ry**2/2.0_rk - (rx**2*ry**2)/2.0_rk
    Wg( 0, 0) = 1 - rx**2 - ry**2 + rx**2*ry**2
    Wg( 0, 1) = ry/2.0_rk - (rx**2*ry)/2.0_rk + ry**2/2.0_rk - (rx**2*ry**2)/2.0_rk
    Wg( 1,-1) = -(rx*ry)/4.0_rk - (rx**2*ry)/4.0_rk + (rx*ry**2)/4.0_rk + (rx**2*ry**2)/4.0_rk
    Wg( 1, 0) = rx/2.0_rk + rx**2/2.0_rk - (rx*ry**2)/2.0_rk - (rx**2*ry**2)/2.0_rk
    Wg( 1, 1) = (rx*ry)/4.0_rk + (rx**2*ry)/4.0_rk + (rx*ry**2)/4.0_rk + (rx**2*ry**2)/4.0_rk
  end subroutine GetMatrixWg
  !
  !  2-nd order Chebyshev interpolation in the (3,3,3) cube, followed by line integral
  !  from the intersection point to the origin.
  !
  subroutine GetMatrixWu(rr,Wu)
    real(rk), intent(in)   :: rr(3)                ! Reduced direction of the integration line:
                                                   ! R/gr%step. The integration limits are [-1:0],
                                                   ! so that the lower plane intersection point
                                                   ! is at rr(1:2).
    real(rk), intent(out)  :: Wu(-1:1,-1:1,-1:1)   ! Weight matrix for line integral from the lower
                                                   ! intersection point to the origin. Contraction
                                                   ! with the function values on grid should give
                                                   ! the value of the integral.
    !
    real(rk)     :: rx, ry ! rr(1), rr(2), and their absolute values
    !
    rx  = rr(1)
    ry  = rr(2)
    !
    Wu(-1,-1,-1)= (rx*ry*(-189.0_rk + 154.0_rk*rx + 154.0_rk*ry - 130.0_rk*rx*ry))/3360.0_rk
    Wu(-1,-1, 0)= (rx*ry*(-56.0_rk + 35.0_rk*rx + 35.0_rk*ry - 24.0_rk*rx*ry))/1680.0_rk
    Wu(-1,-1, 1)= (rx*ry*(21.0_rk - 14.0_rk*ry + 2.0_rk*rx*(-7.0_rk + 5.0_rk*ry)))/3360.0_rk
    Wu(-1, 0,-1)= (rx*(245.0_rk - 154.0_rk*ry**2 + rx*(-189.0_rk + 130.0_rk*ry**2)))/1680.0_rk
    Wu(-1, 0, 0)= (rx*(-35.0_rk*(-3.0_rk + ry**2) + 8.0_rk*rx*(-7.0_rk + 3.0_rk*ry**2)))/840.0_rk
    Wu(-1, 0, 1)= (rx*(-35.0_rk + 21.0_rk*rx + 14.0_rk*ry**2 - 10.0_rk*rx*ry**2))/1680.0_rk
    Wu(-1, 1,-1)= (rx*ry*(7.0_rk*(27.0_rk + 22.0_rk*ry) - 2.0_rk*rx*(77.0_rk + 65.0_rk*ry)))/3360.0_rk
    Wu(-1, 1, 0)= (rx*ry*(56.0_rk - 35.0_rk*rx + 35.0_rk*ry - 24.0_rk*rx*ry))/1680.0_rk
    Wu(-1, 1, 1)= (rx*ry*(-7.0_rk*(3.0_rk + 2.0_rk*ry) + 2.0_rk*rx*(7.0_rk + 5.0_rk*ry)))/3360.0_rk
    Wu( 0,-1,-1)= (ry*(245.0_rk - 189.0_rk*ry + 2.0_rk*rx**2*(-77.0_rk + 65.0_rk*ry)))/1680.0_rk
    Wu( 0,-1, 0)= (ry*(105.0_rk - 56.0_rk*ry + rx**2*(-35.0_rk + 24.0_rk*ry)))/840.0_rk
    Wu( 0,-1, 1)= -(ry*(35.0_rk - 21.0_rk*ry + 2.0_rk*rx**2*(-7.0_rk + 5.0_rk*ry)))/1680.0_rk
    Wu( 0, 0,-1)= -0.4166666666666667_rk + (9.0_rk*ry**2)/40.0_rk + rx**2*(0.225_rk -(13.0_rk*ry**2)/84.0_rk)
    Wu( 0, 0, 0)= (-2.0_rk*(-7.0_rk*(-5.0_rk + ry**2) + rx**2*(-7.0_rk + 3.0_rk*ry**2)))/105.0_rk
    Wu( 0, 0, 1)= (70.0_rk - 21.0_rk*ry**2 + rx**2*(-21.0_rk + 10.0_rk*ry**2))/840.0_rk
    Wu( 0, 1,-1)= (ry*(-7.0_rk*(35.0_rk + 27.0_rk*ry) + 2.0_rk*rx**2*(77.0_rk + 65.0_rk*ry)))/1680.0_rk
    Wu( 0, 1, 0)= (ry*(-7.0_rk*(15.0_rk + 8.0_rk*ry) + rx**2*(35.0_rk + 24.0_rk*ry)))/840.0_rk
    Wu( 0, 1, 1)= (ry*(7.0_rk*(5.0_rk + 3*ry) - 2.0_rk*rx**2*(7.0_rk + 5.0_rk*ry)))/1680.0_rk
    Wu( 1,-1,-1)= (rx*ry*(189.0_rk - 154.0_rk*ry - 2.0_rk*rx*(-77.0_rk + 65.0_rk*ry)))/3360.0_rk
    Wu( 1,-1, 0)= (rx*(56.0_rk + rx*(35.0_rk - 24.0_rk*ry) - 35.0_rk*ry)*ry)/1680.0_rk
    Wu( 1,-1, 1)= (rx*ry*(7.0_rk*(-3.0_rk + 2.0_rk*ry) + 2.0_rk*rx*(-7.0_rk + 5.0_rk*ry)))/3360.0_rk
    Wu( 1, 0,-1)= (rx*(-245.0_rk - 189.0_rk*rx + 154.0_rk*ry**2 + 130.0_rk*rx*ry**2))/1680.0_rk
    Wu( 1, 0, 0)= (rx*(35.0_rk*(-3.0_rk + ry**2) + 8.0_rk*rx*(-7.0_rk + 3.0_rk*ry**2)))/840.0_rk
    Wu( 1, 0, 1)= (rx*(35.0_rk - 14.0_rk*ry**2 + rx*(21.0_rk - 10.0_rk*ry**2)))/1680.0_rk
    Wu( 1, 1,-1)= -(rx*ry*(7.0_rk*(27.0_rk + 22.0_rk*ry) + 2.0_rk*rx*(77.0_rk + 65.0_rk*ry)))/3360.0_rk
    Wu( 1, 1, 0)= -(rx*ry*(56.0_rk + 35.0_rk*rx + 35.0_rk*ry + 24.0_rk*rx*ry))/1680.0_rk
    Wu( 1, 1, 1)= (rx*ry*(7.0_rk*(3.0_rk + 2.0_rk*ry) + 2.0_rk*rx*(7.0_rk + 5.0_rk*ry)))/3360.0_rk
  end subroutine GetMatrixWu
  
end module phase_interpolation
