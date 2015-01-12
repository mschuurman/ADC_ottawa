!
!  Evaluation of Coulomb attraction integrals over cuboidal basis functions.
!  The expression used in this module trades accuracy for simplicity; a more
!  numerically stable, but unwieldy expression can be found in: 
!  Michael E. Mura and Nicholas !  C. Handy, Theor Chim Acta 90: 145-165 (1995). 
!  The implementation used presently should still give 9-10 significant digits
!  in double precision as long as the volume elements are not too far from cubes; 
!  this is more than sufficient for our purposes. The accuracy _will_ suffer
!  greatly in single precision - so beware.
!
  module cubeint
    use accuracy
    implicit none
    private
    public cube_potential
    !
    !  Volume elements where the lower boundary is smaller than (zero_fraction_tol)
    !  times the upper boundary are treated as starting on-axis.
    !
    real(rk), parameter :: zero_fraction_tol = 1e-6_rk
    !
    contains
    !
    !  Out primitive integral: 1/r integrated over a cube, extending from
    !  (0,0,0) to (x,y,z), where x, y, and z are assumed to be strictly positive.
    !
    function corner_integral(x,y,z) result(v)
      real(rk), intent(in) :: x, y, z 
      real(rk)             :: v
      !
      real(rk) :: x2, y2, z2, r
      !
      if (x<=0 .or. y<=0 .or. z<=0) then
        write (out,"('cubeint%corner_integral called for volume not in Ist octant: ',3(g14.7,1x))") x, y, z
        stop 'cubeint%corner_integral - bad volume'
      endif
      !
      x2 = x**2
      y2 = y**2
      z2 = z**2
      r   = sqrt(x2+y2+z2)
      !
      v = y*z*log(r+x) + x*z*log(r+y) + x*y*log(r+z) &
          -0.5_rk * ( z2*atan2(x*y,r*z) + y2*atan2(x*z,r*y) + x2*atan2(y*z,r*x) &
                    + x*y*log(x2+y2) + x*z*log(x2+z2) + y*z*log(y2+z2) )
      ! write (out,"('Corner integral for: ',3(g15.8,1x),' is ',g25.15)") x, y, z, v
    end function corner_integral
    !
    !  Integral for a volume entirely contained in a single octant
    !
    function octant_integral(vol) result(v)
      real(rk), intent(in) :: vol(2,3) ! Corners of the volume element
      real(rk)             :: v        ! The integral
      !
      real(rk) :: vr(2,3) ! Volume "flipped" into the Ist octant
      logical  :: zr(3)   ! True if the lower boundary of the volume is zero
      !
      if (any(vol(1,:)*vol(2,:)<0._rk)) then
        write (out,"('cubeint%octant_integral called for volume spanning several octants:'/3(2x,g14.7,' : ',g14.7,5x))") vol
        stop 'cubeint%octant_integral - bad octant'
      end if
      !
      !  Flip the cube into the first octant
      !
      vr(1,:) = minval(abs(vol),dim=1)
      vr(2,:) = maxval(abs(vol),dim=1)
      ! write (out,"('Original octant: ',3(g15.8,1x,g15.8,3x))") vol
      ! write (out,"(' Reduced octant: ',3(g15.8,1x,g15.8,3x))") vr
      !
      !  Check for lower boundary being on-axis
      !
      zr = vr(1,:) < zero_fraction_tol * vr(2,:)
      !
      v = corner_integral(vr(2,1),vr(2,2),vr(2,3))
      if (.not.    zr(1))        v = v - corner_integral(vr(1,1),vr(2,2),vr(2,3))
      if (.not.    zr(2))        v = v - corner_integral(vr(2,1),vr(1,2),vr(2,3))
      if (.not.    zr(3))        v = v - corner_integral(vr(2,1),vr(2,2),vr(1,3))
      if (.not.any(zr((/1,2/)))) v = v + corner_integral(vr(1,1),vr(1,2),vr(2,3))
      if (.not.any(zr((/1,3/)))) v = v + corner_integral(vr(1,1),vr(2,2),vr(1,3))
      if (.not.any(zr((/2,3/)))) v = v + corner_integral(vr(2,1),vr(1,2),vr(1,3))
      if (.not.any(zr))          v = v - corner_integral(vr(1,1),vr(1,2),vr(1,3))
    end function octant_integral
    !
    !  Volume integral for a volume element which could span multiple octants along Z
    !
    function volume_integral3(vol) result(v)
      real(rk), intent(in) :: vol(2,3) ! Corners of the volume element
      real(rk)             :: v        ! The integral
      !
      real(rk) :: tv(2,3) ! Part of the volume element falling into a given octant
      !
      !  See whether we need to split the volume element along the Z axis
      ! 
      if (vol(1,3)*vol(2,3)<0._rk) then
        tv = vol
        tv(1,3) = 0.0_rk 
        v  = octant_integral(tv)
        tv = vol
        tv(2,3) = 0.0_rk
        v  = v + octant_integral(tv)
      else
        v = octant_integral(vol)
      end if
    end function volume_integral3
    !
    !  Volume integral for a volume element which could span multiple octants along Y/Z
    !
    function volume_integral2(vol) result(v)
      real(rk), intent(in) :: vol(2,3) ! Corners of the volume element
      real(rk)             :: v        ! The integral
      !
      real(rk) :: tv(2,3) ! Part of the volume element falling into a given octant
      !
      !  See whether we need to split the volume element along the Y axis
      ! 
      if (vol(1,2)*vol(2,2)<0._rk) then
        tv = vol
        tv(1,2) = 0.0_rk 
        v  = volume_integral3(tv)
        tv = vol
        tv(2,2) = 0.0_rk
        v  = v + volume_integral3(tv)
      else
        v = volume_integral3(vol)
      end if
    end function volume_integral2
    !
    !  Volume integral for a volume element which could span multiple octants
    !
    function volume_integral(vol) result(v)
      real(rk), intent(in) :: vol(2,3) ! Corners of the volume element
      real(rk)             :: v        ! The integral
      !
      real(rk) :: tv(2,3) ! Part of the volume element falling into a given octant
      !
      !  See whether we need to split the volume element along the X axis
      ! 
      if (vol(1,1)*vol(2,1)<0._rk) then
        tv = vol
        tv(1,1) = 0.0_rk 
        v  = volume_integral2(tv)
        tv = vol
        tv(2,1) = 0.0_rk
        v  = v + volume_integral2(tv)
      else
        v = volume_integral2(vol)
      end if
    end function volume_integral
    !
    !  Asymptotic form of the integral, including the 2nd-order correction
    !
    function asymptotic_potential(r,dr) result(v)
      real(rk), intent(in) :: r (:) ! Coordinates of the centre of the cuboid
      real(rk), intent(in) :: dr(:) ! Extent of the cuboid
      real(rk)             :: v
      !
      real(rk)    :: r12   ! Distance to the centre
      real(rk)    :: dc(3) ! Direction cosines
      !
      r12 = sqrt(sum(r**2))
      dc  = r/r12
      v   = 1._rk/r12 + sum(dr**2 * (3*dc**2-1)) / (24._rk*r12**3)
    end function asymptotic_potential
    !
    !  Finally, the externally-visible interface.
    !  Returns 1/r averaged over cuboid volume element with extent dr,
    !  centered around point r. Switch between the "exact" expression and
    !  the second-order approximant, based on the expected accuracy of each expression
    !
    function cube_potential(r,dr) result(v)
      real(rk), intent(in) :: r (:) ! Coordinates of the centre of the cuboid
      real(rk), intent(in) :: dr(:) ! Extent of the cuboid
      real(rk)             :: v
      !
      real(rk) :: vol(2,3) ! Volume element
      real(rk) :: vbox     ! Volume of the volume element
      real(rk) :: rc       ! Distance to the centre of the volume
      real(rk) :: vint     ! Characteristic magnitude of the elementary volume integral
      real(rk) :: eps_v    ! Expected absolute error in the volume integral expression
      real(rk) :: eps_a    ! Expected absolute error in the asymptotic expression
      logical  :: use_integral
      !
      !  Decide on which expression to use. One special case is where 
      !  the volume element encloses the charge; then we must use the
      !  exact expression. Otherwise, go with the result expected to
      !  yield more significant digits.
      !
      vbox = product(dr)
      if (all(abs(r)-0.5_rk*abs(dr)<=0._rk)) then
        use_integral = .true.
        ! write (out,"('Overlapping origin test is true')")
      else
        rc    = sqrt(sum(r**2))
        vint  = (0.75_rk * sqrtpi)** (2._rk/3._rk) * maxval(abs(r)+abs(dr)) ** 2
        ! write (out,"('Expected primitive integral magnitude = ',g14.7)") vint
        eps_v = spacing(vint*rc/vbox)
        eps_a = maxval(dr)**4/(20._rk*rc**5)
        ! write (out,"('Integral/asymptotic error estimates are: ',g14.7,' / ',g14.7)") eps_v, eps_a
        use_integral = eps_v < eps_a
      end if
      !
      if (use_integral) then
        vol(1,:) = r - 0.5_rk * dr
        vol(2,:) = r + 0.5_rk * dr
        v = volume_integral(vol)/vbox
        ! write (out,"('Volume integral result is: ',g25.15)") v
      else
        v = asymptotic_potential(r,dr)
        ! write (out,"('     Asymptotic result is: ',g25.15)") v
      end if
      !
    end function cube_potential
  end module cubeint
