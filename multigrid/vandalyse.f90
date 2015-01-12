module vandalyse
  use accuracy
  use multigrid
  use timer
  implicit none
  private
  public AnalyzeSlice1D, AnalyzePolar
  !
  !  divideThreshold controls calculation of the angular distribution 
  !  function. For probabilities above divideThreshold, Cartesian boxes 
  !  falling into several angular bins will be subdivided; for smalled
  !  probabilities, they will be simply assigned to whatever boxes their
  !  corners fall in.
  !
  real(rk), parameter :: divideThreshold = 1.0e-6_rk
  !
  integer(ik)         :: angNphi        ! Number of angular bins
  integer(ik)         :: angNtheta      ! along phi and theta
  real(rk), pointer   :: angData(:,:)

  contains

  subroutine AnalyzeSlice1D(typ,src,grd)
 !
 ! AnalyzeSlice1D prints out three reduced 1D distributions
 ! Picking the direction to make statistics about we integrate over 
 ! the remaining two dimensions: 
 !
 !          / 
 ! rho(x) = | |psi(x,y,z)|^2 dy dz 
 !          /
 !
 ! AnalyzeSlice works only at the specified grid level - it does not
 ! traverse the entire multigrid.
 !
    character(len=*), intent(in)  :: typ            ! Field type, either 'COORDINATE', or 'MOMENTUM'
    integer(ik), intent(in)       :: src            ! Field to return
    integer(ik), intent(in)       :: grd            ! Grid level


    integer(ik)                   :: npts(3)        ! Number of points at this grid level
    real(rk)                      :: step(3)        ! Grid spacing
    real(rk), pointer             :: coord(:,:,:,:) ! Coordinates and integration weights
    complex(rk), pointer          :: field(:,:,:)   ! Field values
    real(rk), allocatable         :: v(:)           ! Field slice value 
    integer(ik)                   :: alpha          ! Surviving coordinate
    integer(ik)                   :: ih, iv         ! Two dimensions, which are integrated out
    real(rk)                      :: wgt            ! Partial integration weight
    integer(ik)                   :: mn(3),mx(3)    ! 2D slice of the grid
    integer(ik)                   :: i, alloc

    !
    call TimerStart('AnalyzeSlice1D')
    write(out,"('   Axis     Coord           Field ')")
    !
    ! Take the Field values from the grid 
    !
    call FieldCheckOut(typ,src,grd,npts,step,coord,field)
    allocate (v(maxval(npts)),stat=alloc)
    if (alloc/=0) then
      write (out,"(' AnalyzeSlice1D: error ',i8,' allocating ',i8,' words')") &
             alloc, maxval(npts)
      stop 'AnalyzeSlice1D - no memory'
    end if
    !
    ! first, we choose the surving coordinate. All remaining coordinates will be
    ! integrated out.
    !
    do alpha=1,3
      ih = 1+mod(alpha+0,3)
      iv = 1+mod(alpha+1,3)
      !
      ! Partial integration weight
      !
      wgt = step(ih) * step(iv)
      !
      ! second,  we pick up a point 
      !
      !$omp parallel do private(i,mn,mx)
      do i=1,npts(alpha) 
        mn(ih   ) = 1 ; mx(ih   ) = npts(ih)
        mn(iv   ) = 1 ; mx(iv   ) = npts(iv)
        mn(alpha) = i ; mx(alpha) = i 
        !
        ! third,  we calculate 2D intergral 
        !
        v(i) = wgt*sum(real( field(mn(1):mx(1),mn(2):mx(2),mn(3):mx(3))*  & 
                       conjg(field(mn(1):mx(1),mn(2):mx(2),mn(3):mx(3))),kind=rk))
      end do
      !$omp end parallel do
      do i=1,npts(alpha)
        mn(ih   ) = 1 ; mx(ih   ) = npts(ih)
        mn(iv   ) = 1 ; mx(iv   ) = npts(iv)
        mn(alpha) = i ; mx(alpha) = i 
        write(out,"('#  ',i2,1x,f12.6,f20.9)") alpha,coord(alpha,mn(1),mn(2),mn(3)),v(i)
      end do
    end do
    !
    deallocate (v)
    call TimerStop('AnalyzeSlice1D')
    !
  end subroutine AnalyzeSlice1D
  !
  ! Transform the field to the spherical coordinates (R,theta,phi) 
  ! phi = (0,2pi) ; theta = (0,pi)
  ! and integrates over the Radius R to get the angle distribution of the Field
  ! 
  ! Results is 
  !                  / 
  ! rho(phi,theta) = | |psi(R,phi,theta)|^2  R^2dR 
  !                  /
  !
  subroutine AnalyzePolar(typ,src,grd,Nphi,Ntheta)
    character(len=*), intent(in) :: typ                ! Field type, either 'COORDINATE', or 'MOMENTUM'
    integer(ik), intent(in)      :: src                ! Field to return
    integer(ik), intent(in)      :: grd                ! Grid level
    integer(ik), intent(in)      :: Nphi,Ntheta        ! Sizes of the spherical grid 

    integer(ik)            :: npts(3)            ! Number of points at this grid level
    real(rk)               :: step(3)            ! Grid spacing
    real(rk), pointer      :: coord(:,:,:,:)     ! Coordinates and integration weights
    complex(rk), pointer   :: field(:,:,:)       ! Field values
    integer(ik)            :: ix,iy,iz
    real(rk)               :: x,y,z
    real(rk), target       :: rho(Nphi,Ntheta)   ! spherical angle distribution of the field 
    integer(ik)            :: itheta,iphi
    real(rk)               :: density,norm
    !
    call TimerStart('AnalyzePolar')
    !
    !  Set up module variables, which will be accessed from the recursive
    !  subdivision routines
    !
    angNphi   = Nphi
    angNtheta = Ntheta
    angData   => rho
    rho       = 0
    !
    !  Take the Field values from the grid 
    !
    call FieldCheckOut(typ,src,grd,npts,step,coord,field)
    !
    write(out,"('          phi        theta          Field')") 
    !
    !$omp parallel do private(density,ix,iy,iz,x,y,z)
    do iz=1,npts(3)
      do iy=1,npts(2)
        do ix=1,npts(1)
          !
          !  Number of particles inside the brick
          !
          density = product(step)*real(field(ix,iy,iz)*conjg(field(ix,iy,iz)),kind=rk)
          !
          x = coord(1,ix,iy,iz)
          y = coord(2,ix,iy,iz)
          z = coord(3,ix,iy,iz)
          !
          !  Assign the density to one, or more sperical patches. If necessary,
          !  DividingBricks will recursively subdivide the brick.
          !
          call DividingBricks(x,y,z,density,step)
        end do
      end do
    end do
    !$omp end parallel do
    !
    ! Print out of the spherical angle distribution 
    !
    norm = 0.0_rk
    do iphi=1,Nphi
      do itheta=1,Ntheta
        write(out,"('! ',2f12.4,2x,f16.8)") ((iphi-0.5_rk)*2.0_rk*pi)/Nphi, &
             ((itheta-0.5_rk)*pi)/Ntheta, rho(iphi,itheta)
        norm = norm + rho(iphi,itheta)
      end do
    end do

    write(out,"('! norm ',f12.4)") norm
    call TimerStop('AnalyzePolar')

  end subroutine AnalyzePolar

  !
  ! DividingBricks does the following:
  ! a) If all corners of a brick fall into the same spherical patch,
  !    the total density is assigned to that patch
  ! b) If the total density is smaller than the threshold, one-eighth
  !    of it is assigned to each of the patch, containing the corners
  ! c) Otherwise, the brick is subdivided by a factor of 2 in each direction
  !    (giving the total of 8 new smaller bricks), and the process is
  !    repeated for each small brick
  !
  recursive subroutine DividingBricks(x,y,z,density,h)
    real(rk), intent(in)         :: x,y,z      ! Centre of the brick
    real(rk), intent(in)         :: density    ! Total probability of finding an electron inside
                                               ! the (sibdivided) brick
    real(rk)                     :: h(3)       ! Size of the current brick
    !
    integer(ik)                  :: ix, iy, iz, ic
    real(rk)                     :: x0, y0, z0
    integer(ik)                  :: bins(2,8)
    integer(ik)                  :: itheta, iphi
    real(rk)                     :: sdd        ! Subdivided density - density of 1/8th of the brick
    real(rk)                     :: sdh(3)     ! Subdivided size - each dimension halved
    !
    !  Find the bins, corresponding to the corners of the brick
    !
    call CornerBins(x,y,z,h,bins)
    !
    !  If all corners are inside the same bin, we are done with subdivision
    !
    if ( all(bins(:,2:8)==spread(bins(:,1),2,7)) ) then
      itheta = bins(1,1)
      iphi   = bins(2,1)
      !$omp atomic
      angData(iphi,itheta) = angData(iphi,itheta) + density
      return
    end if
    !
    ! If the density is smaller than threshold, assign the corners
    ! anyway, and be done with the subdivision
    !
    sdd = density / 8.0_rk
    if ( density < divideThreshold ) then
      do ic=1,8
        itheta = bins(1,ic)
        iphi   = bins(2,ic)
        !$omp atomic
        angData(iphi,itheta) = angData(iphi,itheta) + sdd
      end do
      return
    end if
    !
    ! This brick must be subdivided further, twice along each direction
    !
    sdh = h/2
    do iz = -1,1,2
      do iy = -1,1,2
        do ix = -1,1,2
          x0 = x+ix*0.25_rk*h(1); y0 = y+iy*0.25_rk*h(2); z0 = z+iz*0.25_rk*h(3)
          call DividingBricks(x0,y0,z0,sdd,sdh)
        end do
      end do
    end do
  end subroutine DividingBricks
  !
  !  CornerBins calculates the angular bins, corresponding to the eight corners
  !  of a (subdivided) brick.
  !
  subroutine CornerBins(x,y,z,h,bins)
    real(rk), intent(in)     :: x,y,z      ! Coordinates of the centre of the box
    real(rk), intent(in)     :: h(3)       ! Size of the box
    integer(ik), intent(out) :: bins(2,8)  ! Indices of the corner bins. The first
                                           ! index is itheta,iphi; the second index
                                           ! corresponds to the eight corners of the
                                           ! brick.
    integer(ik)           :: ix,iy,iz,ic
    real(rk)              :: x0,y0,z0
    real(rk)              :: phi,theta
    integer(ik)           :: iphi, itheta
    !
    ! looping over the 8 corners of the block to check if the it 
    ! within the spherical cone 
    !
    ic = 0
    do iz =-1,1,2
      do iy = -1,1,2
        do ix = -1,1,2
          !
          ! coordinates of the corners 
          !
          x0    = x+ix*0.5_rk*h(1) 
          y0    = y+iy*0.5_rk*h(2)
          z0    = z+iz*0.5_rk*h(3)
          phi   = modulo( atan2(y0,x0), 2.0_rk*pi )
          theta = modulo( atan2(sqrt(x0**2+y0**2),z0), pi )
          !
          iphi   = 1 + modulo(int((angNphi * phi)/(2.0_rk*pi),kind=ik),angNphi)
          itheta = 1 + int((angNtheta * theta)/pi,kind=ik)
          itheta = min(max(1,itheta),angNtheta)
          !
          ic = ic + 1
          bins(1,ic) = itheta
          bins(2,ic) = iphi
        end do 
      end do 
    end do 
  end subroutine CornerBins

end module vandalyse
