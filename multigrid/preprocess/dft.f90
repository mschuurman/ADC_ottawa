!
!  Exchange-correlation functionals. All functionals assume spin-unpolarized
!  inputs (anything else does not map well on our data structures elsewhere)
!
module dft
  use accuracy
  implicit none
  private
  public :: dft_set_functional, dft_epsilon_xc, dft_v_xc
  !
  real(rk), parameter    :: rho_cut    = 1e-20_rk ! Do not allow densities below rho_cut
  integer(ik), parameter :: fnc_none   = 0        ! Do not use DFT
  integer(ik), parameter :: fnc_slater = 1        ! Use Slater-Dirac exchange
  integer(ik), parameter :: fnc_svwn   = 2        ! Use Slater-Dirac exchange + Vosko-Wilk-Nusair 
                                                  !     correlation.
  integer(ik), save      :: functional = fnc_none ! Functional to use
  !
  contains
  !
  subroutine dft_set_functional(vxc)
    character(len=*), intent(in) :: vxc
    !
    select case (vxc)
      case default
        write (out,"('dft%dft_set_functional: Functional ',a,' is not recognized')") trim(vxc)
        stop 'dft%dft_set_functional - bad functional'
      case ('None','none',' ')
        functional = fnc_none
      case ('Slater','slater')
        functional = fnc_slater
      case ('VWN','SVWN','vwn','svwn')
        functional = fnc_svwn
    end select
  end subroutine dft_set_functional
  !
  function dft_epsilon_xc(r,rho) result(eps)
    real(rk), intent(in)    :: r(*) ! Coordinates of a grid point; not used
    complex(rk), intent(in) :: rho  ! Electron density at a point - use only the real part
    complex(rk)             :: eps  ! Exchange-correlation energy density
    !
    real(rk) :: x, resX(3), resC(3)
    !
    x = max(rho_cut,0.5_rk*real(rho,kind=rk))
    select case (functional)
      case default
        write (out,"()") functional
        stop 'dft%dft_epsilon_xc - bad functional'
      case (fnc_none)
        eps = 0
      case (fnc_slater)
        call dft_slater(x,x,resX)
        eps = resX(1)
      case (fnc_svwn)
        call dft_slater(x,x,resX)
        call dft_vwn(x,x,resC)
        eps = resX(1) + resC(1)
    end select
  end function dft_epsilon_xc
  !
  function dft_v_xc(r,rho) result(v)
    real(rk), intent(in)    :: r(*) ! Coordinates of a grid point; not used
    complex(rk), intent(in) :: rho  ! Electron density at a point - use only the real part
    complex(rk)             :: v    ! Exchange-correlation potential
    !
    real(rk) :: x, resX(3), resC(3)
    !
    x = max(rho_cut,0.5_rk*real(rho,kind=rk))
    select case (functional)
      case default
        write (out,"()") functional
        stop 'dft%dft_v_xc - bad functional'
      case (fnc_none)
        v = 0
      case (fnc_slater)
        call dft_slater(x,x,resX)
        v = resX(2)
      case (fnc_svwn)
        call dft_slater(x,x,resX)
        call dft_vwn(x,x,resC)
        v = resX(2) + resC(2)
    end select
  end function dft_v_xc
  !
  !  Internal routines, one per functional. All internal subroutines 
  !  share the same calling convention: The first two arguments are
  !  sanitized spin-densities; the last argument is a 3-element array,
  !  filled with the total energy density (including the overall density
  !  factor), and the alpha/beta exchange-correlation potentials.
  !
  !  The routines below are computer-generated - don't look for any meaning!
  !
  !  Slater-Dirac exchange
  !
  subroutine dft_slater(x,y,res)
    real(rk), intent(in)  :: x, y   ! alpha-/beta-spin density
    real(rk), intent(out) :: res(3) ! XC energy density and XC potentials
    !
    real(rk), parameter :: CTHIRD=0.3333333333333333333333333333333333333333_rk
    real(rk), parameter :: CN1=-0.9305257363491000250020102180716672510262_rk
    real(rk), parameter :: CN2=-1.240700981798800033336013624095556334702_rk
    !
    real(rk) :: T1, T2, T3, T4
    !
    T1 = x**CTHIRD
    T2 = T1**2
    T3 = y**CTHIRD
    T4 = T3**2
    res(1) = CN1*(T2**2 + T4**2)
    res(2) = CN2*T1
    res(3) = CN2*T3
  end subroutine dft_slater
  !
  !  Vosko-Wilk-Nusair correlation
  !
  subroutine dft_vwn(x,y,res)
    real(rk), intent(in)  :: x, y   ! alpha-/beta-spin density
    real(rk), intent(out) :: res(3) ! XC energy density and XC potentials
    !
    real(rk), parameter :: CTHIRD=0.3333333333333333333333333333333333333333_rk
    real(rk), parameter :: CSIXTH=0.1666666666666666666666666666666666666667_rk
    real(rk), parameter :: CN1=0.0310907_rk
    real(rk), parameter :: CN2=-1.0311676_rk
    real(rk), parameter :: CN3=20.85143832359383809660210955464847854224_rk
    real(rk), parameter :: CN4=4.732516058487829936045617011065299207383_rk
    real(rk), parameter :: CN5=0.0623352_rk
    real(rk), parameter :: CN6=0.1332870645322399251727912116148442659818_rk
    real(rk), parameter :: CN7=-0.3333333333333333333333333333333333333333_rk
    real(rk), parameter :: CN8=1.24742_rk
    real(rk), parameter :: CN9=6.15199_rk
    real(rk), parameter :: CN10=3.72744_rk
    real(rk), parameter :: CN11=1.575246635799486501056023072482332358073_rk
    real(rk), parameter :: CN12=-0.03205972250132_rk
    real(rk), parameter :: CN13=0.00193804500264_rk
    real(rk), parameter :: CN14=-0.01036356666666666666666666666666666666667_rk
    real(rk), parameter :: CN15=0.038783160994_rk
    real(rk), parameter :: CN16=-2.240915430552281694725956725367747673801e-7_rk
    real(rk), parameter :: CN17=-11.85741659049244755676157582286297427344_rk
    real(rk), parameter :: CN18=-1.490274674851275566921548161864495100751_rk
    real(rk), parameter :: CN19=-2.685909448833100455289563871445913715253_rk
    real(rk), parameter :: CN20=-0.5524416463222842320849333333333333333333_rk
    real(rk), parameter :: CN21=2.481401963597600066672027248191112669403_rk
    real(rk), parameter :: CN22=143.7994006102733840268840677628192772197_rk
    real(rk), parameter :: CN23=1144.144335362144148070427053976150871037_rk
    real(rk), parameter :: CN24=510.9285034195520059062446733810175407441_rk
    real(rk), parameter :: CN25=23.817288064230212039746773169472570201_rk
    real(rk), parameter :: CN26=162.1872789377_rk
    real(rk), parameter :: CN27=-0.3039633_rk
    real(rk), parameter :: CN28=-1.000414034_rk
    real(rk), parameter :: CN29=20.96314936600717944270379535712058098851_rk
    real(rk), parameter :: CN30=1.436054487335498306012468619939054142922_rk
    real(rk), parameter :: CN31=0.000828068_rk
    real(rk), parameter :: CN32=0.006041466640028676509260903994552057108474_rk
    real(rk), parameter :: CN33=0.317708_rk
    real(rk), parameter :: CN34=7.12311_rk
    real(rk), parameter :: CN35=1.13107_rk
    real(rk), parameter :: CN36=1.25992104989487316476721060727822835057_rk
    real(rk), parameter :: CN37=-1.333333333333333333333333333333333333333_rk
    real(rk), parameter :: CN38=1.67989473319316421968961414303763780076_rk
    real(rk), parameter :: CN39=1.3171743_rk
    real(rk), parameter :: CN40=-0.6079266_rk
    real(rk), parameter :: CN41=8.127348033960001478852419210706210737855e-8_rk
    real(rk), parameter :: CN42=-363.7880102214421262725848788029313016209_rk
    real(rk), parameter :: CN43=-2.1959954443046128310989266602001058474_rk
    real(rk), parameter :: CN44=-25.0152897236158611454578158774250633911_rk
    real(rk), parameter :: CN45=-17.34650933101205134588266666666666666667_rk
    real(rk), parameter :: CN46=6.587986332913838493296779980600317542199_rk
    real(rk), parameter :: CN47=1091.364030404065324423644050018394001981_rk
    real(rk), parameter :: CN48=150.0608327611678743565263090033866936343_rk
    real(rk), parameter :: CN49=7.141848156598477340130994026746566501567_rk
    real(rk), parameter :: CN50=109.196349464504_rk
    real(rk), parameter :: CN51=3.847322101863072639518916246550536610962_rk
    real(rk), parameter :: CN52=-0.01779322315535_rk
    real(rk), parameter :: CN53=29.10902830723860547813884392324287957047_rk
    real(rk), parameter :: CN54=8.964208955655528791088574264714028885721_rk
    real(rk), parameter :: CN55=0.03205972250132_rk
    real(rk), parameter :: CN56=-0.00193804500264_rk
    real(rk), parameter :: CN57=0.0044957463107_rk
    real(rk), parameter :: CN58=0.4126337966562962057644993691638825151846_rk
    real(rk), parameter :: CN59=0.005181783333333333333333333333333333333333_rk
    real(rk), parameter :: CN60=0.052491361531_rk
    real(rk), parameter :: CN61=4.73093_rk
    real(rk), parameter :: CN62=7.06042_rk
    real(rk), parameter :: CN63=-0.038783160994_rk
    real(rk), parameter :: CN64=-16.6717291080733147712486370683856586475_rk
    real(rk), parameter :: CN65=-15.38928840745229055807566498620214644385_rk
    real(rk), parameter :: CN66=11.74327464028887592699252528258708984955_rk
    real(rk), parameter :: CN67=51.7407899137_rk
    real(rk), parameter :: CN68=22.24380570466282096357193284283141805517_rk
    real(rk), parameter :: CN69=72.2312292413_rk
    real(rk), parameter :: CN70=5.13617152307287475613758698209863668722e-7_rk
    real(rk), parameter :: CN71=11.99334575086991277683014547842639079457_rk
    real(rk), parameter :: CN72=155.698622793677988267704867374003108_rk
    real(rk), parameter :: CN73=0.3240652369125459347303265919135767361909_rk
    real(rk), parameter :: CN74=4862.006894165809045584186219603674714029_rk
    real(rk), parameter :: CN75=25729.75744455533360441321195042812453183_rk
    real(rk), parameter :: CN76=1292.956254647208626639271378066978551807_rk
    real(rk), parameter :: CN77=14217.20776918912179046700939803932558453_rk
    real(rk), parameter :: CN78=34484.58633208389866489456559342658119745_rk
    real(rk), parameter :: CN79=1119.528112838756315483215793386861015614_rk
    real(rk), parameter :: CN80=12864.87872227766680220660597521406226591_rk
    real(rk), parameter :: CN81=646.4781273236043133196356890334892759033_rk
    real(rk), parameter :: CN82=7108.603884594560895233504699019662792265_rk
    real(rk), parameter :: CN83=17242.29316604194933244728279671329059872_rk
    real(rk), parameter :: CN84=15.38928840745229055807566498620214644385_rk
    !
    real(rk) :: T1, T2, T3, T4, T5, T6, T7, T8, T9, T10
    real(rk) :: T11, T12, T13, T14, T15, T16, T17, T18, T19, T20
    real(rk) :: T21, T22, T23, T24, T25, T26, T27, T28, T29, T30
    real(rk) :: T31, T32, T33, T34, T35, T36, T37, T38, T39, T40
    real(rk) :: T41, T42, T43, T44, T45, T46, T47, T48, T49, T50
    real(rk) :: T51, T52, T53, T54, T55, T56, T57, T58, T59, T60
    real(rk) :: T61, T62, T63, T64, T65
    !
    T1 = x + y
    T8 = T1**CTHIRD
    T14 = 1/T8
    T7 = T1**CSIXTH
    T2 = 1/T7
    T3 = LOG(CN3 + T14 + CN4*T2)
    T4 = LOG(CN6 + T2)
    T5 = LOG(T1)
    T15 = CN11*T2
    T6 = ATAN(CN9/(CN10 + T15))
    res(1) = CN1*T1*(CN2*T3 + CN5*T4 + CN7*T5 + CN8*T6)
    T10 = T8**2
    T13 = T10*T8
    T25 = T10**2
    T29 = 1/(T13**3*T25)
    T23 = T29*y
    T27 = x**2
    T40 = y**2
    T21 = T27 + T40
    T22 = CN7*T5 + CN33*ATAN(CN34/(CN35+T15)) + CN31*LOG(CN32+T2) + CN28*LOG(CN29+T14+CN30*T2)
    T16 = x**CTHIRD
    T17 = T16**2
    T33 = T17**2
    T18 = y**CTHIRD
    T19 = T18**2
    T34 = T19**2
    T35 = -(T8*x) - T8*y
    T28 = T22*(CN36*T33 + CN36*T34 + T35)
    T24 = T21*T28
    res(1) = res(1) + CN27*T23*T24*x
    T51 = LOG(CN53 + T14 + CN54*T2)
    T52 = LOG(CN58 + T2)
    T53 = ATAN(CN61/(CN62 + T15))
    T48 = CN55*T3 + CN56*T4 + CN59*T5 + CN52*T51 + CN57*T52 + CN60*T53 + CN63*T6
    T45 = (CN36*T33 + CN36*T34 + T35)*T48
    T43 = -x + y
    T44 = T43**2
    T46 = T44**2
    res(1) = res(1) + CN51*T29*T45*T46
    T11 = SQRT(T1)
    T9 = T7**2
    T30 = T9**2
    T32 = T30*T7
    T12 = CN12*T3 + CN13*T4 + CN14*T5 + CN15*T6 + (T7*(CN16+CN18*T10+CN17*T11+CN20*T7+CN19*T8)) / &
                                                  (CN21+CN23*T10+CN24*T11+CN22*T32+CN25*T7+CN26*T8)
    res(2) = T12
    T26 = T25**2
    T47 = T26**(-2)
    T39 = CN39*T24*T47*x*y
    T31 = T30**2
    T42 = (CN27*T21*(CN36*T33+CN36*T34+T35)*(CN41+CN43*T10+CN42*T11+CN45*T7+CN44*T8)*x*y) / &
          (T30*T31**3*T7*(CN21+CN47*T10+CN48*T11+CN46*T32+CN49*T7+CN50*T8)*T9)
    T20 = T1**2
    T50 = 1/(T1*T20**2)
    T36 = T50*x*y
    T37 = CN37*y
    T38 = CN37*x
    T49 = T10*T16
    T41 = T28*T29
    res(2) = res(2) + CN27*T23*T24 + T39 + T42 + CN27*T21*T22*T36*(T37+T38+CN38*T49) + CN40*T27*T41*y
    T59 = CN64*T45*T46*T47
    T58 = T30*T32
    T57 = CN80 + CN81*T11 + CN83*T7 + CN82*T8
    T54 = CN72*T11 + CN73*T7 + CN71*T8
    T55 = CN74*T32
    T56 = CN79*T10
    T60 = (CN51*(CN36*T33 + CN36*T34 + T35)*T46* (T27*T57 + y*(CN70 + T54 + T55 + T56 + T57*y) + &
                   x*(CN70 + T54 + T55 + T56 + (CN75 + CN76*T11 + CN78*T7 + CN77*T8)*y))) / &
          (T32**2*T58**3*(1 + CN58*T7)*(1 + CN6*T7)*(1 + CN4*T7 + CN3*T8)* (1 + CN54*T7 + CN53*T8)* &
                   (CN21 + CN66*T7 + CN67*T8)* (CN21 + CN68*T7 + CN69*T8))
    T61 = CN37*y
    T63 = CN37*x
    T64 = T46*T50
    T65 = T29*T43*T44
    res(2) = res(2) + T59 + T60 + CN51*T48*(CN38*T49 + T61 + T63)*T64 + CN65*(CN36*T33 + CN36*T34 + T35)* &
                         (CN55*T3 + CN56*T4 + CN59*T5 + CN52*T51 + CN57*T52 + CN60*T53 + CN63*T6)* T65
    res(3) = T12
    T62 = T10*T18
    res(3) = res(3) + T39 + T42 + CN27*T21*T22*T36*(T37 + T38 + CN38*T62) + CN27*T24*T29*x + CN40*T40*T41*x
    res(3) = res(3) + T59 + T60 + CN51*T48*(T61 + CN38*T62 + T63)*T64 + CN84*T45*T65
  end subroutine dft_vwn
end module dft
