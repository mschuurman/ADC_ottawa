!
!  Equation of state and Helmholtz free energy for normal hydrogen.
!  Experimental reference data is from:
!
!    R.L. Mills, D.H. Liebenberg, J.C. Bronson, and L. C. Schmidt
!    "Equation of state for fluid n-H2 from P-V-T and sound velocity
!     measurements to 20 kbar"
!    J. Chem. Phys. 66, 3076-3084 (1977)
!
!  Pressure is calculated from the real solution of the analytical
!  inversion of eq. 11 of the paper. The free energy is calculated
!  from the heat capacity at constant pressure, calculated according
!  to eqs. 4 and 5. The free energy is then calculated from 
!  thermodynamic relations:
!
!    C_p                              = T \frac{\partial S}{\partial T}_P
!    -\frac{\partial S}{\partial P}_T = \frac{\partial V}{\partial T}_P
!    -S                               = \frac{\partial G}{\partial T}_P
!     V                               = \frac{\oartial G}{\partial P}_T
!     G                               = F - PV
!
!  The referemnce values for the enthropy and free energy at 100MPa and
!  300K are taken from:
!
!    R.D. McCarty, J. Hord, H.M. Roder,
!    "Selected Properties of Hydrogen (Engineering Design Data)"
!    U.S. NBS Monograph 168 (1981).
!
!  The nominal range of applicability is:
!
!    75     < T < 307
!    2 kbar < P < 20 kbar
!
!  Comparison to the NBS tabulated data suggests that the equation
!  of state and the free energy are within 10% of the experimental
!  data for pressures down to 1 bar, and behave sensiblly for
!  pressures between 10^-3 and 10^12 Pa.
!
module h2_thermo
  use accuracy
  implicit none
  private
  public h2_tp_f, h2_tv_p, h2_tp_v
  !
  contains
  !
  !  Equation of state: V(T,P) and it's derivatives
  !
  subroutine h2_tp_v(t,p,v,dvdp)
    real(rk), intent(in)            :: t    ! Temperature, Kelvin
    real(rk), intent(in)            :: p    ! Pressure, Pascal
    real(rk), intent(out)           :: v    ! Molar volume, m^3/mole
    real(rk), intent(out), optional :: dvdp ! d V/d P
    !
    real(ark), parameter :: CTHIRD=0.3333333333333333333333333333333333333333_ark
    real(ark), parameter :: CN1=-3181.098375_ark
    real(ark), parameter :: CN2=0.01711699663629408294496037107364086908782_ark
    real(ark), parameter :: CN3=-364.1736837540716913849743736954763694908_ark
    real(ark), parameter :: CN4=146.1075109460260656440459557793639796158_ark
    real(ark), parameter :: CN5=11.14676325_ark
    real(ark), parameter :: CN6=-0.004649605426126382183357521602562151018429_ark
    real(ark), parameter :: CN7=1.538599629555544229852181954851752921084e-6_ark
    real(ark), parameter :: CN8=7015.033725_ark
    real(ark), parameter :: CN9=-0.01047970823039695747139569074420887613643_ark
    real(ark), parameter :: CN10=3181.098375_ark
    real(ark), parameter :: CN11=-0.005705665545431360981653457024546956362607_ark
    real(ark), parameter :: CN12=-673.2346989378807615490114710376187195639_ark
    real(ark), parameter :: CN13=237.1023533217895100671546352053939398138_ark
    real(ark), parameter :: CN14=-11.14676325_ark
    real(ark), parameter :: CN15=-5.128665431851814099507273182839176403615e-7_ark
    real(ark), parameter :: CN16=0.003099736950750921455571681068374767345619_ark
    real(ark), parameter :: CN17=-7015.033725_ark
    real(ark), parameter :: CN18=0.003493236076798985823798563581402958712143_ark
    !
    real(ark) :: t1, t2, t3, t4, t5
    !
    T3 = 1/t
    T1 = p**CTHIRD
    T2 = T1**2
    T4 = t**2
    T5 = SQRT(t)
    v = (T3*(CN1 + CN2*t*(CN3 + T1)*(CN4 + T1) + (CN5 + CN6*T1 + CN7*T2)*T4 + (CN8 + CN9*T2)*T5))/p
    if (.not.present(dvdp)) return
    dvdp = (T3*(CN10 + CN11*t*(CN12 + T1)*(CN13 + T1) + (CN14 + CN16*T1 + CN15*T2)*T4 + (CN17 + CN18*T2)*T5))/p**2
  end subroutine h2_tp_v
  !
  !  Equation of state: P(T,V) and its derivatives
  !
  subroutine h2_tv_p(t,v,p,dpdv)
    real(rk), intent(in)            :: t    ! Temperature, Kelvin
    real(rk), intent(in)            :: v    ! Molar volume, m^3/mole
    real(rk), intent(out)           :: p    ! Pressure, Pascal
    real(rk), intent(out), optional :: dpdv ! d P/d T
    !
    real(ark), parameter :: CTHIRD=0.3333333333333333333333333333333333333333_ark
    real(ark), parameter :: CN1=-3.1395e6_ark
    real(ark), parameter :: CN2=-898860._ark
    real(ark), parameter :: CN3=11001._ark
    real(ark), parameter :: CN4=6.9233e6_ark
    real(ark), parameter :: CN5=261.1071719661484455550136753470364767483_ark
    real(ark), parameter :: CN6=0.000105406271063469140625_ark
    real(ark), parameter :: CN7=-0.0027_ark
    real(ark), parameter :: CN8=3.1395e6_ark
    real(ark), parameter :: CN9=-11001._ark
    real(ark), parameter :: CN10=898860._ark
    real(ark), parameter :: CN11=-6.9233e6_ark
    real(ark), parameter :: CN12=9.e-18_ark
    real(ark), parameter :: CN13=-2.2479e8_ark
    real(ark), parameter :: CN14=33003._ark
    real(ark), parameter :: CN15=3.6716e8_ark
    real(ark), parameter :: CN16=1.7174e7_ark
    real(ark), parameter :: CN17=21393._ark
    real(ark), parameter :: CN18=2._ark
    real(ark), parameter :: CN19=-17.174_ark
    real(ark), parameter :: CN20=-0.021393_ark
    real(ark), parameter :: CN21=4._ark
    real(ark), parameter :: CN22=100.0109484665696610516028224252161896145_ark
    real(ark), parameter :: CN23=-60.71620503704277071126578458763311203747_ark
    real(ark), parameter :: CN24=388.8627343168314709241605218979822354986_ark
    real(ark), parameter :: CN25=-383.8990705469709002326804327561547655907_ark
    real(ark), parameter :: CN26=0.0000298323868029161667182038731013234360111_ark
    real(ark), parameter :: CN27=0.5334828091991671342480745046574689293102_ark
    real(ark), parameter :: CN28=-0.3180644948481200667624726108533735301737_ark
    real(ark), parameter :: CN29=0.0000277202401875_ark
    real(ark), parameter :: CN30=-9.2400800625e-20_ark
    real(ark), parameter :: CN31=0.02053351125_ark
    real(ark), parameter :: CN32=17.174_ark
    real(ark), parameter :: CN33=0.021393_ark
    real(ark), parameter :: CN34=-41452.71788567563212414161443556312221067_ark
    real(ark), parameter :: CN35=-161176.5258608284993572219082417499530222_ark
    real(ark), parameter :: CN36=25165.76191986292099237035328670933373697_ark
    real(ark), parameter :: CN37=159119.1775696045913661293113167018735401_ark
    real(ark), parameter :: CN38=-221.1189147354489532116517890437598733494_ark
    real(ark), parameter :: CN39=-0.0123649813641253706461200655735429625852_ark
    real(ark), parameter :: CN40=131.8319441675551974004921074239149252633_ark
    real(ark), parameter :: CN41=0.0000715000822090270645027032973382854651126_ark
    real(ark), parameter :: CN42=-3._ark
    real(ark), parameter :: CN43=0.119165673313929925434131141880234375_ark
    real(ark), parameter :: CN44=41168.4451789018309118248229296875_ark
    real(ark), parameter :: CN45=0.017553688831267157994525633085702734375_ark
    real(ark), parameter :: CN46=567.75945411599664772145378671875_ark
    real(ark), parameter :: CN47=-0.000051629391700037902094252766792046875_ark
    real(ark), parameter :: CN48=117.04942623211445142658063828125_ark
    real(ark), parameter :: CN49=7573.8195349076262370886806640625_ark
    real(ark), parameter :: CN50=-2.2224979752794988604356739636167890625e-7_ark
    real(ark), parameter :: CN51=0.09299472346854036948736972265625_ark
    real(ark), parameter :: CN52=-0.1207451092193781029766260836283865234375_ark
    real(ark), parameter :: CN53=-9563.771226524715553976572359375_ark
    real(ark), parameter :: CN54=-0.0001475082248812365704618694539903823632812_ark
    real(ark), parameter :: CN55=-15.19666160111484347194203234375_ark
    real(ark), parameter :: CN56=-0.031044284623568379438039088578515625_ark
    real(ark), parameter :: CN57=-33403.933611101110831174430859375_ark
    real(ark), parameter :: CN58=-1.43215377623407358212099822632421875e-11_ark
    real(ark), parameter :: CN59=1.2298426081724144087426434190944921875e-7_ark
    real(ark), parameter :: CN60=13817.57262855854404138053814518770740356_ark
    real(ark), parameter :: CN61=-8388.587306620973664123451095569777912323_ark
    real(ark), parameter :: CN62=53725.50862027616645240730274724998434073_ark
    real(ark), parameter :: CN63=-53039.72585653486378870977043890062451335_ark
    real(ark), parameter :: CN64=-43.94398138918506580016403580797164175445_ark
    real(ark), parameter :: CN65=73.70630491181631773721726301458662444981_ark
    real(ark), parameter :: CN66=0.004121660454708456882040021857847654195067_ark
    real(ark), parameter :: CN67=87.03572398871614851833789178234549224944_ark
    !
    real(ark) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
    real(ark) :: t11, t12, t13, t14, t15, t16, t17, t18, t19, t20
    real(ark) :: t21, t22, t23, t24
    !
    T1 = SQRT(t)
    T8 = T1**2
    T5 = T1*T8
    T3 = t**2
    T4 = CN8 + CN10*t + CN11*T1 + CN9*T3
    T18 = T4**2
    T20 = 1/T3
    T12 = T18*T20*v
    T2 = CN1 + CN2*t + CN4*T1 + CN3*T3
    T13 = ((CN16 + CN17*t)*T2*(CN13 + CN15*T1 + CN14*T5))/T5
    T14 = CN32 + CN33*t
    T6 = CN19 + CN20*t
    T7 = CN7*T12 + CN12*T13 + CN18*T6**3
    T9 = T8**2
    T11 = T1*T9
    T15 = T5*T9
    T10 = CN22 + CN24*t + CN25*T1 + CN27*T11 + CN26*T15 + CN28*T3 + CN23*T5
    T22 = SQRT(CN6*T7**2 + (CN21*T10**3)/(T11*T9))
    T16 = (CN29*T12 + CN30*T13 + CN31*T14**3 + T22)**CTHIRD
    T17 = CN41*(CN16 + CN17*t)*T5 + CN5*T16*T5 + (CN34 + CN35*t + CN37*T1 + CN38*T11 + CN39*T15 + CN40*T3 + CN36*T5)/T16
    T24 = T17**2
    p = (T2**3*T5)/(T17*T24)
    if (.not.present(dpdv)) return
    T19 = T18*T4
    T23 = T16**2
    T21 = T3**2
    dpdv = (CN42*T19**2*(CN60 + CN62*t + CN63*T1 + CN65*T11 + CN66*T15 + CN64*T3 + CN61*T5 + CN67*T23*T5)* (CN29 &
         + (T20*(CN59*T15 + CN58*t*T21 + CN49*v + t*(CN43 + CN44*v) + T3*(CN45 + CN46*v) + T11*(CN47 + CN48*v) &
         + T21*(CN50 + CN51*v) + T5*(CN52 + CN53*v) + t*T3*(CN54 + CN55*v) + T1*(CN56 + CN57*v)))/ T22))/(T1*T2*T23**2*T24**2)
  end subroutine h2_tv_p
  !
  !  Helmholtz free energy, from temperature and pressure and its derivatives
  !
  subroutine h2_tp_f(t,p,f,dfdp)
    real(rk), intent(in)            :: t ! Temperature, Kelvin
    real(rk), intent(in)            :: p ! Pressure, Pascal
    real(rk), intent(out)           :: f ! Helmholtz free energy, Joule/mole
    real(rk), intent(out), optional :: dfdp ! d F/d P
    !
    real(ark), parameter :: CTHIRD=0.3333333333333333333333333333333333333333_ark
    real(ark), parameter :: CN1=27399.86305764357620203130174921905866518_ark
    real(ark), parameter :: CN2=-42048._ark
    real(ark), parameter :: CN3=-0.005239854115198478735697845372104438068214_ark
    real(ark), parameter :: CN4=-135768.121449983925468109389793622537945_ark
    real(ark), parameter :: CN5=0.00855849831814704147248018553682043454391_ark
    real(ark), parameter :: CN6=64025.9380807446334888233334692566914324_ark
    real(ark), parameter :: CN7=-7.465275892889682383675227971991060776001_ark
    real(ark), parameter :: CN8=-0.5615_ark
    real(ark), parameter :: CN9=81.95733333333333333333333333333333333333_ark
    real(ark), parameter :: CN10=-910.769895_ark
    real(ark), parameter :: CN11=11.14676325_ark
    real(ark), parameter :: CN12=7015.033725_ark
    real(ark), parameter :: CN13=-3181.098375_ark
    real(ark), parameter :: CN14=31638._ark
    real(ark), parameter :: CN15=-1259.3_ark
    real(ark), parameter :: CN16=-0.009299210852252764366715043205124302036858_ark
    real(ark), parameter :: CN17=7.692998147777721149260909774258764605422e-7_ark
    real(ark), parameter :: CN18=7463.027933482398908813233535242537181767_ark
    real(ark), parameter :: CN19=0.005705665545431360981653457024546956362607_ark
    real(ark), parameter :: CN20=-673.2346989378807615490114710376187195639_ark
    real(ark), parameter :: CN21=237.1023533217895100671546352053939398138_ark
    real(ark), parameter :: CN22=-0.003099736950750921455571681068374767345619_ark
    real(ark), parameter :: CN23=5.128665431851814099507273182839176403615e-7_ark
    real(ark), parameter :: CN24=-0.003493236076798985823798563581402958712143_ark
    !
    real(ark) :: t1, t2, t3, t4, t5, t6
    !
    T1 = SQRT(t)
    T2 = p**CTHIRD
    T3 = 1/T1
    T4 = T2**2
    T5 = 1/t
    T6 = t**2
    f = CN1 + CN18*t + CN2*T1 + CN9*T1**3 + CN7*T2 + CN16*t*T2 + CN4*T3 + CN5*T4 + CN17*t*T4 + CN3*T3*T4 + CN6*T5 + CN8*T6 &
      + (CN10 + CN11*t + CN12*T3 + CN13*T5)*LOG(p) + (CN14 + CN15*t)*LOG(t)
    if (.not.present(dfdp)) return
    dfdp = (T5*(CN13 + CN19*t*(CN20 + T2)*(CN21 + T2) + T1*(CN12 + CN24*T4) + (CN11 + CN22*T2 + CN23*T4)*T6))/p
  end subroutine h2_tp_f
end module h2_thermo
!
!program test_thermo
!  use accuracy
!  use h2_thermo
!  real(rk) :: t_ref, v_ref, p_ref, f_ref, dpdv_ref, dvdp_ref, dfdp_ref
!  real(rk) :: v, p, dpdv, dvdp, f, dfdp
!  !
!  do 
!    read(5,*,end=100) t_ref, v_ref, p_ref, f_ref, dvdp_ref, dpdv_ref, dfdp_ref
!    call h2_tv_p(t_ref,v_ref,p,dpdv)
!    call h2_tp_v(t_ref,p_ref,v,dvdp)
!    call h2_tp_f(t_ref,p_ref,f,dfdp)
!    call compare('P    ',p_ref,   p)
!    call compare('dP/dV',dpdv_ref,dpdv)
!    call compare('V    ',v_ref,   v)
!    call compare('dV/dP',dvdp_ref,dvdp)
!    call compare('F    ',f_ref,   f)
!    call compare('dF/dP',dfdp_ref,dfdp)
!  end do
!  100 continue
!  !
!  contains
!    subroutine compare(lbl,ref,val)
!      character(len=*), intent(in) :: lbl
!      real(rk), intent(in)         :: ref, val
!      !
!      if (abs(ref-val)<=max(1e-10_rk*abs(ref),100*spacing(abs(ref)))) return
!      write (6,"('ERROR AT T,V=',2g14.7,1x,a,': REF=',g20.12,' VAL=',g20.12,' ERR=',2g20.12)") &
!             t_ref, v_ref, lbl, ref, val, val-ref, abs((val/ref-1.0_rk))
!    end subroutine compare
!end program test_thermo
