!
!  Collection of miscellaneous atomic data. The values are lifted from gamess_dx.c,
!  and do not represent anything in particular.
!
  module atoms
    use accuracy
    implicit none
    private
    public AtomCovalentR, AtomNuclearCharge, AtomElementSymbol
    !
    type element
      character(len=5) :: symbol ! Element symbol
      real(rk)         :: radius ! Element covalent radius, in Angstrom (!)
      real(rk)         :: rgb(3) ! An arbitrary color assignement
    end type element
    !
    !  Index within this table must coincied with the nuclear charge!
    !
    type(element), parameter :: elements(0:103) = (/ &
                  element("Xx",1.000,(/0.5,0.5,0.5/)),        &
                  element("H ",0.354,(/1.0,1.0,1.0/)),        &
                  element("He",0.849,(/1.0,0.8,0.8/)),        &
                  element("Li",1.336,(/0.7,0.1,0.1/)),        &
                  element("Be",1.074,(/0.1,0.1,0.1/)),        &
                  element("B ",0.838,(/0.0,1.0,0.0/)),        &
                  element("C ",0.757,(/0.4,0.5,0.6/)),        &
                  element("N ",0.700,(/0.0,0.0,0.9/)),        &
                  element("O ",0.680,(/1.0,0.0,0.0/)),        &
                  element("F ",0.668,(/0.9,0.6,0.1/)),        &
                  element("Ne",0.920,(/0.1,0.1,0.1/)),        &
                  element("Na",1.539,(/0.0,0.0,1.0/)),        &
                  element("Mg",1.421,(/0.1,0.1,0.1/)),        &
                  element("Al",1.244,(/0.1,0.5,0.1/)),        &
                  element("Si",1.117,(/0.9,0.6,0.1/)),        &
                  element("P ",1.101,(/1.0,0.6,0.0/)),        &
                  element("S ",1.077,(/1.0,1.0,0.0/)),        &
                  element("Cl",1.044,(/0.0,1.0,0.0/)),        &
                  element("Ar",1.032,(/0.0,0.0,0.5/)),        &
                  element("K ",1.953,(/1.0,0.1,0.6/)),        &
                  element("Ca",1.761,(/0.7,0.7,0.7/)),        &
                  element("Sc",1.513,(/0.7,0.7,0.7/)),        &
                  element("Ti",1.412,(/0.7,0.7,0.7/)),        &
                  element("V ",1.402,(/0.7,0.7,0.7/)),        &
                  element("Cr",1.345,(/0.7,0.7,0.7/)),        &
                  element("Mn",1.382,(/0.7,0.7,0.7/)),        &
                  element("Fe",1.335,(/1.0,0.6,0.0/)),        &
                  element("Co",1.241,(/0.6,0.2,0.2/)),        &
                  element("Ni",1.164,(/0.6,0.2,0.2/)),        &
                  element("Cu",1.302,(/0.6,0.2,0.2/)),        &
                  element("Zn",1.193,(/0.6,0.2,0.2/)),        &
                  element("Ga",1.260,(/0.6,0.2,0.2/)),        &
                  element("Ge",1.197,(/0.7,0.6,0.6/)),        &
                  element("As",1.211,(/0.4,0.6,0.6/)),        &
                  element("Se",1.190,(/0.6,0.2,0.2/)),        &
                  element("Br",1.192,(/0.6,0.2,0.2/)),        &
                  element("Kr",1.147,(/0.7,0.7,0.7/)),        &
                  element("Rb",2.260,(/0.0,0.5,0.0/)),        &
                  element("Sr",2.052,(/0.6,0.1,0.9/)),        &
                  element("Y ",1.698,(/0.5,0.5,0.5/)),        &
                  element("Zr",1.564,(/0.5,0.5,0.5/)),        &
                  element("Nb",1.473,(/0.5,0.5,0.5/)),        &
                  element("Mo",1.484,(/0.5,0.5,0.5/)),        &
                  element("Tc",1.322,(/0.5,0.5,0.5/)),        &
                  element("Ru",1.478,(/0.5,0.5,0.5/)),        &
                  element("Rh",1.332,(/0.5,0.5,0.5/)),        &
                  element("Pd",1.338,(/0.5,0.5,0.5/)),        &
                  element("Ag",1.386,(/0.5,0.5,0.5/)),        &
                  element("Cd",1.403,(/0.5,0.5,0.5/)),        &
                  element("In",1.459,(/0.5,0.5,0.5/)),        &
                  element("Sn",1.398,(/0.5,0.5,0.5/)),        &
                  element("Sb",1.407,(/0.5,0.5,0.5/)),        &
                  element("Te",1.386,(/0.5,0.5,0.5/)),        &
                  element("I ",1.382,(/0.6,0.1,0.9/)),        &
                  element("Xe",1.267,(/0.5,0.5,0.5/)),        &
                  element("Cs",2.570,(/0.5,0.5,0.5/)),        &
                  element("Ba",2.277,(/0.5,0.5,0.5/)),        &
                  element("La",1.943,(/0.5,0.5,0.5/)),        &
                  element("Ce",1.841,(/0.5,0.5,0.5/)),        &
                  element("Pr",1.823,(/0.5,0.5,0.5/)),        &
                  element("Nd",1.816,(/0.5,0.5,0.5/)),        &
                  element("Pm",1.801,(/0.5,0.5,0.5/)),        &
                  element("Sm",1.780,(/0.5,0.5,0.5/)),        &
                  element("Eu",1.771,(/0.5,0.5,0.5/)),        &
                  element("Gd",1.735,(/0.5,0.5,0.5/)),        &
                  element("Tb",1.732,(/0.5,0.5,0.5/)),        &
                  element("Dy",1.710,(/0.5,0.5,0.5/)),        &
                  element("Ho",1.696,(/0.5,0.5,0.5/)),        &
                  element("Er",1.673,(/0.5,0.5,0.5/)),        &
                  element("Tm",1.660,(/0.5,0.5,0.5/)),        &
                  element("Yb",1.637,(/0.5,0.5,0.5/)),        &
                  element("Lu",1.671,(/0.5,0.5,0.5/)),        &
                  element("Hf",1.611,(/0.5,0.5,0.5/)),        &
                  element("Ta",1.511,(/0.5,0.5,0.5/)),        &
                  element("W ",1.526,(/0.5,0.5,0.5/)),        &
                  element("Re",1.372,(/0.5,0.5,0.5/)),        &
                  element("Os",1.372,(/0.5,0.5,0.5/)),        &
                  element("Ir",1.371,(/0.5,0.5,0.5/)),        &
                  element("Pt",1.364,(/0.5,0.5,0.5/)),        &
                  element("Au",1.262,(/0.9,0.6,0.1/)),        &
                  element("Hg",1.340,(/0.5,0.5,0.5/)),        &
                  element("Tl",1.518,(/0.5,0.5,0.5/)),        &
                  element("Pb",1.459,(/0.5,0.5,0.5/)),        &
                  element("Bi",1.512,(/0.5,0.5,0.5/)),        &
                  element("Po",1.500,(/0.5,0.5,0.5/)),        &
                  element("At",1.545,(/0.5,0.5,0.5/)),        &
                  element("Rn",1.420,(/0.5,0.5,0.5/)),        &
                  element("Fr",2.880,(/0.5,0.5,0.5/)),        &
                  element("Ra",2.512,(/0.5,0.5,0.5/)),        &
                  element("Ac",1.983,(/0.5,0.5,0.5/)),        &
                  element("Th",1.721,(/0.5,0.5,0.5/)),        &
                  element("Pa",1.711,(/0.5,0.5,0.5/)),        &
                  element("U ",1.684,(/0.5,0.5,0.5/)),        &
                  element("Np",1.666,(/0.5,0.5,0.5/)),        &
                  element("Pu",1.657,(/0.5,0.5,0.5/)),        &
                  element("Am",1.660,(/0.5,0.5,0.5/)),        &
                  element("Cm",1.801,(/0.5,0.5,0.5/)),        &
                  element("Bk",1.761,(/0.5,0.5,0.5/)),        &
                  element("Cf",1.750,(/0.5,0.5,0.5/)),        &
                  element("Es",1.724,(/0.5,0.5,0.5/)),        &
                  element("Fm",1.712,(/0.5,0.5,0.5/)),        &
                  element("Md",1.689,(/0.5,0.5,0.5/)),        &
                  element("No",1.679,(/0.5,0.5,0.5/)),        &
                  element("Lr",1.698,(/0.5,0.5,0.5/)) /)
    !
    contains
    !
    !  Externally visible interfaces
    !
    function AtomCovalentR(atom) result(r)
      character(len=*), intent(in) :: atom ! Element symbol
      real(rk)                     :: r    ! Covalent radius, or -1 if not found
      !
      integer(ik) :: ind
      !
      ind = find_atom(atom)
      r   = -1._rk
      if (ind>=0) r = elements(ind)%radius
    end function AtomCovalentR
    !
    function AtomNuclearCharge(atom) result(q)
      character(len=*), intent(in) :: atom ! Element symbol
      real(rk)                     :: q    ! Nuclear charge, in |e| units
      !
      q = find_atom(atom)
    end function AtomNuclearCharge
    !
    function AtomElementSymbol(q) result(atom)
      real(rk), intent(in) :: q    ! Nuclear charge, in |e| units
      character(len=7)     :: atom ! Element symbol
      !
      integer(ik) :: iq
      !
      iq = nint(q)
      if (iq<lbound(elements,1) .or. iq>ubound(elements,1)) then
        atom = 'Unknown'
      else
        atom = elements(iq)%symbol
      end if
    end function AtomElementSymbol
    !
    !  Internal routines beyond this point
    !
    function find_atom(atom) result(ind)
      character(len=*), intent(in) :: atom ! Element to look up; case is not significant
      integer(ik)                  :: ind  ! Index in elements, or -1 if not found
      !
      character(len=3) :: la
      integer          :: ia    ! Must be of default integer kind - iachar has problems otherwise
      !
      la = atom
      !
      !  Capitalize first character; lowecase the remaining two.
      !
      ia = iachar(la(1:1))
      if (ia>=iachar('a') .and. ia<=iachar('z')) then
        la(1:1) = achar(ia-iachar('a')+iachar('A'))
      end if
      ia = iachar(la(2:2))
      if (ia>=iachar('A') .and. ia<=iachar('Z')) then
        la(2:2) = achar(ia-iachar('A')+iachar('a'))
      end if
      ia = iachar(la(3:3))
      if (ia>=iachar('A') .and. ia<=iachar('Z')) then
        la(3:3) = achar(ia-iachar('A')+iachar('a'))
      end if
      !
      !  Scan the table
      !
      scan_atoms: do ind=lbound(elements,dim=1),ubound(elements,dim=1)
        if (elements(ind)%symbol==la) return
      end do scan_atoms
      !
      !  Not found!
      !
      ind = -1
    end function find_atom
  end module atoms
