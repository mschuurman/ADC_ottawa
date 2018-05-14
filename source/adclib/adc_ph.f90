module adc_ph

!!$The module contains an assortment of functions needed to calculate
!!$matrix elements of the ADC matrix of the polarization propagator.
!!$For the spin-orbital expressions see A.B. Trofimov et al, JCP 111,9982 (1999).
!!$Spin free expressions were taken from the Ph.D. thesis of A.B. Trofimov.
  
  use constants
  use parameters
  use misc
  use channels
  use vpqrsmod
  
  implicit none
  
contains

!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************PH-PH BLOCK***********************************
!!$*******************************************************************************
!!$*******************************************************************************

!#######################################################################

!!$Indices are supplied in the order: PH,PH

!!$Zeroth order contribution K_ak,a'k'. The condition that a=a', k=k'
!!$is checked by the calling procedure

  REAL(DP) FUNCTION K_PH_PH(E_A,E_K)
    
    REAL(DP), INTENT(IN) :: E_A,E_K
    
    K_PH_PH=E_A-E_K
    
  end function K_ph_ph

!#######################################################################

!!$First order contribution C_ak,a'k'
  
  real(dp) function C1_ph_ph(a,k,a1,k1)

    implicit none
    
    integer, intent(in) :: a,k,a1,k1

    C1_ph_ph=2.0d0*vpqrs(a,k,a1,k1)-vpqrs(a,a1,k,k1)

    return
    
  end function C1_ph_ph

!!#######################################################################
!
!!!$First order contribution C_ak,a'k'
!
!  real(dp) function C1a_ph_ph(a,k)
!
!    implicit none
!
!!    integer, intent(in) :: a,k,a1,k1,i,c,r,s,nsym1,nsym2
!    integer :: a,k,a1,k1,i,c,r,s,nsym1,nsym2
!    real(dp)             :: term
!
!    C1a_ph_ph=0._dp
!    do c=nOcc+1,nBas
!       r=roccnum(c)
!       do i=1,nOcc
!          s=roccnum(i)
!
!          nsym1=MT(orbSym(r),orbSym(s))
!          nsym2=MT(orbSym(k),orbSym(a))
!
!          if(MT(nsym1,nsym2) .eq. 1) then
!
!             term=0._dp
!
!             term=term+2.0d0*vpqrs(a,k,r,s)-vpqrs(s,k,r,a)
!
!
!             C1a_ph_ph=C1a_ph_ph+term
!
!          end if
!       end do
!    end do
!
!    return
!
!  end function C1a_ph_ph
!
!!#######################################################################

!!$Second order contribution CA_ak,a'k'. The condition that k=k' is
!!$checked in the calling procedure.

  real(dp) function Ca1_ph_ph(a,a1)

    integer, intent(in) :: a,a1

    integer  :: c,i,j, nsym1,nsym2,nsym3,u,v,r,cnt
    real(dp) :: DA,eija,eijc,term

    Ca1_ph_ph=0._dp
    cnt=0
    do c=nOcc+1,nBas
       r=roccnum(c)
       do i=1,nOcc
          u=roccnum(i)
          do j=1,nOcc
             v=roccnum(j)

             nsym1=MT(orbSym(u),orbSym(v))
             nsym2=MT(orbSym(a),orbSym(r))
             nsym3=MT(orbSym(a1),orbSym(r))


             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

                cnt=cnt+1
                term=0._dp

!                eijc=e(u)+e(v)-e(r)
!                DA=(eijc-0.5_dp*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))
                eija=e(r)+e(a)-e(u)-e(v)
                eijc=e(r)+e(a1)-e(u)-e(v)
                DA=1._dp/((e(r)+e(a)-e(u)-e(v))*(e(r)+e(a1)-e(u)-e(v)))

!                term=term+eijc*vpqrs(a1,u,r,v)*(2._dp*vpqrs(v,a,u,r)-vpqrs(v,r,u,a))
!                term=term+eija*vpqrs(a,u,r,v)*(2._dp*vpqrs(v,a1,u,r)-vpqrs(v,r,u,a1))
                term=term+eijc*vpqrs(a1,u,r,v)*(2._dp*vpqrs(u,a,v,r)-vpqrs(u,r,v,a))
                term=term+eija*vpqrs(a,v,r,u)*(2._dp*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r))

                term=DA*term

                Ca1_ph_ph=Ca1_ph_ph+term

             end if
          end do
       end do
    end do

    Ca1_ph_ph=0.5_dp*Ca1_ph_ph

  end function Ca1_ph_ph

!#######################################################################
             
!!$Second order contribution CB_ak,a'k'. The condition that a=a' is
!!$checked in the calling procedure.

  real(dp) function Cb1_ph_ph(k,k1)

    integer, intent(in) :: k,k1
    
    integer  :: c,dd,i, nsym1, nsym2, nsym3,u,r,s
    real(dp) :: DB,eicd,ejcd,term

!!$    CB_ph_ph=0._dp             
  
!!$! Faster version (taking symmetries into account)

    Cb1_ph_ph=0._dp
    do c=nOcc+1,nBas
       r=roccnum(c)
       do dd=nOcc+1,nBas
          s=roccnum(dd)
          do i=1,nOcc
             u=roccnum(i)

             nsym1=MT(orbSym(roccnum(c)),orbSym(roccnum(dd)))
             nsym2=MT(orbSym(k),orbSym(roccnum(i)))
             nsym3=MT(orbSym(k1),orbSym(roccnum(i)))
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
            
                term=0._dp
                
!                eicd=e(roccnum(i))-e(roccnum(c))-e(roccnum(dd))
!                DB=(eicd+0.5_dp*(e(k)+e(k1)))/((eicd+e(k))*(eicd+e(k1)))

                eicd=e(r)+e(s)-e(u)-e(k)
                ejcd=e(r)+e(s)-e(u)-e(k1)
                DB=1._dp/((e(r)+e(s)-e(u)-e(k))*(e(r)+e(s)-e(u)-e(k1)))

!                term=term+ejcd*vpqrs(k1,r,u,s)*(2._dp*vpqrs(r,k,s,u)-vpqrs(s,k,r,u))
!                term=term+eicd*vpqrs(k,r,u,s)*(2._dp*vpqrs(r,k1,s,u)-vpqrs(s,k1,r,u))
                term=term+ejcd*vpqrs(r,k1,s,u)*(2._dp*vpqrs(k,r,u,s)-vpqrs(k,s,u,r))
                term=term+eicd*vpqrs(r,u,s,k)*(2._dp*vpqrs(k1,s,u,r)-vpqrs(k1,r,u,s))
                
                term=DB*term
                
                Cb1_ph_ph=Cb1_ph_ph+term
                
             end if
             
          end do
       end do
    end do
    
    Cb1_ph_ph=0.5_dp*Cb1_ph_ph   
   
  end function Cb1_ph_ph

!#######################################################################

!!$Second order contribution CC_ak,a'k'.

  real(dp) function Cc1_ph_ph(a,k,a1,k1)

    integer, intent(in) :: a,a1,k,k1

    integer  :: c,i, nsym1, nsym2, nsym3,r,s
    real(dp) :: DC,eic,term

    real(dp), dimension(nocc) :: tau

    Cc1_ph_ph=0._dp

    do c=nOcc+1,nBas
       r=roccnum(c)
       do i=1,nOcc
          s=roccnum(i)


          nsym1=MT(orbSym(r),orbSym(s))
          nsym2=MT(orbSym(k),orbSym(a))
          nsym3=MT(orbSym(k1),orbSym(a1))


          if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

             term=0._dp

!             eic=e(s)-e(r)
!             DC=(0.5_dp*(e(k)+e(k1)-e(a)-e(a1))+eic)/((e(k)-e(a)+eic)*(e(k1)-e(a1)+eic))

             DC=1._dp/(e(r)+e(a)-e(s)-e(k))

!             term=(2._dp*vpqrs(k,a,s,r)-vpqrs(k,r,s,a))*(2._dp*vpqrs(a1,k1,r,s)-vpqrs(r,k1,a1,s))
             term=(2._dp*vpqrs(a,k,r,s)-vpqrs(a,s,r,k))*(2._dp*vpqrs(a1,k1,r,s)-vpqrs(a1,s,r,k1))
             term=term*DC

             Cc1_ph_ph=Cc1_ph_ph+term

          end if

       end do
    end do

    Cc1_ph_ph=-0.5_dp*Cc1_ph_ph

  end function Cc1_ph_ph

!#######################################################################

!!$Second order contribution CC_ak,a'k'.

  real(dp) function Cc2_ph_ph(a,k,a1,k1)

    integer, intent(in) :: a,a1,k,k1

    integer  :: c,i, nsym1, nsym2, nsym3,r,s
    real(dp) :: DC,eic,term

    real(dp), dimension(nocc) :: tau

    Cc2_ph_ph=0._dp

    do c=nOcc+1,nBas
       r=roccnum(c)
       do i=1,nOcc
          s=roccnum(i)


          nsym1=MT(orbSym(r),orbSym(s))
          nsym2=MT(orbSym(k),orbSym(a))
          nsym3=MT(orbSym(k1),orbSym(a1))


          if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

             term=0._dp

!             eic=e(s)-e(r)
!             DC=(0.5_dp*(e(k)+e(k1)-e(a)-e(a1))+eic)/((e(k)-e(a)+eic)*(e(k1)-e(a1)+eic))

             DC=1._dp/(e(r)+e(a1)-e(s)-e(k1))

             term=(2._dp*vpqrs(a,k,r,s)-vpqrs(a,s,r,k))*(2._dp*vpqrs(a1,k1,r,s)-vpqrs(a1,s,r,k1))
             term=term*DC

             Cc2_ph_ph=Cc2_ph_ph+term

          end if

       end do
    end do

    Cc2_ph_ph=-0.5_dp*Cc2_ph_ph

  end function Cc2_ph_ph

!#######################################################################

!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************PH-2PH BLOCK**********************************
!!$*******************************************************************************
!!$*******************************************************************************

!!$  We distinguish here between five different types of coupling. 
!!$ Calculating Cak,a'b'k'l'

!#######################################################################
  
!!$ a'|=b' and k'|=l'; spin case 1
  
  function C11_ph_2p2h(a,k,apr,bpr,kpr,lpr) result(func)
    
    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
    real(dp)            :: func
 
    func=0.0d0

!    if (a.eq.apr) func=func+vpqrs(kpr,k,lpr,bpr)+vpqrs(lpr,k,kpr,bpr)
    if (a.eq.apr) func=func+vpqrs(k,kpr,lpr,bpr)+vpqrs(k,lpr,kpr,bpr)
    
!    if (a.eq.bpr) func=func+vpqrs(kpr,k,lpr,apr)+vpqrs(lpr,k,kpr,apr)
    if (a.eq.bpr) func=func+vpqrs(k,kpr,lpr,apr)+vpqrs(k,lpr,kpr,apr)

!    if (k.eq.kpr) func=func-(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))
    if (k.eq.kpr) func=func-(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))

!    if (k.eq.lpr) func=func-(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
    if (k.eq.lpr) func=func-(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
    func=func/sqrt(2.0d0)
    
    return

  end function C11_ph_2p2h

!#######################################################################

!!$ a'|=b' and k'|=l'; spin case 2
  
  function C22_ph_2p2h(a,k,apr,bpr,kpr,lpr) result(func)
    
    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
    real(dp)            :: func

    func=0.0d0

!    if (a.eq.apr) func=func+(vpqrs(kpr,k,lpr,bpr)+vpqrs(kpr,bpr,lpr,k))
    if (a.eq.apr) func=func+(vpqrs(k,kpr,lpr,bpr)+vpqrs(k,lpr,kpr,bpr)) 

    if (a.eq.bpr) then
!       func=func-(vpqrs(kpr,k,lpr,apr)+vpqrs(kpr,apr,lpr,k))
       func=func-(vpqrs(k,kpr,lpr,apr)+vpqrs(k,lpr,kpr,apr))
!       func=sqrt(2.0d0)*func/3.0d0
       func=func/3.0d0
    endif

!    if (k.eq.kpr) func=func-(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))
    if (k.eq.kpr) func=func-(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))

    if (k.eq.lpr) then
!       func=func+(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
       func=func+(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
!       func=sqrt(2.0d0)*func/3.0d0
       func=func/3.0d0
    endif

    func=sqrt(3.0d0)*func/sqrt(2.0d0)
   
  end function C22_ph_2p2h
  
!#######################################################################

!!$ a'=b' and k'|=l' 
  
  function C33_ph_2p2h(j,k,ipr,kpr,lpr) result(func)
    
    integer, intent(in) :: j,k,ipr,kpr,lpr
    real(dp)            :: func

    func=0.0d0
    
    if (j.eq.ipr) then
!       func=func+(vpqrs(kpr,k,lpr,ipr)+vpqrs(kpr,ipr,lpr,k))
       func=func+(vpqrs(k,kpr,lpr,ipr)+vpqrs(k,lpr,kpr,ipr))
       func=func/sqrt(2.0d0)
    endif

!    if (k.eq.kpr) func=func-vpqrs(j,ipr,lpr,ipr)
    if (k.eq.kpr) func=func-vpqrs(j,ipr,lpr,ipr)
    
!    if (k.eq.lpr) func=func-vpqrs(j,ipr,kpr,ipr)
    if (k.eq.lpr) func=func-vpqrs(j,ipr,kpr,ipr)

!    func=func/sqrt(2.0d0)
    func=sqrt(2.0d0)*func

  end function C33_ph_2p2h

!#######################################################################

!!$ a'|=b' and k'=l'

  function C44_ph_2p2h(a,k,apr,bpr,kpr) result(func)
    
    integer, intent(in) :: a,k,apr,bpr,kpr
    real(dp)            :: func
  
    func=0.0d0

!    if (a.eq.apr) func=func+vpqrs(kpr,k,kpr,bpr)
    if (a.eq.apr) func=func+vpqrs(k,kpr,kpr,bpr)

!    if (a.eq.bpr) func=func+vpqrs(kpr,k,kpr,apr)
    if (a.eq.bpr) func=func+vpqrs(k,kpr,kpr,apr)

    if (k.eq.kpr) then
!       func=func-(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
       func=func-(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
       func=func/sqrt(2.0d0)
    endif

!    func=func/sqrt(2.0d0)
    func=sqrt(2.0d0)*func

  end function C44_ph_2p2h

!#######################################################################

!!$ a'=b' and k'=l'
  
  function C55_ph_2p2h(a,k,apr,kpr) result(func)
    
    integer, intent(in) :: a,k,apr,kpr
    real(dp)            :: func

    func=0.0d0
     
!    if (a.eq.apr) func=func+vpqrs(kpr,apr,kpr,k)
    if (a.eq.apr) func=func+vpqrs(k,kpr,kpr,apr)

!    if (k.eq.kpr) func=func-vpqrs(a,apr,kpr,apr)
    if (k.eq.kpr) func=func-vpqrs(a,apr,kpr,apr)
    
    func=sqrt(2.0d0)*func
   
  end function C55_ph_2p2h

!#######################################################################
  
!!$ a'|=b' and k'|=l'; spin case 1
  
  function C11_1_ph_2p2h(a,k,apr,bpr,kpr,lpr) result(func)
    
    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
    real(dp)            :: func
 
    func=0.0d0

!    if (a.eq.apr) func=func+vpqrs(kpr,k,lpr,bpr)+vpqrs(lpr,k,kpr,bpr)
    !if (a.eq.apr) 
    func=func+vpqrs(k,kpr,lpr,bpr)!+vpqrs(k,lpr,kpr,bpr)
    
!    if (a.eq.bpr) func=func+vpqrs(kpr,k,lpr,apr)+vpqrs(lpr,k,kpr,apr)
    !if (a.eq.bpr) 
    func=func+vpqrs(k,kpr,lpr,apr)!+vpqrs(k,lpr,kpr,apr)

!    if (k.eq.kpr) func=func-(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))
    !if (k.eq.kpr) 
    func=func-vpqrs(a,apr,lpr,bpr)!+vpqrs(a,bpr,lpr,apr)

!    if (k.eq.lpr) func=func-(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
    !if (k.eq.lpr) 
    func=func-vpqrs(a,apr,kpr,bpr)!+vpqrs(a,bpr,kpr,apr)

    !func=func/sqrt(2.0d0)
    
    return

  end function C11_1_ph_2p2h

!#######################################################################

!!$ a'|=b' and k'|=l'; spin case 2
  
  function C22_1_ph_2p2h(a,k,apr,bpr,kpr,lpr) result(func)
    
    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
    real(dp)            :: func

    func=0.0d0

!    if (a.eq.apr) func=func+(vpqrs(kpr,k,lpr,bpr)+vpqrs(kpr,bpr,lpr,k))
    !if (a.eq.apr) 
    func=func+3.0d0*(vpqrs(k,kpr,lpr,bpr))!+vpqrs(k,lpr,kpr,bpr)) 

    !if (a.eq.bpr) then
!       func=func-(vpqrs(kpr,k,lpr,apr)+vpqrs(kpr,apr,lpr,k))
       func=func+(vpqrs(k,kpr,lpr,apr))!+vpqrs(k,lpr,kpr,apr))
!       func=func/3.0d0
    !endif

!    if (k.eq.kpr) func=func-(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))
    !if (k.eq.kpr) 
    func=func-3.0d0*(vpqrs(a,apr,lpr,bpr))!+vpqrs(a,bpr,lpr,apr))

    !if (k.eq.lpr) then
!       func=func+(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
       func=func+(vpqrs(a,apr,kpr,bpr))!+vpqrs(a,bpr,kpr,apr))
!       func=func/3.0d0
    !endif

!    func=sqrt(3.0d0)*func!/sqrt(2.0d0)
    func=2.0d0*func/sqrt(12.0d0)
   
  end function C22_1_ph_2p2h
  
!#######################################################################

!!$ a'=b' and k'|=l' 
  
  function C33_1_ph_2p2h(j,k,ipr,kpr,lpr) result(func)
    
    integer, intent(in) :: j,k,ipr,kpr,lpr
    real(dp)            :: func

    func=0.0d0
    
    !if (j.eq.ipr) then
!       func=func+(vpqrs(kpr,k,lpr,ipr)+vpqrs(kpr,ipr,lpr,k))
       func=func+(vpqrs(k,kpr,lpr,ipr))!+vpqrs(k,lpr,kpr,ipr))
!       func=func/sqrt(2.0d0)
    !endif

!    if (k.eq.kpr) func=func-vpqrs(j,ipr,lpr,ipr)
    !if (k.eq.kpr) 
    func=func-vpqrs(j,ipr,lpr,ipr)
    
!    if (k.eq.lpr) func=func-vpqrs(j,ipr,kpr,ipr)
    !if (k.eq.lpr) 
    func=func-vpqrs(j,ipr,kpr,ipr)

!    func=func/sqrt(2.0d0)
    func=sqrt(2.0d0)*func

  end function C33_1_ph_2p2h

!#######################################################################

!!$ a'|=b' and k'=l'

  function C44_1_ph_2p2h(a,k,apr,bpr,kpr) result(func)
    
    integer, intent(in) :: a,k,apr,bpr,kpr
    real(dp)            :: func
  
    func=0.0d0

!    if (a.eq.apr) func=func+vpqrs(kpr,k,kpr,bpr)
    !if (a.eq.apr) 
    func=func+vpqrs(k,kpr,kpr,bpr)

!    if (a.eq.bpr) func=func+vpqrs(kpr,k,kpr,apr)
    !if (a.eq.bpr) 
    func=func+vpqrs(k,kpr,kpr,apr)

    !if (k.eq.kpr) then
!       func=func-(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
       func=func-(vpqrs(a,apr,kpr,bpr))!+vpqrs(a,bpr,kpr,apr))
!       func=func/sqrt(2.0d0)
    !endif

!    func=func/sqrt(2.0d0)
    func=sqrt(2.0d0)*func

  end function C44_1_ph_2p2h

!#######################################################################

!!$ a'=b' and k'=l'
  
  function C55_1_ph_2p2h(a,k,apr,kpr) result(func)
    
    integer, intent(in) :: a,k,apr,kpr
    real(dp)            :: func

    func=0.0d0
     
!    if (a.eq.apr) func=func+vpqrs(kpr,apr,kpr,k)
    !if (a.eq.apr) 
    func=func+vpqrs(k,kpr,kpr,apr)

!    if (k.eq.kpr) func=func-vpqrs(a,apr,kpr,apr)
    !if (k.eq.kpr) 
    func=func+vpqrs(a,apr,kpr,apr)
    
    func=sqrt(2.0d0)*func
   
  end function C55_1_ph_2p2h

!!#######################################################################
! IS
!!!$ a'|=b' and k'|=l'; spin case 1
!
!  function C11_ph_2p2h(a,k,apr,bpr,kpr,lpr) result(func)
!
!    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
!    real(dp)             :: func
!
!    func=0.0d0
!
!    if (a.eq.apr) func=func-2._dp*(vpqrs(kpr,k,lpr,bpr)+6._dp*vpqrs(kpr,bpr,lpr,k))
!
!    if (a.eq.bpr) func=func-2._dp*(vpqrs(kpr,k,lpr,apr)+6._dp*vpqrs(kpr,apr,lpr,k))
!
!    if (k.eq.kpr) func=func-2._dp*(vpqrs(a,apr,lpr,bpr)+6._dp*vpqrs(a,bpr,lpr,apr))
!
!    if (k.eq.lpr) func=func-2._dp*(vpqrs(a,apr,kpr,bpr)+6._dp*vpqrs(a,bpr,kpr,apr))
!
!    func=func/sqrt(12.0d0)
!
!    return
!
!  end function C11_ph_2p2h
!
!!#######################################################################
!
!!!$ a'|=b' and k'|=l'; spin case 2
!
!  function C22_ph_2p2h(a,k,apr,bpr,kpr,lpr) result(func)
!
!    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
!    real(dp)            :: func
!
!    func=0.0d0
!
!    if (a.eq.apr) func=func+(vpqrs(kpr,k,lpr,bpr)+vpqrs(kpr,bpr,lpr,k))
!
!    if (a.eq.bpr) func=func+(vpqrs(kpr,k,lpr,apr)+vpqrs(kpr,apr,lpr,k))
!
!    if (k.eq.kpr) func=func+(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))
!
!    if (k.eq.lpr) func=func+(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
!
!!    func=sqrt(3.0d0)*func/sqrt(2.0d0)
!
!  end function C22_ph_2p2h
!
!!#######################################################################
!
!!!$ a'=b' and k'|=l' 
!
!  function C33_ph_2p2h(j,k,ipr,kpr,lpr) result(func)
!
!    integer, intent(in) :: j,k,ipr,kpr,lpr
!    real(dp)            :: func
!
!    func=0.0d0
!
!    if (j.eq.ipr) func=func+vpqrs(kpr,k,lpr,ipr)+vpqrs(kpr,ipr,lpr,k)
!
!    if (k.eq.kpr) func=func+vpqrs(j,ipr,lpr,ipr)+vpqrs(ipr,lpr,ipr,j)
!
!    if (k.eq.lpr) func=func+vpqrs(j,ipr,kpr,ipr)+vpqrs(ipr,kpr,ipr,j)
!
!    func=func/sqrt(2.0d0)
!
!  end function C33_ph_2p2h
!
!!#######################################################################
!
!!!$ a'|=b' and k'=l'
!
!  function C44_ph_2p2h(a,k,apr,bpr,kpr) result(func)
!
!    integer, intent(in) :: a,k,apr,bpr,kpr
!    real(dp)            :: func
!
!    func=0.0d0
!
!    if (a.eq.apr) func=func+vpqrs(kpr,k,kpr,bpr)+vpqrs(kpr,bpr,kpr,k)
!
!    if (a.eq.bpr) func=func+vpqrs(kpr,k,kpr,apr)+vpqrs(kpr,apr,kpr,k)
!
!    if (k.eq.kpr) func=func+vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr)
!
!    func=func/sqrt(2.0d0)
!
!  end function C44_ph_2p2h
!
!!#######################################################################
!
!!!$ a'=b' and k'=l'
!
!  function C55_ph_2p2h(a,k,apr,kpr) result(func)
!
!    integer, intent(in) :: a,k,apr,kpr
!    real(dp)            :: func
!
!    func=0.0d0
!
!    if (a.eq.apr) func=func+vpqrs(kpr,apr,kpr,k)
!
!    if (k.eq.kpr) func=func-vpqrs(a,apr,kpr,apr)
!
!    func=sqrt(2.0d0)*func
!
!  end function C55_ph_2p2h
!
!!#######################################################################
!
!  function C1a_ph_ph(a,k,apr,bpr,kpr,lpr) result(func)
!    
!    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
!    real(dp)            :: func
!
!    func=0.0d0
!
!    if (a.eq.apr) func=func-(vpqrs(kpr,k,lpr,bpr)+vpqrs(kpr,bpr,lpr,k))
!    
!    if (a.eq.bpr) func=func-(vpqrs(kpr,k,lpr,apr)+vpqrs(kpr,apr,lpr,k))
!
!    if (k.eq.kpr) func=func+(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))
!    
!    if (k.eq.lpr) func=func+(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
!
!    func=func/sqrt(2.0d0)
!    
!    return
!
!  end function C1a_ph_ph
!
!!#######################################################################
!
!  function C1b_ph_ph(a,k,apr,bpr,kpr,lpr) result(func)
!    
!    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
!    real(dp)            :: func
!
!    func=0.0d0
!
!    if (a.eq.apr) func=func-(vpqrs(kpr,k,lpr,bpr)+vpqrs(kpr,bpr,lpr,k))
!    
!    if (a.eq.bpr) func=func-(vpqrs(kpr,k,lpr,apr)+vpqrs(kpr,apr,lpr,k))
!
!    if (k.eq.kpr) func=func+(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))
!    
!    if (k.eq.lpr) func=func+(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
!
!    func=func/sqrt(2.0d0)
!    
!    return
!
!  end function C1b_ph_ph
!
!
!!#######################################################################

!!$Second order contribution CA_ak,a'k'. The condition that k=k' is
!!$checked in the calling procedure.

  real(dp) function CA_ph_ph(a,a1)

    integer, intent(in) :: a,a1
    
    integer  :: c,i,j, nsym1,nsym2,nsym3,u,v,r,cnt
    real(dp) :: DA,eijc,term

    CA_ph_ph=0._dp
    cnt=0
    do c=nOcc+1,nBas
       r=roccnum(c)
       do i=1,nOcc
          u=roccnum(i)
          do j=1,nOcc
             v=roccnum(j)

             nsym1=MT(orbSym(u),orbSym(v))
             nsym2=MT(orbSym(a),orbSym(r))
             nsym3=MT(orbSym(a1),orbSym(r))
             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1
                term=0._dp

                eijc=e(u)+e(v)-e(r)
                DA=(eijc-0.5_dp*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

                term=term+vpqrs(a,u,r,v)*(2._dp*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
                term=term+vpqrs(a,v,r,u)*(2._dp*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r))

                term=DA*term
                
                CA_ph_ph=CA_ph_ph+term
                
             end if
          end do
       end do
    end do 
    
    CA_ph_ph=-0.5_dp*CA_ph_ph            
    
  end function CA_ph_ph

!#######################################################################
             
!!$Second order contribution CB_ak,a'k'. The condition that a=a' is
!!$checked in the calling procedure.

  real(dp) function CB_ph_ph(k,k1)

    integer, intent(in) :: k,k1
    
    integer  :: c,dd,i, nsym1, nsym2, nsym3,u,r,s
    real(dp) :: DB,eicd,term

!!$    CB_ph_ph=0._dp
  
!!$! Faster version (taking symmetries into account)

!!$! b|=c
!!$    
!!$    do c=nOcc+1,nBas
!!$       do dd=c+1,nBas
!!$          do i=1,nOcc
!!$            
!!$             term=0._dp
!!$ 
!!$             eicd=e(roccnum(i))-e(roccnum(c))-e(roccnum(dd))
!!$             DB=(eicd+0.5_dp*(e(k)+e(k1)))/((eicd+e(k))*(eicd+e(k1)))
!!$             
!!$             term=term+vpqrs(c,k,dd,i)*(2._dp*vpqrs(k1,c,i,dd)-vpqrs(k1,dd,i,c))
!!$             term=term+vpqrs(c,i,dd,k)*(2._dp*vpqrs(k1,dd,i,c)-vpqrs(k1,c,i,dd))
!!$             term=DB*term
!!$             
!!$             CB_ph_ph=CB_ph_ph+term
!!$             
!!$          end do
!!$       end do
!!$    end do
!!$    
!!$    CB_ph_ph=2._dp*CB_ph_ph
       
!!$! b=c
    
!!$    do c=nOcc+1,nBas
!!$       do i=1,nOcc
!!$          
!!$          term=0._dp
!!$          
!!$          eicd=e(roccnum(i))-2._dp*e(roccnum(c))
!!$          DB=(eicd+0.5_dp*(e(k)+e(k1)))/((eicd+e(k))*(eicd+e(k1)))
!!$
!!$          term=term+2._dp*vpqrs(c,k,c,i)*vpqrs(k1,c,i,c)
!!$          term=DB*term
!!$          
!!$          CB_ph_ph=CB_ph_ph+term
!!$          
!!$       end do
!!$    end do
!!$    
!!$    CB_ph_ph=-0.5_dp*CB_ph_ph

!!$ Dumbed down version

    CB_ph_ph=0._dp   
    do c=nOcc+1,nBas
       r=roccnum(c)
       do dd=nOcc+1,nBas
          s=roccnum(dd)
          do i=1,nOcc
             u=roccnum(i)

             nsym1=MT(orbSym(roccnum(c)),orbSym(roccnum(dd)))
             nsym2=MT(orbSym(k),orbSym(roccnum(i)))
             nsym3=MT(orbSym(k1),orbSym(roccnum(i)))
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
            
                term=0._dp
                
                eicd=e(roccnum(i))-e(roccnum(c))-e(roccnum(dd))
                DB=(eicd+0.5_dp*(e(k)+e(k1)))/((eicd+e(k))*(eicd+e(k1)))
                
                term=term+vpqrs(r,k,s,u)*(2._dp*vpqrs(k1,r,u,s)-vpqrs(k1,s,u,r))
                term=term+vpqrs(r,u,s,k)*(2._dp*vpqrs(k1,s,u,r)-vpqrs(k1,r,u,s))
!!$
!!$                term=term+vpqrs(c,k,dd,i)*(vpqrs(k1,c,i,dd)-2._dp*vpqrs(k1,dd,i,c))
!!$                term=term+vpqrs(c,i,dd,k)*(vpqrs(k1,dd,i,c)-2._dp*vpqrs(k1,c,i,dd))
                
                term=DB*term
                
                CB_ph_ph=CB_ph_ph+term
                
             end if
             
          end do
       end do
    end do
    
    CB_ph_ph=-0.5_dp*CB_ph_ph   
   
  end function CB_ph_ph

!#######################################################################

!!$Second order contribution CC_ak,a'k'.

  real(dp) function CC_ph_ph(a,k,a1,k1)

    integer, intent(in) :: a,a1,k,k1
    
    integer  :: c,i, nsym1, nsym2, nsym3,r,s
    real(dp) :: DC,eic,term
    
    real(dp), dimension(nocc) :: tau

    CC_ph_ph=0._dp

    do c=nOcc+1,nBas
       r=roccnum(c)
       do i=1,nOcc
          s=roccnum(i)

    
          nsym1=MT(orbSym(r),orbSym(s))
          nsym2=MT(orbSym(k),orbSym(a))
          nsym3=MT(orbSym(k1),orbSym(a1))
          
             
          if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
             
             term=0._dp
             
             eic=e(s)-e(r)
             DC=(0.5_dp*(e(k)+e(k1)-e(a)-e(a1))+eic)/((e(k)-e(a)+eic)*(e(k1)-e(a1)+eic))
             
             term=(2._dp*vpqrs(a,k,r,s)-vpqrs(a,s,r,k))*(2._dp*vpqrs(a1,k1,r,s)-vpqrs(a1,s,r,k1))
             term=term*DC

             CC_ph_ph=CC_ph_ph+term
             
          end if
          
       end do
    end do

  end function CC_ph_ph

!#######################################################################

!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************PH-2PH BLOCK**********************************
!!$*******************************************************************************
!!$*******************************************************************************

!!$  We distinguish here between five different types of coupling. 
!!$ Calculating Cak,a'b'k'l'

!#######################################################################
  
!!$ a'|=b' and k'|=l'; spin case 1
  
  function C1_ph_2p2h(a,k,apr,bpr,kpr,lpr) result(func)
    
    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
    real(dp)            :: func
 
    func=0.0d0

    if (a.eq.apr) func=func-(vpqrs(kpr,k,lpr,bpr)+vpqrs(kpr,bpr,lpr,k))
    
    if (a.eq.bpr) func=func-(vpqrs(kpr,k,lpr,apr)+vpqrs(kpr,apr,lpr,k))

    if (k.eq.kpr) func=func+(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))
    
    if (k.eq.lpr) func=func+(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
     
    func=func/sqrt(2.0d0)
    
    return

  end function C1_ph_2p2h

!#######################################################################

!!$ a'|=b' and k'|=l'; spin case 2
  
  function C2_ph_2p2h(a,k,apr,bpr,kpr,lpr) result(func)
    
    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
    real(dp)            :: func

    func=0.0d0

    if (a.eq.apr) func=func-(vpqrs(kpr,k,lpr,bpr)-vpqrs(kpr,bpr,lpr,k))
    
    if (a.eq.bpr) func=func+(vpqrs(kpr,k,lpr,apr)-vpqrs(kpr,apr,lpr,k))
    
    if (k.eq.kpr) func=func+(vpqrs(a,apr,lpr,bpr)-vpqrs(a,bpr,lpr,apr))
    
    if (k.eq.lpr) func=func-(vpqrs(a,apr,kpr,bpr)-vpqrs(a,bpr,kpr,apr))
    
    func=sqrt(3.0d0)*func/sqrt(2.0d0)
   
  end function C2_ph_2p2h
  
!#######################################################################

!!$ a'=b' and k'|=l' 
  
  function C3_ph_2p2h(j,k,ipr,kpr,lpr) result(func)
    
    integer, intent(in) :: j,k,ipr,kpr,lpr
    real(dp)            :: func

    func=0.0d0
    
    if (j.eq.ipr) func=func-(vpqrs(kpr,k,lpr,ipr)+vpqrs(kpr,ipr,lpr,k))
    
    if (k.eq.kpr) func=func+vpqrs(j,ipr,lpr,ipr)
    
    if (k.eq.lpr) func=func+vpqrs(j,ipr,kpr,ipr)

  end function C3_ph_2p2h

!#######################################################################

!!$ a'|=b' and k'=l'

  function C4_ph_2p2h(a,k,apr,bpr,kpr) result(func)
    
    integer, intent(in) :: a,k,apr,bpr,kpr
    real(dp)            :: func
  
    func=0.0d0

    if (a.eq.apr) func=func+vpqrs(kpr,k,kpr,bpr)

    if (a.eq.bpr) func=func+vpqrs(kpr,k,kpr,apr)

    if (k.eq.kpr) func=func-(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))

  end function C4_ph_2p2h

!#######################################################################

!!$ a'=b' and k'=l'
  
  function C5_ph_2p2h(a,k,apr,kpr) result(func)
    
    integer, intent(in) :: a,k,apr,kpr
    real(dp)            :: func

    func=0.0d0
     
    if (a.eq.apr) func=func+vpqrs(kpr,apr,kpr,k)

    if (k.eq.kpr) func=func-vpqrs(a,apr,kpr,apr)
    
    func=sqrt(2.0d0)*func
   
  end function C5_ph_2p2h

!#######################################################################
  
!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************2PH-2PH BLOCK*********************************
!!$*******************************************************************************
!!$*******************************************************************************  

!!$ 2P2H-2P2H block is organized in subblocks according to the different
!!$ classes of double excitations. 1- doubles k=l,a=b; 2-doubles
!!$  k=l,a|=b; 3-doubles k|=l,a=b; 4I-doubles k|=l,a|=b (spin case 1);
!!$ 4II-doubles i|=j,a|=b (spin case 2).

!#######################################################################

  real(dp) function K_2p2h_2p2h(e_a,e_b,e_k,e_l)
    
    real(dp), intent(in) :: e_a,e_b,e_k,e_l
    
    K_2p2h_2p2h=e_a+e_b-e_k-e_l
  
  end function K_2p2h_2p2h

!#######################################################################

!!$ (1,1)-Caakk,a'a'k'k' (case 25 in Trofimov's Ph.D)
  
  function C_1_1(a,k,apr,kpr) result(func)
    
    integer, intent(in) :: a,k,apr,kpr
    real(dp)            :: func
    
    func=0.0d0
    
    if (k.eq.kpr) func=func+vpqrs(a,apr,a,apr)

    if (a.eq.apr) func=func+vpqrs(kpr,k,kpr,k)

    if (a.eq.apr.and.k.eq.kpr) func=func-2.0d0*(2.0d0*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
  
  end function C_1_1
  
!#######################################################################

!!$ (2,1)-Cabkk,a'a'k'k' (case 20 in Trofimov's Ph.D)
  
  function C_2_1(a,b,k,apr,kpr) result(func)
    
    integer, intent(in) :: a,b,k,apr,kpr
    real(dp)            :: func

    func=0.0d0
    
    if (k.eq.kpr) func=func+vpqrs(a,apr,b,apr)

    if (a.eq.apr.and.b.eq.apr) func=func+0.5_dp*vpqrs(kpr,k,kpr,k)

    if (b.eq.apr.and.k.eq.kpr) func=func-(2._dp*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))

    if (a.eq.apr.and.k.eq.kpr) func=func-(2._dp*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
    
    func=func*sqrt(2._dp)

  end function C_2_1

!#######################################################################

!!$ (3,1)-Caakl,a'a'k'k' (case 15 in Trofimov's Ph.D)
  
  function C_3_1(a,k,l,apr,kpr) result(func)
    
    integer, intent(in) :: a,k,l,apr,kpr
    real(dp)            :: func
    
    func=0.0d0
    
    if (k.eq.kpr.and.l.eq.kpr) func=func-0.5_dp*vpqrs(a,apr,a,apr)
 
    if (a.eq.apr) func=func-vpqrs(kpr,k,kpr,l)
    
    if (a.eq.apr.and.l.eq.kpr) func=func+(2._dp*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))

    if (a.eq.apr.and.k.eq.kpr) func=func+(2._dp*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
    
    func=func*sqrt(2._dp)
   
  end function C_3_1
  
!#######################################################################

!!$ (4i,1)-Cabkl,a'a'k'k' (spin case I) (case 5 in Trofimov's Ph.D)

  real(dp) function C_4i_1(a,b,k,l,apr,kpr)
    
    integer, intent(in) :: a,b,k,l,apr,kpr
        
    C_4i_1=0._dp
    
    if(kdelta(k,kpr)*kdelta(l,kpr) .eq. 1)&
         C_4i_1=C_4i_1-vpqrs(a,apr,b,apr)
    if(kdelta(a,apr)*kdelta(b,apr) .eq. 1)&
         C_4i_1=C_4i_1-vpqrs(kpr,k,kpr,l)
    if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
         C_4i_1=C_4i_1+(2._dp*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_4i_1=C_4i_1+(2._dp*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_4i_1=C_4i_1+(2._dp*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_4i_1=C_4i_1+(2._dp*vpqrs(b,apr,kpr,l)-vpqrs(b,l,kpr,apr))
    
  end function C_4i_1

!#######################################################################

!!$ (4ii,1)-Cabkl,a'a'k'k' (spin case II) (case 10 in Trofimov's Ph.D) 
  
  real(dp) function C_4ii_1(a,b,k,l,apr,kpr)
    
    integer, intent(in) :: a,b,k,l,apr,kpr
    
    C_4ii_1=0._dp
    
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_4ii_1=C_4ii_1+vpqrs(a,l,kpr,apr)
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_4ii_1=C_4ii_1+vpqrs(b,k,kpr,apr)
    if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
         C_4ii_1=C_4ii_1-vpqrs(a,k,kpr,apr)
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_4ii_1=C_4ii_1-vpqrs(b,l,kpr,apr)
    
    C_4ii_1=sqrt(3._dp)*C_4ii_1

  end function C_4ii_1

!#######################################################################

!!$ (2,2)-Cabkk,a'b'k'k' (case 19 in Trofimov's Ph.D)

  real(dp) function C_2_2(a,b,k,apr,bpr,kpr)
    
    integer, intent(in) :: a,b,k,apr,bpr,kpr
    
    C_2_2=0._dp
    
    if(kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2+(vpqrs(a,apr,b,bpr)+vpqrs(a,bpr,b,apr))
    if(kdelta(a,apr)*kdelta(b,bpr) .eq. 1)&
         C_2_2=C_2_2+vpqrs(kpr,k,kpr,k)
    
    if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2-(2._dp*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2-(2._dp*vpqrs(b,bpr,kpr,k)-vpqrs(b,k,kpr,bpr))
    
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2-(2._dp*vpqrs(a,bpr,kpr,k)-vpqrs(a,k,kpr,bpr))
    if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2-(2._dp*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
    
  end function C_2_2
  
!#######################################################################

!!$ (3,2)-Caakl,a'b'k'k' (case 14 in Trofimov's Ph.D)  
  
  real(dp) function C_3_2(a,k,l,apr,bpr,kpr)
    
    integer, intent(in) ::a,k,l,apr,bpr,kpr 
    
    C_3_2=0._dp
    
    if(kdelta(k,kpr)*kdelta(l,kpr) .eq. 1)&
         C_3_2=C_3_2-vpqrs(a,apr,a,bpr)
    if(kdelta(a,apr)*kdelta(a,bpr) .eq. 1)&
         C_3_2=C_3_2-vpqrs(kpr,k,kpr,l)

    if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
         C_3_2=C_3_2+(2._dp*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
         C_3_2=C_3_2+(2._dp*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_3_2=C_3_2+(2._dp*vpqrs(a,bpr,kpr,k)-vpqrs(a,k,kpr,bpr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_3_2=C_3_2+(2._dp*vpqrs(a,bpr,kpr,l)-vpqrs(a,l,kpr,bpr))
    
  end function C_3_2

!#######################################################################

!!$ (4i,2)-Cabkl,a'b'k'k' (spin case I) (case 4 in Trofimov's Ph.D)

  real(dp) function C_4i_2(a,b,k,l,apr,bpr,kpr)
    
    integer, intent(in) ::a,b,k,l,apr,bpr,kpr 
    
    C_4i_2=0._dp  
  
    if(kdelta(k,kpr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2-(vpqrs(a,apr,b,bpr)+vpqrs(a,bpr,b,apr))
    if(kdelta(a,apr)*kdelta(b,bpr) .eq. 1)&
         C_4i_2=C_4i_2-2._dp*vpqrs(kpr,k,kpr,l)
    
    if(kdelta(b,bpr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._dp*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._dp*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._dp*vpqrs(b,bpr,kpr,k)-vpqrs(b,k,kpr,bpr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._dp*vpqrs(b,bpr,kpr,l)-vpqrs(b,l,kpr,bpr))
    
    if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._dp*vpqrs(a,bpr,kpr,k)-vpqrs(a,k,kpr,bpr))
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._dp*vpqrs(a,bpr,kpr,l)-vpqrs(a,l,kpr,bpr))
    if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._dp*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
    if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._dp*vpqrs(b,apr,kpr,l)-vpqrs(b,l,kpr,apr))
    
    
    C_4i_2=C_4i_2/sqrt(2._dp)
    
  end function C_4i_2

!#######################################################################

!!$ (4ii,2)-Cabkl,a'b'k'k' (spin case II) (case 9 in Trofimov's Ph.D)

  real(dp) function C_4ii_2(a,b,k,l,apr,bpr,kpr)
    
    integer, intent(in) ::a,b,k,l,apr,bpr,kpr 
    
    C_4ii_2=0._dp  
    
    if(kdelta(b,bpr)*kdelta(l,kpr) .eq. 1)&
         C_4ii_2=C_4ii_2-vpqrs(a,k,kpr,apr)
    if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
         C_4ii_2=C_4ii_2+vpqrs(a,l,kpr,apr)
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_4ii_2=C_4ii_2+vpqrs(b,k,kpr,bpr)
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_4ii_2=C_4ii_2-vpqrs(b,l,kpr,bpr)
    
    if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
         C_4ii_2=C_4ii_2-vpqrs(a,k,kpr,bpr)
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_4ii_2=C_4ii_2+vpqrs(a,l,kpr,bpr)
    if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
         C_4ii_2=C_4ii_2+vpqrs(b,k,kpr,apr)
    if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
         C_4ii_2=C_4ii_2-vpqrs(b,l,kpr,apr)

    
    C_4ii_2=sqrt(3._dp)*C_4ii_2/sqrt(2._dp)
    
  end function C_4ii_2
    
!#######################################################################

!!$ (3,3)-Caakl,a'a'k'l' (case 13  in Trofimov's Ph.D)  
    
    real(dp) function C_3_3(a,k,l,apr,kpr,lpr)
    
      integer, intent(in) ::a,k,l,apr,kpr,lpr 
    
      C_3_3=0._dp
      
      if(kdelta(k,kpr)*kdelta(l,lpr) .eq. 1)&
           C_3_3=C_3_3+vpqrs(a,apr,a,apr)
      if(kdelta(a,apr) .eq. 1)&
           C_3_3=C_3_3+(vpqrs(k,kpr,lpr,l)+vpqrs(kpr,l,lpr,k))
      
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_3_3=C_3_3-(2._dp*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_3_3=C_3_3-(2._dp*vpqrs(a,apr,lpr,l)-vpqrs(a,l,lpr,apr))
      
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_3_3=C_3_3-(2._dp*vpqrs(a,apr,lpr,k)-vpqrs(a,k,lpr,apr))
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_3_3=C_3_3-(2._dp*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
      
    end function C_3_3

!#######################################################################

!!$ (4i,3)-Cabkl,a'a'k'l' (spin case I) (case 3 in Trofimov's Ph.D)

    real(dp) function C_4i_3(a,b,k,l,apr,kpr,lpr)
    
      integer, intent(in) ::a,b,k,l,apr,kpr,lpr 
      
      C_4i_3=0._dp 
    
      if(kdelta(k,kpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_3=C_4i_3+2._dp*vpqrs(a,apr,b,apr)
      if(kdelta(a,apr)*kdelta(b,apr) .eq. 1)&
           C_4i_3=C_4i_3+(vpqrs(kpr,k,lpr,l)+vpqrs(kpr,l,lpr,k))
      
      if(kdelta(b,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._dp*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
      if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._dp*vpqrs(a,apr,lpr,l)-vpqrs(a,l,lpr,apr))
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._dp*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._dp*vpqrs(b,apr,lpr,l)-vpqrs(b,l,lpr,apr))
      
      if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._dp*vpqrs(a,apr,lpr,k)-vpqrs(a,k,lpr,apr))
      if(kdelta(b,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._dp*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._dp*vpqrs(b,apr,lpr,k)-vpqrs(b,k,lpr,apr))
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._dp*vpqrs(b,apr,kpr,l)-vpqrs(b,l,kpr,apr))
      
      C_4i_3=C_4i_3/sqrt(2._dp)
      
    end function C_4i_3

!#######################################################################

!!$ (4ii,3)-Cabkl,a'a'k'l' (spin case II) (case 8 in Trofimov's Ph.D)    

    real(dp) function C_4ii_3(a,b,k,l,apr,kpr,lpr)
    
      integer, intent(in) ::a,b,k,l,apr,kpr,lpr 
    
      C_4ii_3=0._dp
      
      if(kdelta(b,apr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_3=C_4ii_3+vpqrs(a,k,kpr,apr)
      if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_3=C_4ii_3-vpqrs(a,l,lpr,apr)
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_3=C_4ii_3-vpqrs(b,k,kpr,apr)
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_3=C_4ii_3+vpqrs(b,l,lpr,apr)
      
      if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_3=C_4ii_3+vpqrs(a,k,lpr,apr)
      if(kdelta(b,apr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_3=C_4ii_3-vpqrs(a,l,kpr,apr)
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_3=C_4ii_3-vpqrs(b,k,lpr,apr)
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_3=C_4ii_3+vpqrs(b,l,kpr,apr) 
      
      
      C_4ii_3=sqrt(3._dp)*C_4ii_3/sqrt(2._dp)
      
    end function C_4ii_3

!#######################################################################

!!$ (4ii,4ii)-Cabkl,a'b'k'l' (case 7 in Trofimov's Ph.D)

    real(dp) function C_4ii_4ii(a,b,k,l,apr,bpr,kpr,lpr)
    
      integer, intent(in) ::a,b,k,l,apr,bpr,kpr,lpr 
    
      C_4ii_4ii=0._dp    

      if(kdelta(k,kpr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,apr,b,bpr)-vpqrs(a,bpr,b,apr))
      if(kdelta(a,apr)*kdelta(b,bpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(kpr,k,lpr,l)-vpqrs(kpr,l,lpr,k))
      
      if(kdelta(b,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(a,apr,kpr,k)-1.5_dp*vpqrs(a,k,kpr,apr))
      if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(a,apr,lpr,l)-1.5_dp*vpqrs(a,l,lpr,apr))
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(b,bpr,kpr,k)-1.5_dp*vpqrs(b,k,kpr,bpr))
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(b,bpr,lpr,l)-1.5_dp*vpqrs(b,l,lpr,bpr))
      
      if(kdelta(b,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,apr,lpr,k)-1.5_dp*vpqrs(a,k,lpr,apr))
      if(kdelta(b,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,apr,kpr,l)-1.5_dp*vpqrs(a,l,kpr,apr))
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(b,bpr,lpr,k)-1.5_dp*vpqrs(b,k,lpr,bpr))
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(b,bpr,kpr,l)-1.5_dp*vpqrs(b,l,kpr,bpr))      

      if(kdelta(b,apr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,bpr,kpr,k)-1.5_dp*vpqrs(a,k,kpr,bpr))
      if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,bpr,lpr,l)-1.5_dp*vpqrs(a,l,lpr,bpr))
      if(kdelta(a,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(b,apr,kpr,k)-1.5_dp*vpqrs(b,k,kpr,apr))
      if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(b,apr,lpr,l)-1.5_dp*vpqrs(b,l,lpr,apr))


      if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(a,bpr,lpr,k)-1.5_dp*vpqrs(a,k,lpr,bpr))
      if(kdelta(b,apr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(a,bpr,kpr,l)-1.5_dp*vpqrs(a,l,kpr,bpr))
      if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(b,apr,lpr,k)-1.5_dp*vpqrs(b,k,lpr,apr))
      if(kdelta(a,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(b,apr,kpr,l)-1.5_dp*vpqrs(b,l,kpr,apr))

    end function C_4ii_4ii

!#######################################################################

!!$ (4i,4ii)-Cabkl,a'b'k'l' (case 2 in Trofimov's Ph.D)
    
    real(dp) function C_4i_4ii(a,b,k,l,apr,bpr,kpr,lpr)
      
      integer, intent(in) ::a,b,k,l,apr,bpr,kpr,lpr 
      
      C_4i_4ii=0._dp

      if(kdelta(b,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii+vpqrs(a,k,kpr,apr)
      if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii-vpqrs(a,l,lpr,apr)
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii-vpqrs(b,k,kpr,bpr)
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii+vpqrs(b,l,lpr,bpr)
      
      
      if(kdelta(b,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii-vpqrs(a,k,lpr,apr)
      if(kdelta(b,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii+vpqrs(a,l,kpr,apr)
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii+vpqrs(b,k,lpr,bpr)
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii-vpqrs(b,l,kpr,bpr)

      
      if(kdelta(b,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii-vpqrs(a,k,kpr,bpr)
      if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii+vpqrs(a,l,lpr,bpr)
      if(kdelta(a,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii+vpqrs(b,k,kpr,apr)
      if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii-vpqrs(b,l,lpr,apr)

      if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii+vpqrs(a,k,lpr,bpr)
      if(kdelta(b,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii-vpqrs(a,l,kpr,bpr)
      if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii-vpqrs(b,k,lpr,apr)
      if(kdelta(a,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4ii=C_4i_4ii+vpqrs(b,l,kpr,apr) 
      
      
      C_4i_4ii=sqrt(3._dp)*C_4i_4ii/2._dp
      
    end function C_4i_4ii

!#######################################################################

!!$ (4ii,4i)-Cabkl,a'b'k'l' (case 6 in Trofimov's Ph.D)
    
    real(dp) function C_4ii_4i(a,b,k,l,apr,bpr,kpr,lpr)
      
      integer, intent(in) ::a,b,k,l,apr,bpr,kpr,lpr 
      
      C_4ii_4i=0._dp

      if(kdelta(b,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i+vpqrs(a,k,kpr,apr)
      if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i-vpqrs(a,l,lpr,apr)
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i-vpqrs(b,k,kpr,bpr)
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i+vpqrs(b,l,lpr,bpr)
      
      
      if(kdelta(b,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i+vpqrs(a,k,lpr,apr)
      if(kdelta(b,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i-vpqrs(a,l,kpr,apr)
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i-vpqrs(b,k,lpr,bpr)
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i+vpqrs(b,l,kpr,bpr)

      
      if(kdelta(b,apr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i+vpqrs(a,k,kpr,bpr)
      if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i-vpqrs(a,l,lpr,bpr)
      if(kdelta(a,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i-vpqrs(b,k,kpr,apr)
      if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i+vpqrs(b,l,lpr,apr)

      if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i+vpqrs(a,k,lpr,bpr)
      if(kdelta(b,apr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i-vpqrs(a,l,kpr,bpr)
      if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i-vpqrs(b,k,lpr,apr)
      if(kdelta(a,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4i=C_4ii_4i+vpqrs(b,l,kpr,apr) 
      
      
      C_4ii_4i=sqrt(3._dp)*C_4ii_4i/2._dp
      
    end function C_4ii_4i

!#######################################################################

!!$ (4i,4i)-Cabkl,a'b'k'l' (case 1 in Trofimov's Ph.D)

    real(dp) function C_4i_4i(a,b,k,l,apr,bpr,kpr,lpr)
      
      integer, intent(in) ::a,b,k,l,apr,bpr,kpr,lpr 
      
      C_4i_4i=0._dp    

      if(kdelta(k,kpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i+2._dp*(vpqrs(a,apr,b,bpr)+vpqrs(a,bpr,b,apr))
      if(kdelta(a,apr)*kdelta(b,bpr) .eq. 1)&
           C_4i_4i=C_4i_4i+2._dp*(vpqrs(kpr,k,lpr,l)+vpqrs(kpr,l,lpr,k))

      if(kdelta(b,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
      if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(a,apr,lpr,l)-vpqrs(a,l,lpr,apr))
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(b,bpr,kpr,k)-vpqrs(b,k,kpr,bpr))
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(b,bpr,lpr,l)-vpqrs(b,l,lpr,bpr))
      
      if(kdelta(b,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(a,apr,lpr,k)-vpqrs(a,k,lpr,apr))
      if(kdelta(b,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(b,bpr,lpr,k)-vpqrs(b,k,lpr,bpr))
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(b,bpr,kpr,l)-vpqrs(b,l,kpr,bpr))
      
      if(kdelta(b,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(a,bpr,kpr,k)-vpqrs(a,k,kpr,bpr))
      if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(a,bpr,lpr,l)-vpqrs(a,l,lpr,bpr))
      if(kdelta(a,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
      if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(b,apr,lpr,l)-vpqrs(b,l,lpr,apr)) 

      if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(a,bpr,lpr,k)-vpqrs(a,k,lpr,bpr))
      if(kdelta(b,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(a,bpr,kpr,l)-vpqrs(a,l,kpr,bpr))
      if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(b,apr,lpr,k)-vpqrs(b,k,lpr,apr))
      if(kdelta(a,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._dp*vpqrs(b,apr,kpr,l)-vpqrs(b,l,kpr,apr))
      
      C_4i_4i=0.5_dp*C_4i_4i

    end function C_4i_4i

!#######################################################################

  end module adc_ph
