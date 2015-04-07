module adc_ph

!!$The module contains an assortment of functions needed to calculate
!!$matrix elements of the ADC matrix of the polarization propagator.
!!$For the spin-orbital expressions see A.B. Trofimov et al, JCP 111,9982 (1999).
!!$Spin free expressions were taken from the Ph.D. thesis of A.B. Trofimov.
  
  use constants
  use parameters
  use misc
  
  implicit none
  
contains

!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************PH-PH BLOCK***********************************
!!$*******************************************************************************
!!$*******************************************************************************


  subroutine MP2(E_MP2)

    real(d), intent(out) :: E_MP2
    
    integer :: c,dd,i,j, nsym1,nsym2,nsym3,u,v,r,s,cnt
    real(d) :: DA,eijc,term
    real(d) :: vpqrs

    external vpqrs


    E_MP2 = 0._d

    cnt=0

    do c=nOcc+1,nBas
       r=roccnum(c)

       do dd=nOcc+1,nBas
          s=roccnum(dd)

          do i=1,nOcc
             u=roccnum(i)

             do j=1,nOcc
                v=roccnum(j)

             nsym1=MT(orbSym(u),orbSym(v))
             nsym2=MT(orbSym(s),orbSym(r))
             
             
             if  (MT(nsym1,nsym2) .eq. 1)  then
                
                cnt=cnt+1
                term=0._d

                eijc=e(u)+e(v)-e(r)-e(s)

                term=term + vpqrs(r,u,s,v)*(2._d*vpqrs(r,u,s,v)-vpqrs(r,v,s,u))

                term=term/eijc
                
                E_MP2 = E_MP2 + term
                
             end if

          end do
       end do
    end do 
  end do
  
  end subroutine MP2


!!$Indices are supplied in the order: PH,PH

!!$Zeroth order contribution K_ak,a'k'. The condition that a=a', k=k'
!!$is checked by the calling procedure

  REAL(D) FUNCTION K_PH_PH(E_A,E_K)
    
    REAL(D), INTENT(IN) :: E_A,E_K
    
    K_PH_PH=E_A-E_K
    
  end function K_ph_ph

!!$First order contribution C_ak,a'k'
  
  real(d) function C1_ph_ph(a,k,a1,k1)
    
    integer, intent(in) :: a,k,a1,k1
    real(d) :: vpqrs
      
    external vpqrs

    C1_ph_ph=2._d*vpqrs(a,k,a1,k1)-vpqrs(a,a1,k,k1)

  end function C1_ph_ph
    
!!$Second order contribution CA_ak,a'k'. The condition that k=k' is
!!$checked in the calling procedure.

  real(d) function CA_ph_ph(a,a1)

    integer, intent(in) :: a,a1
    
    integer :: c,i,j, nsym1,nsym2,nsym3,u,v,r,cnt
    real(d) :: DA,eijc,term
    real(d) :: vpqrs
    
    external vpqrs


!!$ Faster version (taking symmetries into account)

!!$ i|=j 

!!$    CA_ph_ph=0._d
!!$   
!!$    do c=nOcc+1,nBas
!!$       do i=1,nOcc
!!$          do j=i+1,nOcc
!!$             
!!$             term=0._d
!!$
!!$             eijc=e(roccnum(i))+e(roccnum(j))-e(roccnum(c))
!!$             DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))
!!$
!!$             term=term+vpqrs(a,i,c,j)*(2._d*vpqrs(i,a1,j,c)-vpqrs(i,c,j,a1))
!!$             term=term+vpqrs(a,j,c,i)*(2._d*vpqrs(i,c,j,a1)-vpqrs(i,a1,j,c))
!!$             term=DA*term
!!$             
!!$             CA_ph_ph=CA_ph_ph+term
!!$             
!!$          end do
!!$       end do
!!$    end do
!!$    
!!$    CA_ph_ph=2._d*CA_ph_ph
    
!!$ i=j
    
!!$    do c=nOcc+1,nBas
!!$       do i=1,nOcc
!!$          
!!$          term=0._d
!!$          
!!$          eijc=2._d*e(roccnum(i))-e(roccnum(c))
!!$          DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))
!!$          
!!$          term=term+2.0_d*vpqrs(a,i,c,i)*vpqrs(a1,i,c,i)
!!$          term=DA*term
!!$          
!!$          CA_ph_ph=CA_ph_ph+term
!!$        
!!$       end do
!!$    end do
    
    
!!$    CA_ph_ph=-0.5_d*CA_ph_ph

!!$ Dumbed down version


    CA_ph_ph=0._d
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
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
                term=0._d

                eijc=e(u)+e(v)-e(r)
                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!!$                term=term+vpqrs(a,i,c,j)*(vpqrs(i,a1,j,c)-2._d*vpqrs(i,c,j,a1))
!!$                term=term+vpqrs(a,j,c,i)*(vpqrs(i,c,j,a1)-2._d*vpqrs(i,a1,j,c))
                term=DA*term
                
                CA_ph_ph=CA_ph_ph+term
                
             end if
          end do
       end do
    end do 
    
    CA_ph_ph=-0.5_d*CA_ph_ph            
    
  end function CA_ph_ph
    
             
!!$Second order contribution CB_ak,a'k'. The condition that a=a' is
!!$checked in the calling procedure.

  real(d) function CB_ph_ph(k,k1)

    integer, intent(in) :: k,k1
    
    integer :: c,dd,i, nsym1, nsym2, nsym3,u,r,s
    real(d) :: DB,eicd,term
    real(d) :: vpqrs
    
    external vpqrs

!!$    CB_ph_ph=0._d             
  
!!$! Faster version (taking symmetries into account)

!!$! b|=c
!!$    
!!$    do c=nOcc+1,nBas
!!$       do dd=c+1,nBas
!!$          do i=1,nOcc
!!$            
!!$             term=0._d
!!$ 
!!$             eicd=e(roccnum(i))-e(roccnum(c))-e(roccnum(dd))
!!$             DB=(eicd+0.5_d*(e(k)+e(k1)))/((eicd+e(k))*(eicd+e(k1)))
!!$             
!!$             term=term+vpqrs(c,k,dd,i)*(2._d*vpqrs(k1,c,i,dd)-vpqrs(k1,dd,i,c))
!!$             term=term+vpqrs(c,i,dd,k)*(2._d*vpqrs(k1,dd,i,c)-vpqrs(k1,c,i,dd))
!!$             term=DB*term
!!$             
!!$             CB_ph_ph=CB_ph_ph+term
!!$             
!!$          end do
!!$       end do
!!$    end do
!!$    
!!$    CB_ph_ph=2._d*CB_ph_ph
       
!!$! b=c
    
!!$    do c=nOcc+1,nBas
!!$       do i=1,nOcc
!!$          
!!$          term=0._d
!!$          
!!$          eicd=e(roccnum(i))-2._d*e(roccnum(c))
!!$          DB=(eicd+0.5_d*(e(k)+e(k1)))/((eicd+e(k))*(eicd+e(k1)))
!!$
!!$          term=term+2._d*vpqrs(c,k,c,i)*vpqrs(k1,c,i,c)
!!$          term=DB*term
!!$          
!!$          CB_ph_ph=CB_ph_ph+term
!!$          
!!$       end do
!!$    end do
!!$    
!!$    CB_ph_ph=-0.5_d*CB_ph_ph

!!$ Dumbed down version

    CB_ph_ph=0._d   
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
            
                term=0._d
                
                eicd=e(roccnum(i))-e(roccnum(c))-e(roccnum(dd))
                DB=(eicd+0.5_d*(e(k)+e(k1)))/((eicd+e(k))*(eicd+e(k1)))
                
                term=term+vpqrs(r,k,s,u)*(2._d*vpqrs(k1,r,u,s)-vpqrs(k1,s,u,r))
                term=term+vpqrs(r,u,s,k)*(2._d*vpqrs(k1,s,u,r)-vpqrs(k1,r,u,s))
!!$
!!$                term=term+vpqrs(c,k,dd,i)*(vpqrs(k1,c,i,dd)-2._d*vpqrs(k1,dd,i,c))
!!$                term=term+vpqrs(c,i,dd,k)*(vpqrs(k1,dd,i,c)-2._d*vpqrs(k1,c,i,dd))
                
                term=DB*term
                
                CB_ph_ph=CB_ph_ph+term
                
             end if
             
          end do
       end do
    end do
    
    CB_ph_ph=-0.5_d*CB_ph_ph   
   
  end function CB_ph_ph
  
!!$Second order contribution CC_ak,a'k'.

  real(d) function CC_ph_ph(a,k,a1,k1)

    integer, intent(in) :: a,a1,k,k1
    
    integer :: c,i, nsym1, nsym2, nsym3,r,s
    real(d) :: DC,eic,term
    real(d) :: vpqrs
    
    real(d), dimension(nocc) :: tau

    external vpqrs
    
    CC_ph_ph=0._d

    do c=nOcc+1,nBas
       r=roccnum(c)
       do i=1,nOcc
          s=roccnum(i)

    
          nsym1=MT(orbSym(r),orbSym(s))
          nsym2=MT(orbSym(k),orbSym(a))
          nsym3=MT(orbSym(k1),orbSym(a1))
          
             
          if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
             
             term=0._d
             
             eic=e(s)-e(r)
             DC=(0.5_d*(e(k)+e(k1)-e(a)-e(a1))+eic)/((e(k)-e(a)+eic)*(e(k1)-e(a1)+eic))
             
             term=(2._d*vpqrs(a,k,r,s)-vpqrs(a,s,r,k))*(2._d*vpqrs(a1,k1,r,s)-vpqrs(a1,s,r,k1))
             term=term*DC

             CC_ph_ph=CC_ph_ph+term
             
          end if
          
       end do
    end do

  end function CC_ph_ph


!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************PH-2PH BLOCK**********************************
!!$*******************************************************************************
!!$*******************************************************************************

!!$  We distinguish here between five different types of coupling. 
!!$ Calculating Cak,a'b'k'l'
  
!!$ a'|=b' and k'|=l'; spin case 1
  
  real(d) function C1_ph_2p2h(a,k,apr,bpr,kpr,lpr)
    
    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
    
    real(d) :: vpqrs
    
    external vpqrs

    C1_ph_2p2h=0._d
    
    if(kdelta(a,apr) .eq. 1)& 
         C1_ph_2p2h=C1_ph_2p2h-(vpqrs(kpr,k,lpr,bpr)+vpqrs(kpr,bpr,lpr,k))
    if(kdelta(a,bpr) .eq. 1)&
         C1_ph_2p2h=C1_ph_2p2h-(vpqrs(kpr,k,lpr,apr)+vpqrs(kpr,apr,lpr,k))
    if(kdelta(k,kpr) .eq. 1)&
         C1_ph_2p2h=C1_ph_2p2h+(vpqrs(a,apr,lpr,bpr)+vpqrs(a,bpr,lpr,apr))
    if(kdelta(k,lpr) .eq. 1)&
         C1_ph_2p2h=C1_ph_2p2h+(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
     
    C1_ph_2p2h=C1_ph_2p2h/sqrt(2._d)
    
    
  end function C1_ph_2p2h

!!$ a'|=b' and k'|=l'; spin case 2
  
  real(d) function C2_ph_2p2h(a,k,apr,bpr,kpr,lpr)
    
    integer, intent(in) :: a,k,apr,bpr,kpr,lpr
    
    real(d) :: vpqrs
    
    external vpqrs

    C2_ph_2p2h=0._d 
    
    if(kdelta(a,apr) .eq. 1)&
         C2_ph_2p2h=C2_ph_2p2h-(vpqrs(kpr,k,lpr,bpr)-vpqrs(kpr,bpr,lpr,k))
    if(kdelta(a,bpr) .eq. 1)&
         C2_ph_2p2h=C2_ph_2p2h+(vpqrs(kpr,k,lpr,apr)-vpqrs(kpr,apr,lpr,k))
    if(kdelta(k,kpr) .eq. 1)&
         C2_ph_2p2h=C2_ph_2p2h+(vpqrs(a,apr,lpr,bpr)-vpqrs(a,bpr,lpr,apr))
    if(kdelta(k,lpr) .eq. 1)&
         C2_ph_2p2h=C2_ph_2p2h-(vpqrs(a,apr,kpr,bpr)-vpqrs(a,bpr,kpr,apr))

    C2_ph_2p2h=sqrt(3._d)*C2_ph_2p2h/sqrt(2._d)
    
  end function C2_ph_2p2h
  
!!$ a'=b' and k'|=l' 
  
  real(d) function C3_ph_2p2h(j,k,ipr,kpr,lpr)
    
    integer, intent(in) :: j,k,ipr,kpr,lpr
    
    real(d) :: vpqrs
    
    external vpqrs

    C3_ph_2p2h=0._d
    
    if(kdelta(j,ipr) .eq. 1)&
         C3_ph_2p2h=C3_ph_2p2h-(vpqrs(kpr,k,lpr,ipr)+vpqrs(kpr,ipr,lpr,k))
    if(kdelta(k,kpr) .eq. 1)&
         C3_ph_2p2h=C3_ph_2p2h+vpqrs(j,ipr,lpr,ipr)
    if(kdelta(k,lpr) .eq. 1)&
         C3_ph_2p2h=C3_ph_2p2h+vpqrs(j,ipr,kpr,ipr)

 !   C3_ph_2p2h=-C3_ph_2p2h
    
  end function C3_ph_2p2h

!!$ a'|=b' and k'=l'

  real(d) function C4_ph_2p2h(a,k,apr,bpr,kpr)
    
    integer, intent(in) :: a,k,apr,bpr,kpr
    
    real(d) :: vpqrs
    
    external vpqrs

    C4_ph_2p2h=0._d

    if(kdelta(a,apr) .eq. 1)&
         C4_ph_2p2h=C4_ph_2p2h+vpqrs(kpr,k,kpr,bpr)
    if(kdelta(a,bpr) .eq. 1)&
         C4_ph_2p2h=C4_ph_2p2h+vpqrs(kpr,k,kpr,apr)
    if(kdelta(k,kpr) .eq. 1)&
         C4_ph_2p2h=C4_ph_2p2h-(vpqrs(a,apr,kpr,bpr)+vpqrs(a,bpr,kpr,apr))
    
 !   C4_ph_2p2h=-C4_ph_2p2h

  end function C4_ph_2p2h

!!$ a'=b' and k'=l'
  
  real(d) function C5_ph_2p2h(a,k,apr,kpr)
    
    integer, intent(in) :: a,k,apr,kpr
    
    real(d) :: vpqrs
    
    external vpqrs

    C5_ph_2p2h=0._d
     
    if(kdelta(a,apr) .eq. 1)&
         C5_ph_2p2h=C5_ph_2p2h+vpqrs(kpr,apr,kpr,k)
    if(kdelta(k,kpr) .eq. 1)&
         C5_ph_2p2h=C5_ph_2p2h-vpqrs(a,apr,kpr,apr)
    
    C5_ph_2p2h=sqrt(2._d)*C5_ph_2p2h
    
  end function C5_ph_2p2h
  
!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************2PH-2PH BLOCK*********************************
!!$*******************************************************************************
!!$*******************************************************************************  

!!$ 2P2H-2P2H block is organized in subblocks according to the different
!!$ classes of double excitations. 1- doubles k=l,a=b; 2-doubles
!!$  k=l,a|=b; 3-doubles k|=l,a=b; 4I-doubles k|=l,a|=b (spin case 1);
!!$ 4II-doubles i|=j,a|=b (spin case 2).



  real(d) function K_2p2h_2p2h(e_a,e_b,e_k,e_l)
    
    real(d), intent(in) :: e_a,e_b,e_k,e_l
    
    K_2p2h_2p2h=e_a+e_b-e_k-e_l
  
  end function K_2p2h_2p2h
    
!!$ (1,1)-Caakk,a'a'k'k' (case 25 in Trofimov's Ph.D)
  
  real(d) function C_1_1(a,k,apr,kpr)
    
    integer, intent(in) :: a,k,apr,kpr
    
    real(d) :: vpqrs
    
    external vpqrs
    
    C_1_1=0._d
    
    if(kdelta(k,kpr) .eq. 1)&
         C_1_1=C_1_1+vpqrs(a,apr,a,apr)
    if(kdelta(a,apr) .eq. 1)&
         C_1_1=C_1_1+vpqrs(kpr,k,kpr,k)
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_1_1=C_1_1-2._d*(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))

  
  end function C_1_1
  
!!$ (2,1)-Cabkk,a'a'k'k' (case 20 in Trofimov's Ph.D)
  
  real(d) function C_2_1(a,b,k,apr,kpr)
    
    integer, intent(in) :: a,b,k,apr,kpr
    
    real(d) :: vpqrs
    
    external vpqrs
    
    C_2_1=0._d
    
    if(kdelta(k,kpr) .eq. 1)&
         C_2_1=C_2_1+vpqrs(a,apr,b,apr)
    if(kdelta(a,apr)*kdelta(b,apr) .eq. 1)&
         C_2_1=C_2_1+0.5_d*vpqrs(kpr,k,kpr,k)
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_2_1=C_2_1-(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq.1)&
         C_2_1=C_2_1-(2._d*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
    
    C_2_1=C_2_1*sqrt(2._d)

  end function C_2_1
    
!!$ (3,1)-Caakl,a'a'k'k' (case 15 in Trofimov's Ph.D)
  
  real(d) function C_3_1(a,k,l,apr,kpr)
    
    integer, intent(in) :: a,k,l,apr,kpr
    
    real(d) :: vpqrs
    
    external vpqrs
    
    C_3_1=0._d
    
    if(kdelta(k,kpr)*kdelta(l,kpr) .eq. 1)&
         C_3_1=C_3_1-0.5_d*vpqrs(a,apr,a,apr)
    if(kdelta(a,apr) .eq. 1)&
         C_3_1=C_3_1-vpqrs(kpr,k,kpr,l)
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_3_1=C_3_1+(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_3_1=C_3_1+(2._d*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
    
    C_3_1=C_3_1*sqrt(2._d)
    
  end function C_3_1
  


!!$ (4i,1)-Cabkl,a'a'k'k' (spin case I) (case 5 in Trofimov's Ph.D)

  real(d) function C_4i_1(a,b,k,l,apr,kpr)
    
    integer, intent(in) :: a,b,k,l,apr,kpr
    
    real(d) :: vpqrs
    
    external vpqrs
    
    C_4i_1=0._d
    
    if(kdelta(k,kpr)*kdelta(l,kpr) .eq. 1)&
         C_4i_1=C_4i_1-vpqrs(a,apr,b,apr)
    if(kdelta(a,apr)*kdelta(b,apr) .eq. 1)&
         C_4i_1=C_4i_1-vpqrs(kpr,k,kpr,l)
    if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
         C_4i_1=C_4i_1+(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_4i_1=C_4i_1+(2._d*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_4i_1=C_4i_1+(2._d*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_4i_1=C_4i_1+(2._d*vpqrs(b,apr,kpr,l)-vpqrs(b,l,kpr,apr))
    
  end function C_4i_1
  
!!$ (4ii,1)-Cabkl,a'a'k'k' (spin case II) (case 10 in Trofimov's Ph.D) 
  
  real(d) function C_4ii_1(a,b,k,l,apr,kpr)
    
    integer, intent(in) :: a,b,k,l,apr,kpr
    
    real(d) :: vpqrs
    
    external vpqrs
    
    C_4ii_1=0._d
    
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_4ii_1=C_4ii_1+vpqrs(a,l,kpr,apr)
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_4ii_1=C_4ii_1+vpqrs(b,k,kpr,apr)
    if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
         C_4ii_1=C_4ii_1-vpqrs(a,k,kpr,apr)
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_4ii_1=C_4ii_1-vpqrs(b,l,kpr,apr)
    
    C_4ii_1=sqrt(3._d)*C_4ii_1

  end function C_4ii_1

!!$ (2,2)-Cabkk,a'b'k'k' (case 19 in Trofimov's Ph.D)

  real(d) function C_2_2(a,b,k,apr,bpr,kpr)
    
    integer, intent(in) :: a,b,k,apr,bpr,kpr
    
    real(d) :: vpqrs
    
    external vpqrs

    C_2_2=0._d
    
    if(kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2+(vpqrs(a,apr,b,bpr)+vpqrs(a,bpr,b,apr))
    if(kdelta(a,apr)*kdelta(b,bpr) .eq. 1)&
         C_2_2=C_2_2+vpqrs(kpr,k,kpr,k)
    
    if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2-(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2-(2._d*vpqrs(b,bpr,kpr,k)-vpqrs(b,k,kpr,bpr))
    
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2-(2._d*vpqrs(a,bpr,kpr,k)-vpqrs(a,k,kpr,bpr))
    if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
         C_2_2=C_2_2-(2._d*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
    
  end function C_2_2
  
!!$ (3,2)-Caakl,a'b'k'k' (case 14 in Trofimov's Ph.D)  
  
  real(d) function C_3_2(a,k,l,apr,bpr,kpr)
    
    integer, intent(in) ::a,k,l,apr,bpr,kpr 
    
    real(d) :: vpqrs
    
    external vpqrs

    C_3_2=0._d
    
    if(kdelta(k,kpr)*kdelta(l,kpr) .eq. 1)&
         C_3_2=C_3_2-vpqrs(a,apr,a,bpr)
    if(kdelta(a,apr)*kdelta(a,bpr) .eq. 1)&
         C_3_2=C_3_2-vpqrs(kpr,k,kpr,l)

    if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
         C_3_2=C_3_2+(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
         C_3_2=C_3_2+(2._d*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_3_2=C_3_2+(2._d*vpqrs(a,bpr,kpr,k)-vpqrs(a,k,kpr,bpr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_3_2=C_3_2+(2._d*vpqrs(a,bpr,kpr,l)-vpqrs(a,l,kpr,bpr))
    
  end function C_3_2


!!$ (4i,2)-Cabkl,a'b'k'k' (spin case I) (case 4 in Trofimov's Ph.D)

  real(d) function C_4i_2(a,b,k,l,apr,bpr,kpr)
    
    integer, intent(in) ::a,b,k,l,apr,bpr,kpr 
    
    real(d) :: vpqrs
    
    external vpqrs

    C_4i_2=0._d  
  
    if(kdelta(k,kpr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2-(vpqrs(a,apr,b,bpr)+vpqrs(a,bpr,b,apr))
    if(kdelta(a,apr)*kdelta(b,bpr) .eq. 1)&
         C_4i_2=C_4i_2-2._d*vpqrs(kpr,k,kpr,l)
    
    if(kdelta(b,bpr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
    if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._d*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
    if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._d*vpqrs(b,bpr,kpr,k)-vpqrs(b,k,kpr,bpr))
    if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._d*vpqrs(b,bpr,kpr,l)-vpqrs(b,l,kpr,bpr))
    
    if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._d*vpqrs(a,bpr,kpr,k)-vpqrs(a,k,kpr,bpr))
    if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._d*vpqrs(a,bpr,kpr,l)-vpqrs(a,l,kpr,bpr))
    if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._d*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
    if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
         C_4i_2=C_4i_2+(2._d*vpqrs(b,apr,kpr,l)-vpqrs(b,l,kpr,apr))
    
    
    C_4i_2=C_4i_2/sqrt(2._d)
    
  end function C_4i_2

!!$ (4ii,2)-Cabkl,a'b'k'k' (spin case II) (case 9 in Trofimov's Ph.D)

  real(d) function C_4ii_2(a,b,k,l,apr,bpr,kpr)
    
    integer, intent(in) ::a,b,k,l,apr,bpr,kpr 
    
    real(d) :: vpqrs
    
    external vpqrs

    C_4ii_2=0._d  
    
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

    
    C_4ii_2=sqrt(3._d)*C_4ii_2/sqrt(2._d)
    
  end function C_4ii_2
    
!!$ (3,3)-Caakl,a'a'k'l' (case 13  in Trofimov's Ph.D)  
    
    real(d) function C_3_3(a,k,l,apr,kpr,lpr)
    
      integer, intent(in) ::a,k,l,apr,kpr,lpr 
    
      real(d) :: vpqrs
    
      external vpqrs

      C_3_3=0._d
      
      if(kdelta(k,kpr)*kdelta(l,lpr) .eq. 1)&
           C_3_3=C_3_3+vpqrs(a,apr,a,apr)
      if(kdelta(a,apr) .eq. 1)&
           C_3_3=C_3_3+(vpqrs(k,kpr,lpr,l)+vpqrs(kpr,l,lpr,k))
      
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_3_3=C_3_3-(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_3_3=C_3_3-(2._d*vpqrs(a,apr,lpr,l)-vpqrs(a,l,lpr,apr))
      
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_3_3=C_3_3-(2._d*vpqrs(a,apr,lpr,k)-vpqrs(a,k,lpr,apr))
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_3_3=C_3_3-(2._d*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
      
    end function C_3_3

!!$ (4i,3)-Cabkl,a'a'k'l' (spin case I) (case 3 in Trofimov's Ph.D)

    real(d) function C_4i_3(a,b,k,l,apr,kpr,lpr)
    
      integer, intent(in) ::a,b,k,l,apr,kpr,lpr 
      
      real(d) :: vpqrs
    
      external vpqrs

      C_4i_3=0._d 
    
      if(kdelta(k,kpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_3=C_4i_3+2._d*vpqrs(a,apr,b,apr)
      if(kdelta(a,apr)*kdelta(b,apr) .eq. 1)&
           C_4i_3=C_4i_3+(vpqrs(kpr,k,lpr,l)+vpqrs(kpr,l,lpr,k))
      
      if(kdelta(b,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
      if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._d*vpqrs(a,apr,lpr,l)-vpqrs(a,l,lpr,apr))
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._d*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._d*vpqrs(b,apr,lpr,l)-vpqrs(b,l,lpr,apr))
      
      if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._d*vpqrs(a,apr,lpr,k)-vpqrs(a,k,lpr,apr))
      if(kdelta(b,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._d*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._d*vpqrs(b,apr,lpr,k)-vpqrs(b,k,lpr,apr))
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_3=C_4i_3-(2._d*vpqrs(b,apr,kpr,l)-vpqrs(b,l,kpr,apr))
      
      C_4i_3=C_4i_3/sqrt(2._d)
      
    end function C_4i_3

!!$ (4ii,3)-Cabkl,a'a'k'l' (spin case II) (case 8 in Trofimov's Ph.D)    

    real(d) function C_4ii_3(a,b,k,l,apr,kpr,lpr)
    
      integer, intent(in) ::a,b,k,l,apr,kpr,lpr 
    
      real(d) :: vpqrs
    
      external vpqrs

      C_4ii_3=0._d
      
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
      
      
      C_4ii_3=sqrt(3._d)*C_4ii_3/sqrt(2._d)
      
    end function C_4ii_3

!!$ (4ii,4ii)-Cabkl,a'b'k'l' (case 7 in Trofimov's Ph.D)

    real(d) function C_4ii_4ii(a,b,k,l,apr,bpr,kpr,lpr)
    
      integer, intent(in) ::a,b,k,l,apr,bpr,kpr,lpr 
    
      real(d) :: vpqrs
    
      external vpqrs

      C_4ii_4ii=0._d    

      if(kdelta(k,kpr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,apr,b,bpr)-vpqrs(a,bpr,b,apr))
      if(kdelta(a,apr)*kdelta(b,bpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(kpr,k,lpr,l)-vpqrs(kpr,l,lpr,k))
      
      if(kdelta(b,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(a,apr,kpr,k)-1.5_d*vpqrs(a,k,kpr,apr))
      if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(a,apr,lpr,l)-1.5_d*vpqrs(a,l,lpr,apr))
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(b,bpr,kpr,k)-1.5_d*vpqrs(b,k,kpr,bpr))
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(b,bpr,lpr,l)-1.5_d*vpqrs(b,l,lpr,bpr))
      
      if(kdelta(b,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,apr,lpr,k)-1.5_d*vpqrs(a,k,lpr,apr))
      if(kdelta(b,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,apr,kpr,l)-1.5_d*vpqrs(a,l,kpr,apr))
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(b,bpr,lpr,k)-1.5_d*vpqrs(b,k,lpr,bpr))
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(b,bpr,kpr,l)-1.5_d*vpqrs(b,l,kpr,bpr))      

      if(kdelta(b,apr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,bpr,kpr,k)-1.5_d*vpqrs(a,k,kpr,bpr))
      if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(a,bpr,lpr,l)-1.5_d*vpqrs(a,l,lpr,bpr))
      if(kdelta(a,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(b,apr,kpr,k)-1.5_d*vpqrs(b,k,kpr,apr))
      if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii+(vpqrs(b,apr,lpr,l)-1.5_d*vpqrs(b,l,lpr,apr))


      if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(a,bpr,lpr,k)-1.5_d*vpqrs(a,k,lpr,bpr))
      if(kdelta(b,apr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(a,bpr,kpr,l)-1.5_d*vpqrs(a,l,kpr,bpr))
      if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(b,apr,lpr,k)-1.5_d*vpqrs(b,k,lpr,apr))
      if(kdelta(a,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4ii_4ii=C_4ii_4ii-(vpqrs(b,apr,kpr,l)-1.5_d*vpqrs(b,l,kpr,apr))

    end function C_4ii_4ii


!!$ (4i,4ii)-Cabkl,a'b'k'l' (case 2 in Trofimov's Ph.D)
    
    real(d) function C_4i_4ii(a,b,k,l,apr,bpr,kpr,lpr)
      
      integer, intent(in) ::a,b,k,l,apr,bpr,kpr,lpr 
      
      real(d) :: vpqrs
      
      external vpqrs
      
      C_4i_4ii=0._d

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
      
      
      C_4i_4ii=sqrt(3._d)*C_4i_4ii/2._d
      
    end function C_4i_4ii

!!$ (4ii,4i)-Cabkl,a'b'k'l' (case 6 in Trofimov's Ph.D)
    
    real(d) function C_4ii_4i(a,b,k,l,apr,bpr,kpr,lpr)
      
      integer, intent(in) ::a,b,k,l,apr,bpr,kpr,lpr 
      
      real(d) :: vpqrs
      
      external vpqrs
      
      C_4ii_4i=0._d

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
      
      
      C_4ii_4i=sqrt(3._d)*C_4ii_4i/2._d
      
    end function C_4ii_4i


!!$ (4i,4i)-Cabkl,a'b'k'l' (case 1 in Trofimov's Ph.D)

    real(d) function C_4i_4i(a,b,k,l,apr,bpr,kpr,lpr)
      
      integer, intent(in) ::a,b,k,l,apr,bpr,kpr,lpr 
      
      real(d) :: vpqrs
      
      external vpqrs
      
      C_4i_4i=0._d    

      if(kdelta(k,kpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i+2._d*(vpqrs(a,apr,b,bpr)+vpqrs(a,bpr,b,apr))
      if(kdelta(a,apr)*kdelta(b,bpr) .eq. 1)&
           C_4i_4i=C_4i_4i+2._d*(vpqrs(kpr,k,lpr,l)+vpqrs(kpr,l,lpr,k))

      if(kdelta(b,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(a,apr,kpr,k)-vpqrs(a,k,kpr,apr))
      if(kdelta(b,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(a,apr,lpr,l)-vpqrs(a,l,lpr,apr))
      if(kdelta(a,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(b,bpr,kpr,k)-vpqrs(b,k,kpr,bpr))
      if(kdelta(a,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(b,bpr,lpr,l)-vpqrs(b,l,lpr,bpr))
      
      if(kdelta(b,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(a,apr,lpr,k)-vpqrs(a,k,lpr,apr))
      if(kdelta(b,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(a,apr,kpr,l)-vpqrs(a,l,kpr,apr))
      if(kdelta(a,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(b,bpr,lpr,k)-vpqrs(b,k,lpr,bpr))
      if(kdelta(a,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(b,bpr,kpr,l)-vpqrs(b,l,kpr,bpr))
      
      if(kdelta(b,apr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(a,bpr,kpr,k)-vpqrs(a,k,kpr,bpr))
      if(kdelta(b,apr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(a,bpr,lpr,l)-vpqrs(a,l,lpr,bpr))
      if(kdelta(a,bpr)*kdelta(l,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(b,apr,kpr,k)-vpqrs(b,k,kpr,apr))
      if(kdelta(a,bpr)*kdelta(k,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(b,apr,lpr,l)-vpqrs(b,l,lpr,apr)) 

      if(kdelta(b,apr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(a,bpr,lpr,k)-vpqrs(a,k,lpr,bpr))
      if(kdelta(b,apr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(a,bpr,kpr,l)-vpqrs(a,l,kpr,bpr))
      if(kdelta(a,bpr)*kdelta(l,kpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(b,apr,lpr,k)-vpqrs(b,k,lpr,apr))
      if(kdelta(a,bpr)*kdelta(k,lpr) .eq. 1)&
           C_4i_4i=C_4i_4i-(2._d*vpqrs(b,apr,kpr,l)-vpqrs(b,l,kpr,apr))
      
      C_4i_4i=0.5_d*C_4i_4i

    end function C_4i_4i
    

  end module adc_ph
