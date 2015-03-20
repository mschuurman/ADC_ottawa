module dipole_ph
  
!!$ The module contains an assortment of functions needed to calculate
!!$ the modified transition moments between the ground and an excited state.  
!!$For the spin-orbital expressions see A.B. Trofimov et al, JCP 111,9982 (1999).
!!$Spin free expressions were taken from the Ph.D. thesis of A.B. Trofimov.
!!$ We follow the usual convention dy denoting virtuals as a,b,c... and occupied as j,k,l...

  use constants
  use parameters
  use misc
  
  implicit none

  real(d) :: vpqrs
  external vpqrs
  
contains

!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************PH BLOCK, MAIN STATES*************************
!!$*******************************************************************************
!!$*******************************************************************************

  real(d) function F0_ph(a,k)
        
    integer, intent(in) :: a,k
    integer :: b,l,b1,l1,sym1
    real(d) :: e_abkl

    F0_ph=0._d

    do b1= nOcc+1,nBas
       b=roccnum(b1)
       do l1= 1,nOcc
          l=roccnum(l1)
          sym1=MT(orbSym(b),orbSym(l))
          if (sym1 .eq. CHECK_dip) then
             e_abkl=e(a)+e(b)-e(k)-e(l)
             F0_ph=F0_ph+dpl(b,l)/e_abkl*(2._d*vpqrs(a,k,b,l)-vpqrs(a,l,b,k))
          end if
       end do
    end do
    
    F0_ph=-F0_ph
    
  end function F0_ph

!!$------------------------------------------------------------------
  
  real(d) function FA_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,l1,l,m1,m,sym,sym1,sym2
    real(d) :: e_lmba,e_lmbc
    
    FA_ph=0._d
    
    do c1= nOcc+1,nBas
       c=roccnum(c1)
       sym=MT(orbSym(k),orbSym(c))
       if (sym .eq. CHECK_dip) then
          do b1= nOcc+1,nBas
             b=roccnum(b1)
             do l1= 1,nOcc
                l=roccnum(l1)
                
                sym1=MT(orbSym(a),orbSym(b))
                
                if (sym1 .eq. 1) then
                   e_lmba=2._d*e(l)-e(b)-e(a)
                   e_lmbc=2._d*e(l)-e(b)-e(c)
                   FA_ph=FA_ph-0.5_d*dpl(c,k)/e_lmba/e_lmbc*vpqrs(a,l,b,l)*vpqrs(l,c,l,b)
                end if
                
                do m1= l1+1,nOcc
                   m=roccnum(m1)
                   
                   sym2=MT(orbSym(l),orbSym(m))
                   
                   if (MT(sym1,sym2) .eq. 1) then
                   
                      e_lmba=e(l)+e(m)-e(b)-e(a)
                      e_lmbc=e(l)+e(m)-e(b)-e(c)
                      
                      FA_ph=FA_ph-0.5_d*dpl(c,k)/e_lmba/e_lmbc*(&
                           vpqrs(a,l,b,m)*(2._d*vpqrs(l,c,m,b)-vpqrs(l,b,m,c))+&
                           vpqrs(a,m,b,l)*(2._d*vpqrs(l,b,m,c)-vpqrs(l,c,m,b)))                      
                   end if
                   
                end do
                      
             end do
          end do
       end if

    end do
                   
  end function FA_ph
!!$------------------------------------------------------------------
  
  real(d) function FB_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,l1,l,m1,m,sym,sym1,sym2
    real(d) :: e_klbc,e_lmbc
    
    FB_ph=0._d
    
    do m1= 1,nOcc
       m=roccnum(m1)
       sym=MT(orbSym(a),orbSym(m))
       if (sym .eq. CHECK_dip) then
          do l1= 1,nOcc
             l=roccnum(l1)
             do b1= nOcc+1,nBas
                b=roccnum(b1)                

                sym1=MT(orbSym(k),orbSym(l))
                if (sym1 .eq. 1) then
                   e_klbc=e(k)+e(l)-e(b)-e(b)
                   e_lmbc=e(l)+e(m)-e(b)-e(b)                   
                   FB_ph=FB_ph-0.5_d*dpl(a,m)/e_klbc/e_lmbc*vpqrs(b,k,b,l)*vpqrs(m,b,l,b)
                end if

                do c1= b1+1,nBas
                   c=roccnum(c1)
                   sym2=MT(orbSym(b),orbSym(c))
                   
                   if (MT(sym1,sym2) .eq. 1) then
                      e_klbc=e(k)+e(l)-e(b)-e(c)
                      e_lmbc=e(l)+e(m)-e(b)-e(c)                      
                      FB_ph=FB_ph-0.5_d*dpl(a,m)/e_klbc/e_lmbc*(&
                           vpqrs(b,k,c,l)*(2._d*vpqrs(m,b,l,c)-vpqrs(m,c,l,b))+&
                           vpqrs(b,l,c,k)*(2._d*vpqrs(m,c,l,b)-vpqrs(m,b,l,c)))
                      
                   end if
                   
                end do
             end do
          end do
       end if
    end do

  end function FB_ph
!!$------------------------------------------------------------------ 

  real(d) function FC_ph(a,k)

    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,l1,l,m1,m,sym,sym1
    real(d) :: e_klab,e_mlcb
    
    FC_ph=0._d

    do c1= nOcc+1,nBas
       c=roccnum(c1)
       do m1= 1,nOcc
          m=roccnum(m1)
          sym=MT(orbSym(c),orbSym(m))
          if (sym .eq. CHECK_dip) then
             do b1= nOcc+1,nBas
                b=roccnum(b1)
                do l1= 1,nOcc
                   l=roccnum(l1)
                   
                   sym1=MT(orbSym(l),orbSym(b))

                   if (sym1 .eq. CHECK_dip) then
                      
                      e_klab=e(k)+e(l)-e(a)-e(b)
                      e_mlcb=e(m)+e(l)-e(c)-e(b)
                      
                      FC_ph=FC_ph+0.5_d*dpl(c,m)/e_klab/e_mlcb*(&
                           2._d*vpqrs(a,k,b,l)-vpqrs(a,l,b,k))*(&
                           2._d*vpqrs(m,c,l,b)-vpqrs(m,b,l,c))

                   end if
                   
                end do
             end do
          end if
       end do
    end do

  end function FC_ph

!!$------------------------------------------------------------------

  real(d) function F21_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,l1,l,m1,m,sym,sym1,sym2
    real(d) :: e_ma,e_lmbc

    F21_ph=0._d
    
    do m1= 1,nOcc
       m=roccnum(m1)
       sym=MT(orbSym(m),orbSym(k))
          if (sym .eq. CHECK_dip) then
          e_ma=e(m)-e(a)
          
          do l1= 1,nOcc
             l=roccnum(l1)
             do b1= nOcc+1,nBas
                b=roccnum(b1)
                
                sym1=MT(orbSym(l),orbSym(m))

                if (sym1 .eq. 1) then
                   e_lmbc=e(l)+e(m)-e(b)-e(b)
                   F21_ph=F21_ph-dpl(m,k)/e_ma/e_lmbc*vpqrs(b,l,b,m)*vpqrs(l,b,a,b)
                end if

                do c1= b1+1,nBas
                   c=roccnum(c1)
                   
                   sym2=MT(orbSym(c),orbSym(b))
                   
                   if (MT(sym1,sym2) .eq. 1) then
                      e_lmbc=e(l)+e(m)-e(b)-e(c)
                      F21_ph= F21_ph-dpl(m,k)/e_ma/e_lmbc*(&
                           vpqrs(b,l,c,m)*(2._d*vpqrs(l,b,a,c)-vpqrs(l,c,a,b))+&
                           vpqrs(b,m,c,l)*(2._d*vpqrs(l,c,a,b)-vpqrs(l,b,a,c)))
                   end if
           
                end do
             end do
          end do
       end if
    end do
    
    
  end function F21_ph
  
!!$--------------------------------------
!!$--------------------------------------

  real(d) function F22_ph(a,k)

    integer, intent(in) :: a,k
    integer :: b1,b,l1,l,m1,m,n1,n,sym,sym1,sym2
    real(d) :: e_na,e_lmba

    F22_ph=0._d

    do n1= 1,nOcc
       n=roccnum(n1)
       sym=MT(orbSym(k),orbSym(n))
          if (sym .eq. CHECK_dip) then
          e_na=e(n)-e(a)
          do l1= 1,nOcc
             l=roccnum(l1)
             do b1= nOcc+1,nBas
                b=roccnum(b1)
                sym1=MT(orbSym(n),orbSym(b))
             
                if (sym1 .eq. 1) then
                   e_lmba=e(l)+e(l)-e(b)-e(a)
                   
                   F22_ph=F22_ph+dpl(n,k)/e_na/e_lmba*vpqrs(l,b,l,n)*vpqrs(b,l,a,l)
                   
                end if
                
                do m1= l1+1,nOcc
                   m=roccnum(m1)
                   
                   sym2=MT(orbSym(l),orbSym(m))
                   
                   if (MT(sym1,sym2) .eq. 1) then
                      e_lmba=e(l)+e(m)-e(b)-e(a)
                      F22_ph=F22_ph+dpl(n,k)/e_na/e_lmba*(&
                           vpqrs(l,b,m,n)*(2._d*vpqrs(b,l,a,m)-vpqrs(b,m,a,l))+&
                           vpqrs(l,n,m,b)*(2._d*vpqrs(b,m,a,l)-vpqrs(b,l,a,m)))
                   end if
                   
                end do
                
             end do
          end do
       end if
    end do
    
  end function F22_ph

!!$-----------------------------------------

  real(d) function F23_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,l1,l,m1,m,sym,sym1,sym2
    real(d) :: e_kb,e_lmbc

    F23_ph=0._d

    do b1= nOcc+1,nBas
       b=roccnum(b1)
       sym=MT(orbSym(a),orbSym(b))
          if (sym .eq. CHECK_dip) then
          e_kb=e(k)-e(b)
          do c1=nOcc+1,nBas
             c=roccnum(c1)
             do l1= 1,nOcc
                l=roccnum(l1)
                
                sym1=MT(orbSym(b),orbSym(c))   
                
                if (sym1 .eq. 1) then
                   e_lmbc=e(l)+e(l)-e(b)-e(c)
                   F23_ph=F23_ph-dpl(a,b)/e_kb/e_lmbc*vpqrs(b,l,c,l)*vpqrs(l,k,l,c)
                end if
                
                do m1= l1+1,nOcc
                   m=roccnum(m1)
                   
                   sym2=MT(orbSym(l),orbSym(m))
                   
                   if (MT(sym1,sym2) .eq. 1) then

                      e_lmbc=e(l)+e(m)-e(b)-e(c)
                      F23_ph=F23_ph-dpl(a,b)/e_kb/e_lmbc*(&
                           vpqrs(b,l,c,m)*(2._d*vpqrs(l,k,m,c)-vpqrs(l,c,m,k))+&
                           vpqrs(b,m,c,l)*(2._d*vpqrs(l,c,m,k)-vpqrs(l,k,m,c)))
                      
                   end if
                   
                end do
                
             end do
          end do
       end if
    end do

  end function F23_ph

!!$----------------------------------------------

  real(d) function F24_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,d1,dd,l1,l,sym,sym1,sym2
    real(d) :: e_kd,e_klbc
    
    F24_ph=0._d
    
    do d1= nOcc+1,nBas
       dd=roccnum(d1)
       sym=MT(orbSym(a),orbSym(dd))
          if (sym .eq. CHECK_dip) then
          e_kd=e(k)-e(dd)
          do l1= 1,nOcc
             l=roccnum(l1)
             do b1= nOcc+1,nBas
                b=roccnum(b1)
                
                sym1=MT(orbSym(l),orbSym(dd))
                
                if (sym1 .eq. 1) then
                
                   e_klbc=e(k)+e(l)-e(b)-e(b)
                   F24_ph=F24_ph+dpl(a,dd)/e_kd/e_klbc*vpqrs(b,k,b,l)*vpqrs(dd,b,l,b)

                end if

                do c1= b1+1,nBas
                   c=roccnum(c1)
                   
                   sym2=MT(orbSym(b),orbSym(c))  
                   
                   if (MT(sym1,sym2) .eq. 1) then
                      
                      e_klbc=e(k)+e(l)-e(b)-e(c)
                      F24_ph=F24_ph+dpl(a,dd)/e_kd/e_klbc*(&
                           vpqrs(b,k,c,l)*(2._d*vpqrs(dd,b,l,c)-vpqrs(dd,c,l,b))+&
                           vpqrs(b,l,c,k)*(2._d*vpqrs(dd,c,l,b)-vpqrs(dd,b,l,c)))
                      
                   end if
                   
                end do

             end do
          end do
       end if
    end do
       
  end function F24_ph

!!$----------------------------------------------------
  
  real(d) function F25_ph(a,k)

    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,l1,l,m1,m,sym1,sym2
    real(d) :: e_kmac,e_lmbc

    F25_ph=0._d
    
    do m1= 1,nOcc
       m=roccnum(m1)
       do c1= nOcc+1,nBas
          c=roccnum(c1)
          sym1=MT(orbSym(m),orbSym(c))
          if (sym1 .eq. CHECK_dip) then
             e_kmac=e(k)+e(m)-e(a)-e(c)
             do l1= 1,nOcc
                l=roccnum(l1)
                do b1= nOcc+1,nBas
                   b=roccnum(b1)
   
                   sym2=MT(orbSym(l),orbSym(b))
                   
                   if (MT(sym1,sym2) .eq. 1) then
                
                      e_lmbc=e(l)+e(m)-e(b)-e(c)
                
                      F25_ph=F25_ph+dpl(m,c)/e_kmac/e_lmbc*&
                           (2._d*vpqrs(a,k,l,b)-vpqrs(a,b,l,k))*&
                           (2._d*vpqrs(b,l,c,m)-vpqrs(b,m,c,l))

                   end if
                
                end do
             end do
          end if
       end do
    end do

  end function F25_ph

!!$----------------------------------------

  real(d) function F26_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,l1,l,m1,m,sym1,sym2
    real(d) :: e_klab,e_kmac

    F26_ph=0._d

    do m1= 1,nOcc
       m=roccnum(m1)
       do c1= nOcc+1,nBas
          c=roccnum(c1)
          sym1=MT(orbSym(m),orbSym(c))
          if (sym1 .eq. CHECK_dip) then
             e_kmac=e(k)+e(m)-e(a)-e(c)
             do l1= 1,nOcc
                l=roccnum(l1)
                do b1= nOcc+1,nBas
                   b=roccnum(b1)
                   
                   sym2=MT(orbSym(l),orbSym(b))
                   
                   if (MT(sym1,sym2) .eq. 1) then
                      
                      e_klab=e(k)+e(l)-e(a)-e(b)
                   
                      F26_ph=F26_ph+dpl(m,c)/e_klab/e_kmac*&
                           (2._d*vpqrs(l,b,c,m)-vpqrs(l,m,c,b))*&
                           (2._d*vpqrs(a,k,b,l)-vpqrs(a,l,b,k))
                      
                   end if
                end do
             end do
          end if
       end do
    end do

  end function F26_ph
  
!!$----------------------------------

  real(d) function F27_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,l1,l,m1,m,sym1,sym2
    real(d) :: e_klab,e_kmbc

    F27_ph=0._d
    
    do l1= 1,nOcc
       l=roccnum(l1)
       do b1= nOcc+1,nBas
          b=roccnum(b1)  
          sym1=MT(orbSym(l),orbSym(b))
          if (sym1 .eq. CHECK_dip) then
             e_klab=e(k)+e(l)-e(a)-e(b)
             do m1= 1,nOcc
                m=roccnum(m1)
                do c1= nOcc+1,nBas
                   c=roccnum(c1)
                   
                   sym1=MT(orbSym(a),orbSym(l))
                   sym2=MT(orbSym(m),orbSym(c))
                   
                   if (MT(sym1,sym2) .eq. 1) then
                 
                      e_kmbc=e(k)+e(m)-e(b)-e(c)
                      
                      F27_ph=F27_ph-dpl(l,b)/e_klab/e_kmbc*(&
                           vpqrs(a,l,m,c)*(2._d*vpqrs(c,m,b,k)-vpqrs(c,k,b,m))+&
                           vpqrs(a,c,m,l)*(2._d*vpqrs(c,k,b,m)-vpqrs(c,m,b,k)))
                      
                   end if
    
                end do

             end do
          end if
       end do
    end do

  end function F27_ph

!!$--------------------------------------------------------

    real(d) function F28_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,l1,l,m1,m,sym1,sym2
    real(d) :: e_kmab,e_lmac

    F28_ph=0._d
    
    do m1= 1,nOcc
       m=roccnum(m1)
       do b1= nOcc+1,nBas
          b=roccnum(b1)
          sym1=MT(orbSym(m),orbSym(b))
          if (sym1 .eq. CHECK_dip) then
             e_kmab=e(k)+e(m)-e(a)-e(b)
             do l1= 1,nOcc
                l=roccnum(l1)
                do c1= nOcc+1,nBas
                   c=roccnum(c1)
                   
                   sym1=MT(orbSym(a),orbSym(m))
                   sym2=MT(orbSym(c),orbSym(l))
                   
                   if (MT(sym1,sym2) .eq. 1) then
                      
                      e_lmac=e(l)+e(m)-e(a)-e(c)
                      
                      F28_ph=F28_ph-dpl(m,b)/e_kmab/e_lmac*(&
                           vpqrs(a,m,c,l)*(2._d*vpqrs(l,c,b,k)-vpqrs(l,k,b,c))+&
                           vpqrs(a,l,c,m)*(2._d*vpqrs(l,k,b,c)-vpqrs(l,c,b,k)))
                      
                   end if
             
                end do
             
             end do
          end if
       end do
    end do

  end function F28_ph

!!$-------------------------------------------------------------------

  real(d) function F29_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,c1,c,d1,dd,l1,l,sym,sym1,sym2
    real(d) :: e_klab,e_klcd

    F29_ph=0._d

    do b1= nOcc+1,nBas
       b=roccnum(b1)
       do l1= 1,nOcc
          l=roccnum(l1)
          sym=MT(orbSym(l),orbSym(b))
          if (sym .eq. CHECK_dip) then
             e_klab=e(k)+e(l)-e(a)-e(b)
             do c1= nOcc+1,nBas
                c=roccnum(c1)

                sym1=MT(orbSym(a),orbSym(b))
                if (sym1 .eq. 1) then
                   e_klcd=e(k)+e(l)-e(c)-e(c)
                   F29_ph=F29_ph+dpl(l,b)/e_klab/e_klcd*vpqrs(a,c,b,c)*vpqrs(c,k,c,l)
                end if
                   
                do d1= c1+1,nBas
                   dd=roccnum(d1)
                   sym2=MT(orbSym(c),orbSym(dd))
                   
                   if (MT(sym1,sym2) .eq. 1) then
                      
                      e_klcd=e(k)+e(l)-e(c)-e(dd)
                      
                      F29_ph=F29_ph+dpl(l,b)/e_klab/e_klcd*(&
                           vpqrs(a,c,b,dd)*(2._d*vpqrs(c,k,dd,l)-vpqrs(c,l,dd,k))+& 
                           vpqrs(a,dd,b,c)*(2._d*vpqrs(dd,k,c,l)-vpqrs(dd,l,c,k)))
                      
                   end if
                   
                end do

             end do
          end if
       end do
    end do

  end function F29_ph

!!$--------------------------------------------

    real(d) function F210_ph(a,k)
    
    integer, intent(in) :: a,k
    integer :: b1,b,l1,l,m1,m,n1,n,sym,sym1,sym2
    real(d) :: e_klab,e_mnab

    F210_ph=0._d
    
    do b1= nOcc+1,nBas
       b=roccnum(b1)
       do l1= 1,nOcc
          l=roccnum(l1)
          sym=MT(orbSym(l),orbSym(b))
          if (sym .eq. CHECK_dip) then
             e_klab=e(k)+e(l)-e(a)-e(b)
             do m1= 1,nOcc
                m=roccnum(m1)
                
                sym1=MT(orbSym(a),orbSym(b))
                if (sym1 .eq. 1) then
                   e_mnab=e(m)+e(m)-e(a)-e(b)
                   F210_ph=F210_ph+dpl(l,b)/e_klab/e_mnab*vpqrs(a,m,b,m)*vpqrs(m,k,m,l)
                end if

                do n1= m1+1,nOcc
                   n=roccnum(n1)
                   
                   sym2=MT(orbSym(m),orbSym(n))

                   if (MT(sym1,sym2) .eq. 1) then
                      
                      e_mnab=e(m)+e(n)-e(a)-e(b)
                      
                      F210_ph=F210_ph+dpl(l,b)/e_klab/e_mnab*(&
                           vpqrs(a,m,b,n)*(2._d*vpqrs(m,k,n,l)-vpqrs(m,l,n,k))+& 
                           vpqrs(a,n,b,m)*(2._d*vpqrs(n,k,m,l)-vpqrs(n,l,m,k)))
                      
                   end if
                   
                end do
                
             end do
          end if
       end do
    end do
    
  end function F210_ph

!!$-----------------------------------------------------------------
!!$ Transition moments to the doubly excited states are numerated in accordance with our 
!!$ choice for the doubly excited configs order: I-a=b,i=j;II-a|=b,i=j;III a=b,i|=j;IV,V-a|=b,i|=j (1,2)

  real(d) function FI_2p2h(a,k)
    
    integer, intent(in) :: a,k
    integer :: c1,c,m1,m,sym
    real(d) :: e_mkaa,e_kkac
    
    FI_2p2h=0._d
    
    do m1= 1,nOcc
       m=roccnum(m1)
       sym=MT(orbSym(m),orbSym(k))
          if (sym .eq. CHECK_dip) then
          e_mkaa=e(m)+e(k)-e(a)-e(a)
          FI_2p2h=FI_2p2h-dpl(m,k)/e_mkaa*vpqrs(a,k,a,m)
       end if
    end do

    do c1= nOcc+1,nBas
       c=roccnum(c1)
       sym=MT(orbsym(a),orbsym(c))
          if (sym .eq. CHECK_dip) then
          e_kkac=e(k)+e(k)-e(a)-e(c)
          FI_2p2h=FI_2p2h+dpl(a,c)/e_kkac*vpqrs(c,k,a,k)
       end if
    end do

    FI_2p2h=2._d*FI_2p2h
    
  end function FI_2p2h

!!$-------------------------------------------------------------------------------

  real(d) function FII_2p2h(a,b,k)
    
    integer, intent(in) :: a,b,k
    integer :: c1,c,m1,m,sym
    real(d) :: e_mkab,e_kkac,e_kkbc
    
    FII_2p2h=0._d
    
    do m1= 1,nOcc
       m=roccnum(m1)
       sym=MT(orbSym(m),orbSym(k))
          if (sym .eq. CHECK_dip) then
          e_mkab=e(m)+e(k)-e(a)-e(b)
          FII_2p2h=FII_2p2h-dpl(m,k)/e_mkab*(vpqrs(a,m,b,k)+vpqrs(a,k,b,m))
       end if
    end do

    do c1= nOcc+1,nBas
       c=roccnum(c1)
       sym=MT(orbsym(a),orbsym(c))
          if (sym .eq. CHECK_dip) then
          e_kkbc=e(k)+e(k)-e(b)-e(c)
          FII_2p2h=FII_2p2h+dpl(a,c)/e_kkbc*vpqrs(c,k,b,k)
       end if
       sym=MT(orbsym(b),orbsym(c))
          if (sym .eq. CHECK_dip) then
          e_kkac=e(k)+e(k)-e(a)-e(c)
          FII_2p2h=FII_2p2h+dpl(b,c)/e_kkac*vpqrs(c,k,a,k)
       end if
    end do

    FII_2p2h=sqrt(2._d)*FII_2p2h  

  end function FII_2p2h

!!$--------------------------------------------------------------------------------

  real(d) function FIII_2p2h(a,k,l)
    
    integer, intent(in) :: a,k,l
    integer :: c1,c,m1,m,sym
    real(d) :: e_mlaa,e_kmaa,e_klac
    
    FIII_2p2h=0._d
    
    do m1= 1,nOcc
       m=roccnum(m1)
       sym=MT(orbSym(m),orbSym(k))
          if (sym .eq. CHECK_dip) then
          e_mlaa=e(m)+e(l)-e(a)-e(a)
          FIII_2p2h=FIII_2p2h+dpl(m,k)/e_mlaa*vpqrs(a,m,a,l)
       end if
       sym=MT(orbSym(m),orbSym(l))
          if (sym .eq. CHECK_dip) then
          e_kmaa=e(k)+e(m)-e(a)-e(a)
          FIII_2p2h=FIII_2p2h+dpl(m,l)/e_kmaa*vpqrs(a,k,a,m)
       end if
    end do

    do c1= nOcc+1,nBas
       c=roccnum(c1)
       sym=MT(orbsym(a),orbsym(c))
          if (sym .eq. CHECK_dip) then
          e_klac=e(k)+e(l)-e(a)-e(c)
          FIII_2p2h=FIII_2p2h-dpl(a,c)/e_klac*(vpqrs(c,k,a,l)+vpqrs(c,l,a,k))
       end if
    end do

    FIII_2p2h=sqrt(2._d)*FIII_2p2h
    
  end function FIII_2p2h

!!$---------------------------------------------------------------------------------

  real(d) function FIV2_2p2h(a,b,k,l)
    
    integer, intent(in) :: a,b,k,l
    integer :: c1,c,m1,m,sym
    real(d) :: e_mlab,e_kmab,e_klbc,e_klac
    
    FIV2_2p2h=0._d
    
    do m1= 1,nOcc
       m=roccnum(m1)
       sym=MT(orbSym(m),orbSym(k))
          if (sym .eq. CHECK_dip) then
          e_mlab=e(m)+e(l)-e(a)-e(b)
          FIV2_2p2h=FIV2_2p2h+dpl(m,k)/e_mlab*(vpqrs(a,m,b,l)-vpqrs(a,l,m,b))
       end if
       sym=MT(orbSym(m),orbSym(l))
          if (sym .eq. CHECK_dip) then
          e_kmab=e(k)+e(m)-e(a)-e(b)
          FIV2_2p2h=FIV2_2p2h+dpl(m,l)/e_kmab*(vpqrs(a,k,m,b)-vpqrs(a,m,b,k))
       end if
    end do

    do c1= nOcc+1,nBas
       c=roccnum(c1)
       sym=MT(orbsym(a),orbsym(c))
          if (sym .eq. CHECK_dip) then
          e_klbc=e(k)+e(l)-e(b)-e(c)
          FIV2_2p2h=FIV2_2p2h-dpl(a,c)/e_klbc*(vpqrs(c,k,b,l)-vpqrs(c,l,b,k))
       end if
       sym=MT(orbsym(b),orbsym(c))
          if (sym .eq. CHECK_dip) then
          e_klac=e(k)+e(l)-e(a)-e(c)
          FIV2_2p2h=FIV2_2p2h-dpl(b,c)/e_klac*(vpqrs(c,l,a,k)-vpqrs(c,k,a,l))
       end if
    end do

    FIV2_2p2h=sqrt(3._d)*FIV2_2p2h
    
  end function FIV2_2p2h

!!$------------------------------------------------------------------------------------    
    
  real(d) function FIV1_2p2h(a,b,k,l)
    
    integer, intent(in) :: a,b,k,l
    integer :: c1,c,m1,m,sym
    real(d) :: e_mlab,e_kmab,e_klbc,e_klac
    
    FIV1_2p2h=0._d
    
    do m1= 1,nOcc
       m=roccnum(m1)
       sym=MT(orbSym(m),orbSym(k))
          if (sym .eq. CHECK_dip) then
          e_mlab=e(m)+e(l)-e(a)-e(b)
          FIV1_2p2h=FIV1_2p2h+dpl(m,k)/e_mlab*(vpqrs(a,m,b,l)+vpqrs(a,l,m,b))
       end if
       sym=MT(orbSym(m),orbSym(l))
          if (sym .eq. CHECK_dip) then
          e_kmab=e(k)+e(m)-e(a)-e(b)
          FIV1_2p2h=FIV1_2p2h+dpl(m,l)/e_kmab*(vpqrs(a,k,m,b)+vpqrs(a,m,b,k))
       end if
    end do

    do c1= nOcc+1,nBas
       c=roccnum(c1)
       sym=MT(orbsym(a),orbsym(c))
          if (sym .eq. CHECK_dip) then
          e_klbc=e(k)+e(l)-e(b)-e(c)
          FIV1_2p2h=FIV1_2p2h-dpl(a,c)/e_klbc*(vpqrs(c,k,b,l)+vpqrs(c,l,b,k))
       end if
       sym=MT(orbsym(b),orbsym(c))
          if (sym .eq. CHECK_dip) then
          e_klac=e(k)+e(l)-e(a)-e(c)
          FIV1_2p2h=FIV1_2p2h-dpl(b,c)/e_klac*(vpqrs(c,l,a,k)+vpqrs(c,k,a,l))
       end if
    end do

  end function FIV1_2p2h



  subroutine get_Dground_0(Dground_0)
    
    real*8, intent(out) :: Dground_0
    integer :: i
    integer :: cah
    
    Dground_0 = 0.d0

       do i = 1 , hcentre(0)
           cah = hcentre(i)
           Dground_0 = Dground_0 + dpl(cah,cah)
    end do
       
       
        
  end subroutine get_Dground_0


  subroutine get_Dground_2(Dground_2)
    
    real*8, intent(out) :: Dground_2
    integer :: i,ap,a
    integer :: cah
    
    Dground_2 = 0.d0

    do i = 1 , hcentre(0)
           cah = hcentre(i)
           Dground_2 = Dground_2 + dpl(cah,cah)
    end do
       
       
   ! do ap=nOcc+1,nBas
   !    a=roccnum(ap)
   !    do i=1,hcentre(0)
   !        cah=hcentre(i)
   !       Dground_2 = Dground_2 + 
   !    end do       
   ! end do
 


  end subroutine get_Dground_2

!#######################################################################

  real(d) function tauA(c1,a)

    implicit none

    integer, intent(in) :: c1,a
    integer             :: c,b1,b,l1,l,m1,m
    real(d)             :: e_lmba,e_lmbc

    tauA=0.0d0

    c=roccnum(c1)

    do b1= nOcc+1,nBas
       b=roccnum(b1)
       
       do l1=1,nOcc
          l=roccnum(l1)

          e_lmba=2._d*e(l)-e(b)-e(a)
          e_lmbc=2._d*e(l)-e(b)-e(c)

          tauA=tauA-0.5d0/e_lmba/e_lmbc*vpqrs(a,l,b,l)*vpqrs(l,c,l,b)
          
          do m1=l1+1,nOcc
             m=roccnum(m1)
          
             e_lmba=e(l)+e(m)-e(b)-e(a)
             e_lmbc=e(l)+e(m)-e(b)-e(c)

             tauA=tauA-0.5d0/e_lmba/e_lmbc*(&
                  vpqrs(a,l,b,m)*(2._d*vpqrs(l,c,m,b)-vpqrs(l,b,m,c))+&
                  vpqrs(a,m,b,l)*(2._d*vpqrs(l,b,m,c)-vpqrs(l,c,m,b)))

          enddo

       enddo
       
    enddo

  end function tauA

!#######################################################################

  real(d) function tauB(m1,k)

    integer, intent(in) :: k,m1
    integer             :: b1,b,c1,c,l1,l,m,sym,sym1,sym2
    real(d)             :: e_klbc,e_lmbc
    
    tauB=0.0d0
    
    m=roccnum(m1)
    
    do l1=1,nOcc
       l=roccnum(l1)
       do b1=nOcc+1,nBas
          b=roccnum(b1)                
          
          sym1=MT(orbSym(k),orbSym(l))
          if (sym1 .eq. 1) then
             e_klbc=e(k)+e(l)-e(b)-e(b)
             e_lmbc=e(l)+e(m)-e(b)-e(b)                   
             tauB=tauB-0.5d0/e_klbc/e_lmbc*vpqrs(b,k,b&
                  &,l)*vpqrs(m,b,l,b)
          end if
          
          do c1=b1+1,nBas
             c=roccnum(c1)
             sym2=MT(orbSym(b),orbSym(c))
             
             if (MT(sym1,sym2) .eq. 1) then
                e_klbc=e(k)+e(l)-e(b)-e(c)
                e_lmbc=e(l)+e(m)-e(b)-e(c)                      
                tauB=tauB-0.5d0/e_klbc/e_lmbc*(&
                     vpqrs(b,k,c,l)*(2._d*vpqrs(m,b,l,c)-vpqrs(m,c,l,b))+&
                     vpqrs(b,l,c,k)*(2._d*vpqrs(m,c,l,b)-vpqrs(m,b,l,c)))
                
             end if
             
          end do
       end do
    end do

    return

  end function tauB

!#######################################################################
  
  real(d) function tau21(a,m1)

    implicit none

    integer, intent(in) :: a,m1
    integer             :: b1,b,c1,c,l1,l,m,sym,sym1,sym2
    real(d)             :: e_ma,e_lmbc

    tau21=0.0d0

    m=roccnum(m1)

    e_ma=e(m)-e(a)

    do l1= 1,nOcc
       l=roccnum(l1)
       do b1= nOcc+1,nBas
          b=roccnum(b1)
          
          sym1=MT(orbSym(l),orbSym(m))
          
          if (sym1 .eq. 1) then
             e_lmbc=e(l)+e(m)-e(b)-e(b)
             tau21=tau21-1.0d0/e_ma/e_lmbc*vpqrs(b,l,b,m)*vpqrs(l,b,a,b)
          end if
          
          do c1= b1+1,nBas
             c=roccnum(c1)
             
             sym2=MT(orbSym(c),orbSym(b))
             
             if (MT(sym1,sym2) .eq. 1) then
                e_lmbc=e(l)+e(m)-e(b)-e(c)
                tau21= tau21-1.0d0/e_ma/e_lmbc*(&
                     vpqrs(b,l,c,m)*(2._d*vpqrs(l,b,a,c)-vpqrs(l,c,a,b))+&
                     vpqrs(b,m,c,l)*(2._d*vpqrs(l,c,a,b)-vpqrs(l,b,a,c)))
             end if
             
          end do
       end do
    end do

    return

  end function tau21

!#######################################################################

  real(d) function tau22(a,n1)

    implicit none

    integer, intent(in) :: a,n1
    integer             :: b1,b,l1,l,m1,m,n,sym,sym1,sym2
    real(d)             :: e_na,e_lmba

    tau22=0.0d0
    
    n=roccnum(n1)
    
    e_na=e(n)-e(a)
    do l1=1,nOcc
       l=roccnum(l1)
       do b1=nOcc+1,nBas
          b=roccnum(b1)
          sym1=MT(orbSym(n),orbSym(b))
          
          if (sym1 .eq. 1) then
             e_lmba=e(l)+e(l)-e(b)-e(a)
             tau22=tau22+1.0d0/e_na/e_lmba*vpqrs(l,b,l,n)*vpqrs(b,l,a,l)
          end if
                
          do m1= l1+1,nOcc
             m=roccnum(m1)

             sym2=MT(orbSym(l),orbSym(m))

             if (MT(sym1,sym2) .eq. 1) then
                e_lmba=e(l)+e(m)-e(b)-e(a)
                tau22=tau22+1.0d0/e_na/e_lmba*(&
                     vpqrs(l,b,m,n)*(2._d*vpqrs(b,l,a,m)-vpqrs(b,m,a,l))+&
                     vpqrs(l,n,m,b)*(2._d*vpqrs(b,m,a,l)-vpqrs(b,l,a,m)))
             end if
             
          end do
          
       end do
    end do

    return
    
  end function tau22

!#######################################################################

  real(d) function tau23(b1,k)

    implicit none

    integer, intent(in) :: b1,k
    integer             :: b,c1,c,l1,l,m1,m,sym,sym1,sym2
    real(d)             :: e_kb,e_lmbc

    tau23=0.0d0

    b=roccnum(b1)
    e_kb=e(k)-e(b)
    do c1=nOcc+1,nBas
       c=roccnum(c1)
       do l1= 1,nOcc
          l=roccnum(l1)
          
          sym1=MT(orbSym(b),orbSym(c))   
                
          if (sym1 .eq. 1) then
             e_lmbc=e(l)+e(l)-e(b)-e(c)
             tau23=tau23-1.0d0/e_kb/e_lmbc*vpqrs(b,l,c,l)*vpqrs(l,k,l,c)
          end if
          
          do m1= l1+1,nOcc
             m=roccnum(m1)
             
             sym2=MT(orbSym(l),orbSym(m))
             
             if (MT(sym1,sym2) .eq. 1) then
                
                e_lmbc=e(l)+e(m)-e(b)-e(c)
                tau23=tau23-1.0d0/e_kb/e_lmbc*(&
                     vpqrs(b,l,c,m)*(2._d*vpqrs(l,k,m,c)-vpqrs(l,c,m,k))+&
                     vpqrs(b,m,c,l)*(2._d*vpqrs(l,c,m,k)-vpqrs(l,k,m,c)))
                
             end if
             
          end do
          
       end do
    end do

    return

  end function tau23

!#######################################################################

  real(d) function tau24(d1,k)

    implicit none

    integer, intent(in) :: d1,k
    integer             :: b1,b,c1,c,dd,l1,l,sym,sym1,sym2
    real(d)             :: e_kd,e_klbc

    tau24=0.0d0

    dd=roccnum(d1)
    
    e_kd=e(k)-e(dd)
    do l1= 1,nOcc
       l=roccnum(l1)
       do b1= nOcc+1,nBas
          b=roccnum(b1)
          
          sym1=MT(orbSym(l),orbSym(dd))
          
          if (sym1 .eq. 1) then
                
             e_klbc=e(k)+e(l)-e(b)-e(b)
             tau24=tau24+1.0d0/e_kd/e_klbc*vpqrs(b,k,b,l)*vpqrs(dd,b,l,b)
             
          end if
          
          do c1= b1+1,nBas
             c=roccnum(c1)
             
             sym2=MT(orbSym(b),orbSym(c))  
                
             if (MT(sym1,sym2) .eq. 1) then
                
                e_klbc=e(k)+e(l)-e(b)-e(c)
                tau24=tau24+1.0d0/e_kd/e_klbc*(&
                     vpqrs(b,k,c,l)*(2._d*vpqrs(dd,b,l,c)-vpqrs(dd,c,l,b))+&
                     vpqrs(b,l,c,k)*(2._d*vpqrs(dd,c,l,b)-vpqrs(dd,b,l,c)))
                
             end if
                   
          end do
          
       end do
    end do

    return

  end function tau24

!#######################################################################

  real(d) function tauC(b1,l1)

    implicit none

    integer, intent(in) :: b1,l1
    integer             :: b,c1,c,l,m1,m,sym
    real(d)             :: e_klab,e_mlcb
    
    tauC=0.0d0

    b=roccnum(b1)
    l=roccnum(l1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
       do m1=1,nOcc
          m=roccnum(m1)
          sym=MT(orbSym(c),orbSym(m))
          if (sym .eq. CHECK_dip) then
             e_mlcb=e(m)+e(l)-e(c)-e(b)
             tauC=tauC+0.5d0*dpl(c,m)/e_mlcb&
                  *(2.0d0*vpqrs(m,c,l,b)-vpqrs(m,b,l,c))
          endif
       enddo
    enddo
    
    return

  end function tauC

!#######################################################################

  real(d) function FC_ph_new(a,k,tau,nvirt)

    implicit none
    
    integer, intent(in)                        :: a,k,nvirt
    real*8, intent(in), dimension(nvirt,nocc)  :: tau
    integer                                    :: itmp
    integer                                    :: b1,b,l1,l,sym
    real(d) :: e_klab,e_mlcb

    FC_ph_new=0.0d0

    itmp=0
    do b1=nOcc+1,nBas
       itmp=itmp+1
       b=roccnum(b1)
       do l1=1,nOcc
          l=roccnum(l1)
          sym=MT(orbSym(l),orbSym(b))
          if (sym .eq. CHECK_dip) then             
             e_klab=e(k)+e(l)-e(a)-e(b)
             FC_ph_new=FC_ph_new &
                  +tau(itmp,l1)/e_klab*(2.0d0*vpqrs(a,k,b,l)-vpqrs(a,l,b,k))
          endif
       enddo
    enddo

    return

  end function FC_ph_new

!#######################################################################

end module dipole_ph
