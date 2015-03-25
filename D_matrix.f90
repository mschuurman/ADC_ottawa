module D_matrix

!!$The module contains an assortment of functions needed to calculate
!!$matrix elements of the ADC matrix of the polarization propagator.
!!$For the spin-orbital expressions see A.B. Trofimov et al, JCP 111,9982 (1999).
!!$Spin free expressions were taken from the Ph.D. thesis of A.B. Trofimov.
  
  use constants
  use parameters
  use misc
  
  implicit none

  real(d) :: vpqrs
  external vpqrs

  
contains
 

 real(d) function t2_1h1p(apr,kpr)

    integer, intent(in) :: kpr,apr

    integer :: b,c,j,i,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*8 :: DA,DB,eaprbij,ebckprj,eaprkpr,term



  do b1=nOcc+1,nBas
     b=roccnum(b1)
   do c1=nOcc+1,nBas
      c=roccnum(c1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(b),orbSym(c))
             nsym2=MT(orbSym(apr),orbSym(j))
             nsym3=MT(orbSym(kpr),orbSym(j))


             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eaprkpr=e(apr)-e(kpr)
                ebckprj=e(b)+e(c)-e(kpr)-e(j)
                DA=eaprkpr*ebckprj

                term=term+vpqrs(b,apr,c,j)*(+2.0*vpqrs(kpr,b,j,c)-1.0*vpqrs(kpr,c,j,b))
                term=term+vpqrs(b,j,c,apr)*(-1.0*vpqrs(kpr,b,j,c)+2.0*vpqrs(kpr,c,j,b))
                term=term/DA

                t2_1h1p=t2_1h1p+term

             end if
          end do
       end do
    end do


do b1=nOcc+1,nBas
     b=roccnum(b1)
   do i1=1,nOcc
      i=roccnum(i1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(i),orbSym(j))
             nsym2=MT(orbSym(b),orbSym(apr))
             nsym3=MT(orbSym(b),orbSym(kpr))


             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eaprkpr=e(apr)-e(kpr)
                eaprbij=e(apr)+e(b)-e(i)-e(j)
                DB=eaprkpr*eaprbij

                term=term+vpqrs(kpr,i,b,j)*(+2.0*vpqrs(i,apr,j,b)-1.0*vpqrs(i,b,j,apr))
                term=term+vpqrs(kpr,j,b,i)*(-1.0*vpqrs(i,apr,j,b)+2.0*vpqrs(i,b,j,apr))
                term=term/DB

                t2_1h1p=t2_1h1p-term  ! this part has minus sign

             end if
          end do
       end do
    end do


   t2_1h1p = 0.5*t2_1h1p  ! EXPRESSION FACTOR IN THE DENSITY


  end function t2_1h1p



 real(d) function t2_1h1p_hc(apr,kpr)

 integer, intent(in) :: kpr,apr


!!! I only have to write the occupied-unoccupied part of the density matrix !!!
!!!!!!!   THE APR OR KPR HAVE ONLY ONE POSSIBLITY FOR THE SPIN
!!!!!!!   THE FACTOR OF TWO IS TAKEN INTO ACCOUNT IN THE CALLING
!!!!!!!   BY THE FUNCTIONS D2_1,D2_2,D2_3,D2_4.
!!!!!!!   THE FACTOR OF 2 TAKES INTO ACCOUHT THAT OR APR OR KPR
!!!!!!!   CAN TAKE ONE TIME ALPHA AND ONE TIME BETA VALUES AND
!!!!!!!   THE EXPRESSION DOES NOT DEPEND OF WHICH VALUE IT TAKES


    integer :: b,c,j,i,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*8 :: DA,DB,eaprbij,ebckprj,eaprkpr,term

  do b1=nOcc+1,nBas
     b=roccnum(b1)
   do c1=nOcc+1,nBas
      c=roccnum(c1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(b),orbSym(c))
             nsym2=MT(orbSym(apr),orbSym(j))
             nsym3=MT(orbSym(kpr),orbSym(j))


             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eaprkpr=e(apr)-e(kpr)
                ebckprj=e(b)+e(c)-e(kpr)-e(j)
                DA=eaprkpr*ebckprj

                term=term+vpqrs(apr,b,j,c)*(+2.0*vpqrs(b,kpr,c,j)-1.0*vpqrs(c,kpr,b,j))
                term=term+vpqrs(j,b,apr,c)*(-1.0*vpqrs(b,kpr,c,j)+2.0*vpqrs(c,kpr,b,j))
                term=term/DA

                t2_1h1p_hc=t2_1h1p_hc+term

             end if
          end do
       end do
    end do


do b1=nOcc+1,nBas
     b=roccnum(b1)
   do i1=1,nOcc
      i=roccnum(i1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(i),orbSym(j))
             nsym2=MT(orbSym(b),orbSym(apr))
             nsym3=MT(orbSym(b),orbSym(kpr))


             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eaprkpr=e(apr)-e(kpr)
                eaprbij=e(apr)+e(b)-e(i)-e(j)
                DB=eaprkpr*eaprbij

                term=term+vpqrs(i,kpr,j,b)*(+2.0*vpqrs(apr,i,b,j)-1.0*vpqrs(b,i,apr,j))
                term=term+vpqrs(j,kpr,i,b)*(-1.0*vpqrs(apr,i,b,j)+2.0*vpqrs(b,i,apr,j))
                term=term/DB

                t2_1h1p_hc=t2_1h1p_hc-term  ! this part has minus sign

             end if
          end do
       end do
    end do


   t2_1h1p_hc = 0.5*t2_1h1p_hc  ! EXPRESSION FACTOR IN THE DENSITY


  end function t2_1h1p_hc






 real(d) function t2_2h2p(apr,bpr,kpr,lpr)

 integer, intent(in) :: kpr,apr,bpr,lpr


!!! I only have to write the occupied-unoccupied part of the density matrix !!!
!!!!!!!   THE APR OR KPR HAVE ONLY ONE POSSIBLITY FOR THE SPIN
!!!!!!!   THE FACTOR OF TWO IS TAKEN INTO ACCOUNT IN THE CALLING
!!!!!!!   BY THE FUNCTIONS D2_1,D2_2,D2_3,D2_4.
!!!!!!!   THE FACTOR OF 2 TAKES INTO ACCOUHT THAT OR APR OR KPR
!!!!!!!   CAN TAKE ONE TIME ALPHA AND ONE TIME BETA VALUES AND
!!!!!!!   THE EXPRESSION DOES NOT DEPEND OF WHICH VALUE IT TAKES


    integer :: b,c,j,i,nsym1,nsym2,nsym3,nsym4,b1,c1,i1,j1,cnt
    integer :: cpr,mpr,npr,dpr,cpr1,dpr1,mpr1,npr1
    real*8 :: DA,DB,eklab,eikac,term


!!!!!!! virtual-occupied summation


   do cpr1=nOcc+1,nBas
      cpr=roccnum(cpr1)
     do mpr1=1,nOcc
        mpr=roccnum(mpr1)

! first term 
             nsym1=MT(orbSym(bpr),orbSym(mpr))
             nsym2=MT(orbSym(cpr),orbSym(lpr))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(apr),orbSym(cpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(kpr)+e(mpr)-e(apr)-e(cpr)
                DA=eklab*eikac

                term=term+vpqrs(apr,kpr,cpr,mpr)*(+2.0*vpqrs(bpr,lpr,mpr,cpr)-1.0*vpqrs(bpr,cpr,mpr,lpr))
                term=term+vpqrs(apr,mpr,cpr,kpr)*(-1.0*vpqrs(bpr,lpr,mpr,cpr)+2.0*vpqrs(bpr,cpr,mpr,lpr))
                term=term/DA

                t2_2h2p=t2_2h2p+term

             end if

! second term:  exchange   (kpr and lpr)  with respect to the first term
             nsym1=MT(orbSym(bpr),orbSym(mpr))
             nsym2=MT(orbSym(cpr),orbSym(kpr))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(apr),orbSym(cpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(lpr)+e(mpr)-e(apr)-e(cpr)
                DA=eklab*eikac

                term=term+vpqrs(apr,lpr,cpr,mpr)*(+2.0*vpqrs(bpr,kpr,mpr,cpr)-1.0*vpqrs(bpr,cpr,mpr,kpr))
                term=term+vpqrs(apr,mpr,cpr,lpr)*(-1.0*vpqrs(bpr,kpr,mpr,cpr)+2.0*vpqrs(bpr,cpr,mpr,kpr))
                term=term/DA

                t2_2h2p=t2_2h2p-term
             end if
 
! third term:   exchange (apr and bpr)  with respect to the first term
             nsym1=MT(orbSym(apr),orbSym(mpr))
             nsym2=MT(orbSym(cpr),orbSym(lpr))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(bpr),orbSym(cpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(kpr)+e(mpr)-e(bpr)-e(cpr)
                DA=eklab*eikac

                term=term+vpqrs(bpr,kpr,cpr,mpr)*(+2.0*vpqrs(apr,lpr,mpr,cpr)-1.0*vpqrs(apr,cpr,mpr,lpr))
                term=term+vpqrs(bpr,mpr,cpr,kpr)*(-1.0*vpqrs(apr,lpr,mpr,cpr)+2.0*vpqrs(apr,cpr,mpr,lpr))
                term=term/DA

                t2_2h2p=t2_2h2p-term

             end if

! fourth term:  exchange BOTH (abr and bpr) AND  (kpr and lpr) with respect to the first term
             nsym1=MT(orbSym(apr),orbSym(mpr))
             nsym2=MT(orbSym(cpr),orbSym(kpr))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(bpr),orbSym(cpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(lpr)+e(mpr)-e(bpr)-e(cpr)
                DA=eklab*eikac

                term=term+vpqrs(bpr,lpr,cpr,mpr)*(+2.0*vpqrs(apr,kpr,mpr,cpr)-1.0*vpqrs(apr,cpr,mpr,kpr))
                term=term+vpqrs(bpr,mpr,cpr,lpr)*(-1.0*vpqrs(apr,kpr,mpr,cpr)+2.0*vpqrs(apr,cpr,mpr,kpr))
                term=term/DA

                t2_2h2p=t2_2h2p+term

             end if

          end do
       end do


!!!!! occupied-occupied summation


   do mpr1=1,nOcc
      mpr=roccnum(mpr1)
     do npr1=1,nOcc
        npr=roccnum(npr1)

             nsym1=MT(orbSym(mpr),orbSym(npr))
             nsym2=MT(orbSym(bpr),orbSym(apr))
             nsym3=MT(orbSym(lpr),orbSym(kpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(mpr)+e(npr)-e(apr)-e(bpr)
                DB=eklab*eikac

                term=term+vpqrs(apr,mpr,bpr,npr)*(+2.0*vpqrs(mpr,kpr,npr,lpr)-1.0*vpqrs(mpr,lpr,npr,kpr))
                term=term+vpqrs(apr,npr,bpr,mpr)*(-1.0*vpqrs(mpr,kpr,npr,lpr)+2.0*vpqrs(mpr,lpr,npr,kpr))
                term=term/DB

                t2_2h2p=t2_2h2p+term*(0.5) 

             end if
          end do
       end do


!!!!!! virtual-virtual summation

   do cpr1=nOcc+1,nBas
      cpr=roccnum(cpr1)
     do dpr1=nOcc+1,nBas
        dpr=roccnum(dpr1)

             nsym1=MT(orbSym(cpr),orbSym(dpr))
             nsym2=MT(orbSym(apr),orbSym(bpr))
             nsym3=MT(orbSym(kpr),orbSym(lpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(kpr)+e(lpr)-e(cpr)-e(dpr)
                DA=eklab*eikac

                term=term+vpqrs(cpr,kpr,dpr,lpr)*(+2.0*vpqrs(apr,cpr,bpr,dpr)-1.0*vpqrs(apr,dpr,bpr,cpr))
                term=term+vpqrs(cpr,lpr,dpr,kpr)*(-1.0*vpqrs(apr,cpr,bpr,dpr)+2.0*vpqrs(apr,dpr,bpr,cpr))
                term=term/DA

                t2_2h2p=t2_2h2p+term*0.5

             end if
          end do
       end do


  end function t2_2h2p



 real(d) function t2_2h2p_hc(apr,bpr,kpr,lpr,P)

 integer, intent(in) :: kpr,apr,bpr,lpr,P


!!! I only have to write the occupied-unoccupied part of the density matrix !!!
!!!!!!!   THE APR OR KPR HAVE ONLY ONE POSSIBLITY FOR THE SPIN
!!!!!!!   THE FACTOR OF TWO IS TAKEN INTO ACCOUNT IN THE CALLING
!!!!!!!   BY THE FUNCTIONS D2_1,D2_2,D2_3,D2_4.
!!!!!!!   THE FACTOR OF 2 TAKES INTO ACCOUHT THAT OR APR OR KPR
!!!!!!!   CAN TAKE ONE TIME ALPHA AND ONE TIME BETA VALUES AND
!!!!!!!   THE EXPRESSION DOES NOT DEPEND OF WHICH VALUE IT TAKES

    integer :: b,c,j,i,nsym1,nsym2,nsym3,nsym4,b1,c1,i1,j1,cnt
    integer :: cpr,mpr,npr,dpr,cpr1,dpr1,mpr1,npr1
    real*8 :: DA,DB,eklab,eikac,term

!!!!!!! virtual-occupied summation

   do cpr1=nOcc+1,nBas
      cpr=roccnum(cpr1)
     do mpr1=1,nOcc
        mpr=roccnum(mpr1)

! first term 
             nsym1=MT(orbSym(bpr),orbSym(mpr))
             nsym2=MT(orbSym(cpr),orbSym(lpr))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(apr),orbSym(cpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(kpr)+e(mpr)-e(apr)-e(cpr)
                DA=eklab*eikac

              if (P .eq. 1) then 
! from (apr a , bpr lpr)  and  (k kpr, brp lpr)
                term=term+vpqrs(kpr,apr,mpr,cpr)*(+4.0*vpqrs(lpr,bpr,cpr,mpr)-2.0*vpqrs(cpr,bpr,lpr,mpr))
                term=term+vpqrs(mpr,apr,kpr,cpr)*(-2.0*vpqrs(lpr,bpr,cpr,mpr)+1.0*vpqrs(cpr,bpr,lpr,mpr))
                term=term/DA
              else 
! from (apr lpr, bpr a) and   (k lpr, bpr kpr)
                term=term+vpqrs(kpr,apr,mpr,cpr)*(+2.0*vpqrs(lpr,bpr,cpr,mpr)-1.0*vpqrs(cpr,bpr,lpr,mpr))
                term=term+vpqrs(mpr,apr,kpr,cpr)*(-1.0*vpqrs(lpr,bpr,cpr,mpr)+2.0*vpqrs(cpr,bpr,lpr,mpr))
                term=term/DA
              end if

                t2_2h2p_hc=t2_2h2p_hc+term

             end if

! second term:  exchange   (kpr and lpr)  with respect to the first term
             nsym1=MT(orbSym(bpr),orbSym(mpr))
             nsym2=MT(orbSym(cpr),orbSym(kpr))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(apr),orbSym(cpr))
             
         if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(lpr)+e(mpr)-e(apr)-e(cpr)
                DA=eklab*eikac

                if (P .eq. 1) then 
! from  (apr a , bpr lpr)  and (k kpr , bpr lpr)
                term=term+vpqrs(lpr,apr,mpr,cpr)*(+2.0*vpqrs(kpr,bpr,cpr,mpr)-1.0*vpqrs(cpr,bpr,kpr,mpr))
                term=term+vpqrs(mpr,apr,lpr,cpr)*(-1.0*vpqrs(kpr,bpr,cpr,mpr)+2.0*vpqrs(cpr,bpr,kpr,mpr))
                term=term/DA
                else
! from  (apr lpr , bpr a)  and (k lpr , bpr kpr)
                term=term+vpqrs(lpr,apr,mpr,cpr)*(+4.0*vpqrs(kpr,bpr,cpr,mpr)-2.0*vpqrs(cpr,bpr,kpr,mpr))
                term=term+vpqrs(mpr,apr,lpr,cpr)*(-2.0*vpqrs(kpr,bpr,cpr,mpr)+1.0*vpqrs(cpr,bpr,kpr,mpr))
                term=term/DA
                end if 

                t2_2h2p_hc=t2_2h2p_hc-term

             end if

! third term:   exchange (apr and bpr)  with respect to the first term
             nsym1=MT(orbSym(apr),orbSym(mpr))
             nsym2=MT(orbSym(cpr),orbSym(lpr))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(bpr),orbSym(cpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(kpr)+e(mpr)-e(bpr)-e(cpr)
                DA=eklab*eikac

                if ( P .eq. 1) then 
! from  (apr a , bpr lpr)  and (k kpr , bpr lpr)
                term=term+vpqrs(kpr,bpr,mpr,cpr)*(+2.0*vpqrs(lpr,apr,cpr,mpr)-1.0*vpqrs(cpr,apr,lpr,mpr))
                term=term+vpqrs(mpr,bpr,kpr,cpr)*(-1.0*vpqrs(lpr,apr,cpr,mpr)+2.0*vpqrs(cpr,apr,lpr,mpr))
                term=term/DA
                else 
! from  (apr lpr , bpr a)  and (k lpr , bpr kpr)
                term=term+vpqrs(kpr,bpr,mpr,cpr)*(+4.0*vpqrs(lpr,apr,cpr,mpr)-2.0*vpqrs(cpr,apr,lpr,mpr))
                term=term+vpqrs(mpr,bpr,kpr,cpr)*(-2.0*vpqrs(lpr,apr,cpr,mpr)+1.0*vpqrs(cpr,apr,lpr,mpr))
                term=term/DA
                end if

                t2_2h2p_hc=t2_2h2p_hc-term

             end if

! fourth term:  exchange BOTH (abr and bpr) AND  (kpr and lpr) with respect to the first term
             nsym1=MT(orbSym(apr),orbSym(mpr))
             nsym2=MT(orbSym(cpr),orbSym(kpr))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(bpr),orbSym(cpr))
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(lpr)+e(mpr)-e(bpr)-e(cpr)
                DA=eklab*eikac

                if ( P .eq. 1) then
! from  (apr a , bpr lpr)  and (k kpr , bpr lpr)
                term=term+vpqrs(lpr,bpr,mpr,cpr)*(+4.0*vpqrs(kpr,apr,cpr,mpr)-2.0*vpqrs(cpr,apr,kpr,mpr))
                term=term+vpqrs(mpr,bpr,lpr,cpr)*(-2.0*vpqrs(kpr,apr,cpr,mpr)+1.0*vpqrs(cpr,apr,kpr,mpr))
                term=term/DA
                else 
! from  (apr lpr , bpr a)  and (k lpr , bpr kpr)
                term=term+vpqrs(lpr,bpr,mpr,cpr)*(+2.0*vpqrs(kpr,apr,cpr,mpr)-1.0*vpqrs(cpr,apr,kpr,mpr))
                term=term+vpqrs(mpr,bpr,lpr,cpr)*(-1.0*vpqrs(kpr,apr,cpr,mpr)+2.0*vpqrs(cpr,apr,kpr,mpr))
                term=term/DA
                end if 

                t2_2h2p_hc=t2_2h2p_hc+term

             end if

          end do
       end do


!!!!! occupied-occupied summation

   do mpr1=1,nOcc
      mpr=roccnum(mpr1)
     do npr1=1,nOcc
        npr=roccnum(npr1)

             nsym1=MT(orbSym(mpr),orbSym(npr))
             nsym2=MT(orbSym(bpr),orbSym(apr))
             nsym3=MT(orbSym(lpr),orbSym(kpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(mpr)+e(npr)-e(apr)-e(bpr)
                DB=eklab*eikac

                if ( P .eq. 1 ) then
 ! from  (apr a , bpr lpr)  and (k kpr , bpr lpr)
                term=term+vpqrs(mpr,apr,npr,bpr)*(+2.0*vpqrs(kpr,mpr,lpr,npr)-1.0*vpqrs(lpr,mpr,kpr,npr))
                term=term+vpqrs(npr,apr,mpr,bpr)*(-1.0*vpqrs(kpr,mpr,lpr,npr)+2.0*vpqrs(lpr,mpr,kpr,npr))
                term=term/DB
                else
! from  (apr lpr , bpr a)  and (k lpr , bpr kpr)
                term=term+vpqrs(mpr,apr,npr,bpr)*(+1.0*vpqrs(kpr,mpr,lpr,npr)-2.0*vpqrs(lpr,mpr,kpr,npr))
                term=term+vpqrs(npr,apr,mpr,bpr)*(-2.0*vpqrs(kpr,mpr,lpr,npr)+1.0*vpqrs(lpr,mpr,kpr,npr))
                term=term/DB
                end if 

                t2_2h2p_hc=t2_2h2p_hc+term*(0.5) 

             end if
          end do
       end do

!!!!!! virtual-virtual summation

   do cpr1=nOcc+1,nBas
      cpr=roccnum(cpr1)
     do dpr1=nOcc+1,nBas
        dpr=roccnum(dpr1)

             nsym1=MT(orbSym(cpr),orbSym(dpr))
             nsym2=MT(orbSym(apr),orbSym(bpr))
             nsym3=MT(orbSym(kpr),orbSym(lpr))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                cnt=cnt+1
                term=0.0
                
                eklab=e(kpr)+e(lpr)-e(apr)-e(bpr)
                eikac=e(kpr)+e(lpr)-e(cpr)-e(dpr)
                DA=eklab*eikac

                if ( P .eq. 1 ) then
 ! from  (apr a , bpr lpr)  and (k kpr , bpr lpr)
                term=term+vpqrs(kpr,cpr,lpr,dpr)*(+2.0*vpqrs(cpr,apr,dpr,bpr)-1.0*vpqrs(dpr,apr,cpr,bpr))
                term=term+vpqrs(lpr,cpr,kpr,dpr)*(-1.0*vpqrs(cpr,apr,dpr,bpr)+2.0*vpqrs(dpr,apr,cpr,bpr))
                term=term/DA
                else
! from  (apr lpr , bpr a)  and (k lpr , bpr kpr)
                term=term+vpqrs(kpr,cpr,lpr,dpr)*(+1.0*vpqrs(cpr,apr,dpr,bpr)-2.0*vpqrs(dpr,apr,cpr,bpr))
                term=term+vpqrs(lpr,cpr,kpr,dpr)*(-2.0*vpqrs(cpr,apr,dpr,bpr)+1.0*vpqrs(dpr,apr,cpr,bpr))
                term=term/DA
                end if

                t2_2h2p_hc=t2_2h2p_hc+term*0.5

             end if
          end do
       end do

  end function t2_2h2p_hc



 real(d) function t2_3h3p_hc(apr,bpr,cpr,kpr,lpr,mpr,I)

 integer, intent(in) :: kpr,apr,bpr,lpr,cpr,mpr,I


!!! I only have to write the occupied-unoccupied part of the density matrix !!!
!!!!!!!   THE APR OR KPR HAVE ONLY ONE POSSIBLITY FOR THE SPIN
!!!!!!!   THE FACTOR OF TWO IS TAKEN INTO ACCOUNT IN THE CALLING
!!!!!!!   BY THE FUNCTIONS D2_1,D2_2,D2_3,D2_4.
!!!!!!!   THE FACTOR OF 2 TAKES INTO ACCOUHT THAT OR APR OR KPR
!!!!!!!   CAN TAKE ONE TIME ALPHA AND ONE TIME BETA VALUES AND
!!!!!!!   THE EXPRESSION DOES NOT DEPEND OF WHICH VALUE IT TAKES

    integer :: d,f,nsym1,nsym2,nsym3,nsym4,d1,f1,cnt
    real*8 :: DA,eabcklm,eterm,term


!!! part with summation over OCCUPIED index

      do f1=1,nOcc
        f=roccnum(f1)

!!! first term
             nsym1=MT(orbSym(apr),orbSym(bpr))
             nsym2=MT(orbSym(kpr),orbSym(f))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(cpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(kpr)+e(f)-e(apr)-e(bpr)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)    bpr=lpr  cpr=mpr  apr=alpha  kpr=alpha
                term=term+vpqrs(cpr,lpr,f,mpr)*(+2.0*vpqrs(apr,kpr,bpr,f)-1.0*vpqrs(apr,f,bpr,kpr))
                term=term+vpqrs(cpr,mpr,f,lpr)*(-4.0*vpqrs(apr,kpr,bpr,f)+2.0*vpqrs(apr,f,bpr,kpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)    bpr=mpr  cpr=lpr  apr=alpha  kpr=alpha
                term=term+vpqrs(cpr,lpr,f,mpr)*(+4.0*vpqrs(apr,kpr,bpr,f)-2.0*vpqrs(apr,f,bpr,kpr))
                term=term+vpqrs(cpr,mpr,f,lpr)*(-2.0*vpqrs(apr,kpr,bpr,f)+1.0*vpqrs(apr,f,bpr,kpr))
                term=term/DA
              end if

                t2_3h3p_hc=t2_3h3p_hc+term

              end if

!!! second term (cpr and bpr) with respect to the first
             nsym1=MT(orbSym(apr),orbSym(cpr))
             nsym2=MT(orbSym(kpr),orbSym(f))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(bpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(kpr)+e(f)-e(apr)-e(cpr)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(bpr,lpr,f,mpr)*(+4.0*vpqrs(apr,kpr,cpr,f)-2.0*vpqrs(apr,f,cpr,kpr))
                term=term+vpqrs(bpr,mpr,f,lpr)*(-2.0*vpqrs(apr,kpr,cpr,f)+1.0*vpqrs(apr,f,cpr,kpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(bpr,lpr,f,mpr)*(+2.0*vpqrs(apr,kpr,cpr,f)-1.0*vpqrs(apr,f,cpr,kpr))
                term=term+vpqrs(bpr,mpr,f,lpr)*(-4.0*vpqrs(apr,kpr,cpr,f)+2.0*vpqrs(apr,f,cpr,kpr))
                term=term/DA
              end if

                t2_3h3p_hc=t2_3h3p_hc-term

              end if

!!! third term (apr and bpr) with respect to the second
             nsym1=MT(orbSym(bpr),orbSym(cpr))
             nsym2=MT(orbSym(kpr),orbSym(f))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(apr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(kpr)+e(f)-e(bpr)-e(cpr)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(apr,lpr,f,mpr)*(+2.0*vpqrs(bpr,kpr,cpr,f)-1.0*vpqrs(bpr,f,cpr,kpr))
                term=term+vpqrs(apr,mpr,f,lpr)*(-1.0*vpqrs(bpr,kpr,cpr,f)+2.0*vpqrs(bpr,f,cpr,kpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(apr,lpr,f,mpr)*(+1.0*vpqrs(bpr,kpr,cpr,f)-2.0*vpqrs(bpr,f,cpr,kpr))
                term=term+vpqrs(apr,mpr,f,lpr)*(-2.0*vpqrs(bpr,kpr,cpr,f)+1.0*vpqrs(bpr,f,cpr,kpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc+term
               
             end if

! fourth term (kpr and lpr) with respect to the first
             nsym1=MT(orbSym(apr),orbSym(bpr))
             nsym2=MT(orbSym(lpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(cpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(lpr)+e(f)-e(apr)-e(bpr)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(cpr,kpr,f,mpr)*(+1.0*vpqrs(apr,lpr,bpr,f)-2.0*vpqrs(apr,f,bpr,lpr))
                term=term+vpqrs(cpr,mpr,f,kpr)*(-2.0*vpqrs(apr,lpr,bpr,f)+4.0*vpqrs(apr,f,bpr,lpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(cpr,kpr,f,mpr)*(+2.0*vpqrs(apr,lpr,bpr,f)-1.0*vpqrs(apr,f,bpr,lpr))
                term=term+vpqrs(cpr,mpr,f,kpr)*(-1.0*vpqrs(apr,lpr,bpr,f)+2.0*vpqrs(apr,f,bpr,lpr))
                term=term/DA
              end if

                t2_3h3p_hc=t2_3h3p_hc-term

              end if

! fifth term (bpr and cpr) with respect to the fourth
             nsym1=MT(orbSym(apr),orbSym(cpr))
             nsym2=MT(orbSym(lpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(bpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(lpr)+e(f)-e(apr)-e(cpr)
                DA = eabcklm*eterm

               if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(bpr,kpr,f,mpr)*(+2.0*vpqrs(apr,lpr,cpr,f)-1.0*vpqrs(apr,f,cpr,lpr))
                term=term+vpqrs(bpr,mpr,f,kpr)*(-1.0*vpqrs(apr,lpr,cpr,f)+2.0*vpqrs(apr,f,cpr,lpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(bpr,kpr,f,mpr)*(+1.0*vpqrs(apr,lpr,cpr,f)-2.0*vpqrs(apr,f,cpr,lpr))
                term=term+vpqrs(bpr,mpr,f,kpr)*(-2.0*vpqrs(apr,lpr,cpr,f)+4.0*vpqrs(apr,f,cpr,lpr))
                term=term/DA
               
             end if

               t2_3h3p_hc=t2_3h3p_hc+term

              end if

! sixth term (apr and bpr) with respect to the fifth
             nsym1=MT(orbSym(bpr),orbSym(cpr))
             nsym2=MT(orbSym(lpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(apr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(lpr)+e(f)-e(bpr)-e(cpr)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(apr,kpr,f,mpr)*(+4.0*vpqrs(bpr,lpr,cpr,f)-2.0*vpqrs(bpr,f,cpr,lpr))
                term=term+vpqrs(apr,mpr,f,kpr)*(-2.0*vpqrs(bpr,lpr,cpr,f)+1.0*vpqrs(bpr,f,cpr,lpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(apr,kpr,f,mpr)*(+2.0*vpqrs(bpr,lpr,cpr,f)-4.0*vpqrs(bpr,f,cpr,lpr))
                term=term+vpqrs(apr,mpr,f,kpr)*(-1.0*vpqrs(bpr,lpr,cpr,f)+2.0*vpqrs(bpr,f,cpr,lpr))
                term=term/DA

              end if

               t2_3h3p_hc=t2_3h3p_hc-term

              end if

! seventh term (lpr and mpr) with respect to the fourth
             nsym1=MT(orbSym(apr),orbSym(bpr))
             nsym2=MT(orbSym(mpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(lpr))
             nsym4=MT(orbSym(cpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(mpr)+e(f)-e(apr)-e(bpr)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(cpr,kpr,f,lpr)*(+2.0*vpqrs(apr,mpr,bpr,f)-1.0*vpqrs(apr,f,bpr,mpr))
                term=term+vpqrs(cpr,lpr,f,kpr)*(-1.0*vpqrs(apr,mpr,bpr,f)+2.0*vpqrs(apr,f,bpr,mpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(cpr,kpr,f,lpr)*(+1.0*vpqrs(apr,mpr,bpr,f)-2.0*vpqrs(apr,f,bpr,mpr))
                term=term+vpqrs(cpr,lpr,f,kpr)*(-2.0*vpqrs(apr,mpr,bpr,f)+4.0*vpqrs(apr,f,bpr,mpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc+term

              end if

! eighth term (lpr and mpr) with respect to the fifth
             nsym1=MT(orbSym(apr),orbSym(cpr))
             nsym2=MT(orbSym(mpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(lpr))
             nsym4=MT(orbSym(bpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(mpr)+e(f)-e(apr)-e(cpr)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(bpr,kpr,f,lpr)*(+1.0*vpqrs(apr,mpr,cpr,f)-2.0*vpqrs(apr,f,cpr,mpr))
                term=term+vpqrs(bpr,lpr,f,kpr)*(-2.0*vpqrs(apr,mpr,cpr,f)+4.0*vpqrs(apr,f,cpr,mpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(bpr,kpr,f,lpr)*(+2.0*vpqrs(apr,mpr,cpr,f)-1.0*vpqrs(apr,f,cpr,mpr))
                term=term+vpqrs(bpr,lpr,f,kpr)*(-1.0*vpqrs(apr,mpr,cpr,f)+2.0*vpqrs(apr,f,cpr,mpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc-term

              end if

! nineth term (lpr and mpr) with respect to the sixth
             nsym1=MT(orbSym(bpr),orbSym(cpr))
             nsym2=MT(orbSym(mpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(lpr))
             nsym4=MT(orbSym(apr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(mpr)+e(f)-e(bpr)-e(cpr)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(apr,kpr,f,lpr)*(+2.0*vpqrs(bpr,mpr,cpr,f)-4.0*vpqrs(bpr,f,cpr,mpr))
                term=term+vpqrs(apr,lpr,f,kpr)*(-1.0*vpqrs(bpr,mpr,cpr,f)+2.0*vpqrs(bpr,f,cpr,mpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(apr,kpr,f,lpr)*(+4.0*vpqrs(bpr,mpr,cpr,f)-2.0*vpqrs(bpr,f,cpr,mpr))
                term=term+vpqrs(apr,lpr,f,kpr)*(-2.0*vpqrs(bpr,mpr,cpr,f)+1.0*vpqrs(bpr,f,cpr,mpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc+term

             end if

      end do
!!! spin-free expressions written!

!!! part with summation over  VIRTUAL  index

      do d1=nOcc+1,nBas
         d=roccnum(d1)

!!! first term
             nsym1=MT(orbSym(mpr),orbSym(d))
             nsym2=MT(orbSym(bpr),orbSym(cpr))
             nsym3=MT(orbSym(lpr),orbSym(kpr))
             nsym4=MT(orbSym(apr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(kpr)+e(lpr)-e(apr)-e(d)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)     bpr=lpr  cpr=mpr  apr=alpha  kpr=alpha
                term=term+vpqrs(bpr,d,cpr,mpr)*(+4.0*vpqrs(apr,kpr,d,lpr)-2.0*vpqrs(apr,lpr,d,kpr))
                term=term+vpqrs(bpr,mpr,cpr,d)*(-2.0*vpqrs(apr,kpr,d,lpr)+1.0*vpqrs(apr,lpr,d,kpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)     bpr=mpr  cpr=lpr  apr=alpha  kpr=alpha
                term=term+vpqrs(bpr,d,cpr,mpr)*(+2.0*vpqrs(apr,kpr,d,lpr)-1.0*vpqrs(apr,lpr,d,kpr))
                term=term+vpqrs(bpr,mpr,cpr,d)*(-4.0*vpqrs(apr,kpr,d,lpr)+2.0*vpqrs(apr,lpr,d,kpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc+term  

             end if

!!! second term (lpr and mpr) with respect to the first
             nsym1=MT(orbSym(lpr),orbSym(d))
             nsym2=MT(orbSym(bpr),orbSym(cpr))
             nsym3=MT(orbSym(mpr),orbSym(kpr))
             nsym4=MT(orbSym(apr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(kpr)+e(mpr)-e(apr)-e(d)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(bpr,d,cpr,lpr)*(+2.0*vpqrs(apr,kpr,d,mpr)-1.0*vpqrs(apr,mpr,d,kpr))
                term=term+vpqrs(bpr,lpr,cpr,d)*(-4.0*vpqrs(apr,kpr,d,mpr)+2.0*vpqrs(apr,mpr,d,kpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(bpr,d,cpr,lpr)*(+4.0*vpqrs(apr,kpr,d,mpr)-2.0*vpqrs(apr,mpr,d,kpr))
                term=term+vpqrs(bpr,lpr,cpr,d)*(-2.0*vpqrs(apr,kpr,d,mpr)+1.0*vpqrs(apr,mpr,d,kpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc-term  ! this part has minus sign

             end if

!!! third  term  (kpr and lpr) with respect to the second
             nsym1=MT(orbSym(kpr),orbSym(d))
             nsym2=MT(orbSym(bpr),orbSym(cpr))
             nsym3=MT(orbSym(mpr),orbSym(lpr))
             nsym4=MT(orbSym(apr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(lpr)+e(mpr)-e(apr)-e(d)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(bpr,d,cpr,kpr)*(+1.0*vpqrs(apr,lpr,d,mpr)-2.0*vpqrs(apr,mpr,d,lpr))
                term=term+vpqrs(bpr,kpr,cpr,d)*(-2.0*vpqrs(apr,lpr,d,mpr)+1.0*vpqrs(apr,mpr,d,lpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(bpr,d,cpr,kpr)*(+2.0*vpqrs(apr,lpr,d,mpr)-1.0*vpqrs(apr,mpr,d,lpr))
                term=term+vpqrs(bpr,kpr,cpr,d)*(-1.0*vpqrs(apr,lpr,d,mpr)+2.0*vpqrs(apr,mpr,d,lpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc+term  

             end if

!!! fourth term (apr and bpr) with respect to the first
             nsym1=MT(orbSym(mpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(cpr))
             nsym3=MT(orbSym(lpr),orbSym(kpr))
             nsym4=MT(orbSym(bpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(kpr)+e(lpr)-e(bpr)-e(d)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(apr,d,cpr,mpr)*(+2.0*vpqrs(bpr,kpr,d,lpr)-4.0*vpqrs(bpr,lpr,d,kpr))
                term=term+vpqrs(apr,mpr,cpr,d)*(-1.0*vpqrs(bpr,kpr,d,lpr)+2.0*vpqrs(bpr,lpr,d,kpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(apr,d,cpr,mpr)*(+1.0*vpqrs(bpr,kpr,d,lpr)-2.0*vpqrs(bpr,lpr,d,kpr))
                term=term+vpqrs(apr,mpr,cpr,d)*(-2.0*vpqrs(bpr,kpr,d,lpr)+1.0*vpqrs(bpr,lpr,d,kpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc-term   ! this part as a minus sign  

             end if

!!! fifth term (apr and bpr) with respect to the second
             nsym1=MT(orbSym(lpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(cpr))
             nsym3=MT(orbSym(mpr),orbSym(kpr))
             nsym4=MT(orbSym(bpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(kpr)+e(mpr)-e(bpr)-e(d)
                DA = eabcklm*eterm

             if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(apr,d,cpr,lpr)*(+1.0*vpqrs(bpr,kpr,d,mpr)-2.0*vpqrs(bpr,mpr,d,kpr))
                term=term+vpqrs(apr,lpr,cpr,d)*(-2.0*vpqrs(bpr,kpr,d,mpr)+1.0*vpqrs(bpr,mpr,d,kpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(apr,d,cpr,lpr)*(+2.0*vpqrs(bpr,kpr,d,mpr)-4.0*vpqrs(bpr,mpr,d,kpr))
                term=term+vpqrs(apr,lpr,cpr,d)*(-1.0*vpqrs(bpr,kpr,d,mpr)+2.0*vpqrs(bpr,mpr,d,kpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc+term  

             end if

!!! sixth term (apr and bpr) with respect to the third
             nsym1=MT(orbSym(kpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(cpr))
             nsym3=MT(orbSym(mpr),orbSym(lpr))
             nsym4=MT(orbSym(bpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(lpr)+e(mpr)-e(bpr)-e(d)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(apr,d,cpr,kpr)*(+2.0*vpqrs(bpr,lpr,d,mpr)-1.0*vpqrs(bpr,mpr,d,lpr))
                term=term+vpqrs(apr,kpr,cpr,d)*(-4.0*vpqrs(bpr,lpr,d,mpr)+2.0*vpqrs(bpr,mpr,d,lpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(apr,d,cpr,kpr)*(+1.0*vpqrs(bpr,lpr,d,mpr)-2.0*vpqrs(bpr,mpr,d,lpr))
                term=term+vpqrs(apr,kpr,cpr,d)*(-2.0*vpqrs(bpr,lpr,d,mpr)+4.0*vpqrs(bpr,mpr,d,lpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc-term  ! this part has minus sign

             end if

!!! seventh term (bpr and cpr) with respect to the fourth
             nsym1=MT(orbSym(mpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(bpr))
             nsym3=MT(orbSym(lpr),orbSym(kpr))
             nsym4=MT(orbSym(cpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(kpr)+e(lpr)-e(cpr)-e(d)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(apr,d,bpr,mpr)*(+1.0*vpqrs(cpr,kpr,d,lpr)-2.0*vpqrs(cpr,lpr,d,kpr))
                term=term+vpqrs(apr,mpr,bpr,d)*(-2.0*vpqrs(cpr,kpr,d,lpr)+1.0*vpqrs(cpr,lpr,d,kpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(apr,d,bpr,mpr)*(+2.0*vpqrs(cpr,kpr,d,lpr)-4.0*vpqrs(cpr,lpr,d,kpr))
                term=term+vpqrs(apr,mpr,bpr,d)*(-1.0*vpqrs(cpr,kpr,d,lpr)+2.0*vpqrs(cpr,lpr,d,kpr))
                term=term/DA
               end if

               t2_3h3p_hc=t2_3h3p_hc+term  

             end if

!!! eighth term (lpr and mpr) with respect to the seventh
             nsym1=MT(orbSym(lpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(bpr))
             nsym3=MT(orbSym(mpr),orbSym(kpr))
             nsym4=MT(orbSym(cpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(kpr)+e(mpr)-e(cpr)-e(d)
                DA = eabcklm*eterm

              if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(apr,d,bpr,lpr)*(+2.0*vpqrs(cpr,kpr,d,mpr)-4.0*vpqrs(cpr,mpr,d,kpr))
                term=term+vpqrs(apr,lpr,bpr,d)*(-1.0*vpqrs(cpr,kpr,d,mpr)+2.0*vpqrs(cpr,mpr,d,kpr))
                term=term/DA
              else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(apr,d,bpr,lpr)*(+1.0*vpqrs(cpr,kpr,d,mpr)-2.0*vpqrs(cpr,mpr,d,kpr))
                term=term+vpqrs(apr,lpr,bpr,d)*(-2.0*vpqrs(cpr,kpr,d,mpr)+1.0*vpqrs(cpr,mpr,d,kpr))
                term=term/DA
              end if

               t2_3h3p_hc=t2_3h3p_hc-term   ! this part as a minus sign  

             end if

!!! nineth term  (kpr and lpr) with respect to the eighth
             nsym1=MT(orbSym(kpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(bpr))
             nsym3=MT(orbSym(mpr),orbSym(lpr))
             nsym4=MT(orbSym(cpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                cnt=cnt+1
                term=0.0

                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)
                eterm = e(lpr)+e(mpr)-e(cpr)-e(d)
                DA = eabcklm*eterm

               if (I .eq. 1) then
!!! from ( bpr lpr, cpr mpr)
                term=term+vpqrs(apr,d,bpr,kpr)*(+1.0*vpqrs(cpr,lpr,d,mpr)-2.0*vpqrs(cpr,mpr,d,lpr))
                term=term+vpqrs(apr,kpr,bpr,d)*(-2.0*vpqrs(cpr,lpr,d,mpr)+4.0*vpqrs(cpr,mpr,d,lpr))
                term=term/DA
               else
!!! from ( bpr mpr, cpr lpr)
                term=term+vpqrs(apr,d,bpr,kpr)*(+2.0*vpqrs(cpr,lpr,d,mpr)-1.0*vpqrs(cpr,mpr,d,lpr))
                term=term+vpqrs(apr,kpr,bpr,d)*(-4.0*vpqrs(cpr,lpr,d,mpr)+2.0*vpqrs(cpr,mpr,d,lpr))
                term=term/DA
               end if
               t2_3h3p_hc=t2_3h3p_hc+term  
             end if
      end do

!!!!! completed the spin-free expressions !!!!!

  end function t2_3h3p_hc


!!!*****************************************************************************
!!!*****************************************************************************
!!!*****************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!  DENSITY FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!*****************************************************************************
!!!*****************************************************************************
!!!*****************************************************************************

 
 real(d) function density(k,a)
 
 integer, intent(in) :: k,a

!!! I only have to write the occupied-unoccupied part of the density matrix !!!
!!!!!!!   THE APR OR KPR HAVE ONLY ONE POSSIBLITY FOR THE SPIN
!!!!!!!   THE FACTOR OF TWO IS TAKEN INTO ACCOUNT IN THE CALLING
!!!!!!!   BY THE FUNCTIONS D2_1,D2_2,D2_3,D2_4.
!!!!!!!   THE FACTOR OF 2 TAKES INTO ACCOUHT THAT OR APR OR KPR
!!!!!!!   CAN TAKE ONE TIME ALPHA AND ONE TIME BETA VALUES AND
!!!!!!!   THE EXPRESSION DOES NOT DEPEND OF WHICH VALUE IT TAKES

    integer :: b,c,j,i,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    integer :: apr,bpr,cpr,kpr,lpr,mpr,apr1,bpr1,cpr1,kpr1,lpr1,mpr1
    real*8 :: DA,DB,eabij,ebckj,eak,term,term1,term2,eka
    real*8 :: ekkpraapr,ebprcprlprmpr

!!!!!!!!!!!!!!!!!!!!!!!!!! ZEROTH ORDER PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! IT DOES NOT GIVE ANY CONTRIBUTION TO THE HOLE-PARTICLE PART !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! SECOND ORDER PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    density=0.0
    cnt=0
  do b1=nOcc+1,nBas
     b=roccnum(b1)
   do c1=nOcc+1,nBas
      c=roccnum(c1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(b),orbSym(c))
             nsym2=MT(orbSym(a),orbSym(j))
             nsym3=MT(orbSym(k),orbSym(j))
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1
                term=0.0

                eak=e(a)-e(k)
                ebckj=e(b)+e(c)-e(k)-e(j)
                DA=eak*ebckj

                term=term+vpqrs(b,a,c,j)*(+2.0*vpqrs(k,b,j,c)-1.0*vpqrs(k,c,j,b))
                term=term+vpqrs(b,j,c,a)*(-1.0*vpqrs(k,b,j,c)+2.0*vpqrs(k,c,j,b))
                term=term/DA
                
                density = density + term
                
             end if
          end do
       end do
    end do 

do b1=nOcc+1,nBas
     b=roccnum(b1)
   do i1=1,nOcc
      i=roccnum(i1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(i),orbSym(j))
             nsym2=MT(orbSym(b),orbSym(a))
             nsym3=MT(orbSym(b),orbSym(k))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eak=e(a)-e(k)
                eabij=e(a)+e(b)-e(i)-e(j)
                DB=eak*eabij

                term=term+vpqrs(k,i,b,j)*(+2.0*vpqrs(i,a,j,b)-1.0*vpqrs(i,b,j,a))
                term=term+vpqrs(k,j,b,i)*(-1.0*vpqrs(i,a,j,b)+2.0*vpqrs(i,b,j,a))
                term=term/DB

                density = density - term  ! this part has minus sign

             end if
          end do
       end do
    end do

   density = 0.5*density  ! EXPRESSION FACTOR IN THE DENSITY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  END SECOND ORDER PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  THIRD ORDER PART  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(denord .gt. 2) then

!!! t2_1h1p part !!!!!!!!!!!
 do apr1=nOcc+1,nBas
     apr=roccnum(apr1)
   do kpr1=1,nOcc
      kpr=roccnum(kpr1)

             nsym1=MT(orbSym(kpr),orbSym(k))
             nsym2=MT(orbSym(apr),orbSym(a))

      if(MT(nsym1,nsym2) .eq. 1)  then

         cnt=cnt+1
         term=0.0

         ekkpraapr=e(k)+e(kpr)-e(a)-e(apr)

         term=term+(+2.0*vpqrs(k,a,kpr,apr)-1.0*vpqrs(k,apr,kpr,a))*t2_1h1p(apr,kpr)
         term=term/ekkpraapr 

         density=density+term

     end if
   end do
 end do


!!!  end t2_1h1p part !!!!!!!!!!!!!!!!!!!!!!!!!!

!!! t2_1h1p_hc part !!!!!!!!!!!
 do apr1=nOcc+1,nBas
     apr=roccnum(apr1)
   do kpr1=1,nOcc
      kpr=roccnum(kpr1)

             nsym1=MT(orbSym(kpr),orbSym(k))
             nsym2=MT(orbSym(apr),orbSym(a))

      if(MT(nsym1,nsym2) .eq. 1)  then

         cnt=cnt+1
         term=0.0

         eka=e(k)-e(a)

         term=term+(+2.0*vpqrs(k,a,apr,kpr)-1.0*vpqrs(k,kpr,apr,a))*t2_1h1p_hc(apr,kpr)
         term=term/eka 

         density=density+term

     end if
   end do
 end do
!!! end t2_1h1p_hc_part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! t2_2h2p_hc part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 do apr1=nOcc+1,nBas
     apr=roccnum(apr1)
  do bpr1=nOcc+1,nBas
     bpr=roccnum(bpr1)
   do lpr1=1,nOcc
      lpr=roccnum(lpr1)

             nsym1=MT(orbSym(lpr),orbSym(a))
             nsym2=MT(orbSym(apr),orbSym(bpr))

      if(MT(nsym1,nsym2) .eq. 1)  then

         cnt=cnt+1
         term=0.0

         eka=e(k)-e(a)

         term=term+(+1.0*vpqrs(apr,a,bpr,lpr))*t2_2h2p_hc(apr,bpr,k,lpr,1)
         term=term+(-1.0*vpqrs(apr,lpr,bpr,a))*t2_2h2p_hc(apr,bpr,k,lpr,2)
         term=term/eka 

         density=density+term*(0.5)   !!! PLUS PART

         end if
      end do
    end do
  end do


 do bpr1=nOcc+1,nBas
     bpr=roccnum(bpr1)
  do kpr1=1,nOcc
     kpr=roccnum(kpr1)
   do lpr1=1,nOcc
      lpr=roccnum(lpr1)

             nsym1=MT(orbSym(bpr),orbSym(k))
             nsym2=MT(orbSym(kpr),orbSym(lpr))

      if(MT(nsym1,nsym2) .eq. 1)  then

         cnt=cnt+1
         term=0.0

         eka=e(k)-e(a)

         term=term+(+1.0*vpqrs(k,kpr,bpr,lpr))*t2_2h2p_hc(a,bpr,kpr,lpr,1)
         term=term+(-1.0*vpqrs(k,lpr,bpr,kpr))*t2_2h2p_hc(a,bpr,kpr,lpr,2)
         term=term/eka 

         density=density-term*(0.5)  !!! MINUS PART

         end if
       end do
     end do
   end do
!!! end t2_2h2p_hc_part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! t2_3h3p_hc part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do bpr1=nOcc+1,nBas
     bpr=roccnum(bpr1)
  do cpr1=nOcc+1,nBas
     cpr=roccnum(cpr1)
   do lpr1=1,nOcc
      lpr=roccnum(lpr1)
    do mpr1=1,nOcc
       mpr=roccnum(mpr1)

             nsym1=MT(orbSym(bpr),orbSym(cpr))
             nsym2=MT(orbSym(mpr),orbSym(lpr))

      if(MT(nsym1,nsym2) .eq. 1)  then

         cnt=cnt+1
         term=0.0
         term1=0.0
         term2=0.0

         eka=e(k)-e(a)
         ebprcprlprmpr=e(bpr)+e(cpr)-e(lpr)-e(mpr)

         term=term+(+1.0*vpqrs(bpr,lpr,cpr,mpr))*t2_3h3p_hc(a,bpr,cpr,k,lpr,mpr,1)
         term=term+(-1.0*vpqrs(bpr,mpr,cpr,lpr))*t2_3h3p_hc(a,bpr,cpr,k,lpr,mpr,2)
         term1=term/eka 
         term2=term/ebprcprlprmpr

         density=density+(term1+term2)*(0.25)

         end if
       end do
     end do
   end do
 end do
!!! end t2_3h3p_hc_part !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 end if

 end function density




 real(d) function densityhc(k,a)
 
 integer, intent(in) :: k,a

!!! I only have to write the occupied-unoccupied part of the density matrix !!!
!!!!!!!   THE APR OR KPR HAVE ONLY ONE POSSIBLITY FOR THE SPIN
!!!!!!!   THE FACTOR OF TWO IS TAKEN INTO ACCOUNT IN THE CALLING
!!!!!!!   BY THE FUNCTIONS D2_1,D2_2,D2_3,D2_4.
!!!!!!!   THE FACTOR OF 2 TAKES INTO ACCOUHT THAT OR APR OR KPR
!!!!!!!   CAN TAKE ONE TIME ALPHA AND ONE TIME BETA VALUES AND
!!!!!!!   THE EXPRESSION DOES NOT DEPEND OF WHICH VALUE IT TAKES

    integer :: b,c,j,i,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*8 :: DA,DB,eabij,ebckj,eak,term

!!!!!!!!!!!!!!!!!!!!!!!!!! ZEROTH ORDER PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! IT DOES NOT GIVE ANY CONTRIBUTION TO THE HOLE-PARTICLE PART !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! SECOND ORDER PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    densityhc=0.0
    cnt=0
  do b1=nOcc+1,nBas
     b=roccnum(b1)
   do c1=nOcc+1,nBas
      c=roccnum(c1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(b),orbSym(c))
             nsym2=MT(orbSym(a),orbSym(j))
             nsym3=MT(orbSym(k),orbSym(j))
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1
                term=0.0

                eak=e(a)-e(k)
                ebckj=e(b)+e(c)-e(k)-e(j)
                DA=eak*ebckj

                term=term+vpqrs(a,b,j,c)*(+2.0*vpqrs(b,k,c,j)-1.0*vpqrs(c,k,b,j))
                term=term+vpqrs(j,b,a,c)*(-1.0*vpqrs(b,k,c,j)+2.0*vpqrs(c,k,b,j))
                term=term/DA
                
                densityhc=densityhc+term
                
             end if
          end do
       end do
    end do 

do b1=nOcc+1,nBas
     b=roccnum(b1)
   do i1=1,nOcc
      i=roccnum(i1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(i),orbSym(j))
             nsym2=MT(orbSym(b),orbSym(a))
             nsym3=MT(orbSym(b),orbSym(k))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eak=e(a)-e(k)
                eabij=e(a)+e(b)-e(i)-e(j)
                DB=eak*eabij

                term=term+vpqrs(i,k,j,b)*(+2.0*vpqrs(a,i,b,j)-1.0*vpqrs(b,i,a,j))
                term=term+vpqrs(j,k,i,b)*(-1.0*vpqrs(a,i,b,j)+2.0*vpqrs(b,i,a,j))
                term=term/DB

                densityhc=densityhc-term  ! this part has minus sign

             end if
          end do
       end do
    end do

   densityhc = 0.5*densityhc  ! EXPRESSION FACTOR IN THE DENSITY

 end function densityhc


 real(d) function proper_density(k,a)
 
 integer, intent(in) :: k,a

!!! I only have to write the occupied-unoccupied part of the density matrix !!!
!!!!!!!   THE APR OR KPR HAVE ONLY ONE POSSIBLITY FOR THE SPIN
!!!!!!!   THE FACTOR OF TWO IS TAKEN INTO ACCOUNT IN THE CALLING
!!!!!!!   BY THE FUNCTIONS D2_1,D2_2,D2_3,D2_4.
!!!!!!!   THE FACTOR OF 2 TAKES INTO ACCOUHT THAT OR APR OR KPR
!!!!!!!   CAN TAKE ONE TIME ALPHA AND ONE TIME BETA VALUES AND
!!!!!!!   THE EXPRESSION DOES NOT DEPEND OF WHICH VALUE IT TAKES

    integer :: b,c,j,i,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*16 :: DA,DB,eabij,ebckj,eak,term

!!!!!!!!!!!!!!!!!!!!!!!!!! ZEROTH ORDER PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! IT DOES NOT GIVE ANY CONTRIBUTION TO THE HOLE-PARTICLE PART !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! SECOND ORDER PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    proper_density=0.0
    cnt=0
  do b1=nOcc+1,nBas
     b=roccnum(b1)
   do c1=nOcc+1,nBas
      c=roccnum(c1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(b),orbSym(c))
             nsym2=MT(orbSym(a),orbSym(j))
             nsym3=MT(orbSym(k),orbSym(j))
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1
                term=0.0

                eak=e(a)-e(k)
                ebckj=e(b)+e(c)-e(k)-e(j)
                DA=eak*ebckj

                term=term+vpqrs(a,b,j,c)*(+2.0*vpqrs(b,k,c,j)-1.0*vpqrs(c,k,b,j))
                term=term+vpqrs(j,b,a,c)*(-1.0*vpqrs(b,k,c,j)+2.0*vpqrs(c,k,b,j))
                term=term/DA
                
                proper_density = proper_density + term
                
             end if
          end do
       end do
    end do 

do b1=nOcc+1,nBas
     b=roccnum(b1)
   do i1=1,nOcc
      i=roccnum(i1)
     do j1=1,nOcc
        j=roccnum(j1)

             nsym1=MT(orbSym(i),orbSym(j))
             nsym2=MT(orbSym(b),orbSym(a))
             nsym3=MT(orbSym(b),orbSym(k))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then

                cnt=cnt+1
                term=0.0

                eak=e(a)-e(k)
                eabij=e(a)+e(b)-e(i)-e(j)
                DB=eak*eabij

                term=term+vpqrs(i,k,j,b)*(+2.0*vpqrs(a,i,b,j)-1.0*vpqrs(b,i,a,j))
                term=term+vpqrs(j,k,i,b)*(-1.0*vpqrs(a,i,b,j)+2.0*vpqrs(b,i,a,j))
                term=term/DB

                proper_density = proper_density - term  ! this part has minus sign

             end if
          end do
       end do
    end do

   proper_density = 0.5*proper_density  ! EXPRESSION FACTOR IN THE DENSITY

 end function proper_density




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! MATRIX BLOCKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************PH-PH BLOCK***********************************
!!$*******************************************************************************
!!$*******************************************************************************

!!$Indices are supplied in the order: PH,PH
!!$********************ZEROTH-ORDER******************ZEROTH-ORDER**************************ZEROTH-ORDER

!!$Zeroth order contribution D0_1_ak,a'k'. The condition that k=k'
!!$is checked by the calling procedure

  real(d) function D0_1_ph_ph(a,apr)
    
    integer, intent(in) :: a,apr

    D0_1_ph_ph = 0.0
 
    D0_1_ph_ph = D0_1_ph_ph + 2._d*dpl(a,apr)
    
    D0_1_ph_ph = D0_1_ph_ph/2._d

  end function D0_1_ph_ph

!!$Zeroth order contribution D0_2_ak,a'k'. The condition that a=a'
!!$is checked by the calling procedure


 real(d) function D0_2_ph_ph(k,kpr)

    integer, intent(in) :: k,kpr

    D0_2_ph_ph = 0.0
 
    D0_2_ph_ph = D0_2_ph_ph + 2._d*dpl(kpr,k)

    D0_2_ph_ph = D0_2_ph_ph/2._d

    D0_2_ph_ph = -D0_2_ph_ph   ! EXPRESSION FACTOR IN THE PAPER

  end function D0_2_ph_ph

!!$********************SECOND-ORDER******************SECOND-ORDER**************************SECOND-ORDER
!!$********************SECOND-ORDER******************SECOND-ORDER**************************SECOND-ORDER
!!$********************SECOND-ORDER******************SECOND-ORDER**************************SECOND-ORDER
!!$********************SECOND-ORDER******************SECOND-ORDER**************************SECOND-ORDER
!!$********************SECOND-ORDER******************SECOND-ORDER**************************SECOND-ORDER
!!$********************SECOND-ORDER******************SECOND-ORDER**************************SECOND-ORDER
!!$********************SECOND-ORDER******************SECOND-ORDER**************************SECOND-ORDER




!!!!!!!  SECOND OREDR PART TERMS INVOLVING THE DENSITY MATRIX HOLE-PARTICLE ELEMENTS
!!!!!!!  SECOND OREDR PART TERMS INVOLVING THE DENSITY MATRIX HOLE-PARTICLE ELEMENTS
!!!!!!!  SECOND OREDR PART TERMS INVOLVING THE DENSITY MATRIX HOLE-PARTICLE ELEMENTS
!!!!!!!  SECOND OREDR PART TERMS INVOLVING THE DENSITY MATRIX HOLE-PARTICLE ELEMENTS
!!!!!!!  SECOND OREDR PART TERMS INVOLVING THE DENSITY MATRIX HOLE-PARTICLE ELEMENTS
!!!!!!!  SECOND OREDR PART TERMS INVOLVING THE DENSITY MATRIX HOLE-PARTICLE ELEMENTS
!!!!!!!  SECOND OREDR PART TERMS INVOLVING THE DENSITY MATRIX HOLE-PARTICLE ELEMENTS
!!!!!!!  SECOND OREDR PART TERMS INVOLVING THE DENSITY MATRIX HOLE-PARTICLE ELEMENTS
!!!!!!!  SECOND OREDR PART TERMS INVOLVING THE DENSITY MATRIX HOLE-PARTICLE ELEMENTS



!!$SECOND order contribution D2_1_ak,a'k'. The condition that k=k' is checked by the calling procedure

 real(d) function D2_1_ph_ph(a,apr)

    integer :: l1,l,sym
    integer, intent(in) :: a,apr

    D2_1_ph_ph = 0.0
    
    do l1=1,nOcc
       l=roccnum(l1)
       
       sym=MT(orbSym(a),orbSym(l))
       if (sym .eq. CHECK_dip) then
          
          D2_1_ph_ph = D2_1_ph_ph  + 2.0*density(l,apr)*dpl(a,l)
          
       end if
       
    end do
    
    D2_1_ph_ph = D2_1_ph_ph/2.0   !!! NORMALIZATION FACTOR
    D2_1_ph_ph = -D2_1_ph_ph   ! EXPRESSION FACTOR IN THE PAPER
    
  end function D2_1_ph_ph

!!$SECOND order contribution D2_2_ak,a'k'. The condition that k'=k
!!$is checked by the calling procedure

 real(d) function D2_2_ph_ph(a,apr)

    integer :: l1,l,sym
    integer, intent(in) :: a,apr


    D2_2_ph_ph = 0.0
 

   do l1=1,nOcc
       l=roccnum(l1)

     sym=MT(orbSym(apr),orbSym(l))
     if (sym .eq. CHECK_dip) then

     D2_2_ph_ph = D2_2_ph_ph + 2.0*density(l,a)*dpl(apr,l)
!     D2_2_ph_ph=+2.0*density_matrix(l,a)*dpl(apr,l)

!     write(6,*) 'd2',l,a,densityhc(l,a)
     end if

   end do

    D2_2_ph_ph = D2_2_ph_ph/2.0   !!! NORMALIZATION FACTOR
    D2_2_ph_ph = -D2_2_ph_ph   ! EXPRESSION FACTOR IN THE PAPER

  end function D2_2_ph_ph


!!$SECOND order contribution D2_3_ak,a'k'. The condition that a=a'
!!$is checked by the calling procedure

 real(d) function D2_3_ph_ph(k,kpr)

    integer :: b1,b,sym
    integer, intent(in) :: k,kpr


    D2_3_ph_ph = 0.0
 
   do b1=nOcc+1,nBas
       b=roccnum(b1)

     sym=MT(orbSym(b),orbSym(k))
     if (sym .eq. CHECK_dip) then

!     D2_3_ph_ph=+2.0*density_matrix(kpr,b)*dpl(b,k)
     D2_3_ph_ph = D2_3_ph_ph + 2.0*density(kpr,b)*dpl(b,k)


!     write(6,*) 'd3',kpr,b,density(kpr,b)

     end if

   end do

    D2_3_ph_ph= D2_3_ph_ph/2.0   !!! NORMALIZATION FACTOR
    D2_3_ph_ph= -D2_3_ph_ph   ! EXPRESSION FACTOR IN THE PAPER

  end function D2_3_ph_ph


!!$SECOND order contribution D2_4_ak,a'k'. The condition that a'=a
!!$is checked by the calling procedure

 real(d) function D2_4_ph_ph(k,kpr)

    integer :: b1,b,sym
    integer, intent(in) :: k,kpr

    D2_4_ph_ph = 0.0

   do b1=nOcc+1,nBas
       b=roccnum(b1)

     sym=MT(orbSym(b),orbSym(kpr))
     if (sym .eq. CHECK_dip) then

     D2_4_ph_ph = D2_4_ph_ph + 2.0*density(k,b)*dpl(b,kpr)
!     D2_4_ph_ph=+2.0*density_matrix(k,b)*dpl(b,kpr)

!     write(6,*) 'd4',k,b,densityhc(k,b)

     end if

   end do

    D2_4_ph_ph = D2_4_ph_ph/2.0   !!! NORMALIZATION FACTOR
    D2_4_ph_ph = -D2_4_ph_ph   ! EXPRESSION FACTOR IN THE PAPER

  end function D2_4_ph_ph


!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
!!!!    SECOND ORDER PART  NOT INCLUDING  THE DENSITY MATRIX TERMS
 

 !!$Second order contribution D2_2_1_ak,a'k'.The condition that k=k' is
!!$checked in the calling procedure.


    real(d) function D2_2_1_ph_ph(a,apr)

    integer, intent(in) :: a,apr

    integer :: b,c,i,j,sym,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_2_1_ph_ph = 0.0
    cnt=0
  do b1=nOcc+1,nBas
     b=roccnum(b1)
     sym=MT(orbSym(b),orbSym(apr))
     if (sym .eq. CHECK_dip) then
   do c1=nOcc+1,nBas
      c=roccnum(c1)
    do i1=1,nOcc
        i=roccnum(i1)
       do j1=1,nOcc
          j=roccnum(j1)

             nsym1=MT(orbSym(i),orbSym(j))
             nsym2=MT(orbSym(b),orbSym(c))
             nsym3=MT(orbSym(a),orbSym(c))
             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(a)+e(c)-e(i)-e(j)
                ebcij=e(b)+e(c)-e(i)-e(j)
                DA=eacij*ebcij

                term = term + vpqrs(i,b,j,c)*(+4.0*vpqrs(a,i,c,j)-2.0*vpqrs(a,j,c,i))
                term = term + vpqrs(i,c,j,b)*(-2.0*vpqrs(a,i,c,j)+4.0*vpqrs(a,j,c,i))
                term = term*dpl(b,apr)
                term = term/DA
                
                D2_2_1_ph_ph = D2_2_1_ph_ph + term
                
             end if
          end do
       end do
    end do 
   end if
   end do
 
    D2_2_1_ph_ph = + 0.5*D2_2_1_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             
    D2_2_1_ph_ph = -0.25*D2_2_1_ph_ph  ! EXPRESSION FACTOR


  end function D2_2_1_ph_ph


 !!$Second order contribution D2_2_2_ak,a'k'. (hc)The condition that k=k' is
!!$checked in the calling procedure.


    real(d) function D2_2_2_ph_ph(a,apr)

    integer, intent(in) :: a,apr

    integer :: b,c,i,j,sym,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_2_2_ph_ph = 0.0
    cnt=0
  do b1=nOcc+1,nBas
     b=roccnum(b1)
     sym=MT(orbSym(b),orbSym(a))
     if (sym .eq. CHECK_dip) then
 

   do c1=nOcc+1,nBas
     c=roccnum(c1)
    do i1=1,nOcc
       i=roccnum(i1)
       do j1=1,nOcc
          j=roccnum(j1)
    
             nsym1=MT(orbSym(i),orbSym(j))
             nsym2=MT(orbSym(b),orbSym(c))
             nsym3=MT(orbSym(apr),orbSym(c))
             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(apr)+e(c)-e(i)-e(j)
                ebcij=e(b)+e(c)-e(i)-e(j)
                DA=eacij*ebcij

                term=term+vpqrs(b,i,c,j)*(+4.0*vpqrs(i,apr,j,c)-2.0*vpqrs(i,c,j,apr))
                term=term+vpqrs(b,j,c,i)*(-2.0*vpqrs(i,apr,j,c)+4.0*vpqrs(i,c,j,apr))
                term=term*dpl(b,a)
                term=term/DA
                
                D2_2_2_ph_ph=D2_2_2_ph_ph+term
                
             end if
          end do
       end do
    end do 
    end if
  end do
    
    D2_2_2_ph_ph=+0.5*D2_2_2_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             
    D2_2_2_ph_ph=-0.25*D2_2_2_ph_ph  ! EXPRESSION FACTOR


  end function D2_2_2_ph_ph

 !!$Second order contribution D2_3_1_ak,a'k'.The condition that k=k' is
!!$checked in the calling procedure.


    real(d) function D2_3_1_ph_ph(a,apr)

    integer, intent(in) :: a,apr

    integer :: b,c,i,j,sym,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_3_1_ph_ph=0.0
    cnt=0
  do b1=nOcc+1,nBas
     b=roccnum(b1)
   do c1=nOcc+1,nBas
      c=roccnum(c1)
     sym=MT(orbSym(b),orbSym(c))
     if (sym .eq. CHECK_dip) then
 

    do i1=1,nOcc
       i=roccnum(i1)
       do j1=1,nOcc
          j=roccnum(j1)

             nsym1=MT(orbSym(i),orbSym(j))
             nsym2=MT(orbSym(apr),orbSym(b))
             nsym3=MT(orbSym(a),orbSym(c))
             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(apr)+e(b)-e(i)-e(j)
                ebcij=e(a)+e(c)-e(i)-e(j)
                DA=eacij*ebcij

                term=term+vpqrs(i,apr,j,b)*(+4.0*vpqrs(a,i,c,j)-2.0*vpqrs(a,j,c,i))
                term=term+vpqrs(i,b,j,apr)*(-2.0*vpqrs(a,i,c,j)+4.0*vpqrs(a,j,c,i))
                term=term*dpl(b,c)
                term=term/DA
                
                D2_3_1_ph_ph=D2_3_1_ph_ph+term
                
             end if
          end do
       end do
     end if
   end do 
  end do  

    D2_3_1_ph_ph=+0.5*D2_3_1_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             
    D2_3_1_ph_ph=-0.5*D2_3_1_ph_ph  ! EXPRESSION FACTOR


  end function D2_3_1_ph_ph


 !!$Second order contribution D2_3_2_ak,a'k'.(hc)The condition that k=k' is
!!$checked in the calling procedure.


    real(d) function D2_3_2_ph_ph(a,apr)

    integer, intent(in) :: a,apr

    integer :: b,i,j,l,sym,nsym1,nsym2,nsym3,nsym4,b1,i1,j1,l1,cnt
    real*8   :: DA,eacij,ebcij,term

    D2_3_2_ph_ph=0.0
    cnt=0
 do l1=1,nOcc
      l=roccnum(l1)
  do i1=1,nOcc
       i=roccnum(i1)
     sym=MT(orbSym(l),orbSym(i))
     if (sym .eq. CHECK_dip) then
    do b1=nOcc+1,nBas
       b=roccnum(b1)
      do j1=1,nOcc
          j=roccnum(j1)

             nsym1=MT(orbSym(j),orbSym(l))
             nsym2=MT(orbSym(a),orbSym(b))
             nsym3=MT(orbSym(apr),orbSym(b))
             nsym4=MT(orbSym(i),orbSym(j))             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(apr)+e(b)-e(i)-e(j)
                ebcij=e(a)+e(b)-e(l)-e(j)
                DA=eacij*ebcij

                term=term+vpqrs(i,apr,j,b)*(+4.0*vpqrs(a,l,b,j)-2.0*vpqrs(a,j,b,l))
                term=term+vpqrs(i,b,j,apr)*(-2.0*vpqrs(a,l,b,j)+4.0*vpqrs(a,j,b,l))
                term=term*dpl(l,i)
                term=term/DA
                
                D2_3_2_ph_ph=D2_3_2_ph_ph+term
                
             end if
           end do
        end do
      end if
    end do 
 end do
  
    D2_3_2_ph_ph=+0.5*D2_3_2_ph_ph   ! NORMALIZATION OF WAVE FUNCTION
             
! EXPRESSION FACTOR +1


  end function D2_3_2_ph_ph

!!$Second order contribution D2_4_1_ak,a'k'. The condition that a=a' is
!!$checked in the calling procedure.

  real(d) function D2_4_1_ph_ph(k,kpr)

    integer, intent(in) :: k,kpr

    integer :: b,c,i,j,sym,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_4_1_ph_ph=0.0
    cnt=0
do j1=1,nOcc
   j=roccnum(j1)
   sym=MT(orbSym(j),orbSym(kpr))
     if (sym .eq. CHECK_dip) then
  do b1=nOcc+1,nBas
     b=roccnum(b1)
    do c1=nOcc+1,nBas
       c=roccnum(c1)
       do i1=1,nOcc
          i=roccnum(i1)
 
             nsym1=MT(orbSym(b),orbSym(c))
             nsym2=MT(orbSym(i),orbSym(k))
             nsym3=MT(orbSym(i),orbSym(j))
                         
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(b)+e(c)-e(j)-e(i)
                ebcij=e(b)+e(c)-e(k)-e(i)
                DA=eacij*ebcij

                term=term+vpqrs(j,b,i,c)*(+4.0*vpqrs(b,k,c,i)-2.0*vpqrs(b,i,c,k))
                term=term+vpqrs(j,c,i,b)*(-2.0*vpqrs(b,k,c,i)+4.0*vpqrs(b,i,c,k))
                term=term*dpl(kpr,j)
                term=term/DA
                
                D2_4_1_ph_ph=D2_4_1_ph_ph+term
                
             end if
     
          end do
       end do
    end do 
  end if
end do
  
    D2_4_1_ph_ph=+0.5*D2_4_1_ph_ph   ! NORMALIZATION OF WAVE FUNCTION  
    D2_4_1_ph_ph=+0.25*D2_4_1_ph_ph  ! EXPRESSION FACTOR            


  end function D2_4_1_ph_ph

 !!$Second order contribution D2_4_2_ak,a'k'.(hc) The condition that a=a' is
!!$checked in the calling procedure.


  real(d) function D2_4_2_ph_ph(k,kpr)

    integer, intent(in) :: k,kpr

    integer :: b,c,i,j,sym,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_4_2_ph_ph=0.0
    cnt=0
do j1=1,nOcc
   j=roccnum(j1)
   sym=MT(orbSym(j),orbSym(k))
     if (sym .eq. CHECK_dip) then
  do b1=nOcc+1,nBas
     b=roccnum(b1)
    do c1=nOcc+1,nBas
       c=roccnum(c1)
       do i1=1,nOcc
          i=roccnum(i1)
 
             nsym1=MT(orbSym(b),orbSym(c))
             nsym2=MT(orbSym(i),orbSym(j))
             nsym3=MT(orbSym(i),orbSym(kpr))
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(b)+e(c)-e(j)-e(i)
                ebcij=e(b)+e(c)-e(kpr)-e(i)
                DA=eacij*ebcij

                term=term+vpqrs(b,j,c,i)*(+4.0*vpqrs(kpr,b,i,c)-2.0*vpqrs(kpr,c,i,b))
                term=term+vpqrs(b,i,c,j)*(-2.0*vpqrs(kpr,b,i,c)+4.0*vpqrs(kpr,c,i,b))
                term=term*dpl(k,j)
                term=term/DA
                
                D2_4_2_ph_ph=D2_4_2_ph_ph+term
                
             end if
     
          end do
       end do
    end do 
  end if
end do
 
    D2_4_2_ph_ph=+0.5*D2_4_2_ph_ph   ! NORMALIZATION OF WAVE FUNCTION  
    D2_4_2_ph_ph=+0.25*D2_4_2_ph_ph  ! EXPRESSION FACTOR            


  end function D2_4_2_ph_ph
 

!!$Second order contribution D2_5_1_ak,a'k'. The condition that a=a' is
!!$checked in the calling procedure.

  real(d) function D2_5_1_ph_ph(k,kpr)

    integer, intent(in) :: k,kpr

    integer :: b,c,d,i,sym,nsym1,nsym2,nsym3,nsym4,b1,c1,d1,i1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_5_1_ph_ph=0.0
    cnt=0
  do b1=nOcc+1,nBas
     b=roccnum(b1)
    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(b),orbSym(c))
     if (sym .eq. CHECK_dip) then
 
       do d1=nOcc+1,nBas
          d=roccnum(d1)
          do i1=1,nOcc
             i=roccnum(i1)

             nsym1=MT(orbSym(b),orbSym(d))
             nsym2=MT(orbSym(i),orbSym(k))
             nsym3=MT(orbSym(c),orbSym(d))
             nsym4=MT(orbSym(i),orbSym(kpr))
                          
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(b)+e(d)-e(k)-e(i)
                ebcij=e(c)+e(d)-e(kpr)-e(i)
                DA=eacij*ebcij

                term=term+vpqrs(kpr,c,i,d)*(+4.0*vpqrs(b,k,d,i)-2.0*vpqrs(b,i,d,k))
                term=term+vpqrs(kpr,d,i,c)*(-2.0*vpqrs(b,k,d,i)+4.0*vpqrs(b,i,d,k))
                term=term*dpl(c,b)
                term=term/DA
                
                D2_5_1_ph_ph=D2_5_1_ph_ph+term
                
             end if
         end do
       end do
     end if
   end do 
  end do
    
    D2_5_1_ph_ph=+0.5*D2_5_1_ph_ph   ! NORMALIZATION OF WAVE FUNCTION
    D2_5_1_ph_ph=-1.0*D2_5_1_ph_ph    ! EXPRESSION FACTOR            


  end function D2_5_1_ph_ph


!!$Second order contribution D2_5_2_ak,a'k'. The condition that a=a' is
!!$checked in the calling procedure.

  real(d) function D2_5_2_ph_ph(k,kpr)

    integer, intent(in) :: k,kpr

    integer :: b,c,i,j,sym,nsym1,nsym2,nsym3,b1,c1,i1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_5_2_ph_ph=0.0
    cnt=0
do i1=1,nOcc
   i=roccnum(i1)
  do j1=1,nOcc
     j=roccnum(j1)
     sym=MT(orbSym(j),orbSym(i))
     if (sym .eq. CHECK_dip) then
   do b1=nOcc+1,nBas
      b=roccnum(b1)
     do c1=nOcc+1,nBas
        c=roccnum(c1)

             nsym1=MT(orbSym(b),orbSym(c))
             nsym2=MT(orbSym(k),orbSym(j))
             nsym3=MT(orbSym(kpr),orbSym(i))             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(b)+e(c)-e(k)-e(j)
                ebcij=e(b)+e(c)-e(kpr)-e(i)
                DA=eacij*ebcij

                term=term+vpqrs(kpr,b,i,c)*(+4.0*vpqrs(b,k,c,j)-2.0*vpqrs(b,j,c,k))
                term=term+vpqrs(kpr,c,i,b)*(-2.0*vpqrs(b,k,c,j)+4.0*vpqrs(b,j,c,k))
                term=term*dpl(j,i)
                term=term/DA
                
                D2_5_2_ph_ph=D2_5_2_ph_ph+term
                
             end if
       
          end do
       end do
     end if
    end do 
  end do
 
    D2_5_2_ph_ph=+0.5*D2_5_2_ph_ph   ! NORMALIZATION OF WAVE FUNCTION
    D2_5_2_ph_ph=+0.5*D2_5_2_ph_ph    ! EXPRESSION FACTOR            


  end function D2_5_2_ph_ph
 
 !!$Second order contribution D2_6_1_ak,a'k'. NO DELTA FUNCTIONS(hc)

    real(d) function D2_6_1_ph_ph(a,apr,k,kpr)

    integer, intent(in) :: a,apr,k,kpr

    integer :: b,c,i,sym,nsym1,nsym2,nsym3,nsym4,b1,c1,i1,cnt
    real*8 :: DA,eacij,ebcij,term,ftmp

    D2_6_1_ph_ph=0.0
    cnt=0
  do b1=nOcc+1,nBas
     b=roccnum(b1)
     sym=MT(orbSym(b),orbSym(apr))
     if (sym .eq. CHECK_dip) then
 
    do c1=nOcc+1,nBas
       c=roccnum(c1)
       do i1=1,nOcc
          i=roccnum(i1)

             nsym1=MT(orbSym(a),orbSym(c))
             nsym2=MT(orbSym(i),orbSym(k))
             nsym3=MT(orbSym(b),orbSym(c))
             nsym4=MT(orbSym(i),orbSym(kpr))             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(a)+e(c)-e(k)-e(i)
                ebcij=e(b)+e(c)-e(kpr)-e(i)
                DA=eacij*ebcij

                term=term+vpqrs(kpr,b,i,c)*(+8.0*vpqrs(a,k,c,i)-4.0*vpqrs(a,i,c,k))
                term=term+vpqrs(kpr,c,i,b)*(-4.0*vpqrs(a,k,c,i)+2.0*vpqrs(a,i,c,k))
                term=term*dpl(b,apr)
                term=term/DA
                
                D2_6_1_ph_ph=D2_6_1_ph_ph+term
                
             end if
          end do
   end do
        end if
    end do 

    D2_6_1_ph_ph=+0.5*D2_6_1_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             
    D2_6_1_ph_ph=+0.5*D2_6_1_ph_ph   ! EXPRESSION FACTOR

  end function D2_6_1_ph_ph


  !!$Second order contribution D2_6_2_ak,a'k'. NO DELTA FUNCTIONS (hc)

    real(d) function D2_6_2_ph_ph(a,apr,k,kpr)

    integer, intent(in) :: a,apr,k,kpr

    integer :: b,c,i,sym,nsym1,nsym2,nsym3,nsym4,b1,c1,i1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_6_2_ph_ph=0.0
    cnt=0
    do b1=nOcc+1,nBas
       b=roccnum(b1)
       sym=MT(orbSym(b),orbSym(a))
       if (sym .eq. CHECK_dip) then
          
          do c1=nOcc+1,nBas
             c=roccnum(c1)
             do i1=1,nOcc
                i=roccnum(i1)
                
                nsym1=MT(orbSym(apr),orbSym(c))
                nsym2=MT(orbSym(i),orbSym(kpr))
                nsym3=MT(orbSym(b),orbSym(c))
                nsym4=MT(orbSym(i),orbSym(k))
                
                if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                   
                   cnt=cnt+1
                   
                   term=0.0

                   eacij=e(b)+e(c)-e(k)-e(i)
                   ebcij=e(apr)+e(c)-e(kpr)-e(i)
                   DA=eacij*ebcij

                   term=term+vpqrs(b,k,c,i)*(+8.0*vpqrs(kpr,apr,i,c)-4.0*vpqrs(kpr,c,i,apr))
                   term=term+vpqrs(b,i,c,k)*(-4.0*vpqrs(kpr,apr,i,c)+2.0*vpqrs(kpr,c,i,apr))
                   term=term*dpl(b,a)
                   term=term/DA
                   
                   D2_6_2_ph_ph=D2_6_2_ph_ph+term
                
                end if
             end do
          end do
       end if
    end do
    
    D2_6_2_ph_ph=+0.5*D2_6_2_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             
    D2_6_2_ph_ph=+0.5*D2_6_2_ph_ph   ! EXPRESSION FACTOR
    

  end function D2_6_2_ph_ph


!!$Second order contribution D2_6_3_ak,a'k'. 

  real(d) function D2_6_3_ph_ph(a,apr,k,kpr)

    integer, intent(in) :: a,apr,k,kpr

    integer :: b,i,j,sym,nsym1,nsym2,nsym3,nsym4,b1,i1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_6_3_ph_ph=0.0
    cnt=0
do j1=1,nOcc
   j=roccnum(j1)
   sym=MT(orbSym(j),orbSym(kpr))
     if (sym .eq. CHECK_dip) then
  do b1=nOcc+1,nBas
     b=roccnum(b1)
       do i1=1,nOcc
          i=roccnum(i1)

             nsym1=MT(orbSym(a),orbSym(b))
             nsym2=MT(orbSym(i),orbSym(k))
             nsym3=MT(orbSym(apr),orbSym(b))
             nsym4=MT(orbSym(i),orbSym(j))
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(apr)+e(b)-e(j)-e(i)
                ebcij=e(a)+e(b)-e(k)-e(i)
                DA=eacij*ebcij

                term=term+vpqrs(j,apr,i,b)*(+8.0*vpqrs(a,k,b,i)-4.0*vpqrs(a,i,b,k))
                term=term+vpqrs(j,b,i,apr)*(-4.0*vpqrs(a,k,b,i)+2.0*vpqrs(a,i,b,k))
                term=term*dpl(kpr,j)
                term=term/DA
                
                D2_6_3_ph_ph=D2_6_3_ph_ph+term
                
             end if
       
          end do
       end do
     end if
    end do 
 
    D2_6_3_ph_ph=+0.5*D2_6_3_ph_ph   ! NORMALIZATION OF WAVE FUNCTION 
    D2_6_3_ph_ph=-0.5*D2_6_3_ph_ph   ! EXPRESSION FACTOR            


  end function D2_6_3_ph_ph


!!$Second order contribution D2_6_4_ak,a'k'. (hc)

  real(d) function D2_6_4_ph_ph(a,apr,k,kpr)

    integer, intent(in) :: a,apr,k,kpr

    integer :: b,i,j,sym,nsym1,nsym2,nsym3,nsym4,b1,i1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_6_4_ph_ph=0.0
    cnt=0
do j1=1,nOcc
   j=roccnum(j1)
   sym=MT(orbSym(j),orbSym(k))
     if (sym .eq. CHECK_dip) then
  do b1=nOcc+1,nBas
     b=roccnum(b1)
       do i1=1,nOcc
          i=roccnum(i1)
 
             nsym1=MT(orbSym(apr),orbSym(b))
             nsym2=MT(orbSym(i),orbSym(kpr))
             nsym3=MT(orbSym(a),orbSym(b))
             nsym4=MT(orbSym(i),orbSym(j))
             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(apr)+e(b)-e(kpr)-e(i)
                ebcij=e(a)+e(b)-e(j)-e(i)
                DA=eacij*ebcij

                term=term+vpqrs(a,j,b,i)*(+8.0*vpqrs(kpr,apr,i,b)-4.0*vpqrs(kpr,b,i,apr))
                term=term+vpqrs(a,i,b,j)*(-4.0*vpqrs(kpr,apr,i,b)+2.0*vpqrs(kpr,b,i,apr))
                term=term*dpl(k,j)
                term=term/DA
                
                D2_6_4_ph_ph=D2_6_4_ph_ph+term
                
             end if
       
          end do
       end do
     end if
    end do 
 
    D2_6_4_ph_ph=+0.5*D2_6_4_ph_ph   ! NORMALIZATION OF WAVE FUNCTION 
    D2_6_4_ph_ph=-0.5*D2_6_4_ph_ph   ! EXPRESSION FACTOR            


  end function D2_6_4_ph_ph
 
 
!!$Second order contribution D2_7_1_ak,a'k'. 

  real(d) function D2_7_1_ph_ph(a,apr,k,kpr)

    integer, intent(in) :: a,apr,k,kpr

    integer :: b,i,j,sym,nsym1,nsym2,nsym3,nsym4,b1,i1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_7_1_ph_ph=0.0
    cnt=0
  do i1=1,nOcc
     i=roccnum(i1)
      do j1=1,nOcc
         j=roccnum(j1)
         sym=MT(orbSym(j),orbSym(i))
     if (sym .eq. CHECK_dip) then
            do b1=nOcc+1,nBas
               b=roccnum(b1)
 
             nsym1=MT(orbSym(a),orbSym(b))
             nsym2=MT(orbSym(k),orbSym(j))
             nsym3=MT(orbSym(apr),orbSym(b))
             nsym4=MT(orbSym(kpr),orbSym(i))
             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(a)+e(b)-e(k)-e(j)
                ebcij=e(apr)+e(b)-e(kpr)-e(i)
                DA=eacij*ebcij

                term=term+vpqrs(kpr,apr,i,b)*(+8.0*vpqrs(a,k,b,j)-4.0*vpqrs(a,j,b,k))
                term=term+vpqrs(kpr,b,i,apr)*(-4.0*vpqrs(a,k,b,j)+2.0*vpqrs(a,j,b,k))
                term=term*dpl(j,i)
                term=term/DA
                
                D2_7_1_ph_ph=D2_7_1_ph_ph+term
                
             end if
          end do
        end if
       end do
    end do 
 
    D2_7_1_ph_ph=+0.5*D2_7_1_ph_ph   ! NORMALIZATION OF WAVE FUNCTION 
    D2_7_1_ph_ph=-1.0*D2_7_1_ph_ph    ! EXPRESSION FACTOR            


  end function D2_7_1_ph_ph
 
!!$Second order contribution D2_7_2_ak,a'k'. 

  real(d) function D2_7_2_ph_ph(a,apr,k,kpr)

    integer, intent(in) :: a,apr,k,kpr

    integer :: b,c,j,sym,nsym1,nsym2,nsym3,nsym4,b1,c1,j1,cnt
    real*8 :: DA,eacij,ebcij,term

    D2_7_2_ph_ph=0.0
    cnt=0
  do b1=nOcc+1,nBas
     b=roccnum(b1)
       do c1=nOcc+1,nBas
          c=roccnum(c1)
     sym=MT(orbSym(b),orbSym(c))
     if (sym .eq. CHECK_dip) then
 
          do j1=1,nOcc
             j=roccnum(j1)

             nsym1=MT(orbSym(a),orbSym(c))
             nsym2=MT(orbSym(j),orbSym(k))
             nsym3=MT(orbSym(apr),orbSym(b))
             nsym4=MT(orbSym(j),orbSym(kpr))             
             
             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
                cnt=cnt+1

                term=0.0

                eacij=e(a)+e(c)-e(k)-e(j)
                ebcij=e(apr)+e(b)-e(kpr)-e(j)
                DA=eacij*ebcij

                term=term+vpqrs(kpr,apr,j,b)*(+8.0*vpqrs(a,k,c,j)-4.0*vpqrs(a,j,c,k))
                term=term+vpqrs(kpr,b,j,apr)*(-4.0*vpqrs(a,k,c,j)+2.0*vpqrs(a,j,c,k))
                term=term*dpl(b,c)
                term=term/DA
                
                D2_7_2_ph_ph=D2_7_2_ph_ph+term
                
             end if
          end do
       end if
       end do
 end do 
    
    D2_7_2_ph_ph=+0.5*D2_7_2_ph_ph   ! NORMALIZATION OF WAVE FUNCTION 
    D2_7_2_ph_ph=+1.0*D2_7_2_ph_ph    ! EXPRESSION FACTOR            


  end function D2_7_2_ph_ph

  
!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************PH-2P2H BLOCK**********************************
!!$*******************************************************************************
!!$*******************************************************************************


!!$  We distinguish here between five different types of coupling. 
!!$ Calculating Dak,a'b'k'l'
  
!!$ a'=b' and k'=l'; EMPTY AND FILL THE SAME SPATIAL ORBITALS  SPIN CASE 5   (5 IN TROFIMOV FILE)

! the condition a=a' k=k' is checked in the call  
 real(d) function D5_1_ph_2p2h(a,k,bpr,lpr)

    integer, intent(in) :: a,k,bpr,lpr

    D5_1_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
!      D5_1_ph_2p2h=dpl(lpr,bpr)
!      D5_1_ph_2p2h=D5_1_ph_2p2h/sqrt(2._d)
!     end if

  end function D5_1_ph_2p2h


! the condition a=a' k=l' is checked in the call  
 real(d) function D5_2_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    D5_2_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
    D5_2_ph_2p2h=-1.0*dpl(kpr,bpr)
    D5_2_ph_2p2h=D5_2_ph_2p2h/sqrt(2.0)
 

!    end if

  end function D5_2_ph_2p2h


! the condition a=b' k=k' is checked in the call  
 real(d) function D5_3_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    D5_3_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
    
    D5_3_ph_2p2h=-1.0*dpl(lpr,apr)
    D5_3_ph_2p2h=D5_3_ph_2p2h/sqrt(2.0)


!    end if

   end function D5_3_ph_2p2h

! the condition a=b' k=l' is checked in the call  
 real(d) function D5_4_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    D5_4_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
!    D5_4_ph_2p2h=dpl(kpr,apr)
!    D5_4_ph_2p2h=D5_4_ph_2p2h/sqrt(2._d)

!    D5_4_ph_2p2h=-D5_4_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

!    end if  

! MAYBE A MINUS FROM EXPRESSION ?

  end function D5_4_ph_2p2h


! the condition a=a' k=k' is checked in the call  
  real(d) function D5_5_ph_2p2h(a,k,bpr,lpr)
    
    integer, intent(in) :: a,k,bpr,lpr
    
    integer :: j,d
    
    real*8 :: E

    D5_5_ph_2p2h=0.0
    
!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
  
!    do d=nOcc+1,nBas
!       r=roccnum(d)
!       do j=1,nOcc
!          u=roccnum(j)

!    E=e(bpr)+e(d)-e(lpr)-e(j)

!    D5_5_ph_2p2h=D5_5_ph_2p2h+(vpqrs(lpr,bpr,j,d)-vpqrs(lpr,d,j,bpr))
!    D5_5_ph_2p2h=D5_5_ph_2p2h/E      
!    D5_5_ph_2p2h=D5_5_ph_2p2h*dpl(d,j) 
      
!      end do
!   end do
    
!    D5_5_ph_2p2h=D5_5_ph_2p2h/sqrt(2._d)


! end if

  end function D5_5_ph_2p2h

! the condition a=a' k=l' is checked in the call  
 real(d) function D5_6_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    integer :: j,d,j1,d1,sym

    real*8 :: En

    D5_6_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(bpr)+e(d)-e(kpr)-e(j)

    D5_6_ph_2p2h=D5_6_ph_2p2h+(-2.0*vpqrs(kpr,bpr,j,d)+1.0*vpqrs(kpr,d,j,bpr))
    D5_6_ph_2p2h=D5_6_ph_2p2h/En
    D5_6_ph_2p2h=D5_6_ph_2p2h*dpl(d,j)  
     
       end if
      end do
   end do

    D5_6_ph_2p2h=D5_6_ph_2p2h/sqrt(2.0)

    D5_6_ph_2p2h=-D5_6_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER
 
  end function D5_6_ph_2p2h


! the condition a=b' k=k' is checked in the call  
 real(d) function D5_7_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    integer :: j,d,j1,d1
    integer :: sym
 
    real*8 :: En

    D5_7_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(d)-e(lpr)-e(j)
 
    D5_7_ph_2p2h=D5_7_ph_2p2h+(-2.0*vpqrs(lpr,apr,j,d)+1.0*vpqrs(lpr,d,j,apr))
    D5_7_ph_2p2h=D5_7_ph_2p2h/En
    D5_7_ph_2p2h=D5_7_ph_2p2h*dpl(d,j)  
     
       end if
      end do
   end do

    D5_7_ph_2p2h=D5_7_ph_2p2h/sqrt(2.0)

    D5_7_ph_2p2h=-D5_7_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER
 
  end function D5_7_ph_2p2h


! the condition a=b' k=l' is checked in the call  
 real(d) function D5_8_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    integer :: j,d,j1,d1
    integer :: sym
 
    real*8 :: En

    D5_8_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
!    do d=nOcc+1,nBas
!       r=roccnum(d)
!       do j=1,nOcc
!          u=roccnum(j)

!    E=e(apr)+e(d)-e(kpr)-e(j)

!    D5_8_ph_2p2h=D5_8_ph_2p2h+(vpqrs(kpr,apr,j,d)-vpqrs(kpr,d,j,apr))
!    D5_8_ph_2p2h=D5_8_ph_2p2h/E
!    D5_8_ph_2p2h=D5_8_ph_2p2h*dpl(d,j)  

!      end do
!   end do

!    D5_8_ph_2p2h=D5_8_ph_2p2h/sqrt(2._d)



!    end if

  end function D5_8_ph_2p2h

! the condition a=a'  is checked in the call  
 real(d) function D5_9_ph_2p2h(a,k,bpr,kpr,lpr)

    integer, intent(in) :: a,k,bpr,kpr,lpr

    integer :: c,c1
    integer :: sym
 
     real*8 :: En
    
    D5_9_ph_2p2h=0.0

!    if(kdelta(a,apr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)

     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then
 
    En=e(c)+e(bpr)-e(kpr)-e(lpr)

    D5_9_ph_2p2h=D5_9_ph_2p2h+(0.0*vpqrs(kpr,c,lpr,bpr)+1.0*vpqrs(kpr,bpr,lpr,c))
    D5_9_ph_2p2h=D5_9_ph_2p2h/En
    D5_9_ph_2p2h=D5_9_ph_2p2h*dpl(c,k)  

       end if
      end do
!   end do

    D5_9_ph_2p2h=D5_9_ph_2p2h/sqrt(2.0)

    D5_9_ph_2p2h=-D5_9_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER


  end function D5_9_ph_2p2h


! the condition a=b'  is checked in the call  
 real(d) function D5_10_ph_2p2h(a,k,apr,kpr,lpr)

    integer, intent(in) :: a,k,apr,kpr,lpr

    integer :: c,c1
     real*8 :: En
     integer :: sym

    D5_10_ph_2p2h=0.0

!    if(kdelta(a,bpr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then
 

    En=e(c)+e(apr)-e(kpr)-e(lpr)

    D5_10_ph_2p2h=D5_10_ph_2p2h+(-1.0*vpqrs(kpr,c,lpr,apr)-0.0*vpqrs(kpr,apr,lpr,c))
    D5_10_ph_2p2h=D5_10_ph_2p2h/En
    D5_10_ph_2p2h=D5_10_ph_2p2h*dpl(c,k)  

       end if
      end do
!   end do

    D5_10_ph_2p2h=D5_10_ph_2p2h/sqrt(2.0)



  end function D5_10_ph_2p2h


! the condition  k=k' is checked in the call  
 real(d) function D5_11_ph_2p2h(a,k,apr,bpr,lpr)

    integer, intent(in) :: a,k,apr,bpr,lpr

    integer :: j,j1
     real*8 :: En
     integer :: sym

    D5_11_ph_2p2h=0.0

!    if(kdelta(k,kpr) .eq. 1)

       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(lpr)

    D5_11_ph_2p2h=D5_11_ph_2p2h+(0.0*vpqrs(j,apr,lpr,bpr)+1.0*vpqrs(j,bpr,lpr,apr))
    D5_11_ph_2p2h=D5_11_ph_2p2h/En
    D5_11_ph_2p2h=D5_11_ph_2p2h*dpl(a,j)  
    
       end if
      end do
!   end do

    D5_11_ph_2p2h=D5_11_ph_2p2h/sqrt(2.0)

    D5_11_ph_2p2h=-D5_11_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER
 
  end function D5_11_ph_2p2h


! the condition  k=l' is checked in the call  
 real(d) function D5_12_ph_2p2h(a,k,apr,bpr,kpr)

    integer, intent(in) :: a,k,apr,bpr,kpr

    integer :: j,j1
     real*8 :: En
     integer :: sym

    D5_12_ph_2p2h=0.0

!    if(kdelta(k,lpr) .eq. 1)


       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(kpr)

    D5_12_ph_2p2h=D5_12_ph_2p2h+(-1.0*vpqrs(j,apr,kpr,bpr)-0.0*vpqrs(j,bpr,kpr,apr))
    D5_12_ph_2p2h=D5_12_ph_2p2h/En
    D5_12_ph_2p2h=D5_12_ph_2p2h*dpl(a,j)  

       end if
      end do
!   end do

    D5_12_ph_2p2h=D5_12_ph_2p2h/sqrt(2.0)


  end function D5_12_ph_2p2h


!!$ a'|=b' and k'=l';  EMPTY THE SAME SPATIAL ORBITAL  SPIN CASE 4 (4 IN TROFIMOV FILE)
 
! the condition a=a' k=k' is checked in the call  
real(d) function D4_1_ph_2p2h(a,k,bpr,lpr)

    integer, intent(in) :: a,k,bpr,lpr

    D4_1_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
    
      D4_1_ph_2p2h=1.0*dpl(lpr,bpr)

      D4_1_ph_2p2h=D4_1_ph_2p2h/sqrt(2.0)
      D4_1_ph_2p2h=D4_1_ph_2p2h/sqrt(2.0)

      D4_1_ph_2p2h=-D4_1_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER


!  MAYBE A MINUS FROM EXPRESSION ?

  end function D4_1_ph_2p2h


! the condition a=a' k=l' is checked in the call  
 real(d) function D4_2_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    D4_2_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
    D4_2_ph_2p2h=-1.0*dpl(kpr,bpr)

    D4_2_ph_2p2h=D4_2_ph_2p2h/sqrt(2.0)
    D4_2_ph_2p2h=D4_2_ph_2p2h/sqrt(2.0)


  end function D4_2_ph_2p2h


! the condition a=b' k=k' is checked in the call  
 real(d) function D4_3_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    D4_3_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
    
    D4_3_ph_2p2h=-1.0*dpl(lpr,apr)

    D4_3_ph_2p2h=D4_3_ph_2p2h/sqrt(2.0)
    D4_3_ph_2p2h=D4_3_ph_2p2h/sqrt(2.0)


   end function D4_3_ph_2p2h

! the condition a=b' k=l' is checked in the call  
 real(d) function D4_4_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    D4_4_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
    D4_4_ph_2p2h=1.0*dpl(kpr,apr)

    D4_4_ph_2p2h=D4_4_ph_2p2h/sqrt(2.0)
    D4_4_ph_2p2h=D4_4_ph_2p2h/sqrt(2.0)
 
    D4_4_ph_2p2h=-D4_4_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

! MAYBE A MINUS FROM EXPRESSION ?

  end function D4_4_ph_2p2h


! the condition a=a' k=k' is checked in the call  
  real(d) function D4_5_ph_2p2h(a,k,bpr,lpr)
    
    integer, intent(in) :: a,k,bpr,lpr
    
    integer :: j,d,j1,d1
     real*8 :: En
      integer :: sym

    D4_5_ph_2p2h=0.0
    
!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
  
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(bpr)+e(d)-e(lpr)-e(j)

    D4_5_ph_2p2h=D4_5_ph_2p2h+(+2.0*vpqrs(lpr,bpr,j,d)-1.0*vpqrs(lpr,d,j,bpr))
    D4_5_ph_2p2h=D4_5_ph_2p2h/En    
    D4_5_ph_2p2h=D4_5_ph_2p2h*dpl(d,j) 
    
       end if      
      end do
   end do
    
    D4_5_ph_2p2h=D4_5_ph_2p2h/sqrt(2.0)
    D4_5_ph_2p2h=D4_5_ph_2p2h/sqrt(2.0)


  end function D4_5_ph_2p2h

! the condition a=a' k=l' is checked in the call  
 real(d) function D4_6_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D4_6_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(bpr)+e(d)-e(kpr)-e(j)

    D4_6_ph_2p2h=D4_6_ph_2p2h+(-2.0*vpqrs(kpr,bpr,j,d)+1.0*vpqrs(kpr,d,j,bpr))
    D4_6_ph_2p2h=D4_6_ph_2p2h/En
    D4_6_ph_2p2h=D4_6_ph_2p2h*dpl(d,j)  

        end if
      end do
   end do

    D4_6_ph_2p2h=D4_6_ph_2p2h/sqrt(2.0)
    D4_6_ph_2p2h=D4_6_ph_2p2h/sqrt(2.0)

    D4_6_ph_2p2h=-D4_6_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER


  end function D4_6_ph_2p2h

! the condition a=b' k=k' is checked in the call  
 real(d) function D4_7_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D4_7_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(d)-e(lpr)-e(j)
 
    D4_7_ph_2p2h=D4_7_ph_2p2h+(-2.0*vpqrs(lpr,apr,j,d)+1.0*vpqrs(lpr,d,j,apr))
    D4_7_ph_2p2h=D4_7_ph_2p2h/En
    D4_7_ph_2p2h=D4_7_ph_2p2h*dpl(d,j)  

       end if
      end do
   end do

    D4_7_ph_2p2h=D4_7_ph_2p2h/sqrt(2.0)
    D4_7_ph_2p2h=D4_7_ph_2p2h/sqrt(2.0)

    D4_7_ph_2p2h=-D4_7_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

  end function D4_7_ph_2p2h

! the condition a=b' k=l' is checked in the call  
 real(d) function D4_8_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D4_8_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(d)-e(kpr)-e(j)

    D4_8_ph_2p2h=D4_8_ph_2p2h+(+2.0*vpqrs(kpr,apr,j,d)-1.0*vpqrs(kpr,d,j,apr))
    D4_8_ph_2p2h=D4_8_ph_2p2h/En
    D4_8_ph_2p2h=D4_8_ph_2p2h*dpl(d,j)  

       end if
      end do
   end do

    D4_8_ph_2p2h=D4_8_ph_2p2h/sqrt(2.0)
    D4_8_ph_2p2h=D4_8_ph_2p2h/sqrt(2.0)


  end function D4_8_ph_2p2h

! the condition a=a'  is checked in the call  
 real(d) function D4_9_ph_2p2h(a,k,bpr,kpr,lpr)

    integer, intent(in) :: a,k,bpr,kpr,lpr

    integer :: c,c1
     real*8 :: En
     integer :: sym

    D4_9_ph_2p2h=0.0

!    if(kdelta(a,apr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then
 

    En=e(c)+e(bpr)-e(kpr)-e(lpr)

    D4_9_ph_2p2h=D4_9_ph_2p2h+(1.0*vpqrs(kpr,c,lpr,bpr)+1.0*vpqrs(kpr,bpr,lpr,c))
    D4_9_ph_2p2h=D4_9_ph_2p2h/En
    D4_9_ph_2p2h=D4_9_ph_2p2h*dpl(c,k)  

        end if
      end do
!   end do

    D4_9_ph_2p2h=D4_9_ph_2p2h/sqrt(2.0)
    D4_9_ph_2p2h=D4_9_ph_2p2h/sqrt(2.0)

    D4_9_ph_2p2h=-D4_9_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER

  end function D4_9_ph_2p2h

! the condition a=b'  is checked in the call  
 real(d) function D4_10_ph_2p2h(a,k,apr,kpr,lpr)

    integer, intent(in) :: a,k,apr,kpr,lpr

    integer :: c,c1
     real*8 :: En
     integer :: sym

    D4_10_ph_2p2h=0.0

!    if(kdelta(a,bpr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then
 

    En=e(c)+e(apr)-e(kpr)-e(lpr)

    D4_10_ph_2p2h=D4_10_ph_2p2h+(-1.0*vpqrs(kpr,c,lpr,apr)-1.0*vpqrs(kpr,apr,lpr,c))
    D4_10_ph_2p2h=D4_10_ph_2p2h/En
    D4_10_ph_2p2h=D4_10_ph_2p2h*dpl(c,k)  

        end if
      end do
!   end do

    D4_10_ph_2p2h=D4_10_ph_2p2h/sqrt(2.0)
    D4_10_ph_2p2h=D4_10_ph_2p2h/sqrt(2.0)


  end function D4_10_ph_2p2h

! the condition  k=k' is checked in the call  
 real(d) function D4_11_ph_2p2h(a,k,apr,bpr,lpr)

    integer, intent(in) :: a,k,apr,bpr,lpr

    integer :: j,j1
     real*8 :: En
     integer :: sym

    D4_11_ph_2p2h=0.0

!    if(kdelta(k,kpr) .eq. 1)

       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(lpr)

    D4_11_ph_2p2h=D4_11_ph_2p2h+(1.0*vpqrs(j,apr,lpr,bpr)+1.0*vpqrs(j,bpr,lpr,apr))
    D4_11_ph_2p2h=D4_11_ph_2p2h/En
    D4_11_ph_2p2h=D4_11_ph_2p2h*dpl(a,j)  

        end if
      end do
!   end do

    D4_11_ph_2p2h=D4_11_ph_2p2h/sqrt(2.0)
    D4_11_ph_2p2h=D4_11_ph_2p2h/sqrt(2.0)

    D4_11_ph_2p2h=-D4_11_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER

  end function D4_11_ph_2p2h

! the condition  k=l' is checked in the call  
 real(d) function D4_12_ph_2p2h(a,k,apr,bpr,kpr)

    integer, intent(in) :: a,k,apr,bpr,kpr

    integer :: j,j1
     real*8 :: En
     integer :: sym
 
    D4_12_ph_2p2h=0.0

!    if(kdelta(k,lpr) .eq. 1)


       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(kpr)

    D4_12_ph_2p2h=D4_12_ph_2p2h+(-1.0*vpqrs(j,apr,kpr,bpr)-1.0*vpqrs(j,bpr,kpr,apr))
    D4_12_ph_2p2h=D4_12_ph_2p2h/En
    D4_12_ph_2p2h=D4_12_ph_2p2h*dpl(a,j)  

        end if
      end do
!   end do

    D4_12_ph_2p2h=D4_12_ph_2p2h/sqrt(2.0)
    D4_12_ph_2p2h=D4_12_ph_2p2h/sqrt(2.0)


  end function D4_12_ph_2p2h

!!$ a'=b' and k'/=l';  FILL THE SAME SPATIAL ORBITAL  SPIN CASE 3 (3 IN TROFIMOV FILE)

! the condition a=a' k=k' is checked in the call  
real(d) function D3_1_ph_2p2h(a,k,bpr,lpr)

    integer, intent(in) :: a,k,bpr,lpr

    D3_1_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
    
      D3_1_ph_2p2h=-1.0*dpl(lpr,bpr)

      D3_1_ph_2p2h=D3_1_ph_2p2h/sqrt(2.0)
      D3_1_ph_2p2h=D3_1_ph_2p2h/sqrt(2.0)

      D3_1_ph_2p2h=-D3_1_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

!  MAYBE A MINUS FROM EXPRESSION ?

  end function D3_1_ph_2p2h


! the condition a=a' k=l' is checked in the call  
 real(d) function D3_2_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    D3_2_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
    D3_2_ph_2p2h=+1.0*dpl(kpr,bpr)

    D3_2_ph_2p2h=D3_2_ph_2p2h/sqrt(2.0)
    D3_2_ph_2p2h=D3_2_ph_2p2h/sqrt(2.0)


  end function D3_2_ph_2p2h


! the condition a=b' k=k' is checked in the call  
 real(d) function D3_3_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    D3_3_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
    
    D3_3_ph_2p2h=+1.0*dpl(lpr,apr)

    D3_3_ph_2p2h=D3_3_ph_2p2h/sqrt(2.0)
    D3_3_ph_2p2h=D3_3_ph_2p2h/sqrt(2.0)


   end function D3_3_ph_2p2h

! the condition a=b' k=l' is checked in the call  
 real(d) function D3_4_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    D3_4_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
    D3_4_ph_2p2h=-1.0*dpl(kpr,apr)

    D3_4_ph_2p2h=D3_4_ph_2p2h/sqrt(2.0)
    D3_4_ph_2p2h=D3_4_ph_2p2h/sqrt(2.0)

    D3_4_ph_2p2h=-D3_4_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

! MAYBE A MINUS FROM EXPRESSION ?

  end function D3_4_ph_2p2h

! the condition a=a' k=k' is checked in the call  
  real(d) function D3_5_ph_2p2h(a,k,bpr,lpr)
    
    integer, intent(in) :: a,k,bpr,lpr
    
    integer :: j,d,j1,d1
     real*8 :: En
      integer :: sym

    D3_5_ph_2p2h=0.0
    
!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
  
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(bpr)+e(d)-e(lpr)-e(j)

    D3_5_ph_2p2h=D3_5_ph_2p2h+(-2.0*vpqrs(lpr,bpr,j,d)+1.0*vpqrs(lpr,d,j,bpr))
    D3_5_ph_2p2h=D3_5_ph_2p2h/En     
    D3_5_ph_2p2h=D3_5_ph_2p2h*dpl(d,j) 
     
        end if
      end do
   end do
    
    D3_5_ph_2p2h=D3_5_ph_2p2h/sqrt(2.0)
    D3_5_ph_2p2h=D3_5_ph_2p2h/sqrt(2.0)
 
  end function D3_5_ph_2p2h

! the condition a=a' k=l' is checked in the call  
 real(d) function D3_6_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D3_6_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(bpr)+e(d)-e(kpr)-e(j)

    D3_6_ph_2p2h=D3_6_ph_2p2h+(+2.0*vpqrs(kpr,bpr,j,d)-1.0*vpqrs(kpr,d,j,bpr))
    D3_6_ph_2p2h=D3_6_ph_2p2h/En
    D3_6_ph_2p2h=D3_6_ph_2p2h*dpl(d,j)  

       end if 
      end do
   end do

    D3_6_ph_2p2h=D3_6_ph_2p2h/sqrt(2.0)
    D3_6_ph_2p2h=D3_6_ph_2p2h/sqrt(2.0)

    D3_6_ph_2p2h=-D3_6_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

  end function D3_6_ph_2p2h

! the condition a=b' k=k' is checked in the call  
 real(d) function D3_7_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D3_7_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(d)-e(lpr)-e(j)
 
    D3_7_ph_2p2h=D3_7_ph_2p2h+(+2.0*vpqrs(lpr,apr,j,d)-1.0*vpqrs(lpr,d,j,apr))
    D3_7_ph_2p2h=D3_7_ph_2p2h/En
    D3_7_ph_2p2h=D3_7_ph_2p2h*dpl(d,j)  

        end if
      end do
   end do

    D3_7_ph_2p2h=D3_7_ph_2p2h/sqrt(2.0)
    D3_7_ph_2p2h=D3_7_ph_2p2h/sqrt(2.0)

    D3_7_ph_2p2h=-D3_7_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

  end function D3_7_ph_2p2h

! the condition a=b' k=l' is checked in the call  
 real(d) function D3_8_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D3_8_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(d)-e(kpr)-e(j)

    D3_8_ph_2p2h=D3_8_ph_2p2h+(-2.0*vpqrs(kpr,apr,j,d)+1.0*vpqrs(kpr,d,j,apr))
    D3_8_ph_2p2h=D3_8_ph_2p2h/En
    D3_8_ph_2p2h=D3_8_ph_2p2h*dpl(d,j)  
 
       end if
      end do
   end do

    D3_8_ph_2p2h=D3_8_ph_2p2h/sqrt(2.0)
    D3_8_ph_2p2h=D3_8_ph_2p2h/sqrt(2.0)


  end function D3_8_ph_2p2h

! the condition a=a'  is checked in the call  
 real(d) function D3_9_ph_2p2h(a,k,bpr,kpr,lpr)

    integer, intent(in) :: a,k,bpr,kpr,lpr
    
     integer :: c,c1
     real*8 :: En
     integer :: sym

    D3_9_ph_2p2h=0.0

!    if(kdelta(a,apr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then
 
    En=e(c)+e(bpr)-e(kpr)-e(lpr)

    D3_9_ph_2p2h=D3_9_ph_2p2h+(-1.0*vpqrs(kpr,c,lpr,bpr)-1.0*vpqrs(kpr,bpr,lpr,c))
    D3_9_ph_2p2h=D3_9_ph_2p2h/En
    D3_9_ph_2p2h=D3_9_ph_2p2h*dpl(c,k)  
       end if
 
      end do
!   end do

    D3_9_ph_2p2h=D3_9_ph_2p2h/sqrt(2.0)
    D3_9_ph_2p2h=D3_9_ph_2p2h/sqrt(2.0)

    D3_9_ph_2p2h=-D3_9_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER

  end function D3_9_ph_2p2h

! the condition a=b'  is checked in the call  
 real(d) function D3_10_ph_2p2h(a,k,apr,kpr,lpr)

    integer, intent(in) :: a,k,apr,kpr,lpr

    integer :: c,c1
     real*8 :: En
     integer :: sym

    D3_10_ph_2p2h=0.0

!    if(kdelta(a,bpr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then

    En=e(c)+e(apr)-e(kpr)-e(lpr)

    D3_10_ph_2p2h=D3_10_ph_2p2h+(+1.0*vpqrs(kpr,c,lpr,apr)+1.0*vpqrs(kpr,apr,lpr,c))
    D3_10_ph_2p2h=D3_10_ph_2p2h/En
    D3_10_ph_2p2h=D3_10_ph_2p2h*dpl(c,k)  
       end if
 
      end do
!   end do

    D3_10_ph_2p2h=D3_10_ph_2p2h/sqrt(2.0)
    D3_10_ph_2p2h=D3_10_ph_2p2h/sqrt(2.0)

  end function D3_10_ph_2p2h

! the condition  k=k' is checked in the call  
 real(d) function D3_11_ph_2p2h(a,k,apr,bpr,lpr)

    integer, intent(in) :: a,k,apr,bpr,lpr

    integer :: j,j1
     real*8 :: En
     integer :: sym

    D3_11_ph_2p2h=0.0

!    if(kdelta(k,kpr) .eq. 1)

       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(lpr)

    D3_11_ph_2p2h=D3_11_ph_2p2h+(-1.0*vpqrs(j,apr,lpr,bpr)-1.0*vpqrs(j,bpr,lpr,apr))
    D3_11_ph_2p2h=D3_11_ph_2p2h/En
    D3_11_ph_2p2h=D3_11_ph_2p2h*dpl(a,j)  
       end if
 
      end do
!   end do

    D3_11_ph_2p2h=D3_11_ph_2p2h/sqrt(2.0)
    D3_11_ph_2p2h=D3_11_ph_2p2h/sqrt(2.0)

    D3_11_ph_2p2h=-D3_11_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER

  end function D3_11_ph_2p2h

! the condition  k=l' is checked in the call  
 real(d) function D3_12_ph_2p2h(a,k,apr,bpr,kpr)

    integer, intent(in) :: a,k,apr,bpr,kpr

    integer :: j,j1
     real*8 :: En
     integer :: sym

    D3_12_ph_2p2h=0.0

!    if(kdelta(k,lpr) .eq. 1)


       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(kpr)

    D3_12_ph_2p2h=D3_12_ph_2p2h+(+1.0*vpqrs(j,apr,kpr,bpr)+1.0*vpqrs(j,bpr,kpr,apr))
    D3_12_ph_2p2h=D3_12_ph_2p2h/En
    D3_12_ph_2p2h=D3_12_ph_2p2h*dpl(a,j)  
       end if
 
      end do
!   end do

    D3_12_ph_2p2h=D3_12_ph_2p2h/sqrt(2.0)
    D3_12_ph_2p2h=D3_12_ph_2p2h/sqrt(2.0)


  end function D3_12_ph_2p2h


!!$ a'/=b' and k'/=l'; EMPTY AND  FILL DIFFERENT  SPATIAL ORBITALS  SPIN CASE 2 (2 IN TROFIMOV FILE) sqrt(1/12)  IS THE NORMALIZATION OF THE STATE, THE STATE IS THE A STATE IN SZABO AND OSTLUND

! the condition a=a' k=k' is checked in the call  
real(d) function D2_1_ph_2p2h(a,k,bpr,lpr)

    integer, intent(in) :: a,k,bpr,lpr

    D2_1_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
    
      D2_1_ph_2p2h=-6.0*dpl(lpr,bpr)

      D2_1_ph_2p2h=D2_1_ph_2p2h/sqrt(2.0)
      D2_1_ph_2p2h=D2_1_ph_2p2h/sqrt(12.0)

      D2_1_ph_2p2h=-D2_1_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

!  MAYBE A MINUS FROM EXPRESSION ?

  end function D2_1_ph_2p2h

! the condition a=a' k=l' is checked in the call  
 real(d) function D2_2_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    D2_2_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
    D2_2_ph_2p2h=-6.0*dpl(kpr,bpr)

    D2_2_ph_2p2h=D2_2_ph_2p2h/sqrt(2.0)
    D2_2_ph_2p2h=D2_2_ph_2p2h/sqrt(12.0)


  end function D2_2_ph_2p2h

! the condition a=b' k=k' is checked in the call  
 real(d) function D2_3_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    D2_3_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
    
    D2_3_ph_2p2h=-6.0*dpl(lpr,apr)

    D2_3_ph_2p2h=D2_3_ph_2p2h/sqrt(2.0)
    D2_3_ph_2p2h=D2_3_ph_2p2h/sqrt(12.0)


   end function D2_3_ph_2p2h

! the condition a=b' k=l' is checked in the call  
 real(d) function D2_4_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    D2_4_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
    D2_4_ph_2p2h=-6.0*dpl(kpr,apr)

    D2_4_ph_2p2h=D2_4_ph_2p2h/sqrt(2.0)
    D2_4_ph_2p2h=D2_4_ph_2p2h/sqrt(12.0)

    D2_4_ph_2p2h=-D2_4_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

! MAYBE A MINUS FROM EXPRESSION ?

  end function D2_4_ph_2p2h

! the condition a=a' k=k' is checked in the call  
  real(d) function D2_5_ph_2p2h(a,k,bpr,lpr)
    
    integer, intent(in) :: a,k,bpr,lpr
    
    integer :: j,d,j1,d1
     real*8 :: En
      integer :: sym

    D2_5_ph_2p2h=0.0
    
!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
  
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(bpr)+e(d)-e(lpr)-e(j)

    D2_5_ph_2p2h=D2_5_ph_2p2h+(-12.0*vpqrs(lpr,bpr,j,d)+6.0*vpqrs(lpr,d,j,bpr))
    D2_5_ph_2p2h=D2_5_ph_2p2h/En      
    D2_5_ph_2p2h=D2_5_ph_2p2h*dpl(d,j) 
        end if
      
      end do
   end do
    
    D2_5_ph_2p2h=D2_5_ph_2p2h/sqrt(2.0)
    D2_5_ph_2p2h=D2_5_ph_2p2h/sqrt(12.0)


  end function D2_5_ph_2p2h

! the condition a=a' k=l' is checked in the call  
 real(d) function D2_6_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D2_6_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(bpr)+e(d)-e(kpr)-e(j)

    D2_6_ph_2p2h=D2_6_ph_2p2h+(-12.0*vpqrs(kpr,bpr,j,d)+6.0*vpqrs(kpr,d,j,bpr))
    D2_6_ph_2p2h=D2_6_ph_2p2h/En
    D2_6_ph_2p2h=D2_6_ph_2p2h*dpl(d,j)  
       end if
 
      end do
   end do

    D2_6_ph_2p2h=D2_6_ph_2p2h/sqrt(2.0)
    D2_6_ph_2p2h=D2_6_ph_2p2h/sqrt(12.0)

    D2_6_ph_2p2h=-D2_6_ph_2p2h  !!! FACTOR FROM EXPRESSION IN THE PAPER

  end function D2_6_ph_2p2h

! the condition a=b' k=k' is checked in the call  
 real(d) function D2_7_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D2_7_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(d)-e(lpr)-e(j)
 
    D2_7_ph_2p2h=D2_7_ph_2p2h+(-12.0*vpqrs(lpr,apr,j,d)+6.0*vpqrs(lpr,d,j,apr))
    D2_7_ph_2p2h=D2_7_ph_2p2h/En
    D2_7_ph_2p2h=D2_7_ph_2p2h*dpl(d,j)  
       end if
 
      end do
   end do

    D2_7_ph_2p2h=D2_7_ph_2p2h/sqrt(2.0)
    D2_7_ph_2p2h=D2_7_ph_2p2h/sqrt(12.0)

    D2_7_ph_2p2h=-D2_7_ph_2p2h  !!! FACTOR FROM EXPRESSION IN THE PAPER

  end function D2_7_ph_2p2h

! the condition a=b' k=l' is checked in the call  
 real(d) function D2_8_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D2_8_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(d)-e(kpr)-e(j)

    D2_8_ph_2p2h=D2_8_ph_2p2h+(-12.0*vpqrs(kpr,apr,j,d)+6.0*vpqrs(kpr,d,j,apr))
    D2_8_ph_2p2h=D2_8_ph_2p2h/En
    D2_8_ph_2p2h=D2_8_ph_2p2h*dpl(d,j)  
       end if
 
      end do
   end do

    D2_8_ph_2p2h=D2_8_ph_2p2h/sqrt(2.0)
    D2_8_ph_2p2h=D2_8_ph_2p2h/sqrt(12.0)

  end function D2_8_ph_2p2h

! the condition a=a'  is checked in the call  
 real(d) function D2_9_ph_2p2h(a,k,bpr,kpr,lpr)

    integer, intent(in) :: a,k,bpr,kpr,lpr

    integer :: c,c1
     real*8 :: En
     integer :: sym

    D2_9_ph_2p2h=0.0

!    if(kdelta(a,apr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then
 
    En=e(c)+e(bpr)-e(kpr)-e(lpr)

    D2_9_ph_2p2h=D2_9_ph_2p2h+(-6.0*vpqrs(kpr,c,lpr,bpr)+6.0*vpqrs(kpr,bpr,lpr,c))
    D2_9_ph_2p2h=D2_9_ph_2p2h/En
    D2_9_ph_2p2h=D2_9_ph_2p2h*dpl(c,k)  
       end if
 
      end do
!   end do

    D2_9_ph_2p2h=D2_9_ph_2p2h/sqrt(2.0)
    D2_9_ph_2p2h=D2_9_ph_2p2h/sqrt(12.0)

    D2_9_ph_2p2h=-D2_9_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER

  end function D2_9_ph_2p2h

! the condition a=b' is checked in the call  
 real(d) function D2_10_ph_2p2h(a,k,apr,kpr,lpr)

    integer, intent(in) :: a,k,apr,kpr,lpr

    integer :: c,c1
     real*8 :: En
     integer :: sym

    D2_10_ph_2p2h=0.0

!    if(kdelta(a,bpr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then
 

    En=e(c)+e(apr)-e(kpr)-e(lpr)

    D2_10_ph_2p2h=D2_10_ph_2p2h+(-6.0*vpqrs(kpr,c,lpr,apr)+6.0*vpqrs(kpr,apr,lpr,c))
    D2_10_ph_2p2h=D2_10_ph_2p2h/En
    D2_10_ph_2p2h=D2_10_ph_2p2h*dpl(c,k)  
       end if
 
      end do
!   end do

    D2_10_ph_2p2h=D2_10_ph_2p2h/sqrt(2.0)
    D2_10_ph_2p2h=D2_10_ph_2p2h/sqrt(12.0)


  end function D2_10_ph_2p2h

! the condition  k=k' is checked in the call  
 real(d) function D2_11_ph_2p2h(a,k,apr,bpr,lpr)

    integer, intent(in) :: a,k,apr,bpr,lpr

    integer :: j,j1
     real*8 :: En
     integer :: sym

    D2_11_ph_2p2h=0.0

!    if(kdelta(k,kpr) .eq. 1)

       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(lpr)

    D2_11_ph_2p2h=D2_11_ph_2p2h+(-6.0*vpqrs(j,apr,lpr,bpr)+6.0*vpqrs(j,bpr,lpr,apr))
    D2_11_ph_2p2h=D2_11_ph_2p2h/En
    D2_11_ph_2p2h=D2_11_ph_2p2h*dpl(a,j)  
       end if
 
      end do
!   end do

    D2_11_ph_2p2h=D2_11_ph_2p2h/sqrt(2.0)
    D2_11_ph_2p2h=D2_11_ph_2p2h/sqrt(12.0)

    D2_11_ph_2p2h=-D2_11_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER


  end function D2_11_ph_2p2h

! the condition  k=l' is checked in the call  
 real(d) function D2_12_ph_2p2h(a,k,apr,bpr,kpr)

    integer, intent(in) :: a,k,apr,bpr,kpr

    integer :: j,j1
     real*8 :: En
     integer :: sym

    D2_12_ph_2p2h=0.0

!    if(kdelta(k,lpr) .eq. 1)

       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(kpr)

    D2_12_ph_2p2h=D2_12_ph_2p2h+(-6.0*vpqrs(j,apr,kpr,bpr)+6.0*vpqrs(j,bpr,kpr,apr))
    D2_12_ph_2p2h=D2_12_ph_2p2h/En
    D2_12_ph_2p2h=D2_12_ph_2p2h*dpl(a,j)  
       end if
 
      end do
!   end do

    D2_12_ph_2p2h=D2_12_ph_2p2h/sqrt(2.0)
    D2_12_ph_2p2h=D2_12_ph_2p2h/sqrt(12.0)

  end function D2_12_ph_2p2h


!!$ a'/=b' and k'/=l'; EMPTY AND  FILL DIFFERENT  SPATIAL ORBITALS  SPIN CASE 1 (1 IN TROFIMOV FILE) 1/2 IS THE NORMALIZATION OF THE STATE, THE STATE IS THE B STATE IN SZABO AND OSTLUND

! the condition a=a' k=k' is checked in the call  
real(d) function D1_1_ph_2p2h(a,k,bpr,lpr)

    integer, intent(in) :: a,k,bpr,lpr

    D1_1_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
    
      D1_1_ph_2p2h=-2.0*dpl(lpr,bpr)

      D1_1_ph_2p2h=D1_1_ph_2p2h/sqrt(2.0)
      D1_1_ph_2p2h=D1_1_ph_2p2h/2.0

      D1_1_ph_2p2h=-D1_1_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

!  MAYBE A MINUS FROM EXPRESSION ?

  end function D1_1_ph_2p2h

! the condition a=a' k=l' is checked in the call  
 real(d) function D1_2_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    D1_2_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
    D1_2_ph_2p2h=+2.0*dpl(kpr,bpr)

    D1_2_ph_2p2h=D1_2_ph_2p2h/sqrt(2.0)
    D1_2_ph_2p2h=D1_2_ph_2p2h/2.0

  end function D1_2_ph_2p2h

! the condition a=b' k=k' is checked in the call  
 real(d) function D1_3_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    D1_3_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
    
    D1_3_ph_2p2h=+2.0*dpl(lpr,apr)

    D1_3_ph_2p2h=D1_3_ph_2p2h/sqrt(2.0)
    D1_3_ph_2p2h=D1_3_ph_2p2h/2.0

   end function D1_3_ph_2p2h

! the condition a=b' k=l' is checked in the call  
 real(d) function D1_4_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    D1_4_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
    
    D1_4_ph_2p2h=-2.0*dpl(kpr,apr)

    D1_4_ph_2p2h=D1_4_ph_2p2h/sqrt(2.0)
    D1_4_ph_2p2h=D1_4_ph_2p2h/2.0

    D1_4_ph_2p2h=-D1_4_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

! MAYBE A MINUS FROM EXPRESSION ?

  end function D1_4_ph_2p2h

! the condition a=a' k=k' is checked in the call  
  real(d) function D1_5_ph_2p2h(a,k,bpr,lpr)
    
    integer, intent(in) :: a,k,bpr,lpr
    
    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D1_5_ph_2p2h=0.0
    
!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
  
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(bpr)+e(d)-e(lpr)-e(j)

    D1_5_ph_2p2h=D1_5_ph_2p2h+(-4.0*vpqrs(lpr,bpr,j,d)+2.0*vpqrs(lpr,d,j,bpr))
    D1_5_ph_2p2h=D1_5_ph_2p2h/En      
    D1_5_ph_2p2h=D1_5_ph_2p2h*dpl(d,j) 
        end if
      
      end do
   end do
    
    D1_5_ph_2p2h=D1_5_ph_2p2h/sqrt(2.0)
    D1_5_ph_2p2h=D1_5_ph_2p2h/2.0


  end function D1_5_ph_2p2h

! the condition a=a' k=l' is checked in the call  
 real(d) function D1_6_ph_2p2h(a,k,bpr,kpr)

    integer, intent(in) :: a,k,bpr,kpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D1_6_ph_2p2h=0.0

!    if((kdelta(a,apr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(bpr)+e(d)-e(kpr)-e(j)

    D1_6_ph_2p2h=D1_6_ph_2p2h+(+4.0*vpqrs(kpr,bpr,j,d)-2.0*vpqrs(kpr,d,j,bpr))
    D1_6_ph_2p2h=D1_6_ph_2p2h/En
    D1_6_ph_2p2h=D1_6_ph_2p2h*dpl(d,j)  
       end if
 
      end do
   end do

    D1_6_ph_2p2h=D1_6_ph_2p2h/sqrt(2.0)
    D1_6_ph_2p2h=D1_6_ph_2p2h/2.0

    D1_6_ph_2p2h=-D1_6_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

  end function D1_6_ph_2p2h

! the condition a=b' k=k' is checked in the call  
 real(d) function D1_7_ph_2p2h(a,k,apr,lpr)

    integer, intent(in) :: a,k,apr,lpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D1_7_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,kpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(d)-e(lpr)-e(j)
 
    D1_7_ph_2p2h=D1_7_ph_2p2h+(+4.0*vpqrs(lpr,apr,j,d)-2.0*vpqrs(lpr,d,j,apr))
    D1_7_ph_2p2h=D1_7_ph_2p2h/En
    D1_7_ph_2p2h=D1_7_ph_2p2h*dpl(d,j)  
       end if
 
      end do
   end do

    D1_7_ph_2p2h=D1_7_ph_2p2h/sqrt(2.0)
    D1_7_ph_2p2h=D1_7_ph_2p2h/2.0

    D1_7_ph_2p2h=-D1_7_ph_2p2h  !!! FACTOR FORM EXPRESSION IN THE PAPER

  end function D1_7_ph_2p2h

! the condition a=b' k=l' is checked in the call  
 real(d) function D1_8_ph_2p2h(a,k,apr,kpr)

    integer, intent(in) :: a,k,apr,kpr

    integer :: j,d,j1,d1
     real*8 :: En
     integer :: sym

    D1_8_ph_2p2h=0.0

!    if((kdelta(a,bpr) .eq. 1).and.(kdelta(k,lpr) .eq. 1))
   
    do d1=nOcc+1,nBas
       d=roccnum(d1)
       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(d),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(d)-e(kpr)-e(j)

    D1_8_ph_2p2h=D1_8_ph_2p2h+(-4.0*vpqrs(kpr,apr,j,d)+2.0*vpqrs(kpr,d,j,apr))
    D1_8_ph_2p2h=D1_8_ph_2p2h/En
    D1_8_ph_2p2h=D1_8_ph_2p2h*dpl(d,j)  
       end if
 
      end do
   end do

    D1_8_ph_2p2h=D1_8_ph_2p2h/sqrt(2.0)
    D1_8_ph_2p2h=D1_8_ph_2p2h/2.0

  end function D1_8_ph_2p2h

! the condition a=a'  is checked in the call  
 real(d) function D1_9_ph_2p2h(a,k,bpr,kpr,lpr)

    integer, intent(in) :: a,k,bpr,kpr,lpr

    integer :: c,c1
     real*8 :: En
     integer :: sym

    D1_9_ph_2p2h=0.0

!    if(kdelta(a,apr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then
 
    En=e(c)+e(bpr)-e(kpr)-e(lpr)

    D1_9_ph_2p2h=D1_9_ph_2p2h+(-2.0*vpqrs(kpr,c,lpr,bpr)-2.0*vpqrs(kpr,bpr,lpr,c))
    D1_9_ph_2p2h=D1_9_ph_2p2h/En
    D1_9_ph_2p2h=D1_9_ph_2p2h*dpl(c,k)  
       end if
 
      end do

    D1_9_ph_2p2h=D1_9_ph_2p2h/sqrt(2.0)
    D1_9_ph_2p2h=D1_9_ph_2p2h/2.0

    D1_9_ph_2p2h=-D1_9_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER

  end function D1_9_ph_2p2h

! the condition a=b'  is checked in the call  
 real(d) function D1_10_ph_2p2h(a,k,apr,kpr,lpr)

    integer, intent(in) :: a,k,apr,kpr,lpr

    integer :: c,c1
     real*8 :: En
     integer :: sym

    D1_10_ph_2p2h=0.0

!    if(kdelta(a,bpr) .eq. 1)

    do c1=nOcc+1,nBas
       c=roccnum(c1)
     sym=MT(orbSym(c),orbSym(k))
     if (sym .eq. CHECK_dip) then
 
    En=e(c)+e(apr)-e(kpr)-e(lpr)

    D1_10_ph_2p2h=D1_10_ph_2p2h+(+2.0*vpqrs(kpr,c,lpr,apr)+2.0*vpqrs(kpr,apr,lpr,c))
    D1_10_ph_2p2h=D1_10_ph_2p2h/En
    D1_10_ph_2p2h=D1_10_ph_2p2h*dpl(c,k)  
       end if
 
      end do
!   end do

    D1_10_ph_2p2h=D1_10_ph_2p2h/sqrt(2.0)
    D1_10_ph_2p2h=D1_10_ph_2p2h/2.0

  end function D1_10_ph_2p2h

! the condition  k=k' is checked in the call  
 real(d) function D1_11_ph_2p2h(a,k,apr,bpr,lpr)

    integer, intent(in) :: a,k,apr,bpr,lpr

    integer :: j,j1
     real*8 :: En
     integer :: sym
 
    D1_11_ph_2p2h=0._d

!    if(kdelta(k,kpr) .eq. 1)

       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(lpr)

    D1_11_ph_2p2h=D1_11_ph_2p2h+(-2.0*vpqrs(j,apr,lpr,bpr)-2.0*vpqrs(j,bpr,lpr,apr))
    D1_11_ph_2p2h=D1_11_ph_2p2h/En
    D1_11_ph_2p2h=D1_11_ph_2p2h*dpl(a,j)  
       end if
 
      end do
!   end do

    D1_11_ph_2p2h=D1_11_ph_2p2h/sqrt(2.0)
    D1_11_ph_2p2h=D1_11_ph_2p2h/2.0

    D1_11_ph_2p2h=-D1_11_ph_2p2h   !!! EXPRESSION FACTOR FROM THE PAPER

  end function D1_11_ph_2p2h

! the condition  k=l' is checked in the call  
 real(d) function D1_12_ph_2p2h(a,k,apr,bpr,kpr)

    integer, intent(in) :: a,k,apr,bpr,kpr

    integer :: j,j1
     real*8 :: En
     integer :: sym
 
    D1_12_ph_2p2h=0._d

!    if(kdelta(k,lpr) .eq. 1)

       do j1=1,nOcc
          j=roccnum(j1)
     sym=MT(orbSym(a),orbSym(j))
     if (sym .eq. CHECK_dip) then
 
    En=e(apr)+e(bpr)-e(j)-e(kpr)

    D1_12_ph_2p2h=D1_12_ph_2p2h+(+2.0*vpqrs(j,apr,kpr,bpr)+2.0*vpqrs(j,bpr,kpr,apr))
    D1_12_ph_2p2h=D1_12_ph_2p2h/En
    D1_12_ph_2p2h=D1_12_ph_2p2h*dpl(a,j)  

       end if
      end do
!   end do

    D1_12_ph_2p2h=D1_12_ph_2p2h/sqrt(2.0)
    D1_12_ph_2p2h=D1_12_ph_2p2h/2.0

  end function D1_12_ph_2p2h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! END OF THE 1PH-2P2H PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!$*******************************************************************************
!!$*******************************************************************************
!!$*******************************************************************************
!!$*******************************************************************************
!!$*********************************2PH-2PH BLOCK*********************************
!!$*******************************************************************************
!!$*******************************************************************************  
!!$*******************************************************************************
!!$*******************************************************************************

!!$ 2P2H-2P2H block is organized in subblocks according to the different
!!$ classes of double excitations. 1- doubles k=l,a=b; 2-doubles
!!$  k=l,a|=b; 3-doubles k|=l,a=b; 4ii-doubles k|=l,a|=b (spin case 2, A in Stzabo, sqrt(1/12) normalization);
!!$ 4i-doubles k|=l,a|=b (spin case 1, B in Stzabo, 1/2 normalization).


!!!!!  FIRST TERM

!!$ (1,1)-Dabkl,a'b'k'l' TERM WITH   b=b' k=k' l=l'(THE LOWER PART OF THE MATRIX)

  real(d) function D_1_1_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)
    
    integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr
    
    D_1_1_2p2h_2p2h=0.0

    if((kdelta(k,kpr).eq.1).and.(kdelta(l,lpr).eq.1)) then

    if(kdelta(b,bpr) .eq. 1) &   
    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 1.0*dpl(a,apr)
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&   
    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h  - 0.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&   
    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h  + 1.0*dpl(b,bpr)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h  - 0.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
   end  if


    if((kdelta(a,apr).eq.1).and.(kdelta(b,bpr).eq.1)) then



       if(kdelta(l,lpr) .eq. 1)&  
    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h  + 1.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h  - 0.0*dpl(lpr,k)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)& 
    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h  + 1.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)& 
    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h  - 0.0*dpl(kpr,l)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if

!    D_1_1_2p2h_2p2h= +0.25*D_1_1_2p2h_2p2h    
 

  
 
!    D_1_1_2p2h_2p2h=+4._d*dpl(a,apr)
!    D_1_1_2p2h_2p2h=+0.25*D_1_1_2p2h_2p2h    
  
  end function D_1_1_2p2h_2p2h

  !!$ (2,1)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_2_1_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

     integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

    D_2_1_2p2h_2p2h=0.0

    if((kdelta(k,kpr).eq.1).and.(kdelta(l,lpr).eq.1)) then



    if(kdelta(b,bpr) .eq. 1) &  
    D_2_1_2p2h_2p2h = D_2_1_2p2h_2p2h + 1.0*dpl(a,apr)
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)& 
    D_2_1_2p2h_2p2h = D_2_1_2p2h_2p2h  - 1.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_2_1_2p2h_2p2h = D_2_1_2p2h_2p2h  + 1.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_2_1_2p2h_2p2h = D_2_1_2p2h_2p2h  - 1.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if



    if((kdelta(a,apr).eq.1).and.(kdelta(b,bpr).eq.1))  then



       if(kdelta(l,lpr) .eq. 1)& 
    D_2_1_2p2h_2p2h = D_2_1_2p2h_2p2h  + 1.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_2_1_2p2h_2p2h = D_2_1_2p2h_2p2h  - 0.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)&  
    D_2_1_2p2h_2p2h = D_2_1_2p2h_2p2h  + 1.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)& 
    D_2_1_2p2h_2p2h = D_2_1_2p2h_2p2h  - 0.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if

    D_2_1_2p2h_2p2h= D_2_1_2p2h_2p2h*(1/sqrt(2.0))    
 
  end function D_2_1_2p2h_2p2h

   !!$ (3,1)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_3_1_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

      integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

    D_3_1_2p2h_2p2h=0.0

    if((kdelta(k,kpr).eq.1).and.(kdelta(l,lpr).eq.1))  then



    if(kdelta(b,bpr) .eq. 1) & 
    D_3_1_2p2h_2p2h = D_3_1_2p2h_2p2h - 1.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&  
    D_3_1_2p2h_2p2h = D_3_1_2p2h_2p2h  - 0.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_3_1_2p2h_2p2h = D_3_1_2p2h_2p2h  - 1.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_3_1_2p2h_2p2h = D_3_1_2p2h_2p2h  - 0.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if



    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1)) then



       if(kdelta(l,lpr) .eq. 1)&  
    D_3_1_2p2h_2p2h = D_3_1_2p2h_2p2h  -  1.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_3_1_2p2h_2p2h = D_3_1_2p2h_2p2h  + 1.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)& 
    D_3_1_2p2h_2p2h = D_3_1_2p2h_2p2h  - 1.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)& 
    D_3_1_2p2h_2p2h = D_3_1_2p2h_2p2h  + 1.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if


    D_3_1_2p2h_2p2h= D_3_1_2p2h_2p2h*(1/sqrt(2.0))    
 
  end function D_3_1_2p2h_2p2h

 !!$ (4i,1)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_4i_1_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

       integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

    D_4i_1_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) & 
    D_4i_1_2p2h_2p2h = D_4i_1_2p2h_2p2h - 1.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)& 
    D_4i_1_2p2h_2p2h = D_4i_1_2p2h_2p2h  + 1.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_4i_1_2p2h_2p2h = D_4i_1_2p2h_2p2h  - 1.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_4i_1_2p2h_2p2h = D_4i_1_2p2h_2p2h  + 1.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1)) then



       if(kdelta(l,lpr) .eq. 1)& 
    D_4i_1_2p2h_2p2h = D_4i_1_2p2h_2p2h  - 1.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_4i_1_2p2h_2p2h = D_4i_1_2p2h_2p2h  + 1.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)& 
    D_4i_1_2p2h_2p2h = D_4i_1_2p2h_2p2h  - 1.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)&   
    D_4i_1_2p2h_2p2h = D_4i_1_2p2h_2p2h  + 1.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if


    D_4i_1_2p2h_2p2h= D_4i_1_2p2h_2p2h*0.5    
 
  end function D_4i_1_2p2h_2p2h

 !!$ (4ii,1)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_4ii_1_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

        integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

    D_4ii_1_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) &  
    D_4ii_1_2p2h_2p2h = D_4ii_1_2p2h_2p2h + 1.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&  
    D_4ii_1_2p2h_2p2h = D_4ii_1_2p2h_2p2h  + 1.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_4ii_1_2p2h_2p2h = D_4ii_1_2p2h_2p2h  + 1.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_4ii_1_2p2h_2p2h = D_4ii_1_2p2h_2p2h  + 1.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1))  then



       if(kdelta(l,lpr) .eq. 1)&  
    D_4ii_1_2p2h_2p2h = D_4ii_1_2p2h_2p2h  + 1.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_4ii_1_2p2h_2p2h = D_4ii_1_2p2h_2p2h  + 1.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)& 
    D_4ii_1_2p2h_2p2h = D_4ii_1_2p2h_2p2h  + 1.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)& 
    D_4ii_1_2p2h_2p2h = D_4ii_1_2p2h_2p2h  + 1.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
     end if



    D_4ii_1_2p2h_2p2h= D_4ii_1_2p2h_2p2h*(1/sqrt(12.0))    
 
  end function D_4ii_1_2p2h_2p2h

 !!$ (2,2)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_2_2_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

        integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

    D_2_2_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) & 
    D_2_2_2p2h_2p2h = D_2_2_2p2h_2p2h  + 2.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&  
    D_2_2_2p2h_2p2h = D_2_2_2p2h_2p2h  - 2.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_2_2_2p2h_2p2h = D_2_2_2p2h_2p2h  + 2.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)& 
    D_2_2_2p2h_2p2h = D_2_2_2p2h_2p2h  - 2.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if
    

    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1)) then



       if(kdelta(l,lpr) .eq. 1)& 
    D_2_2_2p2h_2p2h = D_2_2_2p2h_2p2h  + 2.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_2_2_2p2h_2p2h = D_2_2_2p2h_2p2h  - 0.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)&  
    D_2_2_2p2h_2p2h = D_2_2_2p2h_2p2h  + 2.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)&  
    D_2_2_2p2h_2p2h = D_2_2_2p2h_2p2h  - 0.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if



    D_2_2_2p2h_2p2h= D_2_2_2p2h_2p2h*0.5    
 
  end function D_2_2_2p2h_2p2h

 !!$ (3,2)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_3_2_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

     integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

    D_3_2_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) & 
    D_3_2_2p2h_2p2h = D_3_2_2p2h_2p2h  - 1.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&  
    D_3_2_2p2h_2p2h = D_3_2_2p2h_2p2h  + 1.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)& 
    D_3_2_2p2h_2p2h = D_3_2_2p2h_2p2h  - 1.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_3_2_2p2h_2p2h = D_3_2_2p2h_2p2h  + 1.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1))  then



       if(kdelta(l,lpr) .eq. 1)&  
    D_3_2_2p2h_2p2h = D_3_2_2p2h_2p2h  - 1.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_3_2_2p2h_2p2h = D_3_2_2p2h_2p2h  + 1.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)&  
    D_3_2_2p2h_2p2h = D_3_2_2p2h_2p2h  - 1.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)&  
    D_3_2_2p2h_2p2h = D_3_2_2p2h_2p2h  + 1.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if



    D_3_2_2p2h_2p2h= D_3_2_2p2h_2p2h*0.5    
 
  end function D_3_2_2p2h_2p2h

 !!$ (4i,2)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_4i_2_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

      integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

    D_4i_2_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) & 
    D_4i_2_2p2h_2p2h = D_4i_2_2p2h_2p2h - 2.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&  
    D_4i_2_2p2h_2p2h = D_4i_2_2p2h_2p2h  + 2.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_4i_2_2p2h_2p2h = D_4i_2_2p2h_2p2h  - 2.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_4i_2_2p2h_2p2h = D_4i_2_2p2h_2p2h  + 2.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1))  then



       if(kdelta(l,lpr) .eq. 1)&  
    D_4i_2_2p2h_2p2h = D_4i_2_2p2h_2p2h  - 2.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_4i_2_2p2h_2p2h = D_4i_2_2p2h_2p2h  + 2.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)& 
    D_4i_2_2p2h_2p2h = D_4i_2_2p2h_2p2h  - 2.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)& 
    D_4i_2_2p2h_2p2h = D_4i_2_2p2h_2p2h  + 2.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if


    D_4i_2_2p2h_2p2h= D_4i_2_2p2h_2p2h*0.5    
    D_4i_2_2p2h_2p2h= D_4i_2_2p2h_2p2h*(1/sqrt(2.0))    
 

  end function D_4i_2_2p2h_2p2h

 !!$ (4ii,2)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_4ii_2_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

       integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

     D_4ii_2_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) & 
    D_4ii_2_2p2h_2p2h = D_4ii_2_2p2h_2p2h + 0.0*dpl(a,apr)
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)& 
    D_4ii_2_2p2h_2p2h = D_4ii_2_2p2h_2p2h  - 0.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_4ii_2_2p2h_2p2h = D_4ii_2_2p2h_2p2h  + 0.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)& 
    D_4ii_2_2p2h_2p2h = D_4ii_2_2p2h_2p2h  - 0.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1))  then



       if(kdelta(l,lpr) .eq. 1)& 
    D_4ii_2_2p2h_2p2h = D_4ii_2_2p2h_2p2h  + 0.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)& 
    D_4ii_2_2p2h_2p2h = D_4ii_2_2p2h_2p2h  - 0.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)& 
    D_4ii_2_2p2h_2p2h = D_4ii_2_2p2h_2p2h  + 0.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)&  
    D_4ii_2_2p2h_2p2h = D_4ii_2_2p2h_2p2h  - 0.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
     end if



    D_4ii_2_2p2h_2p2h= D_4ii_2_2p2h_2p2h*(1/sqrt(2.0))    
    D_4ii_2_2p2h_2p2h= D_4ii_2_2p2h_2p2h*(1/sqrt(12.0))    
 

  end function D_4ii_2_2p2h_2p2h

 !!$ (3,3)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_3_3_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

        integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr


    D_3_3_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1)) then



    if(kdelta(b,bpr) .eq. 1) & 
    D_3_3_2p2h_2p2h = D_3_3_2p2h_2p2h + 2.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)& 
    D_3_3_2p2h_2p2h = D_3_3_2p2h_2p2h  - 0.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_3_3_2p2h_2p2h = D_3_3_2p2h_2p2h  + 2.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_3_3_2p2h_2p2h = D_3_3_2p2h_2p2h  - 0.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1)) then



       if(kdelta(l,lpr) .eq. 1)& 
    D_3_3_2p2h_2p2h = D_3_3_2p2h_2p2h  + 2.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_3_3_2p2h_2p2h = D_3_3_2p2h_2p2h  - 2.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)& 
    D_3_3_2p2h_2p2h = D_3_3_2p2h_2p2h  + 2.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)& 
    D_3_3_2p2h_2p2h = D_3_3_2p2h_2p2h  - 2.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if



    D_3_3_2p2h_2p2h= D_3_3_2p2h_2p2h*0.5    
 
  end function D_3_3_2p2h_2p2h

 !!$ (4i,3)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_4i_3_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

        integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

     D_4i_3_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) & 
    D_4i_3_2p2h_2p2h = D_4i_3_2p2h_2p2h  + 2.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&
    D_4i_3_2p2h_2p2h = D_4i_3_2p2h_2p2h  - 2.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_4i_3_2p2h_2p2h = D_4i_3_2p2h_2p2h  + 2.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)& 
    D_4i_3_2p2h_2p2h = D_4i_3_2p2h_2p2h  - 2.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1)) then



       if(kdelta(l,lpr) .eq. 1)& 
    D_4i_3_2p2h_2p2h = D_4i_3_2p2h_2p2h  + 2.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)& 
    D_4i_3_2p2h_2p2h = D_4i_3_2p2h_2p2h  - 2.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)& 
    D_4i_3_2p2h_2p2h = D_4i_3_2p2h_2p2h  + 2.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)&  
    D_4i_3_2p2h_2p2h = D_4i_3_2p2h_2p2h  - 2.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if




    D_4i_3_2p2h_2p2h= D_4i_3_2p2h_2p2h*0.5    
    D_4i_3_2p2h_2p2h= D_4i_3_2p2h_2p2h*(1/sqrt(2.0))    
 

  end function D_4i_3_2p2h_2p2h

 !!$ (4ii,3)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_4ii_3_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

        integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

     D_4ii_3_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) &  
    D_4ii_3_2p2h_2p2h = D_4ii_3_2p2h_2p2h  + 0.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)& 
    D_4ii_3_2p2h_2p2h = D_4ii_3_2p2h_2p2h  - 0.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)& 
    D_4ii_3_2p2h_2p2h = D_4ii_3_2p2h_2p2h  + 0.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_4ii_3_2p2h_2p2h = D_4ii_3_2p2h_2p2h  - 0.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1))  then



       if(kdelta(l,lpr) .eq. 1)& 
    D_4ii_3_2p2h_2p2h = D_4ii_3_2p2h_2p2h  + 0.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_4ii_3_2p2h_2p2h = D_4ii_3_2p2h_2p2h  - 0.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)&  
    D_4ii_3_2p2h_2p2h = D_4ii_3_2p2h_2p2h  + 0.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)&  
    D_4ii_3_2p2h_2p2h = D_4ii_3_2p2h_2p2h  - 0.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if



    D_4ii_3_2p2h_2p2h= D_4ii_3_2p2h_2p2h*(1/sqrt(12.0))
    D_4ii_3_2p2h_2p2h= D_4ii_3_2p2h_2p2h*(1/sqrt(2.0))    
 
 end function D_4ii_3_2p2h_2p2h

 !!$ (4i,4i)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_4i_4i_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

        integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr


     D_4i_4i_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1)) then



    if(kdelta(b,bpr) .eq. 1) & 
    D_4i_4i_2p2h_2p2h = D_4i_4i_2p2h_2p2h  + 4.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&  
    D_4i_4i_2p2h_2p2h = D_4i_4i_2p2h_2p2h  - 4.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_4i_4i_2p2h_2p2h = D_4i_4i_2p2h_2p2h  + 4.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&  
    D_4i_4i_2p2h_2p2h = D_4i_4i_2p2h_2p2h  - 4.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1))  then



       if(kdelta(l,lpr) .eq. 1)&  
    D_4i_4i_2p2h_2p2h = D_4i_4i_2p2h_2p2h  + 4.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)& 
    D_4i_4i_2p2h_2p2h = D_4i_4i_2p2h_2p2h  - 4.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)&  
    D_4i_4i_2p2h_2p2h = D_4i_4i_2p2h_2p2h  + 4.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)&  
    D_4i_4i_2p2h_2p2h = D_4i_4i_2p2h_2p2h  - 4.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if



    D_4i_4i_2p2h_2p2h= D_4i_4i_2p2h_2p2h*0.25    
 
 end function D_4i_4i_2p2h_2p2h

 !!$ (4ii,4i)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_4ii_4i_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)

        integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

 


     D_4ii_4i_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) & 
    D_4ii_4i_2p2h_2p2h = D_4ii_4i_2p2h_2p2h  + 0.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&  
    D_4ii_4i_2p2h_2p2h = D_4ii_4i_2p2h_2p2h  - 0.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)& 
    D_4ii_4i_2p2h_2p2h = D_4ii_4i_2p2h_2p2h  + 0.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)& 
    D_4ii_4i_2p2h_2p2h = D_4ii_4i_2p2h_2p2h  - 0.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1))  then



       if(kdelta(l,lpr) .eq. 1)& 
    D_4ii_4i_2p2h_2p2h = D_4ii_4i_2p2h_2p2h  + 0.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_4ii_4i_2p2h_2p2h = D_4ii_4i_2p2h_2p2h  - 0.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)& 
    D_4ii_4i_2p2h_2p2h = D_4ii_4i_2p2h_2p2h  + 0.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)&  
    D_4ii_4i_2p2h_2p2h = D_4ii_4i_2p2h_2p2h  - 0.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
    end if



    D_4ii_4i_2p2h_2p2h= D_4ii_4i_2p2h_2p2h*(1/sqrt(12.0))    
    D_4ii_4i_2p2h_2p2h= D_4ii_4i_2p2h_2p2h*0.5    
 

  end function D_4ii_4i_2p2h_2p2h

 !!$ (4ii,4ii)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

  real(d) function D_4ii_4ii_2p2h_2p2h(a,b,k,l,apr,bpr,kpr,lpr)
    
    integer , intent(in) :: a,b,k,l,apr,bpr,kpr,lpr

 


     D_4ii_4ii_2p2h_2p2h=0.0

    if((kdelta(k,kpr) .eq. 1).and.(kdelta(l,lpr) .eq. 1))  then



    if(kdelta(b,bpr) .eq. 1) & 
    D_4ii_4ii_2p2h_2p2h = D_4ii_4ii_2p2h_2p2h  + 12.0*dpl(a,apr)!*0.25
!    D_1_1_2p2h_2p2h = D_1_1_2p2h_2p2h + 0.25*D_1_1_2p2h_2p2h    
 
    if(kdelta(b,apr) .eq. 1)&  
    D_4ii_4ii_2p2h_2p2h = D_4ii_4ii_2p2h_2p2h  + 12.0*dpl(a,bpr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h    
  
    if(kdelta(a,apr) .eq. 1)&  
    D_4ii_4ii_2p2h_2p2h = D_4ii_4ii_2p2h_2p2h  + 12.0*dpl(b,bpr)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h
 

    if(kdelta(a,bpr) .eq. 1)&
    D_4ii_4ii_2p2h_2p2h = D_4ii_4ii_2p2h_2p2h  + 12.0*dpl(b,apr)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)
 
    end if


    if((kdelta(a,apr) .eq. 1).and.(kdelta(b,bpr) .eq. 1))  then



       if(kdelta(l,lpr) .eq. 1)&  
    D_4ii_4ii_2p2h_2p2h = D_4ii_4ii_2p2h_2p2h  + 12.0*dpl(kpr,k)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  - 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(l,kpr) .eq. 1)&  
    D_4ii_4ii_2p2h_2p2h = D_4ii_4ii_2p2h_2p2h  + 12.0*dpl(lpr,k)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    


       if(kdelta(k,kpr) .eq. 1)&  
    D_4ii_4ii_2p2h_2p2h = D_4ii_4ii_2p2h_2p2h  + 12.0*dpl(lpr,l)*(-1.0)
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h*(-1._d)

 
       if(kdelta(k,lpr) .eq. 1)&  
    D_4ii_4ii_2p2h_2p2h = D_4ii_4ii_2p2h_2p2h  + 12.0*dpl(kpr,l)!*0.25
!    D_1_1_2p2h_2p2h= D_1_1_2p2h_2p2h  + 0.25*D_1_1_2p2h_2p2h    
  
     end if    



    D_4ii_4ii_2p2h_2p2h= D_4ii_4ii_2p2h_2p2h*(1/12)    
 
  end function D_4ii_4ii_2p2h_2p2h


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!$******************DIAGONAL**************************
!!$****

!  real(d) function DD_1_1_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: term
!    integer :: j
    
!    DD_1_1_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+4._d*dpl(j,j)

!    end do

!    DD_1_1_2p2h_2p2h=DD_1_1_2p2h_2p2h+term
!    DD_1_1_2p2h_2p2h=DD_1_1_2p2h_2p2h+4._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))
!    DD_1_1_2p2h_2p2h=+0.25*DD_1_1_2p2h_2p2h

!  end function DD_1_1_2p2h_2p2h

  !!$ (2,1)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_2_1_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: term
!    integer :: j

!    DD_2_1_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+0._d*dpl(j,j)

!    end do

!    DD_2_1_2p2h_2p2h=DD_2_1_2p2h_2p2h+term
!    DD_2_1_2p2h_2p2h=DD_2_1_2p2h_2p2h+0._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))

!    DD_2_1_2p2h_2p2h=+0.5*DD_2_1_2p2h_2p2h
!    DD_2_1_2p2h_2p2h=+(1/sqrt(12))*DD_2_1_2p2h_2p2h

!  end function DD_2_1_2p2h_2p2h

   !!$ (3,1)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_3_1_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_3_1_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+2._d*dpl(j,j)

!    end do

!    DD_3_1_2p2h_2p2h=DD_3_1_2p2h_2p2h+term
!    DD_3_1_2p2h_2p2h=DD_3_1_2p2h_2p2h+2._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))

!    DD_3_1_2p2h_2p2h=+0.5*DD_3_1_2p2h_2p2h
!    DD_3_1_2p2h_2p2h=+(1/sqrt(2))*DD_3_1_2p2h_2p2h

!  end function DD_3_1_2p2h_2p2h

 !!$ (4i,1)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_4i_1_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_4i_1_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term-2._d*dpl(j,j)

!    end do

!    DD_4i_1_2p2h_2p2h=DD_4i_1_2p2h_2p2h+term
!    DD_4i_1_2p2h_2p2h=DD_4i_1_2p2h_2p2h-2._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))
 
!    DD_4i_1_2p2h_2p2h=+0.5*DD_4i_1_2p2h_2p2h
!    DD_4i_1_2p2h_2p2h=+(1/sqrt(2))*DD_4i_1_2p2h_2p2h

!  end function DD_4i_1_2p2h_2p2h

 !!$ (4ii,1)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_4ii_1_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_4ii_1_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term-1._d*dpl(j,j)

!    end do

!    DD_4ii_1_2p2h_2p2h=DD_4ii_1_2p2h_2p2h+term
!    DD_4ii_1_2p2h_2p2h=DD_4ii_1_2p2h_2p2h-1._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))

!    DD_4ii_1_2p2h_2p2h=+0.5*DD_4ii_1_2p2h_2p2h

!  end function DD_4ii_1_2p2h_2p2h

 !!$ (2,2)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_2_2_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_2_2_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+12._d*dpl(j,j)

!    end do

!    DD_2_2_2p2h_2p2h=DD_2_2_2p2h_2p2h+term
!    DD_2_2_2p2h_2p2h=DD_2_2_2p2h_2p2h+12._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))

!    DD_2_2_2p2h_2p2h=+(1/12)*DD_2_2_2p2h_2p2h

!  end function DD_2_2_2p2h_2p2h

 !!$ (3,2)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_3_2_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_3_2_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+0._d*dpl(j,j)

!    end do

!    DD_3_2_2p2h_2p2h=DD_3_2_2p2h_2p2h+term
!    DD_3_2_2p2h_2p2h=DD_3_2_2p2h_2p2h+0._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))
 
!    DD_3_2_2p2h_2p2h=+(1/sqrt(12))*DD_3_2_2p2h_2p2h
!    DD_3_2_2p2h_2p2h=+(1/sqrt(2))*DD_3_2_2p2h_2p2h

!  end function DD_3_2_2p2h_2p2h

 !!$ (4i,2)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_4i_2_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_4i_2_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+0._d*dpl(j,j)

!    end do

!    DD_4i_2_2p2h_2p2h=DD_4i_2_2p2h_2p2h+term
!    DD_4i_2_2p2h_2p2h=DD_4i_2_2p2h_2p2h+0._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))

!    DD_4i_2_2p2h_2p2h=+(1/sqrt(12))*DD_4i_2_2p2h_2p2h
!    DD_4i_2_2p2h_2p2h=+(1/sqrt(2))*DD_4i_2_2p2h_2p2h

!  end function DD_4i_2_2p2h_2p2h

 !!$ (4ii,2)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_4ii_2_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_4ii_2_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+1._d*dpl(j,j)

!    end do

!    DD_4ii_2_2p2h_2p2h=DD_4ii_2_2p2h_2p2h+term
!    DD_4ii_2_2p2h_2p2h=DD_4ii_2_2p2h_2p2h+1._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))
 
!    DD_4ii_2_2p2h_2p2h=+(1/sqrt(12)*DD_4ii_2_2p2h_2p2h

!  end function DD_4ii_2_2p2h_2p2h

 !!$ (3,3)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_3_3_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_3_3_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+2._d*dpl(j,j)

!    end do

!    DD_3_3_2p2h_2p2h=DD_3_3_2p2h_2p2h+term
!    DD_3_3_2p2h_2p2h=DD_3_3_2p2h_2p2h+2._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))

!    DD_3_3_2p2h_2p2h=+0.5*DD_3_3_2p2h_2p2h

!  end function DD_3_3_2p2h_2p2h

 !!$ (4i,3)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_4i_3_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_4i_3_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc
!
!    term=term-1._d*dpl(j,j)

!    end do

!    DD_4i_3_2p2h_2p2h=DD_4i_3_2p2h_2p2h+term
!    DD_4i_3_2p2h_2p2h=DD_4i_3_2p2h_2p2h-1._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))
 
!    DD_4i_3_2p2h_2p2h=+0.5*DD_4i_3_2p2h_2p2h

!  end function DD_4i_3_2p2h_2p2h

 !!$ (4ii,3)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_4ii_3_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_4ii_3_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term-1._d*dpl(j,j)

!    end do

!    DD_4ii_3_2p2h_2p2h=DD_4ii_3_2p2h_2p2h+term
!    DD_4ii_3_2p2h_2p2h=DD_4ii_3_2p2h_2p2h-1._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))

!    DD_4ii_3_2p2h_2p2h=+(1/sqrt(2))*DD_4ii_3_2p2h_2p2h
 
! end function DD_4ii_3_2p2h_2p2h

 !!$ (4i,4i)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_4i_4i_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_4i_4i_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+2._d*dpl(j,j)

!    end do

!    DD_4i_4i_2p2h_2p2h=DD_4i_4i_2p2h_2p2h+term
!    DD_4i_4i_2p2h_2p2h=DD_4i_4i_2p2h_2p2h+2._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))

!    DD_4i_4i_2p2h_2p2h=+0.5*DD_4i_4i_2p2h_2p2h
 
! end function DD_4i_4i_2p2h_2p2h

 !!$ (4ii,4i)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_4ii_4i_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_4ii_4i_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+1._d*dpl(j,j)

!    end do

!    DD_4ii_4i_2p2h_2p2h=DD_4ii_4i_2p2h_2p2h+term
!    DD_4ii_4i_2p2h_2p2h=DD_4ii_4i_2p2h_2p2h+1._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))
 
!    DD_4ii_4i_2p2h_2p2h=+(1/sqrt(2))*DD_4ii_4i_2p2h_2p2h

!  end function DD_4ii_4i_2p2h_2p2h

 !!$ (4ii,4ii)-Dabkl,a'b'k'l'  b=b' k=k' l=l'

!  real(d) function DD_4ii_4ii_2p2h_2p2h(a,b,k,l)
    
!    real(d), intent(in) :: a,b,k,l
!    real(d) :: dpl
!    external dpl
!    real(d) :: term
!    integer :: j

!    DD_4ii_4ii_2p2h_2p2h=+0._d

!    term=+0._d

!    do j=1,nOcc

!    term=term+1._d*dpl(j,j)

!    end do

!    DD_4ii_4ii_2p2h_2p2h=DD_4ii_4ii_2p2h_2p2h+term
!    DD_4ii_4ii_2p2h_2p2h=DD_4ii_4ii_2p2h_2p2h+1._d*(dpl(a,a)+dpl(b,b)-dpl(k,k)-dpl(l,l))

!    DD_4ii_4ii_2p2h_2p2h=+1._d**DD_4ii_4ii_2p2h_2p2h

!  end function DD_4ii_4ii_2p2h_2p2h
    


!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF
!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF
!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF
!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF
!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF
!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF
!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF
!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF
!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF
!!!!!!          SECOND ORDER PART THAT CANCELS ITSELF


!!$Second order contribution D2_1_1_ak,a'k'. The condition that k=k' is
!!$checked in the calling procedure.

!  real(d) function D2_1_1_ph_ph(a,apr)

!    integer, intent(in) :: a,apr

!    integer :: b,c,i,j, nsym1,nsym2,nsym3,u,v,r,s,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_1_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do c=nOcc+1,nBas
!       s=roccnum(c)
!       do i=1,nOcc
!          u=roccnum(i)
!          do j=1,nOcc
!             v=roccnum(j)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(a)+e(c)-e(i)-e(j)
!                ebcij=e(b)+e(c)-e(i)-e(j)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(a,i,c,j)*(+4._d*vpqrs(i,b,j,c)-2._d*vpqrs(i,c,j,b))
!                term=term+vpqrs(a,j,c,i)*(-2._d*vpqrs(i,b,j,c)+4._d*vpqrs(i,c,j,b))
!                term=term*dpl(b,apr)
!                term=term/DA
!                term=DA*term
                
!                D2_1_1_ph_ph=D2_1_1_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_1_ph_ph=+0.5_d*D2_1_1_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_1_ph_ph


  !!$Second order contribution D2_1_2_ak,a'k'. The condition that k'=k is
!!$checked in the calling procedure. (hc)

!  real(d) function D2_1_2_ph_ph(a,apr)

!    integer, intent(in) :: a,apr

!    integer :: b,c,i,j, nsym1,nsym2,nsym3,u,v,r,s,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_2_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do c=nOcc+1,nBas
!       s=roccnum(c)
!       do i=1,nOcc
!          u=roccnum(i)
!          do j=1,nOcc
!             v=roccnum(j)
!
!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(apr)+e(c)-e(i)-e(j)
!                ebcij=e(b)+e(c)-e(i)-e(j)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(i,apr,j,c)*(+4._d*vpqrs(b,i,c,j)-2._d*vpqrs(b,j,c,i))
!                term=term+vpqrs(i,c,j,apr)*(-2._d*vpqrs(b,i,c,j)+4._d*vpqrs(b,j,c,i))
!                term=term*dpl(b,a)
!                term=term/DA
!                term=DA*term
                
!                D2_1_2_ph_ph=D2_1_2_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_2_ph_ph=+0.5_d*D2_1_2_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_2_ph_ph

 !!$Second order contribution D2_1_3_ak,a'k'. NO DELTA FUNCTIONS

!    real(d) function D2_1_3_ph_ph(a,apr,k,kpr)
!
!    integer, intent(in) :: a,apr,k,kpr

!    integer :: b,c,i, nsym1,nsym2,nsym3,u,v,r,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_3_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do c=nOcc+1,nBas
!       s=roccnum(c)
!       do i=1,nOcc
!          u=roccnum(i)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(b)+e(c)-e(kpr)-e(i)
!                ebcij=e(b)+e(c)-e(k)-e(i)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(b,k,c,i)*(+4._d*vpqrs(kpr,b,i,c)-2._d*vpqrs(kpr,c,i,b))
!                term=term+vpqrs(b,i,c,k)*(-2._d*vpqrs(kpr,b,i,c)+4._d*vpqrs(kpr,c,i,b))
!                term=term*dpl(a,apr)
!                term=term/DA
!                term=DA*term
                
!                D2_1_3_ph_ph=D2_1_3_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_3_ph_ph=+0.5_d*D2_1_3_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_3_ph_ph

 !!$Second order contribution D2_1_4_ak,a'k'. NO DELTA FUNCTIONS(hc)

!    real(d) function D2_1_4_ph_ph(a,apr,k,kpr)

!    integer, intent(in) :: a,apr,k,kpr

!    integer :: b,c,i, nsym1,nsym2,nsym3,u,v,r,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_4_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do c=nOcc+1,nBas
!       s=roccnum(c)
!       do i=1,nOcc
!          u=roccnum(i)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(b)+e(c)-e(kpr)-e(i)
!                ebcij=e(b)+e(c)-e(k)-e(i)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(b,k,c,i)*(+4._d*vpqrs(kpr,b,i,c)-2._d*vpqrs(kpr,c,i,b))
!                term=term+vpqrs(b,i,c,k)*(-2._d*vpqrs(kpr,b,i,c)+4._d*vpqrs(kpr,c,i,b))
!                term=term*dpl(apr,a)
!                term=term/DA
!                term=DA*term
                
!                D2_1_4_ph_ph=D2_1_4_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_4_ph_ph=+0.5_d*D2_1_4_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_4_ph_ph

 !!$Second order contribution D2_1_5_ak,a'k'. NO DELTA FUNCTIONS(hc)

!    real(d) function D2_1_5_ph_ph(a,apr,k,kpr)

!    integer, intent(in) :: a,apr,k,kpr

!    integer :: b,c,i, nsym1,nsym2,nsym3,u,v,r,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_5_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do c=nOcc+1,nBas
!       s=roccnum(c)
!       do i=1,nOcc
!          u=roccnum(i)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(b)+e(c)-e(kpr)-e(i)
!                ebcij=e(a)+e(c)-e(k)-e(i)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(a,k,c,i)*(+8._d*vpqrs(kpr,b,i,c)-4._d*vpqrs(kpr,c,i,b))
!                term=term+vpqrs(a,i,c,k)*(-4._d*vpqrs(kpr,b,i,c)+2._d*vpqrs(kpr,c,i,b))
!                term=term*dpl(b,apr)
!                term=term/DA
!                term=DA*term
                
!                D2_1_5_ph_ph=D2_1_5_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_5_ph_ph=+0.5_d*D2_1_5_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_5_ph_ph


 !!$Second order contribution D2_1_6_ak,a'k'. NO DELTA FUNCTIONS(hc)

 !   real(d) function D2_1_6_ph_ph(a,apr,k,kpr)

 !   integer, intent(in) :: a,apr,k,kpr

!    integer :: b,c,i, nsym1,nsym2,nsym3,u,v,r,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_6_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do c=nOcc+1,nBas
!       s=roccnum(c)
!       do i=1,nOcc
!          u=roccnum(i)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(apr)+e(c)-e(kpr)-e(i)
!                ebcij=e(b)+e(c)-e(k)-e(i)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(kpr,apr,i,c)*(+8._d*vpqrs(b,k,c,i)-4._d*vpqrs(b,i,c,k))
!                term=term+vpqrs(kpr,c,i,apr)*(-4._d*vpqrs(b,k,c,i)+2._d*vpqrs(b,i,c,k))
!                term=term*dpl(b,a)
!                term=term/DA
!                term=DA*term
                
!                D2_1_6_ph_ph=D2_1_6_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_6_ph_ph=+0.5_d*D2_1_6_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_6_ph_ph

 !!$Second order contribution D2_1_7_ak,a'k'. NO DELTA FUNCTIONS(hc)

!    real(d) function D2_1_7_ph_ph(a,apr,k,kpr)

!    integer, intent(in) :: a,apr,k,kpr

!    integer :: b,i,j, nsym1,nsym2,nsym3,u,v,r,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_7_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do i=1,nOcc
!       u=roccnum(i)
!       do j=1,nOcc
!          v=roccnum(j)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(apr)+e(c)-e(i)-e(j)
!                ebcij=e(a)+e(c)-e(i)-e(j)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(a,i,c,j)*(+4._d*vpqrs(i,apr,j,c)-2._d*vpqrs(i,c,j,apr))
!                term=term+vpqrs(a,j,c,i)*(-2._d*vpqrs(i,apr,j,c)+4._d*vpqrs(i,c,j,apr))
!                term=term*dpl(kpr,k)
!                term=term/DA
!                term=DA*term
                
!                D2_1_7_ph_ph=D2_1_7_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_7_ph_ph=+0.5_d*D2_1_7_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_7_ph_ph

 !!$Second order contribution D2_1_8_ak,a'k'. NO DELTA FUNCTIONS(hc)

!    real(d) function D2_1_8_ph_ph(a,apr,k,kpr)

!    integer, intent(in) :: a,apr,k,kpr

!    integer :: b,i,j, nsym1,nsym2,nsym3,u,v,r,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_7_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do i=1,nOcc
!       u=roccnum(i)
!       do j=1,nOcc
!          v=roccnum(j)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(apr)+e(c)-e(i)-e(j)
!                ebcij=e(a)+e(c)-e(i)-e(j)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(a,i,c,j)*(+4._d*vpqrs(i,apr,j,c)-2._d*vpqrs(i,c,j,apr))
!                term=term+vpqrs(a,j,c,i)*(-2._d*vpqrs(i,apr,j,c)+4._d*vpqrs(i,c,j,apr))
!                term=term*dpl(k,kpr)
!                term=term/DA
!                term=DA*term
                
!                D2_1_8_ph_ph=D2_1_8_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_8_ph_ph=+0.5_d*D2_1_8_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_8_ph_ph

 !!$Second order contribution D2_1_9_ak,a'k'. The condition that a=a' is
!!$checked in the calling procedure. (hc)

!  real(d) function D2_1_9_ph_ph(k,kpr)

!    integer, intent(in) :: k,kpr

!    integer :: b,c,i,j, nsym1,nsym2,nsym3,u,v,r,s,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_9_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do c=nOcc+1,nBas
!       s=roccnum(c)
!       do i=1,nOcc
!          u=roccnum(i)
!          do j=1,nOcc
!             v=roccnum(j)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(b)+e(c)-e(k)-e(i)
!                ebcij=e(b)+e(c)-e(j)-e(i)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(b,k,c,i)*(+4._d*vpqrs(j,b,i,c)-2._d*vpqrs(j,c,i,b))
!                term=term+vpqrs(b,i,c,k)*(-2._d*vpqrs(j,b,i,c)+4._d*vpqrs(j,c,i,b))
!                term=term*dpl(kpr,j)
!                term=term/DA
!                term=DA*term
                
!                D2_1_9_ph_ph=D2_1_9_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_9_ph_ph=+0.5_d*D2_1_9_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_9_ph_ph


!!$Second order contribution D2_1_10_ak,a'k'. The condition that a'=a is
!!$checked in the calling procedure.


!  real(d) function D2_1_10_ph_ph(k,kpr)

!    integer, intent(in) :: k,kpr

!    integer :: b,c,i,j, nsym1,nsym2,nsym3,u,v,r,s,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_10_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do c=nOcc+1,nBas
!       s=roccnum(c)
!       do i=1,nOcc
!          u=roccnum(i)
!          do j=1,nOcc
!             v=roccnum(j)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(b)+e(c)-e(j)-e(i)
!                ebcij=e(b)+e(c)-e(kpr)-e(i)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(kpr,b,i,c)*(+4._d*vpqrs(b,j,c,i)-2._d*vpqrs(b,i,c,j))
!                term=term+vpqrs(kpr,c,i,b)*(-2._d*vpqrs(b,j,c,i)+4._d*vpqrs(b,i,c,j))
!                term=term*dpl(k,j)
!                term=term/DA
!                term=DA*term
                
!                D2_1_10_ph_ph=D2_1_10_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_10_ph_ph=+0.5_d*D2_1_10_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_10_ph_ph

 !!$Second order contribution D2_1_11_ak,a'k'. NO DELTA FUNCTIONS

!   real(d) function D2_1_11_ph_ph(a,apr,k,kpr)

!    integer, intent(in) :: a,apr,k,kpr

!    integer :: b,i,j, nsym1,nsym2,nsym3,u,v,r,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_11_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do i=1,nOcc
!       u=roccnum(i)
!       do j=1,nOcc
!          v=roccnum(j)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(apr)+e(b)-e(j)-e(i)
!                ebcij=e(a)+e(b)-e(k)-e(i)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(a,k,b,i)*(+8._d*vpqrs(j,apr,i,b)-4._d*vpqrs(j,b,i,apr))
!                term=term+vpqrs(a,i,b,k)*(-4._d*vpqrs(j,apr,i,b)+2._d*vpqrs(j,b,i,apr))
!                term=term*dpl(kpr,j)
!                term=term/DA
!                term=DA*term
                
!                D2_1_11_ph_ph=D2_1_11_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_11_ph_ph=+0.5_d*D2_1_11_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_11_ph_ph


 !!$Second order contribution D2_1_12_ak,a'k'. NO DELTA FUNCTIONS (hc)

!    real(d) function D2_1_12_ph_ph(a,apr,k,kpr)

!    integer, intent(in) :: a,apr,k,kpr

!    integer :: b,i,j, nsym1,nsym2,nsym3,u,v,r,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_12_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!    do i=1,nOcc
 !      u=roccnum(i)
!       do j=1,nOcc
!          v=roccnum(j)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(apr)+e(b)-e(kpr)-e(i)
!                ebcij=e(a)+e(b)-e(j)-e(i)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(kpr,apr,i,b)*(+8._d*vpqrs(a,j,b,i)-4._d*vpqrs(a,i,b,j))
!                term=term+vpqrs(kpr,b,i,apr)*(-4._d*vpqrs(a,j,b,i)+2._d*vpqrs(a,i,b,j))
!                term=term*dpl(k,j)
!                term=term/DA
!                term=DA*term
                
!                D2_1_12_ph_ph=D2_1_12_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_12_ph_ph=+0.5_d*D2_1_12_ph_ph   ! NORMALIZATION OF WAVE FUNCTION             

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_12_ph_ph

!!$Second order contribution D2_1_13_ak,a'k'. 

!  real(d) function D2_1_13_ph_ph(a,apr,k,kpr)

!    integer, intent(in) :: a,apr,k,kpr

!    integer :: b,i,j, nsym1,nsym2,nsym3,u,v,r,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_13_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!       do i=1,nOcc
!          u=roccnum(i)
!          do j=1,nOcc
!             v=roccnum(j)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(a)+e(b)-e(i)-e(j)
!                ebcij=e(apr)+e(b)-e(i)-e(j)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(i,apr,j,b)*(+4._d*vpqrs(a,i,b,j)-2._d*vpqrs(a,j,b,i))
!                term=term+vpqrs(i,b,j,apr)*(-2._d*vpqrs(a,i,b,j)+4._d*vpqrs(a,j,b,i))
!                term=term*dpl(k,kpr)
!                term=term/DA
!                term=DA*term
                
!                D2_1_13_ph_ph=D2_1_13_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_13_ph_ph=+0.5_d*D2_1_13_ph_ph   ! NORMALIZATION OF WAVE FUNCTION 

!    D2_1_13_ph_ph=+1._d*D2_1_13_ph_ph    ! EXPRESSION FACTOR            

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_13_ph_ph
 
!!$Second order contribution D2_1_14_ak,a'k'. 

!  real(d) function D2_1_14_ph_ph(a,apr,k,kpr)

!    integer, intent(in) :: a,apr,k,kpr

!    integer :: b,c,j, nsym1,nsym2,nsym3,u,r,t,cnt
!    real(d) :: DA,eacij,ebcij,term

!    D2_1_14_ph_ph=0._d
!    cnt=0
!  do b=nOcc+1,nBas
!     r=roccnum(b)
!       do c=nOcc+1,nBas
!          t=roccnum(i)
!          do j=1,nOcc
!             u=roccnum(j)

!             nsym1=MT(orbSym(u),orbSym(v))
!             nsym2=MT(orbSym(a),orbSym(r))
!             nsym3=MT(orbSym(a1),orbSym(r))
             
             
!             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym1,nsym3) .eq. 1)) then
                
!                cnt=cnt+1
!!$                write(6,*) a,a1,orbSym(a),orbSym(a1)
!!$                write(6,*) u,v,r,orbSym(u),orbSym(v),orbSym(r)
!                term=0._d

!                eacij=e(b)+e(c)-e(k)-e(j)
!                ebcij=e(b)+e(c)-e(kpr)-e(j)
!                DA=eacij*ebcij
!                DA=(eijc-0.5_d*(e(a)+e(a1)))/((eijc-e(a))*(eijc-e(a1)))

!                term=term+vpqrs(a,u,r,v)*(2._d*vpqrs(u,a1,v,r)-vpqrs(u,r,v,a1))
!                term=term+vpqrs(a,v,r,u)*(2._d*vpqrs(u,r,v,a1)-vpqrs(u,a1,v,r)) 
!!$
!                term=term+vpqrs(kpr,b,j,c)*(+4._d*vpqrs(b,k,c,j)-2._d*vpqrs(b,j,c,k))
!                term=term+vpqrs(kpr,c,j,b)*(-2._d*vpqrs(b,k,c,j)+4._d*vpqrs(b,j,c,k))
!                term=term*dpl(a,apr)
!                term=term/DA
!                term=DA*term
                
!                D2_1_14_ph_ph=D2_1_14_ph_ph+term
                
!             end if
!          end do
!       end do
!    end do 
    
!    D2_1_14_ph_ph=+0.5_d*D2_1_14_ph_ph   ! NORMALIZATION OF WAVE FUNCTION 

!    D2_1_14_ph_ph=+1._d*D2_1_14_ph_ph    ! EXPRESSION FACTOR            

! FACTOR FOR THE ENTIRE EXPRESSION?

!  end function D2_1_14_ph_ph

!####################################################################### 
  
  function tau_D2_6_1(a,k,kpr,b1) result(func)

    implicit none

    real(d)             :: func
    integer, intent(in) :: a,k,kpr,b1
    integer             :: b,c,i,sym,nsym1,nsym2,nsym3,nsym4,c1,i1,cnt
    real*8              :: DA,eacij,ebcij,term

    func=0.0d0

    b=roccnum(b1)

    cnt=0
    
    do c1=nOcc+1,nBas
       c=roccnum(c1)
       do i1=1,nOcc
          i=roccnum(i1)
             
          nsym1=MT(orbSym(a),orbSym(c))
          nsym2=MT(orbSym(i),orbSym(k))
          nsym3=MT(orbSym(b),orbSym(c))
          nsym4=MT(orbSym(i),orbSym(kpr))             
                
          if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
             cnt=cnt+1
                   
             term=0.0

             eacij=e(a)+e(c)-e(k)-e(i)
             ebcij=e(b)+e(c)-e(kpr)-e(i)
             DA=eacij*ebcij
                   
             term=term+vpqrs(kpr,b,i,c)*(+8.0*vpqrs(a,k,c,i)-4.0*vpqrs(a,i,c,k))
             term=term+vpqrs(kpr,c,i,b)*(-4.0*vpqrs(a,k,c,i)+2.0*vpqrs(a,i,c,k))
!             term=term*dpl(b,apr)
             term=term/DA

             func=func+term
                
          end if
       end do
    end do
    
    return

  end function tau_D2_6_1

!####################################################################### 

      function tau_D2_6_2(apr,k,kpr,b1) result(func)

        implicit none

        integer, intent(in) :: apr,k,kpr,b1
        integer             :: b,c,i,sym,nsym1,nsym2,nsym3,nsym4,c1,i1,cnt
        real*8              :: DA,eacij,ebcij,term,func

        func=0.0d0
        cnt=0

        b=roccnum(b1)
          
        do c1=nOcc+1,nBas
           c=roccnum(c1)
           do i1=1,nOcc
              i=roccnum(i1)
                
              nsym1=MT(orbSym(apr),orbSym(c))
              nsym2=MT(orbSym(i),orbSym(kpr))
              nsym3=MT(orbSym(b),orbSym(c))
              nsym4=MT(orbSym(i),orbSym(k))
                
              if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                   
                 cnt=cnt+1
                   
                 term=0.0

                 eacij=e(b)+e(c)-e(k)-e(i)
                 ebcij=e(apr)+e(c)-e(kpr)-e(i)
                 DA=eacij*ebcij

                 term=term+vpqrs(b,k,c,i)*(+8.0*vpqrs(kpr,apr,i,c)-4.0*vpqrs(kpr,c,i,apr))
                 term=term+vpqrs(b,i,c,k)*(-4.0*vpqrs(kpr,apr,i,c)+2.0*vpqrs(kpr,c,i,apr))
!                 term=term*dpl(b,a)
                 term=term/DA
                 
                 func=func+term
                
              end if
           end do
        end do
       
        return
        
      end function tau_D2_6_2

!####################################################################### 

      function tau_D2_6_3(a,apr,k,j1) result(func)

        implicit none

        integer, intent(in) :: a,apr,k,j1
        integer             :: b,b1,i,j,sym,nsym1,nsym2,nsym3,nsym4,i1,cnt
        real*8              :: DA,eacij,ebcij,term, func

        func=0.0d0

        cnt=0

        j=roccnum(j1)
  
        do b1=nOcc+1,nBas
           b=roccnum(b1)
           do i1=1,nOcc
              i=roccnum(i1)

              nsym1=MT(orbSym(a),orbSym(b))
              nsym2=MT(orbSym(i),orbSym(k))
              nsym3=MT(orbSym(apr),orbSym(b))
              nsym4=MT(orbSym(i),orbSym(j))
             
              if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
                 cnt=cnt+1

                 term=0.0

                 eacij=e(apr)+e(b)-e(j)-e(i)
                 ebcij=e(a)+e(b)-e(k)-e(i)
                 DA=eacij*ebcij

                 term=term+vpqrs(j,apr,i,b)*(+8.0*vpqrs(a,k,b,i)-4.0*vpqrs(a,i,b,k))
                 term=term+vpqrs(j,b,i,apr)*(-4.0*vpqrs(a,k,b,i)+2.0*vpqrs(a,i,b,k))
!                 term=term*dpl(kpr,j)
                 term=term/DA
                
                 func=func+term
                
              end if
       
           end do
        end do

        return

      end function tau_D2_6_3

!####################################################################### 

      function tau_D2_6_4(a,apr,kpr,j1) result(func)

        implicit none

        integer, intent(in) :: a,apr,kpr,j1
        integer             :: b,i,j,sym,nsym1,nsym2,nsym3,nsym4,b1,i1,cnt
        real*8              :: DA,eacij,ebcij,term,func


        func=0.0d0

        cnt=0

        j=roccnum(j1)
   
        do b1=nOcc+1,nBas
           b=roccnum(b1)
           do i1=1,nOcc
              i=roccnum(i1)
 
              nsym1=MT(orbSym(apr),orbSym(b))
              nsym2=MT(orbSym(i),orbSym(kpr))
              nsym3=MT(orbSym(a),orbSym(b))
              nsym4=MT(orbSym(i),orbSym(j))            
             
              if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then
                
                 cnt=cnt+1

                 term=0.0

                 eacij=e(apr)+e(b)-e(kpr)-e(i)
                 ebcij=e(a)+e(b)-e(j)-e(i)
                 DA=eacij*ebcij

                 term=term+vpqrs(a,j,b,i)*(+8.0*vpqrs(kpr,apr,i,b)-4.0*vpqrs(kpr,b,i,apr))
                 term=term+vpqrs(a,i,b,j)*(-4.0*vpqrs(kpr,apr,i,b)+2.0*vpqrs(kpr,b,i,apr))
!                 term=term*dpl(k,j)
                 term=term/DA
                
                 func=func+term
                
             end if
       
          end do
       end do

        return

      end function tau_D2_6_4

!####################################################################### 

  end module D_matrix
