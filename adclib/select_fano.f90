module select_fano
  
  use constants
  use parameters
  use adc_ph
  use misc
  
  implicit none
  

contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! KEY SUBROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine select_atom_is(kpq)
  
    implicit none
  
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    integer               :: i,ap,n2pne
    integer               :: isym,a,cah,ic
    integer, dimension(7) :: col
    
    kpq(1,0)=0

    if (lifrzcore) then
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do i=1,hcentre(0)
             cah=hcentre(i)
             call iscore(cah,ic)
             if (ic.eq.0) then
                isym=MT(orbSym(cah),orbSym(a))
                if(isym .eq. nirrep) then
                   kpq(1,0)=kpq(1,0)+1
                   call fill_indices(col(:),1,1,a,-1,cah,-1,0)
                   kpq(:,kpq(1,0))=col(:)
                end if
             endif
          end do
       end do
    else
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do i=1,hcentre(0)
             cah=hcentre(i)
             isym=MT(orbSym(cah),orbSym(a))
             if(isym .eq. nirrep) then
                kpq(1,0)=kpq(1,0)+1
                call fill_indices(col(:),1,1,a,-1,cah,-1,0)
                kpq(:,kpq(1,0))=col(:)
             end if
          end do
       end do
    endif

  end subroutine select_atom_is

!#######################################################################

  subroutine select_atom_is_cvs(kpq)

    implicit none

    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq

    integer :: i,ap,n2pne,ic
    integer :: isym,a,cah
    integer, dimension(7) :: col

    kpq(1,0)=0

    do ap=nOcc+1,nBas
       a=roccnum(ap)
       do i=1,hcentre(0)
          cah=hcentre(i)
          call iscore(cah,ic)
          if (ic.eq.1) then
             isym=MT(orbSym(cah),orbSym(a))
             if(isym .eq. nirrep) then
                kpq(1,0)=kpq(1,0)+1
                call fill_indices(col(:),1,1,a,-1,cah,-1,0)
                kpq(:,kpq(1,0))=col(:)
             end if
          endif
       enddo
    enddo

    return

  end subroutine select_atom_is_cvs

!#######################################################################

  subroutine iscore(cah,ic)

    implicit none

    integer :: cah,k,ic

    ic=0

    do k=1,ncore       
       if (cah.eq.icore(k)) ic=1
    enddo

    return

  end subroutine iscore

!#######################################################################

  subroutine select_atom_isf(kpq)
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    
    integer :: i,ap,n2pne
    integer :: isym,isym2,a,cah
    integer, dimension(7) :: col
    
    
100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))
    
!    write(ilog,*) 'Selecting FINAL  1h1p subspace'
    
    kpq(1,0)=0

    do ap=nOcc+1,nBas
       a=roccnum(ap)
       do i=1,hcentre(0)
          cah=hcentre(i)
          isym=MT(orbSym(cah),orbSym(a))
          if(nirrep .eq. nirrep2) then
             
             if(isym .eq. 1) then
                kpq(1,0)=kpq(1,0)+1
                call fill_indices(col(:),1,1,a,-1,cah,-1,0)
                kpq(:,kpq(1,0))=col(:)
             end if

          else

             isym2=MT(nirrep,nirrep2)
             if(isym .eq. isym2) then
                kpq(1,0)=kpq(1,0)+1
                call fill_indices(col(:),1,1,a,-1,cah,-1,0)
                kpq(:,kpq(1,0))=col(:)
                
             end if
          end if
       end do
    end do
       
!    write(ilog,100) "Number of 1h-1p FINAL configs in the IS ADC",kpq(1,0) 
!    write(ilog,103)
!    write(ilog,*) "FINAL Singly excited configurations allowed in the Fin. St. Manif."
!    write(ilog,101) "CNF","SPN","HL1","PT1"
!    do i=1,kpq(1,0)
!       write(ilog,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
!    enddo
       
        
  end subroutine select_atom_isf

  subroutine select_atom_ist(kpq)
    
    integer, dimension(7,0:2*nBas**2*nOcc**2), intent(inout) :: kpq
    
    integer :: i,ap,n2pne
    integer :: isym,isym2,a,cah
    integer, dimension(7) :: col
    
    
100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))
    
!    write(ilog,'(a)') 'Selecting TOTAL 1h1p subspace'
    
    kpq(1,0)=0

    do ap=nOcc+1,nBas
       a=roccnum(ap)
       do i=1,hcentre(0)
           cah=hcentre(i)
          isym=MT(orbSym(cah),orbSym(a))
   if(nirrep .eq. nirrep2) then 
          if((isym .eq. nirrep).or.(isym .eq. 1)) then
             kpq(1,0)=kpq(1,0)+1
             call fill_indices(col(:),1,1,a,-1,cah,-1,0)
             kpq(:,kpq(1,0))=col(:)
          end if
   else
          
          isym2=MT(nirrep,nirrep2)
          if((isym .eq. nirrep).or.(isym .eq. isym2)) then
             kpq(1,0)=kpq(1,0)+1
             call fill_indices(col(:),1,1,a,-1,cah,-1,0)
             kpq(:,kpq(1,0))=col(:)
          end if
   end if
       end do
    end do
       
!    write(ilog,100) "Number of 1h-1p TOTAL configs in the IS ADC",kpq(1,0) 
!    write(ilog,103)
!    write(ilog,*) "TOTAL Singly excited configurations allowed in the Fin. St. Manif."
!    write(ilog,101) "CNF","SPN","HL1","PT1"
!    do i=1,kpq(1,0)
!       write(ilog,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
!    enddo
       
        
  end subroutine select_atom_ist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    subroutine select_atom_d(kpq,flag)

!!$Includes 2p2pnxnx configurations in the both, initial and final states.
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    integer, intent(in) :: flag

    integer :: i,j,a,b,k
    integer :: ih,jh,ap,bp
    integer :: cnti,cntf
    integer :: isym1, isym2,ic
    integer, dimension(7) :: col
    real(d) :: einit,ei,ej

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I5,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-")) 

    einit=abs(e(hinit))

! Filling the kpq arrays in the order: a=b,i=j;a|=b,i=j;a=b,i|=j;a|=b,i|=j(I,II).

    if (flag .eq. -1) then
       
       kpq(2:5,0)=0
       cnti=kpq(1,0)
    
!!$a=b,i=j i=hcentre(ih)

       if(nirrep .eq. 1) then
          do ih=1,hcentre(0)
             i=hcentre(ih)
             ei=abs(e(i))

             call iscore(i,ic)
             if (lifrzcore.and.ic.eq.1) cycle

             do ap=nOcc+1,nBas
                a=roccnum(ap)
                if(einit .le. 2._d*ei) then
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                else
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                end if
             end do
          end do
       end if

!!$a|=b,i=j
  
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do bp=ap+1,nBas
             b=roccnum(bp)
             isym1=MT(orbSym(a),orbSym(b))
             if(isym1 .eq. nirrep) then
                do ih=1,hcentre(0)
                   i=hcentre(ih)
                   ei=abs(e(i))

                   call iscore(i,ic)
                   if (lifrzcore.and.ic.eq.1) cycle

                   if(einit .le. 2._d*ei) then
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   else
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   end if
                end do
             end if
          end do
       end do

!a=b,i|=j
  
       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))

          call iscore(i,ic)
          if (lifrzcore.and.ic.eq.1) cycle

          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))

             call iscore(j,ic)
             if (lifrzcore.and.ic.eq.1) cycle

             isym1=MT(orbSym(i),orbSym(j))
             if (isym1 .eq. nirrep) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                      cnti=cnti+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cnti)=col(:)
                end do
             end if
          end do
       end do

!a|=b,i|=j spin I

       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))

          call iscore(i,ic)
          if (lifrzcore.and.ic.eq.1) cycle

          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))

             call iscore(j,ic)
             if (lifrzcore.and.ic.eq.1) cycle

             isym1=MT(orbSym(i),orbSym(j))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                do bp=ap+1,nBas
                   b=roccnum(bp)
                   isym2=MT(orbSym(a),orbSym(b))
                   if(MT(isym1,isym2) .eq. nirrep) then 
                      if(einit .le. (ei+ej)) then
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      else
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      end if
                   end if
                end do
             end do
          end do
       end do
    
       kpq(:,cnti+1:cnti+kpq(5,0))=kpq(:,cnti+1-kpq(5,0):cnti)
       kpq(2,cnti+1:cnti+kpq(5,0))=12
       kpq(7,cnti+1:cnti+kpq(5,0))=5

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

    elseif(flag .eq. 1) then

       kpq(2:5,0)=0
       cntf=kpq(1,0)

       kpq(:,cntf+1:cntf+kpq(5,0))=kpq(:,cntf+1-kpq(5,0):cntf)
       kpq(2,cntf+1:cntf+kpq(5,0))=12 
       kpq(7,cntf+1:cntf+kpq(5,0))=5

    end if
    
  end subroutine select_atom_d

!#######################################################################

  subroutine select_atom_d_cvs(kpq,flag)

!!$Includes 2p2pnxnx configurations in the both, initial and final states.
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    integer, intent(in)                                    :: flag

    integer               :: i,j,a,b,k,ic1,ic2
    integer               :: ih,jh,ap,bp
    integer               :: cnti,cntf
    integer               :: isym1, isym2
    integer, dimension(7) :: col
    real(d)               :: einit,ei,ej

!-----------------------------------------------------------------------
! Filling the kpq arrays in the order:
!
! a=b,i=j;
! a|=b,i=j;
! a=b,i|=j;
! a|=b,i|=j(I,II).
!
! Note that in the CVS approximation, the following configurations are
! ignored: a=b,i=j
!          a|=b,i=j;
!
! i.e., we only need to consider a=b,i|=j and a|=b,i|=j(I,II) for which
!       one, and only one, of i or j corresponds to a core orbital
!-----------------------------------------------------------------------

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I5,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-")) 
    
    einit=abs(e(hinit))

    kpq(2:5,0)=0
    cnti=kpq(1,0)

!-----------------------------------------------------------------------
! a=b,i|=j, i or j a core orbital, but not both
!-----------------------------------------------------------------------
    do ih=1,hcentre(0)
       i=hcentre(ih)
       call iscore(i,ic1)
       ei=abs(e(i))

       do jh=ih+1,hcentre(0)
          j=hcentre(jh)
          call iscore(j,ic2)
          ej=abs(e(j))

          if ((ic1.eq.0.and.ic2.eq.1) &
               .or.(ic1.eq.1.and.ic2.eq.0)) then
             isym1=MT(orbSym(i),orbSym(j))
             if (isym1 .eq. nirrep) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                   cnti=cnti+1
                   kpq(4,0)=kpq(4,0)+1
                   call fill_indices(col(:),2,1,a,a,i,j,3) 
                   kpq(:,cnti)=col(:)
                end do
             end if
             
          endif

       end do
    end do

!-----------------------------------------------------------------------
! a|=b,i|=j spin I, i or j a core orbital, but not both
!-----------------------------------------------------------------------
       do ih=1,hcentre(0)
          i=hcentre(ih)          
          call iscore(i,ic1)
          ei=abs(e(i))
          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             call iscore(j,ic2)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))

             if ((ic1.eq.0.and.ic2.eq.1) &
                  .or.(ic1.eq.1.and.ic2.eq.0)) then
                
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                   do bp=ap+1,nBas
                      b=roccnum(bp)
                      isym2=MT(orbSym(a),orbSym(b))
                      if(MT(isym1,isym2) .eq. nirrep) then 
                         if(einit .le. (ei+ej)) then
                            cnti=cnti+1
                            kpq(5,0)=kpq(5,0)+1
                            call fill_indices(col(:),2,11,a,b,i,j,4)
                            kpq(:,cnti)=col(:)
                         else
                            cnti=cnti+1
                            kpq(5,0)=kpq(5,0)+1
                            call fill_indices(col(:),2,11,a,b,i,j,4)
                            kpq(:,cnti)=col(:)
                         end if
                      end if
                   end do
                end do
                
             endif

          end do
       end do

       kpq(:,cnti+1:cnti+kpq(5,0))=kpq(:,cnti+1-kpq(5,0):cnti)
       kpq(2,cnti+1:cnti+kpq(5,0))=12
       kpq(7,cnti+1:cnti+kpq(5,0))=5

    return

  end subroutine select_atom_d_cvs

!#######################################################################

    subroutine select_atom_df(kpq,flag)

!!$Includes 2p2pnxnx configurations in the both, initial and final states.
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    integer, intent(in) :: flag

    integer :: i,j,a,b,k
    integer :: ih,jh,ap,bp
    integer :: cnti,cntf
    integer :: isym1,isym2,isym3
    integer, dimension(7) :: col
    real(d) :: einit,ei,ej

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I5,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-")) 

    einit=abs(e(hinit))

! Filling the kpq arrays in the order: a=b,i=j;a|=b,i=j;a=b,i|=j;a|=b,i|=j(I,II).

    if (flag .eq. -1) then
       
!       write(ilog,'(a)') 'Selecting FINAL  2h2p subspace'
       
       kpq(2:5,0)=0
       cnti=kpq(1,0)
    
!!$a=b,i=j i=hcentre(ih)
       
       if(nirrep .eq. nirrep2) then
          do ih=1,hcentre(0)
             i=hcentre(ih)
             ei=abs(e(i))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                if(einit .le. 2._d*ei) then
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                else
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                end if
             end do
          end do
       end if

!        if(nirrep .eq. 1) then
!           do ih=1,hneighb(0)
!              i=hneighb(ih)
!              ei=abs(e(i))
!              do ap=nOcc+1,nBas
!                 a=roccnum(ap)
!                 if(einit .le. 2._d*ei) then
!                    cnti=cnti+1
!                    kpq(2,0)=kpq(2,0)+1
!                    call fill_indices(col(:),2,1,a,a,i,i,1)
!                    kpq(:,cnti)=col(:)
!                 else
!                    cnti=cnti+1
!                    kpq(2,0)=kpq(2,0)+1
!                    call fill_indices(col(:),2,1,a,a,i,i,1)
!                    kpq(:,cnti)=col(:)
!                 end if
!              end do
!           end do
!        end if

!!$a|=b,i=j
  
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do bp=ap+1,nBas
             b=roccnum(bp)
             isym1=MT(orbSym(a),orbSym(b))

     if(nirrep .eq. nirrep2) then
             if(isym1 .eq. 1) then
                do ih=1,hcentre(0)
                   i=hcentre(ih)
                   ei=abs(e(i))
                   if(einit .le. 2._d*ei) then
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   else
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   end if
                end do
             end if
      else
             isym2=MT(nirrep,nirrep2)
             if(isym1 .eq. isym2) then
                do ih=1,hcentre(0)
                   i=hcentre(ih)
                   ei=abs(e(i))
                   if(einit .le. 2._d*ei) then
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   else
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   end if
                end do
             end if
       end if
 
          end do
       end do


!        do ap=nOcc+1,nBas
!           a=roccnum(ap)
!           do bp=ap+1,nBas
!              b=roccnum(bp)
!              isym1=MT(orbSym(a),orbSym(b))
!              if(isym1 .eq. nirrep) then
!                 do ih=1,hneighb(0)
!                    i=hneighb(ih)
!                    ei=abs(e(i))
!                    if(einit .le. 2._d*ei) then
!                       cnti=cnti+1
!                       kpq(3,0)=kpq(3,0)+1
!                       call fill_indices(col(:),2,1,a,b,i,i,2)
!                       kpq(:,cnti)=col(:)
!                    else
!                       cnti=cnti+1
!                       kpq(3,0)=kpq(3,0)+1
!                       call fill_indices(col(:),2,1,a,b,i,i,2)
!                       kpq(:,cnti)=col(:)
!                    end if
!                 end do
!              end if
!           end do
!        end do
!a=b,i|=j
  
       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))
          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))

     if(nirrep .eq. nirrep2) then
             if (isym1 .eq. 1) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                      cnti=cnti+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cnti)=col(:)
                end do
             end if
      else 
             isym2=MT(nirrep,nirrep2)
             if (isym1 .eq. isym2) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                      cnti=cnti+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cnti)=col(:)
                end do
             end if
       end if

          end do
       end do

!        do ih=1,hneighb(0)
!          i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=1,hcentre(0)
!              j=hcentre(jh)
!             ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              if (isym1 .eq. nirrep) then
!                 do ap=nOcc+1,nBas
!                    a=roccnum(ap)
!                       cnti=cnti+1
!                       kpq(4,0)=kpq(4,0)+1
!                       call fill_indices(col(:),2,1,a,a,i,j,3) 
!                       kpq(:,cnti)=col(:)
!                 end do
!              end if
!           end do
!        end do

!        do ih=1,hneighb(0)
!           i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=ih+1,hneighb(0)
!              j=hneighb(jh)
!              ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              if (isym1 .eq. nirrep) then
!                 do ap=nOcc+1,nBas
!                    a=roccnum(ap)
!                       cnti=cnti+1
!                       kpq(4,0)=kpq(4,0)+1
!                       call fill_indices(col(:),2,1,a,a,i,j,3) 
!                       kpq(:,cnti)=col(:)
!                 end do
!              end if
!           end do
!        end do
    
!a|=b,i|=j spin I

       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))
          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                do bp=ap+1,nBas
                   b=roccnum(bp)
                   isym2=MT(orbSym(a),orbSym(b))

           if(nirrep .eq. nirrep2) then
                   if(MT(isym1,isym2) .eq. 1) then 
                      if(einit .le. (ei+ej)) then
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      else
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      end if
                   end if
            else
                   isym3=MT(nirrep,nirrep2)
                   if(MT(isym1,isym2) .eq. isym3) then 
                      if(einit .le. (ei+ej)) then
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      else
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      end if
                   end if
            end if 

                end do
             end do
          end do
       end do

!        do ih=1,hneighb(0)
!           i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=1,hcentre(0)
!              j=hcentre(jh)
!              ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              do ap=nOcc+1,nBas
!                a=roccnum(ap)
!                 do bp=ap+1,nBas
!                    b=roccnum(bp)
!                    isym2=MT(orbSym(a),orbSym(b))
!                    if(MT(isym1,isym2) .eq. nirrep) then 
!                       if(einit .le. (ei+ej)) then
!                         cnti=cnti+1
!                          kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       else
!                          cnti=cnti+1
!                         kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       end if
!                    end if
!                 end do
!              end do
!           end do
!        end do

! aukommentier -> Ne sonst HF

!        do ih=1,hneighb(0)
!           i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=ih+1,hneighb(0)
!              j=hneighb(jh)
!              ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              do ap=nOcc+1,nBas
!                 a=roccnum(ap)
!                 do bp=ap+1,nBas
!                    b=roccnum(bp)
!                    isym2=MT(orbSym(a),orbSym(b))
!                    if(MT(isym1,isym2) .eq. nirrep) then 
!                       if(einit .le. (ei+ej)) then
!                          cnti=cnti+1
!                          kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       else
!                          cnti=cnti+1
!                          kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       end if
!                    end if
!                 end do
!              end do
!           end do
!        end do
    
       kpq(:,cnti+1:cnti+kpq(5,0))=kpq(:,cnti+1-kpq(5,0):cnti)
       kpq(2,cnti+1:cnti+kpq(5,0))=12
       kpq(7,cnti+1:cnti+kpq(5,0))=5

       
!       write(ilog,100) "Number of 2h-2p |abij> FINAL configs in the IS ADC", cnti+kpq(5,0)-kpq(1,0)
!       write(ilog,103)
!       write(ilog,*) " FINAL Doubly excited |abij> configs allowed in the Init. St. Manif."
!       write(ilog,101) "CNF","SPN","HL1","HL2","PT1","PT2"
!    
!       do k=kpq(1,0)+1,cnti+kpq(5,0) 
!          write(ilog,102) k,kpq(2,k),kpq(3,k),&
!               kpq(4,k),kpq(5,k),kpq(6,k)
!       end do

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
    elseif(flag .eq. 1) then
       
!       write(ilog,'(a)') 'Selecting final 2h2p subspace'
       kpq(2:5,0)=0
       cntf=kpq(1,0)

       kpq(:,cntf+1:cntf+kpq(5,0))=kpq(:,cntf+1-kpq(5,0):cntf)
       kpq(2,cntf+1:cntf+kpq(5,0))=12 
       kpq(7,cntf+1:cntf+kpq(5,0))=5
    
!       write(ilog,100) "Number of 2h-2p |abij> configs in the FS ADC", cntf+kpq(5,0)-kpq(1,0)
!       write(ilog,103)
!       write(ilog,*) " Doubly excited |abij> configs allowed in the Fin. St. Manif."
!       write(ilog,101) "CNF","SPN","HL1","HL2","PT1","PT2"
!       
!       do k=kpq(1,0)+1, cntf+kpq(5,0)
!          write(ilog,102) k,kpq(2,k),kpq(3,k),&
!               kpq(4,k),kpq(5,k),kpq(6,k)
!       end do
       
    end if
    
  end subroutine select_atom_df
!!$--------------------------------------------



    subroutine select_atom_dt(kpq,flag)

!!$Includes 2p2pnxnx configurations in the both, initial and final states.
    
    integer, dimension(7,0:2*nBas**2*nOcc**2), intent(inout) :: kpq
    integer, intent(in) :: flag

    integer :: i,j,a,b,k
    integer :: ih,jh,ap,bp
    integer :: cnti,cntf
    integer :: isym1,isym2,isym3
    integer, dimension(7) :: col
    real(d) :: einit,ei,ej

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I5,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-")) 

    einit=abs(e(hinit))

! Filling the kpq arrays in the order: a=b,i=j;a|=b,i=j;a=b,i|=j;a|=b,i|=j(I,II).

    if (flag .eq. -1) then
       
!       write(ilog,'(a)') 'Selecting TOTAL  2h2p subspace'
       
       kpq(2:5,0)=0
       cnti=kpq(1,0)
    
!!$a=b,i=j i=hcentre(ih)
       
       if( ( nirrep .eq. 1) .or. (nirrep .eq. nirrep2) ) then
          do ih=1,hcentre(0)
             i=hcentre(ih)
             ei=abs(e(i))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                if(einit .le. 2._d*ei) then
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                else
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                end if
             end do
          end do
       end if

!        if(nirrep .eq. 1) then
!           do ih=1,hneighb(0)
!              i=hneighb(ih)
!              ei=abs(e(i))
!              do ap=nOcc+1,nBas
!                 a=roccnum(ap)
!                 if(einit .le. 2._d*ei) then
!                    cnti=cnti+1
!                    kpq(2,0)=kpq(2,0)+1
!                    call fill_indices(col(:),2,1,a,a,i,i,1)
!                    kpq(:,cnti)=col(:)
!                 else
!                    cnti=cnti+1
!                    kpq(2,0)=kpq(2,0)+1
!                    call fill_indices(col(:),2,1,a,a,i,i,1)
!                    kpq(:,cnti)=col(:)
!                 end if
!              end do
!           end do
!        end if

!!$a|=b,i=j
  
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do bp=ap+1,nBas
             b=roccnum(bp)
             isym1=MT(orbSym(a),orbSym(b))

       if(nirrep .eq. nirrep2) then
             if((isym1 .eq. nirrep).or.(isym1 .eq. 1)) then
                do ih=1,hcentre(0)
                   i=hcentre(ih)
                   ei=abs(e(i))
                   if(einit .le. 2._d*ei) then
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   else
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   end if
                end do
             end if
        else 
              isym2=MT(nirrep,nirrep2)
              if((isym1 .eq. nirrep).or.(isym1 .eq. isym2)) then
                do ih=1,hcentre(0)
                   i=hcentre(ih)
                   ei=abs(e(i))
                   if(einit .le. 2._d*ei) then
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   else
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   end if
                end do
             end if
        end if

          end do
       end do


!        do ap=nOcc+1,nBas
!           a=roccnum(ap)
!           do bp=ap+1,nBas
!              b=roccnum(bp)
!              isym1=MT(orbSym(a),orbSym(b))
!              if(isym1 .eq. nirrep) then
!                 do ih=1,hneighb(0)
!                    i=hneighb(ih)
!                    ei=abs(e(i))
!                    if(einit .le. 2._d*ei) then
!                       cnti=cnti+1
!                       kpq(3,0)=kpq(3,0)+1
!                       call fill_indices(col(:),2,1,a,b,i,i,2)
!                       kpq(:,cnti)=col(:)
!                    else
!                       cnti=cnti+1
!                       kpq(3,0)=kpq(3,0)+1
!                       call fill_indices(col(:),2,1,a,b,i,i,2)
!                       kpq(:,cnti)=col(:)
!                    end if
!                 end do
!              end if
!           end do
!        end do
!a=b,i|=j
  
       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))
          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
       if(nirrep .eq. nirrep2) then
             if ((isym1 .eq. nirrep).or.(isym1 .eq. 1)) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                      cnti=cnti+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cnti)=col(:)
                end do
             end if
        else
             isym2=MT(nirrep,nirrep2)
             if ((isym1 .eq. nirrep).or.(isym1 .eq. isym2)) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                      cnti=cnti+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cnti)=col(:)
                end do
             end if
        end if 


          end do
       end do

!        do ih=1,hneighb(0)
!          i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=1,hcentre(0)
!              j=hcentre(jh)
!             ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              if (isym1 .eq. nirrep) then
!                 do ap=nOcc+1,nBas
!                    a=roccnum(ap)
!                       cnti=cnti+1
!                       kpq(4,0)=kpq(4,0)+1
!                       call fill_indices(col(:),2,1,a,a,i,j,3) 
!                       kpq(:,cnti)=col(:)
!                 end do
!              end if
!           end do
!        end do

!        do ih=1,hneighb(0)
!           i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=ih+1,hneighb(0)
!              j=hneighb(jh)
!              ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              if (isym1 .eq. nirrep) then
!                 do ap=nOcc+1,nBas
!                    a=roccnum(ap)
!                       cnti=cnti+1
!                       kpq(4,0)=kpq(4,0)+1
!                       call fill_indices(col(:),2,1,a,a,i,j,3) 
!                       kpq(:,cnti)=col(:)
!                 end do
!              end if
!           end do
!        end do
    
!a|=b,i|=j spin I

       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))
          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                do bp=ap+1,nBas
                   b=roccnum(bp)
                   isym2=MT(orbSym(a),orbSym(b))
        if(nirrep .eq. nirrep2) then
                   if((MT(isym1,isym2) .eq. nirrep).or. (MT(isym1,isym2) .eq. 1)) then 
                      if(einit .le. (ei+ej)) then
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      else
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      end if
                   end if
         else
                   isym3=MT(nirrep,nirrep2)
                   if((MT(isym1,isym2) .eq. nirrep).or. (MT(isym1,isym2) .eq. isym3)) then 
                      if(einit .le. (ei+ej)) then
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      else
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      end if
                   end if
          end if


                end do
             end do
          end do
       end do

!        do ih=1,hneighb(0)
!           i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=1,hcentre(0)
!              j=hcentre(jh)
!              ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              do ap=nOcc+1,nBas
!                a=roccnum(ap)
!                 do bp=ap+1,nBas
!                    b=roccnum(bp)
!                    isym2=MT(orbSym(a),orbSym(b))
!                    if(MT(isym1,isym2) .eq. nirrep) then 
!                       if(einit .le. (ei+ej)) then
!                         cnti=cnti+1
!                          kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       else
!                          cnti=cnti+1
!                         kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       end if
!                    end if
!                 end do
!              end do
!           end do
!        end do

! aukommentier -> Ne sonst HF

!        do ih=1,hneighb(0)
!           i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=ih+1,hneighb(0)
!              j=hneighb(jh)
!              ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              do ap=nOcc+1,nBas
!                 a=roccnum(ap)
!                 do bp=ap+1,nBas
!                    b=roccnum(bp)
!                    isym2=MT(orbSym(a),orbSym(b))
!                    if(MT(isym1,isym2) .eq. nirrep) then 
!                       if(einit .le. (ei+ej)) then
!                          cnti=cnti+1
!                          kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       else
!                          cnti=cnti+1
!                          kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       end if
!                    end if
!                 end do
!              end do
!           end do
!        end do
    
       kpq(:,cnti+1:cnti+kpq(5,0))=kpq(:,cnti+1-kpq(5,0):cnti)
       kpq(2,cnti+1:cnti+kpq(5,0))=12
       kpq(7,cnti+1:cnti+kpq(5,0))=5

       
!       write(ilog,100) "Number of 2h-2p |abij> TOTAL configs in the IS ADC", cnti+kpq(5,0)-kpq(1,0)
!       write(ilog,103)
!       write(ilog,*) " TOTAL Doubly excited |abij> configs allowed in the Init. St. Manif."
!       write(ilog,101) "CNF","SPN","HL1","HL2","PT1","PT2"
!    
!       do k=kpq(1,0)+1,cnti+kpq(5,0) 
!          write(ilog,102) k,kpq(2,k),kpq(3,k),&
!               kpq(4,k),kpq(5,k),kpq(6,k)
!       end do

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
    elseif(flag .eq. 1) then
       
!       write(ilog,'(a)') 'Selecting final 2h2p subspace'
       kpq(2:5,0)=0
       cntf=kpq(1,0)

       kpq(:,cntf+1:cntf+kpq(5,0))=kpq(:,cntf+1-kpq(5,0):cntf)
       kpq(2,cntf+1:cntf+kpq(5,0))=12 
       kpq(7,cntf+1:cntf+kpq(5,0))=5
    
!       write(ilog,100) "Number of 2h-2p |abij> configs in the FS ADC", cntf+kpq(5,0)-kpq(1,0)
!       write(ilog,103)
!       write(ilog,*) " Doubly excited |abij> configs allowed in the Fin. St. Manif."
!       write(ilog,101) "CNF","SPN","HL1","HL2","PT1","PT2"
!       
!       do k=kpq(1,0)+1, cntf+kpq(5,0)
!          write(ilog,102) k,kpq(2,k),kpq(3,k),&
!               kpq(4,k),kpq(5,k),kpq(6,k)
!       end do
       
    end if
    
  end subroutine select_atom_dt
!!$--------------------------------------------

  subroutine get_ncnfi_ryd(kpq,ncnfi,nryd)

!!$ Subroutine finds the maximum number of Rydberg states in a given basis.
!!$ It is used to identify the Rydberg series of Fano states (initial states).

    integer, dimension(nBas), intent(out) :: ncnfi
    integer, dimension(7,0:nBas**2*4*nOcc**2),intent(in) :: kpq
    integer, intent(out) :: nryd

    integer :: i,j,a,cnt

    cnt=0

    do i=1,kpq(1,0)
       j=kpq(3,i)
       if((j .eq. hinit)) then
          cnt=cnt+1
          ncnfi(cnt)=i
       end if
    end do
    nryd=cnt

    write(ilog,'(a)') 'Number of rydberg states', nryd
    write(ilog,'(a)') 'Cnfs contributing to the Rydbergs', ncnfi(1:nryd)

  end subroutine get_ncnfi_ryd


!!$---------------------------------------------
!!$---------------------------------------------

  subroutine select_fstate_ryd(ndim,stsel,vec,nstate,nisri)

    integer, intent(in) :: ndim
    real(d), intent(in) :: stsel
    real(d), dimension(ndim),intent(in) :: vec 
    integer, intent(out) :: nstate
    integer, dimension(ndim),intent(out) :: nisri

    integer :: i
    integer, dimension(ndim) :: indx
    
    nstate=0
    
!!$    call dsortqx("D",ndim,vec(:),1,indx(:))
    
    call dsortindxa1("D",ndim,vec(:),indx(:))
    
    nisri(:)=-1
    
    do i=1,ndim
    if(vec(indx(i)) .ge. stsel) then
       nstate=nstate+1
       nisri(i)=indx(i)
       write(ilog,'(a)') "selecting states nr.",indx(i)
    else
       exit
    end if
    end do
    
  end subroutine select_fstate_ryd

!!$---------------------------------------------------
!!$---------------------------------------------------

  subroutine select_fstate(ndim,stsel,vec,nstate,nisri)

    integer, intent(in) :: ndim
    real(d), intent(in) :: stsel
    real(d), dimension(ndim),intent(in) :: vec 
    integer,intent(out) :: nstate
    integer, dimension(ndim),intent(out) :: nisri

    integer :: i
    integer, dimension(ndim) :: indx
    real(d), dimension(ndim) :: coeff
    
    nstate=0
    
    do i=1,ndim
       coeff(i)=abs(vec(i))**2 
    end do

!!$    call dsortqx("D",ndim,coeff(:),1,indx(:))
   
    call dsortindxa1("D",ndim,coeff(:),indx(:)) 
    
    nisri(:)=-1
    
    do i=1,ndim
    if(coeff(indx(i)) .ge. stsel) then
       nstate=nstate+1
       nisri(i)=indx(i)
    else
       exit
    end if
    end do
    
    
  end subroutine select_fstate
    
!!$-----------------------------------------------






  subroutine select_atom_is_ALL(kpq,Simmetry,dimensione)
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    integer, intent(in) :: Simmetry    
    integer, intent(inout) :: dimensione   

    integer :: i,ap,n2pne
    integer :: isym,a,cah
    integer, dimension(7) :: col
    
    
100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))
    
!    write(ilog,'(a)') 'Selecting INITIAL  1h1p subspace'
    
    dimensione = 0
    kpq(1,0)=0

    do ap=nOcc+1,nBas
       a=roccnum(ap)
       do i=1,hcentre(0)
           cah=hcentre(i)
          isym=MT(orbSym(cah),orbSym(a))
          if(isym .eq. Simmetry) then
             kpq(1,0)=kpq(1,0)+1
             call fill_indices(col(:),1,1,a,-1,cah,-1,0)
             kpq(:,kpq(1,0))=col(:)
          end if
       end do
    end do
       
!    write(ilog,100) "Number of 1h-1p INITIAL configs in the IS ADC",kpq(1,0) 
!    write(ilog,103)
!    write(ilog,*) "INITIAL Singly excited configurations allowed in the Fin. St. Manif."
!    write(ilog,101) "CNF","SPN","HL1","PT1"
!    do i=1,kpq(1,0)
!       write(ilog,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
!    enddo
       
       
   dimensione = kpq(1,0)
 
  end subroutine select_atom_is_ALL



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    subroutine select_atom_d_ALL(kpq,flag,Simmetry,dimensione)

!!$Includes 2p2pnxnx configurations in the both, initial and final states.
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    integer, intent(in) :: flag
    integer, intent(in) :: Simmetry    
    integer, intent(inout) :: dimensione   

    integer :: i,j,a,b,k
    integer :: ih,jh,ap,bp
    integer :: cnti,cntf
    integer :: isym1, isym2
    integer, dimension(7) :: col
    real(d) :: einit,ei,ej

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I5,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-")) 

    einit=abs(e(hinit))

! Filling the kpq arrays in the order: a=b,i=j;a|=b,i=j;a=b,i|=j;a|=b,i|=j(I,II).

    if (flag .eq. -1) then
       
!       write(ilog,'(a)') 'Selecting INITIAL 2h2p subspace'
       
       kpq(2:5,0)=0
       cnti=kpq(1,0)
    
!!$a=b,i=j i=hcentre(ih)
       
       if(Simmetry .eq. 1) then
          do ih=1,hcentre(0)
             i=hcentre(ih)
             ei=abs(e(i))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                if(einit .le. 2._d*ei) then
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                else
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                end if
             end do
          end do
       end if

!        if(nirrep .eq. 1) then
!           do ih=1,hneighb(0)
!              i=hneighb(ih)
!              ei=abs(e(i))
!              do ap=nOcc+1,nBas
!                 a=roccnum(ap)
!                 if(einit .le. 2._d*ei) then
!                    cnti=cnti+1
!                    kpq(2,0)=kpq(2,0)+1
!                    call fill_indices(col(:),2,1,a,a,i,i,1)
!                    kpq(:,cnti)=col(:)
!                 else
!                    cnti=cnti+1
!                    kpq(2,0)=kpq(2,0)+1
!                    call fill_indices(col(:),2,1,a,a,i,i,1)
!                    kpq(:,cnti)=col(:)
!                 end if
!              end do
!           end do
!        end if

!!$a|=b,i=j
  
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do bp=ap+1,nBas
             b=roccnum(bp)
             isym1=MT(orbSym(a),orbSym(b))
             if(isym1 .eq. Simmetry) then
                do ih=1,hcentre(0)
                   i=hcentre(ih)
                   ei=abs(e(i))
                   if(einit .le. 2._d*ei) then
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   else
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   end if
                end do
             end if
          end do
       end do


!        do ap=nOcc+1,nBas
!           a=roccnum(ap)
!           do bp=ap+1,nBas
!              b=roccnum(bp)
!              isym1=MT(orbSym(a),orbSym(b))
!              if(isym1 .eq. nirrep) then
!                 do ih=1,hneighb(0)
!                    i=hneighb(ih)
!                    ei=abs(e(i))
!                    if(einit .le. 2._d*ei) then
!                       cnti=cnti+1
!                       kpq(3,0)=kpq(3,0)+1
!                       call fill_indices(col(:),2,1,a,b,i,i,2)
!                       kpq(:,cnti)=col(:)
!                    else
!                       cnti=cnti+1
!                       kpq(3,0)=kpq(3,0)+1
!                       call fill_indices(col(:),2,1,a,b,i,i,2)
!                       kpq(:,cnti)=col(:)
!                    end if
!                 end do
!              end if
!           end do
!        end do
!a=b,i|=j
  
       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))
          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             if (isym1 .eq. Simmetry) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                      cnti=cnti+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cnti)=col(:)
                end do
             end if
          end do
       end do

!        do ih=1,hneighb(0)
!          i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=1,hcentre(0)
!              j=hcentre(jh)
!             ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              if (isym1 .eq. nirrep) then
!                 do ap=nOcc+1,nBas
!                    a=roccnum(ap)
!                       cnti=cnti+1
!                       kpq(4,0)=kpq(4,0)+1
!                       call fill_indices(col(:),2,1,a,a,i,j,3) 
!                       kpq(:,cnti)=col(:)
!                 end do
!              end if
!           end do
!        end do

!        do ih=1,hneighb(0)
!           i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=ih+1,hneighb(0)
!              j=hneighb(jh)
!              ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              if (isym1 .eq. nirrep) then
!                 do ap=nOcc+1,nBas
!                    a=roccnum(ap)
!                       cnti=cnti+1
!                       kpq(4,0)=kpq(4,0)+1
!                       call fill_indices(col(:),2,1,a,a,i,j,3) 
!                       kpq(:,cnti)=col(:)
!                 end do
!              end if
!           end do
!        end do
    
!a|=b,i|=j spin I

       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))
          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                do bp=ap+1,nBas
                   b=roccnum(bp)
                   isym2=MT(orbSym(a),orbSym(b))
                   if(MT(isym1,isym2) .eq. Simmetry) then 
                      if(einit .le. (ei+ej)) then
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      else
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      end if
                   end if
                end do
             end do
          end do
       end do

!        do ih=1,hneighb(0)
!           i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=1,hcentre(0)
!              j=hcentre(jh)
!              ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              do ap=nOcc+1,nBas
!                a=roccnum(ap)
!                 do bp=ap+1,nBas
!                    b=roccnum(bp)
!                    isym2=MT(orbSym(a),orbSym(b))
!                    if(MT(isym1,isym2) .eq. nirrep) then 
!                       if(einit .le. (ei+ej)) then
!                         cnti=cnti+1
!                          kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       else
!                          cnti=cnti+1
!                         kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       end if
!                    end if
!                 end do
!              end do
!           end do
!        end do

! aukommentier -> Ne sonst HF

!        do ih=1,hneighb(0)
!           i=hneighb(ih)
!           ei=abs(e(i))
!           do jh=ih+1,hneighb(0)
!              j=hneighb(jh)
!              ej=abs(e(j))
!              isym1=MT(orbSym(i),orbSym(j))
!              do ap=nOcc+1,nBas
!                 a=roccnum(ap)
!                 do bp=ap+1,nBas
!                    b=roccnum(bp)
!                    isym2=MT(orbSym(a),orbSym(b))
!                    if(MT(isym1,isym2) .eq. nirrep) then 
!                       if(einit .le. (ei+ej)) then
!                          cnti=cnti+1
!                          kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       else
!                          cnti=cnti+1
!                          kpq(5,0)=kpq(5,0)+1
!                          call fill_indices(col(:),2,11,a,b,i,j,4)
!                          kpq(:,cnti)=col(:)
!                       end if
!                    end if
!                 end do
!              end do
!           end do
!        end do
    
       kpq(:,cnti+1:cnti+kpq(5,0))=kpq(:,cnti+1-kpq(5,0):cnti)
       kpq(2,cnti+1:cnti+kpq(5,0))=12
       kpq(7,cnti+1:cnti+kpq(5,0))=5

       
!       write(ilog,100) "Number of 2h-2p |abij> INITIAL configs in the IS ADC", cnti+kpq(5,0)-kpq(1,0)
!       write(ilog,103)
!       write(ilog,*) "INITIAL  Doubly excited |abij> configs allowed in the Init. St. Manif."
!       write(ilog,101) "CNF","SPN","HL1","HL2","PT1","PT2"
!    
!       do k = kpq(1,0) + 1 , cnti + kpq(5,0) 
!          write(ilog,102) k,kpq(2,k),kpq(3,k),&
!               kpq(4,k),kpq(5,k),kpq(6,k)
!       end do

       dimensione = cnti + kpq(5,0)


!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
    elseif(flag .eq. 1) then
       
!       write(ilog,'(a)') 'Selecting final 2h2p subspace'
       kpq(2:5,0)=0
       cntf=kpq(1,0)

       kpq(:,cntf+1:cntf+kpq(5,0))=kpq(:,cntf+1-kpq(5,0):cntf)
       kpq(2,cntf+1:cntf+kpq(5,0))=12 
       kpq(7,cntf+1:cntf+kpq(5,0))=5
    
!       write(ilog,100) "Number of 2h-2p |abij> configs in the FS ADC", cntf+kpq(5,0)-kpq(1,0)
!       write(ilog,103)
!       write(ilog,*) " Doubly excited |abij> configs allowed in the Fin. St. Manif."
!       write(ilog,101) "CNF","SPN","HL1","HL2","PT1","PT2"
!       
!       do k = kpq(1,0) + 1 , cntf + kpq(5,0)
!          write(ilog,102) k,kpq(2,k),kpq(3,k),&
!               kpq(4,k),kpq(5,k),kpq(6,k)
!       end do


       dimensione = cnti + kpq(5,0)

       
    end if
    
  end subroutine select_atom_d_ALL
!!$--------------------------------------------



  
 subroutine BUILD_ALL_SYMMETRY_SPACES(NSYMA,DIMEN)

  INTEGER, INTENT(IN) :: NSYMA
  INTEGER, DIMENSION(NSYMA), INTENT(INOUT) :: DIMEN

    integer, dimension(7,0:nBas**2*nOcc**2) :: kpq

  INTEGER :: dimensione
  INTEGER :: Simmetry

  do Simmetry = 1 , NSYMA

  dimensione = 0

  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )

  DIMEN(Simmetry) = dimensione


  end do


 end subroutine BUILD_ALL_SYMMETRY_SPACES




 subroutine BUILD_ALL_SYMMETRY_SPACES_ADC1(NSYMA,DIMEN)

  INTEGER, INTENT(IN) :: NSYMA
  INTEGER, DIMENSION(NSYMA), INTENT(INOUT) :: DIMEN

    integer, dimension(7,0:nBas**2*nOcc**2) :: kpq

  INTEGER :: dimensione
  INTEGER :: Simmetry

  do Simmetry = 1 , NSYMA

  dimensione = 0

  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )

  DIMEN(Simmetry) = dimensione


  end do


 end subroutine BUILD_ALL_SYMMETRY_SPACES_ADC1













end module select_fano

    
    
