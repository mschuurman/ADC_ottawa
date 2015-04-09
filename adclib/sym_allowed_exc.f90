module sym_allowed_exc
  
  
  use constants
  use parameters
  implicit none

contains

  subroutine get_symallowed_adc(kpq,kpqf,kpqd)
    
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(inout) :: kpq
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(inout) :: kpqf
    integer, dimension(7,0:2*nBas**2*4*nOcc**2), intent(inout) :: kpqd
    integer :: a,b,i,j,isym,isym1,bookmark,bookmarkf,bookmarkd
    integer :: a1,b1,isym2

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A5,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I6,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))

200 FORMAT(/,3("*"),A50,3x,I7)
201 FORMAT(/,("*"),3x,A5,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
202 FORMAT(("*"),3x,I7,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
203 FORMAT(/,60("-"))

    kpq(:,0)=0
    bookmark=0
    kpqf(:,0)=0
    bookmarkf=0
    kpqd(:,0)=0
    bookmarkd=0



!  if(statenumber.eq.0) then



!!!!!!!!!!!!!!!! from ground state !!!!!!!!!!!!!!!!!!!!!!!!!!



!!$  Allowed singles

    do i=1,nOcc
       do a=nOcc+1,nBas
          
          isym=MT(orbsym(roccnum(a)),orbsym(roccnum(i)))

         if ((isym .eq. nirrep).or.(isym .eq.1)) then
             kpqd(1,0)=kpqd(1,0)+1
             bookmark=bookmarkd+1
             kpqd(1,bookmarkd)=1
             kpqd(2,bookmarkd)=1
             kpqd(3,bookmarkd)=roccnum(i)
             kpqd(5,bookmarkd)=roccnum(a)
          

          if (isym .eq. nirrep) then
             kpq(1,0)=kpq(1,0)+1
             bookmark=bookmark+1
             kpq(1,bookmark)=1
             kpq(2,bookmark)=1
             kpq(3,bookmark)=roccnum(i)
             kpq(5,bookmark)=roccnum(a)
          end if
         end if

       end do
    end do
    
    write(ilog,100) "Number of 1h-1p TOTAL configs in the ADC D  matrix", kpqd(1,0)
    write(ilog,103)
    write(ilog,*) "Singly excited TOTAL D  configurations"
    write(ilog,101) "Nr","SPN","HL1","PT1"       
    do i=1,kpqd(1,0)
       write(ilog,102) i,kpqd(2,i),kpqd(3,i),kpqd(5,i)
    end do


    write(ilog,100) "Number of 1h-1p INITIAL configs in the ADC matrix", kpq(1,0)
    write(ilog,103)
    write(ilog,*) "Singly excited INITIAL configurations"
    write(ilog,101) "Nr","SPN","HL1","PT1"       
    do i=1,kpq(1,0)
       write(ilog,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
    end do

!!$Class 2, doubles i=j,a=b


!    if (nirrep .eq. 1) then

       do i=1,nOcc
          do a=nOcc+1,nBas

!     if (nirrep .eq. 1) then

!             kpq(2,0)=kpq(2,0)+1
!             bookmark=bookmark+1
!             kpq(1,bookmark)=2
!             kpq(2,bookmark)=1
!             kpq(3,bookmark)=roccnum(i)
!             kpq(4,bookmark)=roccnum(i)
!             kpq(5,bookmark)=roccnum(a)
!             kpq(6,bookmark)=roccnum(a)
     
!      end if           

             kpqd(2,0)=kpqd(2,0)+1
             kpqd(1,bookmarkd)=2
             kpqd(2,bookmarkd)=1
             kpqd(3,bookmarkd)=roccnum(i)
             kpqd(4,bookmarkd)=roccnum(i)
             kpqd(5,bookmarkd)=roccnum(a)
             kpqd(6,bookmarkd)=roccnum(a)
             
          end do
       end do
       

       write(ilog,200) "Number of 2h-2p TOTAL  configs i=j,a=b ", kpqd(2,0)
       write(ilog,203)
       write(ilog,*) "Doubly excited TOTAL  configs i=j,a=b"
       write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
       do i=kpqd(1,0)+1,bookmarkd
          write(ilog,202) i,kpqd(2,i),kpqd(3,i),kpqd(4,i),kpqd(5,i),kpqd(6,i)
       enddo
!    end if
    
    write(ilog,200) "Number of 2h-2p TOTAL configs i=j,a=b", kpqd(2,0)



       write(ilog,200) "Number of 2h-2p INITIAL configs i=j,a=b ", kpq(2,0)
       write(ilog,203)
       write(ilog,*) "Doubly excited INITIAL configs i=j,a=b"
       write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
       do i=kpq(1,0)+1,bookmark
          write(ilog,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
       enddo
!    end if
    
    write(ilog,200) "Number of 2h-2p INITIAL configs i=j,a=b", kpq(2,0)

!!$Class 3, doubles i=j,a|=b 
    
    do i=1, nOcc
       do a=nOcc+1,nBas
          do b=a+1,nBas
            
             isym=MT(orbsym(roccnum(a)),orbsym(roccnum(b)))
             
            if ((isym .eq. nirrep).or.(isym .eq.1)) then

                kpqd(3,0)=kpqd(3,0)+1
                bookmarkd=bookmarkd+1
                kpqd(1,bookmarkd)=2
                kpqd(2,bookmarkd)=1
                kpqd(3,bookmarkd)=roccnum(i)
                kpqd(4,bookmarkd)=roccnum(i)
                kpqd(5,bookmarkd)=roccnum(a)
                kpqd(6,bookmarkd)=roccnum(b)
 

             if(isym .eq. nirrep) then
                kpq(3,0)=kpq(3,0)+1
                bookmark=bookmark+1
                kpq(1,bookmark)=2
                kpq(2,bookmark)=1
                kpq(3,bookmark)=roccnum(i)
                kpq(4,bookmark)=roccnum(i)
                kpq(5,bookmark)=roccnum(a)
                kpq(6,bookmark)=roccnum(b)
             end if
            end if
          end do
       end do
    end do

    write(ilog,200) "Number of 2h-2p TOTAL configs i=j,a|=b ", kpqd(3,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited TOTAL configs i=j,a|=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpqd(1,0)+kpqd(2,0)+1,bookmarkd
       write(ilog,202) i,kpqd(2,i),kpqd(3,i),kpqd(4,i),kpqd(5,i),kpqd(6,i)
    enddo



    
    write(ilog,200) "Number of 2h-2p INITIAL configs i=j,a|=b ", kpq(3,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited INITIAL configs i=j,a|=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpq(1,0)+kpq(2,0)+1,bookmark
       write(ilog,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
    enddo

!!$Class 4, doubles i|=j,a=b
    
    do i=1,nOcc
       do j=i+1,nOcc
          do a=nOcc+1,nBas
             
             isym=MT(orbsym(roccnum(i)),orbsym(roccnum(j)))
 

            if ((isym .eq. nirrep).or.(isym .eq.1)) then

                kpqd(4,0)=kpqd(4,0)+1
                bookmarkd=bookmarkd+1
                kpqd(1,bookmarkd)=2
                kpqd(2,bookmarkd)=1
                kpqd(3,bookmarkd)=roccnum(i)
                kpqd(4,bookmarkd)=roccnum(j)
                kpqd(5,bookmarkd)=roccnum(a)
                kpqd(6,bookmarkd)=roccnum(a)

             if(isym .eq. nirrep) then
                kpq(4,0)=kpq(4,0)+1
                bookmark=bookmark+1
                kpq(1,bookmark)=2
                kpq(2,bookmark)=1
                kpq(3,bookmark)=roccnum(i)
                kpq(4,bookmark)=roccnum(j)
                kpq(5,bookmark)=roccnum(a)
                kpq(6,bookmark)=roccnum(a)
             end if
            end if
          end do
       end do
    end do


    write(ilog,200) "Number of 2h-2p TOTAL configs i|=j,a=b ", kpqd(4,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited TOTAL configs i|=j,a=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpqd(1,0)+kpqd(2,0)+kpqd(3,0)+1,bookmarkd
       write(ilog,202) i,kpqd(2,i),kpqd(3,i),kpqd(4,i),kpqd(5,i),kpqd(6,i)
    enddo


    
    write(ilog,200) "Number of 2h-2p INITIAL configs i|=j,a=b ", kpq(4,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited INITIAL configs i|=j,a=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1,bookmark
       write(ilog,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
    enddo

    
!!$Class 5a, doubles i|=j,a|=b spin case I
    
    do i=1,nOcc
       do j=i+1,nOcc
          do a=nOcc+1,nBas
             do b=a+1,nBas
                
                isym=MT(orbsym(roccnum(i)),orbsym(roccnum(j)))
                isym1=MT(orbsym(roccnum(a)),orbsym(roccnum(b)))
 
               if ((MT(isym,isym1) .eq. nirrep).or.(MT(isym,isym1) .eq.1)) then

                   kpqd(5,0)=kpqd(5,0)+1
                   bookmarkd=bookmarkd+1
                   kpqd(1,bookmarkd)=2
                   kpqd(2,bookmarkd)=11
                   kpqd(3,bookmarkd)=roccnum(i)
                   kpqd(4,bookmarkd)=roccnum(j)
                   kpqd(5,bookmarkd)=roccnum(a)
                   kpqd(6,bookmarkd)=roccnum(b)

               
                if(nirrep .eq. MT(isym,isym1)) then
                   kpq(5,0)=kpq(5,0)+1
                   bookmark=bookmark+1
                   kpq(1,bookmark)=2
                   kpq(2,bookmark)=11
                   kpq(3,bookmark)=roccnum(i)
                   kpq(4,bookmark)=roccnum(j)
                   kpq(5,bookmark)=roccnum(a)
                   kpq(6,bookmark)=roccnum(b)
                end if
               end if
             end do
          end do
       end do
    end do

!!!!!!!!!!!!!! TOTAL CONFIG D MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    kpqd(1,bookmarkd+1:bookmarkd+kpqd(5,0))=2
    kpqd(2,bookmarkd+1:bookmarkd+kpqd(5,0))=12
    kpqd(3,bookmarkd+1:bookmarkd+kpqd(5,0))=kpqd(3,bookmarkd-kpqd(5,0)+1:bookmarkd)
    kpqd(4,bookmarkd+1:bookmarkd+kpqd(5,0))=kpqd(4,bookmarkd-kpqd(5,0)+1:bookmarkd)
    kpqd(5,bookmarkd+1:bookmarkd+kpqd(5,0))=kpqd(5,bookmarkd-kpqd(5,0)+1:bookmarkd)
    kpqd(6,bookmarkd+1:bookmarkd+kpqd(5,0))=kpqd(6,bookmarkd-kpqd(5,0)+1:bookmarkd)
    
    write(ilog,200) "Number of 2h-2p INITIAL configs i|=j,a|=b. Spin1 ", kpqd(5,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited INITIAL configs i|=j,a|=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpqd(1,0)+kpqd(2,0)+kpqd(3,0)+kpqd(4,0)+1,bookmarkd
       write(ilog,202) i,kpqd(2,i),kpqd(3,i),kpqd(4,i),kpqd(5,i),kpqd(6,i)
    enddo

    write(ilog,200) "Number of 2h-2p INITIAL configs i|=j,a|=b. Spin2 ", kpqd(5,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited INITIAL configs i|=j,a|=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=bookmarkd+1,bookmarkd+kpqd(5,0)
       write(ilog,202) i,kpqd(2,i),kpqd(3,i),kpqd(4,i),kpqd(5,i),kpqd(6,i)
    enddo
    write(ilog,203)
    write(ilog,203)
    
    write(ilog,200) "Number of 1h-1p INITIAL configs", kpqd(1,0)
    write(ilog,200) "Number of 2h-2p INITIAL configs", kpqd(2,0)+kpqd(3,0)+kpqd(4,0)+2*kpqd(5,0)
    write(ilog,203)


!!!!!!!!!!!!!! INITIAL CONFIG H MATRIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    kpq(1,bookmark+1:bookmark+kpq(5,0))=2
    kpq(2,bookmark+1:bookmark+kpq(5,0))=12
    kpq(3,bookmark+1:bookmark+kpq(5,0))=kpq(3,bookmark-kpq(5,0)+1:bookmark)
    kpq(4,bookmark+1:bookmark+kpq(5,0))=kpq(4,bookmark-kpq(5,0)+1:bookmark)
    kpq(5,bookmark+1:bookmark+kpq(5,0))=kpq(5,bookmark-kpq(5,0)+1:bookmark)
    kpq(6,bookmark+1:bookmark+kpq(5,0))=kpq(6,bookmark-kpq(5,0)+1:bookmark)
    
    write(ilog,200) "Number of 2h-2p INITIAL configs i|=j,a|=b. Spin1 ", kpq(5,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited INITIAL configs i|=j,a|=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1,bookmark
       write(ilog,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
    enddo

    write(ilog,200) "Number of 2h-2p INITIAL configs i|=j,a|=b. Spin2 ", kpq(5,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited INITIAL configs i|=j,a|=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=bookmark+1,bookmark+kpq(5,0)
       write(ilog,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
    enddo
    write(ilog,203)
    write(ilog,203)
    
    write(ilog,200) "Number of 1h-1p INITIAL configs", kpq(1,0)
    write(ilog,200) "Number of 2h-2p INITIAL configs", kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)
    write(ilog,203)


!!!!!!!!!!!!!!!!!!!!!!!!! from excited state !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (statenumber.gt.0) then




    do i=1,nOcc
       do a=nOcc+1,nBas
          
          isym=MT(orbsym(roccnum(a)),orbsym(roccnum(i)))
          if (isym .eq. 1) then
             kpqf(1,0)=kpqf(1,0)+1
             bookmarkf=bookmarkf+1
             kpqf(1,bookmarkf)=1
             kpqf(2,bookmarkf)=1
             kpqf(3,bookmarkf)=roccnum(i)
             kpqf(5,bookmarkf)=roccnum(a)
          end if
       end do
    end do
    write(ilog,100) "Number of 1h-1p configs in the ADC matrix for final states", kpqf(1,0)
    write(ilog,103)
    write(ilog,*) "Singly excited  FINAL  configurations"
    write(ilog,101) "Nr","SPN","HL1","PT1"       
    do i=1,kpqf(1,0)
       write(ilog,102) i,kpqf(2,i),kpqf(3,i),kpqf(5,i)
    end do

!!$Class 2, doubles i=j,a=b

!    if (nirrep .eq. 1) then
!!!!!!!!!!!!!!!!!!  IT IS ALWAYS GOOD (POSSIBLE FINAL STATE BY SYMMETRY) !!!!!!!!!!!!!!!!!


       do i=1,nOcc
          do a=nOcc+1,nBas

             kpqf(2,0)=kpqf(2,0)+1
             bookmarkf=bookmarkf+1
             kpqf(1,bookmarkf)=2
             kpqf(2,bookmarkf)=1
             kpqf(3,bookmarkf)=roccnum(i)
             kpqf(4,bookmarkf)=roccnum(i)
             kpqf(5,bookmarkf)=roccnum(a)
             kpqf(6,bookmarkf)=roccnum(a)
             
          end do
       end do
       
       write(ilog,200) "Number of 2h-2p FINAL configs i=j,a=b ", kpqf(2,0)
       write(ilog,203)
       write(ilog,*) "Doubly excited FINAL configs i=j,a=b"
       write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
       do i=kpqf(1,0)+1,bookmarkf
          write(ilog,202) i,kpqf(2,i),kpqf(3,i),kpqf(4,i),kpqf(5,i),kpqf(6,i)
       enddo
!    end if
    
    write(ilog,200) "Number of 2h-2p FINAL configs i=j,a=b", kpqf(2,0)

!!$Class 3, doubles i=j,a|=b 
    
    do i=1, nOcc
       do a=nOcc+1,nBas
          do b=a+1,nBas
            
             isym=MT(orbsym(roccnum(a)),orbsym(roccnum(b)))
             
             if(isym .eq. 1) then
                kpqf(3,0)=kpqf(3,0)+1
                bookmarkf=bookmarkf+1
                kpqf(1,bookmarkf)=2
                kpqf(2,bookmarkf)=1
                kpqf(3,bookmarkf)=roccnum(i)
                kpqf(4,bookmarkf)=roccnum(i)
                kpqf(5,bookmarkf)=roccnum(a)
                kpqf(6,bookmarkf)=roccnum(b)
             end if
          end do
       end do
    end do
    
    write(ilog,200) "Number of 2h-2p FINAL configs i=j,a|=b ", kpqf(3,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited FINAL configs i=j,a|=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpqf(1,0)+kpqf(2,0)+1,bookmarkf
       write(ilog,202) i,kpqf(2,i),kpqf(3,i),kpqf(4,i),kpqf(5,i),kpqf(6,i)
    enddo

!!$Class 4, doubles i|=j,a=b
    
    do i=1,nOcc
       do j=i+1,nOcc
          do a=nOcc+1,nBas
             
             isym=MT(orbsym(roccnum(i)),orbsym(roccnum(j)))
             
             if(isym .eq. 1) then
                kpqf(4,0)=kpqf(4,0)+1
                bookmarkf=bookmarkf+1
                kpqf(1,bookmarkf)=2
                kpqf(2,bookmarkf)=1
                kpqf(3,bookmarkf)=roccnum(i)
                kpqf(4,bookmarkf)=roccnum(j)
                kpqf(5,bookmarkf)=roccnum(a)
                kpqf(6,bookmarkf)=roccnum(a)
             end if
          end do
       end do
    end do
    
    write(ilog,200) "Number of 2h-2p FINAL configs i|=j,a=b ", kpqf(4,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited FINAL configs i|=j,a=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+1,bookmarkf
       write(ilog,202) i,kpqf(2,i),kpqf(3,i),kpqf(4,i),kpqf(5,i),kpqf(6,i)
    enddo

    
!!$Class 5a, doubles i|=j,a|=b spin case I
    
    do i=1,nOcc
       do j=i+1,nOcc
          do a=nOcc+1,nBas
             do b=a+1,nBas
                
                isym=MT(orbsym(roccnum(i)),orbsym(roccnum(j)))
                isym1=MT(orbsym(roccnum(a)),orbsym(roccnum(b)))
                
                if(1 .eq. MT(isym,isym1)) then
                   kpqf(5,0)=kpqf(5,0)+1
                   bookmarkf=bookmarkf+1
                   kpqf(1,bookmarkf)=2
                   kpqf(2,bookmarkf)=11
                   kpqf(3,bookmarkf)=roccnum(i)
                   kpqf(4,bookmarkf)=roccnum(j)
                   kpqf(5,bookmarkf)=roccnum(a)
                   kpqf(6,bookmarkf)=roccnum(b)
                end if
             end do
          end do
       end do
    end do

    kpqf(1,bookmarkf+1:bookmarkf+kpqf(5,0))=2
    kpqf(2,bookmarkf+1:bookmarkf+kpqf(5,0))=12
    kpqf(3,bookmarkf+1:bookmarkf+kpqf(5,0))=kpqf(3,bookmarkf-kpqf(5,0)+1:bookmarkf)
    kpqf(4,bookmarkf+1:bookmarkf+kpqf(5,0))=kpqf(4,bookmarkf-kpqf(5,0)+1:bookmarkf)
    kpqf(5,bookmarkf+1:bookmarkf+kpqf(5,0))=kpqf(5,bookmarkf-kpqf(5,0)+1:bookmarkf)
    kpqf(6,bookmarkf+1:bookmarkf+kpqf(5,0))=kpqf(6,bookmarkf-kpqf(5,0)+1:bookmarkf)
    
    write(ilog,200) "Number of 2h-2p FINAL configs i|=j,a|=b. Spin1 ", kpqf(5,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited FINAL configs i|=j,a|=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+1,bookmarkf
       write(ilog,202) i,kpqf(2,i),kpqf(3,i),kpqf(4,i),kpqf(5,i),kpqf(6,i)
    enddo

    write(ilog,200) "Number of 2h-2p FINAL configs i|=j,a|=b. Spin2 ", kpqf(5,0)
    write(ilog,203)
    write(ilog,*) "Doubly excited FINAL configs i|=j,a|=b "
    write(ilog,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=bookmarkf+1,bookmarkf+kpqf(5,0)
       write(ilog,202) i,kpqf(2,i),kpqf(3,i),kpqf(4,i),kpqf(5,i),kpqf(6,i)
    enddo
    write(ilog,203)
    write(ilog,203)
    
    write(ilog,200) "Number of 1h-1p FINAL configs", kpqf(1,0)
    write(ilog,200) "Number of 2h-2p FINAL configs", kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+2*kpqf(5,0)
    write(ilog,203)


  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  END  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  end subroutine get_symallowed_adc

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

  subroutine get_symallowed_tda(kpq,kpqf,kpqd)
    
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(inout) :: kpq
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(inout) :: kpqf
    integer, dimension(7,0:2*nBas**2*4*nOcc**2), intent(inout) :: kpqd
 

    integer :: a,i,isym,bookmark,bookmarkf,bookmarkd

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I4,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))

200 FORMAT(/,3("*"),A50,3x,I4)

    kpq(:,0)=0
    kpqf(:,0)=0
    kpqd(:,0)=0
    
    do i=1,nOcc
       do a=nOcc+1,nBas
          
          isym=MT(orbsym(roccnum(a)),orbsym(roccnum(i)))

        if ((isym .eq. nirrep).or.(isym .eq. 1)) then

             kpqd(1,0)=kpqd(1,0)+1
             bookmarkd=bookmarkd+1
             kpqd(1,bookmarkd)=1
             kpqd(2,bookmarkd)=1
             kpqd(3,bookmarkd)=roccnum(i)
             kpqd(5,bookmarkd)=roccnum(a)

          if (isym .eq. nirrep) then
             kpq(1,0)=kpq(1,0)+1
             bookmark=bookmark+1
             kpq(1,bookmark)=1
             kpq(2,bookmark)=1
             kpq(3,bookmark)=roccnum(i)
             kpq(5,bookmark)=roccnum(a)
          end if
         end if
       end do
    end do


    write(ilog,100) "Number of 1h-1p TOTAL configs in the ADC matrix", kpqd(1,0)
    write(ilog,103)
    write(ilog,*) "Singly excited TOTAL configurations"
    write(ilog,101) "Nr","SPN","HL1","PT1"       
    do i=1,kpqd(1,0)
       write(ilog,102) i,kpqd(2,i),kpqd(3,i),kpqd(5,i)
    end do

    write(ilog,103)
    write(ilog,103)
    
    write(ilog,200) "Number of 1h-1p TOTAL configs", kpqd(1,0)
    write(ilog,103)  






    write(ilog,100) "Number of 1h-1p INITIAL  configs in the ADC matrix", kpq(1,0)
    write(ilog,103)
    write(ilog,*) "Singly excited INITIAL configurations"
    write(ilog,101) "Nr","SPN","HL1","PT1"       
    do i=1,kpq(1,0)
       write(ilog,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
    end do

    write(ilog,103)
    write(ilog,103)
    
    write(ilog,200) "Number of 1h-1p INITIAL configs", kpq(1,0)
    write(ilog,103)  


   if (statenumber.gt.0) then

    do i=1,nOcc
       do a=nOcc+1,nBas
          
          isym=MT(orbsym(roccnum(a)),orbsym(roccnum(i)))
          if (isym .eq. 1) then
             kpqf(1,0)=kpqf(1,0)+1
             bookmarkf=bookmarkf+1
             kpqf(1,bookmarkf)=1
             kpqf(2,bookmarkf)=1
             kpqf(3,bookmarkf)=roccnum(i)
             kpqf(5,bookmarkf)=roccnum(a)
          end if
       end do
    end do
    write(ilog,100) "Number of 1h-1p FINAL configs in the ADC matrix", kpqf(1,0)
    write(ilog,103)
    write(ilog,*) "Singly excited FINAL configurations"
    write(ilog,101) "Nr","SPN","HL1","PT1"       
    do i=1,kpqf(1,0)
       write(ilog,102) i,kpqf(2,i),kpqf(3,i),kpqf(5,i)
    end do

    write(ilog,103)
    write(ilog,103)
    
    write(ilog,200) "Number of 1h-1p FINAL configs", kpqf(1,0)
    write(ilog,103)  
 
 
      
   end if


  end subroutine get_symallowed_tda
  

end module sym_allowed_exc
