module sym_allowed
  
  
  use constants
  use parameters
  implicit none

contains

  subroutine get_symallowed_adc(kpq)
    
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(inout) :: kpq
    integer :: a,b,i,j,isym,isym1,bookmark

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

          write(6,*) 'MT 1 50', MT(orbsym(1),orbsym(50))    
          write(6,*) 'MT 1 51', MT(orbsym(1),orbsym(51))
          write(6,*) 'MT 1 64', MT(orbsym(1),orbsym(64)) 



!!$  Allowed singles

    do i=1,nOcc
       do a=nOcc+1,nBas
!          write(6,*) 'MT 1 50', MT(orbsym(1),orbsym(50))  
!          write(6,*) 'MT 1 51', MT(orbsym(1),orbsym(51))
!          write(6,*) 'MT 1 64', MT(orbsym(1),orbsym(64))    
          isym=MT(orbsym(roccnum(a)),orbsym(roccnum(i)))
          if (isym .eq. nirrep) then
             kpq(1,0)=kpq(1,0)+1
             bookmark=bookmark+1
             kpq(1,bookmark)=1
             kpq(2,bookmark)=1
             kpq(3,bookmark)=roccnum(i)
             kpq(5,bookmark)=roccnum(a)
          end if
       end do
    end do
    write(6,100) "Number of 1h-1p configs in the ADC matrix", kpq(1,0)
    write(6,103)
    write(6,*) "Singly excited configurations"
    write(6,101) "Nr","SPN","HL1","PT1"       
    do i=1,kpq(1,0)
       write(6,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
    end do

!!$Class 2, doubles i=j,a=b

    if (nirrep .eq. 1) then

       do i=1,nOcc
          do a=nOcc+1,nBas

             kpq(2,0)=kpq(2,0)+1
             bookmark=bookmark+1
             kpq(1,bookmark)=2
             kpq(2,bookmark)=1
             kpq(3,bookmark)=roccnum(i)
             kpq(4,bookmark)=roccnum(i)
             kpq(5,bookmark)=roccnum(a)
             kpq(6,bookmark)=roccnum(a)
             
          end do
       end do
       
       write(6,200) "Number of 2h-2p configs i=j,a=b ", kpq(2,0)
       write(6,203)
       write(6,*) "Doubly excited configs i=j,a=b"
       write(6,201) "Nr","SPN","HL1","HL2","PT1","PT2"
       do i=kpq(1,0)+1,bookmark
          write(6,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
       enddo
    end if
    
    write(6,200) "Number of 2h-2p configs i=j,a=b", kpq(2,0)

!!$Class 3, doubles i=j,a|=b 
    
    do i=1, nOcc
       do a=nOcc+1,nBas
          do b=a+1,nBas
            
             isym=MT(orbsym(roccnum(a)),orbsym(roccnum(b)))
             
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
          end do
       end do
    end do
    
    write(6,200) "Number of 2h-2p configs i=j,a|=b ", kpq(3,0)
    write(6,203)
    write(6,*) "Doubly excited configs i=j,a|=b "
    write(6,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpq(1,0)+kpq(2,0)+1,bookmark
       write(6,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
    enddo

!!$Class 4, doubles i|=j,a=b
    
    do i=1,nOcc
       do j=i+1,nOcc
          do a=nOcc+1,nBas
             
             isym=MT(orbsym(roccnum(i)),orbsym(roccnum(j)))
             
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
          end do
       end do
    end do
    
    write(6,200) "Number of 2h-2p configs i|=j,a=b ", kpq(4,0)
    write(6,203)
    write(6,*) "Doubly excited configs i|=j,a=b "
    write(6,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpq(1,0)+kpq(2,0)+kpq(3,0)+1,bookmark
       write(6,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
    enddo

    
!!$Class 5a, doubles i|=j,a|=b spin case I
    
    do i=1,nOcc
       do j=i+1,nOcc
          do a=nOcc+1,nBas
             do b=a+1,nBas
                
                isym=MT(orbsym(roccnum(i)),orbsym(roccnum(j)))
                isym1=MT(orbsym(roccnum(a)),orbsym(roccnum(b)))
                
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
             end do
          end do
       end do
    end do

    kpq(1,bookmark+1:bookmark+kpq(5,0))=2
    kpq(2,bookmark+1:bookmark+kpq(5,0))=12
    kpq(3,bookmark+1:bookmark+kpq(5,0))=kpq(3,bookmark-kpq(5,0)+1:bookmark)
    kpq(4,bookmark+1:bookmark+kpq(5,0))=kpq(4,bookmark-kpq(5,0)+1:bookmark)
    kpq(5,bookmark+1:bookmark+kpq(5,0))=kpq(5,bookmark-kpq(5,0)+1:bookmark)
    kpq(6,bookmark+1:bookmark+kpq(5,0))=kpq(6,bookmark-kpq(5,0)+1:bookmark)
    
    write(6,200) "Number of 2h-2p configs i|=j,a|=b. Spin1 ", kpq(5,0)
    write(6,203)
    write(6,*) "Doubly excited configs i|=j,a|=b "
    write(6,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+1,bookmark
       write(6,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
    enddo

    write(6,200) "Number of 2h-2p configs i|=j,a|=b. Spin2 ", kpq(5,0)
    write(6,203)
    write(6,*) "Doubly excited configs i|=j,a|=b "
    write(6,201) "Nr","SPN","HL1","HL2","PT1","PT2"
    do i=bookmark+1,bookmark+kpq(5,0)
       write(6,202) i,kpq(2,i),kpq(3,i),kpq(4,i),kpq(5,i),kpq(6,i)
    enddo
    write(6,203)
    write(6,203)
    
    write(6,200) "Number of 1h-1p configs", kpq(1,0)
    write(6,200) "Number of 2h-2p configs", kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)
    write(6,203)


  end subroutine get_symallowed_adc

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

  subroutine get_symallowed_tda(kpq)
    
    integer, dimension(7,0:nBas**2*4*nOcc**2), intent(inout) :: kpq

    integer :: a,i,isym,bookmark

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I4,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))

200 FORMAT(/,3("*"),A50,3x,I4)

    kpq(:,0)=0
    
    do i=1,nOcc
       do a=nOcc+1,nBas
          
          isym=MT(orbsym(roccnum(a)),orbsym(roccnum(i)))
          if (isym .eq. nirrep) then
             kpq(1,0)=kpq(1,0)+1
             bookmark=bookmark+1
             kpq(1,bookmark)=1
             kpq(2,bookmark)=1
             kpq(3,bookmark)=roccnum(i)
             kpq(5,bookmark)=roccnum(a)
          end if
       end do
    end do
    write(6,100) "Number of 1h-1p configs in the ADC matrix", kpq(1,0)
    write(6,103)
    write(6,*) "Singly excited configurations"
    write(6,101) "Nr","SPN","HL1","PT1"       
    do i=1,kpq(1,0)
       write(6,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
    end do

    write(6,103)
    write(6,103)
    
    write(6,200) "Number of 1h-1p configs", kpq(1,0)
    write(6,103)  
    
  end subroutine get_symallowed_tda
  

end module sym_allowed
