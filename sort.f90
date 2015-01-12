program main

  use constants
  implicit none
  
    integer :: ndim,i
    real(d), dimension(:), allocatable:: arrin
    integer, dimension(:), allocatable :: indx 
    
    ndim=8
    
    allocate(arrin(8),indx(8))
    arrin(1:8)=(/ 4._d,7._d,1._d,5._d,12._d,105._d,22._d,0.5_d /)
    call indexx('D',ndim,arrin,indx)
    
    write(6,*) arrin(:)
    write(6,*) (arrin(indx(i)), i=1,ndim)

  end program main
  
  subroutine indexx(order,ndim,arrin,indx)

    use constants
    implicit none

    character(1), intent(in) :: order
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: arrin
    integer, dimension(ndim), intent(inout) :: indx
    
    integer :: i,l,ir,indxt,j
    real(d) :: q
!!$ The subroutine is taken from the NR p233, employs heapsort.

    do i= 1,ndim
       indx(i)=i
    end do
    
    l=ndim/2+1
    ir=ndim
    
    if(order .eq. 'D') then
       
10     continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
20     if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .gt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .gt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 20
       end if
       indx(i)=indxt
       go to 10
       
    elseif(order .eq. 'A') then
       
100    continue
       if(l .gt. 1) then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
       else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir .eq. 1) then
             indx(1)=indxt
             return
          end if
       end if
       
       i=l
       j=l+l
       
200    if(j .le. ir) then
          if(j .lt. ir) then
             if(arrin(indx(j)) .lt. arrin(indx(j+1))) j=j+1 !
          end if
          if(q .lt. arrin(indx(j))) then !
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          go to 200
       end if
       indx(i)=indxt
       go to 100
       
    end if
       
  end subroutine indexx
  
