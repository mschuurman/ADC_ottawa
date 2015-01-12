module photoionisation 
  
  use constants 
  use parameters
  use misc
  
  implicit none
  
contains

!!$------------------------------------------------------
  
  subroutine get_dominant(ndim,ene,fosc,indx)
    
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: ene, fosc 
    integer, dimension(norder), intent(out) :: indx
    
    integer :: i,j,numfosc,nbound,ndelta
    real(d) :: maxfosc,fsum
    
!!$ the subroutine subdivides the tda spectrum into norder regions; inside each region a state
!!$ with the largest fosc is found and its number is saved. Inefficient procedure.
    

    ndelta=ndim/norder
    nbound=1
    
    do i= 1,norder-1
       fsum=0._d
       numfosc=nbound
       maxfosc=fosc(nbound)
       do j= nbound,nbound+ndelta
          fsum=fsum+fosc(j)
          if (fosc(j) .gt. maxfosc) then
             maxfosc=fosc(j)
             numfosc=j
          end if
       end do
       write(6,*) 'The region ',i,' lies between ',ene(nbound),' and ', ene(nbound+ndelta)
       write(6,*) 'The dominant state in region', i,' is ',numfosc,' and carries',&
            maxfosc/fsum,' osc.strength','(',maxfosc,')'   
       indx(i)=numfosc
       nbound=nbound+ndelta+1
    end do
    
    fsum=0._d
    numfosc=nbound
    maxfosc=fosc(nbound)
    do j= nbound,ndim
       fsum=fsum+fosc(j)
       if (fosc(j) .gt. maxfosc) then
          maxfosc=fosc(j)
          numfosc=j
       end if
    end do
    write(6,*) 'The region ',i,' lies between ',ene(nbound),' and ', ene(ndim)
    write(6,*) 'The dominant state in region', i,' is ',numfosc,' and carries',&
         maxfosc/fsum,' osc.strength','(',maxfosc,')' 
    indx(i)=numfosc

  end subroutine get_dominant
    
!!$--------------------------
!!$--------------------------

  subroutine get_lanc_conf(ndim,arr)
    
    integer, intent(in) :: ndim
    real(d), dimension(ndim,norder), intent(in) :: arr
    
    integer :: i,j
    integer, dimension(ndim) :: indx
    real(d), dimension(ndim) :: coeff
    
102 FORMAT(5(F8.6,"(",I4,")",1x))
    
    do i= 1,norder
       coeff(:)=arr(:,i)**2
       call dsortindxa1("D",ndim,coeff(:),indx(:))
       write(6,*) 'Region ',i,' Conf.: ',(indx(j),j=1,5)
       write(6,102) (coeff(indx(j)),indx(j),j=1,5)
       
    end do
    
  end subroutine get_lanc_conf

end module photoionisation
       
    
    
