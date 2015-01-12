module misc

!!$Contains miscellaneous ancillary functions and subroutines
  
  use constants
  use parameters  
  
  implicit none

contains
!!$------------------------------------  
  integer function kdelta(k,kpr)
    
    integer, intent(in) :: k,kpr
    
    kdelta=0

    if (k .eq. kpr) kdelta=1
    
  end function kdelta
!!$------------------------------------------------------------  

  logical function lkdelta(k,kpr)

    integer, intent(in) :: k,kpr

    lkdelta=.false.

    if (k .eq. kpr) lkdelta=.true.

  end function lkdelta

  subroutine diagonalise(nsout,ndim,arr)
    
    integer, intent(in) :: ndim,nsout
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: info,lwork,i,j
    real(d), dimension(ndim) :: w
    real(d), dimension(:), allocatable :: work
    
    integer, dimension(ndim) :: indx
    real(d), dimension(ndim) :: coeff
    
    real(d) :: dnrm2, nrm
    
    external dsyev 

100 FORMAT(60("-"),/)
101 FORMAT(3(A14,2x),/)
102 FORMAT(I10,2x,2(F10.6,2x),5(F8.6,"(",I4,")",1x))

    write(6,*) "Starting diagonalisation"
    write(6,100)
    lwork=3*ndim
    allocate(work(lwork))

    call dsyev("V","L", ndim, arr, ndim, w, work,lwork, info)
    
    write(6,101) "Eigenval","Ev","a.u."
    write(6,100)
    do i=1,nsout
       coeff(:)=arr(:,i)**2
!!$       call dsortqx("D",ndim,coeff(:),1,indx(:))
       call dsortindxa1("D",ndim,coeff(:),indx(:))
       
       write(6,102) i,w(i),w(i)*27.211396,(coeff(indx(j)),indx(j),j=1,5)
    end do
    
    deallocate(work)
    
  end subroutine diagonalise
!!$---------------------------------------------------
  subroutine vdiagonalise(ndim,arr,evector)

    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(inout) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: info,lwork, i,j
    integer*8 :: lworkl,ndiml
    
    real(d), dimension(ndim) :: w
    real(d), dimension(:), allocatable :: work

    integer, dimension(ndim) :: indx
    real(d), dimension(ndim) :: coeff

    external dsyev

    write(6,*) "Starting diagonalisation"
    lwork=3*ndim
    allocate(work(lwork))
    lworkl=int(lwork,lng)
    ndiml=int(ndim,lng)
    call dsyev("V","L",ndiml,arr(:,:),ndiml,w,work(:),lworkl,info)
    evector(:)=w(:)

    deallocate(work)

  end subroutine vdiagonalise

!!$----------------------------------------------------

  subroutine get_indices(col,a,b,j,k,spin)
    
    integer, dimension(7),intent(in) :: col
    integer, intent(out) :: a,b,j,k,spin 
    
    spin=col(2)
    j=col(3)
    k=col(4)
    a=col(5)
    b=col(6)
    
  end subroutine get_indices

!!$---------------------------------------------------------------------

  subroutine get_indices1(col,spin,a,b,j,k,type)
    
    integer, dimension(7),intent(in) :: col
    integer, intent(out) :: a,b,j,k,spin,type 
    
    spin=col(2)
    j=col(3)
    k=col(4)
    a=col(5)
    b=col(6)
    type=col(7)
    
  end subroutine get_indices1

!----------------------------------------------------------
!----------------------------------------------------------

  subroutine fill_indices(col,cnf,spin,a,b,j,k,type)
    
    integer, dimension(7),intent(out) :: col
    integer, intent(in) :: a,b,j,k,spin,cnf,type
    
    col(1)=cnf
    col(2)=spin
    col(3)=j
    col(4)=k
    col(5)=a
    col(6)=b
    col(7)=type
    
  end subroutine fill_indices

!-------------------------------------

  real(d) function dsp(ndim,vec1,vec2)
    
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: vec1,vec2
    
    integer :: i
    
    dsp=0._d
    
    do i=1, ndim
       dsp=dsp+vec1(i)*vec2(i)
    end do
    
  end function dsp
!!$---------------------------------------------
  subroutine dsortindxa1(order,ndim,arrin,indx)

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
       
  end subroutine dsortindxa1


!!$-------------------------------------------------

  subroutine dsortindxa(order,ndim,arr,indarr)
    
    character(1), intent(in) :: order
    integer, intent(in) :: ndim
    integer, dimension(ndim), intent(out) :: indarr
    real(d), dimension(ndim), intent(in) :: arr
    
    integer :: i, ifail
    integer, dimension(ndim) :: irank
    
!    external M01DAF,M01EBF
    
    ifail=0
!    call M01DAF (arr(:),1,ndim,order,irank(:),ifail)
    do i=1,ndim
       indarr(i)=i
    end do
!    call M01EBF (indarr(:),1,ndim,irank(:),ifail)
    
  end subroutine dsortindxa
!!$------------------------------------------------------
  
  subroutine table1(ndim1,ndim2,en,vspace)
    
    integer, intent(in) :: ndim1,ndim2
    real(d), dimension(ndim2), intent(in) :: en
    real(d), dimension(ndim1,ndim2), intent(in) :: vspace

    real(d), dimension(:), allocatable :: coeff
    integer, dimension(:), allocatable :: indx
    integer :: i,j

    allocate(coeff(ndim1),indx(ndim1))
    
100 FORMAT(60("-"),/)    
101 FORMAT(3(A10,2x),/)
102 FORMAT(I10,x,"|",2(F10.6,2x),"|",x,5(F8.6,"(",I4,")",1x))    
  
    write(6,101) "Eigenval","a.u.","Ev"
    write(6,100)
    do i=1,ndim2
       coeff(:)=vspace(:,i)**2
       call dsortindxa1("D",ndim1,coeff(:),indx(:))
       write(6,102) i,en(i),en(i)*27.211396,(coeff(indx(j)),indx(j),j=1,5)
    end do 

    deallocate(coeff,indx)

  end subroutine table1

!!$-----------------------------------------------------------

  subroutine table2(ndim1,ndim2,en,vspace,tmvec,osc_str)
    
    integer, intent(in) :: ndim1,ndim2
    real(d), dimension(ndim2), intent(in) :: en,tmvec,osc_str
    real(d), dimension(ndim1,ndim2), intent(in) :: vspace

    real(d), dimension(:), allocatable :: coeff
    integer, dimension(:), allocatable :: indx
    integer :: i,j,nlim

100 FORMAT(60("-"),/)    
101 FORMAT(4(A10,2x),/)
102 FORMAT(I10,x,"|",4(F10.5,2x),"|",x,5(F8.6,"(",I7,")",1x)) 

    allocate(coeff(ndim1),indx(ndim1))
    nlim=5
    if (nlim .gt. ndim1) nlim=ndim1
  
    write(6,101) "Eigenval","a.u.","Ev","Trans.Mom."
    write(6,100)
    do i=1,ndim2
       coeff(:)=vspace(:,i)**2
       call dsortindxa1("D",ndim1,coeff(:),indx(:))
       
       coeff(:)=vspace(:,i)
       write(6,102) i,en(i),en(i)*27.211396,tmvec(i),osc_str(i),(coeff(indx(j)),indx(j),j=1,nlim)

!!$       write(6,102) i,en(i),en(i)*27.211396,tmvec(i),(coeff(indx(j)),indx(j),j=1,5)
    end do 

    deallocate(coeff,indx)

  end subroutine table2


 subroutine mat_vec_multiply_SYM(ndim,A,vec,fin)

    integer, intent(in) :: ndim

    real(d), dimension(ndim,ndim), intent(in) :: A
    real(d), dimension(ndim), intent(in) :: vec
    real(d), dimension(ndim), intent(inout) :: fin

    integer :: LDA
    real(d) :: alpha    
    real(d) :: beta


    external dsymv

    write(6,*) "Starting matrix-vector multiplication"

    LDA=max(1,ndim)

    alpha=1._d
    beta=0._d


    call dsymv("L",ndim,alpha,A,LDA,vec,1,beta,fin,1)


  end subroutine mat_vec_multiply_SYM

 subroutine mat_vec_multiply(ndimiA,ndimjA,A,vec,fin)

    integer, intent(in) :: ndimiA
    integer, intent(in) :: ndimjA

    real(d), dimension(ndimiA,ndimjA), intent(in) :: A
    real(d), dimension(ndimjA), intent(in) :: vec
    real(d), dimension(ndimiA), intent(inout) :: fin

    integer :: LDA
    real(d) :: alpha    
    real(d) :: beta


    external dgemv

    write(6,*) "Starting matrix-vector multiplication"

    LDA=max(1,ndimiA)

    alpha=1._d
    beta=0._d


    call dgemv("N",ndimiA,ndimjA,alpha,A,LDA,vec,1,beta,fin,1)

  end subroutine mat_vec_multiply




  subroutine read_density_matrix


  integer ::  k,a
  integer :: i
  real*8, dimension(nBas,nOcc) :: density_matrix

  OPEN(unit=77,file="density.dat",status='OLD',access='SEQUENTIAL',form='FORMATTED')
  write(6,*) "nBas from read_density_matrix",nBas

  do i=1,nBas
  read(77) k,a,density_matrix(a,k)
  end do

  close(77)

  end subroutine read_density_matrix
 


end module misc
