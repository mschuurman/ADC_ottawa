!######################################################################
! Routines to compute the unitary transformation matrix that
! maximally approximately diagonalises a set of matrices
!######################################################################

module simdiag

contains

!######################################################################
! simdiag_jacobi: Simulataneous diagonalisation of a set of matrices
!                 using a Jacobi rotation scheme.
!
!                 For a set of real symmetric matrices A, returns an
!                 orthogonal transformation matrix Q.
!
!                 Uses the algorithms reported in
!
!                 SIAM J, Matrix Anal. Appl. 14, 927 (1993)
!
!                 and
!
!                 SIAM J, Matrix Anal. Appl. 17, 161 (1996)
!######################################################################
    
  subroutine simdiag_jacobi(A1,Q,dim,nmat)

    use constants
    use iomod
    use channels
    
    implicit none

    integer                          :: dim,nmat
    integer                          :: i,j,k,l,m,n,error
    integer, parameter               :: maxit=100
    real(d), dimension(dim,dim,nmat) :: A1,A
    real(d), dimension(dim,dim)      :: Q,R,tmp
    real(d), dimension(dim)          :: lambda
    real(d), dimension(3*dim)        :: work
    real(d)                          :: off,c,s
    real(d)                          :: last
    real(d), parameter               :: thrsh=1e-6_d
    logical                          :: converged

    A=A1

!----------------------------------------------------------------------
! Ouput what we are doing
!----------------------------------------------------------------------
    write(ilog,'(/,2x,a)') 'Approximate simultaneous diagonalisation &
         of X, Y, and Z'
    write(ilog,'(2x,a)') 'Alorithm: Jacobi'
    
!----------------------------------------------------------------------
! Initialisation: set the orthogonal transformation matrix Q as the
! matrix that diagonalises one of the sum matrices in the set A
!----------------------------------------------------------------------
    ! Diagonalise sum_k A_k
    Q=A(:,:,1)+A(:,:,2)+A(:,:,3)
    call dsyev('V','U',dim,Q,dim,lambda,work,3*dim,error)    
    if (error.ne.0) then
       errmsg='Initialisation failure in subroutine simdiag_jacobi'
       call error_control
    endif
    
    ! Similarity transform all matrices in the set
    do k=1,nmat
       tmp=A(:,:,k)
       A(:,:,k)=matmul(transpose(Q),matmul(tmp,Q))
    enddo
    
    ! Compute and save the initial objective function value
    off=0.0d0
    do k=1,nmat
       do i=1,dim-1
          do j=i+1,dim
             off=off+A(i,j,k)**2
          enddo
       enddo
    enddo
    last=off

    ! Output our progress
    write(ilog,'(2x,a10,i3,a2,F11.7)') 'Iteration ',0,': ',off
    
!----------------------------------------------------------------------
! Perform the simultaneous diagonalisation
!----------------------------------------------------------------------
    converged=.false.

    do n=1,maxit

       ! Loop over unique pairs of off-diagonal elements
       do i=1,dim-1
          do j=i+1,dim

             ! Compute the elements of the rotation matrix R(i,j,c,s)
             call get_cs(i,j,A,c,s,dim,nmat)
             
             R=0.0d0
             do k=1,dim
                R(k,k)=1.0d0
             enddo
             R(i,i)=c
             R(j,j)=c
             R(i,j)=s
             R(j,i)=-s
             
             ! Update the transformation matrix Q
             tmp=Q
             Q=matmul(tmp,transpose(R))
             
             ! Similarity transform the matrices in the set A
             do k=1,nmat
                tmp=A(:,:,k)
                A(:,:,k)=matmul(R,matmul(tmp,transpose(R)))
             enddo
             
          enddo
       enddo
       
       ! Compute the objective function value
       off=0.0d0
       do k=1,nmat
          do i=1,dim-1
             do j=i+1,dim
                off=off+A(i,j,k)**2
             enddo
          enddo
       enddo
       
       ! Output our progress
       write(ilog,'(2x,a10,i3,a2,F11.7)') 'Iteration ',n,': ',off
       
       ! Exit if we have reached convergence
       if (abs(last-off).lt.thrsh) then
          converged=.true.
          exit
       else
          last=off
       endif
          
    enddo

!----------------------------------------------------------------------
! Exit if convergence was not reached
!----------------------------------------------------------------------
    if (.not.converged) then
       errmsg='Convergence not reaced in simdiag_jacobi'
       call error_control
    endif

    return
      
  end subroutine simdiag_jacobi

!######################################################################

  subroutine get_cs(i,j,A,c,s,dim,nmat)

    use constants
    use iomod

    integer                          :: i,j,dim,nmat
    integer                          :: k,p,q,error
    real(d)                          :: c,s
    real(d)                          :: x,y,r
    real(d), dimension(dim,dim,nmat) :: A
    real(d), dimension(2,nmat)       :: h
    real(d), dimension(2,2)          :: G,eigvec
    real(d), dimension(2)            :: eigval
    real(d), dimension(6)            :: work
    real(d)                          :: tiny

    real(d) :: theta,val
    real(d) :: t_on,t_off
  
!----------------------------------------------------------------------
! Cutoff for taking double precision values to be zero
!----------------------------------------------------------------------
    tiny=epsilon(tiny)

!----------------------------------------------------------------------
! TEST: Jacobi angle for a single matrix
!----------------------------------------------------------------------
    !! TEST: Jacobi angle for a single matrix
    !if (A(j,j,1)-A(i,i,1).ne.0.0d0) then
    !   val=2.0d0*A(i,j,1)/(A(j,j,1)-A(i,i,1))
    !   theta=atan(val)/2.0d0
    !else
    !   theta=pi/4.0d0
    !endif
    !c=cos(theta)
    !s=-sin(theta)
    !return
    !! TEST: Jacobi angle for a single matrix

!----------------------------------------------------------------------
! Calculate the h-vectors
!----------------------------------------------------------------------
    do k=1,nmat
       h(1,k)=A(i,i,k)-A(j,j,k)
       h(2,k)=2.0d0*A(i,j,k)
    enddo

!----------------------------------------------------------------------
! Calculate the G-matrix
!----------------------------------------------------------------------
    G=0.0d0
    do k=1,nmat
       do p=1,2
          do q=1,2
             G(p,q)=G(p,q)+h(p,k)*h(q,k)
          enddo
       enddo
    enddo

!----------------------------------------------------------------------
! Diagonalise the G-matrix
!----------------------------------------------------------------------
    eigvec=G

    call dsyev('V','U',2,eigvec,2,eigval,work,6,error)

    if (error.ne.0) then
       errmsg='Diagonalisation of the G-matrix failed in subroutine &
            get_cs'
       call error_control
    endif

!----------------------------------------------------------------------
! Calculate the cosine and sine of the rotation angle
! Note that: (1) We only consider inner rotations, which can be
!                achieved by ensuring that x>=0
!            (2) The constraint c**2 + s**2 = 1 has to hold
!----------------------------------------------------------------------    
!    print*,

    if (abs(eigvec(1,2)).ge.tiny) then
       !print*,"CASE I"
       x=abs(eigvec(1,2))
       y=eigvec(2,2)

       ! CHECK
       if (y.lt.1.0d0/sqrt(2.0d0)) y=-y
       ! CHECK

       r=sqrt(x**2+y**2)
       c=sqrt((x+r)/(2.0d0*r))
       s=y/(sqrt(2.0d0*r*(x+r)))
    else if (abs(eigvec(1,1)).ge.tiny) then
       !print*,"CASE II"
       x=abs(eigvec(1,1))
       y=eigvec(2,1)

       ! CHECK
       if (y.lt.1.0d0/sqrt(2.0d0)) y=-y
       ! CHECK

       r=sqrt(x**2+y**2)
       c=sqrt((x+r)/(2.0d0*r))
       s=y/(sqrt(2.0d0*r*(x+r)))
    else
       !print*,"CASE III"
       c=cos(pi/4.0d0)
       s=sin(pi/4.0d0)
    endif

    !! TEST: Single matrix case
    !print*,'ACTUAL: ',c,s
    !if (A(j,j,1)-A(i,i,1).ne.0.0d0) then
    !   val=2.0d0*A(i,j,1)/(A(j,j,1)-A(i,i,1))
    !   theta=atan(val)/2.0d0
    !else
    !   theta=pi/4.0d0
    !endif
    !!c=cos(theta)
    !!s=-sin(theta)
    !print*,'CORRECT:',cos(theta),-sin(theta)
    !! TEST: Single matrix case

    return
    
  end subroutine get_cs

!######################################################################

  function arctan2(x,y)

    use constants
    
    implicit none

    real(d) :: x,y,arctan2

    if (x.gt.0.0d0) then
       arctan2=atan(x/y)
    else if (x.lt.0.0d0.and.y.ge.0.0d0) then
       arctan2=atan(x/y)+pi
    else if (x.lt.0.0d0.and.y.lt.0.0d0) then
       arctan2=atan(x/y)-pi
    else if (x.eq.0.0d0.and.y.gt.0) then
       arctan2=0.5d0*pi
    else if (x.eq.0.0d0.and.y.lt.0) then
       arctan2=-0.5d0*pi
    else
       ! Undefined: x=0 and y=0, but compute this just so that NaN
       ! is returned
       arctan2=atan(x/y)
    endif
       
    return
    
  end function arctan2
    
!######################################################################
    
end module simdiag
