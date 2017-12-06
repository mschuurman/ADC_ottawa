!######################################################################
! Routines to compute the unitary transformation matrix that
! maximally approximately diagonalises a set of matrices
!######################################################################

module simdiag

contains

!######################################################################
! simdiag_jacobi: Simultaneous diagonalisation of a set of matrices
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
    use timingmod
    
    implicit none

    integer                          :: dim,nmat
    integer                          :: i,j,k,l,m,n,noff,error
    integer                          :: ks
    integer, parameter               :: maxit=1500
    real(d)                          :: thrsh
    real(d), dimension(dim,dim,nmat) :: A1,A
    real(d), dimension(dim,dim)      :: Q,tmp
    real(d), dimension(dim)          :: lambda
    real(d), dimension(3*dim)        :: work
    real(d)                          :: func,last,c,s,cs,sc
    real(d), dimension(nmat)         :: off
    real(d)                          :: tw1,tw2,tc1,tc2
    logical                          :: converged
    
!----------------------------------------------------------------------
! Ouput what we are doing
!----------------------------------------------------------------------
    write(ilog,'(2x,a)') 'Approximate simultaneous diagonalisation &
         of X, Y, and Z'
    write(ilog,'(2x,a,/)') 'Alorithm: Jacobi'

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call times(tw1,tc1)

!----------------------------------------------------------------------
! Convergence threshold
!----------------------------------------------------------------------
    thrsh=sqrt(epsilon(thrsh))

!----------------------------------------------------------------------
! Initialisation: set the orthogonal transformation matrix Q as the
! matrix that diagonalises one of the sum matrices in the set A
!----------------------------------------------------------------------
    ! Number of off-diagonal elements in the lower triangle
    noff=dim*(dim-1)/2

    ! Make a copy of the set of matrices to diagonalise
    A=A1
    
    ! Diagonalise sum_k A_k/nmat
    Q=0.0d0
    do k=1,nmat
       Q=Q+A(:,:,k)/nmat
    enddo
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
             off(k)=off(k)+A(j,i,k)**2/noff
          enddo
       enddo
    enddo
    func=sum(off(1:nmat))/nmat
    last=func

    ! Output our progress
    write(ilog,'(2x,a10,i3,a2,ES15.7)') &
         'Iteration ',0,': ',func
    
!----------------------------------------------------------------------
! Perform the simultaneous diagonalisation.
!
! We terminate when all Givens rotations in a sweep have sin(theta)'s
! smaller in magnitude than threshold.
!----------------------------------------------------------------------
    converged=.false.

    ! Loop over sweeps
    do n=1,maxit

       converged=.true.
       
       ! Loop over unique pairs of off-diagonal elements
       do i=1,dim-1
          do j=i+1,dim

             ! Compute cosine and sine of the rotation angle
             call get_cs(i,j,A,c,s,dim,nmat)

             if (abs(s).gt.thrsh) then
                converged=.false.
                
                ! Update the matrices in the set A
                do l=1,nmat
                   
                   ! Recalculate rows i and j
                   do k=1,dim
                      cs=c*A(i,k,l)+s*A(j,k,l)
                      sc=-s*A(i,k,l)+c*A(j,k,l)
                      A(i,k,l)=cs
                      A(j,k,l)=sc
                   enddo
                   ! New matrix A_{n+1} from A_{n}
                   do k=1,dim
                      cs=c*A(k,i,l)+s*A(k,j,l)
                      sc=-s*A(k,i,l)+c*A(k,j,l)
                      A(k,i,l)=cs
                      A(k,j,l)=sc
                   enddo
                enddo
                
                ! Update the matrix of eigenvectors
                tmp=Q
                Q(:,i)=c*tmp(:,i)+s*tmp(:,j)
                Q(:,j)=c*tmp(:,j)-s*tmp(:,i)

             endif
                
          enddo
       enddo

       ! Compute the objective function value
       off=0.0d0
       do k=1,nmat
          do i=1,dim-1
             do j=i+1,dim
                off(k)=off(k)+A(j,i,k)**2/noff
             enddo
          enddo
       enddo
       func=sum(off(1:nmat))/nmat
       
       ! Output our progress
       write(ilog,'(2x,a10,i3,a2,ES15.7,2x,ES15.7)') &
            'Iteration ',n,': ',func,abs(last-func)
       
       ! Exit if we have reached convergence
       if (converged) then
          exit
       else
          last=func
       endif
       
    enddo

!----------------------------------------------------------------------
! Exit if convergence was not reached
!----------------------------------------------------------------------
    if (.not.converged) then
       errmsg='Convergence not reaced in simdiag_jacobi'
       call error_control
    endif

!-----------------------------------------------------------------------    
! Output timings
!-----------------------------------------------------------------------    
    call times(tw2,tc2)
    write(ilog,'(/,a,1x,F11.4,1x,a)') 'Time taken:',tw2-tw1," s"

    return
      
  end subroutine simdiag_jacobi

!######################################################################

  subroutine get_cs(i,j,A,c,s,dim,nmat)

    use constants
    use iomod

    integer                          :: i,j,dim,nmat
    integer                          :: k,p,q
    real(d)                          :: theta,ton,toff
    real(d)                          :: c,s
    real(d), dimension(dim,dim,nmat) :: A
    real(d), dimension(2,nmat)       :: h
    real(d), dimension(2,2)          :: G

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
! Calculate the cosine and sine of the rotation angle
!----------------------------------------------------------------------
    ton=g(1,1)-g(2,2)
    toff=g(1,2)+g(2,1)

    theta=0.5d0*atan2(toff,ton+sqrt(ton*ton+toff*toff))

    c=cos(theta)
    s=sin(theta)

    return
    
  end subroutine get_cs
  
!######################################################################

end module simdiag
