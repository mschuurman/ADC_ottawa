! **********************************************************************
! * Taken and adapted from the Quantics code                           *
! **********************************************************************

module eigchessmod

contains
    
! **********************************************************************
! *                                                                    *
! *   EIGENVECTORS OF A COMPLEX UPPER HESSENBERG MATRIX (eigchess.f)   *
! *                                                                    *
! * Library module storing a routine that computes the eigenvalues and *
! * eigenvectors of a complex upper hessenberg matrix.                 *
! *                                                                    *
! * Contains:                                                          *
! *   eigchess: The eigenvalue routine. Taken from:                    *
! *             J. H. Wilkinson and C. Reinsch: Linear Algebra         *
! *             Routine: comlr2                                        *
! *                                                                    *
! **********************************************************************


! **********************************************************************
! *                                                                    *
! *                         SUBROUTINE EIGCHESS                        *
! *                                                                    *
! * Determines all eigenvalues and eigenvectors of a complex upper     *
! * Hessenberg matrix using an LR algorithm. The eigenvectors have the *
! * Euklidian norm one. All floating point parameters are double       *
! * precision (i. e. complex*16).                                      *
! *                                                                    *
! * Input parameters:                                                  *
! *   dim:       Leading dimension of h and v (for v see below).       *
! *   n:         True dimension of h and v.                            *
! *   h:         Upper Hessenberg matrix (destroyed on output).        *
! *   macheps:   Machine precision, i. e. the smallest positive number *
! *              for which 1.0+macheps > 1.0.                          *
! *                                                                    *
! * Output parameters:                                                 *
! *   w:         Vector of length n containing the eigenvalues.        *
! *   v:         Matrix giving the eigenvectors (stored columnwise).   *
! *   err:       Error flag having the following meaning:              *
! *              0: everything was o. k.,                              *
! *              1: the algorithm did not converge.                    *
! *                                                                    *
! * V7.0 MB                                                            *
! *                                                                    *
! **********************************************************************

  subroutine eigchess (dim,n,h,macheps,w,v,err)

    implicit none

    integer, parameter :: maxits=30

    integer            :: dim,n,err
    real*8             :: macheps
    complex*16         :: h(dim,n),w(n),v(dim,n)

    integer            :: i,j,k,m,its,en
    real*8             :: norm,tmp
    complex*16         :: s,t,x,y,z

!----------------------------------------------------------------------
! Initialise variables
!----------------------------------------------------------------------
    err = 0
    do j = 1,n
       do i = j+2,n
          h(i,j) = (0.0d0,0.0d0)
       enddo
    enddo
    do j = 1,n
       do i = 1,n
          v(i,j) = (0.0d0,0.0d0)
       enddo
       v(j,j) = (1.0d0,0.0d0)
       w(j)   = (0.0d0,0.0d0)
    enddo
    do j = n-1,2,-1
       do i = j+1,n
          v(i,j) = h(i,j-1)
       enddo
    enddo
    t = (0.0d0,0.0d0)
    en = n

!----------------------------------------------------------------------
! Begin search for next eigenvalue
!----------------------------------------------------------------------
100 continue
    if (en .lt. 1) goto 600
    its = 0

!----------------------------------------------------------------------
! Look for single small sub-diagonal element
!----------------------------------------------------------------------
200 continue
    do k = en,2,-1
       if (abs(dble(h(k,k-1)))+abs(dimag(h(k,k-1))) .le. &
            macheps*(abs(dble(h(k-1,k-1)))+abs(dimag(h(k-1,k-1))) &
            +abs(dble(h(k,k)))+abs(dimag(h(k,k))))) goto 300
    enddo
    k = 1

!----------------------------------------------------------------------
! One root found
!----------------------------------------------------------------------
300 continue
    if (k .eq. en) goto 500
      
! ----------------------------------------------------------------------
! Stop if algorithm does not converge
! ----------------------------------------------------------------------
    if (its .ge. maxits) then
       err = 1
       return
    endif

! ----------------------------------------------------------------------
! Form shift to accelerate convergence
! ----------------------------------------------------------------------
    if ((mod(its,10) .eq. 0) .and. (its .ne. 0)) then
       s = dcmplx(abs(dble(h(en,en-1)))+abs(dble(h(en-1,en-2))), &
            abs(dimag(h(en,en-1)))+abs(dimag(h(en-1,en-2))))
    else
       s = h(en,en)
       x = h(en-1,en)*h(en,en-1)
       if (x .ne. (0.0d0,0.0d0)) then
          y = 0.5d0*(h(en-1,en-1)-s)
          z = sqrt(y**2+x)
          if (dble(y)*dble(z)+dimag(y)*dimag(z) .lt. 0.0d0) z = -z
          x = x/(y+z)
          s = s-x
       endif
    endif
    do i = 1,en
       h(i,i) = h(i,i)-s
    enddo
    t = t+s
    its = its+1
    j = k+1

!----------------------------------------------------------------------
! Look for two consecutive small sub-diagonal elements
!----------------------------------------------------------------------
    x = dcmplx(abs(dble(h(en-1,en-1)))+abs(dimag(h(en-1,en-1))), &
         dimag(x))
    y = dcmplx(abs(dble(h(en,en-1)))+abs(dimag(h(en,en-1))),dimag(y))
    z = dcmplx(abs(dble(h(en,en)))+abs(dimag(h(en,en))),dimag(z))
    do m = en-1,j,-1
       y = dcmplx(abs(dble(h(m,m-1)))+abs(dimag(h(m,m-1))),dble(y))
       tmp = dble(x)
       x = dcmplx(abs(dble(h(m-1,m-1)))+abs(dimag(h(m-1,m-1))), &
            dble(z))
       z = dcmplx(tmp,dimag(z))
       if (dble(y) .le. macheps*dble(z)/dimag(y)*(dble(z)+dble(x) &
            +dimag(x))) goto 400
    enddo
    m = k

! ----------------------------------------------------------------------
! Triangular decomposition H = LR
! ----------------------------------------------------------------------
400 continue
    do i = m+1,en
       x = h(i-1,i-1)
       y = h(i,i-1)
       if (abs(dble(x))+abs(dimag(x)) .lt. abs(dble(y)+abs(dimag(y)))) &
            then
          do j = i-1,n
             z = h(i-1,j)
             h(i-1,j) = h(i,j)
             h(i,j) = z
          enddo
          z = x/y
          w(i) = dcmplx(1.0d0,dimag(w(i)))
       else
          z = y/x
          w(i) = dcmplx(-1.0d0,dimag(w(i)))
       endif
       h(i,i-1) = z
       do j = i,n
          h(i,j) = h(i,j)-z*h(i-1,j)
       enddo
    enddo
      
!----------------------------------------------------------------------
! Composition H = RL
!----------------------------------------------------------------------
    do j = m+1,en
       x = h(j,j-1)
       h(j,j-1) = (0.0d0,0.0d0)
       if (dble(w(j)) .gt. 0.0d0) then
          do i = 1,j
             z = h(i,j-1)
             h(i,j-1) = h(i,j)
             h(i,j) = z
          enddo
          do i = 1,n
             z = v(i,j-1)
             v(i,j-1) = v(i,j)
             v(i,j) = z
          enddo
       endif
       do i = 1,j
          h(i,j-1) = h(i,j-1)+x*h(i,j)
       enddo
       do i = 1,n
          v(i,j-1) = v(i,j-1)+x*v(i,j)
       enddo
    enddo
    goto 200

!----------------------------------------------------------------------
! One root found
!----------------------------------------------------------------------
500 continue
    w(en) = h(en,en)+t
    en = en-1
    goto 100

!----------------------------------------------------------------------
! All roots found
!----------------------------------------------------------------------
600 continue

!----------------------------------------------------------------------
! Compute norm
!----------------------------------------------------------------------
    norm = 0.0d0
    do i = 1,n
       norm = norm+abs(dble(w(i)))+abs(dimag(w(i)))
       do j = i+1,n
          norm = norm+abs(dble(h(i,j)))+abs(dimag(h(i,j)))
       enddo
    enddo
      
!----------------------------------------------------------------------
! Backsubstitute to find vectors of upper triangular form
!----------------------------------------------------------------------
    do en = n,2,-1
       x = w(en)
       do i = en-1,1,-1
          z = h(i,en)
          do j = i+1,en-1
             z = z+h(i,j)*h(j,en)
          enddo
          y = x-w(i)
          if (y .eq. (0.0d0,0.0d0)) y = dcmplx(macheps*norm,dimag(y))
          h(i,en) = z/y
       enddo
    enddo

!----------------------------------------------------------------------
! Multiply by transformation matrix to give vectors of original matrix
!----------------------------------------------------------------------
    do j = n,1,-1
       do i = 1,n
          z = v(i,j)
          do k = 1,j-1
             z = z+v(i,k)*h(k,j)
          enddo
          v(i,j) = z
       enddo
    enddo

!----------------------------------------------------------------------
! Normalise eigenvectors (Euklidian norm)
!----------------------------------------------------------------------
    do j = 1,n
       norm = 0.0d0
       do i = 1,n
          norm = norm+dble(dconjg(v(i,j))*v(i,j))
       enddo
       norm = 1.0d0/dsqrt(norm)
       do i = 1,n
          v(i,j) = norm*v(i,j)
       enddo
    enddo
    
  end subroutine eigchess

end module eigchessmod
