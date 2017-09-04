! **********************************************************************
! * Linear equation solver library taken and adapted from the          *
! * Quantics code                                                      *
! **********************************************************************

     module lineqmod

      implicit none

      private
      public :: lineqd,lineqz,choleskyzhp,choleskydsp,matinvzhp,&
                lineqchd,precond

      integer, parameter :: dop=selected_real_kind(8)
      
      contains

! **********************************************************************
! *                                                                    *
! *                        LINEAR EQUATION (lineq.f)                   *
! *                                                                    *
! * Library module storing routines that solve linear equations.       *
! * There is little good reason to use these routines. LAPACK          *
! * equivalents should be used instead.                                *
! *                                                                    *
! * Contains:                                                          *
! *   lineqd:      Solves the real linear equation Ax = b employing    *
! *                the Gauss algorithm.                                *
! *   lineqz:      Solves the complex linear equation Ax = b employing *
! *                the Gauss algorithm.                                *
! *   choleskyzhp: Decomposition of a hermitian positive definite      *
! *                matrix.                                             *
! *   choleskydsp: Decomposition of a real symmetric positive definite *
! *                matrix.                                             *
! *   matinvzhp:   Inverts a hermitian positive definite matrix.       *
! *   lineqchd:    Solves the real linear equation Ax = b employing    *
! *                the Cholesky decomposition.                         *
! *                                                                    *
! *                                                                    *
! * Updated to F90 and OpenMP (Matinvzhp). GWR 09/13
! **********************************************************************

! ----------------------------------------------------------------------
!                                 LINEQD                         
!                                                                       
! Solves the linear equation Ax = b employing the Gauss algorithm. 
! All floating point parameters are double precision
!
!   Input parameters:
!     dim:       Leading dimension of A.
!     n:         True dimension of A, x and b.
!     a:         Coefficient matrix A (destroyed on output).
!     b:         Constant vector b (destroyed on output).
!
!   Output parameters:
!     x:         Solution vector.
!     err:       Error flag having the following meaning:
!                0: everything was o. k.,
!                1: the coefficient matrix is (numerically) singular.
! ----------------------------------------------------------------------

      subroutine lineqd (dim,n,a,b,x,err)

      implicit none

      integer,intent(in)                         :: dim,n
      integer, intent(out)                       :: err
      real(dop), dimension(dim,n), intent(inout) :: a
      real(dop), dimension(n), intent(inout)     :: b
      real(dop), dimension(n), intent(out)       :: x

      integer   :: pivrow,i,j,k
      real(dop) :: pivabs,pivabs2,pivot,swap,tmp

! ----------------------------------------------------------------------
! Initialise variables
! ----------------------------------------------------------------------
      err = 0
      k   = 1

! ----------------------------------------------------------------------
! Loop over all columns
! ----------------------------------------------------------------------
 100  continue

! ----------------------------------------------------------------------
! Find pivot element in current column
! ----------------------------------------------------------------------
      pivot  = a(k,k)
      pivrow = k
      pivabs = dble(abs(pivot))
      do i = k+1,n
         pivabs2 = dble(abs(a(i,k)))
         if (pivabs2 .gt. pivabs) then
            pivot  = a(i,k)
            pivrow = i
            pivabs = pivabs2
         endif
      enddo

! ----------------------------------------------------------------------
! Abort if matrix is singular
! ----------------------------------------------------------------------
      if (pivabs .eq. 0.0_dop) then
         err = 1
         return
      endif

! ----------------------------------------------------------------------
! Swap current and pivot row
! ----------------------------------------------------------------------
      if (k .ne. pivrow) then
         do j = k,n
            swap        = a(k,j)
            a(k,j)      = a(pivrow,j)
            a(pivrow,j) = swap
         enddo
         swap      = b(k)
         b(k)      = b(pivrow)
         b(pivrow) = swap
      endif

! ----------------------------------------------------------------------
! Eliminate subdiagonal elements in current column
! ----------------------------------------------------------------------
      do j = k+1,n
         tmp = a(k,j)/a(k,k)
         do i = k+1,n
            a(i,j) = a(i,j)-tmp*a(i,k)
         enddo
      enddo
      tmp = b(k)/a(k,k)
      do i = k+1,n
         b(i) = b(i)-tmp*a(i,k)
      enddo
      
! ----------------------------------------------------------------------
! Continue with next column
! ----------------------------------------------------------------------
      k = k+1
      if (k .le. n) goto 100

! ----------------------------------------------------------------------
! Determine solution vector
! ----------------------------------------------------------------------
      do i = n,1,-1
         x(i) = b(i)
         do j = i+1,n
            x(i) = x(i)-a(i,j)*x(j)
         enddo
         x(i) = x(i)/a(i,i)
      enddo

      return
      end subroutine lineqd 

! #######################################################################

! ----------------------------------------------------------------------
!                                 LINEQZ                         
!                                                                       
! Solves the complex linear equation Ax = b employing the Gauss
! algorithm. All floating point parameters are double precision
! (i. e. complex*16).
!
!   Input parameters:
!     dim:       Leading dimension of A.
!     n:         True dimension of A.
!     a:         Coefficient matrix A (destroyed on output).
!     b:         Constant vector b (destroyed on output).
!
!   Output parameters:
!     x:         Solution vector.
!     err:       Error flag having the following meaning:
!                0: everything was o. k.,
!                1: the coefficient matrix is (numerically) singular.
! ----------------------------------------------------------------------

      subroutine lineqz (dim,n,a,b,x,err)

      implicit none

      integer, intent(in)                           :: dim,n
      integer, intent(out)                          :: err
      complex(dop), dimension(dim,n), intent(inout) :: a
      complex(dop), dimension(n), intent(inout)     :: b
      complex(dop), dimension(n), intent(out)       :: x

      integer      :: pivrow,i,j,k
      real(dop)    :: pivabs,pivabs2
      complex(dop) :: pivot,swap,tmp

! ----------------------------------------------------------------------
! Initialise variables
! ----------------------------------------------------------------------
      err = 0
      k   = 1

! ----------------------------------------------------------------------
! Loop over all columns
! ----------------------------------------------------------------------
 100  continue

! ----------------------------------------------------------------------
! Find pivot element in current column
! ----------------------------------------------------------------------
      pivot  = a(k,k)
      pivrow = k
      pivabs = dble(abs(pivot))
      do i = k+1,n
         pivabs2 = dble(abs(a(i,k)))
         if (pivabs2 .gt. pivabs) then
            pivot  = a(i,k)
            pivrow = i
            pivabs = pivabs2
         endif
      enddo

! ----------------------------------------------------------------------
! Abort if matrix is singular
! ----------------------------------------------------------------------
      if (pivabs .eq. 0.0_dop) then
         err = 1
         return
      endif

! ----------------------------------------------------------------------
! Swap current and pivot row
! ----------------------------------------------------------------------
      if (k .ne. pivrow) then
         do j = k,n
            swap        = a(k,j)
            a(k,j)      = a(pivrow,j)
            a(pivrow,j) = swap
         enddo
         swap      = b(k)
         b(k)      = b(pivrow)
         b(pivrow) = swap
      endif

! ----------------------------------------------------------------------
! Eliminate subdiagonal elements in current column
! ----------------------------------------------------------------------
      do j = k+1,n
         tmp = a(k,j)/a(k,k)
         do i = k+1,n
            a(i,j) = a(i,j)-tmp*a(i,k)
         enddo
      enddo
      tmp = b(k)/a(k,k)
      do i = k+1,n
         b(i) = b(i)-tmp*a(i,k)
      enddo
      
! ----------------------------------------------------------------------
! Continue with next column
! ----------------------------------------------------------------------
      k = k+1
      if (k .le. n) goto 100

! ----------------------------------------------------------------------
! Determine solution vector
! ----------------------------------------------------------------------
      do i = n,1,-1
         x(i) = b(i)
         do j = i+1,n
            x(i) = x(i)-a(i,j)*x(j)
         enddo
         x(i) = x(i)/a(i,i)
      enddo

      return
      end subroutine lineqz

! #######################################################################

! ----------------------------------------------------------------------
!                                 LINEQZOMP                         
!                                                                       
! Solves the complex linear equation Ax = b employing the Gauss
! algorithm. All floating point parameters are double precision
! (i. e. complex*16).
!
!   Input parameters:
!     dim:       Leading dimension of A.
!     n:         True dimension of A.
!     a:         Coefficient matrix A (destroyed on output).
!     b:         Constant vector b (destroyed on output).
!     nthread:   Number of openmp threads.
!
!   Output parameters:
!     x:         Solution vector.
!     err:       Error flag having the following meaning:
!                0: everything was o. k.,
!                1: the coefficient matrix is (numerically) singular.
! ----------------------------------------------------------------------

      subroutine lineqzomp (dim,n,a,b,x,err,nthread)

!$use omp_lib

      implicit none

      integer, intent(in)                           :: dim,n,nthread
      integer, intent(out)                          :: err
      complex(dop), dimension(dim,n), intent(inout) :: a
      complex(dop), dimension(n), intent(inout)     :: b
      complex(dop), dimension(n), intent(out)       :: x

      integer      :: pivrow,i,j,k
      real(dop)    :: pivabs,pivabs2
      complex(dop) :: pivot,swap,tmp

! ----------------------------------------------------------------------
! Initialise variables
! ----------------------------------------------------------------------
      err = 0
      k   = 1

! ----------------------------------------------------------------------
! Loop over all columns
! ----------------------------------------------------------------------
 100  continue

! ----------------------------------------------------------------------
! Find pivot element in current column
! ----------------------------------------------------------------------
      pivot  = a(k,k)
      pivrow = k
      pivabs = dble(abs(pivot))
      do i = k+1,n
         pivabs2 = dble(abs(a(i,k)))
         if (pivabs2 .gt. pivabs) then
            pivot  = a(i,k)
            pivrow = i
            pivabs = pivabs2
         endif
      enddo

! ----------------------------------------------------------------------
! Abort if matrix is singular
! ----------------------------------------------------------------------
      if (pivabs .eq. 0.0_dop) then
         err = 1
         return
      endif

! ----------------------------------------------------------------------
! Swap current and pivot row
! ----------------------------------------------------------------------
      if (k .ne. pivrow) then
         do j = k,n
            swap        = a(k,j)
            a(k,j)      = a(pivrow,j)
            a(pivrow,j) = swap
         enddo
         swap      = b(k)
         b(k)      = b(pivrow)
         b(pivrow) = swap
      endif

! ----------------------------------------------------------------------
! Eliminate subdiagonal elements in current column
! ----------------------------------------------------------------------
!$omp parallel private(i,j,tmp) num_threads(nthread) 
!$omp do
      do j = k+1,n
         tmp = a(k,j)/a(k,k)
         do i = k+1,n
            a(i,j) = a(i,j)-tmp*a(i,k)
         enddo
      enddo
!$omp end do
      tmp = b(k)/a(k,k)
!$omp do
      do i = k+1,n
         b(i) = b(i)-tmp*a(i,k)
      enddo
!$omp end do
!$omp end parallel
      
! ----------------------------------------------------------------------
! Continue with next column
! ----------------------------------------------------------------------
      k = k+1
      if (k .le. n) goto 100

! ----------------------------------------------------------------------
! Determine solution vector
! ----------------------------------------------------------------------
      do i = n,1,-1
         x(i) = b(i)
         do j = i+1,n
            x(i) = x(i)-a(i,j)*x(j)
         enddo
         x(i) = x(i)/a(i,i)
      enddo

      return
      end subroutine lineqzomp

!-----------------------------------------------------------------------
!                               CHOLESKYZHP
!
! Decomposites a complex hermitian (H), positive definite (P) matrix
! A into a product of a left lower triangular matrix L and its adjoint,
! i. e. A = L L^dagger where L_ij = 0 for j > i and L_ii is real and
! positive. All calculations are performed with double precision (Z).
!
! Input parameters:
!   a:     The matrix to be decomposited (a(i,j) is needed for j <= i
!          only).
!   dim:   The number of allocated rows (and columns) of "a".
!   n:     The number of used rows (and columns) of "a".
!   tol:   A (relative) tolerance for the check wether "a" is hermitian
!          and positive definite.
!
! Output parameters:
!   a:     The left lower triangular matrix L. (The part of "a" right of
!          the diagonal still contains the corresponding part of A.)
!   error: An error flag having the following meaning:
!          0: Everything was o. k.
!          1: The matrix is not hermitian and positive definite.
!
! Needs OpenMP parallelisation once we can deal with tasks.
!-----------------------------------------------------------------------

      subroutine choleskyzhp (a,dim,n,tol,error)

      implicit none

      integer, intent(in)                             :: dim,n
      integer, intent(out)                            :: error
      real(dop), intent(in)                           :: tol
      complex(dop), dimension(dim,dim), intent(inout) :: a

      integer      :: i,j,k
      real(dop)    :: y
      complex(dop) :: x

      error = 0
      y = 0.0_dop
      do j = 1,n
         do i = j,n
            x = a(i,j)
            do k = 1,j-1
               x = x-a(i,k)*dconjg(a(j,k))
            enddo
            if (i .eq. j) then
               if (dabs(dimag(x)) .gt. tol*dble(x)) then
                  error = 1
                  return
               else
                  a(i,j) = dcmplx(dsqrt(dble(x)))
                  y = 1.0_dop/dble(a(i,j))
               endif
            else
               a(i,j) = x*y 
            endif
         enddo
      enddo

      return
      end subroutine choleskyzhp 

!-----------------------------------------------------------------------
!                               CHOLESKYDSP
!
! Decomposites a real symmetric (S), positive definite (P) matrix
! A into a product of a left lower triangular matrix L and its adjoint,
! i. e. A = L L^dagger where L_ij = 0 for j > i and L_ii is real and
! positive. All calculations are performed with double precision (D).
!
! Input parameters:
!   a:     The matrix to be decomposited (a(i,j) is needed for j <= i
!          only).
!   dim:   The number of allocated rows (and columns) of "a".
!   n:     The number of used rows (and columns) of "a".
!   tol:   A (relative) tolerance for the check wether "a" is symmetric
!          and positive definite.
!
! Output parameters:
!   a:     The left lower triangular matrix L. (The part of "a" right of
!          the diagonal still contains the corresponding part of A.)
!   error: An error flag having the following meaning:
!          0: Everything was o. k.
!          1: The matrix is not symmetric and positive definite.!
!
! Needs OpenMP parallelisation once we can deal with tasks.
!-----------------------------------------------------------------------

      subroutine choleskydsp (a,dim,n,error)

      implicit none

      integer,intent(in)                          :: dim,n
      integer, intent(out)                        :: error
      real(dop),dimension(dim,dim), intent(inout) :: a

      integer   :: i,j,k
      real(dop) :: x,y

      error = 0
      y = 0.0_dop
      do j = 1,n
         do i = j,n
            x = a(i,j)
            do k = 1,j-1
               x = x-a(i,k)*a(j,k)
            enddo
            if (i .eq. j) then
               if (x.le.0.0_dop) then
                  error=1
               else
                  a(i,j) = sqrt(x)
                  y = 1.0_dop/a(i,j) 
               endif
            else
               a(i,j) = x*y 
            endif
         enddo
      enddo

      return
      end subroutine choleskydsp

!-----------------------------------------------------------------------
!                               MATINVZHP
!
! Inverts a double precision complex (Z), Hermitian (H), and positive
! definite (P) matrix A employing a Cholesky-decomposition of A.
!
! Input parameters:
!   a:     The matrix to be inverted (a(i,j) is needed for j <= i only).
!   dim:   The number of allocated rows (and columns) of "a".
!   n:     The number of used rows (and columns) of "a".
!   tol:   A (relative) tolerance for the check whether "a" is Hermitian
!          and positive definite.
! Output parameters:
!   b:     The inverse of "a".
!   a:     The matrix "a" is overwritten.
!   error: An error flag having the following meaning:
!          0: Everything was o. k.
!          1: The matrix "a" is not hermitian and positive definite.
! Remarks:
!   Library routine "choleskyzhp" is used.
!-----------------------------------------------------------------------

      subroutine matinvzhp (b,a,dim,n,tol,error)

!!$ use omp_lib

      implicit none

      complex(dop), parameter :: zeroz=(0.0_dop,0.0_dop), &
                                 onez=(1.0_dop,0.0_dop)

      integer, intent(in)                             :: dim,n
      integer, intent(out)                            :: error
      real(dop),intent(in)                            :: tol
      complex(dop), dimension(dim,dim), intent(inout) :: a
      complex(dop), dimension(dim,dim), intent(out)   :: b

      integer      :: i,j,k
      complex(dop) :: x

! 1. Construct the A = L L^dagger Cholesky-decomposition.
! 2. For every column j of the identity matrix do:
!    a. Solve L b = 1 (first do-loop over i).
!    b. Solve L_dagger u = b (u is stored directly in the j-th column of
!       the matrix b) (second do-loop over i).

      call choleskyzhp(a,dim,n,tol,error)

      if (error .ne. 0) return
!!$omp parallel do if(lompthread) num_threads(ompthread) private(i,j,k,x)
      do j = 1,n
         do i = 1,n
            if (i .lt. j) then
               b(i,j) = zeroz
            elseif (i .eq. j) then
               b(i,j) = onez/a(i,i)
            else
               x = zeroz
               do k = j,i-1
                  x = x-a(i,k)*b(k,j)
               enddo
               b(i,j) = x/a(i,i)
            endif
         enddo

         do i = n,1,-1
            x = b(i,j)
            do k = i+1,n
               x = x-dconjg(a(k,i))*b(k,j)
            enddo
            b(i,j) = x/a(i,i)
         enddo
      enddo
!!$omp end parallel do
      return
      end subroutine matinvzhp

! ----------------------------------------------------------------------
!                                 LINEQCHD
!
! Solves the linear equation Ax = b employing the Cholesky
! decomposition. All floating point parameters are double precision
!
!   Input parameters:
!     dim:       Leading dimension of A.
!     n:         True dimension of A, x and b.
!     a:         Coefficient matrix A (destroyed on output).
!     b:         Constant vector b (destroyed on output), holds solution
!
!   Output parameters:
!     x:         auxiliary vector.
!     err:       Error flag having the following meaning:
!                0: everything was o. k.,
!                1: the coefficient matrix is (numerically) singular.
! ----------------------------------------------------------------------

      subroutine lineqchd (dim,n,a,c,b,err)

      implicit none

      integer, intent(in)                        :: dim,n
      integer, intent(out)                       :: err
      real(dop), dimension(dim,n), intent(inout) :: a
      real(dop), dimension(n), intent(inout)     :: b
      real(dop), dimension(n), intent(out)       :: c

      integer   :: i,j
      real(dop) :: x

! ----------------------------------------------------------------------
! Initialise variables
! ----------------------------------------------------------------------
      err = 0

!-----------------------------------------------------------------------
! Do the cholesky decomposition
!-----------------------------------------------------------------------
      call choleskydsp(a,dim,n,err)
      if (err .ne. 0) return

!-----------------------------------------------------------------------
! Solve system G*x`=b
!-----------------------------------------------------------------------
      do i=1,n
         x=c(i)
         do j=1,i-1
            x=x-a(i,j)*b(j)
         enddo
         b(i)=x/a(i,i)
      enddo

      do i=1,n
         x=0.0_dop
         do j=1,i
            x=x+b(j)*a(i,j)
         enddo
         c(i)=x
      enddo

! ----------------------------------------------------------------------
! Determine solution vector: Solve system G^T*x=x`
! ----------------------------------------------------------------------
      do i=0,n-1
         x=b(n-i)
         do j=n-i+1,n
            x=x-a(j,n-i)*c(j)
         enddo
         c(n-i)=x/a(n-i,n-i)
      enddo

      do i=1,n
         x=0.0_dop
         do j=i,n
            x=x+a(j,i)*c(j)
         enddo
         b(i)=x
      enddo

      return
      end subroutine lineqchd

!##########################################################################
! Routine to solve the linear equation Ax=b
! A is Hermitian complex and is diagonalised to form a preconditioner, P, 
! from the eigenvalues. The transformed equations PAx=Pb is then solved for
! x.
! GWR 03/16
!##########################################################################
      subroutine precond(dim,a,b,x,err,lomp,nthread,time)

use channels, only: ilog

      implicit none

      logical(kind=4), intent(in)                  :: lomp
      integer, intent(in)                          :: dim,nthread
      integer, intent(out)                         :: err
      complex(dop), dimension(dim,dim), intent(in) :: a
      complex(dop), dimension(dim), intent(in)     :: b
      complex(dop), dimension(dim), intent(out)    :: x

      integer                                      :: i,j
      complex(dop), dimension(dim)                 :: precondmat,btmp,tmpvec
      complex(dop), dimension(dim,dim)             :: atmp

      integer                                      :: maxidx,ninc
      real(dop)                                    :: maxrow
      complex(dop)                                 :: func

! test
      real(dop), intent(in) :: time
      real(dop)                                     :: cond1,cond2
      complex(dop), dimension(dim)                  :: ytmp
! test

      do i = 1,dim
         maxrow=0.0_dop
         maxidx = 1
         ninc = 0
         do j = 1,dim
            if(abs(a(i,j)).gt.maxrow)then
               ninc = ninc + 1
               maxidx = j
               maxrow = abs(a(i,j))
            endif
         enddo
         if(ninc.gt.0)then
            precondmat(i) = 1.0_dop/a(i,maxidx)
!            precondmat(i) = dcmplx(1.0_dop,0.0_dop)
         else
            precondmat(i) = dcmplx(1.0_dop,0.0_dop)
         endif
      enddo

      do i = 1, dim
         btmp(i) = b(i)*precondmat(i)
         do j = 1, dim
            atmp(i,j) = a(i,j)*precondmat(i)
         enddo
      enddo

      if(lomp)then
         call lineqzomp(dim,dim,atmp,btmp,x,err,nthread)
      else
         call lineqz(dim,dim,atmp,btmp,x,err)
      endif

      tmpvec = matmul(a,x)-b
      func = dconjg(tmpvec(1))*tmpvec(1)
      do i = 2,dim
         func = func + dconjg(tmpvec(i))*tmpvec(i)
      enddo        
      func = 0.5_dop*func

! test
!      tmpvec = -1.0_dop * tmpvec
      do i = dim,1,-1
         ytmp(i) = tmpvec(i)
         do j = i+1,dim
            ytmp(i) = ytmp(i)-atmp(i,j)*ytmp(j)
         enddo
         ytmp(i) = ytmp(i)/atmp(i,i)
      enddo
      cond1 = maxval(abs(ytmp(1:dim)))
      cond2 = maxval(abs(x(1:dim)))
      cond1 = cond1/cond2*1.0e16_dop
!      write(*,*)0.0241888*time,cond1
! test

      if(abs(func).gt.1.0e-8_dop)then
         write(ilog,*)'Function (precond)',abs(func)
         flush(ilog)
         err = 1
      endif

      return
      end subroutine precond

      end module lineqmod
