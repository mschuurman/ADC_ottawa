!#######################################################################
! blnczs: block-Lanczos eigensolver
!#######################################################################

 subroutine blnczs(id,nd,main,ncycles,maxmem,memx,mode,nprint,wthr,erange,unit,iparm,fparm,var1,var2)

   implicit none

   integer       :: id,nd,main,ncycles,maxmem,memx,mode,nprint,wthr,unit,iparm(:),var1,var2
   real          :: erange(:),fparm(:)
   
   return
 end subroutine blnczs

!#######################################################################
! dnvini: block-Davidson initialisation
!
! Currently written only for the case of in core storage of vectors
!#######################################################################

 subroutine dnvini(rstr,ndm)

   use parameters, only: dmain

   implicit none

   integer                    :: ndm
   real, dimension(ndm,dmain) :: invec
   logical                    :: rstr

!-----------------------------------------------------------------------
! If rstr=.true., read initial vectors from file
! If rstr=.false., take initial vectors to be a collection of 1h1p unit
! vectors
!-----------------------------------------------------------------------
   if (rstr) then
      write(6,'(/,a,/)') "Write the davidson restart code!"
      STOP
   else
      call dav_invecs(invec,ndm)
   endif

   return
 end subroutine dnvini

!#######################################################################
! dav_invecs: constructs the block-Davidson initial vectors - a set of
!             N=dmain 1h1p unit vectors
!#######################################################################

 subroutine dav_invecs(invec,ndm)
   
   use parameters, only: dmain

   implicit none

   integer                    :: ndm,i
   real, dimension(ndm,dmain) :: invec


   ! N.B., as we are using SLEPc for the Davidson diagonalisation, we
   ! really should be passing this information straight to
   ! EPSSetInitialSpace

   invec=0.0d0
   do i=1,dmain
      invec(i,i)=1.0d0
   enddo

   return

 end subroutine dav_invecs

!#######################################################################
! dinvop: block-Davidson eigensolver
!#######################################################################

 subroutine dinvop(ndm,dmain,mem,str,myb0,transp,rmtxhd,rmtxq1,rmtxhq1)

   implicit none

   integer       :: ndm,dmain,mem,myb0,transp,rmtxhd,rmtxq1,rmtxhq1
   character*4   :: str
   
   return
 end subroutine dinvop

!#######################################################################
! davidson_diag: SLEPc block-Davidson eigensolver
!#######################################################################

 subroutine davidson_diag(blckdim,matdim,davstates,davname)

   use constants

   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"

   integer                             :: blckdim,matdim,davstates,&
                                          maxbl,nrec,unit,num
   real(d), dimension(matdim)          :: hii
   double precision                    :: val
   double precision, dimension(matdim) :: vec
   character(36)                       :: davname

!----------------------------------------------------------------
! PETSc/SLEPc variables
!----------------------------------------------------------------
      Mat         ham           ! Matrix whose eigenpairs are sought
      EPS         eps           ! Eigenproblem solver context 
      EPSType     tname,method  ! Eigenproblem type and method name
      PetscReal   tol,error     ! Error tolerance and calculated error value
      PetscScalar kr,ki         ! Real and imaginary parts of an eigenvalue
      Vec         xr,xi         ! Real and imaginary parts of an eigenvector
      PetscInt    n,i,j,&       ! Integers...
                  Istart,&
                  Iend
      PetscInt    neig,maxit,&  ! No. eigenvalues sought, max no. iterations
                  its,nconv,&   ! ???, no. converged eigenpairs
                  ncv           ! max. dimension of the working space
      PetscInt    col(3)        ! Integer array used to asign values to
                                ! the matrix elements
      PetscInt    i1,i2,i3,&    ! More integers...
                  dim,nvecs
      PetscMPIInt rank          ! MPI thing (???)
      PetscErrorCode ierr       ! Error code
      PetscBool   flg           ! Logical flag

      PetscScalar elval         ! Real array used for assigning values
                                ! to the matrix elements      

      PetscScalar b1,b2         ! Bounds of a closesd interval in which
                                ! to find eigenvalues
      PetscScalar targ          ! Target value for the eigensolver


      Vec         ivec(blckdim) ! Initial vectors
      PetscInt    nivec         ! No. initial vectors
      PetscInt    indx          ! Integer array used to fill initial
                                ! vector array
      PetscScalar ftmp          ! Scalar array used to fill initial
                                ! vector array

!-----------------------------------------------------------------------
! Initialise SLEPc
!-----------------------------------------------------------------------
      call slepcinitialize(petsc_null_character,ierr)

!-----------------------------------------------------------------------
! Create the matrix ham corresponding to the ADC Hamiltonian matrix
!-----------------------------------------------------------------------
      n=matdim

      call matcreate(petsc_comm_world,ham,ierr)

      call matsetsizes(ham,petsc_decide,petsc_decide,n,n,ierr)

      call matsetfromoptions(ham,ierr)

      call matsetup(ham,ierr)

!-----------------------------------------------------------------------
! Set the values of the non-zero elements of the ADC Hamiltonian matrix
!
! N.B. here nrec and maxbl are read from fill_ondiag (via rdham_diag)
!      and subsequently passed to fill_offdiag
!
! To do: determine whether we actually have to pass the matrix elements
!        one-by-one to matsetvalues
!-----------------------------------------------------------------------
      ! (i) Diagonal elements
      call fill_ondiag(ham,maxbl,nrec,matdim)

      ! (ii) Off-diagonal elements      
      call fill_offdiag(maxbl,nrec,ham)

!-----------------------------------------------------------------------
! Create vectors xr and xi (real and imaginary parts of an eigenvector) 
! of dimensions compatible with the matrix ham
!-----------------------------------------------------------------------
      call matgetvecs(ham,xr,xi,ierr)

!-----------------------------------------------------------------------
! Create the eigensolver context
!-----------------------------------------------------------------------
      call epscreate(petsc_comm_world,eps,ierr)

!-----------------------------------------------------------------------
! Set operators
! N.B. Standard eigenvalue problem so one matrix only
!-----------------------------------------------------------------------
      call epssetoperators(eps,ham,petsc_null_object,ierr)

!-----------------------------------------------------------------------
! Set the problem type (Hermitian eigenvalue problem)
!-----------------------------------------------------------------------
      call epssetproblemtype(eps,eps_hep,ierr)

!-----------------------------------------------------------------------
! Set the eigensolver to be the Davidson method
!-----------------------------------------------------------------------
      call epssettype(eps,EPSGD,ierr)

!-----------------------------------------------------------------------
! Set the number of eigenvalues to be found (neig), the maximum
! dimension of the working subspace (ncv), and the maximum projected
! dimenson (mpd, which for serial execution we take to be equal to ncv)
!-----------------------------------------------------------------------
      neig=davstates
      ncv=blckdim
      call epssetdimensions(eps,neig,ncv,ncv,ierr)

!-----------------------------------------------------------------------
! Set the eigenpairs of interest: eigenvalues with the smallest
! magnitude
!-----------------------------------------------------------------------
      call epssetwhicheigenpairs(eps,EPS_SMALLEST_MAGNITUDE,ierr)

!-----------------------------------------------------------------------
! Set the initial vectors: Davidson diagonalisation is only ever used
! for the initial space, so the initial vectors will always correspond
! to a set of 1h1p unit vectors
!-----------------------------------------------------------------------
      nvecs=blckdim
      dim=matdim
      do i=1,nvecs
         ! Create the ith initial vector
         call veccreateseq(PETSC_COMM_SELF,dim,ivec(i),ierr)
         ! Assign the components of the ith initial vector
         ftmp=1.0d0
         indx=i
         call vecsetvalues(ivec(i),1,indx,ftmp,INSERT_VALUES,ierr)
         call vecassemblybegin(ivec(i),ierr)
         call vecassemblyend(ivec(i),ierr)
      enddo

      ! Set the initial vector space
      call epssetinitialspace(eps,nivec,ivec,ierr)

!-----------------------------------------------------------------------
! Solve the eigenproblem that we have set up
!-----------------------------------------------------------------------
      call epssolve(eps,ierr)

!-----------------------------------------------------------------------
! Get the number of converged eigenpairs
! Exit if we have not converged enough eigenpairs
!-----------------------------------------------------------------------
      call epsgetconverged(eps,nconv,ierr)
      
      if (nconv.lt.davstates) then
         write(6,'(/,2x,i3,1x,a,/)') davstates,'Davidson states requested'
         write(6,'(/,2x,i3,1x,a,/)') nconv,'Davidson states converged'
         STOP
      endif

!-----------------------------------------------------------------------
! Output the converged eigenpairs
!-----------------------------------------------------------------------
      ! Open file
      unit=77

!      open(unit=unit,file='dav_vecs.init',status='old',&
!           access='sequential',form='unformatted')

      open(unit=unit,file='dav_vecs.init',status='unknown',&
           access='sequential',form='unformatted')

      do i=0,nconv-1
         
         ! Get the ith converged eigenpair
         call epsgeteigenpair(eps,i,kr,ki,xr,xi,ierr)
         
         ! Write the ith converged eigenpair to file
         val=PetscRealPart(kr)
         vec=PetscRealPart(xr)
         num=i

         write(unit) num,val,vec(:)
         
      enddo

      ! Close file
      close(unit)

   return

 end subroutine davidson_diag

!#######################################################################
! fill_ondiag: Assembles the on-diagonal elements of the PETSc version
!              of the ADC Hamiltonian matrix
!#######################################################################

 subroutine fill_ondiag(ham,maxbl,nrec,matdim)

   use constants

   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"

   integer                    :: matdim,maxbl,nrec
   real(d), dimension(matdim) :: hii

!----------------------------------------------------------------
! PETSc/SLEPc variables
!----------------------------------------------------------------
   Mat            ham
   PetscErrorCode ierr
   PetscInt       dim1,dim2,indx
   PetscScalar    elval
   
!----------------------------------------------------------------
! Read the diagonal elements of the Hamiltonian matrix
!----------------------------------------------------------------
   call rdham_diag(hii,matdim,maxbl,nrec)

!----------------------------------------------------------------
! Pass the diagonal elements of the Hamiltonian matrix to
! matsetvalues
!----------------------------------------------------------------
   dim1=1
   dim2=1
   do indx=1,matdim
      elval=hii(indx)
      call matsetvalues(ham,dim1,indx,dim2,indx,elval,insert_values,ierr)
   enddo

 end subroutine fill_ondiag

!#######################################################################
! rdham_diag: reads the on-diagonal Hamiltonian matrix elements from 
!             file
!#######################################################################

 subroutine rdham_diag(hii,matdim,maxbl,nrec)

   use constants

   implicit none

   integer                    :: matdim,maxbl,nrec
   integer*8                  :: unit
   real(d), dimension(matdim) :: hii

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
   unit=77
   open(unit,file='hmlt.diai',status='old',access='sequential',&
        form='unformatted')

!-----------------------------------------------------------------------
! Read the on-diagonal Hamiltonian matrix elements
!-----------------------------------------------------------------------
   read(77) maxbl,nrec
   read(77) hii

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
   close(unit)

   return

 end subroutine rdham_diag

!#######################################################################
! fill_offdiag: Assembles the off-diagonal elements of the PETSc version
!               of the ADC Hamiltonian matrix
!#######################################################################

 subroutine fill_offdiag(maxbl,nrec,ham)

   use constants

   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"

   integer                   :: maxbl,nrec,nlim,k,l
   integer*8                 :: unit
   real(d), dimension(maxbl) :: hij
   integer, dimension(maxbl) :: indxi,indxj

!----------------------------------------------------------------
! PETSc/SLEPc variables
!----------------------------------------------------------------
   Mat            ham
   PetscErrorCode ierr
   PetscInt       dim1,dim2,i1,i2
   PetscScalar    elval

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
   unit=78
   open(unit,file='hmlt.offi',status='old',access='sequential',&
        form='unformatted')

!-----------------------------------------------------------------------
! Read the on-diagonal Hamiltonian matrix elements and pass to
! matsetvalues
!-----------------------------------------------------------------------
   dim1=1
   dim2=1

   do k=1,nrec
      read(unit) hij(:),indxi(:),indxj(:),nlim
      do l=1,nlim
         i1=indxi(l)
         i2=indxj(l)
         elval=hij(l)
         call matsetvalues(ham,dim1,i1,dim2,i2,elval,insert_values,ierr)
      enddo
   end do

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
   close(unit)

   return

 end subroutine fill_offdiag

!#######################################################################
