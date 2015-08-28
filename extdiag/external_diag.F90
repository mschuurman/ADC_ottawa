!#######################################################################
! eigensolver: SLEPc block-Davidson and Krylov-Schur eigensolvers
!#######################################################################

 subroutine eigensolver(blckdim,matdim,noffd)

   use constants
   use channels
   use parameters   
   
   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"
#include "finclude/petscvec.h90"
   
   integer                             :: blckdim,matdim,nstates,&
                                          maxbl,nrec,iout,count,kk
   integer*8                           :: noffd
   real(d), dimension(matdim)          :: hii
   double precision                    :: val
   double precision, dimension(matdim) :: vec
   character(36)                       :: vecfile
   character(70)                       :: compiler
   character(120)                      :: atmp
   logical                             :: lrdadc1,lincore
   logical, external                   :: isincore

!-----------------------------------------------------------------------
! PETSc/SLEPc variables
!-----------------------------------------------------------------------
      Mat         ham           ! Matrix whose eigenpairs are sought
      EPS         eps           ! Eigenproblem solver context 
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

      PetscScalar targ          ! Target value for the eigensolver


      Vec         ivec(blckdim)     ! Initial vectors
      PetscInt    nivec             ! No. initial vectors
      PetscInt    indx              ! Integer array used to fill initial
                                    ! vector array
      PetscScalar ftmp              ! Scalar array used to fill initial
                                    ! vector array
      PetscScalar, pointer :: xx(:) ! array to hold a single retrieved 
                                    ! eigenvector

      PetscInt nz,nnz(matdim)       ! Number of non-zero elements per row

      PetscInt  niter               ! Max. no. iterations

      external hxvec

!-----------------------------------------------------------------------
! Write to the log file
!-----------------------------------------------------------------------
      write(ilog,'(/,70a)') ('-',kk=1,70)

      if (solver.eq.1) then
         atmp='Block Davidson diagonalisation in the'
      else if (solver.eq.2) then
         atmp='Krylov-Schur diagonalisation in the'
      endif

      if (hamflag.eq.'i') then
         atmp=trim(atmp)//' initial space'
       else if (hamflag.eq.'f') then
          atmp=trim(atmp)//' final space'
      endif

      write(ilog,'(2x,a)') trim(atmp)
      
      write(ilog,'(70a,/)') ('-',kk=1,70)

!-----------------------------------------------------------------------
! If this is the first call to davidson_diag, then call mpi_init.
! It is crucial that we call mpi_init once before slepcinitialize is
! called for the first time, otherwise we cannot make multiple calls
! to slepcinitialize/slepcfinalize.
!
! Note that this only holds if the code has NOT been compiled using
! ifort
!-----------------------------------------------------------------------
#ifdef __INTEL_COMPILER
#define compiler "intel"
#endif

      ndavcalls=ndavcalls+1
      if (ndavcalls.eq.1.and.compiler.ne.'intel') call mpi_init()

!-----------------------------------------------------------------------
! Set the no. sates sought, the block size, etc
!-----------------------------------------------------------------------
      if (hamflag.eq.'i') then
         nstates=davstates
         lrdadc1=ladc1guess
         vecfile=davname
      else if (hamflag.eq.'f') then
         nstates=davstates_f
         lrdadc1=ladc1guess_f
         vecfile=davname_f
      endif
      
!-----------------------------------------------------------------------
! Initialise SLEPc
!-----------------------------------------------------------------------
      call slepcinitialize(petsc_null_character,ierr)

!-----------------------------------------------------------------------
! MPI
!-----------------------------------------------------------------------
!      call mpi_comm_rank(petsc_comm_world,rank,ierr)

!-----------------------------------------------------------------------
! Determine whether we have enough memory to store the Hamiltonian
! matrix in-core
!-----------------------------------------------------------------------
      lincore=isincore(matdim,noffd)

      if (lincore) then
         write(ilog,'(2x,a,/)') &
              'The Hamiltonian matrix will be stored in-core'
      else
         write(ilog,'(2x,a,/)') &
              'The Hamiltonian matrix will be stored externally'
      endif

!-----------------------------------------------------------------------
! If the Hamiltonian matrix is being stored in external memory then
! register the matrix-vector subroutine for the operator that defines
! the eigensystem, Ax=kx
!-----------------------------------------------------------------------
      if (.not.lincore) then

         n=matdim

         call MatCreateShell(MPI_COMM_WORLD,n,n,n,n,PETSC_NULL_OBJECT,&
              ham,ierr)

         call MatShellSetOperation(ham,MATOP_MULT,hxvec,ierr)

      endif

!-----------------------------------------------------------------------
! If the Hamiltonian matrix is being stored in-core then determine 
! the no. non-zero elements per row: important in order to not have 
! terrible performance from matsetvalues.
!
! We here store the no. non-zero elements per row in the PetscInt
! array nnz: nnz(i)=no. non-zero elements in the ith row.
!-----------------------------------------------------------------------
      if (lincore) call get_nonzeros(nnz,matdim,hamflag)

!-----------------------------------------------------------------------
! If the Hamiltonian matrix is being stored in-core then create the 
! matrix ham corresponding to the ADC Hamiltonian matrix
!-----------------------------------------------------------------------
      n=matdim

      if (lincore) then
         
         call matcreateseqaij(MPI_COMM_WORLD,n,n,nz,nnz,ham,ierr)
         
         call matsetfromoptions(ham,ierr)
         
         call matsetup(ham,ierr)

      endif
         
!-----------------------------------------------------------------------
! If the Hamiltonian matrix is being stored in-core then set the
! values of the non-zero elements of the ADC Hamiltonian matrix
!
! Here nrec and maxbl are read from fill_ondiag (via rdham_diag) and
! subsequently passed to fill_offdiag
!-----------------------------------------------------------------------
      if (lincore) then
         
         ! (i) Diagonal elements
         call fill_ondiag(ham,maxbl,nrec,matdim,hamflag)
         
         ! (ii) Off-diagonal elements      
         call fill_offdiag(maxbl,nrec,ham,hamflag)
         
         ! Assemble the PETSc matrix
         call matassemblybegin(ham,mat_final_assembly,ierr)
         call matassemblyend(ham,mat_final_assembly,ierr)

      endif
         
!-----------------------------------------------------------------------
! Create vectors xr and xi (real and imaginary parts of an eigenvector) 
! of dimensions compatible with the matrix ham
!-----------------------------------------------------------------------
      call matgetvecs(ham,xr,xi,ierr)
      
!-----------------------------------------------------------------------
! Create the eigensolver context
!-----------------------------------------------------------------------
      call epscreate(mpi_comm_world,eps,ierr)
      
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
! Set the eigensolver
!-----------------------------------------------------------------------
      if (solver.eq.1) then
         call epssettype(eps,EPSGD,ierr)
      else if (solver.eq.2) then
         call epssettype(eps,EPSKRYLOVSCHUR,ierr)
      endif

!-----------------------------------------------------------------------
! Set the number of eigenvalues to be found (neig), the maximum
! dimension of the working subspace (ncv), and the maximum projected
! dimenson (mpd, which for serial execution we take to be equal to ncv)
!-----------------------------------------------------------------------
      neig=nstates
      ncv=min(nstates+blckdim+10,(nstates+blckdim)*2)
     
      call epssetdimensions(eps,neig,ncv,ncv,ierr)

!-----------------------------------------------------------------------
! Set the eigenpairs of interest: either eigenvalues with the smallest
! real part or those eigenvalues closest to a target values
!-----------------------------------------------------------------------
      if (eigentype.eq.1) then
         call epssetwhicheigenpairs(eps,EPS_SMALLEST_REAL,ierr)
      else if (eigentype.eq.2) then
         call epssetwhicheigenpairs(eps,EPS_TARGET_MAGNITUDE,ierr)
         targ=davtarg
         call epssettarget(eps,targ,ierr)
      endif
         
!-----------------------------------------------------------------------
! Set the error tolerance and max. no. iterations
!-----------------------------------------------------------------------
      if (hamflag.eq.'i') then
         tol=davtol
         niter=maxiter
      else if (hamflag.eq.'f') then
         tol=davtol_f
         niter=maxiter_f
      endif
         
      call epssettolerances(eps,tol,niter,ierr)

!-----------------------------------------------------------------------
! Set the initial vectors
!-----------------------------------------------------------------------
      if (lrdadc1) then
         call load_adc1_vecs(blckdim,matdim,ivec)
      else
         call guess_vecs_ondiag(blckdim,matdim,ivec,hamflag)
      endif

      ! Set the initial vector space
      nivec=blckdim
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

      ! Output the number of iterations taken
      call epsgetiterationnumber(eps,its,ierr)
      call epsgettolerances(eps,tol,maxit,ierr)
      write(ilog,'(/,x,a,x,i4,/)') 'Number of iterations:',its

      ! If not all states have been converged, then exit
      if (nconv.lt.nstates) then
         write(ilog,'(/,2x,i3,1x,a,/)') nstates,'Davidson states requested'
         write(ilog,'(/,2x,i3,1x,a,/)') nconv,'Davidson states converged'
         STOP
      endif

!-----------------------------------------------------------------------
! Output the converged eigenpairs
!-----------------------------------------------------------------------
      ! Open file
      iout=77

      open(unit=iout,file=vecfile,status='unknown',&
           access='sequential',form='unformatted')
      
      do i=0,nconv-1
         count=i+1
         ! Get the ith converged eigenpair
         call epsgeteigenpair(eps,i,kr,ki,xr,xi,ierr) 
         ! Write the ith converged eigenpair to file
         val=PetscRealPart(kr)
         call VecGetArrayF90(xr,xx,ierr) 
         vec=xx
         write(iout) count,val,vec(:)
         ! Calculate the relative error for the ith eigenpair
         call epscomputerelativeerror(eps,i,error,ierr)
      enddo

      ! Close file
      close(iout)

!-----------------------------------------------------------------------
! Free the workspace
!-----------------------------------------------------------------------
      call epsdestroy(eps,ierr)
      call matdestroy(ham,ierr)

!      call vecdestroy(xr,ierr)
!      call vecdestroy(xi,ierr)
!      call vecdestroy(xx,ierr)

!-----------------------------------------------------------------------
! Finalise SLEPc
!-----------------------------------------------------------------------
      call slepcfinalize(ierr)

! Full diagonalisation for debugging purposes
!      call full_diag(matdim)

   return

 end subroutine eigensolver

!#######################################################################
! fill_ondiag: Assembles the on-diagonal elements of the PETSc version
!              of the ADC Hamiltonian matrix
!#######################################################################

 subroutine fill_ondiag(ham,maxbl,nrec,matdim,flag)

   use constants

   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"

   integer                    :: matdim,maxbl,nrec
   real(d), dimension(matdim) :: hii
   character(1)               :: flag

!-----------------------------------------------------------------------
! PETSc/SLEPc variables
!-----------------------------------------------------------------------
   Mat            ham
   PetscErrorCode ierr
   PetscInt       dim1,dim2,indx,pmatindx
   PetscScalar    elval
   
!-----------------------------------------------------------------------
! Read the diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
   call rdham_diag(hii,matdim,maxbl,nrec,flag)

!-----------------------------------------------------------------------
! Pass the diagonal elements of the Hamiltonian matrix to
! matsetvalues
!-----------------------------------------------------------------------
   dim1=1
   dim2=1
   do indx=1,matdim
      elval=hii(indx)
      ! PETSc matrix indices start from zero
      pmatindx=indx-1
      call matsetvalues(ham,dim1,pmatindx,dim2,pmatindx,elval,&
           insert_values,ierr)
   enddo
   
 end subroutine fill_ondiag

!#######################################################################
! rdham_diag: reads the on-diagonal Hamiltonian matrix elements from 
!             file
!#######################################################################

 subroutine rdham_diag(hii,matdim,maxbl,nrec,flag)

   use constants

   implicit none

   integer                    :: matdim,maxbl,nrec
   integer*8                  :: unit
   real(d), dimension(matdim) :: hii
   character(1)               :: flag
   character(70)              :: filename

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
   unit=77

   if (flag.eq.'i') then
      filename='SCRATCH/hmlt.diai'
   else if (flag.eq.'f') then
      filename='SCRATCH/hmlt.diac'
   endif

   open(unit,file=filename,status='old',access='sequential',&
        form='unformatted')

!-----------------------------------------------------------------------
! Read the on-diagonal Hamiltonian matrix elements
!-----------------------------------------------------------------------
   rewind(unit)
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

 subroutine fill_offdiag(maxbl,nrec,ham,flag)

   use constants

   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"

   integer                   :: maxbl,nrec,nlim,k,l
   integer*8                 :: unit
   integer, dimension(maxbl) :: indxi,indxj
   real(d), dimension(maxbl) :: hij
   character(1)              :: flag
   character(70)             :: filename

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

   if (flag.eq.'i') then
      filename='SCRATCH/hmlt.offi'
   else if (flag.eq.'f') then
      filename='SCRATCH/hmlt.offc'
   endif

   open(unit,file=filename,status='old',access='sequential',&
        form='unformatted')

!-----------------------------------------------------------------------
! Read the off-diagonal Hamiltonian matrix elements and pass to
! matsetvalues
!-----------------------------------------------------------------------
   dim1=1
   dim2=1
   do k=1,nrec     
      read(unit) hij(:),indxi(:),indxj(:),nlim
      do l=1,nlim
         ! N.B., PETSc matrix indices start from zero
         i1=indxi(l)-1
         i2=indxj(l)-1
         elval=hij(l)
         call matsetvalues(ham,dim1,i1,dim2,i2,elval,insert_values,ierr)
         i2=indxi(l)-1
         i1=indxj(l)-1
         elval=hij(l)
         call matsetvalues(ham,dim1,i1,dim2,i2,elval,insert_values,ierr)
      enddo
   enddo

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
   close(unit)

   return

 end subroutine fill_offdiag

!#######################################################################
! full_diag: full diagonalisation of the Hamiltonian matrix using
!            the LAPACK routine dsyev
!#######################################################################

 subroutine full_diag(matdim)

   use constants
   use parameters
   use channels

   implicit none

   integer                            :: matdim,iham,maxbl,nrec,nlim,&
                                         i,k,i1,i2
   real(d), dimension(matdim)         :: hii
   real(d), dimension(matdim,matdim)  :: hmat,umat
   real(d), dimension(:), allocatable :: hij
   integer, dimension(:), allocatable :: indxi,indxj
   character(70)                      :: filename

   integer*4                   :: e2
   real*8, dimension(3*matdim) :: work
   real*8                      :: error
   real*8, dimension(matdim)   :: eigval

   hmat=0.0d0

!-----------------------------------------------------------------------
! Read the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
   iham=77

   if (hamflag.eq.'i') then
      filename='SCRATCH/hmlt.diai'
   else if (hamflag.eq.'f') then
      filename='SCRATCH/hmlt.diac'
   endif

   open(iham,file=filename,status='old',access='sequential',&
        form='unformatted')

   rewind(iham)
   read(iham) maxbl,nrec
   read(iham) hii

   do i=1,matdim
      hmat(i,i)=hii(i)
   enddo

   close(iham)

!-----------------------------------------------------------------------
! Read the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
   if (hamflag.eq.'i') then
      filename='SCRATCH/hmlt.offi'
   else if (hamflag.eq.'f') then
      filename='SCRATCH/hmlt.offc'
   endif

   open(iham,file=filename,status='old',access='sequential',&
        form='unformatted')

   allocate(hij(maxbl))
   allocate(indxi(maxbl))
   allocate(indxj(maxbl))

   do i=1,nrec
      read(iham) hij(:),indxi(:),indxj(:),nlim
      do k=1,nlim
         i1=indxi(k)
         i2=indxj(k)
         hmat(i1,i2)=hij(k)
         hmat(i2,i1)=hij(k)
      enddo
   enddo
   
   close(iham)

!-----------------------------------------------------------------------
! Diagonalise the Hamiltonian matrix
!-----------------------------------------------------------------------
   error=0
   e2=3*matdim
   umat=hmat
   call dsyev('V','U',matdim,umat,matdim,eigval,work,e2,error)

!-----------------------------------------------------------------------
! Exit here if the diagonalisation was not successful
!-----------------------------------------------------------------------
   if (error.ne.0) then
      write(ilog,'(/,a,/)') &
           'Diagonalisation in subroutine full_diag failed.'
      STOP
   endif

   do i=1,15
      print*,i,eigval(i)*eh2ev
   enddo

   return

 end subroutine full_diag

!#######################################################################

 subroutine get_nonzeros(nnz,matdim,flag)

   use constants

   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"
#include "finclude/petscvec.h90"

   integer                            :: matdim,maxbl,nrec,nlim,k,l,&
                                         tot,tot0
   integer*8                          :: unit
   integer, dimension(:), allocatable :: indxi,indxj
   real(d), dimension(:), allocatable :: hij
   character(1)                       :: flag
   character(len=70)                  :: filename

!-----------------------------------------------------------------------
! PETSc/SLEPc variables
!-----------------------------------------------------------------------
   PetscInt nnz(matdim)       ! Number of non-zero elements per row

!-----------------------------------------------------------------------
! Read buffer size and no. records from hmlt.diai
!-----------------------------------------------------------------------
   unit=78

   if (flag.eq.'i') then
      filename='SCRATCH/hmlt.diai'
   else if (flag.eq.'f') then
      filename='SCRATCH/hmlt.diac'
   endif

   open(unit,file=filename,status='old',access='sequential',&
        form='unformatted')

   read(unit) maxbl,nrec

   close(unit)

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
   allocate(hij(maxbl))
   allocate(indxi(maxbl))
   allocate(indxj(maxbl))

!-----------------------------------------------------------------------
! Read the off-diagonal Hamiltonian matrix elements and determine the
! no. non-zero elements per row (nnz)
!-----------------------------------------------------------------------
   if (flag.eq.'i') then
      filename='SCRATCH/hmlt.offi'
   else if (flag.eq.'f') then
      filename='SCRATCH/hmlt.offc'
   endif

   open(unit,file=filename,status='old',access='sequential',&
        form='unformatted')

   nnz=1
   do k=1,nrec
      read(unit) hij(:),indxi(:),indxj(:),nlim
      do l=1,nlim
         nnz(indxi(l))=nnz(indxi(l))+1
         nnz(indxj(l))=nnz(indxj(l))+1
      enddo
   enddo

   close(unit)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
   deallocate(hij)
   deallocate(indxi)
   deallocate(indxj)

   return

 end subroutine get_nonzeros

!#######################################################################

 subroutine load_adc1_vecs(blckdim,matdim,ivec)

   use constants

   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"
#include "finclude/petscvec.h90"

   integer                              :: blckdim,matdim,unit,dim1,&
                                           curr
   integer, dimension(:), allocatable   :: indx1
   real(d), dimension(:,:), allocatable :: vec1

!-----------------------------------------------------------------------
! PETSc/SLEPc variables
!-----------------------------------------------------------------------
   Vec            ivec(blckdim)
   PetscInt       i,nvecs,dim2
   PetscScalar    ftmp
   PetscErrorCode ierr

!-----------------------------------------------------------------------
! Open ADC(1) vector file
!-----------------------------------------------------------------------
   unit=225
   open(unit,file='SCRATCH/adc1_vecs',form='unformatted',status='old')

!-----------------------------------------------------------------------
! Read ADC(1) eigenvectors and set initial Davidson vectors
!-----------------------------------------------------------------------
   read(unit) dim1
   allocate(vec1(dim1,dim1))
   allocate(indx1(dim1))

   rewind(unit)
   read(unit) dim1,vec1

!-----------------------------------------------------------------------
! Set the initial Davidson vectors
!-----------------------------------------------------------------------

   ! Set the indices of the elements to be inserted into the PETSc
   ! vectors, taking into account that PETSc vector indices start from 
   ! zero
   do i=1,dim1
      indx1(i)=i-1
   enddo

   ! Set the PETSc integers corresponding to the no. starting vectors
   ! and the dimension of these vectors
   nvecs=blckdim
   dim2=matdim

   ! Loop over the starting vectors, such that the ith vector
   ! corresponds to a direct sum of the ith ADC(1) eigenvector and
   ! a zero vector in the space spanned by the 2h2p states
   curr=0
   do i=1,nvecs
      curr=curr+1
      ! Create the ith initial vector (of dimension dim2=matdim)
!      call veccreateseq(PETSC_COMM_SELF,dim2,ivec(i),ierr)
      call veccreateseq(MPI_COMM_SELF,dim2,ivec(i),ierr)
      ! Assign the components of the ith initial vector
      call setvec(i,ivec,dim1,vec1(:,curr),blckdim,indx1)
   enddo

!-----------------------------------------------------------------------
! Close ADC(1) vector file
!-----------------------------------------------------------------------
   close(unit)

   return

 end subroutine load_adc1_vecs

!#######################################################################

 subroutine setvec(i,ivec,dim1,vec,blckdim,indx1)

   use constants

   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"
#include "finclude/petscvec.h90"

   integer                  :: dim1,blckdim
   integer, dimension(dim1) :: indx1
   real(d), dimension(dim1) :: vec

   Vec            ivec(blckdim)
   PetscInt       i
   PetscInt       indx(dim1)
   PetscReal      ftmp(dim1)
   PetscErrorCode ierr

!-----------------------------------------------------------------------
! Copy the index and element value arrays to their PETSc counterparts
!-----------------------------------------------------------------------
   indx=indx1
   ftmp=vec

!-----------------------------------------------------------------------
! Insert the element values into the ith PETSc starting vector
!
! Note that the only reason we (superfluously) call this subroutine is
! the inability to dynamicall allocate PETSc arrays, and that with the
! current implementation we do not know in advance the dimension of the
! ADC(1) eigenvectors in subroutine load_adc1_vecs
! 
! However, this could be easily remedied in the future by passing 
! kqp(1,0) to load_adc1_vecs
!-----------------------------------------------------------------------
   call vecsetvalues(ivec(i),dim1,indx,ftmp,INSERT_VALUES,ierr)   

!-----------------------------------------------------------------------
! Assemble the ith PETSc starting vector
!-----------------------------------------------------------------------
   call vecassemblybegin(ivec(i),ierr)
   call vecassemblyend(ivec(i),ierr)

   return

 end subroutine setvec

!#######################################################################

 subroutine guess_vecs_ondiag(blckdim,matdim,ivec,flag)

   use constants
   use misc, only: dsortindxa1

   implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/slepcsys.h"
#include "finclude/slepceps.h"
#include "finclude/petscvec.h90"

   integer                    :: blckdim,matdim,unit,maxbl,nrec,k
   integer, dimension(matdim) :: indx_hii
   real(d), dimension(matdim) :: hii
   character(1)               :: flag
   character(70)              :: filename

!-----------------------------------------------------------------------
! PETSc/SLEPc variables
!-----------------------------------------------------------------------
   Vec            ivec(blckdim)
   PetscInt       i,nvecs,dim
   PetscInt       indx
   PetscScalar    ftmp
   PetscErrorCode ierr

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
   unit=77

   if (flag.eq.'i') then
      filename='SCRATCH/hmlt.diai'
   else if (flag.eq.'f') then
      filename='SCRATCH/hmlt.diac'
   endif

   open(unit,file=filename,status='old',access='sequential',&
        form='unformatted')

!-----------------------------------------------------------------------
! Read the on-diagonal Hamiltonian matrix elements
!-----------------------------------------------------------------------
   rewind(unit)
   read(77) maxbl,nrec
   read(77) hii

!-----------------------------------------------------------------------
! Determine the indices of the on-diagonal elements with the smallest
! values
!-----------------------------------------------------------------------
   hii=abs(hii)   
   call dsortindxa1('A',matdim,hii,indx_hii)

   nvecs=blckdim
   dim=matdim

   do i=1,nvecs
      ! Create the ith initial vector
!      call veccreateseq(PETSC_COMM_SELF,dim,ivec(i),ierr)
      call veccreateseq(MPI_COMM_SELF,dim,ivec(i),ierr)
      ! Assign the components of the ith initial vector
      ftmp=1.0d0
      ! PETSc indices start from zero...
      indx=indx_hii(i)-1
      call vecsetvalues(ivec(i),1,indx,ftmp,INSERT_VALUES,ierr)
      call vecassemblybegin(ivec(i),ierr)
      call vecassemblyend(ivec(i),ierr)
   enddo

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
   close(unit)

   return

 end subroutine guess_vecs_ondiag

!#######################################################################
 
 subroutine hxvec(ham,x,y,ierr)

   implicit none

#include <finclude/petscsys.h>
#include <finclude/petscvec.h>
#include <finclude/petscvec.h90>
#include <finclude/petscmat.h>

   Mat                     ham
   Vec                     x,y
   PetscInt                n
   PetscScalar, pointer :: xvec(:),yvec(:)
   PetscErrorCode          ierr

   external hxvec_ext
   
!-----------------------------------------------------------------------
! Retrieve the vectors:
!
! xvec: input vector
! yvec: output vector
!-----------------------------------------------------------------------
   call MatGetSize(ham,n,PETSC_NULL_INTEGER,ierr)
   call VecGetArrayF90(x,xvec,ierr)
   call VecGetArrayF90(y,yvec,ierr)

!-----------------------------------------------------------------------
! Calculate the matrix-vector product y=Hx
!-----------------------------------------------------------------------
   call hxvec_ext(n,xvec,yvec)

!-----------------------------------------------------------------------
! Restore the Petsc vector arrays
!-----------------------------------------------------------------------
   call VecRestoreArrayF90(x,xvec,ierr)
   call VecRestoreArrayF90(y,yvec,ierr)
   
   return
   
 end subroutine hxvec

!#######################################################################

 subroutine hxvec_ext(dim,xvec,yvec)

   use constants
   use parameters, only: hamflag
   use iomod, only: freeunit
   
   implicit none

   PetscInt     dim
   PetscScalar  xvec(dim), yvec(dim)

   integer                            :: unit
   integer                            :: maxbl,nrec,nlim,k,l
   integer, dimension(:), allocatable :: indxi,indxj
   real(d), dimension(:), allocatable :: hii,hij
   character(70)                      :: filename
   
!-----------------------------------------------------------------------
! Contribution from the on-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
   allocate(hii(dim))
   
   call freeunit(unit)

   ! Note that the i/f flag that we pass to the main subroutine is
   ! going to have to be available globally if it is to be used here...

   if (hamflag.eq.'i') then
      filename='SCRATCH/hmlt.diai'
   else if (hamflag.eq.'f') then
      filename='SCRATCH/hmlt.diac'
   endif
      
   open(unit,file=filename,status='old',access='sequential',&
        form='unformatted')
   
   read(unit) maxbl,nrec
   read(unit) hii
   
   close(unit)

   yvec=0.0d0
   do k=1,dim
      yvec(k)=yvec(k)+hii(k)*xvec(k)
   enddo

   deallocate(hii)

!-----------------------------------------------------------------------
! Contribution from the off-diagonal elements of the Hamiltonian matrix
!-----------------------------------------------------------------------
   allocate(hij(maxbl),indxi(maxbl),indxj(maxbl))

   if (hamflag.eq.'i') then
      filename='SCRATCH/hmlt.offi'
   else if (hamflag.eq.'f') then
      filename='SCRATCH/hmlt.offc'
   endif
   
   open(unit,file=filename,status='old',access='sequential',&
        form='unformatted')

   do k=1,nrec
      read(unit) hij(:),indxi(:),indxj(:),nlim
      do l=1,nlim
         yvec(indxi(l))=yvec(indxi(l))+hij(l)*xvec(indxj(l))
         yvec(indxj(l))=yvec(indxj(l))+hij(l)*xvec(indxi(l))
      enddo
   enddo

   close(unit)

   deallocate(hij,indxi,indxj)
   
   return
   
 end subroutine hxvec_ext

!#######################################################################
   
 function isincore(matdim,noffd)

   use constants

   implicit none

   integer   :: matdim
   integer*8 :: noffd
   real(d)   :: hammem
   logical   :: isincore
   
   ! Non-zero matrix elements
   hammem=8.0d0*(matdim+noffd)
   
   ! Indices (2x 32 bit per off-diagonal element)
   hammem=hammem+8.0d0*noffd
   
   ! Convert to mb
   hammem=hammem/1024.0d0

   ! Temporary fix: if the Hamiltonian will occupy more than 10 gb of
   ! memory, then do not store it in-core
   ! To do: make this decision based on how much memory is actually
   ! available (which is tricky as we do not know how many vectors
   ! are being stored by PETSc/SLEPc...)
   if (hammem.lt.10000.0d0) then
      isincore=.true.
   else
      isincore=.false.
   endif

   return

 end function isincore

!#######################################################################


