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

   integer                             :: blckdim,matdim,davstates,&
                                          maxbl,nrec,unit
   real(d), dimension(matdim)          :: hii
   double precision                    :: val
   double precision, dimension(matdim) :: vec
   character(36)                       :: davname

 end subroutine davidson_diag

!#######################################################################
! fill_ondiag: Assembles the on-diagonal elements of the PETSc version
!              of the ADC Hamiltonian matrix
!#######################################################################

 subroutine fill_ondiag(maxbl,nrec,matdim)

   use constants

   implicit none

   integer                    :: matdim,maxbl,nrec
   real(d), dimension(matdim) :: hii


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
   open(unit,file='hmlt.dia',status='old',access='sequential',&
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

 subroutine fill_offdiag(maxbl,nrec)

   use constants

   implicit none

   integer                   :: maxbl,nrec,nlim,k,l
   integer*8                 :: unit
   real(d), dimension(maxbl) :: hij
   integer, dimension(maxbl) :: indxi,indxj


   return

 end subroutine fill_offdiag

!#######################################################################
