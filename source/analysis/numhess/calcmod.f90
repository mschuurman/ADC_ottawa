  module calcmod

    use constants

    implicit none

    save
    
    integer                                 :: nzero,nmodes
    integer, dimension(:), allocatable      :: iimag
    real(dp), dimension(:), allocatable     :: ener,mass,freq0
    real(dp), dimension(:,:), allocatable   :: grad,gradq0,q0,hess_gs
    real(dp), dimension(:,:,:), allocatable :: hess,hessq0
    real(dp), parameter                     :: eh2ev=27.2113845d0
    real(dp), parameter                     :: b2a=0.529177249d0

  contains

!#######################################################################

    subroutine hesscalc

      use constants
      use iomod
      use hessmod
      
      implicit none

!-----------------------------------------------------------------------
! Initialise things
!-----------------------------------------------------------------------
      call initialise

!-----------------------------------------------------------------------
! Read the names of the ADC output files
!-----------------------------------------------------------------------
      call rdoutfilenames

!-----------------------------------------------------------------------
! Determine the numbers of initial and final space states
!-----------------------------------------------------------------------
      call getnsta

!-----------------------------------------------------------------------
! Read the ADC state energies at the displaced geometries
!-----------------------------------------------------------------------
      call getener

!-----------------------------------------------------------------------
! Calculate the gradients and Hessians
!-----------------------------------------------------------------------
      call calc_hess

!-----------------------------------------------------------------------
! Mass-weight the ground-state Hessian
!-----------------------------------------------------------------------
      call mass_weight

!-----------------------------------------------------------------------
! Project out the rotational and translational degrees of freedom
! from the ground state Hessians
!-----------------------------------------------------------------------
      call proj_hess

!-----------------------------------------------------------------------
! Diagonalise the ground state projected Hessian
!-----------------------------------------------------------------------
      call diaghess_gs

!-----------------------------------------------------------------------
! Transform the gradients and Hessians to be in terms of the ground 
! state mass- and frequency-scaled normal modes
!-----------------------------------------------------------------------
      call transform

!-----------------------------------------------------------------------
! Ouput the transformed gradients and Hessians
!-----------------------------------------------------------------------
      call wrout

      return

    end subroutine hesscalc

!#######################################################################
    
    subroutine initialise

      use hessmod

      implicit none

      ! Number of geometries
      ngeom=1+2*ncoo+4*(ncoo*(ncoo-1)/2)

      ! Names of the ADC output files
      allocate(adcout(ngeom))
      adcout=''

      return

    end subroutine initialise

!#######################################################################

    subroutine rdoutfilenames

      use parsemod
      use iomod
      use hessmod

      implicit none

      integer :: unit,k

!-----------------------------------------------------------------------
! Open the list file
!-----------------------------------------------------------------------
      call freeunit(unit)
      open(unit,file=listfile,form='formatted',status='old')

!-----------------------------------------------------------------------
! Read the names of the ADC output files
!-----------------------------------------------------------------------
      k=0
5     call rdinp(unit)      
      if (.not.lend) then
         k=k+1
         if (k.gt.ngeom) goto 999
         adcout(k)=trim(keyword(1))         
         goto 5
      endif

      if (k.lt.ngeom) goto 888

!-----------------------------------------------------------------------
! Close the list file
!-----------------------------------------------------------------------
      close(unit)

      return

888   continue
      errmsg='There are more ADC output files than displacements'
      call error_control

999   continue
      errmsg='There are more ADC output files than displacements'
      call error_control

    end subroutine rdoutfilenames

!#######################################################################

    subroutine getnsta

      use iomod
      use hessmod

      implicit none

      integer            :: unit
      character(len=120) :: string
      logical            :: lfspace

!-----------------------------------------------------------------------
! Open file
!-----------------------------------------------------------------------
      call freeunit(unit)
      open(unit,file=adcout(1),form='formatted',status='old')

!-----------------------------------------------------------------------
! Number of initial space states (including the ground state)
!-----------------------------------------------------------------------
      nsta_i=1

5     read(unit,'(a)') string
      if (string(3:15).ne.'Initial space') goto 5
      read(unit,*)

10    read(unit,'(a)') string
      if (string(1:1).ne.'*') then
         if (string(3:7).eq.'State') nsta_i=nsta_i+1
         goto 10
      endif

!-----------------------------------------------------------------------
! Number of final space states
!-----------------------------------------------------------------------
      nsta_f=0

15    read(unit,'(a)',end=25) string
      if (string(3:13).ne.'Final space') goto 15
      read(unit,*)

20    read(unit,'(a)') string
      if (string(1:15).ne.'Final wall time') then
         if (string(3:7).eq.'State') nsta_f=nsta_f+1
         goto 20
      endif

25    continue

!-----------------------------------------------------------------------
! Total number of states (including the ground state)
!-----------------------------------------------------------------------
      nsta=nsta_i+nsta_f

!-----------------------------------------------------------------------
! Close file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine getnsta

!#######################################################################

    subroutine getener

      use iomod
      use hessmod

      implicit none

      integer               :: i,itype,unit
      integer, dimension(2) :: xindx

!-----------------------------------------------------------------------
! Allocate and initialise arrays
!-----------------------------------------------------------------------
      allocate(ener(nsta))
      allocate(ref(nsta))
      allocate(pos(nsta,ncoo))
      allocate(neg(nsta,ncoo))
      allocate(pospos(nsta,ncoo,ncoo))
      allocate(negneg(nsta,ncoo,ncoo))
      allocate(posneg(nsta,ncoo,ncoo))
      allocate(negpos(nsta,ncoo,ncoo))
      ener=0.0d0
      ref=0.0d0
      pos=0.0d0
      neg=0.0d0
      pospos=0.0d0
      negneg=0.0d0
      posneg=0.0d0
      negpos=0.0d0

!-----------------------------------------------------------------------
! Fill in the various energy arrays
!-----------------------------------------------------------------------
      call freeunit(unit)

      do i=1,ngeom
         
         ! Open the current ADC output file
         open(unit,file=adcout(i),form='formatted',status='old')

         ! Read the Cartesian coordinates
         call rdgeom(unit)

         ! Determine which type of displacement we are dealing with
         call disptype(itype,xindx)

         ! Read the state energies
         call rdadcener(unit)

         ! Fill in the correct energy array for the current
         ! displacement
         call fillarr(itype,xindx)

         ! Close the current ADC output file
         close(unit)

      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(ener)

      return

    end subroutine getener

!#######################################################################

    subroutine rdgeom(unit)

      use hessmod

      implicit none

      integer            :: unit,i,j
      character(len=120) :: string

      rewind(unit)

5     read(unit,'(a)') string
      if (index(string,'Atom ').eq.0) goto 5
      read(unit,*)
      
      do i=1,natm
         read(unit,'(13x,3(6x,F14.10))') (xcoo(j),j=i*3-2,i*3)
      enddo
      
      return

    end subroutine rdgeom

!#######################################################################

    subroutine disptype(itype,xindx)

      use hessmod

      implicit none

      integer                :: itype,i,k,nneg,npos,nindx,pindx
      integer, dimension(2)  :: xindx,ilbl
      real(dp)               :: diff
      real(dp), parameter    :: tol=1e-6_dp

!-----------------------------------------------------------------------
! Set the displacement type:
!
! itype = 0 <-> reference geometry
!         1 <-> pos.
!         2 <-> neg.
!         3 <-> pos. pos.
!         4 <-> neg. neg.
!         5 <-> pos. neg.
!         6 <-> neg. pos.
!-----------------------------------------------------------------------
      nneg=0
      npos=0

      ilbl=0
      k=0
      do i=1,ncoo
         diff=xcoo(i)-xcoo0(i)
         if (abs(diff).gt.tol) then
            k=k+1
            if (diff.lt.0.0d0) then
               nneg=nneg+1
               ilbl(k)=-i
            else if (diff.gt.0.0d0) then
               npos=npos+1
               ilbl(k)=i
            endif            
         endif
      enddo
      
      if (npos.eq.0.and.nneg.eq.0) then

         ! Reference point
         itype=0     

      else if (nneg.eq.0.and.npos.eq.1) then
         
         ! Pos.
         itype=1

      else if (nneg.eq.1.and.npos.eq.0) then

         ! Neg.
         itype=2

      else if (npos.eq.2) then

         ! Pos. pos
         itype=3

      else if (nneg.eq.2) then
         
         ! Neg. neg.
         itype=4

      else if (nneg.eq.1.and.npos.eq.1) then

         if (ilbl(1).lt.0) then
            nindx=abs(ilbl(1))
            pindx=ilbl(2)
         else
            nindx=abs(ilbl(2))
            pindx=ilbl(1)
         endif

         if (nindx.gt.pindx) then
            
            ! Pos. neg.
            itype=5

         else if (nindx.lt.pindx) then
            
            ! Neg. pos.
            itype=6

         endif
         
      endif

!-----------------------------------------------------------------------
! Fill in the coordinate index array
!-----------------------------------------------------------------------
      xindx=0
      ilbl=abs(ilbl)

      if (ilbl(2).eq.0) then
         xindx(1)=ilbl(1)
      else
         if (ilbl(1).lt.ilbl(2)) then
            xindx(1)=ilbl(1)
            xindx(2)=ilbl(2)
         else
            xindx(1)=ilbl(2)
            xindx(2)=ilbl(1)
         endif
      endif

      return

    end subroutine disptype

!#######################################################################

    subroutine rdadcener(unit)

      use hessmod

      implicit none

      integer            :: unit,k
      character(len=120) :: string

!-----------------------------------------------------------------------
! Ground state MP2 energy
!-----------------------------------------------------------------------
      rewind(unit)
      
5     read(unit,'(a)') string
      if (string(5:14).ne.'MP2 energy') goto 5

      read(string,'(31x,F16.10)') ener(1)
      
!-----------------------------------------------------------------------
! Initial space excited state energies
!-----------------------------------------------------------------------
10    read(unit,'(a)') string
      if (string(3:15).ne.'Initial space') goto 10
      read(unit,*)
      
      k=1
15    read(unit,'(a)') string
      if (string(1:1).ne.'*') then
         if (string(3:7).eq.'State') then
            k=k+1
            read(string,'(27x,F14.8)') ener(k)
         endif
         goto 15
      endif
      
!-----------------------------------------------------------------------
! Final space excited state energies
!-----------------------------------------------------------------------
      if (nsta_f.gt.0) then

20       read(unit,'(a)') string
         if (string(3:13).ne.'Final space') goto 20
         read(unit,*)

25       read(unit,'(a)') string
         if (string(1:15).ne.'Final wall time') then
            if (string(3:7).eq.'State') then
               k=k+1
               read(string,'(27x,F14.8)') ener(k)
            endif
            goto 25
         endif

      endif

      return

    end subroutine rdadcener

!#######################################################################

    subroutine fillarr(itype,xindx)

      use hessmod

      implicit none
      
      integer               :: itype,i
      integer, dimension(2) :: xindx

! itype = 0 <-> reference geometry
!         1 <-> pos.
!         2 <-> neg.
!         3 <-> pos. pos.
!         4 <-> neg. neg.
!         5 <-> pos. neg.
!         6 <-> neg. pos.

      if (itype.eq.0) then
         ref=ener
      else if (itype.eq.1) then
         pos(:,xindx(1))=ener
      else if (itype.eq.2) then
         neg(:,xindx(1))=ener
      else if (itype.eq.3) then
         pospos(:,xindx(1),xindx(2))=ener
      else if (itype.eq.4) then
         negneg(:,xindx(1),xindx(2))=ener
      else if (itype.eq.5) then
         posneg(:,xindx(1),xindx(2))=ener
      else if (itype.eq.6) then
         negpos(:,xindx(1),xindx(2))=ener
      endif

      return

    end subroutine fillarr

!#######################################################################

    subroutine calc_hess

      use constants
      use hessmod

      implicit none

      integer  :: n,i,j
      real(dp) :: diff

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(hess(nsta,ncoo,ncoo))
      hess=0.0d0

      allocate(grad(nsta,ncoo))
      grad=0.0d0
      
!-----------------------------------------------------------------------
! Calculate the Hessian for each state
!-----------------------------------------------------------------------
      diff=dx/0.529177249d0

      do n=1,nsta

         ! (1) Quadratic (on-diagonal) terms
         do i=1,ncoo
            hess(n,i,i)=(pos(n,i)+neg(n,i)-2*ref(n))/(diff**2)
         enddo

         ! (2) Bi-linear (off-diagonal) terms
         do i=1,ncoo-1
            do j=i+1,ncoo
               hess(n,i,j)=(pospos(n,i,j)+negneg(n,i,j)-posneg(n,i,j) &
                    -negpos(n,i,j))/(4.0d0*diff**2)
               hess(n,j,i)=hess(n,i,j)
            enddo
         enddo

      enddo

!-----------------------------------------------------------------------
! Save a copy of the ground-state Hessian to be used in the calculation
! of the ground state normal modes
!-----------------------------------------------------------------------
      allocate(hess_gs(ncoo,ncoo))
      hess_gs=hess(1,:,:)

!-----------------------------------------------------------------------
! Calculate the gradients for each state
!-----------------------------------------------------------------------
      do n=1,nsta
         do i=1,ncoo
            grad(n,i)=(pos(n,i)-neg(n,i))/(2.0d0*diff)
         enddo
      enddo

!-----------------------------------------------------------------------
! Convert the gradients to units of eV/Angtrom 
!-----------------------------------------------------------------------
      grad=grad*eh2ev/b2a

!-----------------------------------------------------------------------
! Convert the Hessians to units of eV/Angtrom^2
!-----------------------------------------------------------------------
      hess=hess*eh2ev/b2a**2

      return

    end subroutine calc_hess

!#######################################################################

    subroutine mass_weight

      use constants
      use iomod
      use hessmod

      implicit none

      integer :: i,j

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(mass(ncoo))

!-----------------------------------------------------------------------
! Set the atomic masses
!-----------------------------------------------------------------------
      do i=1,natm
         if (aatm(i).eq.'C'.or.aatm(i).eq.'c') then
            mass(i*3-2:i*3)=12.0107d0
         else if (aatm(i).eq.'H'.or.aatm(i).eq.'h') then
            mass(i*3-2:i*3)=1.00782504d0
         else if (aatm(i).eq.'N'.or.aatm(i).eq.'n') then
            mass(i*3-2:i*3)=14.0067d0
         else
            errmsg='Unknown atom type: '//trim(aatm(i))
            call error_control
         endif
      enddo

!-----------------------------------------------------------------------
! Mass-weight the ground state Hessian
!-----------------------------------------------------------------------
      do i=1,ncoo
         do j=1,ncoo
            hess_gs(i,j)=hess_gs(i,j)/sqrt(mass(i)*mass(j))
         enddo
      enddo
      
      return

    end subroutine mass_weight

!#######################################################################
! The projector onto the translation-rotation subspace is constructed
! by:
!
! (i)  Constructing six (non-orthogonal) vectors spanning the 
!      translation-rotation subspace;
! (ii) Constructing the contribution to the projector by
!      orthogonalising these vectors
!
! The procedure used is taken from J. Chem. Phys., 88, 922 (1988)
!#######################################################################

    subroutine proj_hess

      use constants
      use iomod
      use hessmod

      implicit none

      integer                        :: i,j,k,l,info
      integer, dimension(6)          :: ipiv
      real(dp), dimension(6,ncoo)    :: vec
      real(dp), dimension(6,6)       :: smat,invsmat
      real(dp), dimension(ncoo,6)    :: bmat
      real(dp), dimension(ncoo,ncoo) :: rmat,pmat,tmpmat
      real(dp), dimension(6)         :: work

!------------------------------------------------------------------
! Initialise arrays
!------------------------------------------------------------------
      vec=0.0d0

!------------------------------------------------------------------
! Vectors 1-3: translation along the three Cartesian axes
!------------------------------------------------------------------
      ! Loop over the translational DOFs
      do i=1,3
         ! Construct the vector for the current DOF
         do j=1,natm
            k=j*3-3+i
            vec(i,k)=sqrt(mass(j*3))
         enddo
      enddo

!------------------------------------------------------------------
! Vectors 4-6: infinitesimal displacements corresponding to
!              rotation about the three Cartesian axes
!------------------------------------------------------------------
      ! Rotation about the x-axis
      do i=1,natm
         j=i*3-1
         k=i*3
         vec(4,j)=sqrt(mass(i*3))*xcoo0(k)
         vec(4,k)=-sqrt(mass(i*3))*xcoo0(j)
      enddo

      ! Rotation about the y-axis
      do i=1,natm
         j=i*3-2
         k=i*3
         vec(5,j)=-sqrt(mass(i*3))*xcoo0(k)
         vec(5,k)=sqrt(mass(i*3))*xcoo0(j)
      enddo

      ! Rotation about the z-axis
      do i=1,natm
         j=i*3-2
         k=i*3-1
         vec(6,j)=sqrt(mass(i*3))*xcoo0(k)
         vec(6,k)=-sqrt(mass(i*3))*xcoo0(j)
      enddo

!------------------------------------------------------------------
! Calculate the projector R onto the translational and rotational
! DOFs using R=b*S^-1*b^T, where S=vec^T*vec.
!
! Here, R <-> rmat, S <-> smat, b <-> bmat (matrix of vectors)
!------------------------------------------------------------------
      ! Construct the b-matrix
      bmat=0.0d0
      do i=1,6
         do j=1,ncoo
            bmat(j,i)=vec(i,j)
         enddo
      enddo

      ! Calculate the S-matrix
      smat=0.0d0
      do i=1,6
         do j=1,6
            do k=1,ncoo
               smat(i,j)=smat(i,j)+bmat(k,i)*bmat(k,j)
            enddo
         enddo
      enddo

      ! Invert the S-matrix
      invsmat=smat
      call dgetrf(6,6,invsmat,6,ipiv,info)
      if (info.ne.0) then
         errmsg='LU factorisation of the S-matrix failed'
         call error_control
      endif
      call dgetri(6,invsmat,6,ipiv,work,6,info)
      if (info.ne.0) then
         errmsg='Diagonalisation of the S-matrix failed'
         call error_control
      endif

      ! Calculate the projection matrix R <-> rmat
      rmat=0.0d0
      do i=1,ncoo
         do j=1,ncoo
            do k=1,6
               do l=1,6
                  rmat(i,j)=rmat(i,j)+&
                       bmat(i,k)*invsmat(k,l)*bmat(j,l)
               enddo
            enddo
         enddo
      enddo

!------------------------------------------------------------------
! Project the Hessians
!------------------------------------------------------------------
      ! Construct the projector P=1-R
      pmat=0.0d0
      do i=1,ncoo
         pmat(i,i)=1.0d0
      enddo
      pmat=pmat-rmat

      ! Similarity transform the ground state Hessian using the 
      ! projection matrix P
      do i=1,nsta
         tmpmat=matmul(hess_gs(:,:),pmat)
         hess_gs(:,:)=matmul(pmat,tmpmat)
      enddo

      return

    end subroutine proj_hess

!#######################################################################

    subroutine diaghess_gs

      use constants
      use iomod
      use hessmod
      
      implicit none
      
      integer                        :: i,j,k,error,e2
      integer, dimension(ncoo)       :: indx
      real(dp), dimension(ncoo,ncoo) :: tmp
      real(dp), dimension(ncoo)      :: lambda
      real(dp), dimension(3*ncoo)    :: work
      real(dp), parameter            :: tol=1e-5_dp
      real(dp)                       :: norm,ftmp
      
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(iimag(ncoo))
      allocate(freq0(ncoo))
      allocate(q0(ncoo,ncoo))

!-----------------------------------------------------------------------
! Diagonalise the ground state projected, mass-weigthe Hessian
!-----------------------------------------------------------------------
      q0=hess_gs(:,:)

      e2=3*ncoo

      call dsyev('V','U',ncoo,q0,ncoo,freq0,work,e2,error)

      if (error.ne.0) then
         errmsg='Diagonalisation of the projected ground state &
              Hessian failed'
         call error_control
      endif

!-----------------------------------------------------------------------
! Flag any imaginary frequencies and convert the eigenvalues to
! frequencies
!-----------------------------------------------------------------------
      iimag=0
      do i=1,ncoo
         if (abs(freq0(i)).gt.tol.and.freq0(i).lt.0) iimag(i)=1
         freq0(i)=sqrt(abs(freq0(i)))
      enddo

!-----------------------------------------------------------------------
! Check the number of zero frequencies
!-----------------------------------------------------------------------
      nzero=0
      do i=1,ncoo
         if (freq0(i).lt.tol) nzero=nzero+1
      enddo

      if (nzero.ne.5.and.nzero.ne.6) then
         errmsg=''
         write(errmsg,'(a,1x,i2)') &
         'Unexpected number of zero frequencies:',nzero
         call error_control
      endif

      nmodes=ncoo-nzero

!-----------------------------------------------------------------------
! Reorder the ground state normal mode and frequency arrays s.t. the
! zero-frequency terms come last
!-----------------------------------------------------------------------
      call dsortindxa1('A',ncoo,freq0,indx)

      k=0
      do i=nzero+1,ncoo
         k=k+1
         tmp(:,k)=q0(:,indx(i))
         lambda(k)=freq0(indx(i))
      enddo
      
      k=ncoo-nzero
      do i=1,nzero
         k=k+1
         tmp(:,k)=q0(:,indx(i))
         lambda(k)=freq0(indx(i))
      enddo

      q0=tmp
      freq0=lambda

!-----------------------------------------------------------------------
! Convert the frequencies to eV
!-----------------------------------------------------------------------
      freq0=freq0*0.6373641d0

!-----------------------------------------------------------------------
! Mass and frequency scale the transformation matrix.
!
! N.B. the scaling used here is for frequencies in eV, mass in amu
!      and length in Angstrom
!-----------------------------------------------------------------------
      do i=1,ncoo-nzero
         do j=1,ncoo
            q0(j,i)=q0(j,i)/(15.4644*sqrt(freq0(i))*sqrt(mass(j)))
         enddo
      enddo

      return
      
    end subroutine diaghess_gs

!#######################################################################

    subroutine dsortindxa1(order,ndim,arrin,indx)

      use constants

      implicit none

      character(1), intent(in) :: order
      integer, intent(in) :: ndim
      real(dp), dimension(ndim), intent(in)   :: arrin
      integer, dimension(ndim), intent(inout) :: indx
    
      integer  :: i,l,ir,indxt,j
      real(dp) :: q

!!$ The subroutine is taken from the NR p233, employs heapsort.

      do i= 1,ndim
         indx(i)=i
      end do
      
      l=ndim/2+1
      ir=ndim
      
      if(order .eq. 'D') then
         
10       continue
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
         
20       if(j .le. ir) then
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
       
100      continue
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
      
200   if(j .le. ir) then
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
    
   return
   
 end subroutine dsortindxa1

!#######################################################################

    subroutine transform

      use constants
      use hessmod

      implicit none

      integer                        :: i,j
      real(dp), dimension(ncoo,ncoo) :: tmpmat

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(hessq0(nsta,ncoo,ncoo))
      hessq0=0.0d0

      allocate(gradq0(nsta,ncoo))
      gradq0=0.0d0

!-----------------------------------------------------------------------
! Hessians
!-----------------------------------------------------------------------
      do i=1,nsta
         tmpmat(:,1:nmodes)=matmul(hess(i,:,:),q0(:,1:nmodes))
         hessq0(i,1:nmodes,1:nmodes)=matmul(transpose(q0(:,1:nmodes)),tmpmat(:,1:nmodes))
      enddo

!-----------------------------------------------------------------------
! Gradients
!-----------------------------------------------------------------------
      do i=1,nsta
         gradq0(i,1:nmodes)=matmul(transpose(q0(:,1:nmodes)),grad(i,:))
      enddo

      return

    end subroutine transform

!#######################################################################

    subroutine wrout
      
      use constants
      use iomod
      use hessmod

      implicit none

      integer :: unit

!-----------------------------------------------------------------------
! Ground state normal modes and frequencies
!-----------------------------------------------------------------------
      call wrmodes

!-----------------------------------------------------------------------
! Gradients and Hessians in terms of the ground state normal modes
!-----------------------------------------------------------------------
      call freeunit(unit)
      open(unit,file='numhess.out',form='formatted',status='unknown')

      call wrnumhessout(unit)

      close(unit)

      return

    end subroutine wrout

!#######################################################################
    
    subroutine wrmodes
      
      use constants
      use iomod
      use hessmod

      implicit none

      integer                        :: unit,i,j,k
      real(dp)                       :: ftmp
      real(dp), dimension(ncoo,ncoo) :: qu
      character(len=60)              :: atmp

!-----------------------------------------------------------------------
! Open the xyz file
!-----------------------------------------------------------------------
      call freeunit(unit)
      open(unit,file='eigvec.xyz',form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Undo the mass-scaling of the normal modes
!-----------------------------------------------------------------------
      do i=1,ncoo
         do j=1,ncoo
            qu(j,i)=q0(j,i)/sqrt(mass(j))
         enddo
         qu(:,i)=qu(:,i)/sqrt(dot_product(qu(:,i),qu(:,i)))
      enddo

!-----------------------------------------------------------------------
! Write the unscaled ground state normal modes and frequencies to file
!-----------------------------------------------------------------------
      do i=1,ncoo-nzero
         
         ftmp=freq0(i)*5140.66d0

         atmp=''
         if (i.lt.10) then
            write(atmp,'(a1,i1,2x,F10.4)') 'Q',i,ftmp
         else if (i.lt.100) then
            write(atmp,'(a1,i2,2x,F10.4)') 'Q',i,ftmp
         else
            write(atmp,'(a1,i3,2x,F10.4)') 'Q',i,ftmp
         endif
         
         if (iimag(i).eq.1) then
            atmp=trim(atmp)//'I cm-1'
         else
            atmp=trim(atmp)//'  cm-1'
         endif

         write(unit,'(i3)') natm
         write(unit,'(a)') trim(atmp)
         
         do j=1,natm
            write(unit,'(a1,6(2x,F10.7))') aatm(j),&
                 (xcoo0(k),k=j*3-2,j*3),(qu(k,i),k=j*3-2,j*3)
         enddo

      enddo

!-----------------------------------------------------------------------
! Close the xyz file
!-----------------------------------------------------------------------
      close(unit)

      return

    end subroutine wrmodes

!#######################################################################

    subroutine wrnumhessout(unit)

      use constants
      use iomod
      use hessmod

      implicit none

      integer :: unit,i,j,n

!-----------------------------------------------------------------------
! Preamble
!-----------------------------------------------------------------------
      write(unit,'(68a)') ('#',i=1,68)
      write(unit,'(a)') '# Gradients and Hessians in terms of the &
           dimensionless mass- and'
      write(unit,'(a)') '# frequency-scaled ground state normal modes.'
      write(unit,'(a)') '# All quantities are given in units of eV.'
      write(unit,'(68a)') ('#',i=1,68)

!-----------------------------------------------------------------------
! System dimensions
!-----------------------------------------------------------------------
      write(unit,'(2/,a,3x,i3)') '# nstates:',nsta
      write(unit,'(a,4x,i3)') '# nmodes:',nmodes

!-----------------------------------------------------------------------
! Frequencies
!-----------------------------------------------------------------------
      write(unit,'(2/,a)') '# Frequencies'
      do i=1,nmodes
         write(unit,'(i3,3x,F6.4)') i,freq0(i)
      enddo

!-----------------------------------------------------------------------
! Energies
!-----------------------------------------------------------------------
      write(unit,'(2/,a)') '# Energies'
      do i=1,nsta
         write(unit,'(i3,3x,F8.4)') i,(ref(i)-ref(1))*eh2ev
      enddo

!-----------------------------------------------------------------------
! Gradients
!-----------------------------------------------------------------------
      write(unit,'(2/,a)') '# Gradients'
      do n=1,nsta
         write(unit,'(/,a,1x,i3)') 'State:',n
         do i=1,nmodes
            write(unit,'(i3,3x,F7.4)') i,gradq0(n,i)
         enddo
      enddo

!-----------------------------------------------------------------------
! Hessians
!-----------------------------------------------------------------------
      write(unit,'(2/,a)') '# Hessians: On-diagonal elements'
      do n=1,nsta
         write(unit,'(/,a,1x,i3)') 'State:',n
         do i=1,nmodes
            write(unit,'(i3,3x,F7.4)') i,hessq0(n,i,i)
         enddo
      enddo

      write(unit,'(2/,a)') '# Hessians: Off-diagonal elements'
      do n=1,nsta
         write(unit,'(/,a,1x,i3)') 'State:',n
         do i=1,nmodes
            do j=i+1,nmodes
               write(unit,'(2(i3,3x),F7.4)') i,j,hessq0(n,i,j)
            enddo
         enddo
      enddo

      return

    end subroutine wrnumhessout

!#######################################################################

  end module calcmod
