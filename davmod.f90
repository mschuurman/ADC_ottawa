  module davmod
  
    use constants
    use parameters
  
    implicit none 

    integer                            :: ndm,ndms,mem,main,nrec
    integer*8                          :: noffdel
 
    real(d), dimension(:), allocatable :: offdiag,diag
    integer, dimension(:), allocatable :: indi,indip
    integer, dimension(:), allocatable :: stbl_index
    character(1)                       :: fl
  
  contains

!#######################################################################

    subroutine master_dav(ndim,noffd,flag,ndims)
    
      integer, intent(in)      :: ndim,ndims
      integer*8, intent(in)    :: noffd
      character(1), intent(in) :: flag

      integer                  :: maxbl, nrec,i
      logical                  :: rstr,prs

!-----------------------------------------------------------------------
! Unknown flag: this appears to never be used here...
!-----------------------------------------------------------------------
      fl=flag
!!$    fl='i'

!-----------------------------------------------------------------------
! ndm:     Hamiltonian matrix dimension
! ndms:    No. 1h1p ISs
! noffdel: No. non-zero off-diagonal Hamiltonian matrix elements ???
!-----------------------------------------------------------------------
      noffdel=noffd
      ndm=ndim

!-----------------------------------------------------------------------
! Check to see if the dav_vecs file holding the Davidson eigenvectors
! already exists:
!
!    - If yes, do nothing
!
!    - If no, read the Hamiltonian matrix from file and perform Davidson
!      diagonalisation
!-----------------------------------------------------------------------
      inquire(file='dav_vecs.'//mtxidd,exist=prs)

      if (prs) then
         write(6,*) 'Older file ','dav_vecs.'//mtxidd,' will be engaged'
      else

!-----------------------------------------------------------------------
! Output restart status
!-----------------------------------------------------------------------
         if (rstr) then
            write(6,*) 'This is the Davidson restart run'
         else
            write(6,*) 'This is the Davidson startup run'
         endif
!-----------------------------------------------------------------------
! Enter block-Davidson routine
!-----------------------------------------------------------------------
         call davidson_diag(dmain,ndm,davstates,davname,ladc1guess,ndms)


      end if

    end subroutine master_dav

!#######################################################################

    subroutine read_matrix()
    
      real(d)                            :: mtrl
      integer                            :: mark1,maxbl,nrec
      integer                            :: irec,nlim,i
      integer*8                          :: count
    
      real(d), dimension(:), allocatable :: dgl,offdgl
      integer, dimension(:), allocatable :: oi,oj
    

!!$ Reading diagonal elements
    
    
      OPEN(UNIT=77,FILE='hmlt.dia'//fl,STATUS='OLD',ACCESS='SEQUENTIAL',&
           FORM='UNFORMATTED')
    
      read(77) maxbl,nrec
      read(77) diag(:)
    
      CLOSE(77)

      allocate(offdgl(maxbl),oi(maxbl),oj(maxbl))
      
      OPEN(UNIT=78,FILE='hmlt.off'//fl,STATUS='OLD',ACCESS='SEQUENTIAL',&
           FORM='UNFORMATTED')
    
      count=0
      do irec= 1,nrec
         read(78) offdgl(:),oi(:),oj(:),nlim
         offdiag(count+1:count+nlim)=offdgl(1:nlim)
         indi(count+1:count+nlim)=oi(1:nlim)
         indip(count+1:count+nlim)=oj(1:nlim)
         count=count+nlim
      end do
      close(78)
      deallocate(offdgl,oi,oj)
      
    end subroutine read_matrix

!#######################################################################

    subroutine read_matrix2()

      implicit none

      integer :: maxbl,rc

!!$ Reading diagonal elements

      OPEN(UNIT=77,FILE='hmlt.dia'//fl,STATUS='OLD',ACCESS='SEQUENTIAL',&
           FORM='UNFORMATTED')
      
      read(77) maxbl,nrec
      read(77) diag
      
      CLOSE(77)

      write(*,*) 'read: nrec',nrec
    
      rc=0
      allocate(offdiag(maxbl), indi(maxbl), indip(maxbl), stat=rc)
      if (rc.ne.0) call errmsg('memory allocation error read_222')

      return

    end subroutine read_matrix2

!#######################################################################

    subroutine rmtxhd(dgl,work2)

      real(d), dimension(ndm), intent(out)   :: dgl
      real(d), dimension(mem), intent(inout) :: work2

      OPEN(UNIT=12,FILE='hmlt.dia'//fl,STATUS='OLD',ACCESS='SEQUENTIAL',&
           FORM='UNFORMATTED')

      read(12) 
      read(12) dgl(:)
      
      CLOSE(12)
    
    end subroutine rmtxhd

!#######################################################################

    subroutine rmtxq1(b,i1,i2,work1)

      integer, intent(in)                         :: i1,i2
      real(d), dimension(ndm,i2-i1+1),intent(out) :: b
      real(d), dimension(mem),intent(inout)       :: work1
      integer                                     :: i


      b(:,:)=rzero
      do i=i1,i2
         b(i,i-i1+1)=rone
      end do
      

    end subroutine rmtxq1

!#######################################################################

    subroutine mtxq1_l(amx,bmx,i1,i2,work)

      implicit none

      real(d), dimension(ndm,i2-i1+1) :: bmx
      real(d) :: work,amx
      integer :: i1,i2,i

      bmx=0.0d0
   
      write(*,*) '---mtx_l----------'
      write(*,*) 'lmain',lmain
      write(*,*) 'main',main
      write(*,*) 'stvc',stvc_lbl(1)
      write(*,*) 'stvc',stvc_lbl(main)
      write(*,*) 'stvc',stvc_lbl(lmain)
      write(*,*) '---mtx_l----------'
      
      do i=i1,i2
         bmx(stvc_lbl(i),i-i1+1)=1.0d0
      end do
      
    end subroutine mtxq1_l

!#######################################################################

    subroutine  rmtxhq1(qv,zv,nq,work3)

      implicit none

      integer, intent(in)        :: nq
      real(d), dimension(nq,ndm) :: qv
      real(d), dimension(nq,ndm) :: zv
      real(d)                    :: work3
      
      real(d)                    :: mtrl
      integer                    :: mark1,ndim,maxbl,nrec,type
      integer*8                  :: irec,nlim,i,j,irow,jcol

      zv=0.0d0

      do i= 1,ndm
         zv(:,i)=diag(i)*qv(:,i)
      end do
      
      do i= 1,noffdel
         irow=indi(i)
         jcol=indip(i)
         mtrl=offdiag(i)
         
         zv(:,irow)=zv(:,irow)+mtrl*qv(:,jcol)
         zv(:,jcol)=zv(:,jcol)+mtrl*qv(:,irow)         
      end do
 
    end subroutine rmtxhq1

!#######################################################################

    subroutine mtxhq_gg(q,z,nq,wrk)

      implicit none

      integer :: nq, i, nw, k
      real(8) :: q(nq,ndm), z(nq,ndm), wrk
      
      z = 0.0d0
      
      do i = 1, ndm
         z(:,i) = q(:,i) * diag(i)
      end do

      write(*,*) 'use: nrec',nrec

!    OPEN(UNIT=78,FILE='hmlt.off'//fl,STATUS='OLD',ACCESS='SEQUENTIAL',&
!         FORM='UNFORMATTED')

      open(unit=78, file='hmlt.off'//fl, status='old', form='unformatted')

      do i = 1, nrec
         read(78) offdiag, indi, indip, nw
         do k = 1, nw
            z(:,indi(k)) = z(:,indi(k)) + q(:,indip(k)) * offdiag(k)
            z(:,indip(k)) = z(:,indip(k)) + q(:,indi(k)) * offdiag(k)
         end do
      end do

      close(unit=78)
      
    end subroutine mtxhq_gg

!#######################################################################

    subroutine readdavvc(nvec,evec,rvec)
! Reads Davidson vectors from disk
      implicit none

      integer, intent(in)                                :: nvec
      double precision, dimension(ndm,nvec), intent(out) :: rvec
      double precision, dimension(nvec), intent(out)     :: evec      
      integer                                            :: i,nr,nfl,count

      nfl=77
      
      open(unit=nfl, file='dav_vecs.init', status='old',access='sequential', form='unformatted')

      count=0
      do i=1,nvec
         read(nfl,end=77) nr,evec(i),rvec(:,i)
         count=count+1
      end do

77    close(nfl)
      write(6,*) nvec,' vectors requested - ', count, 'found on disc'
      
    end subroutine readdavvc

!#######################################################################

  end module davmod
