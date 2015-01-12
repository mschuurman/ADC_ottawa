module lancdiag

  use constants
  use parameters
  use get_matrix
  use fspace

  implicit none

  integer :: main,nsat,mdim,noffdel
  real(d), dimension(:), allocatable :: diag, offdiag
  integer, dimension(:), allocatable :: indi, indj

contains

!!&-------------------------------------------------------------

!!$-----------------------------------------------  
  subroutine master_lancdiag(ndim,noff,flag)
    
    integer, intent(in) :: ndim,noff
    character(2),intent(in) :: flag

    integer :: i,j
    
    external blnczs

!!$flag switches between adc2 and adc2e

!    if (stvc_flg .eq. 0) then

       write(6,*) 'Set stvc_flg to 1'

!!$       allocate(kpq_loc(7,0:nBas**2*nOcc**2))
!!$       kpq_loc(:,0:nBas**2*nOcc**2)=kpq_gl(:,0:nBas**2*nOcc**2)
!!$    
!!$       main=kpq_loc(1,0)
!!$       nsat=ndim-main
!!$
!!$       write(6,*) 'parameter main=',main
!!$       write(6,*) 'parameter nsat=',nsat
!!$
!!$       if (flag .eq. '2s') then
!!$          mtxid='mx2s'   
!!$          call blnczs(mtxid,ndim,main,ncycles,maxmem,memx,mode,nprint,wthr,erange(:),unit,iparm(:),fparm(:),&
!!$               mtxq1_0,mtxhq2s_0)
!!$       end if
!!$
!!$       
!!$       deallocate(kpq_loc) 

!    elseif (stvc_flg .eq. 1) then


    write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
       mdim=ndim
       main=lmain
       noffdel=noff

       write(6,*) 'noffdel=',noffdel
       allocate(diag(ndim),offdiag(noffdel),indi(noffdel),indj(noffdel))
       call read_matrix()
       if (flag .eq. '2s') then       
!       if (flag .eq. 'f') then
          mtxid='mx2s'   
          call blnczs(mtxid,ndim,main,ncycles,maxmem,memx,mode,nprint,wthr,erange(:),unit,iparm(:),fparm(:),&
               MTXQ1_1,mtxhq2s_2)
       end if
       
       deallocate(diag,offdiag,indi,indj)
!    end if

  end subroutine master_lancdiag

!!$----------------------------------------------
!!$----------------------------------------------

!!$  subroutine MTXQ1_0(amx,bmx,i1,i2,work)
!!$
!!$    real(d), dimension(main,main), intent(out) :: amx
!!$    real(d), dimension(nsat,i2-i1+1), intent(out) :: bmx
!!$    real(d), dimension(memx), intent(inout) :: work
!!$    integer, intent(in) :: i1,i2
!!$
!!$    call get_phph_adc2(main,kpq_loc(:,:),amx(:,:))
!!$    call get_ph_2p2h(nsat,i1,i2,kpq_loc(:,:),bmx(:,:))
!!$    
!!$  end subroutine MTXQ1_0
  
!!$-----------------------------------------------------
!!$-----------------------------------------------------

!!$  subroutine MTXHQ2s_0(q,z,nq,work)
!!$    
!!$    integer, intent(in) :: nq
!!$    real(d), dimension(nq,nsat), intent(in) :: q
!!$    real(d), dimension(nq,nsat), intent(out) :: z
!!$    real(d), dimension(memx), intent(inout) :: work
!!$    
!!$    integer :: k
!!$    real(d), dimension(:),allocatable :: ar_diag
!!$!    
!!$    allocate(ar_diag(nsat))
!!$   
!!$    call get_2p2h2p2h_dg2s(nsat,kpq_loc(:,:),ar_diag(:)) 
!!$    do k=1,nsat
!!$       z(:,k)=q(:,k)*ar_diag(k)
!!$    end do
!!$   
!!$   deallocate(ar_diag) 
!!$ end subroutine MTXHQ2s_0

!!$------------------------------------------

  subroutine MTXQ1_1(amx,bmx,i1,i2,work)

    real(d), dimension(main,main), intent(out) :: amx
    real(d), dimension(mdim,1), intent(out) :: bmx
!    real(d), dimension(mdim,i2-i1+1), intent(out) :: bmx
    real(d), dimension(memx), intent(inout) :: work
    integer, intent(in) :: i1,i2

    integer :: i

! Hier kann ich den Startblock fest legen bmx(i,j)=bla
! Frage bleibt, wie man i1, i2 festlegen kann
    write(6,*) 'BBBBBBBBBBBBBBBBBBBBBBBBB'
    amx(:,:)=rzero
    bmx(:,:)=rzero
    write(6,*) 'ratatata'
    bmx(1,1)=rone
    bmx(2:mdim,1)=rzero
!     do i=i1,i2
!        bmx(stvc_lbl(i),i-i1+1)=rone
!     end do
       
  end subroutine MTXQ1_1

!!$-----------------------------------------

  subroutine  MTXHQ2s_1(qv,zv,nq,work)

    integer, intent(in) :: nq
    real(d), dimension(nq,mdim), intent(in) :: qv
    real(d), dimension(nq,mdim), intent(out) :: zv
    real(d), dimension(memx), intent(inout) :: work

    real(d) :: mtrl
    integer :: mark1,ndim,maxbl,nrec,type
    integer :: irec,nlim,i,j, irow,jcol

    real(d), dimension(:), allocatable :: dgl,offdgl
    integer, dimension(:), allocatable :: oi,oj

!!$ Reading diagonal elements

    write(6,*) 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'

    OPEN(UNIT=77,FILE='adc2.diagc',STATUS='OLD',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')

    read(77) ndim,maxbl,nrec,type
    allocate(dgl(ndim),offdgl(maxbl),oi(maxbl),oj(maxbl))
    read(77) dgl(:)

    CLOSE(77)

    do i= 1,mdim
       zv(:,i)=dgl(i)*qv(:,i)
    end do

    OPEN(UNIT=78,FILE='adc2.offc',STATUS='OLD',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')

    do irec= 1,nrec
       read(78) offdgl(:),oi(:),oj(:),nlim

       do i= 1,nlim
          irow=oi(i)
          jcol=oj(i)
          mtrl=offdgl(i)

          zv(:,irow)=zv(:,irow)+mtrl*qv(:,jcol)
          zv(:,jcol)=zv(:,jcol)+mtrl*qv(:,irow)

       end do
    end do

    close(78)
    deallocate(dgl,offdgl,oi,oj)

  end subroutine MTXHQ2s_1

!$$-------------------------------------------

  subroutine  read_matrix()
    
    real(d) :: mtrl
    integer :: mark1,ndim,maxbl,nrec,type
    integer :: irec,nlim,i,count
    
    real(d), dimension(:), allocatable :: dgl,offdgl
    integer, dimension(:), allocatable :: oi,oj
    

!!$ Reading diagonal elements
    
    
    OPEN(UNIT=77,FILE='adc2.diagc',STATUS='OLD',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    
    read(77) ndim,maxbl,nrec,type
    allocate(offdgl(maxbl),oi(maxbl),oj(maxbl))
    read(77) diag(:)
    
    CLOSE(77)

    
    OPEN(UNIT=78,FILE='adc2.offc',STATUS='OLD',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    
    count=0
    do irec= 1,nrec
       read(78) offdgl(:),oi(:),oj(:),nlim
       offdiag(count+1:count+nlim)=offdgl(1:nlim)
       indi(count+1:count+nlim)=oi(1:nlim)
       indj(count+1:count+nlim)=oj(1:nlim)
       count=count+nlim
    end do
    close(78)
    deallocate(offdgl,oi,oj)
    
  end subroutine read_matrix

!!$------------------------------------------

  subroutine  MTXHQ2s_2(qv,zv,nq,work)

    integer, intent(in) :: nq
    real(d), dimension(nq,mdim), intent(in) :: qv
    real(d), dimension(nq,mdim), intent(out) :: zv
    real(d), dimension(memx), intent(inout) :: work

    real(d) :: mtrl
    integer :: mark1,ndim,maxbl,nrec,type
    integer :: irec,nlim,i,j,irow,jcol

    do i= 1,mdim
       zv(:,i)=diag(i)*qv(:,i)
    end do
    
    do i= 1,noffdel
       irow=indi(i)
       jcol=indj(i)
       mtrl=offdiag(i)
       
       zv(:,irow)=zv(:,irow)+mtrl*qv(:,jcol)
       zv(:,jcol)=zv(:,jcol)+mtrl*qv(:,irow)
       
    end do
 
  end subroutine MTXHQ2s_2

end module lancdiag
