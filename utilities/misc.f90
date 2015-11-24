module misc

!!$Contains miscellaneous ancillary functions and subroutines
  
  use constants
  use parameters  
  use channels
  
  implicit none

contains
!!$------------------------------------  
  integer function kdelta(k,kpr)
    
    integer, intent(in) :: k,kpr
    
    kdelta=0

    if (k .eq. kpr) kdelta=1
    
  end function kdelta
!!$------------------------------------------------------------  

  logical function lkdelta(k,kpr)

    integer, intent(in) :: k,kpr

    lkdelta=.false.

    if (k .eq. kpr) lkdelta=.true.

  end function lkdelta

  subroutine diagonalise(nsout,ndim,arr)
    
    integer, intent(in) :: ndim,nsout
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: info,lwork,i,j
    real(d), dimension(ndim) :: w
    real(d), dimension(:), allocatable :: work
    
    integer, dimension(ndim) :: indx
    real(d), dimension(ndim) :: coeff
    
    real(d) :: dnrm2, nrm
    
    external dsyev 

100 FORMAT(60("-"),/)
101 FORMAT(3(A14,2x),/)
102 FORMAT(I10,2x,2(F10.6,2x),5(F8.6,"(",I4,")",1x))

    write(ilog,*) "Starting diagonalisation"
    write(ilog,100)
    lwork=3*ndim
    allocate(work(lwork))

    call dsyev("V","L", ndim, arr, ndim, w, work,lwork, info)
    
    write(ilog,101) "Eigenval","Ev","a.u."
    write(ilog,100)
    do i=1,nsout
       coeff(:)=arr(:,i)**2
!!$       call dsortqx("D",ndim,coeff(:),1,indx(:))
       call dsortindxa1("D",ndim,coeff(:),indx(:))
       
       write(ilog,102) i,w(i),w(i)*eh2ev,(coeff(indx(j)),indx(j),j=1,5)
    end do
    
    deallocate(work)
    
  end subroutine diagonalise
!!$---------------------------------------------------

  subroutine vdiagonalise(ndim,arr,evector)

    implicit none

    integer, intent(in)                          :: ndim
    integer                                      :: info,lwork,i,j
    integer, dimension(ndim)                     :: indx
    integer*8                                    :: lworkl,ndiml
    real(d), dimension(ndim), intent(inout)      :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr
    real(d), dimension(ndim)                     :: w
    real(d), dimension(:), allocatable           :: work
    real(d), dimension(ndim)                     :: coeff

    external dsyev

    lwork=3*ndim
    allocate(work(lwork))
    lworkl=int(lwork,lng)
    ndiml=int(ndim,lng)

    call dsyev("V","L",ndiml,arr(:,:),ndiml,w,work(:),lworkl,info)

    if (info.ne.0) then
       write(ilog,'(/,2x,a,/)') 'In subroutine vdiagonalise: &
            diagonalisation of the Hamiltonian matrix failed'
       STOP
    endif

    evector(:)=w(:)

    deallocate(work)

  end subroutine vdiagonalise

!!$----------------------------------------------------

  subroutine get_indices(col,a,b,j,k,spin)
    
    integer, dimension(7),intent(in) :: col
    integer, intent(out) :: a,b,j,k,spin 
    
    spin=col(2)
    j=col(3)
    k=col(4)
    a=col(5)
    b=col(6)
    
  end subroutine get_indices

!!$---------------------------------------------------------------------

  subroutine get_indices1(col,spin,a,b,j,k,type)
    
    integer, dimension(7),intent(in) :: col
    integer, intent(out) :: a,b,j,k,spin,type 
    
    spin=col(2)
    j=col(3)
    k=col(4)
    a=col(5)
    b=col(6)
    type=col(7)
    
  end subroutine get_indices1

!----------------------------------------------------------
!----------------------------------------------------------

  subroutine fill_indices(col,cnf,spin,a,b,j,k,type)
    
    integer, dimension(7),intent(out) :: col
    integer, intent(in) :: a,b,j,k,spin,cnf,type
    
    col(1)=cnf
    col(2)=spin
    col(3)=j
    col(4)=k
    col(5)=a
    col(6)=b
    col(7)=type
    
  end subroutine fill_indices

!-------------------------------------

  real(d) function dsp(ndim,vec1,vec2)
    
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: vec1,vec2
    
    integer :: i
    
    dsp=0._d
    
    do i=1, ndim
       dsp=dsp+vec1(i)*vec2(i)
    end do
    
  end function dsp
!!$---------------------------------------------
  subroutine dsortindxa1(order,ndim,arrin,indx)

    use constants
    implicit none

    character(1), intent(in) :: order
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: arrin
    integer, dimension(ndim), intent(inout) :: indx
    
    integer :: i,l,ir,indxt,j
    real(d) :: q
!!$ The subroutine is taken from the NR p233, employs heapsort.

    do i= 1,ndim
       indx(i)=i
    end do
    
    l=ndim/2+1
    ir=ndim
    
    if(order .eq. 'D') then
       
10     continue
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
       
20     if(j .le. ir) then
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
       
100    continue
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
       
200    if(j .le. ir) then
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
       
  end subroutine dsortindxa1


!!$-------------------------------------------------

  subroutine dsortindxa(order,ndim,arr,indarr)
    
    character(1), intent(in) :: order
    integer, intent(in) :: ndim
    integer, dimension(ndim), intent(out) :: indarr
    real(d), dimension(ndim), intent(in) :: arr
    
    integer :: i, ifail
    integer, dimension(ndim) :: irank
    
!    external M01DAF,M01EBF
    
    ifail=0
!    call M01DAF (arr(:),1,ndim,order,irank(:),ifail)
    do i=1,ndim
       indarr(i)=i
    end do
!    call M01EBF (indarr(:),1,ndim,irank(:),ifail)
    
  end subroutine dsortindxa
!!$------------------------------------------------------
  
  subroutine table1(ndim1,ndim2,en,vspace)
    
    integer, intent(in) :: ndim1,ndim2
    real(d), dimension(ndim2), intent(in) :: en
    real(d), dimension(ndim1,ndim2), intent(in) :: vspace

    real(d), dimension(:), allocatable :: coeff
    integer, dimension(:), allocatable :: indx
    integer :: i,j

    allocate(coeff(ndim1),indx(ndim1))
    
100 FORMAT(60("-"),/)    
101 FORMAT(3(A10,2x),/)
102 FORMAT(I10,x,"|",2(F10.6,2x),"|",x,5(F8.6,"(",I4,")",1x))    
  
    write(ilog,101) "Eigenval","a.u.","Ev"
    write(ilog,100)
    do i=1,ndim2
       coeff(:)=vspace(:,i)**2
       call dsortindxa1("D",ndim1,coeff(:),indx(:))
       write(ilog,102) i,en(i),en(i)*eh2ev,(coeff(indx(j)),indx(j),j=1,5)
    end do 

    deallocate(coeff,indx)

  end subroutine table1

!!$-----------------------------------------------------------

  subroutine table2(ndim1,ndim2,en,vspace,tmvec,osc_str,&
       kpq,kpqdim2,flag)
    
    integer, intent(in) :: ndim1,ndim2
    real(d), dimension(ndim2), intent(in) :: en,tmvec,osc_str
    real(d), dimension(ndim1,ndim2), intent(in) :: vspace

    real(d), dimension(:), allocatable :: coeff
    integer, dimension(:), allocatable :: indx
    integer :: i,j,nlim,k

    integer                           :: kpqdim2,iout
    integer, dimension(7,0:kpqdim2-1) :: kpq

    character(len=1)  :: flag
    character(len=70) :: filename
    
100 FORMAT(60("-"),/)    
101 FORMAT(4(A10,2x),/)
102 FORMAT(I10,x,"|",4(F10.5,2x),"|",x,5(F8.6,"(",I7,")",1x)) 
    
    allocate(coeff(ndim1),indx(ndim1))
    
    nlim=5
    if (nlim .gt. ndim1) nlim=ndim1
  
    iout=127

    if (flag.eq.'i') then
       filename='davstates.dat'
    else if (flag.eq.'f') then
       filename='davstates_f.dat'
    endif
    
    open(iout,file=filename,form='formatted',status='unknown')

    write(ilog,'(/,2x,a,2x,F14.8)') 'Ground state MP2 energy:',ehf+e_mp2
    write(iout,'(/,2x,a,2x,F14.8)') 'Ground state MP2 energy:',ehf+e_mp2

    do i=1,ndim2
       coeff(:)=vspace(:,i)**2
       call dsortindxa1("D",ndim1,coeff(:),indx(:))
       coeff(:)=vspace(:,i)
       call wrstateinfo(i,indx,coeff,kpq,kpqdim2,en(i),tmvec(i),&
            osc_str(i),ndim1,iout)
    end do

    close(iout)
    
    deallocate(coeff,indx)

  end subroutine table2

!#######################################################################

  subroutine wrstateinfo(i,indx,coeff,kpq,kpqdim2,en,tmvec,osc_str,&
       ndim1,iout)

    implicit none

    integer                           :: i,k,ndim1,kpqdim2,ilbl,iout
    integer, dimension(ndim1)         :: indx
    integer, dimension(7,0:kpqdim2-1) :: kpq
    real(d), dimension(ndim1)         :: coeff
    real(d)                           :: en,tmvec,osc_str
    real(d), parameter                :: tol=0.05d0
    character(len=120)                :: fmat
    character(len=2)                  :: spincase
    
!-----------------------------------------------------------------------
! State energy in a.u.
!-----------------------------------------------------------------------
    if (i.lt.10) then
       fmat='(2/,2x,a,3x,i1,a1,i1,2x,a,5x,F14.8)'
    else
       fmat='(2/,2x,a,2x,i2,a1,i1,2x,a,5x,F14.8)'
    endif
    
    write(ilog,fmat) 'State',i,'.',nirrep,'Energy:',en+ehf+e_mp2
    write(iout,fmat) 'State',i,'.',nirrep,'Energy:',en+ehf+e_mp2
    
!-----------------------------------------------------------------------
! Excitation energy in eV
!-----------------------------------------------------------------------
    if (en*eh2ev.lt.10.0d0) then
       write(ilog,'(2x,a,2x,F10.5)') &
            'Excitation Energy (eV):',en*eh2ev
       write(iout,'(2x,a,2x,F10.5)') &
            'Excitation Energy (eV):',en*eh2ev
    else
       write(ilog,'(2x,a,3x,F10.5)') &
            'Excitation Energy (eV):',en*eh2ev
       write(iout,'(2x,a,3x,F10.5)') &
            'Excitation Energy (eV):',en*eh2ev
    endif

!-----------------------------------------------------------------------
! Transition dipole moment along the chosen direction
!-----------------------------------------------------------------------
    if (ltdm_gs2i) then
       write(ilog,'(2x,a,F10.5)') 'Transition Dipole Moment:',tmvec
       write(iout,'(2x,a,F10.5)') 'Transition Dipole Moment:',tmvec
    endif

!-----------------------------------------------------------------------
! Oscillator strength
!-----------------------------------------------------------------------
    if (ltdm_gs2i) then
       write(ilog,'(2x,a,5x,F10.5)') 'Oscillator Strength:',osc_str
       write(iout,'(2x,a,5x,F10.5)') 'Oscillator Strength:',osc_str
    endif

!-----------------------------------------------------------------------
! Dominant configurations
!-----------------------------------------------------------------------
    write(ilog,'(/,2x,a,/)') 'Dominant Configurations:'
    write(ilog,'(2x,29a)') ('*',k=1,29)
    write(ilog,'(3x,a)') 'j   k -> a  b        C_jkab'
    write(ilog,'(2x,29a)') ('*',k=1,29)
    write(iout,'(/,2x,a,/)') 'Dominant Configurations:'
    write(iout,'(2x,29a)') ('*',k=1,29)
    write(iout,'(3x,a)') 'j   k -> a  b        C_jkab'
    write(iout,'(2x,29a)') ('*',k=1,29)

    do k=1,50
       ilbl=indx(k)
       if (abs(coeff(ilbl)).ge.tol) then
          if (kpq(4,ilbl).eq.-1) then
             ! Single excitations
             write(ilog,'(3x,i2,4x,a2,1x,i2,9x,F8.5)') &
                  kpq(3,ilbl),'->',kpq(5,ilbl),coeff(ilbl)
             write(iout,'(3x,i2,4x,a2,1x,i2,9x,F8.5)') &
                  kpq(3,ilbl),'->',kpq(5,ilbl),coeff(ilbl)
          else
             ! Double excitations
             if (kpq(3,ilbl).ne.kpq(4,ilbl).and.kpq(5,ilbl).ne.kpq(6,ilbl)) then
                ! a|=b, i|=j
                spincase=getspincase(ilbl,kpq,kpqdim2)
                write(ilog,'(3x,2(i2,1x),a2,2(1x,i2),2x,a2,2x,F8.5)') &
                     kpq(3,ilbl),kpq(4,ilbl),'->',kpq(5,ilbl),&
                     kpq(6,ilbl),spincase,coeff(ilbl)
                write(iout,'(3x,2(i2,1x),a2,2(1x,i2),2x,a2,2x,F8.5)') &
                     kpq(3,ilbl),kpq(4,ilbl),'->',kpq(5,ilbl),&
                     kpq(6,ilbl),spincase,coeff(ilbl)
             else
                ! a=b,  i=j
                ! a|=b, i=j
                ! a=b,  i=|j
                write(ilog,'(3x,2(i2,1x),a2,2(1x,i2),6x,F8.5)') &
                     kpq(3,ilbl),kpq(4,ilbl),'->',kpq(5,ilbl),&
                     kpq(6,ilbl),coeff(ilbl)
                write(iout,'(3x,2(i2,1x),a2,2(1x,i2),6x,F8.5)') &
                     kpq(3,ilbl),kpq(4,ilbl),'->',kpq(5,ilbl),&
                     kpq(6,ilbl),coeff(ilbl)
             endif
          endif
       endif
    enddo

    write(ilog,'(2x,29a)') ('*',k=1,29)
    write(iout,'(2x,29a)') ('*',k=1,29)

    return

  end subroutine wrstateinfo

!#######################################################################

  function getspincase(ilbl,kpq,kpqdim2) result(spinlbl)

    implicit none

    integer, dimension(7,0:kpqdim2-1) :: kpq
    integer                           :: ilbl,kpqdim2,lim,i
    character(len=2)                  :: spinlbl

    lim=0
    do i=1,5
       lim=lim+kpq(i,0)
    enddo

    if (ilbl.le.lim) then
       spinlbl='I'
    else
       spinlbl='II'
    endif

    return

  end function getspincase

!#######################################################################

  subroutine wrstateinfo_neutral(ndim,kpq,kpqdim2,filename,nstates)

    use iomod, only: freeunit

    implicit none

    integer                            :: ndim,kpqdim2,ieig,nstates,&
                                          i,k,itmp,ilbl
    integer, dimension(7,0:kpqdim2-1)  :: kpq
    integer, dimension(:), allocatable :: indx
    real(d)                            :: ener
    real(d), dimension(:), allocatable :: coeff,coeffsq
    real(d), parameter                 :: tol=0.05d0
    character(len=36)                  :: filename
    character(len=120)                 :: fmat
    character(len=2)                   :: spincase

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(coeff(ndim))
    allocate(coeffsq(ndim))
    allocate(indx(ndim))

!-----------------------------------------------------------------------
! Open the eigenpair file
!-----------------------------------------------------------------------
    call freeunit(ieig)
    open(unit=ieig,file=filename,status='unknown',access='sequential',&
         form='unformatted')

!-----------------------------------------------------------------------
! Read the eigenpairs and output the state energies and dominant
! configurations
!-----------------------------------------------------------------------
    do i=1,nstates

       ! Read the next eigenpair from file
       read(ieig) itmp,ener,coeff
       coeffsq=coeff**2
       call dsortindxa1("D",ndim,coeffsq,indx)

       ! State energy in a.u.
       if (i.lt.10) then
          fmat='(2/,2x,a,3x,i1,a1,i1,2x,a,5x,F14.8)'
       else
          fmat='(2/,2x,a,2x,i2,a1,i1,2x,a,5x,F14.8)'
       endif
       write(ilog,fmat) 'State',i,'.',nirrep,'Energy:',ener+ehf+e_mp2

       ! Excitation energy in eV
       if (ener*eh2ev.lt.10.0d0) then
          write(ilog,'(2x,a,2x,F10.5)') &
               'Excitation Energy (eV):',ener*eh2ev
       else
          write(ilog,'(2x,a,3x,F10.5)') &
               'Excitation Energy (eV):',ener*eh2ev
       endif
       
       ! Dominant configurations
       write(ilog,'(/,2x,a,/)') 'Dominant Configurations:'
       write(ilog,'(2x,29a)') ('*',k=1,29)
       write(ilog,'(3x,a)') 'j   k -> a  b        C_jkab'
       write(ilog,'(2x,29a)') ('*',k=1,29)
       do k=1,50
          ilbl=indx(k)
          if (abs(coeff(ilbl)).lt.tol) cycle
          if (kpq(4,ilbl).eq.-1) then
             ! Single excitations
             write(ilog,'(3x,i2,4x,a2,1x,i2,9x,F8.5)') &
                  kpq(3,ilbl),'->',kpq(5,ilbl),coeff(ilbl)
          else
             ! Double excitations
             if (kpq(3,ilbl).ne.kpq(4,ilbl).and.kpq(5,ilbl).ne.kpq(6,ilbl)) then
                ! a|=b, i|=j
                spincase=getspincase(ilbl,kpq,kpqdim2)
                write(ilog,'(3x,2(i2,1x),a2,2(1x,i2),2x,a2,2x,F8.5)') &
                     kpq(3,ilbl),kpq(4,ilbl),'->',kpq(5,ilbl),&
                     kpq(6,ilbl),spincase,coeff(ilbl)
             else
                ! a=b,  i=j
                ! a|=b, i=j
                ! a=b,  i=|j
                write(ilog,'(3x,2(i2,1x),a2,2(1x,i2),6x,F8.5)') &
                     kpq(3,ilbl),kpq(4,ilbl),'->',kpq(5,ilbl),&
                     kpq(6,ilbl),coeff(ilbl)
             endif
          endif
       enddo
       write(ilog,'(2x,29a)') ('*',k=1,29)

    enddo

!-----------------------------------------------------------------------
! Close the eigenpair file
!-----------------------------------------------------------------------
    close(ieig)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(coeff)
    deallocate(coeffsq)
    deallocate(indx)

    return

  end subroutine wrstateinfo_neutral

!#######################################################################

  subroutine wrstateinfo_cation(ndim,kpq,kpqdim2,filename,nstates)

    use iomod, only: freeunit

    implicit none

    integer                              :: ndim,kpqdim2,ieig,nstates,&
                                            i,k,itmp,ilbl,iout
    integer, dimension(7,0:kpqdim2-1)    :: kpq
    integer, dimension(:), allocatable   :: indx
    real(d)                              :: ener
    real(d), dimension(:), allocatable   :: coeff,coeffsq
    real(d), parameter                   :: tol=0.05d0
    character(len=36)                    :: filename
    character(len=120)                   :: fmat

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
    allocate(coeff(ndim))
    allocate(coeffsq(ndim))
    allocate(indx(ndim))

!-----------------------------------------------------------------------
! Open files
!-----------------------------------------------------------------------
    call freeunit(ieig)
    open(unit=ieig,file=filename,status='unknown',access='sequential',&
         form='unformatted')

    call freeunit(iout)
    open(unit=iout,file='ipadcstates.dat',form='formatted',&
         status='unknown')

!-----------------------------------------------------------------------
! Read the eigenpairs and output the state energies and dominant
! configurations
!-----------------------------------------------------------------------
    do i=1,nstates

       ! Read the next eigenpair from file
       read(ieig) itmp,ener,coeff
       coeffsq=coeff**2
       call dsortindxa1("D",ndim,coeffsq,indx)

       ! State energy in a.u.
       if (i.lt.10) then
          fmat='(2/,2x,a,3x,i1,a1,i1,2x,a,5x,F14.8)'
       else if (i.lt.100) then
          fmat='(2/,2x,a,2x,i2,a1,i1,2x,a,5x,F14.8)'
       else 
          fmat='(2/,2x,a,1x,i3,a1,i1,2x,a,5x,F14.8)'
       endif
       write(iout,fmat) 'State',i,'.',nirrep,'Energy:',ener+ehf+e_mp2

       ! Excitation energy in eV
       if (ener*eh2ev.lt.10.0d0) then
          write(iout,'(2x,a,2x,F10.5)') &
               'Excitation Energy (eV):',ener*eh2ev
       else
          write(iout,'(2x,a,3x,F10.5)') &
               'Excitation Energy (eV):',ener*eh2ev
       endif
       
       ! Dominant configurations
       write(iout,'(/,2x,a,/)') 'Dominant Configurations:'
       write(iout,'(2x,25a)') ('*',k=1,25)
       write(iout,'(3x,a)') 'j   k -> a  b    C_jkab'
       write(iout,'(2x,25a)') ('*',k=1,25)
       do k=1,50          
          ilbl=indx(k)
          if (abs(coeff(ilbl)).lt.tol) cycle
          if (kpq(4,ilbl).eq.-1) then
             ! Single excitations
             write(iout,'(3x,i2,4x,a2,1x,i2,5x,F8.5)') &
                  kpq(3,ilbl),'->',kpq(5,ilbl),coeff(ilbl)
          else
             ! Double excitations
             write(iout,'(3x,2(i2,1x),a2,2(1x,i2),2x,F8.5)') &
                  kpq(3,ilbl),kpq(4,ilbl),'->',kpq(5,ilbl),&
                  kpq(6,ilbl),coeff(ilbl)
          endif
       enddo
       write(iout,'(2x,25a)') ('*',k=1,25)

    enddo

!-----------------------------------------------------------------------
! Close files
!-----------------------------------------------------------------------
    close(ieig)
    close(iout)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
    deallocate(coeff)
    deallocate(coeffsq)
    deallocate(indx)

    return

  end subroutine wrstateinfo_cation

!#######################################################################

 subroutine mat_vec_multiply_SYM(ndim,A,vec,fin)

    integer, intent(in) :: ndim

    real(d), dimension(ndim,ndim), intent(in) :: A
    real(d), dimension(ndim), intent(in) :: vec
    real(d), dimension(ndim), intent(inout) :: fin

    integer :: LDA
    real(d) :: alpha    
    real(d) :: beta


    external dsymv

    write(ilog,*) "Starting matrix-vector multiplication"

    LDA=max(1,ndim)

    alpha=1._d
    beta=0._d


    call dsymv("L",ndim,alpha,A,LDA,vec,1,beta,fin,1)


  end subroutine mat_vec_multiply_SYM

 subroutine mat_vec_multiply(ndimiA,ndimjA,A,vec,fin)

    integer, intent(in) :: ndimiA
    integer, intent(in) :: ndimjA

    real(d), dimension(ndimiA,ndimjA), intent(in) :: A
    real(d), dimension(ndimjA), intent(in) :: vec
    real(d), dimension(ndimiA), intent(inout) :: fin

    integer :: LDA
    real(d) :: alpha    
    real(d) :: beta


    external dgemv

    write(ilog,*) "Starting matrix-vector multiplication"

    LDA=max(1,ndimiA)

    alpha=1._d
    beta=0._d


    call dgemv("N",ndimiA,ndimjA,alpha,A,LDA,vec,1,beta,fin,1)

  end subroutine mat_vec_multiply




  subroutine read_density_matrix


  integer ::  k,a
  integer :: i
  real*8, dimension(nBas,nOcc) :: density_matrix

  OPEN(unit=77,file="density.dat",status='OLD',access='SEQUENTIAL',form='FORMATTED')
  write(ilog,*) "nBas from read_density_matrix",nBas

  do i=1,nBas
  read(77) k,a,density_matrix(a,k)
  end do

  close(77)

  end subroutine read_density_matrix
 


end module misc
