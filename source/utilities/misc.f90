module misc

  ! Contains miscellaneous ancillary functions and subroutines
  
  use constants
  use parameters  
  use channels
  
  implicit none

contains
  
!#######################################################################

  function atomic_number(lbl) result(num)

    real*8           :: num
    character(len=*) :: lbl
    
    if (lbl.eq.'h') then
       num=1.0d0
    else if (lbl.eq.'he') then
       num=2.0d0
    else if (lbl.eq.'c') then
       num=6.0d0
    else if (lbl.eq.'n') then
       num=7.0d0
    else if (lbl.eq.'o') then
       num=8.0d0
    else if (lbl.eq.'f') then
       num=9.0d0
    else if (lbl.eq.'ne') then
       num=10.0d0
    else if (lbl.eq.'s') then
       num=16.0d0
    else if (lbl.eq.'cl') then
       num=17.0d0
    else
       write(6,'(/,2(2x,a),/)') 'Atomic no. not known for atom:',lbl
    endif
    
    return
    
  end function atomic_number

!#######################################################################

  function uppercase(string1) result(string2)

    implicit none
    
    integer                     :: dim,i,j
    character(len=*)            :: string1
    character(len=len(string1)) :: string2
    
    do i = 1, len(string1)
       j = iachar(string1(i:i))
       if (j>= iachar("a") .and. j<=iachar("z") ) then
          string2(i:i) = achar(iachar(string1(i:i))-32)
       else
          string2(i:i) = string1(i:i)
       endif
    enddo

    return

  end function uppercase

!#######################################################################
  
  subroutine getcom(xcom)

    use parameters

    implicit none

    integer              :: i,j
    real*8, dimension(3) :: xcom
    real*8               :: m,tmass

    xcom=0.0d0
    tmass=0.0d0

    do i=1,natm
       m=mass(atlbl(i))
       tmass=tmass+m
       do j=1,3
          xcom(j)=xcom(j)+xcoo(i*3-3+j)*m
       enddo
    enddo
    
    xcom=xcom/tmass

    return

  end subroutine getcom

!#######################################################################

  function mass(label)

    implicit none

    real*8           :: mass
    character(len=*) :: label

    if (label.eq.'h') then
       mass=1.00794d0
    else if (label.eq.'he') then
       mass=4.002602d0
    else if (label.eq.'c') then
       mass=12.0107d0
    else if (label.eq.'n') then
       mass=14.0067d0      
    else if (label.eq.'o') then
       mass=15.9994d0
    else if (label.eq.'f') then
       mass=18.998403d0
    else if (label.eq.'ne') then
       mass=20.1797d0
    else if (label.eq.'s') then
       mass=32.065d0
    else if (label.eq.'cl') then
       mass=35.453d0
    else
       write(6,'(2(2x,a))') 'Unknown atom type:',trim(label)
       STOP
    endif
    
    return

  end function mass

!#######################################################################

  subroutine get_vdwr(gam,vdwr,natom)

    use constants
    use channels
    use parameters
    use iomod
    use parsemod, only: lowercase
    use import_gamess
    
    implicit none
    
    integer                   :: natom,i
    real(d), dimension(natom) :: vdwr
    real(d), parameter        :: ang2bohr=1.889725989d0
    character(len=20)         :: name
    type(gam_structure)       :: gam
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------   
    vdwr=0.0d0

!----------------------------------------------------------------------
! Fill in the van der Waals radius array (units of Bohr)
! Van der Waals radii taken from Mantina et al., JPCA, 113, 5806 (2009)
!----------------------------------------------------------------------
    ! Van der Waals radii in Angstrom
    do i=1,natom
       
       ! Atom name
       name=gam%atoms(i)%name
       call lowercase(name)

       ! Skip dummy atoms
       if (name.eq.'x') cycle

       ! Assign the Van der Waals radius for the current atom
       if (name.eq.'h') then
          vdwr(i)=1.10d0

       else if (name.eq.'he') then
          vdwr(i)=1.40d0

       else if (name.eq.'li') then
          vdwr(i)=1.81d0

       else if (name.eq.'be') then
          vdwr(i)=1.53d0

       else if (name.eq.'b') then
          vdwr(i)=1.92d0

       else if (name.eq.'c') then
          vdwr(i)=1.70d0

       else if (name.eq.'n') then
          vdwr(i)=1.55d0

       else if (name.eq.'o') then
          vdwr(i)=1.52d0

       else if (name.eq.'f') then
          vdwr(i)=1.47d0
                    
       else if (name.eq.'ne') then
          vdwr(i)=1.54d0

       else if (name.eq.'na') then
          vdwr(i)=2.27d0

       else if (name.eq.'mg') then
          vdwr(i)=1.73d0

       else if (name.eq.'al') then
          vdwr(i)=1.84d0

       else if (name.eq.'si') then
          vdwr(i)=2.10d0

       else if (name.eq.'p') then
          vdwr(i)=1.80d0

       else if (name.eq.'s') then
          vdwr(i)=1.80d0

       else if (name.eq.'cl') then
          vdwr(i)=1.75d0

       else if (name.eq.'ar') then
          vdwr(i)=1.88d0

       else if (name.eq.'k') then
          vdwr(i)=2.75d0

       else if (name.eq.'ca') then
          vdwr(i)=2.31d0

       else if (name.eq.'ga') then
          vdwr(i)=1.87d0

       else if (name.eq.'ge') then
          vdwr(i)=2.11d0

       else if (name.eq.'as') then
          vdwr(i)=1.85d0

       else if (name.eq.'se') then
          vdwr(i)=1.90d0

       else if (name.eq.'br') then
          vdwr(i)=1.83d0

       else if (name.eq.'kr') then
          vdwr(i)=2.02d0

       else if (name.eq.'rb') then
          vdwr(i)=3.03d0

       else if (name.eq.'sr') then
          vdwr(i)=2.49d0

       else if (name.eq.'in') then
          vdwr(i)=1.93d0

       else if (name.eq.'sn') then
          vdwr(i)=2.17d0

       else if (name.eq.'sb') then
          vdwr(i)=2.06d0

       else if (name.eq.'te') then
          vdwr(i)=2.06d0

       else if (name.eq.'i') then
          vdwr(i)=1.98d0

       else if (name.eq.'xe') then
          vdwr(i)=2.16d0

       else if (name.eq.'cs') then
          vdwr(i)=3.43d0

       else if (name.eq.'ba') then
          vdwr(i)=2.68d0

       else if (name.eq.'tl') then
          vdwr(i)=1.96d0

       else if (name.eq.'pb') then
          vdwr(i)=2.02d0

       else if (name.eq.'bi') then
          vdwr(i)=2.07d0

       else if (name.eq.'po') then
          vdwr(i)=1.97d0

       else if (name.eq.'at') then
          vdwr(i)=2.02d0

       else if (name.eq.'rn') then
          vdwr(i)=2.20d0

       else if (name.eq.'fr') then
          vdwr(i)=3.48d0

       else if (name.eq.'ra') then
          vdwr(i)=2.83d0

       else
          errmsg='Currently CAPs are not supported for the element '&
               //trim(name)
          call error_control
          
       endif
       
    enddo

    ! Convert to Bohr
    vdwr=vdwr*ang2bohr

    return
    
  end subroutine get_vdwr
    
!#######################################################################

  integer function kdelta(k,kpr)
    
    integer, intent(in) :: k,kpr
    
    kdelta=0

    if (k .eq. kpr) kdelta=1
    
  end function kdelta

!#######################################################################

  logical function lkdelta(k,kpr)

    integer, intent(in) :: k,kpr

    lkdelta=.false.

    if (k .eq. kpr) lkdelta=.true.

  end function lkdelta

!#######################################################################

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

!#######################################################################

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

!#######################################################################

  subroutine get_indices(col,a,b,j,k,spin)
    
    integer, dimension(7),intent(in) :: col
    integer, intent(out) :: a,b,j,k,spin 
    
    spin=col(2)
    j=col(3)
    k=col(4)
    a=col(5)
    b=col(6)
    
  end subroutine get_indices

!#######################################################################

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

!#######################################################################

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

!#######################################################################

  real(d) function dsp(ndim,vec1,vec2)
    
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: vec1,vec2
    
    integer :: i
    
    dsp=0._d
    
    do i=1, ndim
       dsp=dsp+vec1(i)*vec2(i)
    end do
    
  end function dsp

!#######################################################################

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

    if (ndim.eq.1) return

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

!#######################################################################

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

!#######################################################################

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
       call dsortindxa1('D',ndim1,coeff(:),indx(:))
       write(ilog,102) i,en(i),en(i)*eh2ev,(coeff(indx(j)),indx(j),j=1,5)
    end do 

    deallocate(coeff,indx)

  end subroutine table1

!#######################################################################

  subroutine table2(ndim1,ndim2,en,vspace,tmvec,osc_str,&
       kpq,kpqdim2,flag)

    use iomod, only: freeunit

    implicit none
    
    integer, intent(in)                         :: ndim1,ndim2
    integer, dimension(:), allocatable          :: indx
    integer                                     :: i,j,k
    integer                                     :: kpqdim2,iout
    integer, dimension(7,0:kpqdim2-1)           :: kpq
    real(d), dimension(ndim2), intent(in)       :: en,tmvec,osc_str
    real(d), dimension(ndim1,ndim2), intent(in) :: vspace
    real(d), dimension(:), allocatable          :: coeff
    character(len=1)                            :: flag
    character(len=70)                           :: filename

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(coeff(ndim1),indx(ndim1))

!----------------------------------------------------------------------
! Open the davstates file
!----------------------------------------------------------------------
    call freeunit(iout)

    if (flag.eq.'i') then
       filename='davstates.dat'
    else if (flag.eq.'f') then
       filename='davstates_f.dat'
    endif    

    open(iout,file=filename,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Output the excited state information to the log and davstate files
!----------------------------------------------------------------------    
    ! MP2 energy
    write(ilog,'(/,2x,a,2x,F14.8)') 'Ground state MP2 energy:',ehf+e_mp2
    write(iout,'(/,2x,a,2x,F14.8)') 'Ground state MP2 energy:',ehf+e_mp2

    ! ADC state information
    do i=1,ndim2
       coeff(:)=vspace(:,i)**2
       call dsortindxa1('D',ndim1,coeff(:),indx(:))
       coeff(:)=vspace(:,i)
       call wrstateinfo(i,indx,coeff,kpq,kpqdim2,en(i),tmvec(i),&
            osc_str(i),ndim1,iout,flag)
    enddo

!----------------------------------------------------------------------
! Close the davstates file
!----------------------------------------------------------------------
    close(iout)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(coeff,indx)
    
100 format(60("-"),/)
101 format(4(A10,2x),/)
102 format(I10,x,"|",4(F10.5,2x),"|",x,5(F8.6,"(",I7,")",1x)) 

    return
    
  end subroutine table2

!#######################################################################

  subroutine wrstateinfo(i,indx,coeff,kpq,kpqdim2,en,tmvec,osc_str,&
       ndim1,iout,flag)

    implicit none

    integer                           :: i,k,ndim1,kpqdim2,ilbl,iout
    integer, dimension(ndim1)         :: indx
    integer, dimension(7,0:kpqdim2-1) :: kpq
    real(d), dimension(ndim1)         :: coeff
    real(d)                           :: en,tmvec,osc_str
    real(d), parameter                :: tol=0.05d0
    character(len=120)                :: fmat
    character(len=2)                  :: spincase
    character(len=1)                  :: flag
    
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
! Dipole moment along the chosen direction
!-----------------------------------------------------------------------
    if (ldipole) then
       if (flag.eq.'i'.and.statenumber.gt.0.or.method.lt.0) then
          write(ilog,'(2x,a,11x,F10.5)') 'Dipole Moment:',dipmom(i)
          write(iout,'(2x,a,11x,F10.5)') 'Dipole Moment:',dipmom(i)
       else if (flag.eq.'f') then
          write(ilog,'(2x,a,11x,F10.5)') 'Dipole Moment:',dipmom_f(i)
          write(iout,'(2x,a,11x,F10.5)') 'Dipole Moment:',dipmom_f(i)
       endif
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
       call dsortindxa1('D',ndim,coeffsq,indx)

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

    integer                            :: ndim,kpqdim2,ieig,nstates,&
                                          i,k,itmp,ilbl,iout
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
       call dsortindxa1('D',ndim,coeffsq,indx)

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
       write(iout,'(2x,29a)') ('*',k=1,29)
       write(iout,'(3x,a)') 'j   k -> a  b        C_jkab'
       write(iout,'(2x,29a)') ('*',k=1,29)
       do k=1,50          
          ilbl=indx(k)
          if (abs(coeff(ilbl)).lt.tol) cycle
          if (kpq(4,ilbl).eq.-1) then
             ! Single excitations
             write(iout,'(3x,i2,4x,a2,1x,i2,9x,F8.5)') &
                  kpq(3,ilbl),'->',kpq(5,ilbl),coeff(ilbl)
          else
             ! Double excitations
             if (kpq(3,ilbl).ne.kpq(4,ilbl).and.kpq(5,ilbl).ne.kpq(6,ilbl)) then
                ! a|=b, i|=j
                spincase=getspincase(ilbl,kpq,kpqdim2)
                write(iout,'(3x,2(i2,1x),a2,2(1x,i2),2x,a2,2x,F8.5)') &
                     kpq(3,ilbl),kpq(4,ilbl),'->',kpq(5,ilbl),&
                     kpq(6,ilbl),spincase,coeff(ilbl)
             else
                ! a=b,  i=j
                ! a|=b, i=j
                ! a=b,  i=|j
                write(iout,'(3x,2(i2,1x),a2,2(1x,i2),6x,F8.5)') &
                     kpq(3,ilbl),kpq(4,ilbl),'->',kpq(5,ilbl),&
                     kpq(6,ilbl),coeff(ilbl)
             endif
          endif
       enddo
       write(iout,'(2x,29a)') ('*',k=1,29)

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

!#######################################################################

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

!#######################################################################

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
 
!#######################################################################

end module misc
