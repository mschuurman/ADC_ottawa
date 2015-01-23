module filetools

  use constants
  use parameters
  
  implicit none

contains

!!$---------------------------------------------------
  subroutine wrtdg(untnr,ndim,maxbl,nbuf,type,diag)
  
    integer, intent(in) :: untnr,maxbl,ndim,nbuf,type
    real(d), dimension(ndim), intent(in) :: diag

    integer :: j
    
    write(untnr) maxbl,nbuf
    write(untnr) diag(:)

  end subroutine wrtdg
!!$------------------------------------------------------
  subroutine wrtoffdg(untnr,maxbl,buff,oi,oj,nrec)
    
    integer, intent(in) :: untnr,maxbl,nrec
    integer, dimension(maxbl), intent(in) :: oi,oj
    real(d), dimension(maxbl), intent(in) :: buff

    write(untnr) buff(:),oi(:),oj(:),nrec
    
  end subroutine wrtoffdg
!!$------------------------------------------------------
  subroutine readvct(ndim,fnm,nvecin,evec,rvec,nvecout)

!!$ reads the first nvecin vectors from a file generated by F. Tarantelli's libs

    integer, intent(in) :: ndim,nvecin,fnm
    integer, intent(out) :: nvecout
    real(d), dimension(ndim,nvecin),intent(out) :: rvec
    real(d), dimension(nvecin), intent(out) :: evec

    integer :: i,nr,nfl,count
    character(20) :: fname

    nfl=77
    count=0

    if (fnm .eq. 1) then 
       open(unit=nfl, file=davname,status='OLD',&
            access='SEQUENTIAL',form='UNFORMATTED')
       fname=davname
    elseif(fnm .eq. 2) then
       open(unit=nfl, file=lancname,status='OLD',&
            access='SEQUENTIAL',form='UNFORMATTED')
       fname=lancname
    end if
    
    do i= 1,nvecin
       read(nfl,END=77) nr,evec(i),rvec(:,i)
       count=count+1
    end do

77    close(nfl)
    
    nvecout=count
    write(6,*) nvecin,' Vectors requested - ', count, 'found on disc in the file ', fname

  end subroutine readvct
!!$-------------------------------------------------------------------------------------
  subroutine readvct1(ndim,fnm,i1,i2,rvec,nvecout)

!!$ reads the vectors i1 to i2 from a file generated by F. Tarantelli's libs

    integer, intent(in) :: ndim,i1,i2
    integer, intent(inout) :: fnm
    integer, intent(out) :: nvecout
    real(d), dimension(ndim,i2-i1+1),intent(out) :: rvec

    real(d) :: en
    integer :: i,nr,nfl,count
    character(20) :: fname

    nfl=77
    count=0
    
    if (fnm .eq. 1) then 
       open(unit=nfl, file=davname,status='OLD',&
            access='SEQUENTIAL',form='UNFORMATTED')
       fname=davname
    elseif(fnm .eq. 2) then
       open(unit=nfl, file=lancname,status='OLD',&
            access='SEQUENTIAL',form='UNFORMATTED')
       fname=lancname
    end if
    
    !Rewinding to the right spot
    do i= 1,i1-1
       read(nfl,END=77)
    end do

    do i= i1,i2
       read(nfl,END=77) nr,en,rvec(:,i-i1+1)
       count=count+1
    end do

77    close(nfl)
    
    nvecout=count

  end subroutine readvct1
!!$-------------------------------------------------------------------------------------

   subroutine readen(ndim,fnm,nvecin,evec,nvecout)

!!$ reads the first nvecin energies from a file generated by F. Tarantelli's libs

    integer, intent(in) :: ndim,nvecin,fnm
    integer, intent(out) :: nvecout
    real(d), dimension(nvecin), intent(out) :: evec

    integer :: i,nr,nfl,count
    character(20) :: fname

    nfl=77
    count=0

    if (fnm .eq. 1) then 
       open(unit=nfl, file=davname,status='OLD',&
            access='SEQUENTIAL',form='UNFORMATTED')
       fname=davname
    elseif(fnm .eq. 2) then
       open(unit=nfl, file=lancname,status='OLD',&
            access='SEQUENTIAL',form='UNFORMATTED')
       fname=lancname
    end if
    
    do i= 1,nvecin
       read(nfl,END=77) nr,evec(i)
       count=count+1
    end do

77    close(nfl)
    
    nvecout=count
    write(6,*) nvecin,' Energies requested - ', count, 'found on disc in the file ', fname

  end subroutine readen

!!$---------------------------------------------------------

  subroutine write_vec(ndim,name,vec)
    
    integer, intent(in) :: ndim
    character(4), intent(in) :: name
    real(d), dimension(ndim), intent(in) :: vec

    open(unit=99,file=name,status='UNKNOWN',access='SEQUENTIAL',&
         form='UNFORMATTED')
    write(99) 1,ndim
    write(99) 1
    write(99) 1,0.0_d,vec(:)
    close(99)
    
  end subroutine write_vec
!!$----------------------------------------------------------------------
    subroutine read_vec(ndim,name,vec)

      integer, intent(in) :: ndim
      character(4), intent(in) :: name
      real(d), dimension(ndim), intent(out) :: vec
      integer :: itmp,ndim1
      real(d) :: tmp

      open(unit=99,file=name,status='OLD',access='SEQUENTIAL',&
           form='UNFORMATTED')
      read(99) itmp,ndim1

      if (ndim .ne. ndim1) then
         write(6,*) 'The vector size in ', name,' is wrong', ndim1
         stop
      end if

      read(99) itmp
      read(99) itmp,tmp,vec(:)
      close(99)

    end subroutine read_vec


end module filetools
    
    
