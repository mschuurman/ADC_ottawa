!#######################################################################
! relax2: A program to calculate the 2nd-order relaxation correction
!         to the Hartree-Fock orbital energies following
!         core-excitation
!#######################################################################

  program relax2

    implicit none

!-----------------------------------------------------------------------
! Read the input file
!-----------------------------------------------------------------------
    call read_input

!-----------------------------------------------------------------------
! Run GAMESS if requested
!-----------------------------------------------------------------------
    

  contains

!#######################################################################

    subroutine read_input

      use channels
      use constants
      use relaxmod
      use parsemod
      use iomod

      implicit none
      
      integer :: i,k,n,l

!-----------------------------------------------------------------------
! Set defaults
!-----------------------------------------------------------------------
      ain=''

!-----------------------------------------------------------------------
! I/O channels
!-----------------------------------------------------------------------
      iin=1
      ilog=2

!-----------------------------------------------------------------------
! Determine the input file name
!-----------------------------------------------------------------------
      call getarg(1,ain)

      if (ain.eq.'') then         
         write(6,'(/,a,/)') 'The input file name has not been given'
         STOP
      endif

      k=index(ain,'.inp')
      if (k.eq.0) then
         alog=trim(ain)//'.log'
         ain=trim(ain)//'.inp'
      else
         alog=ain(1:k)//'log'
      endif

!-----------------------------------------------------------------------
! Open the input and log files
!-----------------------------------------------------------------------
      open(iin,file=ain,form='formatted',status='old')
      open(ilog,file=alog,form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Read input file
!-----------------------------------------------------------------------
5     continue
      call rdinp(iin)

      i=0
      if (keyword(1).ne.'end-input') then

10       continue
         i=i+1
         
         if (keyword(i).eq.'basis') then
            lrungamess=.true.
            if (keyword(i+1).eq.'=') then
               i=i+2
               read(keyword(i),*) basname
            else
               goto 100
            endif
            
         else if (keyword(i).eq.'diffuse') then
            if (keyword(i+1).eq.'=') then
               i=i+2
               ! Diffuse function type
               if (keyword(i).eq.'kbj') then
                  difftype=1
               else
                  goto 100
               endif
               ! Numbers of diffuse functions
35             continue
               if (keyword(i+1).eq.',') then
                  i=i+2
                  call getdiffinfo(n,l,keyword(i))
                  ndiff(l)=n
                  goto 35
               endif
            else
               goto 100
            endif

         else if (keyword(i).eq.'geometry') then
            lrungamess=.true.
40          continue
            call rdinp(iin)
            if (keyword(1).ne.'end-geometry') then
               natm=natm+1
               goto 40
            endif
            ncoo=natm*3
            
         else
            ! Exit if the keyword is not recognised
            errmsg='Unknown keyword: '//trim(keyword(i))
            write(6,'(/,a,/)') trim(errmsg)
            STOP
         endif
         
         ! If there are more keywords to be read on the current line,
         ! then read them, else read the next line
         if (i.lt.inkw) then
            goto 10
         else
            goto 5
         endif
         
         ! Exit if a required argument has not been given with a keyword
100      continue
         errmsg='No argument given with the keyword '//trim(keyword(i))
         write(6,'(/,a,/)') trim(errmsg)
         STOP

      endif

!-----------------------------------------------------------------------
! Make sure that all the required information has been given
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Read the geometry section if necessary
!-----------------------------------------------------------------------
      if (lrungamess) call rdgeometry

!-----------------------------------------------------------------------
! Close the input file
!-----------------------------------------------------------------------
      close(iin)

      return

    end subroutine read_input

!#######################################################################

    subroutine getdiffinfo(n,l,string)
      
      implicit none
      
      integer          :: n,l,i,k
      character(len=*) :: string
      character(len=5) :: labels
      
      labels='spdfg'
      
      do i=1,len_trim(string)
         k=index(labels,string(i:i))
         if (k.ne.0) then
            read(string(1:i-1),*) n
            l=k
            exit
         endif
      enddo

      return

    end subroutine getdiffinfo 

!#######################################################################

    subroutine rdgeometry

      use relaxmod
      use parsemod
      use channels

      implicit none

      integer :: i,j

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      allocate(xcoo(ncoo))
      allocate(atlbl(natm))

!-----------------------------------------------------------------------
! Read the Cartesian coordinates from the input file
!-----------------------------------------------------------------------
      rewind(iin)

5     continue
      call rdinp(iin)
      if (keyword(1).ne.'geometry') goto 5

      do i=1,natm
         call rdinp(iin)
         atlbl(i)=keyword(1)
         do j=1,3
            read(keyword(j+1),*) xcoo(i*3-3+j)
         enddo
      enddo

      return 

    end subroutine rdgeometry

!#######################################################################

  end program relax2
