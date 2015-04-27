  module parsemod

    implicit none

    save

    integer                              :: iline,inkw
    integer, parameter                   :: maxkw=60
    integer, dimension(maxkw)            :: ilkw
    character(len=120), dimension(maxkw) :: keyword

  contains

!#######################################################################
!
! rdinp: reads the current line in the file opened using the 
!        io-unit 'unit' and splits this into keywords
! 
!        Keywords are separated by spaces and/or tabs (achar(9))
!
!        Nothing past a # character is read. That is, anything 
!        proceeding a # acts as a comment.
!
!        keyword: array of keywords present on the current line
!        nkw:     number of keyword on the current line
!        ilkw:    array of lengths of the keywords on the current line
!
!#######################################################################

    subroutine rdinp(unit)

      use iomod, only: errmsg,error_control

      implicit none

      integer            :: unit,i,k,istart,iend
      character(len=120) :: string,message

!------------------------------------------------------------------------
! Initialise arrays
!------------------------------------------------------------------------
      keyword=''
      ilkw=0

!------------------------------------------------------------------------
! Read the next line that is not blank and does not start with a comment
!------------------------------------------------------------------------
5     continue
      read(unit,'(a)') string
    
      ! Skip blank lines
      if (string.eq.'') goto 5
    
      ! Read to first keyword
      iend=len_trim(string)
      do i=1,iend
         if (string(i:i).ne.''.and.string(i:i).ne.achar(9)) then
            istart=i
            exit
         endif
      enddo

      ! Skip lines beginning with a comment
      if (string(i:i).eq.'#') goto 5

      ! If line is not blank and does not start with a comment, then
      ! extract the keywords
      inkw=0
      k=0
      do i=istart,iend

         ! Terminate program with an error message if the number of
         ! keywords exceeds the maximum value
         if (inkw.eq.maxkw) goto 100

         ! Break if we have reached a comment
         if (string(i:i).eq.'#') exit

         ! If the current character is a delimiter other than a space
         ! (=, ',', ], etc.) then read in as a separate keyword and
         ! then go to the next keyword...
         if (string(i:i).eq.'='&
              .or.string(i:i).eq.','&
              .or.string(i:i).eq.'('&
              .or.string(i:i).eq.')'&
              .or.string(i:i).eq.'['&
              .or.string(i:i).eq.']'&
              .or.string(i:i).eq.'{'&
              .or.string(i:i).eq.'}') then
            inkw=inkw+1
            read(string(i:i),'(a1)') keyword(inkw)
            k=0

         ! ... otherwise keep adding to the current keyword...
         else if (string(i:i).ne.' '.and.string(i:i).ne.achar(9)) then
            k=k+1
            if (k.eq.1) inkw=inkw+1
            read(string(i:i),'(a1)') keyword(inkw)(k:k)

         ! ... until a blank space or tab is reached, at which point
         ! we move to the next keyword
         else
            k=0
         endif

      enddo

      ! Convert all keywords to lowercase
      do i=1,inkw
         call lowercase(keyword(i))
      enddo

      ! Determine keyword lengths
      do i=1,inkw
         ilkw(i)=len_trim(keyword(i))
      enddo

      return

100   continue
      write(string(1:),'(i4)') iline
      k=len_trim(string)
      errmsg='Number of keywords on line'//string(1:k)//' exceeds maxkw'
      call error_control

    end subroutine rdinp

!#######################################################################
! 
!   lowercase: replaces all uppercase letters in a given character
!              string with the corresponding lowercase letters
!
!#######################################################################
    subroutine lowercase(string)
      
      implicit none
      
      integer*8    ::  i,j,length
      character(*) :: string
      
      length=len(string)
      
      do i=1,length
         do j=65,90
            if (ichar(string(i:i)).eq.j) string(i:i)=char(j+32)
         enddo
      enddo

      return

    end subroutine lowercase

!#######################################################################

  end module parsemod
