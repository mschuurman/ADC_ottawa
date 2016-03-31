  module prepmod

    use hessmod

    implicit none

    save

    integer            :: ityp=1
    character(len=120) :: filename

  contains

!#######################################################################
    
    subroutine hessprep
      
      use constants
      use hessmod

      implicit none

      integer           :: i,j,dir1,dir2,pm1,pm2
      character(len=60) :: acmmnd
      logical           :: exists

!-----------------------------------------------------------------------
! Create the output directory
!-----------------------------------------------------------------------
      inquire(file='adcinp/.',exist=exists)
      
      if (exists) then
         acmmnd='rm -rf adcinp/*'
      else
         acmmnd='mkdir adcinp'
      endif

      call system(acmmnd)

!-----------------------------------------------------------------------
! Reference geometry
!-----------------------------------------------------------------------
      ityp=0
      xcoo=xcoo0
      call wrfile(0,0,0,0)

!-----------------------------------------------------------------------
! Geometries required for the quadratic (on-diagonal) terms
!-----------------------------------------------------------------------
      ityp=1

      ! Loop over the plus/minus directions
      do dir1=1,2
         pm1=(-1)**dir1

         ! Loop over Cartesian coordinates
         do i=1,ncoo
            
            ! Set the current Cartesian coordinates
            xcoo=xcoo0
            xcoo(i)=xcoo(i)+pm1*dx
            
            ! Write the current input file
            call wrfile(i,0,pm1,0)

         enddo

      enddo

!-----------------------------------------------------------------------
! Geometries required for the bilinear (off-diagonal) terms
!-----------------------------------------------------------------------
      ityp=2
      
      ! Loop over the plus/minus directions
      do dir1=1,2
         pm1=(-1)**dir1         
         do dir2=1,2
            pm2=(-1)**dir2

            ! Loop over the unique pairs of Cartesian coordinates
            do i=1,ncoo
               do j=i+1,ncoo
                  
                  ! Set the current Cartesian coordinates
                  xcoo=xcoo0
                  xcoo(i)=xcoo(i)+pm1*dx
                  xcoo(j)=xcoo(j)+pm2*dx

                  ! Write the current input file
                  call wrfile(i,j,pm1,pm2)
                  
               enddo
            enddo

         enddo
      enddo
   
      return

    end subroutine hessprep

!#######################################################################

    subroutine wrfile(indx1,indx2,d1,d2)

      use iomod
      use hessmod

      implicit none
      
      integer          :: indx1,indx2,d1,d2,unit,i,j,k
      character(len=1) :: atmp1,atmp2

!-----------------------------------------------------------------------
! Write the filename
!-----------------------------------------------------------------------
      ! Reference geometry
      if (ityp.eq.0) then
         filename='adcinp/adc_x000.inp'
      endif

      ! Quadratic (on-diagonal) terms
      if (ityp.eq.1) then
         if (d1.eq.1) then
            atmp1='p'
         else if (d1.eq.-1) then
            atmp1='m'
         endif
         filename='adcinp/adc_x'
         k=len_trim(filename)
         if (indx1.lt.10) then
            write(filename(k+1:k+8),'(a2,i1,a5)') '00',indx1,atmp1//'.inp'
         else if (indx1.lt.100) then
            write(filename(k+1:k+8),'(a1,i2,a5)') '0',indx1,atmp1//'.inp'
         else
            write(filename(k+1:k+8),'(i3,a5)') indx1,atmp1//'.inp'
         endif         
      endif

      ! Bilinear (off-diagonal) terms
      if (ityp.eq.2) then
         if (d1.eq.1) then
            atmp1='p'
         else if (d1.eq.-1) then
            atmp1='m'
         endif
         if (d2.eq.1) then
            atmp2='p'
         else if (d2.eq.-1) then
            atmp2='m'
         endif
         filename='adcinp/adc_x'
         k=len_trim(filename)
         if (indx1.lt.10) then
            if (indx2.lt.10) then
               write(filename(k+1:k+13),'(a2,i1,a4,i1,a5)') &
                    '00',indx1,atmp1//'_00',indx2,atmp2//'.inp'
            else if (indx2.lt.100) then               
               write(filename(k+1:k+13),'(a2,i1,a3,i2,a5)') &
                    '00',indx1,atmp1//'_0',indx2,atmp2//'.inp'
            else               
               write(filename(k+1:k+13),'(a2,i1,a2,i3,a5)') &
                    '00',indx1,atmp1//'_',indx2,atmp2//'.inp'               
            endif
         else if (indx1.lt.100) then
            if (indx2.lt.10) then
               write(filename(k+1:k+13),'(a1,i2,a4,i1,a5)') &
                    '0',indx1,atmp1//'_00',indx2,atmp2//'.inp'
            else if (indx2.lt.100) then               
               write(filename(k+1:k+13),'(a1,i2,a3,i2,a5)') &
                    '0',indx1,atmp1//'_0',indx2,atmp2//'.inp'
            else               
               write(filename(k+1:k+13),'(a1,i2,a2,i3,a5)') &
                    '0',indx1,atmp1//'_',indx2,atmp2//'.inp'               
            endif
         else
            if (indx2.lt.10) then
               write(filename(k+1:k+13),'(i3,a4,i1,a5)') &
                    indx1,atmp1//'_00',indx2,atmp2//'.inp'
            else if (indx2.lt.100) then               
               write(filename(k+1:k+13),'(i3,a3,i2,a5)') &
                    indx1,atmp1//'_0',indx2,atmp2//'.inp'
            else               
               write(filename(k+1:k+13),'(i3,i2,a1,i3,a5)') &
                    indx1,atmp1//'_',indx2,atmp2//'.inp'               
            endif
         endif
      endif

!-----------------------------------------------------------------------
! Write the ADC input file
!-----------------------------------------------------------------------
      call freeunit(unit)

      open(unit,file=filename,form='formatted',status='new')

      do i=1,geomline
         write(unit,'(a)') trim(aline(i))
      enddo
      
      k=0
      do i=geomline+1,geomline+natm
         k=k+1
         write(unit,'(a,3(2x,F10.7))') aatm(k),(xcoo(j),j=k*3-2,k*3)
      enddo

      do i=geomline+natm+1,nlines
         write(unit,'(a)') trim(aline(i))
      enddo

      close(unit)

      return

    end subroutine wrfile

!#######################################################################

  end module prepmod
