  module mp2

    contains

!#######################################################################

      subroutine mp2_master(e0)

        use channels
        use constants
        use parameters

        implicit none

        integer           :: i
        real(d)           :: e0,d2
        character(len=60) :: atmp

!-----------------------------------------------------------------------
! Calculation of the MP2 correlation energy
!-----------------------------------------------------------------------
        call mp2_ener

!-----------------------------------------------------------------------
! Calculation of the D2 diagnostic
!-----------------------------------------------------------------------
        call mp2_d2(d2)

!-----------------------------------------------------------------------
! Output results
!-----------------------------------------------------------------------
        e0=Ehf+E_MP2
 
        write(ilog,'(/,2x,90a)') ('*',i=1,90)

        write(ilog,'(2x,1a)') '*'
        atmp='* HF energy:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),Ehf

        write(ilog,'(2x,1a)') '*'
        atmp='* MP2 correlation energy:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),E_MP2

        write(ilog,'(2x,1a)') '*'
        atmp='* MP2 energy:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),e0

        write(ilog,'(2x,1a)') '*'
        atmp='* D2 diagnostic:'
        write(ilog,'(2x,a25,2x,F18.10)') adjustl(atmp),d2
        
        write(ilog,'(2x,1a)') '*'
        write(ilog,'(2x,90a,/)') ('*',i=1,90)

        return

      end subroutine mp2_master

!#######################################################################
! mp2_ener: Calculates the MP2 correlation energy
!#######################################################################

      subroutine mp2_ener

        use channels
        use constants
        use parameters
        use vpqrsmod

        implicit none

        integer :: r,s,u,v,nsym1,nsym2,i,j,a,b,cnt
        real(d) :: DA,eijc,term,etotal

        E_MP2 = 0._d
        etotal = 0._d
        cnt=0
        do r=nOcc+1,nBas
           a=roccnum(r)  !r
           
           do s=nOcc+1,nBas
              b=roccnum(s) !s
          
              do u=1,nOcc
                 i=roccnum(u) !u
                 
                 do v=1,nOcc
                    j=roccnum(v) !v
                    
                    nsym1=MT(orbSym(i),orbSym(a))
                    nsym2=MT(orbSym(j),orbSym(b))


                    if  (MT(nsym1,nsym2) .eq. 1)  then

                       cnt=cnt+1

                       eijc=e(i)+e(j)-e(a)-e(b)
                       
                       term= vpqrs(i,a,j,b)*(2._d*vpqrs(i,a,j,b)-vpqrs(i,b,j,a))
                       
                       term=term/eijc

                       E_MP2 = E_MP2 + term

                    else

                       eijc=e(i)+e(j)-e(a)-e(b)
                       term= vpqrs(i,a,j,b)*(2._d*vpqrs(i,a,j,b)-vpqrs(i,b,j,a)) / eijc
                       etotal = etotal + term
                       if(abs(term).gt.1e-5) write(ilog,100) &
                            labsym(orbsym(i)),i,labsym(orbsym(j)),j,&
                            labsym(orbsym(a)),a,labsym(orbsym(b)),b,&
                            MT(nsym1,nsym2),term

                       
                    endif
                    
                 end do
              end do
           end do
        end do
        
        if(Etotal.gt.1e-5) write(ilog,*) 'E[sym residual]=',Etotal
        
100     format(a3,'(',i3,') ',a3,'(',i3,') ',a3,'(',i3,') ',a3,'(',i3,') sym=',i3,'term=',f15.8)
        
        return
        
      end subroutine mp2_ener

!#######################################################################
! mp2_d2: Calculates the D2(MP1) MP2 diagnostic of Nielsen and
!         Janssen. See Chem. Phys. Lett., 310, 568 (1999).
!#######################################################################

    subroutine mp2_d2(d2)

      use channels
      use constants
      use parameters
      use iomod

      implicit none

      integer                              :: i,j,k,a,b,c,albl,blbl,&
                                              workdim,error
      real(d)                              :: d2
      real(d), dimension(:,:), allocatable :: toto,tvtv
      real(d), dimension(:), allocatable   :: eigval,workarr

!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      nvirt=nbas-nocc
      allocate(toto(nocc,nocc))
      allocate(tvtv(nvirt,nvirt))

!-----------------------------------------------------------------------
! Calculation of T^o * T^o
!-----------------------------------------------------------------------
      toto=0.0d0
      do i=1,nocc
         do j=1,nocc

            do k=1,nocc
               do a=nocc+1,nbas
                  do b=nocc+1,nbas
                     toto(i,j)=toto(i,j)+t2amp(i,k,a,b)*t2amp(j,k,a,b)
                  enddo
               enddo
            enddo

         enddo
      enddo

!-----------------------------------------------------------------------
! Calculation of T^v * T^v
!-----------------------------------------------------------------------
      tvtv=0.0d0
      albl=0
      do a=nocc+1,nbas
         albl=albl+1
        
         blbl=0
         do b=nocc+1,nbas
            blbl=blbl+1
            
            do i=1,nocc
               do j=1,nocc
                  do c=nocc+1,nbas
                     tvtv(albl,blbl)=tvtv(albl,blbl)+&
                          t2amp(i,j,a,c)*t2amp(i,j,b,c)
                  enddo
               enddo
            enddo

         enddo
      enddo

!-----------------------------------------------------------------------
! Eigenvalues of T^o * T^o
!-----------------------------------------------------------------------
      workdim=3*nocc
      allocate(workarr(workdim))
      allocate(eigval(nocc))

      call dsyev('N','U',nocc,toto,nocc,eigval,workarr,workdim,error)

      if (error.ne.0) then
         errmsg='Diagonalisation of T^o * T^o failed in subroutine &
              mp2_d2'
         call error_control
      endif      

      d2=sqrt(eigval(nocc))

      deallocate(workarr)
      deallocate(eigval)

!-----------------------------------------------------------------------
! Eigenvalues of T^v * T^v
!-----------------------------------------------------------------------
      workdim=3*nvirt
      allocate(workarr(workdim))
      allocate(eigval(nvirt))

      call dsyev('N','U',nvirt,tvtv,nvirt,eigval,workarr,workdim,error)

      if (error.ne.0) then
         errmsg='Diagonalisation of T^v * T^v failed in subroutine &
              mp2_d2'
         call error_control
      endif

      if (sqrt(eigval(nvirt)).gt.d2) d2=sqrt(eigval(nvirt))

      deallocate(workarr)
      deallocate(eigval)

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(toto)
      deallocate(tvtv)

      return

    end subroutine mp2_d2

!#######################################################################

    function t2amp(i,j,a,b)

      use constants
      use parameters
      use vpqrsmod

      implicit none

      integer, intent(in) :: i,j,a,b
      real(d)             :: t2amp

      t2amp=vpqrs(i,a,j,b)/(e(i)+e(j)-e(a)-e(b))

      return

    end function t2amp

!#######################################################################

  end module mp2
