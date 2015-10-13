  module dysonmod

    use constants
    use parameters
    use vpqrsmod

  contains

!#######################################################################
! rhogs2_uu: returns an element of the unoccupied-unoccupied block of
!            the 2nd-order correction to the ground state density
!            matrix
!#######################################################################
   
    function rhogs2_uu(a,b) result(fret)

      implicit none

      integer :: a,b,i1,j1,c1,i,j,c
      real(d) :: fret,ftmp,delta_ijac,delta_ijbc
      
      fret=0.0d0

      do i1=1,nocc
         i=roccnum(i1)
         do j1=1,nocc
            j=roccnum(j1)
            do c1=nocc+1,nbas
               c=roccnum(c1)

               delta_ijac=1.0d0/(e(a)+e(c)-e(i)-e(j))

               delta_ijbc=1.0d0/(e(b)+e(c)-e(i)-e(j))

               ftmp=vpqrs(a,i,c,j)*(2.0d0*vpqrs(i,b,j,c)-vpqrs(i,c,j,b))&
                    +vpqrs(c,i,a,j)*(2.0d0*vpqrs(i,c,j,b)-vpqrs(i,b,j,c))

               fret=fret+delta_ijac*delta_ijbc*ftmp

            enddo
         enddo
      enddo

      fret=0.5d0*fret

      return

    end function rhogs2_uu

!#######################################################################
! rhogs2_oo: returns an element of the occupied-occupied block of
!            the 2nd-order correction to the ground state density
!            matrix
!#######################################################################

    function rhogs2_oo(i,j) result(fret)

      implicit none

      integer :: i,j,k1,a1,b1,k,a,b
      real(d) :: fret,ftmp,delta_ikab,delta_jkab

      fret=0.0d0

      do k1=1,nocc
         k=roccnum(k1)
         do a1=nocc+1,nbas
            a=roccnum(a1)
            do b1=nocc+1,nbas
               b=roccnum(b1)
               
               delta_ikab=1.0d0/(e(a)+e(b)-e(i)-e(k))
               delta_jkab=1.0d0/(e(a)+e(b)-e(j)-e(k))

               ftmp=vpqrs(a,i,b,k)*(2.0d0*vpqrs(j,a,k,b)-vpqrs(j,b,k,a))&
                    +vpqrs(b,i,a,k)*(2.0d0*vpqrs(j,b,k,a)-vpqrs(j,a,k,b))

               fret=fret+delta_ikab*delta_jkab*ftmp

            enddo
         enddo
      enddo

      fret=-0.5d0*fret

      return

    end function rhogs2_oo

!#######################################################################
! rhogs2_ou: returns an element of the occupied-unoccupied block of
!            the 2nd-order correction to the ground state density
!            matrix
!#######################################################################

    function rhogs2_ou(i,a) result(fret)

      implicit none

      integer :: i,a,j1,k1,b1,c1,j,k,b,c
      real(d) :: fret,ftmp,term1,term2,delta_ijbc,delta_jkab

      term1=0.0d0
      do j1=1,nocc
         j=roccnum(j1)
         do b1=nocc+1,nbas
            b=roccnum(b1)
            do c1=nocc+1,nbas
               c=roccnum(c1)
               
               delta_ijbc=1.0d0/(e(b)+e(c)-e(i)-e(j))
      
               ftmp=vpqrs(i,b,j,c)*(vpqrs(j,b,a,c)-2.0d0*vpqrs(j,c,a,b))&
                    +vpqrs(i,c,j,b)*(vpqrs(j,c,a,b)-2.0d0*vpqrs(j,b,a,c))

               term1=term1+delta_ijbc*ftmp

            enddo
         enddo
      enddo

      term2=0.0d0
      do j1=1,nocc
         j=roccnum(j1)
         do k1=1,nocc
            k=roccnum(k1)
            do b1=nocc+1,nbas
               b=roccnum(b1)
               
               delta_jkab=1.0d0/(e(a)+e(b)-e(j)-e(k))

               ftmp=vpqrs(j,i,k,b)*(2.0d0*vpqrs(j,a,k,b)-vpqrs(j,b,k,a))&
                    +vpqrs(j,b,k,i)*(2.0d0*vpqrs(j,b,k,a)-vpqrs(j,a,k,b))

               term2=term2+delta_jkab*ftmp

            enddo
         enddo
      enddo

      fret=-0.5d0*(term1+term2)/(e(a)-e(i))

    end function rhogs2_ou
    
!#######################################################################
! dyscoeff_gs_occ: calculates the coefficients for the occupied MOs in
!                  the expansion of the Dyson orbital corresponding to
!                  ionisation from the ground state to a single
!                  IP-ADC state
!
! rhogs2   - 2nd-order correction to the ground state density matrix
! dyscoeff - Dyson orbital coefficients in the MO basis
! vec      - IP-ADC state vector of interest
!
!#######################################################################
    
    subroutine dyscoeff_gs_occ(rhogs2,dyscoeff,kpqf,vec,ndimf,ndimsf)
        
      use constants        
      use parameters
      use vpqrsmod

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
      integer                                   :: ndimf,ndimsf,&
                                                   alpha,i,j,ilbl,&
                                                   jlbl,b,n,k,c
      real(d), dimension(nbas,nbas)             :: rhogs2
      real(d), dimension(nbas)                  :: dyscoeff
      real(d), dimension(ndimf)                 :: vec
      real(d), dimension(:,:), allocatable      :: chi,zeta
      real(d)                                   :: delta_ijab,ftmp
        
!-----------------------------------------------------------------------
! Allocate arrays
!-----------------------------------------------------------------------
      ! Note that these arrays only need to be of
      ! dimension (nocc,nvirt), but for simplicity they are of
      ! dimension (nbas,nbas), as this requires only a negligible
      ! extra amount of memory 
      allocate(chi(nbas,nbas))
      allocate(zeta(nbas,nbas))
        
!-----------------------------------------------------------------------
! Index of the 'continuum' orbital
!-----------------------------------------------------------------------
      alpha=ifakeorb

!-----------------------------------------------------------------------
! Zeroth-order contribution
!-----------------------------------------------------------------------
      do i=1,ndimsf
         ilbl=kpqf(3,i)
         dyscoeff(ilbl)=vec(i)
      enddo

!-----------------------------------------------------------------------
! Second-oder contribution
! Note that the first-order contribution vanishes
!-----------------------------------------------------------------------
      ! Term C
      ! (i) Pre-calculate the chi and zeta terms
      chi=0.0d0
      zeta=0.0d0
      do j=1,nocc
         do b=nocc+1,nbas
            do n=1,ndimsf
               k=kpqf(3,n)
               c=kpqf(5,n)
               chi(j,b)=chi(j,b)+vec(n)*(vpqrs(b,j,c,k)-vpqrs(c,j,b,k))
               zeta(j,b)=zeta(j,b)+vec(n)*(vpqrs(c,j,b,k)-vpqrs(b,j,c,k))
            enddo
         enddo
      enddo
      ! (ii) Contribution of the C terms
      do i=1,nocc
         do j=1,nocc
            do b=nocc+1,nbas
               delta_ijab=1.0d0/(e(alpha)+e(b)-e(j)-e(j))
               ftmp=vpqrs(i,alpha,j,b)*chi(j,b)+vpqrs(i,b,j,alpha)*zeta(j,b)
               ftmp=delta_ijab*ftmp
               dyscoeff(i)=dyscoeff(i)+ftmp
            enddo
         enddo
      enddo

      ! Term D
      do i=1,ndimsf
         ilbl=kpqf(3,i)
         dyscoeff(ilbl)=dyscoeff(ilbl)-0.5d0*rhogs2(alpha,alpha)*vec(i)
      enddo

      ! Term E
      do i=1,nocc
         do j=1,ndimsf
            jlbl=kpqf(3,j)
            dyscoeff(i)=dyscoeff(i)+0.5d0*rhogs2(i,jlbl)*vec(j)
         enddo
      enddo

!-----------------------------------------------------------------------
! Deallocate arrays
!-----------------------------------------------------------------------
      deallocate(chi)
      deallocate(zeta)
      
      return

    end subroutine dyscoeff_gs_occ

!#######################################################################
! dyscoeff_gs_unocc: calculates the coefficients for the unoccupied MOs
!                    in the expansion of the Dyson orbital corresponding
!                    to ionisation from the ground state to a single
!                    IP-ADC state
!
! rhogs2   - 2nd-order correction to the ground state density matrix
! dyscoeff - Dyson orbital coefficients in the MO basis
! vec      - IP-ADC state vector of interest
!
!#######################################################################

    subroutine dyscoeff_gs_unocc(rhogs2,dyscoeff,kpqf,vec,ndimf,ndimsf)

      use constants
      use parameters
      use vpqrsmod
      
      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
      integer                                   :: ndimf,ndimsf,&
                                                   alpha,b,i,ilbl,&
                                                   j,c,n
      real(d), dimension(nbas,nbas)             :: rhogs2
      real(d), dimension(nbas)                  :: dyscoeff
      real(d), dimension(ndimf)                 :: vec
      real(d)                                   :: delta_ijbc,ftmp
      
!-----------------------------------------------------------------------
! Index of the 'continuum' orbital
!-----------------------------------------------------------------------
      alpha=ifakeorb
        
!-----------------------------------------------------------------------
! Second-order contributions
! Note that the zeroth and first-order contributions vanish
!-----------------------------------------------------------------------
      ! Term F
      do b=nocc+1,nbas
         do i=1,nocc
            ilbl=kpqf(3,i)
            dyscoeff(b)=dyscoeff(b)+vec(i)*rhogs2(ilbl,b)
         enddo
      enddo

      ! Term G
      do b=nocc+1,nbas
         do n=ndimsf+1,ndimf
            i=kpqf(3,n)
            j=kpqf(4,n)
            if (kpqf(5,n).eq.alpha) then
               c=kpqf(6,n)
            else
               c=kpqf(5,n)
            endif
            delta_ijbc=1.0d0/(e(b)+e(c)-e(i)-e(j))
            ftmp=vec(n)*(2.0d0*vpqrs(b,i,c,j)-vpqrs(c,i,b,j))
            ftmp=delta_ijbc*ftmp
            dyscoeff(b)=dyscoeff(b)+ftmp
         enddo
      enddo
        
      return
        
    end subroutine dyscoeff_gs_unocc
      
!#######################################################################

  end module dysonmod
