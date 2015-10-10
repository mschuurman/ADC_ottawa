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

  end module dysonmod
