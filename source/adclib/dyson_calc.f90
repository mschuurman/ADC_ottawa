!#######################################################################
! dyson_calc: subroutines and functions required for the calculation of
!             Dyson orbital expansion coefficients in the MO basis.
!
!             Equations adapted from spin-orbtal expressions for the 
!             density and transition density matrices given in 
!             Dreuw et al., Mol. Phys, 112, 774 (2014).
!#######################################################################

  module dyson_calc

    use constants
    use parameters
    use vpqrsmod

    save

    real(d), parameter :: invsqrt2=1.0d0/sqrt(2.0d0)

  contains

!#######################################################################
! dyson_precalc: pre-calculates intermediate quantities required in
!                the calculation of the ground-to-excited state and
!                state-to-state Dyson orbital coefficients
!#######################################################################

    subroutine dyson_precalc(rhogs2,rmat,smat,vec_init,ndim,ndims,kpq)

      use mp2
      
      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,ndims,i,j,k,a,&
                                                   b,c,n
      real(d), dimension(nbas,nbas)             :: rhogs2,rmat,smat
      real(d), dimension(ndim)                  :: vec_init
      real(d)                                   :: delta_ijab,ftmp,&
                                                   delta_ikab,&
                                                   delta_ijac,prei,&
                                                   pref

!-----------------------------------------------------------------------
! Pre-calculation of the 2nd-order correction to the ground state
! density matrix
!
! Note that this matrix is Hermitian, and so we only need to calculate
! the lower triangle
!-----------------------------------------------------------------------
        ! Occupied-occupied block
        do i=1,nocc
           do j=i,nocc
              rhogs2(i,j)=rhogs2_oo(i,j)
              rhogs2(j,i)=rhogs2(i,j)
           enddo
        enddo

        ! Unoccupied-unoccupied block
        do a=nocc+1,nbas
           do b=a,nbas
              rhogs2(a,b)=rhogs2_uu(a,b)
              rhogs2(b,a)=rhogs2(a,b)
           enddo
        enddo

        ! Occupied-unoccupied block
        do i=1,nocc
           do a=nocc+1,nbas
              rhogs2(i,a)=rhogs2_ou(i,a)
              rhogs2(a,i)=rhogs2(i,a)
           enddo
        enddo

!-----------------------------------------------------------------------
! Intermediate terms required for the calculation of the state-to-state
! transition density matrix elements
!-----------------------------------------------------------------------
        if (statenumber.gt.0) then

           rmat=0.0d0
           smat=0.0d0

           ! R-matrix
           do i=1,nocc
              do a=nocc+1,nbas
                 do n=1,ndims
                    j=kpq(3,n)
                    b=kpq(5,n)
                    delta_ijab=1.0d0/(e(a)+e(b)-e(i)-e(j))
                    ftmp=2.0d0*vpqrs(a,i,b,j)-vpqrs(b,i,a,j)
                    ftmp=delta_ijab*0.5d0*vec_init(n)*ftmp
                    rmat(i,a)=rmat(i,a)+ftmp
                 enddo
              enddo
           enddo

           ! S-matrix, occupied-occupied block
           call precalc_smat(smat,vec_init,kpq,ndim,ndims)

        endif

      return

    end subroutine dyson_precalc

!#######################################################################

    subroutine precalc_smat(smat,vec_init,kpq,ndim,ndims)

      use constants
      use parameters

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq
      integer                                   :: ndim,ndims,count,n,&
                                                   i,j,k,a,b
      real(d), dimension(nbas,nbas)             :: smat
      real(d), dimension(ndim)                  :: vec_init
      real(d)                                   :: ftmp,delta

!-----------------------------------------------------------------------
! Initialise arrays
!-----------------------------------------------------------------------
      smat=0.0d0

!-----------------------------------------------------------------------
! j=k, a|=b
!-----------------------------------------------------------------------
      count=ndims+kpq(2,0)
      do n=count+1,count+kpq(3,0)
         j=kpq(3,n)
         k=kpq(4,n)
         a=kpq(5,n)
         b=kpq(6,n)
         do i=1,nocc
            delta=1.0d0/(e(a)+e(b)-e(i)-e(j))
            ftmp=vpqrs(a,i,b,j)-vpqrs(b,i,a,j)
            smat(i,j)=smat(i,j)+delta*ftmp*invsqrt2*vec_init(n)
         enddo

      enddo

!-----------------------------------------------------------------------
! j|=k, a|=b, spin case I
!-----------------------------------------------------------------------
      count=ndims+kpq(2,0)+kpq(3,0)+kpq(4,0)
      do n=count+1,count+kpq(5,0)
         j=kpq(3,n)
         k=kpq(4,n)
         a=kpq(5,n)
         b=kpq(6,n)
         do i=1,nocc
            delta=1.0d0/(e(a)+e(b)-e(i)-e(k))
            ftmp=vpqrs(a,i,b,k)-vpqrs(b,i,a,k)
            smat(i,j)=smat(i,j)+delta*ftmp*0.5d0*vec_init(n)
         enddo
      enddo

!-----------------------------------------------------------------------
! j|=k, a|=b, spin case II
!-----------------------------------------------------------------------
      count=ndims+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
      do n=count+1,count+kpq(5,0)
         j=kpq(3,n)
         k=kpq(4,n)
         a=kpq(5,n)
         b=kpq(6,n)
         do i=1,nocc
            delta=1.0d0/(e(a)+e(b)-e(i)-e(k))
            ftmp=vpqrs(a,i,b,k)-vpqrs(b,i,a,k)
            smat(i,j)=smat(i,j)&
                 +delta*ftmp*(3.0d0/sqrt(12.0d0))*vec_init(n)
         enddo
      enddo

      return

    end subroutine precalc_smat

!#######################################################################
! dyscoeff_gs: calculates the coefficients in the MO expansion of a
!              single dyson orbital corresponding to ionisation from
!              the ground state
!#######################################################################

    subroutine dyscoeff_gs(rhogs2,dyscoeff,kpqf,vec,ndimf,ndimsf)

      use constants        
      use parameters
      
      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpqf
      integer                                   :: ndimf,ndimsf
      real(d), dimension(nbas,nbas)             :: rhogs2
      real(d), dimension(nbas)                  :: dyscoeff
      real(d), dimension(ndimf)                 :: vec

!-----------------------------------------------------------------------
! Coefficients for the occupied MOs
!-----------------------------------------------------------------------
      call dyscoeff_gs_occ(rhogs2,dyscoeff,kpqf,vec,ndimf,ndimsf)

!-----------------------------------------------------------------------
! Coefficients for the unoccupied MOs
!-----------------------------------------------------------------------
      call dyscoeff_gs_unocc(rhogs2,dyscoeff,kpqf,vec,ndimf,ndimsf)

!-----------------------------------------------------------------------
! Account for the other spin component
!-----------------------------------------------------------------------
      dyscoeff=sqrt(2.0d0)*dyscoeff

      return

    end subroutine dyscoeff_gs

!#######################################################################
! dyscoeff_exci: calculates the coefficients in the MO expansion of a
!                single Dyson orbital corresponding to ionisation from
!                an excited state
!#######################################################################

    subroutine dyscoeff_exci(rhogs2,dyscoeff,kpq,kpqf,vec_init,&
         vec_final,ndim,ndims,ndimf,ndimsf,rmat,smat)

      use constants
      use parameters
      
      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
      integer                                   :: ndim,ndims,ndimf,&
                                                   ndimsf
      real(d), dimension(nbas,nbas)             :: rhogs2,rmat,smat,&
                                                   pmat,qmat
      real(d), dimension(nbas)                  :: dyscoeff
      real(d), dimension(ndim)                  :: vec_init
      real(d), dimension(ndimf)                 :: vec_final

!-----------------------------------------------------------------------
! Pre-calculation of the P-matrix
! Note that only the occupied-occupied and occupied-unoccupied blocks
! vanish
!
! Note that the current implementation is very, very stupid
!-----------------------------------------------------------------------
      call precalc_pmat(ndim,ndimf,ndims,ndimsf,pmat,vec_init,&
           vec_final,kpq,kpqf)

!-----------------------------------------------------------------------
! Pre-calculation of the Q-matrix
! Not that only the unoccupied-unoccupied block does not vanish
!
! Note that, again, the current implementation is very, very stupid
!-----------------------------------------------------------------------
      call precalc_qmat(ndim,ndimf,ndims,ndimsf,qmat,vec_init,&
           vec_final,kpq,kpqf)

!-----------------------------------------------------------------------
! Coefficients for the occupied MOs
!-----------------------------------------------------------------------
      call dyscoeff_exci_occ(rhogs2,dyscoeff,kpq,kpqf,vec_init,vec_final,&
           ndim,ndims,ndimf,ndimsf,smat,pmat)

!-----------------------------------------------------------------------      
! Coefficients for the unoccupied MOs
!-----------------------------------------------------------------------
      call dyscoeff_exci_unocc(rhogs2,dyscoeff,kpq,kpqf,vec_init,vec_final,&
           ndim,ndims,ndimf,ndimsf,rmat,pmat,qmat)

!-----------------------------------------------------------------------
! Account for the other spin component
!-----------------------------------------------------------------------
      dyscoeff=sqrt(2.0d0)*dyscoeff

      return

    end subroutine dyscoeff_exci

!#######################################################################

    subroutine precalc_pmat(ndim,ndimf,ndims,ndimsf,pmat,vec_init,&
         vec_final,kpq,kpqf)

      use constants
      use parameters

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
      integer                                   :: ndim,ndimf,ndims,&
                                                   ndimsf,m,n,a,b,&
                                                   count,i,j
      real(d), dimension(nbas,nbas)             :: pmat
      real(d), dimension(ndim)                  :: vec_init
      real(d), dimension(ndimf)                 :: vec_final

!-----------------------------------------------------------------------
! Initialisation of arrays
!-----------------------------------------------------------------------
      pmat=0.0d0

!-----------------------------------------------------------------------
! Unoccupied-unoccupied block
!-----------------------------------------------------------------------
      do m=1,ndimsf
         do n=1,ndims
            if (kpqf(3,m).eq.kpq(3,n)) then
               a=kpqf(5,m)
               b=kpq(5,n)
               pmat(a,b)=pmat(a,b)+invsqrt2*vec_final(m)*invsqrt2*vec_init(n)
            endif
         enddo
      enddo

!-----------------------------------------------------------------------
! Unoccupied-occupied block
!-----------------------------------------------------------------------
      do n=1,ndims
         j=kpq(3,n)
         b=kpq(5,n)
         
         ! i=j, a|=b
         count=kpqf(1,0)+kpqf(2,0)
         do m=count+1,count+kpqf(3,0)
            if (j.eq.kpqf(4,m).and.b.eq.kpqf(6,m)) then
               i=kpqf(3,m)
               a=kpqf(5,m)
               pmat(a,i)=pmat(a,i)-vec_final(m)*vec_init(n)
            endif
         enddo

         ! i|=j, a|=b, spin case I
         count=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
         do m=count+1,count+kpqf(5,0)
            if (j.eq.kpqf(4,m).and.b.eq.kpqf(6,m)) then
               i=kpqf(3,m)
               a=kpqf(5,m)
               pmat(a,i)=pmat(a,i)-invsqrt2*vec_final(m)*vec_init(n)
            endif
         enddo

         ! i|=j, a|=b, spin case II
         count=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
         do m=count+1,count+kpqf(5,0)
            if (j.eq.kpqf(4,m).and.b.eq.kpqf(6,m)) then
               i=kpqf(3,m)
               a=kpqf(5,m)
               pmat(a,i)=pmat(a,i)-(6.0d0/sqrt(24.0d0))*vec_final(m)*vec_init(n)
            endif
         enddo

      enddo

!-----------------------------------------------------------------------
! CHECK: should the prefactor be -1/2 or -1 instead of -2?
!-----------------------------------------------------------------------
!      pmat=pmat/4.0d0
      pmat=pmat/2.0d0

      return

    end subroutine precalc_pmat

!#######################################################################

    subroutine precalc_qmat(ndim,ndimf,ndims,ndimsf,qmat,vec_init,&
         vec_final,kpq,kpqf)

      use constants
      use parameters

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
      integer                                   :: ndim,ndimf,ndims,&
                                                   ndimsf,m,n,a,b,&
                                                   count,countf,i,j,&
                                                   alpha
      real(d), dimension(nbas,nbas)             :: qmat
      real(d), dimension(ndim)                  :: vec_init
      real(d), dimension(ndimf)                 :: vec_final
      logical                                   :: lcontrib

!-----------------------------------------------------------------------
! Initialise arrays
!-----------------------------------------------------------------------
      qmat=0.0d0

!-----------------------------------------------------------------------
! Index of the 'continuum' orbital
!-----------------------------------------------------------------------
      alpha=ifakeorb

!-----------------------------------------------------------------------
! Final:   i=j, a|=c
! Initial: i=j, b=c
!-----------------------------------------------------------------------
      count=kpq(1,0)
      countf=kpqf(1,0)+kpqf(2,0)
      do m=countf+1,countf+kpqf(3,0)
         do n=count+1,count+kpq(2,0)
            call qcontrib(lcontrib,alpha,kpq,kpqf,m,n,a,b)
            if (.not.lcontrib) cycle
            qmat(a,b)=qmat(a,b)+sqrt(2.0d0)*vec_final(m)*vec_init(n)
         enddo
      enddo

!-----------------------------------------------------------------------
! Final:   i=j, a|=c
! Initial: i=j, b|=c
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)
      countf=kpqf(1,0)+kpqf(2,0)
      do m=countf+1,countf+kpqf(3,0)
         do n=count+1,count+kpq(3,0)
            call qcontrib(lcontrib,alpha,kpq,kpqf,m,n,a,b)
            if (.not.lcontrib) cycle
            qmat(a,b)=qmat(a,b)+vec_final(m)*vec_init(n)
         enddo
      enddo

!-----------------------------------------------------------------------
! Final:   i|=j, a|=c, spin case II
! Initial: i|=j, b=c
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)+kpq(3,0)
      countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
      do m=countf+1,countf+kpqf(5,0)
         do n=count+1,count+kpq(4,0)
            call qcontrib(lcontrib,alpha,kpq,kpqf,m,n,a,b)
            if (.not.lcontrib) cycle
            qmat(a,b)=qmat(a,b)+ &
                 (1.0d0/sqrt(6.0d0))*vec_final(m)*vec_init(n)
         enddo
      enddo

!-----------------------------------------------------------------------
! Final:   i|=j, a|=c, spin case I
! Initial: i|=j, b|=c, spin case I
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)
      countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)
      do m=countf+1,countf+kpqf(5,0)
         do n=count+1,count+kpq(5,0)
            call qcontrib(lcontrib,alpha,kpq,kpqf,m,n,a,b)
            if (.not.lcontrib) cycle
            qmat(a,b)=qmat(a,b)+vec_final(m)*vec_init(n)
         enddo
      enddo

!-----------------------------------------------------------------------
! Final:   i|=j, a|=c, spin case II
! Initial: i|=j, b|=c, spin case II
!-----------------------------------------------------------------------
      count=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+kpq(5,0)
      countf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+kpqf(5,0)
      do m=countf+1,countf+kpqf(5,0)
         do n=count+1,count+kpq(5,0)
            call qcontrib(lcontrib,alpha,kpq,kpqf,m,n,a,b)
            if (.not.lcontrib) cycle
            qmat(a,b)=qmat(a,b)+ &
                 ((1.0d0/6.0d0)+(4.0d0/sqrt(24.0d0)))*vec_final(m)*vec_init(n)
         enddo
      enddo

      return
      
    end subroutine precalc_qmat

!#######################################################################

    subroutine qcontrib(lcontrib,alpha,kpq,kpqf,m,n,a,b)

      use constants
      use parameters

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
      integer                                   :: alpha,m,n,a,b
      logical                                   :: lcontrib

      lcontrib=.false.

      if (kpqf(5,m).eq.alpha) then
         if (kpqf(3,m).eq.kpq(3,n).and.kpqf(4,m).eq.kpq(4,n) &
              .and.kpqf(6,m).eq.kpq(6,n)) then
            lcontrib=.true.
            a=kpqf(5,m)
            b=kpq(5,n)
         endif
      else if (kpqf(6,m).eq.alpha) then
         if (kpqf(3,m).eq.kpq(3,n).and.kpqf(4,m).eq.kpq(4,n) &
              .and.kpqf(5,m).eq.kpq(5,n)) then
            lcontrib=.true.
            a=kpqf(6,m)
            b=kpq(6,n)
         endif
      endif

      return

    end subroutine qcontrib

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
         dyscoeff(ilbl)=invsqrt2*vec(i)
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
               chi(j,b)=chi(j,b)+invsqrt2*vec(n)*(vpqrs(b,j,c,k)-vpqrs(c,j,b,k))
               zeta(j,b)=zeta(j,b)+invsqrt2*vec(n)*(vpqrs(c,j,b,k)-vpqrs(b,j,c,k))
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
         dyscoeff(ilbl)=dyscoeff(ilbl)-0.25d0*rhogs2(alpha,alpha)*invsqrt2*vec(i)
      enddo

      ! Term E
      do i=1,nocc
         do j=1,ndimsf
            jlbl=kpqf(3,j)
            dyscoeff(i)=dyscoeff(i)+0.25d0*rhogs2(i,jlbl)*invsqrt2*vec(j)
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
         do n=1,ndimsf 
            i=kpqf(3,n)
            dyscoeff(b)=dyscoeff(b)+0.5d0*invsqrt2*vec(n)*rhogs2(i,b)
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
            ftmp=0.25d0*vec(n)*(2.0d0*vpqrs(b,i,c,j)-vpqrs(c,i,b,j))
            ftmp=delta_ijbc*ftmp
            dyscoeff(b)=dyscoeff(b)+ftmp
         enddo
      enddo
        
      return
        
    end subroutine dyscoeff_gs_unocc

!#######################################################################
! dyscoeff_exci_occ: calculates the coefficients for the occupied MOs in
!                    the expansion of the Dyson orbital corresponding to
!                    ionisation from an excited state to a single
!                    IP-ADC state
!
! rhogs2   - 2nd-order correction to the ground state density matrix
! dyscoeff - Dyson orbital coefficients in the MO basis
! vec      - IP-ADC state vector of interest
!
!#######################################################################

    subroutine dyscoeff_exci_occ(rhogs2,dyscoeff,kpq,kpqf,vec_init,&
         vec_final,ndim,ndims,ndimf,ndimsf,smat,pmat)

      use constants
      use parameters
      use vpqrsmod

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
      integer                                   :: ndim,ndims,ndimf,&
                                                   ndimsf,alpha,n,m,&
                                                   i,j,a,b
      real(d), dimension(nbas,nbas)             :: rhogs2,smat,pmat
      real(d), dimension(nbas)                  :: dyscoeff
      real(d), dimension(ndimf)                 :: vec_init
      real(d), dimension(ndimf)                 :: vec_final

!-----------------------------------------------------------------------
! Index of the continuum orbital
!-----------------------------------------------------------------------
      alpha=ifakeorb

!-----------------------------------------------------------------------
! Zeroth-order term
!-----------------------------------------------------------------------
      do i=1,nocc
         dyscoeff(i)=pmat(alpha,i)
      enddo

!-----------------------------------------------------------------------
! Second-order contributions
!
! Note that all first-order constributions vanish
!-----------------------------------------------------------------------
      ! Term H
      do i=1,nocc
         do b=nocc+1,nbas
            dyscoeff(i)=dyscoeff(i)+0.5d0*rhogs2(i,b)*pmat(alpha,b)
         enddo
      enddo

!      ! Term K
!      ! Won't S_alpha,alpha always be zero?
!      do n=1,ndimsf
!         i=kpqf(3,n)
!         dyscoeff(i)=dyscoeff(i)-invsqrt2*vec_final(n)*smat(alpha,alpha)
!      enddo

      ! Term L
      do i=1,nocc
         do n=1,ndimsf
            j=kpqf(3,n)
            dyscoeff(i)=dyscoeff(i)-invsqrt2*vec_final(n)*smat(i,j)
         enddo
      enddo      

      return

    end subroutine dyscoeff_exci_occ

!#######################################################################

!#######################################################################
! dyscoeff_exci_occ: calculates the coefficients for the unoccupied MOs 
!                    in the expansion of the Dyson orbital corresponding
!                    to ionisation from an excited state to a single
!                    IP-ADC state
!
! rhogs2   - 2nd-order correction to the ground state density matrix
! dyscoeff - Dyson orbital coefficients in the MO basis
! vec      - IP-ADC state vector of interest
!
!#######################################################################

    subroutine dyscoeff_exci_unocc(rhogs2,dyscoeff,kpq,kpqf,vec_init,&
         vec_final,ndim,ndims,ndimf,ndimsf,rmat,pmat,qmat)

      use constants
      use parameters
      use vpqrsmod

      implicit none

      integer, dimension(7,0:nBas**2*4*nOcc**2) :: kpq,kpqf
      integer                                   :: ndim,ndims,ndimf,&
                                                   ndimsf,alpha,n,m,&
                                                   i,j,k,l,a,b,c,&
                                                   clbl
      real(d), dimension(nbas,nbas)             :: rhogs2,rmat,pmat,&
                                                   qmat,rmatf
      real(d), dimension(nbas)                  :: dyscoeff
      real(d), dimension(ndimf)                 :: vec_init
      real(d), dimension(ndimf)                 :: vec_final
      real(d)                                   :: ftmp,ftmp2,&
                                                   delta_ijab,&
                                                   delta_klad,&
                                                   delta_klbc,&
                                                   delta_ikac,&
                                                   delta_ikbc
      real(d), dimension(:,:,:), allocatable    :: tau1,tau2
      real(d), dimension(:), allocatable        :: zeta

!-----------------------------------------------------------------------
! R-matrix for the final state
!-----------------------------------------------------------------------
      rmatf=0.0d0
      do i=1,nocc
         do a=nocc+1,nbas
            do n=1,ndimsf
               j=kpqf(3,n)
               b=kpqf(5,n)
               delta_ijab=1.0d0/(e(a)+e(b)-e(i)-e(j))
               ftmp=2.0d0*vpqrs(a,i,b,j)-vpqrs(b,i,a,j)
               ftmp=delta_ijab*invsqrt2*vec_final(n)*ftmp
               rmatf(i,a)=rmatf(i,a)+ftmp
            enddo
         enddo
      enddo

!-----------------------------------------------------------------------
! Index of the continuum orbital
!-----------------------------------------------------------------------
      alpha=ifakeorb

!-----------------------------------------------------------------------
! Zeroth-order contribution
!-----------------------------------------------------------------------
      do b=nocc+1,nbas
         dyscoeff(b)=pmat(alpha,b)+qmat(alpha,b)
      enddo

!-----------------------------------------------------------------------
! Second-order contributions
!-----------------------------------------------------------------------
      ! Term 1
      do b=nocc+1,nbas
         do c=nocc+1,nbas
            ftmp=0.5d0*rhogs2(alpha,c)*pmat(c,b) &
                 +pmat(alpha,c)*0.5d0*rhogs2(c,b)
            dyscoeff(b)=dyscoeff(b)-0.5d0*ftmp
         enddo
      enddo

      ! Term 2
      do b=nocc+1,nbas
         do i=1,nocc
            dyscoeff(b)=dyscoeff(b)+rmatf(i,b)*rmat(i,alpha)
         enddo
      enddo

      ! Term 3
      !
      ! (i) Intermediate terms
      allocate(tau1(nvirt,nocc,nocc))
      allocate(tau2(nvirt,nocc,nocc))
      tau1=0.0d0
      tau2=0.0d0
      clbl=0
      do c=nocc+1,nbas
         clbl=clbl+1
         do k=1,nocc
            do l=1,nocc
               do a=nocc+1,nbas
                  
                  delta_klad=1.0d0/(e(alpha)+e(a)-e(k)-e(l))

                  ftmp=2.0d0*vpqrs(alpha,k,a,l)-vpqrs(a,k,alpha,l)
                  ftmp=delta_klad*pmat(c,a)*ftmp
                  tau1(clbl,k,l)=tau1(clbl,k,l)+ftmp

                  ftmp=2.0d0*vpqrs(a,k,alpha,l)-vpqrs(alpha,k,a,l)
                  ftmp=delta_klad*pmat(c,a)*ftmp
                  tau2(clbl,k,l)=tau2(clbl,k,l)+ftmp

               enddo
            enddo
         enddo
      enddo
      ! (ii) Calculation of term 3
      do b=nocc+1,nbas
         clbl=0
         do c=nocc+1,nbas
            clbl=clbl+1
            do k=1,nocc
               do l=1,nocc
                  delta_klbc=1.0d0/(e(b)+e(c)-e(k)-e(l))
                  ftmp=delta_klbc*vpqrs(k,b,l,c)*tau1(clbl,k,l)
                  ftmp=ftmp+delta_klbc*vpqrs(k,c,l,b)*tau2(clbl,k,l)
                  dyscoeff(b)=dyscoeff(b)-0.5d0*ftmp
               enddo
            enddo
         enddo
      enddo
      deallocate(tau1)
      deallocate(tau2)

      ! Term 5 (note that term 4 vanishes)
      !
      ! (i) Intermediate terms
      allocate(zeta(nocc))
      zeta=0.0d0
      do i=1,nocc
         do c=nocc+1,nbas
            do k=1,nocc
               delta_ikac=1.0d0/(e(alpha)+e(c)-e(i)-e(k))
               ftmp=2.0d0*vpqrs(alpha,i,c,k)-vpqrs(c,i,alpha,k)
               zeta(i)=zeta(i)+delta_ikac*rmatf(k,c)*ftmp
            enddo
         enddo
      enddo
      ! (ii) Calculation of term 5
      do n=1,ndim
         i=kpq(3,n)
         b=kpq(5,n)
         dyscoeff(b)=dyscoeff(b)+0.5d0*invsqrt2*vec_init(n)*zeta(i)
      enddo
      deallocate(zeta)

      ! Term 6
      do b=nocc+1,nbas
         do m=1,ndimsf
            i=kpqf(3,m)
            ftmp=0.0d0
            do c=nocc+1,nbas
               do k=1,nocc
                  delta_ikbc=1.0d0/(e(b)+e(c)-e(i)-e(k))
                  ftmp2=2.0d0*vpqrs(i,b,k,c)-vpqrs(i,c,k,b)
                  ftmp=ftmp+delta_ikbc*rmat(k,c)*ftmp2
               enddo
            enddo
            dyscoeff(b)=dyscoeff(b)+0.5d0*invsqrt2*vec_final(m)*ftmp
         enddo
      enddo

      return

    end subroutine dyscoeff_exci_unocc

!#######################################################################

  end module dyson_calc
