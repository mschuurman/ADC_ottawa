!#######################################################################
! dyson_io: subroutines to write the ADC/IP-ADC Dyson orbitals to input 
!           files for external programs to be used in further analysis
!#######################################################################

  module dyson_io

    save

    private :: dp
    
    ! Annoyingly, the gamess_internal module contains a variable
    ! named 'd', so we will use 'dp' here instead
    integer, parameter :: dp=selected_real_kind(8)

    contains

!#######################################################################
      
      subroutine pre_dysout(gam)

        use constants
        use parameters
        use gamess_internal

        implicit none

        integer             :: i,j,p1,p2,np
        real(dp)            :: alpha
        real(dp), parameter :: zettol=1e-15_dp
        logical             :: ldir
        type(gam_structure) :: gam

!-----------------------------------------------------------------------
! Determine where the 'continuum' orbital is centred
!-----------------------------------------------------------------------
        ! Determine if there is a dummy atom holding only the 
        ! 'continuum' orbital
        lrmatom=.false.
        ccent=0
        do i=1,gam%natoms
           if (gam%atoms(i)%nshell.eq.1) then
              p1=gam%atoms(i)%sh_p(j)
              alpha=gam%atoms(i)%p_zet(p1)
              if (alpha.lt.zettol) then 
                 ccent=i
                 lrmatom=.true.
              endif
           endif
        enddo
        
        ! If there isn't a dummy atom holding only the 'continuum'
        ! orbital, then determine the index of the 'continuum' orbital
        if (.not.lrmatom) then
           do i=1,gam%natoms
              do j=1,gam%atoms(i)%nshell
                 p1=gam%atoms(i)%sh_p(j)
                 p2=gam%atoms(i)%sh_p(j+1)
                 np=p2-p1
                 if (np.eq.1) then
                    alpha=gam%atoms(i)%p_zet(p1)
                    if (alpha.lt.zettol) then 
                       ccent=i
                       cshell=j
                    endif
                 endif
              enddo
           enddo
        endif

!-----------------------------------------------------------------------
! Transformation matrix for the reordering of the f-functions to
! comply with the molden ordering of:
!  
!  1    2    3    4    5    6    7    8    9   10
! xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
!
! Note that this ordering can also be used for ezdyson input files.
!-----------------------------------------------------------------------
        ftransmat=0.0d0
        ftransmat(1,1)=1.0d0
        ftransmat(2,2)=1.0d0
        ftransmat(3,3)=1.0d0
        ftransmat(6,4)=1.0d0
        ftransmat(4,5)=1.0d0
        ftransmat(5,6)=1.0d0
        ftransmat(8,7)=1.0d0
        ftransmat(9,8)=1.0d0
        ftransmat(7,9)=1.0d0
        ftransmat(10,10)=1.0d0

!-----------------------------------------------------------------------
! Create output directories
!-----------------------------------------------------------------------
        if (dysout(1).eq.1) then
           inquire(file='molden/.',exist=ldir)
           if (ldir) call system('rm -rf molden')
           call system('mkdir molden')
        endif

        if (dysout(2).eq.1) then
           inquire(file='ezdyson/.',exist=ldir)
           if (ldir) call system('rm -rf ezdyson')
           call system('mkdir ezdyson')
        endif

        return
        
      end subroutine pre_dysout

!#######################################################################

      subroutine wrmolden(dyscoeff,gam,n)

        use constants
        use parameters
        use import_gamess
        use iomod, only: freeunit
        use gamess_internal

        implicit none

        integer, intent(in)          :: n
        integer                      :: imolden,i,j,k,iang,p1,p2,np,&
                                        count,pk,norb,atmcnt,irmfunc
        real(dp), dimension(nbas)    :: dyscoeff
        real(dp), dimension(nbas_ao) :: dyscoeff_ao
        real(dp)                     :: nfac,alpha,coeff
        character(len=60)            :: filename
        type(gam_structure)          :: gam

!-----------------------------------------------------------------------
! Set shell labels
!-----------------------------------------------------------------------
        shlbl(0:3)=(/ 's','p','d','f' /)

!-----------------------------------------------------------------------
! Open the molden file
!-----------------------------------------------------------------------
        call freeunit(imolden)
        
        call outfilename(n,filename,1)

        open(imolden,file=filename,form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Preamble
!-----------------------------------------------------------------------
        write(imolden,'(a)') '[Molden Format]'
        write(imolden,'(a)') '[Title]'
        
!-----------------------------------------------------------------------
! Cartesian coordinates
!-----------------------------------------------------------------------
        write(imolden,'(a)') '[Atoms] Angs'
        atmcnt=0
        do i=1,gam%natoms
           
           ! Skip the current atom if it is a dummy atom holding only
           ! the 'continuum' orbital
           if (lrmatom.and.i.eq.ccent) cycle

           atmcnt=atmcnt+1

           write(imolden,'(1x,a2,2x,i3,2x,i3,3(2x,F10.7))') &
                gam%atoms(i)%name,atmcnt,int(gam%atoms(i)%znuc),&
                (gam%atoms(i)%xyz(j),j=1,3)

        enddo

!-----------------------------------------------------------------------
! AO basis set section
!-----------------------------------------------------------------------
        atmcnt=0

        write(imolden,'(a)') '[GTO]'
        ! Loop over atoms
        do i=1,gam%natoms

           ! Skip if we are at a dummy atom holding only the
           ! 'continuum' orbital
           if (lrmatom.and.i.eq.ccent) cycle

           atmcnt=atmcnt+1

           ! Primitive counter for the current atom: pk
           pk=0

           ! Atom number
           if (atmcnt.gt.1) write(imolden,*)
           write(imolden,'(2x,i3,2x,a1)') atmcnt,'0'

           ! Loop over shells for the current atom
           do j=1,gam%atoms(i)%nshell

              ! Skip if we are at the shell containing the 'continuum'
              ! orbital
              if (i.eq.ccent.and.j.eq.cshell) then
                 pk=pk+1
                 cycle
              endif

              ! Angular momentum quantum number and primitive indices
              ! for the current shell
              iang=gam%atoms(i)%sh_l(j)
              p1=gam%atoms(i)%sh_p(j)
              p2=gam%atoms(i)%sh_p(j+1)
              np=p2-p1
              
              ! Label and no. primitives for the current shell
              write(imolden,'(2x,a1,2x,i3,2x,a4)') shlbl(iang),&
                   np,'1.00'
              
              ! Primitive exponents and coefficients for the current 
              ! shell.
              ! Note that we output the original, unscaled,
              ! non-normalised primitive coefficients
              do k=1,np
                 pk=pk+1
                 alpha=gam%atoms(i)%p_zet(p1-1+k)
                 coeff=gam%atoms(i)%p_c_orig(p1-1+k)
                 write(imolden,*) alpha,coeff
              enddo

           enddo
        enddo

!-----------------------------------------------------------------------
! Write the *normalised* Dyson orbital to the molden file
!-----------------------------------------------------------------------
        ! Calculate the *normalised* Dyson orbital expansion 
        ! coefficients in the AO basis
        nfac=sqrt(dot_product(dyscoeff,dyscoeff))
        dyscoeff_ao=matmul(ao2mo(1:nbas_ao,1:nbas),dyscoeff(1:nbas)/nfac)
        
        ! Write the *normalised* Dyson orbital to the molden file
        write(imolden,'(/,a)') '[MO]'
        write(imolden,'(a)') 'Sym= 1'
        write(imolden,'(a,x,F6.1)') 'Ene=',real(n)
        write(imolden,'(a)') 'Spin= Alpha'
        write(imolden,'(a,E14.4)') 'Occup= ',nfac
        count=0
        irmfunc=0
        do i=1,gam%natoms

           ! Skip if we are at a dummy atom holding only the
           ! 'continuum' orbital
           if (lrmatom.and.i.eq.ccent) then
              count=count+1
              irmfunc=1
              cycle
           endif

           do j=1,gam%atoms(i)%nshell

              ! Skip if we are at the shell containing the 'continuum'
              ! orbital
              if (i.eq.ccent.and.j.eq.cshell) then
                 count=count+1
                 irmfunc=1
                 cycle
              endif

              ! Angular momentum quantum number and primitive indices
              ! for the current shell
              iang=gam%atoms(i)%sh_l(j)
              norb=gam_orbcnt(iang)
              p1=count+1
              p2=count+norb

              ! If we are at a shell of f-functions, then rearrange
              ! to correspond to the molden ordering
              if (iang.eq.3) then
                 dyscoeff_ao(p1:p2)=matmul(ftransmat,dyscoeff_ao(p1:p2))
              endif

              do k=ang_loc(iang),ang_loc(iang)+norb-1                
                 count=count+1
                 ! Note that here k tells us what type of function we
                 ! are at (s, px, py,...) and count indexes the current
                 ! AO basis function
                 !
                 ! Also, the coefficient that we output has to
                 ! correspond to a normalised AO - hence the
                 ! multiplication by ang_c(k)
                 write(imolden,*) count-irmfunc,dyscoeff_ao(count)*ang_c(k)
              enddo
              
           enddo
        enddo

!-----------------------------------------------------------------------
! Close the molden file
!-----------------------------------------------------------------------
        close(imolden)

        return

      end subroutine wrmolden

!#######################################################################

      subroutine wrezdyson(dyscoeff,gam,n,de)

        use constants
        use parameters
        use import_gamess
        use iomod, only: freeunit
        use gamess_internal

        implicit none

        integer, intent(in)               :: n
        integer                           :: iezd,i,j,k,m,pk,iang,p1,p2,&
                                             np,count,norb
        real(dp), dimension(nbas)         :: dyscoeff
        real(dp), dimension(nbas_ao)      :: dyscoeff_ao
        real(dp)                          :: alpha,dnorm,si,sf,coeff,de
        character(len=60)                 :: filename
        character(len=60), dimension(2)   :: comment
        type(gam_structure)               :: gam

!-----------------------------------------------------------------------
! Set shell labels
!-----------------------------------------------------------------------
        shlbl(0:3)=(/ 'S','P','D','F' /)

!-----------------------------------------------------------------------
! Open the ezdyson input file
!-----------------------------------------------------------------------
        call freeunit(iezd)

        call outfilename(n,filename,2)
        
        open(iezd,file=filename,form='formatted',status='unknown')

!-----------------------------------------------------------------------
! Preamble
!-----------------------------------------------------------------------
        write(iezd,'(a)') '<?xml version="1.0" encoding="ISO-8859-1"?>'
        write(iezd,'(/,a,2(/,2x,a))') '<root','job = "dyson"','>'

!-----------------------------------------------------------------------
! Geometry
!-----------------------------------------------------------------------
        if (lrmatom) then
           m=gam%natoms-1
        else
           m=gam%natoms
        endif

        write(iezd,'(/,a,/,x,a,i2,a1,/,x,a)') '<geometry',&
             'n_of_atoms="',m,'"','text = "'
        do i=1,gam%natoms

           ! Skip the current atom if it is a dummy atom holding only
           ! the 'continuum' orbital
           if (lrmatom.and.i.eq.ccent) cycle

           write(iezd,'(6x,a2,1x,3(3x,F14.10))') gam%atoms(i)%name,&
                (gam%atoms(i)%xyz(j),j=1,3)

        enddo

        write(iezd,'(/,2x,a,/,a)') '"','/>'

!-----------------------------------------------------------------------
! Photoelectron wavefunction
!-----------------------------------------------------------------------
        write(iezd,'(/,a)') '<free_electron'
        write(iezd,'(3x,a,1x,i2,a1)') 'l_max = "',lmax,'"'
        write(iezd,'(3x,a,F5.3,a)') 'charge_of_ionized_core = "',&
             zcore,'">'
        write(iezd,'(a,i2,a,2(1x,a,F7.3,a),1x,a)') &
             '<k_grid n_points="',nelen,'"','min="',eleni,'"','max="',&
             elenf,'"','/>'
        write(iezd,'(a)') '</free_electron>'

!-----------------------------------------------------------------------
! Averaging over molecular orientations
!-----------------------------------------------------------------------
        write(iezd,'(/,a)') '<averaging'
        write(iezd,'(3x,a)') 'method = "avg"'
        write(iezd,'(3x,a)') 'method_possible_values = "noavg, avg">'
        write(iezd,'(a)') '</averaging>'

!-----------------------------------------------------------------------
! Ionization energy - we use a value of zero here so that the 
! calculated cross-sections are given as a function of the 
! photoelectron kinetic energy
!-----------------------------------------------------------------------
        write(iezd,'(/,a)') '<laser'
        write(iezd,'(3x,a,F6.2,a)') 'ionization_energy = "',de,'" >'
        write(iezd,'(3x,a)') &
             '<laser_polarization x="0.0" y="0.0" z="1.0" />'
        write(iezd,'(a)') '</laser>'

!-----------------------------------------------------------------------
! Cartesian grid
!-----------------------------------------------------------------------
        write(iezd,'(/,a)') '<lab_xyz_grid>'
        write(iezd,'(3x,a,i4,a,2(1x,a,F7.3,a))') '<axis n_points="',&
             ngrdpnts,'"','min="',grdi,'"','max="',grdf,'" />'
        write(iezd,'(a)') '</lab_xyz_grid>'

!-----------------------------------------------------------------------
! Job parameters
!-----------------------------------------------------------------------
        write(iezd,'(/,a)') '<job_parameters'
        write(iezd,'(3x,a)') 'unrestricted = "false"'
        write(iezd,'(3x,a)') 'Dyson_MO_transitions = "1"'
        write(iezd,'(3x,a)') 'spin_degeneracy = "2"'
        write(iezd,'(3x,a)') 'orbital_degeneracy = "1"'
        write(iezd,'(3x,a)') 'number_of_MOs_to_plot = "0"'
        write(iezd,'(3x,a)') 'MOs_to_plot = ""'
        write(iezd,'(a)') '/>'

!-----------------------------------------------------------------------
! AO basis
!-----------------------------------------------------------------------
        write(iezd,'(/,a)') '<basis'
        write(iezd,'(3x,a,i3,a)') 'n_of_basis_functions="',&
             gam%nbasis-1,'"'
        write(iezd,'(3x,a)') 'AO_ordering = "Molden"'
        write(iezd,'(3x,a)') 'purecart = "2222">'
        
        do i=1,gam%natoms

           ! Skip if we are at a dummy atom holding only the
           ! 'continuum' orbital
           if (lrmatom.and.i.eq.ccent) cycle

           write(iezd,'(3x,a,/,3x,a)') '<atom','text = "'
           ! Primitive counter for the current atom: pk
           pk=0
           ! Atomic symbol
           write(iezd,'(a2,4x,a)') gam%atoms(i)%name,'0'           
           ! Loop over shells for the current atom
           do j=1,gam%atoms(i)%nshell

              ! Skip if we are at the shell containing the 'continuum'
              ! orbital
              if (i.eq.ccent.and.j.eq.cshell) then
                 pk=pk+1
                 cycle
              endif

              ! Angular momentum quantum number and primitive indices
              ! for the current shell
              iang=gam%atoms(i)%sh_l(j)
              p1=gam%atoms(i)%sh_p(j)
              p2=gam%atoms(i)%sh_p(j+1)
              np=p2-p1

              ! Label and no. primitives for the current shell
              write(iezd,'(a1,4x,i3,4x,a)') shlbl(iang),np,'1.00'

              ! Primitive exponents and coefficients for the current 
              ! shell.
              ! Note that we output the original, unscaled,
              ! non-normalised primitive coefficients
              do k=1,np
                 pk=pk+1
                 alpha=gam%atoms(i)%p_zet(p1-1+k)
                 coeff=gam%atoms(i)%p_c_orig(p1-1+k)
                 write(iezd,'(2(3x,E14.4))') alpha,coeff
              enddo

           enddo
           write(iezd,'(a,/,7x,a,/,a)') '****','"','/>'
        enddo

        write(iezd,'(a)') '</basis>'

!-----------------------------------------------------------------------
! Dyson orbital
!
! Note that ezdyson appears to assume that the Dyson orbitals are
! from an EOM-CCSD calculation, and as such requires the input of
! both left and right Dyson orbitals. As the ADC Hamiltonian is
! Hermitian, we only have one Dyson orbital for each transition, and
! this will be entered as both the 'left' and 'right' Dyson orbitals.
!-----------------------------------------------------------------------
        write(iezd,'(/,a)') '<dyson_molecular_orbitals>'

        dnorm=sqrt(dot_product(dyscoeff,dyscoeff))
        si=real(statenumber)+real(nirrep)/10.0d0
        sf=real(n)+real(dysirrep)/10.0d0
        comment(1)='comment="dyson right"'
        comment(2)='comment="dyson left"'

        do m=1,2
           write(iezd,'(3x,a,F8.6,a,x,2(F6.1,x,a),x,a)') &
                '<DMO norm="',dnorm,'" transition="',si,'to',sf,'"',&
                trim(comment(m))           
           write(iezd,'(3x,a)') 'text="'
           ! Convert the Dyson orbital to the AO basis and normalise
           dyscoeff_ao=matmul(ao2mo(1:nbas_ao,1:nbas),dyscoeff(1:nbas)/dnorm)
           count=0
           do i=1,gam%natoms

              ! Skip if we are at a dummy atom holding only the
              ! 'continuum' orbital
              if (lrmatom.and.i.eq.ccent) then
                 count=count+1
                 cycle
              endif

              do j=1,gam%atoms(i)%nshell
                 
                 ! Skip if we are at the shell containing the 'continuum'
                 ! orbital
                 if (i.eq.ccent.and.j.eq.cshell) then
                    count=count+1
                    cycle
                 endif

                 ! Angular momentum quantum number and primitive indices
                 ! for the current shell
                 iang=gam%atoms(i)%sh_l(j)
                 norb=gam_orbcnt(iang)
                 p1=count+1
                 p2=count+norb

                 ! If we are at a shell of f-functions, then rearrange
                 ! to correspond to the molden ordering
                 if (iang.eq.3) then
                    dyscoeff_ao(p1:p2)=matmul(ftransmat,dyscoeff_ao(p1:p2))
                 endif

                 do k=ang_loc(iang),ang_loc(iang)+norb-1
                    count=count+1
                    ! Note that here k tells us what type of function we
                    ! are at (s, px, py,...) and count indexes the current
                    ! AO basis function
                    !
                    ! Also, the coefficient that we output has to
                    ! correspond to a normalised AO - hence the
                    ! multiplication by ang_c(k)
                    write(iezd,'(3x,E14.4)') dyscoeff_ao(count)*ang_c(k)
                 enddo

              enddo
           enddo
           write(iezd,'(3x,a,/)') '"  />'
        enddo

        write(iezd,'(a)') '</dyson_molecular_orbitals>'

!-----------------------------------------------------------------------
! End
!-----------------------------------------------------------------------
        write(iezd,'(/,a)') '</root>'

!-----------------------------------------------------------------------
! Close the ezdyson input file
!-----------------------------------------------------------------------
        close(iezd)

        return

      end subroutine wrezdyson

!#######################################################################

      subroutine outfilename(n,filename,flag)

        use parameters

        implicit none

        integer           :: n,flag
        character(len=60) :: filename
        
        filename=''

        if (flag.eq.1) then
           if (n.lt.10) then
              write(filename,'(a14,i1,a7)') 'molden/dysorb_',n,'.molden'
           else if (n.lt.100) then
              write(filename,'(a14,i2,a7)') 'molden/dysorb_',n,'.molden'
           else if (n.lt.1000) then
              write(filename,'(a14,i3,a7)') 'molden/dysorb_',n,'.molden'
           else
              write(filename,'(a14,i4,a7)') 'molden/dysorb_',n,'.molden'
           endif
        else if (flag.eq.2) then
           if (n.lt.10) then
              write(filename,'(a15,i1,a4)') 'ezdyson/dysorb_',n,'.xml'
           else if (n.lt.100) then
              write(filename,'(a15,i2,a4)') 'ezdyson/dysorb_',n,'.xml'
           else if (n.lt.1000) then
              write(filename,'(a15,i3,a4)') 'ezdyson/dysorb_',n,'.xml'
           else
              write(filename,'(a15,i4,a4)') 'ezdyson/dysorb_',n,'.xml'
           endif
        endif

        return

      end subroutine outfilename

!#######################################################################

    end module dyson_io
