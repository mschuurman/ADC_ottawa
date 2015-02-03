      program dipole_calc

        use stieltjesmod

        implicit none


        integer           :: scf,iv,ov(maxbas),nov
        integer           :: nfin,nch,conv_max,ninit,noffd,nham
        integer           :: nhole_fin,hole_fin(maxbas)
        integer           :: sym_fin,nhf,ncontim,nhole_in,np
        integer           :: i,j,k
        real*8            :: ecutoff
        real*8            :: e_init,d_init(2)
        real*8            :: xpqrs
        real*16           :: overmax
        character(LEN=20) :: filename

!  Allocatable arrays
        integer, dimension(:,:), allocatable            :: kpq_fin,kpq_in
        double precision, dimension(:,:), allocatable   :: ham2h1p,evec
        double precision, dimension(:,:,:), allocatable :: d_matr,d_matr0
        double precision, dimension(:), allocatable     :: e_fin
        double precision, dimension(:), allocatable     :: vec_init
        double precision, dimension(:,:), allocatable   :: coupling_1p
        double precision, dimension(:), allocatable     :: ediff,dijlen,dijvel
        double precision, dimension(:), allocatable     :: idiff


      namelist /Input/ scf,iv,adc2,eps_min,eps_max,coup_threshold,&
           sym_fin,overmax,ecutoff,filename
! here we read and print the values of the input list variables
      read(5,nml=Input)
      write(6,*) 'SCF = ',scf,' (0 or 131072 = MOLCAS, 1 = GUK)'
      write(6,*) 'IV =',iv
      write(6,*) 'ADC2 =',adc2
      write(6,*) 'EPS_MIN =',eps_min
      write(6,*) 'EPS_MAX =',eps_max
      write(6,*) 'COUP_THRESHOLD =',coup_threshold
      write(6,*) 'OVERMAX =',overmax
      write(6,*)'SYM_FIN =',sym_fin ! symmetry for the final states
      write(6,*) 'ECUTOFF =',ecutoff ! energy to remove discrete states below ioniation threshold
      write(6,*) 'FILENAME =',filename  
      write(6,*) ''

      open(unit=4,file=filename)

      READ(4,*) np ! number of series
       
      write(6,*) 'numberpoints=',np
      allocate(ediff(np))
      allocate(dijlen(np))
 
      do i=1,np
         READ(4,*) ediff(i), dijlen(i)
      end do

      call stieltjes (np,ediff,dijlen,1,overmax,'grace.dat','int.dat')
      
      deallocate(ediff,dijlen)
      
      close(4)
      
      stop
    
    end program dipole_calc
