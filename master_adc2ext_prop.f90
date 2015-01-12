subroutine master_adc2ext_prop()

  use constants
  use parameters
  use select_fano
  use davmod
  use fspace
  use get_moment
  use misc
  use fspace2
  use get_matrix
  use get_matrix_DIPOLE
  use propagate_prepare
 
  implicit none

  integer, dimension(:,:), allocatable :: kpq,kpqd,kpqf
  integer :: i,j,ndim,ndims,ndimsf,nout,nstates,ndimf,ndimd,noutf
  integer*8 :: noffd,noffdf  
  integer :: k,l,m,n,k1,b1,b


  real(d) :: time,enerstate
  real(d), dimension(:), allocatable :: ener,enerdav,mtm,tmvec,osc_str
  real(d), dimension(:), allocatable :: enerf,tmvecf,osc_strf
  real(d), dimension(:), allocatable :: autvec,travec,excit,coeff,coeff_tra
  INTEGER, dimension(:), allocatable :: indx,indx_tra
  real(d), dimension(:,:), allocatable :: arr,arrd,arrf 
  real(d) :: e_init
  real(d), dimension(:,:), allocatable :: rvec
  real(d), dimension(:), allocatable :: vec_init 
  real(d), dimension(:), allocatable :: travec_norm
  real(d) :: norma  
  INTEGER, dimension(:), allocatable :: MAP
  INTEGER :: count
  INTEGER, dimension(:), allocatable :: DIMEN
  real*8, dimension(:), allocatable :: mtmf
  INTEGER :: TOTDIM_PROP, DIV_DIP, DIV_HAM, NBLOCKS_DIP

  real(d) :: E_groundstate
  
  logical :: prs,prs1,prs2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! INFO=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (info .eq. 1) then
     allocate(kpq(7,0:nBas**2*4*nOcc**2))
     allocate(kpqf(7,0:nBas**2*4*nOcc**2))
     allocate(kpqd(7,0:nBas**2*4*nOcc**2))

  kpq(:,:)=-1
       call  select_atom_is(kpq(:,:))
       call  select_atom_d(kpq(:,:),-1)

  kpqf(:,:)=-1
       call  select_atom_isf(kpqf(:,:))
       call  select_atom_df(kpqf(:,:),-1)

  kpqd(:,:)=-1
       call  select_atom_ist(kpqd(:,:))
       call  select_atom_dt(kpqd(:,:),-1)

       ndim=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)

       ndimf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+2*kpqf(5,0)

       ndimd=kpqd(1,0)+kpqd(2,0)+kpqd(3,0)+kpqd(4,0)+2*kpqd(5,0)

       write(6,*)
       write(6,*) 'ADC(2)ext INITIAL Space dim',ndim
       write(6,*)
       write(6,*) 'ADC(2)ext FINAL Space dim',ndimf
       write(6,*)
       write(6,*) 'ADC(2)ext TOTAL Space dim WITHOUT GROUND STATE',ndimd
       write(6,*)


!!!!!!!!!!!!!!!!!!!!!! CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(6,*) 'dimension of various INITIAL configuration spaces'
       write(6,*) '      1p1h       2p2h_1      2p2h_2      2p2h_3      2p2h_4i      2p2h_4ii'
       write(6,*) kpq(1,0),kpq(2,0),kpq(3,0),kpq(4,0),kpq(5,0),kpq(5,0)
       write(6,*)

       write(6,*) 'dimension of various FINAL configuration spaces'
       write(6,*) '      1p1h       2p2h_1      2p2h_2      2p2h_3      2p2h_4i      2p2h_4ii'
       write(6,*) kpqf(1,0),kpqf(2,0),kpqf(3,0),kpqf(4,0),kpqf(5,0),kpqf(5,0)
       write(6,*)

       write(6,*) 'dimension of various TOTAL configuration spaces'
       write(6,*) '      1p1h       2p2h_1      2p2h_2      2p2h_3      2p2h_4i      2p2h_4ii'
       write(6,*) kpqd(1,0),kpqd(2,0),kpqd(3,0),kpqd(4,0),kpqd(5,0),kpqd(5,0)
       write(6,*)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call cpu_time(time)
     write(6,*) 'Time=',time," s"

     deallocate(kpq,kpqd,kpqf)
     stop

  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! INFO=0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  allocate(kpq(7,0:nBas**2*4*nOcc**2))
  allocate(kpqd(7,0:nBas**2*4*nOcc**2))
  allocate(kpqf(7,0:nBas**2*4*nOcc**2))

  kpq(:,:)=-1

  call  select_atom_is(kpq(:,:))
  call  select_atom_d(kpq(:,:),-1)

  kpqf(:,:)=-1

  call  select_atom_isf(kpqf(:,:))
  call  select_atom_df(kpqf(:,:),-1)

  kpqd(:,:)=-1

  call  select_atom_ist(kpqd(:,:))
  call  select_atom_dt(kpqd(:,:),-1)

  ndim=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)
  ndimf=kpqf(1,0)+kpqf(2,0)+kpqf(3,0)+kpqf(4,0)+2*kpqf(5,0)
  ndimd=kpqd(1,0)+kpqd(2,0)+kpqd(3,0)+kpqd(4,0)+2*kpqd(5,0)

  write(6,*) 'ADC(2)ext INITIAL Space dim',ndim
  write(6,*) 'ADC(2)ext FINAL Space dim',ndimf
  write(6,*) 'ADC(2)ext TOTAL Space dim WITHOUT GROUND STATE',ndimd

!!!!!!!!!!!!!!!!!!!!!! CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*) 'dimension of various INITIAL configuration spaces'
   write(6,*) kpq(1,0),kpq(2,0),kpq(3,0),kpq(4,0),kpq(5,0)
   write(6,*) 'dimension of various FINAL configuration spaces'
   write(6,*) kpqf(1,0),kpqf(2,0),kpqf(3,0),kpqf(4,0),kpqf(5,0)
   write(6,*) 'dimension of various TOTAL configuration spaces'
   write(6,*) kpqd(1,0),kpqd(2,0),kpqd(3,0),kpqd(4,0),kpqd(5,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nout=ndim
  noutf=ndimf
  ndims=kpq(1,0)
  ndimsf=kpqf(1,0)






IF ( WHAT .EQ. 'PROP' ) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 ALLOCATE(DIMEN(NSYMA))
 DIMEN(:) = 0

 CALL BUILD_ALL_SYMMETRY_SPACES( NSYMA , DIMEN )


 CALL SYM_TOPROP_CALC( DIPOLESYM , ELECTRIC_FIELD , NSYMA_PROP )

 allocate(SYM_MAP(NSYMA_PROP))
 allocate(DIM_PROP(NSYMA_PROP))

 CALL SYM_TOPROP_VECT( DIMEN , NSYMA_PROP , SYM_MAP , DIM_PROP )

!!! ********************************************************* !!!
 allocate(NBUF_SYM( NSYMA_PROP ))
!!! ********************************************************* !!!

!!! SCRIVE I PEZZI HAMILTONIANI SU FILE ( IL PEZZO I , NEL FILE [100 + I] ) !!!
 CALL  SAVE_HAMPIECES_ADC2EXT( DIMEN , NSYMA_PROP , SYM_MAP , DIM_PROP , TOTDIM_PROP, NBUF_SYM)
!!! SCRIVE I PEZZI HAMILTONIANI SU FILE ( IL PEZZO I , NEL FILE [100 + I] ) !!!

 CALL SYM_TOPROP_CALC_DIP2( DIPOLESYM , ELECTRIC_FIELD , NSYMA_PROP , SYM_MAP , NBLOCKS_DIP )
!!! RESTITUISCE NBLOCKS_DIP

!!! ********************************************************* !!!
 allocate(NREC_VECTOR( NBLOCKS_DIP ))
!!! ********************************************************* !!!

!!! SCRIVE I PEZZI DIPOLO SU FILE ( IL PEZZO [I], NEL FILE [100 + NSYMA_PROP + I] ) !!!
 CALL SAVE_DIPOLE_PIECES_ADC2( DIMEN , NSYMA_PROP , DIM_PROP , SYM_MAP , DIPOLESYM , ELECTRIC_FIELD , NBLOCKS_DIP , NREC_VECTOR )
!!! SCRIVE I PEZZI DIPOLO SU FILE ( IL PEZZO [I], NEL FILE [100 + NSYMA_PROP + I] ) !!!



 allocate(NDIV_DIP(NBLOCKS_DIP))
 NDIV_DIP(:) = DIV_DIP
 allocate(NDIV(NSYMA_PROP))
 NDIV(:) = DIV_HAM


!!! CALCOLA QUANTI PEZZI DI DIPOLO AVRO DA PROPAGARE, IL NUMERO DI PROCESSORI NECESSARI PER I DIPOLI !!! 
 CALL SYM_TOPROP_CALC_DIP( DIPOLESYM, ELECTRIC_FIELD, NSYMA_PROP , DIM_PROP , SYM_MAP , NDIV_DIP , NDIV , NBLOCKS_DIP , NBUF_SYM , NREC_VECTOR , KLPDTOT , HAM_PIECES )
!!! RESTITUISCE KLPDTOT


write(*,*) 'NUMBER OF HAMILTONIAN PIECES IN TOTAL =', HAM_PIECES
write(*,*) 'NUMBER OF DIPOLE PIECES IN TOTAL ='     , KLPDTOT

!!! ********************************************************* !!!
ALLOCATE(NRECTOT_VECT( NSYMA_PROP + NBLOCKS_DIP ))
!!! ********************************************************* !!!

!!! ********************************************************* !!!
ALLOCATE(NREC_VECTOR_BIS( HAM_PIECES + KLPDTOT ))
!!! ********************************************************* !!!

!!! ********************************************************* !!!
ALLOCATE(RECINI_VECT( HAM_PIECES + KLPDTOT ))
!!! ********************************************************* !!!

ALLOCATE(SIMMETRIA(HAM_PIECES))
ALLOCATE(NLPD(NSYMA_PROP))
ALLOCATE(LPSYM(KLPDTOT))
ALLOCATE(MU1D(KLPDTOT))
ALLOCATE(MU2D(KLPDTOT))
ALLOCATE(KAPPAD(KLPDTOT))
ALLOCATE(DIPOLE_BLOCK(KLPDTOT))
ALLOCATE(NUM(KLPDTOT))
ALLOCATE(NUM_DIAG(HAM_PIECES))

!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! TUTTE LE QUANTITA TENGONO CONTO DEL GROUND STATE IN SYMMETRY 1 !!! 
CALL SYM_TOPROP_VECT_DIP(NSYMA_PROP,SYM_MAP,DIM_PROP,NDIV,DIPOLESYM,ELECTRIC_FIELD,NBLOCKS_DIP,NDIV_DIP,HAM_PIECES,KLPDTOT,SIMMETRIA,NLPD,LPSYM,MU1D,MU2D,KAPPAD,DIPOLE_BLOCK,NUM,NUM_DIAG,NBUF_SYM,NREC_VECTOR,NRECTOT_VECT,NREC_VECTOR_BIS, RECINI_VECT )
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! TUTTE LE QUANTITA TENGONO CONTO DEL GROUND STATE IN SYMMETRY 1 !!! 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!














ELSE IF ( WHAT .EQ. 'CROS' ) THEN


  if (tranmom2 .eq. 'x') then
     dpl(:,:)=x_dipole(:,:)
  elseif (tranmom2 .eq. 'y') then
     dpl(:,:)=y_dipole(:,:)
  elseif (tranmom2 .eq. 'z') then
     dpl(:,:)=z_dipole(:,:)
  end if


  CHECK_dip = nirrep2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! FULL DIAGONALIZATION CASE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
IF (chrun .eq. 'dire') THEN
write(*,*) 'I WILL PERFORM FULL ADC2 DIAGONALIZATION IN THE FINAL SPACE'
 
IF (chrun2.eq.'full') THEN
write(*,*) 'I WILL PERFORM FULL ADC2 DIAGONALIZATION IN THE FINAL SPACE'
   



!!! ALLOCATION PART !!!
!!! ALLOCATION PART !!!
!!! ALLOCATION PART !!!
   allocate(arr(ndim,ndim),mtm(ndim),tmvec(nout),osc_str(nout))
   allocate(arrf(ndimf,ndimf),enerf(ndimf),tmvecf(noutf),osc_strf(noutf))
   allocate(arrd(ndimf,ndim),autvec(ndim),travec(ndimf))
   allocate(ener(ndim))
   allocate(mtmf(ndimf))
!!! ALLOCATION PART !!!
!!! ALLOCATION PART !!!
!!! ALLOCATION PART !!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*) "DIAGONALIZING IN THE INITIAL SYMMETRY  SPACE" 
   call get_fspace_adc2e_direct(ndim,kpq(:,:),arr,ener)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*) 'Calculating ADC2 transition moments from ground state in ',tranmom2,' direction.'
   write(6,*) 'I AM Calculating ADC2 transition moments in', tranmom2,' direction, i.e. : < PSI0  D',tranmom2,' PSIn  > IN THE INITIAL SPACE'
   call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:))
   write(100,*) ndim , mtm
   do i = 1 , ndim
   tmvec(i) = tm(ndim,arr(:,i),mtm(:))
   osc_str(i) = 2._d/3._d*ener(i)*tmvec(i)**2
   write(6,*) i,ener(i),tmvec(i),os2cs*osc_str(i)
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nout = 100
   write(6,*) "ADC2 diagonalization result in the INITIAL space"
   call table2(ndim,nout,ener(1:nout),arr(:,1:nout),tmvec(:),osc_str(:))
   write(6,*) ' sums calculated with respect to the ground state'
   call get_sigma(ndim,ener(:),os2cs*osc_str(:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   deallocate(mtm)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*) "DIAGONALIZING IN THE FINAL SYMMETRY  SPACE" 
   call get_fspace_adc2e_direct(ndimf,kpqf(:,:),arrf,enerf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*) 'Calculating ADC2 transition moments from ground state in ',tranmom2,' direction.'
   write(6,*) 'I AM Calculating ADC2 transition moments in', tranmom2,' direction, i.e. : < PSI0  D',tranmom2,' PSIm  > IN THE FINAL SPACE'
   call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtmf(:))
   write(100,*) ndimf , mtmf
   do i = 1 , ndimf
   tmvecf(i) = tm(ndimf,arrf(:,i),mtmf(:))
   osc_strf(i) = 2._d/3._d*enerf(i)*tmvecf(i)**2
   write(6,*) i,enerf(i),tmvecf(i),os2cs*osc_strf(i)
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nout = 100
   write(6,*) "ADC2 diagonalization result in the FINAL space"
   call table2(ndimf,nout,enerf(1:nout),arrf(:,1:nout),tmvecf(:),osc_strf(:))
   write(6,*) ' sums calculated with respect to the ground state'
   call get_sigma(ndimf,enerf(:),os2cs*osc_strf(:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   deallocate(mtmf)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!! INITIAL STATE VECTOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   write(6,*) statenumber,'SELECTED STATE IN THE INITIAL SPACE has been obtained' 
       allocate(coeff(ndim))
       allocate(indx(ndim))
       coeff(:)=arr(:,statenumber)**2
       call dsortindxa1("D",ndim,coeff(:),indx(:))
       write(105,*) 'DESCRIPTION OF THIS SELECTED STATE IN THE INITIAL SPACE'
       write(6,*) statenumber,ener(statenumber),ener(statenumber)*27.211396,(coeff(indx(j)),indx(j),j=1,5)
       write(105,*) 'ENERGY OF THIS SELECTED STATE IN THE INITIAL SPACE', ener(statenumber)
       deallocate(coeff)
       deallocate(indx)
   write(6,*) 'Calculating ADC2 transition moments TO THE FINAL SPACE STATES,  from THIS ', statenumber, 'EXCITE STATE IN THE INITIAL SPACE  in ',tranmom2,' direction.'
   write(6,*) 'I AM Calculating ADC2 transition moments TO THE FINAL SPACE STATES, in', tranmom2,' direction, i.e. : < PSIin(INITIALSPACE)  D',tranmom2,' PSIm(FINAL SPACE)  > '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!! MATRIX VECTOR MULTIPLICATION !!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   write(6,*) "DOING THE MATRIX VECTOR MULTIPLICATION TO OBTAIN CONTRACTING VECTOR"
   if (matvec.eq.1) then
   call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,arr(:,statenumber),travec)
   write(6,*) "first matrix-vector product done"
   else if (matvec.eq. 2) then
   call get_fspace_adc2_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,arr(:,statenumber),arrd,travec) 
   write(6,*) "second matrix-vector product done"
   else if (matvec.eq.3) then
   call get_fspace_adc2_DIPOLE_direct(ndim,kpq,arr(:,statenumber),arrd,travec) 
   write(6,*) "third matrix-vector product done"
   end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TRAVEC VECTOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*)  ' vector to make scalar product with has been calculated'
   write(101,*) ndimf , travec
   allocate(coeff_tra(ndimf))
   allocate(indx_tra(ndimf))
   coeff_tra(:)=travec(:)**2
   call dsortindxa1("D",ndimf,coeff_tra(:),indx_tra(:))
   write(6,*) "BIGGER COEFFICIENTS OF ! D",tranmom2," PSIin >"
   write(6,*) (coeff_tra(indx_tra(j)),indx_tra(j),j=1,30)
   allocate(travec_norm(ndimf))
   norma = 0
   do j = 1 , ndimf
   norma = norma + coeff_tra(j)
   end do
   norma = sqrt(norma)
   write(6,*) 'norma OF ! D PSIin>', norma 
   travec_norm(:) = travec(:)/norma
   write(102,*) ndimf , travec_norm
   allocate(MAP(kpqf(1,0)+100))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do i = 1 , kpqf(1,0)
   MAP(i) = i
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   count = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   DO i = kpqf(1,0) + 1 , kpqf(1,0) + 100
   count = count + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do j = count , ndimf
   IF( indx_tra(j) .GT. kpqf(1,0) ) then
   MAP(i) = indx_tra(j)
   count = j
   exit
   END IF
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(103,*) kpqf(1,0) + 100 
   write(103,*) MAP
   deallocate(coeff_tra,indx_tra,travec_norm,MAP)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM EXCITED STATE !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*)  'energies and  pseudo-cross-section values OF THE FINAL SPACE STATES, from THE CHOOSEN',statenumber,'EXCITED STATE IN THE INITIAL SPACE'
   allocate(excit(ndimf))
   do i=1,ndimf
   excit(i)= enerf(i)-ener(statenumber)
   end do
   do i = 1 , noutf
   tmvecf(i)=tm(ndimf,arrf(:,i),travec(:))
   osc_strf(i)=2._d/3._d*excit(i)*tmvecf(i)**2
!!!write(6,*) i, enerf(i), os2cs*osc_strf(i)
   write(6,*) i, excit(i), os2cs*osc_strf(i)
   end do
    write(6,*) "ADC2 diagonalization results and trans. mom. FROM THE CHOOSEN EXCITED INITIAL SPACE STATE, in the FINAL space"
!!! call table2(ndimf,nout,enerf(1:nout),arrf(:,1:nout),tmvecf(:),osc_strf(:))
    call table2(ndimf,nout,excit(1:nout),arrf(:,1:nout),tmvecf(:),osc_strf(:))
    write(6,*) ' sums calculated with respect to the',statenumber,'initial excited state'
!!! call get_sigma(ndimf,enerf(:),os2cs*osc_strf(:))
    call get_sigma(ndimf,excit(:),os2cs*osc_strf(:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM EXCITED STATE !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call MP2(E_MP2)
    E_groundstate = Ehf + E_MP2
    write(6,*) 'THE ADC2 GROUND STATE ENERGY AT THIS GEOMETRY IS', E_groundstate




!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!
   deallocate(excit,ener)  
   deallocate(arr,tmvec,osc_str,arrd)
   deallocate(autvec,travec,arrf,enerf,tmvecf,osc_strf)
   deallocate(kpq,kpqf,kpqd)
!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!









ELSE IF  (chrun2.eq.'davi') THEN
write(*,*) 'I WILL PERFORM DAVIDSON ADC2 DIAGONALIZATION IN THE INITIAL SPACE'







   


!!! ALLOCATION PART !!!
!!! ALLOCATION PART !!!
!!! ALLOCATION PART !!!
   allocate(arrf(ndim,ndim),mtm(ndim),tmvec(nout),osc_str(nout))
   allocate(arrd(ndim,ndim),enerf(ndimf),tmvecf(noutf),osc_strf(noutf))
   allocate(autvec(ndim),travec(ndimf))
   allocate(mtmf(ndimf))
!!! ALLOCATION PART !!!
!!! ALLOCATION PART !!!
!!! ALLOCATION PART !!!
   

!!! DAVIDSON IN THE INITIAL SPACE PART !!!
!!! DAVIDSON IN THE INITIAL SPACE PART !!!
!!! DAVIDSON IN THE INITIAL SPACE PART !!!
   inquire(file='hmlt.diai',exist=prs)
   inquire(file='hmlt.offi',exist=prs1)
   inquire(file='fort.111',exist=prs2)
   if(prs .and. prs1 .and. prs2) then
      write(6,*) 'Older hmlt-file will be engaged'
      read(111,*) noffd
   else
      write(6,*) 'Saving complete INITIAL SPACE ADC2 matrix in file'
      call  write_fspace_adc2e_1(ndim,kpq(:,:),noffd,'i') 
      write(111,*) noffd
   end if
   call cpu_time(time)
   write(6,*) 'Time=',time," s"
   allocate(enerdav(davstates),rvec(ndim,davstates))
   allocate(vec_init(ndim))
! DAVIDSON DIAGONALIZATION IN THE INITIAL SPACE
!  call initial_diag (ninit,e_init,vec_init,noffd,nvec,davmem)
   call  master_dav(ndim,noffd,'i')
! Reading Davidson eigenvectors
   call readdavvc(davstates,enerdav,rvec)
!!! DAVIDSON IN THE INITIAL SPACE PART !!!
!!! DAVIDSON IN THE INITIAL SPACE PART !!!
!!! DAVIDSON IN THE INITIAL SPACE PART !!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*) 'Calculating ADC2 transition moments from ground state in ',tranmom2,' direction.'
   write(6,*) 'I AM Calculating ADC2 transition moments in', tranmom2,' direction, i.e. : < PSI0  D',tranmom2,' PSIn  > IN THE INITIAL SPACE'
   call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:))
   write(100,*) ndim , mtm
   do i = 1 , davstates
   tmvec(i) = tm(ndim,rvec(:,i),mtm(:))
   osc_str(i) = 2._d/3._d*enerdav(i)*tmvec(i)**2
   write(6,*) i,enerdav(i),tmvec(i),os2cs*osc_str(i)
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(6,*) "ADC2 diagonalization result in the INITIAL space"
    call table2(ndim,davstates,enerdav(1:davstates),rvec(:,1:davstates),tmvec(:),osc_str(:))
    write(6,*) ' sums calculated with respect to the ground state'
    call get_sigma(ndim,enerdav(:),os2cs*osc_str(:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   deallocate(mtm)




   write(6,*) "DIAGONALIZING IN THE FINAL SYMMETRY  SPACE" 
   call get_fspace_adc2e_direct(ndimf,kpqf(:,:),arrf,enerf)
   write(6,*) 'Calculating ADC2 transition moments from ', statenumber, 'excited state  in ',tranmom,' direction.'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*) 'Calculating ADC2 transition moments from ground state in ',tranmom2,' direction.'
   write(6,*) 'I AM Calculating ADC2 transition moments in', tranmom2,' direction, i.e. : < PSI0  D',tranmom2,' PSIm  > IN THE FINAL SPACE'
   call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtmf(:))
   write(100,*) ndimf , mtmf
   do i = 1 , ndimf
   tmvecf(i) = tm(ndimf,arrf(:,i),mtmf(:))
   osc_strf(i) = 2._d/3._d*enerf(i)*tmvecf(i)**2
   write(6,*) i,enerf(i),tmvecf(i),os2cs*osc_strf(i)
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   nout = 100
   write(6,*) "ADC2 diagonalization result in the FINAL space"
   call table2(ndimf,nout,enerf(1:nout),arrf(:,1:nout),tmvecf(:),osc_strf(:))
   write(6,*) ' sums calculated with respect to the ground state'
   call get_sigma(ndimf,enerf(:),os2cs*osc_strf(:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   deallocate(mtmf)




! here we analyse the eigenvectors and determine the initial state
! here we analyse the eigenvectors and determine the initial state
! here we analyse the eigenvectors and determine the initial state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!! INITIAL STATE VECTOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
write(6,*) statenumber,'selected initial state has been obtained' 
write(6,*) 'The initial state vector is the ', statenumber,'st one'
   do i=1,ndim
      vec_init(i)=rvec(i,statenumber)
!      write(6,*) i,vec_init(i)
   end do
   write(6,*)
   e_init=enerdav(statenumber)
       allocate(coeff(ndim))
       allocate(indx(ndim))
       coeff(:)=vec_init(:)**2
       call dsortindxa1("D",ndim,coeff(:),indx(:))
       write(6,*) statenumber,e_init,e_init*27.211396,(coeff(indx(j)),indx(j),j=1,5)
       write(105,*) e_init
! here we analyse the eigenvectors and determine the initial state
! here we analyse the eigenvectors and determine the initial state
! here we analyse the eigenvectors and determine the initial state

   call cpu_time(time)
   write(6,*) 'Time=',time," s"




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!! MATRIX VECTOR MULTIPLICATION !!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
   write(6,*) "DOING THE MATRIX VECTOR MULTIPLICATION TO OBTAIN CONTRACTING VECTOR"
   if (matvec.eq.1) then
   call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,vec_init,travec)
   write(6,*) "first matrix-vector product done"
   else if (matvec.eq. 2) then
   call get_fspace_adc2_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,vec_init,arrd,travec) 
   write(6,*) "second matrix-vector product done"
   else if (matvec.eq.3) then
   call get_fspace_adc2_DIPOLE_direct(ndim,kpq,vec_init,arrd,travec) 
   write(6,*) "third matrix-vector product done"
   end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TRAVEC VECTOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*)  ' vector to make scalar product with has been calculated'
   write(101,*) ndimf , travec
   allocate(coeff_tra(ndimf))
   allocate(indx_tra(ndimf))
   coeff_tra(:)=travec(:)**2
   call dsortindxa1("D",ndimf,coeff_tra(:),indx_tra(:))
   write(6,*) "BIGGER COEFFICIENTS OF ! D PSIin>"
   write(6,*) (coeff_tra(indx_tra(j)),indx_tra(j),j=1,30)
   allocate(travec_norm(ndimf))
   norma = 0
   do j = 1 , ndimf
   norma = norma + coeff_tra(j)
   end do
   norma = sqrt(norma)
   write(6,*) 'norma OF ! D PSIin>', norma 
   travec_norm(:) = travec(:)/norma
   write(102,*) ndimf , travec_norm
   allocate(MAP(kpqf(1,0)+100))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do i = 1 , kpqf(1,0)
   MAP(i) = i
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   count = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   DO i = kpqf(1,0) + 1 , kpqf(1,0) + 100
   count = count + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do j = count , ndimf
   IF( indx_tra(j) .GT. kpqf(1,0) ) then
   MAP(i) = indx_tra(j)
   count = j
   exit
   END IF
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(103,*) kpqf(1,0) + 100 
   write(103,*) MAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM EXCITED STATE !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*)  'energies and  pseudo-cross-section values OF THE FINAL SPACE STATES, from THE CHOOSEN',statenumber,'EXCITED STATE IN THE INITIAL SPACE'
   allocate(excit(ndimf))
   do i = 1 , ndimf
   excit(i) = enerf(i) - e_init
   end do
   do i = 1 , noutf
   tmvecf(i)=tm(ndimf,arrf(:,i),travec(:))
   osc_strf(i)=2._d/3._d*excit(i)*tmvecf(i)**2
!!!write(6,*) i, enerf(i), os2cs*osc_strf(i)
   write(6,*) i, excit(i), os2cs*osc_strf(i)
   end do
    write(6,*) "ADC2 diagonalization results and trans. mom. FROM THE CHOOSEN EXCITED INITIAL SPACE STATE, in the FINAL space"
!!! call table2(ndimf,noutf,enerf(1:noutf),arrf(:,1:noutf),tmvecf(:),osc_strf(:))
    call table2(ndimf,nout,excit(1:nout),arrf(:,1:nout),tmvecf(:),osc_strf(:))
    write(6,*) ' sums calculated with respect to the',statenumber,'initial excited state'
!!! call get_sigma(ndimf,enerf(:),os2cs*osc_strf(:))
    call get_sigma(ndimf,excit(:),os2cs*osc_strf(:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM EXCITED STATE !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call MP2(E_MP2)
    E_groundstate = Ehf + E_MP2
    write(6,*) 'THE ADC2 GROUND STATE ENERGY AT THIS GEOMETRY IS', E_groundstate


!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!
   deallocate(excit,enerdav,rvec,vec_init)  
   deallocate(arrf,arrd,tmvec,osc_str)
   deallocate(autvec,travec,enerf,tmvecf,osc_strf)
   deallocate(kpq,kpqf,kpqd)
!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!







 END IF











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! LANCZOS DIAGONALIZATION CASE !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF (chrun .eq. 'save') THEN
write(*,*) 'I WILL PERFORM LANCZOS ADC2 DIAGONALIZATION IN THE FINAL SPACE'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! LANCZOS DIAGONALIZATION CASE !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!! DAVIDSON IN THE INITIAL SPACE PART !!!
!!! DAVIDSON IN THE INITIAL SPACE PART !!!
!!! DAVIDSON IN THE INITIAL SPACE PART !!!
   inquire(file='hmlt.diai',exist=prs)
   inquire(file='hmlt.offi',exist=prs1)
   inquire(file='fort.111',exist=prs2)
   if(prs .and. prs1 .and. prs2) then
      write(6,*) 'Older hmlt-file will be engaged'
      read(111,*) noffd
   else
      write(6,*) 'Saving complete INITIAL SPACE ADC2 matrix in file'
      call  write_fspace_adc2e_1(ndim,kpq(:,:),noffd,'i') 
      write(111,*) noffd
   end if
   call cpu_time(time)
   write(6,*) 'Time=',time," s"
   allocate(ener(davstates),rvec(ndim,davstates))
   allocate(vec_init(ndim))
! DAVIDSON DIAGONALIZATION IN THE INITIAL SPACE
!  call initial_diag (ninit,e_init,vec_init,noffd,nvec,davmem)
   call  master_dav(ndim,noffd,'i')
! Reading Davidson eigenvectors
   call readdavvc(davstates,ener,rvec)
!!! DAVIDSON IN THE INITIAL SPACE PART !!!
!!! DAVIDSON IN THE INITIAL SPACE PART !!!
!!! DAVIDSON IN THE INITIAL SPACE PART !!!



!!! ALLOCATION PART !!!
   allocate(mtm(ndim),tmvec(davstates),osc_str(davstates))
!!! ALLOCATION PART !!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(6,*) 'Calculating ADC2 transition moments from ground state in ',tranmom2,' direction.'
   write(6,*) 'I AM Calculating ADC2 transition moments in', tranmom2,' direction, i.e. : < PSI0  D',tranmom2,' PSIn  > IN THE INITIAL SPACE'
   call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:))
   write(100,*) ndim , mtm
   do i = 1 , davstates
   tmvec(i) = tm(ndim,rvec(:,i),mtm(:))
   osc_str(i) = 2._d/3._d * ener(i) * tmvec(i)**2
   write(6,*) i,ener(i),tmvec(i),os2cs*osc_str(i)
   end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM GROUND STATE !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(6,*) "ADC2 diagonalization result in the INITIAL space"
    call table2(ndim,davstates,ener(1:davstates),rvec(:,1:davstates),tmvec(:),osc_str(:))
    write(6,*) ' sums calculated with respect to the ground state'
    call get_sigma(davstates,ener(:),os2cs*osc_str(:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   deallocate(mtm)




! here we analyse the eigenvectors and determine the initial state
! here we analyse the eigenvectors and determine the initial state
! here we analyse the eigenvectors and determine the initial state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!!!!!!!!!!!!!!!!!!!!!!!!! INITIAL STATE VECTOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
write(6,*) statenumber,'selected initial state has been obtained' 
write(6,*) 'The initial state vector is the ', statenumber,'st one'
   do i=1,ndim
      vec_init(i)=rvec(i,statenumber)
!      write(6,*) i,vec_init(i)
   end do
   write(6,*)
   write(6,*) " ALL THE FIRST ENERGIES IN THE INITIAL SPACE SYMMETRY"
   write(6,*) ener
   e_init=ener(statenumber)
   WRITE(6,*) " I AM CHOOSING TO IONIZE FROM THE STATE WITH ENERGY", e_init     
       allocate(coeff(ndim))
       allocate(indx(ndim))
       coeff(:)=vec_init(:)**2
       call dsortindxa1("D",ndim,coeff(:),indx(:))
       write(6,*) statenumber,e_init,e_init*27.211396,(coeff(indx(j)),indx(j),j=1,5)
       write(105,*) e_init
! here we analyse the eigenvectors and determine the initial state
! here we analyse the eigenvectors and determine the initial state
! here we analyse the eigenvectors and determine the initial state

   call cpu_time(time)
   write(6,*) 'Time=',time," s"



!!! ALLOCATION PART !!!
     allocate(travec(ndimf))
!!! ALLOCATION PART !!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(6,*) "DOING THE MATRIX VECTOR MULTIPLICATION TO OBTAIN CONTRACTING VECTOR"
     if (matvec.eq.0) then
     write(6,*) "you did not request to continue with calculations"
     stop
     else if (matvec.eq.1) then
     call get_dipole_initial_product(ndim,ndimf,kpq,kpqf,vec_init,travec)
     write(6,*) "first matrix-vector product done"
     else if (matvec.eq. 2) then
     call get_fspace_adc2_DIPOLE_direct_OK(ndim,ndimf,kpq,kpqf,vec_init,arrd,travec) 
     write(6,*) "second matrix-vector product done"
     else if (matvec.eq.3) then
     call get_fspace_adc2_DIPOLE_direct(ndim,kpq,vec_init,arrd,travec) 
     write(6,*) "third matrix-vector product done"
     end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(101,*) ndimf , travec
     allocate(coeff_tra(ndimf))
     allocate(indx_tra(ndimf))
     coeff_tra(:)=travec(:)**2
     call dsortindxa1("D",ndimf,coeff_tra(:),indx_tra(:))
     write(6,*) "BIGGER COEFFICIENTS OF ! D PSIin>"
     write(6,*) (coeff_tra(indx_tra(j)),indx_tra(j),j=1,30)
     allocate(travec_norm(ndimf))
     norma = 0
     do j = 1 , ndimf
     norma = norma + coeff_tra(j)
     end do
     norma = sqrt(norma)
     write(6,*) 'norma OF ! D PSIin>', norma 
     travec_norm(:) = travec(:)/norma
     write(102,*) ndimf , travec_norm
     write(104,*) indx_tra
     allocate(MAP(kpqf(1,0)+100))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do i = 1 , kpqf(1,0)
     MAP(i) = i
     end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     count = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO i = kpqf(1,0) + 1 , kpqf(1,0) + 100
     count = count + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do j = count , ndimf
     IF( indx_tra(j) .GT. kpqf(1,0) ) then
     MAP(i) = indx_tra(j)
     count = j
     exit
     END IF
     end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(103,*) kpqf(1,0) + 100 
     write(103,*) MAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



     inquire(file='hmlt.diac',exist=prs)
     inquire(file='hmlt.offc',exist=prs1)
     inquire(file='fort.111',exist=prs2)
     if(prs .and. prs1 .and. prs2) then
        write(6,*) 'Older hmlt-file will be engaged'
        read(111,*) noffdf
     else
        write(6,*) 'Saving complete FINAL SPACE ADC2 matrix in file'
        call  write_fspace_adc2e_1(ndimf,kpqf(:,:),noffdf,'c') 
        write(111,*) noffdf
     end if

     call cpu_time(time)
     write(6,*) 'Time=',time," s"



IF ( GO .EQ. 'OK' ) THEN


!!! ALLOCATION PART !!!
   allocate(enerf(ndimf),tmvecf(noutf),mtmf(ndimf))
!!! ALLOCATION PART !!!

     inquire(file='mmnt',exist=prs)
     if (prs) then
        write(6,*) 'Older mmnt-file will be engaged'
        call read_vec(ndim,'mmnt',mtm(:))
     else
        write(6,*) 'Calculating ADC2 transition moments from ground state in ',tranmom2,' direction.'
        write(6,*) 'I AM Calculating ADC2 transition moments in', tranmom2,' direction, i.e. : < PSI0  D',tranmom2,' PSIm  > IN THE FINAL SPACE'
        call get_modifiedtm_adc2(ndimf,kpqf(:,:),mtmf(:))
        call write_vec(ndimf,'mmnt',mtmf(:))
        call write_vec(ndimf,'mmnt',travec(:))
     end if 
     write(100,*) ndimf , mtmf

       write(6,*) 'check', mtmf(1:ndimsf)
       call fill_stvc(ndimsf,mtmf(1:ndimsf))
       write(6,*) 'check', travec(1:ndimsf) 
       call fill_stvc(ndimsf,travec(1:ndimsf))


    mtxidl='full'
    call master_lancdiag(ndimf,noffdf,'c')
    call cpu_time(time)
    write(6,*) 'Time=',time," s"
    allocate(enerf(lancstates),tmvec(lancstates))

!!$ ***lancstates is at most ncyclesXmain, usually we expect it to be lower than that***

    write(6,*) 'Calculating ADC2 transition moments from ground state in ',tranmom2,' direction.'
    write(6,*) 'I AM Calculating ADC2 transition moments in', tranmom2,' direction, i.e. : < PSI0  D',tranmom2,' PSIm  > IN THE FINAL SPACE'
    tmvecf(:) = 0.d0
    call get_tranmom_1(ndimf,lancstates,lancname,mtmf(:),nstates,enerf(:),tmvecf(:),ndimsf)
    allocate(osc_strf(nstates))
    osc_strf(:) = 0.d0
    do i = 1 , nstates
       osc_strf(i) = 2._d/3._d * enerf(i) * tmvecf(i)**2
    end do



!!!! FROM NOW ON IT IS THE SAME BUT WITH travec(:) INSTEAD OF mtm(:) !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM EXCITED STATE !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(6,*)  'energies and  pseudo-cross-section values OF THE FINAL SPACE STATES, from THE CHOOSEN',statenumber,'EXCITED STATE IN THE INITIAL SPACE'
    tmvecf(:) = 0.d0
    call get_tranmom_3(ndimf,lancstates,lancname,travec(:),nstates,enerf(:),tmvecf(:),ndimsf)
    allocate(excit(nstates))
    do i = 1 , nstates
    excit(i) = enerf(i) - e_init
    end do
    osc_strf(:) = 0.d0
    do i = 1 , nstates
    osc_strf(i) = 2._d/3._d * excit(i) * tmvecf(i)**2
    end do
    call cpu_time(time)
    write(6,*) 'Time=',time," s"
    call get_sigma(nstates,excit(1:nstates),os2cs*osc_strf(:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! TRANSITION MOMENTS FROM EXCITED STATE !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    call MP2(E_MP2)
    E_groundstate = Ehf + E_MP2
    write(6,*) 'THE ADC2 GROUND STATE ENERGY AT THIS GEOMETRY IS', E_groundstate



!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!
    deallocate(ener,rvec,vec_init)
    deallocate(enerf,tmvec,osc_str,excit)
    deallocate(travec)          
    deallocate(kpq,kpqf,kpqd)
!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!
!!! DEALLOCATION PART !!!



END IF !!! END IF ON THE GO CONDITION 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 END IF




end subroutine master_adc2ext_prop
