
 module propagate_prepare 

  use constants
  use parameters
  use select_fano
  use misc
  use get_matrix
  use get_moment
  use get_matrix_DIPOLE

  implicit none


  contains





!!! CALCULATE HOW MANY HAMILTONIAN BLOCKS MUST BE PROPAGATED !!!
 subroutine SYM_TOPROP_CALC( DIPOLESYM, ELECTRIC_FIELD, NSYMA_PROP )

 INTEGER, DIMENSION(3), INTENT(IN) :: DIPOLESYM
 REAL*8, DIMENSION(3) , INTENT(IN) :: ELECTRIC_FIELD
 INTEGER, INTENT(OUT) :: NSYMA_PROP

 IF ( POLARIZATION .EQ. 'LIN' ) THEN
 NSYMA_PROP = 2
 ELSE IF ( POLARIZATION .EQ. 'CIR') THEN
 NSYMA_PROP = 4
 END IF

 end subroutine SYM_TOPROP_CALC
!!! CALCULATE HOW MANY HAMILTONIAN BLOCKS MUST BE PROPAGATED !!!




!!! CALCULATE THE MAPPING BETWEEN THE HAMILTONIAN BLOCKS PROPAGATED, AND THE SYMMETRY ORDERED FULL LIST !!!
 subroutine SYM_TOPROP_VECT( DIMEN , NSYMA_PROP , SYM_MAP , DIM_PROP )

 INTEGER, INTENT(IN) :: NSYMA_PROP
 INTEGER, DIMENSION(NSYMA), INTENT(IN) :: DIMEN
 INTEGER, DIMENSION(NSYMA_PROP), INTENT(OUT) :: SYM_MAP
 INTEGER, DIMENSION(NSYMA_PROP), INTENT(OUT) :: DIM_PROP

 IF ( POLARIZATION .EQ. 'LIN' ) THEN
 SYM_MAP(1) = NIRREP
 SYM_MAP(2) = MT( NIRREP , NIRREP2 )
 DIM_PROP(1) = DIMEN(NIRREP)
 DIM_PROP(2) = DIMEN( MT(NIRREP,NIRREP2) )
 ELSE IF ( POLARIZATION .EQ. 'CIR') THEN
 SYM_MAP(1) = NIRREP
 SYM_MAP(2) = 2
 SYM_MAP(3) = 3
 SYM_MAP(4) = 4
 DIM_PROP(1) = DIMEN(NIRREP)
 DIM_PROP(2) = DIMEN(2)
 DIM_PROP(3) = DIMEN(3)
 DIM_PROP(4) = DIMEN(4)
 END IF

 end subroutine SYM_TOPROP_VECT
!!! CALCULATE THE MAPPING BETWEEN THE HAMILTONIAN BLOCKS PROPAGATED, AND THE SYMMETRY ORDERED FULL LIST !!!





    subroutine  SAVE_HAMPIECES_ADC1( numofconfig , NSYMA_PROP , SYM_MAP, DIM_PROP , TOTDIM_PROP, NBUF_SYM)


    INTEGER, DIMENSION(NSYMA), INTENT(IN) :: numofconfig
    INTEGER, INTENT(IN)  :: NSYMA_PROP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: SYM_MAP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: DIM_PROP

    INTEGER, DIMENSION(NSYMA_PROP), INTENT(OUT) :: NBUF_SYM
    INTEGER, INTENT(OUT) :: TOTDIM_PROP

    
    integer :: Simmetrie, spaceSym, part, J, i, K, step, dim, dimsingle, dim_local, J1, J2, IA
    integer :: start
    integer :: UNIT_HAM
    integer :: unit_ham_diag
    integer :: Simmetry
    integer :: dimensione
    

    INTEGER :: nbuf
    INTEGER*8 :: count

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: kpq
    integer, dimension(5) :: config


 
!   TOTDIM_PROP = 0
    TOTDIM_PROP = 1  !!! PARTO DA 1 PERCHE' AGGIUNGO ANCHE IL GROUND STATE !!!
    DO J = 1 , NSYMA_PROP
    TOTDIM_PROP = TOTDIM_PROP + DIM_PROP(J)
    END DO

!!! IL VETTORE DIM_PROP NON TIENE CONTO DEL GROUND STATE E QUINDI IL SUO
!!! ELEMENTO NUMERO 1 E' PIU PICCOLO DI 1 RISPETTO ALLA EFFETTIVA DIMENSIONE 1 CHE
!!! PROPAGHERA'



 NBUF_SYM(:)    = 0





DO Simmetrie = 1 , NSYMA_PROP
   spaceSym = SYM_MAP(Simmetrie)
   Simmetry = spaceSym



!!!***************************************************!!!
IF ( spaceSym .eq. 1 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
dimsingle=kpq(1,0)
dim = kpq(1,0)

else if ( spaceSym .eq. 2 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
dimsingle=kpq(1,0)
dim = kpq(1,0)

else if ( spaceSym .eq. 3 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
dimsingle=kpq(1,0)
dim = kpq(1,0)

else if ( spaceSym .eq. 4 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
dimsingle=kpq(1,0)
dim = kpq(1,0)

else if ( spaceSym .eq. 5 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
dimsingle=kpq(1,0)
dim = kpq(1,0)

else if ( spaceSym .eq. 6 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
dimsingle=kpq(1,0)
dim = kpq(1,0)

else if ( spaceSym .eq. 7 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
dimsingle=kpq(1,0)
dim = kpq(1,0)

else if ( spaceSym .eq. 8 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
dimsingle=kpq(1,0)
dim = kpq(1,0)

end if
!!!***************************************************!!!



IF ( dimsingle .NE. 0 ) THEN


unit_ham_diag = 200 + Simmetrie
OPEN( unit_ham_diag, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')

UNIT_HAM = 100 + Simmetrie
OPEN( UNIT_HAM, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')

!!!***************************************************!!!
 IF ( spaceSym .eq. 1 ) then
!!!***************************************************!!!

 CALL get_diag_tda_save_GS( dim , kpq , unit_ham_diag )

 CALL get_offdiag_tda_save_GS( dim , kpq , nbuf , count , UNIT_HAM )

 NBUF_SYM(Simmetrie) = nbuf

!!!***************************************************!!!
 ELSE
!!!***************************************************!!!

 CALL get_diag_tda_save_OK( dim , kpq , unit_ham_diag )

 CALL get_offdiag_tda_save_OK( dim , kpq , nbuf , count , UNIT_HAM )

 NBUF_SYM(Simmetrie) = nbuf

!!!***************************************************!!!
 END IF
!!!***************************************************!!!

CLOSE( unit_ham_diag )

CLOSE( UNIT_HAM )


END IF




 deallocate(kpq)


!!!***************************************************!!!

END DO




!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!



    end  subroutine SAVE_HAMPIECES_ADC1





    subroutine  SAVE_HAMPIECES_ADC2( numofconfig , NSYMA_PROP , SYM_MAP , DIM_PROP , TOTDIM_PROP , NBUF_SYM )


    INTEGER, DIMENSION(NSYMA), INTENT(IN) :: numofconfig
    INTEGER, INTENT(IN)  :: NSYMA_PROP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: SYM_MAP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: DIM_PROP

    INTEGER, DIMENSION(NSYMA_PROP), INTENT(OUT) :: NBUF_SYM
    INTEGER, INTENT(OUT) :: TOTDIM_PROP

    
    integer :: Simmetrie, spaceSym, part, J, i, K, step, dim, dimsingle, dim_local, J1, J2, IA
    integer :: start
    integer :: UNIT_HAM
    integer :: unit_ham_diag
    integer :: Simmetry
    integer :: dimensione

    

    INTEGER :: nbuf
    INTEGER*8 :: count

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: kpq
    integer, dimension(5) :: config


 
!   TOTDIM_PROP = 0
    TOTDIM_PROP = 1  !!! PARTO DA 1 PERCHE' AGGIUNGO ANCHE IL GROUND STATE !!!
    DO J = 1 , NSYMA_PROP
    TOTDIM_PROP = TOTDIM_PROP + DIM_PROP(J)
    END DO

!!! IL VETTORE DIM_PROP NON TIENE CONTO DEL GROUND STATE E QUINDI IL SUO
!!! ELEMENTO NUMERO 1 E' PIU PICCOLO DI 1 RISPETTO ALLA EFFETTIVA DIMENSIONE 1 CHE
!!! PROPAGHERA'


 NBUF_SYM(:)    = 0





DO Simmetrie = 1 , NSYMA_PROP
   spaceSym = SYM_MAP(Simmetrie)
   Simmetry = spaceSym

!!!***************************************************!!!
IF ( spaceSym .eq. 1 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 2 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 3 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 4 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 5 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 6 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 7 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 8 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

end if
!!!***************************************************!!!











IF ( dimsingle .NE. 0 ) THEN


unit_ham_diag = 200 + Simmetrie
OPEN( unit_ham_diag, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')

UNIT_HAM = 100 + Simmetrie
OPEN( UNIT_HAM, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')

!!!***************************************************!!!
 IF ( spaceSym .eq. 1 ) then
!!!***************************************************!!!

 CALL get_diag_adc2_save_GS(dimsingle,dim-dimsingle,kpq, unit_ham_diag )

 CALL get_offdiag_adc2_save_GS(dim,kpq,nbuf,count, UNIT_HAM )

 NBUF_SYM(Simmetrie) = nbuf

!!!***************************************************!!!
 ELSE
!!!***************************************************!!!

 CALL get_diag_adc2_save_OK(dimsingle,dim-dimsingle,kpq, unit_ham_diag )

 CALL get_offdiag_adc2_save_OK(dim,kpq,nbuf,count, UNIT_HAM )

 NBUF_SYM(Simmetrie) = nbuf

!!!***************************************************!!!
 END IF
!!!***************************************************!!!

CLOSE( unit_ham_diag )
CLOSE( UNIT_HAM )


END IF




 deallocate(kpq)


!!!***************************************************!!!

END DO




!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!



    end  subroutine SAVE_HAMPIECES_ADC2






    subroutine  SAVE_HAMPIECES_ADC2EXT( numofconfig , NSYMA_PROP , SYM_MAP , DIM_PROP , TOTDIM_PROP , NBUF_SYM )


    INTEGER, DIMENSION(NSYMA), INTENT(IN) :: numofconfig
    INTEGER, INTENT(IN)  :: NSYMA_PROP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: SYM_MAP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: DIM_PROP

    INTEGER, DIMENSION(NSYMA_PROP), INTENT(OUT) :: NBUF_SYM
    INTEGER, INTENT(OUT) :: TOTDIM_PROP

    
    integer :: Simmetrie, spaceSym, part, J, i, K, step, dim, dimsingle, dim_local, J1, J2, IA
    integer :: start
    integer :: UNIT_HAM
    integer :: unit_ham_diag
    integer :: Simmetry
    integer :: dimensione

    

    INTEGER :: nbuf
    INTEGER*8 :: count

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: kpq
    integer, dimension(5) :: config


 
!   TOTDIM_PROP = 0
    TOTDIM_PROP = 1  !!! PARTO DA 1 PERCHE' AGGIUNGO ANCHE IL GROUND STATE !!!
    DO J = 1 , NSYMA_PROP
    TOTDIM_PROP = TOTDIM_PROP + DIM_PROP(J)
    END DO

!!! IL VETTORE DIM_PROP NON TIENE CONTO DEL GROUND STATE E QUINDI IL SUO
!!! ELEMENTO NUMERO 1 E' PIU PICCOLO DI 1 RISPETTO ALLA EFFETTIVA DIMENSIONE 1 CHE
!!! PROPAGHERA'


 NBUF_SYM(:)    = 0





DO Simmetrie = 1 , NSYMA_PROP
   spaceSym = SYM_MAP(Simmetrie)
   Simmetry = spaceSym

!!!***************************************************!!!
IF ( spaceSym .eq. 1 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 2 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 3 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 4 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 5 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 6 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 7 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

else if ( spaceSym .eq. 8 ) then

dim = numofconfig(spaceSym)
    allocate(kpq(7,0:nBas**2*nOcc**2))
  CALL  select_atom_is_ALL( kpq , Simmetry , dimensione )
  CALL  select_atom_d_ALL( kpq ,-1 , Simmetry , dimensione )
dimsingle=kpq(1,0)

end if
!!!***************************************************!!!




IF ( dimsingle .NE. 0 ) THEN


unit_ham_diag = 200 + Simmetrie
OPEN( unit_ham_diag, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')

UNIT_HAM = 100 + Simmetrie
OPEN( UNIT_HAM, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')

!!!***************************************************!!!
 IF ( spaceSym .eq. 1 ) then
!!!***************************************************!!!

 CALL get_diag_adc2ext_save_GS(dimsingle,dim-dimsingle,kpq, unit_ham_diag )

 CALL get_offdiag_adc2ext_save_GS(dim,kpq,nbuf,count, UNIT_HAM )

 NBUF_SYM(Simmetrie) = nbuf

!!!***************************************************!!!
 ELSE
!!!***************************************************!!!

 CALL get_diag_adc2ext_save_OK(dimsingle,dim-dimsingle,kpq, unit_ham_diag )

 CALL get_offdiag_adc2ext_save_OK(dim,kpq,nbuf,count, UNIT_HAM )

 NBUF_SYM(Simmetrie) = nbuf

!!!***************************************************!!!
 END IF
!!!***************************************************!!!

CLOSE( unit_ham_diag )
CLOSE( UNIT_HAM )


END IF




 deallocate(kpq)


!!!***************************************************!!!

END DO



!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!
!!!*************************************************************************!!!



    end  subroutine SAVE_HAMPIECES_ADC2EXT

















!!! CALCOLA QUANTI BLOCCHI DI DIPOLO AVRO DA PROPAGARE !!! 
 SUBROUTINE SYM_TOPROP_CALC_DIP2( DIPOLE, ELECTRIC_FIELD, NSYMA_PROP , SYM_MAP , NBLOCKS_DIP )

 INTEGER, DIMENSION(3), INTENT(IN) :: DIPOLE
 REAL*8, DIMENSION(3) , INTENT(IN) :: ELECTRIC_FIELD
 INTEGER, INTENT(IN) :: NSYMA_PROP
 INTEGER, DIMENSION(NSYMA_PROP) , INTENT(IN) :: SYM_MAP

 INTEGER, INTENT(OUT) :: NBLOCKS_DIP


 INTEGER :: KLPD , NUMLOC , IA , IB , BLOCK , DIMSINGLE , MT1 , COORD , DIP , STEP , J1 , J2 , DIM_LOCAL
 character(2) :: flag


    NBLOCKS_DIP = 0

 DO IA = 1 , NSYMA_PROP
 DO IB = IA , NSYMA_PROP
    MT1 = MT(SYM_MAP(IA),SYM_MAP(IB))

    flag = 'NO'

    DO COORD = 1 , 3
    DIP = DIPOLE(COORD)
    IF ( MT(MT1,DIP) .EQ. 1 ) THEN
    IF ( ELECTRIC_FIELD(COORD) .NE. 0.d0 ) THEN
    flag = 'si'  
    END IF
    END IF
    END DO

 IF ( flag .eq. 'si') THEN
 NBLOCKS_DIP = NBLOCKS_DIP + 1
 END IF

 END DO
 END DO


 end subroutine SYM_TOPROP_CALC_DIP2
!!! CALCOLA QUANTI BLOCCHI DI DIPOLO AVRO DA PROPAGARE !!! 

 
!!! SCRIVE I PEZZI DIPOLO SU FILE ( IL PEZZO [KLPD], NEL FILE [100 + NSYMA_PROP + KLPD] ) !!!
    subroutine SAVE_DIPOLE_PIECES_ADC1(numofconfig,NSYMA_PROP,DIM_PROP,SYM_MAP,DIPOLE_SYM,ELECTRIC_FIELD ,NBLOCKS_DIP,NREC_VEC)


    INTEGER, DIMENSION(NSYMA), INTENT(IN) :: numofconfig
    INTEGER, INTENT(IN)  :: NSYMA_PROP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: DIM_PROP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: SYM_MAP
    INTEGER, DIMENSION(3), INTENT(IN) :: DIPOLE_SYM
    REAL*8, DIMENSION(3), INTENT(IN) :: ELECTRIC_FIELD

    INTEGER, INTENT(IN)  :: NBLOCKS_DIP
    INTEGER, DIMENSION(NBLOCKS_DIP), INTENT(INOUT) :: NREC_VEC

    
    integer :: part, J, i, K, step, dimA , dimB , dimsingleA , dimsingleB , DIMSINGLE , dim_local, J1, J2, IA , IB
    integer :: start, MT1, DIP, COORD, KLPD, dimrow
    integer :: UNIT_DIPOLE
    integer :: unit_dipole_diag
     integer :: dimensione   
 
    real*8, dimension(:), allocatable :: ar_offdiag
    integer :: nbuf
    integer*8 :: count  
 
    INTEGER :: IA_VERA , IB_VERA
    
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: kpqA , kpqB

    REAL*8, DIMENSION(:), ALLOCATABLE :: mtm, COLONNA_DIPOLE
    real*8, dimension(:,:) , allocatable :: DIPOLE
    integer :: dimSinBound 
    INTEGER :: NBLOCKS_DIP_LOC

    REAL*8 :: GROUND_DIP

 


!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
    KLPD = 0
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!



 DO IA = 1 , NSYMA_PROP
 IA_VERA = SYM_MAP(IA)
 dimrow = DIM_PROP(IA)
 dimA = numofconfig( IA_VERA )
 allocate(kpqA(7,0:nBas**2*nOcc**2))
 CALL  select_atom_is_ALL( kpqA , IA_VERA , dimensione )


dimA = kpqA(1,0)


  DO IB = IA , NSYMA_PROP
  IB_VERA = SYM_MAP(IB)


DIMSINGLE = DIM_PROP(IB)
dimB = numofconfig( IB_VERA )
allocate(kpqB(7,0:nBas**2*nOcc**2))
CALL  select_atom_is_ALL( kpqB , IB_VERA , dimensione )


dimB = kpqB(1,0)    

MT1 = MT( IA_VERA , IB_VERA )

   DO COORD = 1 , 3
    DIP = DIPOLESYM(COORD)
    IF ( MT(MT1,DIP) .EQ. 1 ) THEN
    IF ( ELECTRIC_FIELD(COORD) .NE. 0.d0 ) THEN


  if ( COORD .eq. 1 ) then
     dpl(:,:)=x_dipole(:,:)
  elseif ( COORD .eq. 2 ) then
     dpl(:,:)=y_dipole(:,:)
  elseif ( COORD .eq. 3 ) then
     dpl(:,:)=z_dipole(:,:)
  end if

    CHECK_dip = DIPOLESYM(COORD)

!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
    KLPD = KLPD + 1
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!


 IF ( IA .EQ. IB ) THEN
!!!************* CALCULATE THE OTHER CONTRIBUTIONS ************!!! 
! call   get_diag_ADC1_DIP(GROUND_DIP)
! WRITE(*,*) 'SYMMETRY , GROUND_DIP' , IA , GROUND_DIP
!!!************************************************************!!!
 END IF





unit_dipole_diag = 200 + NSYMA_PROP + IA
OPEN( unit_dipole_diag , STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')


!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
UNIT_DIPOLE = 100 + NSYMA_PROP + KLPD
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
OPEN( UNIT_DIPOLE, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')


!!!************************************************************!!!
IF ( ( IA .EQ. IB ) .AND. ( IA .EQ. 1 ) ) THEN
!!!************************************************************!!!


  CALL  get_diag_tda_DIPOLE_OK_SAME_GS(dimA,kpqA,ar_offdiag, unit_dipole_diag )

  CALL  get_offdiag_tda_DIPOLE_SAVE_OK_SAME_GS(dimA,kpqA,nbuf,count, UNIT_DIPOLE )

   NREC_VEC( KLPD ) = nbuf


!!!************************************************************!!!
ELSE  IF ( IA .EQ. 1 ) THEN
!!!************************************************************!!!


 CALL  get_offdiag_tda_DIPOLE_SAVE_OK_GS( dimA , dimB , kpqA , kpqB , nbuf , count , UNIT_DIPOLE )

   NREC_VEC( KLPD ) = nbuf


!!!************************************************************!!!
ELSE IF ( IA .EQ. IB ) THEN
!!!************************************************************!!!


  CALL get_diag_tda_DIPOLE_OK_SAME(dimA,kpqA,ar_offdiag, unit_dipole_diag )

  CALL get_offdiag_tda_DIPOLE_SAVE_OK_SAME(dimA,kpqA,nbuf,count, UNIT_DIPOLE )

   NREC_VEC( KLPD ) = nbuf


!!!************************************************************!!!
ELSE
!!!************************************************************!!!


  CALL get_offdiag_tda_DIPOLE_SAVE_OK(dimA,dimB,kpqA,kpqB,nbuf,count, UNIT_DIPOLE )

   NREC_VEC( KLPD ) = nbuf


!!!************************************************************!!!
END IF
!!!************************************************************!!!

CLOSE( UNIT_DIPOLE )
CLOSE( unit_dipole_diag )


    END IF
    END IF
   END DO
  deallocate(kpqB)
  END DO
 deallocate(kpqA)
 END DO

IF ( KLPD .NE. NBLOCKS_DIP ) THEN
 WRITE(*,*) 'ERROR IN DIPOLE BLOCKS NUMBER'
 STOP
END IF


    end  subroutine SAVE_DIPOLE_PIECES_ADC1
!!! SCRIVE I PEZZI DIPOLO SU FILE ( IL PEZZO [KLPD], NEL FILE [100 + NSYMA_PROP + KLPD] ) !!!





!!! SCRIVE I PEZZI DIPOLO SU FILE ( IL PEZZO [KLPD], NEL FILE [100 + NSYMA_PROP + KLPD] ) !!!
    subroutine  SAVE_DIPOLE_PIECES_ADC2( numofconfig , NSYMA_PROP , DIM_PROP , SYM_MAP , DIPOLESYM , ELECTRIC_FIELD , NBLOCKS_DIP , NREC_VEC )


    INTEGER, DIMENSION(NSYMA), INTENT(IN) :: numofconfig
    INTEGER, INTENT(IN)  :: NSYMA_PROP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: DIM_PROP
    INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN) :: SYM_MAP
    INTEGER, DIMENSION(3), INTENT(IN) :: DIPOLESYM
    REAL*8, DIMENSION(3), INTENT(IN) :: ELECTRIC_FIELD

    INTEGER, INTENT(IN)  :: NBLOCKS_DIP
    INTEGER, DIMENSION(NBLOCKS_DIP), INTENT(INOUT) :: NREC_VEC


    real*8, dimension(:), allocatable :: ar_offdiag

    
    integer :: part, J, i, K, step, dimA , dimB , dimsingleA , dimsingleB , DIMSINGLE , dim_local, J1, J2, IA , IB
    integer :: start, MT1, DIP, COORD, KLPD, dimrow

    integer :: UNIT_DIPOLE
    integer :: unit_dipole_diag
    integer :: nbuf
    integer*8 :: count  

     integer :: dimensione   
    
    INTEGER :: IA_VERA , IB_VERA
    
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: kpqA , kpqB

    REAL*8, DIMENSION(:), ALLOCATABLE :: mtm, COLONNA_DIPOLE
    real*8, dimension(:,:) , allocatable :: DIPOLE
    integer :: dimSinBound 
    INTEGER :: NBLOCKS_DIP_LOC

    REAL*8 :: GROUND_DIP

 

!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
    KLPD = 0
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!




 DO IA = 1 , NSYMA_PROP
 IA_VERA = SYM_MAP(IA)
 dimrow = DIM_PROP(IA)
 dimA = numofconfig( IA_VERA )
 allocate(kpqA(7,0:nBas**2*nOcc**2))
 CALL  select_atom_is_ALL( kpqA , IA_VERA , dimensione )
 CALL  select_atom_d_ALL( kpqA ,-1 , IA_VERA , dimensione )


dimsingleA=kpqA(1,0)


  DO IB = IA , NSYMA_PROP
  IB_VERA = SYM_MAP(IB)


DIMSINGLE = DIM_PROP(IB)
dimB = numofconfig( IB_VERA )
allocate(kpqB(7,0:nBas**2*nOcc**2))
CALL  select_atom_is_ALL( kpqB , IB_VERA , dimensione )
CALL  select_atom_d_ALL( kpqB ,-1 , IB_VERA , dimensione )

dimsingleB=kpqB(1,0)    

MT1 = MT( IA_VERA , IB_VERA )

   DO COORD = 1 , 3
    DIP = DIPOLESYM(COORD)
    IF ( MT(MT1,DIP) .EQ. 1 ) THEN
    IF ( ELECTRIC_FIELD(COORD) .NE. 0.d0 ) THEN



  if ( COORD .eq. 1 ) then
     dpl(:,:)=x_dipole(:,:)
  elseif ( COORD .eq. 2 ) then
     dpl(:,:)=y_dipole(:,:)
  elseif ( COORD .eq. 3 ) then
     dpl(:,:)=z_dipole(:,:)
  end if

    CHECK_dip = DIPOLESYM(COORD)


!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
    KLPD = KLPD + 1
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!


 IF ( IA .EQ. IB ) THEN
!!!************* CALCULATE THE OTHER CONTRIBUTIONS ************!!! 
! call   get_diag_ADC1_DIP(GROUND_DIP)
! WRITE(*,*) 'SYMMETRY , GROUND_DIP' , IA , GROUND_DIP
!!!************************************************************!!!
 END IF





unit_dipole_diag = 200 + NSYMA_PROP + IA
OPEN( unit_dipole_diag , STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')


!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
UNIT_DIPOLE = 100 + NSYMA_PROP + KLPD
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
!!! INDICE CHE CONTA I BLOCCHI DA SALVARE !!!
OPEN( UNIT_DIPOLE, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')


!!!************************************************************!!!
IF ( ( IA .EQ. IB ) .AND. ( IA .EQ. 1 ) ) THEN
!!!************************************************************!!!


 CALL  get_diag_adc2_DIPOLE_OK_SAME_GS(dimA,kpqA,ar_offdiag, unit_dipole_diag )

  CALL get_offdiag_adc2_DIPOLE_SAVE_OK_SAME_GS(dimA,kpqA,nbuf,count, UNIT_DIPOLE )

   NREC_VEC( KLPD ) = nbuf


!!!************************************************************!!!
ELSE  IF ( IA .EQ. 1 ) THEN
!!!************************************************************!!!


  CALL get_offdiag_adc2_DIPOLE_SAVE_OK_GS(dimA,dimB,kpqA,kpqB,nbuf,count, UNIT_DIPOLE )

   NREC_VEC( KLPD ) = nbuf



!!!************************************************************!!!
ELSE IF ( IA .EQ. IB ) THEN
!!!************************************************************!!!


  CALL get_diag_adc2_DIPOLE_OK_SAME(dimA,kpqA,ar_offdiag, unit_dipole_diag )

  CALL get_offdiag_adc2_DIPOLE_SAVE_OK_SAME(dimA,kpqA,nbuf,count, UNIT_DIPOLE )

   NREC_VEC( KLPD ) = nbuf


!!!************************************************************!!!
ELSE
!!!************************************************************!!!


  CALL get_offdiag_adc2_DIPOLE_SAVE_OK(dimA,dimB,kpqA,kpqB,nbuf,count, UNIT_DIPOLE )

   NREC_VEC( KLPD ) = nbuf



!!!************************************************************!!!
END IF
!!!************************************************************!!!

CLOSE( UNIT_DIPOLE )
CLOSE( unit_dipole_diag )


    END IF
    END IF
   END DO
  deallocate(kpqB)
  END DO
 deallocate(kpqA)
 END DO

IF ( KLPD .NE. NBLOCKS_DIP ) THEN
 WRITE(*,*) 'ERROR IN DIPOLE BLOCKS NUMBER'
 STOP
END IF


    end  subroutine SAVE_DIPOLE_PIECES_ADC2
!!! SCRIVE I PEZZI DIPOLO SU FILE ( IL PEZZO [KLPD], NEL FILE [100 + NSYMA_PROP + KLPD] ) !!!






















!!! CALCOLA QUANTI PEZZI DI DIPOLO AVRO DA PROPAGARE, IL NUMERO DI PROCESSORI NECESSARI PER I DIPOLI !!! 
!!! CALCOLA QUANTI PEZZI DI DIPOLO AVRO DA PROPAGARE, IL NUMERO DI PROCESSORI NECESSARI PER I DIPOLI !!! 
!!! CALCOLA QUANTI PEZZI DI DIPOLO AVRO DA PROPAGARE, IL NUMERO DI PROCESSORI NECESSARI PER I DIPOLI !!! 
 SUBROUTINE SYM_TOPROP_CALC_DIP( DIPOLE, ELECTRIC_FIELD, NSYMA_PROP , DIM_PROP , SYM_MAP , NDIV_DIP , NDIV , NBLOCKS_DIP , NBUF_SYM , NREC_VECTOR , KLPDTOT , HAM_PIECES )

 INTEGER, DIMENSION(3),           INTENT(IN)    :: DIPOLE
 REAL*8, DIMENSION(3) ,           INTENT(IN)    :: ELECTRIC_FIELD
 INTEGER,                         INTENT(IN)    :: NSYMA_PROP
 INTEGER, DIMENSION(NSYMA_PROP) , INTENT(IN)    :: DIM_PROP
 INTEGER, DIMENSION(NSYMA_PROP) , INTENT(IN)    :: SYM_MAP
 INTEGER,                         INTENT(IN)    :: NBLOCKS_DIP
 INTEGER, DIMENSION(NBLOCKS_DIP), INTENT(IN)    :: NDIV_DIP
 INTEGER, DIMENSION(NSYMA_PROP) , INTENT(INOUT) :: NDIV
 INTEGER, DIMENSION(NSYMA_PROP) , INTENT(IN)    :: NBUF_SYM
 INTEGER, DIMENSION(NBLOCKS_DIP), INTENT(IN)    :: NREC_VECTOR

 INTEGER, INTENT(OUT) :: KLPDTOT
 INTEGER, INTENT(OUT) :: HAM_PIECES

 INTEGER :: KLPD , NUMLOC , IA , IB , BLOCK , DIMSINGLE , MT1 , COORD , DIP , STEP , J1 , J2 , DIM_LOCAL
 INTEGER :: PART
 INTEGER :: NBLOCKS_DIP_LOC

  INTEGER :: DIV_DIP 

  DIV_DIP = 1
  NDIV(:) = 1

    KLPD = 0
!!! ****************************************************** !!!
    BLOCK = 0
!!! ****************************************************** !!!

 DO IA = 1 , NSYMA_PROP
  DO IB = IA , NSYMA_PROP
     MT1 = MT(SYM_MAP(IA),SYM_MAP(IB))
     DIMSINGLE = DIM_PROP(IB)
     DO COORD = 1 , 3
     DIP = DIPOLE(COORD)
     IF ( MT(MT1,DIP) .EQ. 1 ) THEN
     IF ( ELECTRIC_FIELD(COORD) .NE. 0.d0 ) THEN
!!! ****************************************************** !!!
     BLOCK = BLOCK + 1
!!! ****************************************************** !!!
     NUMLOC = 0
!!!    STEP  = ( DIMSINGLE + NDIV_DIP(BLOCK) - 1 )/NDIV_DIP(BLOCK)
!!!    DO PART = 1 , NDIV_DIP(BLOCK)
     STEP  = ( NREC_VECTOR(BLOCK) + DIV_DIP - 1 )/DIV_DIP
     DO PART = 1 , DIV_DIP
     J1 = (PART-1)*(STEP) + 1  !!! COLONNA DI PARTENZA
     J2 = MIN ( J1 + STEP -1 , NREC_VECTOR(BLOCK) )  !!! COLONNA DI ARRIVO
     DIM_LOCAL = J2 -J1 + 1  
     NUMLOC = NUMLOC + 1
     KLPD = KLPD + 1
     WRITE(*,*) 'KLPD,IA,IB,DIP,PART',KLPD,IA,IB,DIP,PART
     WRITE(*,*) 'KLPD,VERA_IA,VERA_IB,DIP,PART',KLPD,SYM_MAP(IA),SYM_MAP(IB),DIP,PART
!    NLPD(IA) = NLPD(IA) + 1
!    LPSYM(KLPD) = IB
!    MU1D(KLPD) = IA
!    MU2D(KLPD) = IB
!    KAPPAD(KLPD) = COORD
!    COLONNE_DIPOLO(KLPD) = DIM_LOCAL
!    NUM(KLPD) = NUMLOC

    END DO
    END IF
    END IF
   END DO
  END DO
 END DO


 KLPDTOT = KLPD
WRITE(*,*) 'KLPD, KLPDTOT', KLPD , KLPDTOT







   HAM_PIECES = 0
DO IA = 1 , NSYMA_PROP
   DIMSINGLE = DIM_PROP(IA)
   STEP  = ( NBUF_SYM(IA) + NDIV(IA) - 1 )/NDIV(IA)
   NUMLOC = 0
DO PART = 1 , NDIV(IA)  !!!  --------------------------------------------------------------------------------->
   NUMLOC = NUMLOC + 1
   HAM_PIECES = HAM_PIECES + 1
   J1 = (PART-1)*(STEP) + 1  !!! COLONNA DI PARTENZA
   J2 = MIN ( J1 + STEP -1 , NBUF_SYM(IA) )  !!! COLONNA DI ARRIVO
!   COLONNE_HAM(HAM_PIECES_LOC) = J2 - J1 + 1  
!   SIMMETRIA(HAM_PIECES_LOC) = IA
!   NUM_DIAG(HAM_PIECES_LOC) = NUMLOC
END DO !!!  <---------------------------------------------------------------------------------

END DO


WRITE(*,*) 'NSYMA_PROP , HAM_PIECES , KLPDTOT'
WRITE(*,*) NSYMA_PROP , HAM_PIECES , KLPDTOT










 end subroutine SYM_TOPROP_CALC_DIP
!!! CALCOLA QUANTI PEZZI DI DIPOLO AVRO DA PROPAGARE, IL NUMERO DI PROCESSORI NECESSARI PER I DIPOLI !!! 
!!! CALCOLA QUANTI PEZZI DI DIPOLO AVRO DA PROPAGARE, IL NUMERO DI PROCESSORI NECESSARI PER I DIPOLI !!! 
!!! CALCOLA QUANTI PEZZI DI DIPOLO AVRO DA PROPAGARE, IL NUMERO DI PROCESSORI NECESSARI PER I DIPOLI !!! 




















!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! TUTTE LE QUANTITA TENGONO CONTO DEL GROUND STATE IN SYMMETRY 1 !!! 
!!! TUTTE LE QUANTITA TENGONO CONTO DEL GROUND STATE IN SYMMETRY 1 !!! 
!!! TUTTE LE QUANTITA TENGONO CONTO DEL GROUND STATE IN SYMMETRY 1 !!! 
subroutine SYM_TOPROP_VECT_DIP(NSYMA_PROP,SYM_MAP,DIM_PROP,NDIV,DIPOLE,ELECTRIC_FIELD,NBLOCKS_DIP,NDIV_DIP,HAM_PIECES,KLPDTOT,SIMMETRIA,NLPD,LPSYM,MU1D,MU2D,KAPPAD,DIPOLE_BLOCK,NUM,NUM_DIAG,NBUF_SYM,NREC_VECTOR,NRECTOT_VECT,NREC_VECTOR_BIS, RECINI_VECT )

 INTEGER, INTENT(IN)                                     :: NSYMA_PROP
 INTEGER, INTENT(IN)                                     :: NBLOCKS_DIP
 INTEGER, INTENT(IN)                                     :: HAM_PIECES
 INTEGER, INTENT(IN)                                     :: KLPDTOT
 INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN)              :: SYM_MAP
 INTEGER, DIMENSION(NSYMA_PROP), INTENT(INOUT)           :: DIM_PROP
 INTEGER, DIMENSION(NSYMA_PROP), INTENT(INOUT)           :: NDIV
 INTEGER, DIMENSION(3), INTENT(IN)                       :: DIPOLE
 REAL*8, DIMENSION(3) , INTENT(IN)                       :: ELECTRIC_FIELD
 INTEGER, DIMENSION(NBLOCKS_DIP), INTENT(IN)             :: NDIV_DIP
 INTEGER, DIMENSION(NSYMA_PROP), INTENT(IN)              :: NBUF_SYM
 INTEGER, DIMENSION(NBLOCKS_DIP), INTENT(IN)             :: NREC_VECTOR
 INTEGER, DIMENSION( NSYMA_PROP + NBLOCKS_DIP ), INTENT(OUT) :: NRECTOT_VECT
 INTEGER, DIMENSION( HAM_PIECES + KLPDTOT ), INTENT(INOUT)   :: NREC_VECTOR_BIS
 INTEGER, DIMENSION( HAM_PIECES + KLPDTOT ), INTENT(INOUT)   :: RECINI_VECT

 INTEGER, DIMENSION(HAM_PIECES), INTENT(OUT) :: SIMMETRIA
 INTEGER, DIMENSION(NSYMA_PROP), INTENT(OUT) :: NLPD
 INTEGER, DIMENSION(KLPDTOT), INTENT(OUT)    :: LPSYM
 INTEGER, DIMENSION(KLPDTOT), INTENT(OUT)    :: MU1D
 INTEGER, DIMENSION(KLPDTOT), INTENT(OUT)    :: MU2D
 INTEGER, DIMENSION(KLPDTOT), INTENT(OUT)    :: KAPPAD
 INTEGER, DIMENSION(KLPDTOT), INTENT(OUT)    :: DIPOLE_BLOCK
 INTEGER, DIMENSION(KLPDTOT), INTENT(OUT)    :: NUM
 INTEGER, DIMENSION(HAM_PIECES), INTENT(OUT) :: NUM_DIAG

 INTEGER :: KLPD , NUMLOC , IA , IB , BLOCK , DIMSINGLE , MT1 , COORD , DIP , STEP , J1 , J2 , DIM_LOCAL , PART
 INTEGER :: HAM_PIECES_LOC , I
 INTEGER :: NBLOCKS_DIP_LOC
  INTEGER :: DIV_DIP 

    NREC_VECTOR_BIS(:) = 0
    RECINI_VECT(:) = 0
    NRECTOT_VECT(:) = 0



    KLPD = 0
    NLPD(:) = 0
!!! ****************************************************** !!!
    BLOCK = 0
!!! ****************************************************** !!!
   
    DIV_DIP = 1
    NDIV(:) = 1

 DO IA = 1 , NSYMA_PROP
  DO IB = IA , NSYMA_PROP
     DIMSINGLE = DIM_PROP(IB)
     MT1 = MT(SYM_MAP(IA),SYM_MAP(IB))
     DO COORD = 1 , 3
     DIP = DIPOLE(COORD)
     IF ( MT(MT1,DIP) .EQ. 1 ) THEN
     IF ( ELECTRIC_FIELD(COORD) .NE. 0.d0 ) THEN
!!! ****************************************************** !!!
     BLOCK = BLOCK + 1
!!! ****************************************************** !!!
     NUMLOC = 0
!!!    STEP  = ( DIMSINGLE + NDIV_DIP(BLOCK) - 1 )/NDIV_DIP(BLOCK)
!!!    DO PART = 1 , NDIV_DIP(BLOCK)
!!! ****************************************************** !!!
     NRECTOT_VECT( NSYMA_PROP + BLOCK ) = NREC_VECTOR( BLOCK )
!!! ****************************************************** !!!
     STEP  = ( NREC_VECTOR(BLOCK) + DIV_DIP - 1 )/DIV_DIP
     DO PART = 1 , DIV_DIP
     J1 = (PART-1)*(STEP) + 1  !!! COLONNA DI PARTENZA
     J2 = MIN ( J1 + STEP -1 , NREC_VECTOR(BLOCK) )  !!! COLONNA DI ARRIVO
     DIM_LOCAL = J2 - J1 + 1  
     NUMLOC = NUMLOC + 1

!!! ****************************************************** !!!
     KLPD = KLPD + 1
!!! ****************************************************** !!!

     NLPD(IA) = NLPD(IA) + 1
     LPSYM(KLPD) = IB
     MU1D(KLPD) = IA
     MU2D(KLPD) = IB
     KAPPAD(KLPD) = COORD
     DIPOLE_BLOCK(KLPD) = BLOCK

!!! ****************************************************** !!!
     RECINI_VECT( HAM_PIECES + KLPD ) = J1
!!! ****************************************************** !!!

!!! ****************************************************** !!!
     NREC_VECTOR_BIS( HAM_PIECES + KLPD ) = DIM_LOCAL
!!! ****************************************************** !!!

     NUM(KLPD) = NUMLOC
     END DO
     END IF
     END IF
   END DO
  END DO
 END DO


WRITE(*,*) 'KLPD, KLPDTOT', KLPD , KLPDTOT
IF ( KLPD .NE. KLPDTOT ) THEN
 WRITE(*,*) 'ERROR IN DIPOLE BLOCKS NUMBER'
 STOP
END IF


!!! ****************************************************** !!!
   HAM_PIECES_LOC = 0
!!! ****************************************************** !!!

DO IA = 1 , NSYMA_PROP
   DIMSINGLE = DIM_PROP(IA)
!!! ****************************************************** !!!
   NRECTOT_VECT( IA ) = NBUF_SYM( IA )
!!! ****************************************************** !!!
   STEP  = ( NBUF_SYM(IA) + NDIV(IA) - 1 )/NDIV(IA)
   NUMLOC = 0
DO PART = 1 , NDIV(IA)  !!!  --------------------------------------------------------------------------------->
   NUMLOC = NUMLOC + 1

!!! ****************************************************** !!!
   HAM_PIECES_LOC = HAM_PIECES_LOC + 1
!!! ****************************************************** !!!

   J1 = (PART-1)*(STEP) + 1  !!! COLONNA DI PARTENZA
   J2 = MIN ( J1 + STEP -1 , NBUF_SYM(IA) )  !!! COLONNA DI ARRIVO

!!! ****************************************************** !!!
   RECINI_VECT( HAM_PIECES_LOC ) = J1   
!!! ****************************************************** !!!

!!! ****************************************************** !!!
   NREC_VECTOR_BIS( HAM_PIECES_LOC ) = J2 - J1 + 1  
!!! ****************************************************** !!!

   SIMMETRIA(HAM_PIECES_LOC) = IA
   NUM_DIAG(HAM_PIECES_LOC) = NUMLOC
END DO !!!  <---------------------------------------------------------------------------------

END DO


DIM_PROP(1) = DIM_PROP(1) + 1 !!! GROUND STATE INCLUSION CORRECTION



OPEN( 92, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='FORMATTED')
WRITE(92,*) 'NSYMA_PROP , NBLOCKS_DIP , HAM_PIECES , KLPDTOT'
WRITE(92,*)  NSYMA_PROP , NBLOCKS_DIP , HAM_PIECES , KLPDTOT
WRITE(92,*) 'DIM_PROP'
WRITE(92,*)  DIM_PROP
WRITE(92,*) 'NRECTOT_VECT'
WRITE(92,*)  NRECTOT_VECT
WRITE(92,*) 'NREC_VECTOR_BIS'
WRITE(92,*)  NREC_VECTOR_BIS
WRITE(92,*) 'RECINI_VECT'
WRITE(92,*)  RECINI_VECT
WRITE(92,*) 'SIMMETRIA'
WRITE(92,*)  SIMMETRIA
WRITE(92,*) 'NUM_DIAG'
WRITE(92,*)  NUM_DIAG
KLPD = 0 
DO IA = 1 , NSYMA_PROP
WRITE(92,*) 'IA' , IA
WRITE(92,*) 'NLPD(IA)'
WRITE(92,*) NLPD(IA)
DO I = 1 , NLPD(IA)
KLPD = KLPD + 1
WRITE(92,*) 'KLPD' , KLPD
WRITE(92,*) 'LPSYM(KLPD) , DIPOLE_BLOCK(KLPD) , NUM(KLPD)'
WRITE(92,*) LPSYM(KLPD) ,  DIPOLE_BLOCK(KLPD) , NUM(KLPD)
WRITE(92,*) 'MU1D(KLPD) , MU2D(KLPD) , KAPPAD(KLPD)'
WRITE(92,*) MU1D(KLPD)  , MU2D(KLPD) , KAPPAD(KLPD)
END DO
END DO
CLOSE(92)



OPEN( 91, STATUS='UNKNOWN', ACCESS='SEQUENTIAL', FORM='UNFORMATTED')
WRITE(91) NSYMA_PROP , NBLOCKS_DIP , HAM_PIECES , KLPDTOT
WRITE(91) DIM_PROP
WRITE(91) NRECTOT_VECT
WRITE(91) NREC_VECTOR_BIS
WRITE(91) RECINI_VECT
WRITE(91) SIMMETRIA
WRITE(91) NUM_DIAG
KLPD = 0 
DO IA = 1 , NSYMA_PROP
WRITE(91) NLPD(IA)
DO I = 1 , NLPD(IA)
KLPD = KLPD + 1
WRITE(91) LPSYM(KLPD) , DIPOLE_BLOCK(KLPD) , NUM(KLPD)
WRITE(91) MU1D(KLPD)  , MU2D(KLPD) , KAPPAD(KLPD)
END DO
END DO
CLOSE(91)


end subroutine SYM_TOPROP_VECT_DIP
!!! TUTTE LE QUANTITA TENGONO CONTO DEL GROUND STATE IN SYMMETRY 1 !!! 
!!! TUTTE LE QUANTITA TENGONO CONTO DEL GROUND STATE IN SYMMETRY 1 !!! 
!!! TUTTE LE QUANTITA TENGONO CONTO DEL GROUND STATE IN SYMMETRY 1 !!! 
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!
!!! CALCOLA I VETTORI DA SCRIVERE NEL FILE 91, E CE LI SCRIVE ANCHE !!!


  end  module propagate_prepare 

