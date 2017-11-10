! **********************************************************************
! ******** Bulirsch-Stoer libraries taken and adapted from the *********
! ********                    Quantics code                    *********
! **********************************************************************

! **********************************************************************
! *                                                                    *
! *                        BULIRSCH-STOER (bslib.f)                    *
! *                                                                    *
! * Library module containing a Bulirsch-Stoer integrator.             *
! *                                                                    *
! * Contains:                                                          *
! *   BSStep:       The BS integration routine.                        *
! *   ModMidPoint:  Integrates a system of ODEs with modified midpoint *
! *                 method.                                            *
! *   AbsBSError:   Computes the absolute error of the current BS      *
! *                 integration step.                                  *
! *   RelBSError:   Computes the relative error of the current BS      *
! *                 integration step.                                  *
! *   PolyExtrapol: Extrapolates the estimated values for the solution *
! *                 to stepsize zero using a polynomial extrapolation. *
! *   WriteBSStep:  In the current form WriteBSStep is just a dummy    *
! *                 routine doing absolutely nothing; it is included   *
! *                 here for formal reasons, but can (if desired) be   *
! *                 extended easily such that it writes the size and   *
! *                 error of the current integration step to a file.   *
! *   BSErrorMsg:   Returns for a given error number a corresponding   *
! *                 error message.                                     *
! *                                                                    *
! **********************************************************************
module bslib

  implicit none
      
  private
  public :: bsstep,bserrormsg,polyextrapol,absbserror

  save

  integer, parameter :: dop=selected_real_kind(8)
      
contains

! **********************************************************************
! *                                                                    *
! *                            SUBROUTINE BSSTEP                       *
! *                                                                    *
! * Integrates a system of complex first order differential equations  *
! * employing the Bulirsch-Stoer extrapolation method. The routine     *
! * runs with both variable step size and order. BSStep makes only one *
! * single integration step, so it has to be embedded into a loop that *
! * calls BSStep until the desired time interval is integrated. The    *
! * ODE is of the form dPsi/dt = Func(AbsTime,Psi) =: DtPsi. All       *
! * computations are performed with double precision.                  *
! *                                                                    *
! * Input parameters:                                                  *
! *   Psi:             The (complex) initial-value vector.             *
! *   DtPsi:           Time derivative of the initial-value vector.    *
! *   PsiDim           Length of Psi and DtPsi vectors.                *
! *   NOffD:           No. non-zero off-diagonal Hamiltonian matrix    *
! *                    elements. This shouldn't really be here, but    *
! *                    this number has to be passed to the             *
! *                    matrix-vector multiplication routine.           *
! *   IntPeriod:       Lenght of time interval to be integrated.       *
! *   AbsTime:         Absolute time, i. e. Psi_initial=Psi(AbsTime).  *
! *   IntOrder:        Maximum integration order.                      *
! *   LargeStep:       Suggestion for the size of the large step (can  *
! *                    be equal to IntPeriod if no better value is     *
! *                    known).                                         *
! *   TolError:        Maximum error that is tolerated.                *
! *   Relaxation:      Logical flag. True if a relaxation calculation  *
! *                    is being performed, in which case, the result   *
! *                    of Func is multiplied by -i to convert from the *
! *                    real time derivative to the negative imaginary  *
! *                    time derivative.                                *
! *                                                                    *
! * Output parameters:                                                 *
! *   Psi:             Solution of the ODE at time                     *
! *                    AbsTime+ActualLargeStep.                        *
! *   DtPsi:           DtPsi is undetermined.                          *
! *   IntOrder         Optimal column of the extrapolation tableau.    *
! *   ActualLargeStep: Time interval that actually has been integrated *
! *                    (can be lower than LargeStep).                  *
! *   NextLargeStep:   Suggestion for the next large integration step  *
! *                    (can be used as value for LargeStep when BSStep *
! *                    integrates the next large step).                *
! *   SmallSteps:      Counter for the number of small integration     *
! *                    steps. (Note that SmallSteps isn't set to zero  *
! *                    at the beginning.)                              *
! *   ErrorCode:       Error code having the following meaning:        *
! *                    0: everything was o. k.,                        *
! *                    1: illegal integration order,                   *
! *                    2: stepsize underflow.                          *
! *                                                                    *
! * Other parameters:                                                  *
! *   AuxPsi:          Auxiliary array of (minimum) size               *
! *                    PsiDim*(IntOrder+2).                            *
! *                                                                    *
! * External routines:                                                 *
! *   Func:            Computes the time derivative DtPsi of Psi at    *
! *                    time AbsTime. Called as                         *
! *                    Func(AbsTime,Psi,DtPsi).                        *
! *   CalcError:       Determines the error of the current integration *
! *                    step. Called as                                 *
! *                    CalcError(Psi,Delta_Psi,PsiDim,Error)           *
! *   Extrapol:        Extrapolates the estimated values for Psi to    *
! *                    stepsize zero. Called as                        *
! *                    Extrapol(Column,Psi,Delta_Psi,AuxPsi,Tableau,   *
! *                             PsiDim,StepSize**2,IntOrder).          *
! *   WriteStep:       Writes the stepsize and error to a file.        *
! *                    Called as:                                      *
! *                    WriteStep(Iteration,Steps,Stepsize,Error).      *
! *                                                                    *
! * V6.0 MB                                                            *
! *                                                                    *
! **********************************************************************

  Subroutine BSStep (Psi,DtPsi,PsiDim,NOffD,IntPeriod,AbsTime,&
                     IntOrder,LargeStep,TolError,ActualLargeStep,&
                     NextLargeStep,SmallSteps,ErrorCode,AuxPsi,&
                     Func,CalcError,Extrapol,Relaxation)

!    use wrintegrat, only: writestep

    Implicit None

    Real(dop), Parameter :: RelativeMinStep = 1.0e-10_dop,MaxScale = 0.1_dop,&
                            ErrorSafety = 0.25_dop,ReductionSafety = 0.7_dop,&
                            MinReduction = 0.5_dop,MaxReduction = 1.0e-5_dop,&
                            Tiny = 1.0e-25_dop

    Integer, Parameter   :: MaxOrder = 16

    Integer, intent(out)                                    :: SmallSteps, &
                                                               ErrorCode
    Integer, intent(in)                                     :: PsiDim
    Integer*8, intent(in)                                   :: NOffd
    Integer, intent(inout)                                  :: IntOrder
    Real(dop), intent(in)                                   :: IntPeriod,&
                                                               AbsTime, &
                                                               TolError,&
                                                               LargeStep
    Real(dop), intent(out)                                  :: ActualLargeStep,&
                                                               NextLargeStep
    Complex(dop), dimension(PsiDim), intent(inout)          :: Psi,DtPsi
    Complex(dop), dimension(PsiDim,IntOrder+2), intent(out) :: AuxPsi
    Logical, intent(in)                                     :: Relaxation
    External                                                :: Func,CalcError,Extrapol

    Logical(kind=4)                               :: StepIsReduced
    Integer                                       :: D,P,P1,P2
    Integer, save                                 :: OptimalColumn
    Integer, dimension(MaxOrder+1), save          :: NumberOfSteps
    Real(dop)                                     :: MaxError, &
                                                     MinLargeStep, &
                                                     SquaredStepSize, &
                                                     Reduction,&
                                                     Work,WorkFactor, &
                                                     MinWork, &
                                                     ScaleFactor
    Real(dop), save                               :: IntAccuracy, &
                                                     NewTime

    Real(dop), dimension(MaxOrder-1)              :: Error
    Real(dop), dimension(MaxOrder+1), save        :: Expense
    Real(dop), dimension(MaxOrder,MaxOrder), save :: Correction
    

    real(dop), save                               :: OldAccuracy = -1.0_dop
    integer, save                                 :: OldIntOrder = -1
    logical(kind=4), save                         :: CallIsFirstCall = .True.

! --- CHECK INTEGRATION ORDER ---
    If ((IntOrder .LT. 2) .Or. (IntOrder .GT. MaxOrder)) Then
       ErrorCode = 1
       Return
    EndIf

! --- INITIALISE VARIABLES ---
    ErrorCode = 0
    MinLargeStep = RelativeMinStep*IntPeriod
    NextLargeStep = -Tiny
    NewTime = -Tiny


! --- REINITIALIZE IF TOLERANCE OR ORDER IS NEW ---
    If ((TolError .NE. OldAccuracy) .Or.&
         (IntOrder .NE. OldIntOrder)) Then

!    --- INITIALISE VARIABLES ---
       OldAccuracy = TolError
       IntAccuracy = ErrorSafety*TolError
       OldIntOrder = IntOrder
       Do P = 1,IntOrder+1
          NumberOfSteps(P) = 2*P
       EndDo

!    --- COMPUTE EXPENSE COEFFICIENTS ---
       Expense(1) = NumberOfSteps(1)+1
       Do P = 2,IntOrder+1
          Expense(P) = Expense(P-1)+NumberOfSteps(P)
       EndDo

!    --- COMPUTE CORRECTION FACTOR ---
       Do P = 2,IntOrder
          Do P1 = 1,P-1
             Correction(P1,P) = IntAccuracy**((Expense(P1+1)&
                  -Expense(P+1))/((2*P1+1)*(Expense(P+1)&
                  -Expense(1)+1.0_dop)))
          EndDo
       EndDo

!    --- DETERMINE OPTIMAL ROW NUMBER FOR CONVERGENCE ---
       Do OptimalColumn = 2,IntOrder-1
          If (Expense(OptimalColumn+1) .GT. Expense(OptimalColumn)&
               *Correction(OptimalColumn-1,OptimalColumn)) Then
             Goto 10
          EndIf
       EndDo
10     IntOrder = OptimalColumn
    EndIf

! --- INITIALIZE VARIABLES ---
    ActualLargeStep = LargeStep
    StepIsReduced = .False.
    MinWork = 1.0e35_dop
    Do D = 1,PsiDim
       AuxPsi(D,1) = Psi(D)
    EndDo

! --- RE-ESTABLISH THE ORDER WINDOW ---
    If ((ActualLargeStep .NE. NextLargeStep) .Or.&
         (AbsTime .NE. NewTime)) Then
       CallIsFirstCall = .True.
       OptimalColumn = IntOrder
    EndIf

! --- MAKE ONE LARGE INTEGRATION STEP ---
300 Continue


! --- LOOP OVER EACH COLUMN IN THE EXTRAPOLATION TABLEAU ---
    reduction = 0
    Do P = 1,IntOrder

!    --- INITIALIZE VARIABLES ---
       NewTime = AbsTime+ActualLargeStep
       SquaredStepSize = (ActualLargeStep/NumberOfSteps(P))**2
       SmallSteps = SmallSteps+NumberOfSteps(P)
         
!    --- CHECK IF STEPSIZE IS TOO SMALL ---
       If ((ActualLargeStep .LT. MinLargeStep) .And.&
            (NewTime .LT. IntPeriod)) Then
          ErrorCode = 2
          Return
       EndIf

!    --- INTEGRATE WITH MODIFIED MIDPOINT METHOD ---

       Call ModMidpoint(AuxPsi(1,1),Psi,DtPsi,AuxPsi(1,3),AuxPsi(1,2),&
                        PsiDim,NOffD,ActualLargeStep,AbsTime,&
                        NumberOfSteps(P),Func,Relaxation)

!    --- EXTRAPOLATE TO ZERO STEP SIZE ---
       
       Call Extrapol(P,Psi,AuxPsi(1,2),AuxPsi(1,3),AuxPsi(1,4),PsiDim,&
                     SquaredStepSize,IntOrder)

!    --- COMPUTE ERROR ---

       MaxError = 0.0_dop
       p1=0
       If (P .GE. 2) Then
          Call CalcError(Psi,AuxPsi(1,2),PsiDim,MaxError)
          P1 = P-1
          Error(P1) = (MaxError/IntAccuracy)**(1.0_dop/(2*P1+1.0_dop))
       EndIf

!    --- WRITE STEPSIZE AND ERROR IF DESIRED ---

!       Call WriteStep(P,NumberOfSteps(P),ActualLargeStep,&
!                      MaxError,AbsTime)

!    --- CHECK IF COLUMN LIES IN ORDER WINDOW ---

       If ((P .GE. 2) .And. ((P .GE. OptimalColumn-1) .Or.&
            CallIsFirstCall)) Then
 
!       --- EXIT LOOP WHEN CONVERGED ---

          If (MaxError .LT. TolError) Then
             Goto 500
          EndIf

!       --- CHECK FOR POSSIBLE STEPSIZE REDUCTION ---

          If ((P .Eq. IntOrder) .Or. (P .Eq. OptimalColumn+1)) Then
             Reduction = ReductionSafety/Error(P1)
             Goto 400
          ElseIf (P .Eq. OptimalColumn) Then
             If (Correction(OptimalColumn-1,OptimalColumn)&
                  .LT. Error(P1)) Then
                Reduction = 1.0_dop/Error(P1)
                Goto 400
             EndIf
          ElseIf (OptimalColumn .Eq. IntOrder) Then
             If (Correction(P1,IntOrder-1) .LT. Error(P1)) Then
                Reduction = ReductionSafety*Correction(P1,IntOrder-1)&
                     /Error(P1)
                Goto 400
             EndIf
          ElseIf (Correction(P1,OptimalColumn) .LT. Error(P1)) Then
             Reduction = Correction(P1,OptimalColumn-1)/Error(P1)
             Goto 400
          EndIf
       EndIf
    EndDo

! --- REDUCE STEPSIZE IF STEP WASN'T SUCCESFUL ---

400 Continue
    Reduction = Min(Reduction,MinReduction)
    Reduction = Max(Reduction,MaxReduction)
    ActualLargeStep = ActualLargeStep*Reduction
    StepIsReduced = .True.

! --- REPEAT INTEGRATION STEP ---

    Goto 300

! --- FINISH INTEGRATION STEP ---

500 Continue

! --- UPDATE VARIABLES ---

    CallIsFirstCall = .False.

! --- COMPUTE OPTIMAL ROW AND STEPSIZE FOR CONVERGENCE ---

    scalefactor = 0
    Do P2 = 1,P1
       WorkFactor = Max(Error(P2),MaxScale)
       Work = WorkFactor*Expense(P2+1)
       If (Work .LT. MinWork) Then
          ScaleFactor = WorkFactor
          MinWork = Work
          OptimalColumn = P2+1
       EndIf
    EndDo
    NextLargeStep = ActualLargeStep/ScaleFactor

! --- CHECK FOR POSSIBLE ORDER INCREASE ---

    If ((OptimalColumn .GE. P) .And. (OptimalColumn .LT. IntOrder)&
         .And. (.Not. StepIsReduced)) Then
       WorkFactor = Max(ScaleFactor/Correction(OptimalColumn-1,&
            OptimalColumn),MaxScale)
       If (Expense(OptimalColumn+1)*WorkFactor .LE. MinWork) Then
          NextLargeStep = ActualLargeStep/WorkFactor
          OptimalColumn = OptimalColumn+1
       EndIf
    EndIf
    
    Return
  End Subroutine BSStep

! **********************************************************************
! *                                                                    *
! *                         SUBROUTINE MODMIDPOINT                     *
! *                                                                    *
! * Integrates a system of ordinary differential equations employing a *
! * modified midpoint method.                                          *
! *                                                                    *
! * Input parameters:                                                  *
! *   SavedPsi:      Initial value vector.                             *
! *   DtPsi:         Time derivative of SavedPsi.                      *
! *   PsiDim:        Length of the psi vectors.                        *
! *   LargeStepSize: Time interval to be integrated.                   *
! *   InitialTime:   Absolute time, i. e. SavedPsi=Psi(InitialTime).   *
! *   NumberOfSteps: Number of integration steps to be made.           *
! *   Relaxation:    Logical flag. True if a relaxation calculation    *
! *                  is being performed, in which case, the result     *
! *                  of Func is multiplied by -i to convert from the   *
! *                  real time derivative to the negative imaginary    *
! *                  time derivative.                                  *
! *                                                                    *
! * Output parameters:                                                 *
! *   EstimatedPsi:  Solution of the ODE at time                       *
! *                  InitialTime+LargeStepSize.                        *
! *                                                                    *
! * Other parameters:                                                  *
! *   Psi1,Psi2:     Auxiliary arrays each of size PsiDim.             *
! *                                                                    *
! * External routines:                                                 *
! *   Func:          Computes the time derivative DtPsi of Psi and is  *
! *                  called as                                         *
! *                  Func(Time,Psi,DtPsi,Cdata,Rdata,Idata,Ldata).     *
! * V6.0 MB                                                            *
! * V7.0 GW addition of Cdata,Rdata,Idata,Ldata arrays                 *
! *                                                                    *
! **********************************************************************

  Subroutine ModMidpoint (SavedPsi,EstimatedPsi,DtPsi,Psi1,Psi2,&
                          PsiDim,NOffD,LargeStepSize,InitialTime,&
                          NumberOfSteps,Func,Relaxation)

    Implicit None

    Integer, intent(in)                          :: PsiDim,NumberOfSteps
    Integer*8, intent(in)                        :: NOffD
    Real(dop), intent(in)                        :: LargeStepSize,InitialTime
    Complex(dop), dimension(PsiDim), intent(in)  :: SavedPsi,DtPsi
    Complex(dop), dimension(PsiDim), intent(out) :: EstimatedPsi,Psi2,Psi1
    Logical, intent(in)                          :: relaxation

    External                :: Func
    Integer                 :: D,P
    Real(dop)               :: StepSize,DoubleStepSize,CurrentTime
    Complex(dop)            :: Swap
    Complex(dop), parameter :: ci=(0._dop,1._dop)
    
! --- INITIALIZE VARIABLES ---

    StepSize = LargeStepSize/NumberOfSteps
    DoubleStepSize = 2.0_dop*StepSize
    CurrentTime = InitialTime+StepSize

! --- STORE PSI AND FIRST ESTIMATE OF NEXT PSI ---

    Do D = 1,PsiDim
       Psi1(D) = SavedPsi(D)
       Psi2(D) = SavedPsi(D)+StepSize*DtPsi(D)
    EndDo
         
! --- EVALUATE FUNCTION WITH ESTIMATED PSI ---

    Call Func(CurrentTime,PsiDim,NOffD,Psi2,EstimatedPsi)
    if (Relaxation) EstimatedPsi=-ci*EstimatedPsi
    
! --- LOOP OVER NUMBER OF INTEGRATION STEPS ---

    Do P = 2,NumberOfSteps

!    --- PERFORM NEXT INTEGRATION STEP ---

       CurrentTime = CurrentTime+StepSize
       Do D = 1,PsiDim
          Swap = Psi1(D)+DoubleStepSize*EstimatedPsi(D)
          Psi1(D) = Psi2(D)
          Psi2(D) = Swap
       EndDo

!    --- EVALUATE FUNCTION WITH ESTIMATED PSI ---

       Call Func(CurrentTime,PsiDim,NOffD,Psi2,EstimatedPsi)
       if (Relaxation) EstimatedPsi=-ci*EstimatedPsi
    EndDo

! --- TAKE MEAN VALUE FROM LAST AND NEXT ESTIMATE ---

    Do D = 1,PsiDim
       EstimatedPsi(D) = 0.5_dop*(Psi1(D)+Psi2(D)&
            +StepSize*EstimatedPsi(D))
    EndDo
      
    Return
  End Subroutine ModMidpoint

! **********************************************************************
! *                                                                    *
! *                         SUBROUTINE ABSBSERROR                      *
! *                                                                    *
! * Estimates the absolute error of the current Bulirsch-Stoer         *
! * iteration.                                                         *
! *                                                                    *
! * Input parameters:                                                  *
! *   Psi:      Vector containing the current solution. (Here Psi is a *
! *             dummy parameter but may be used e. g. if a relative    *
! *             error is desired.)                                     *
! *   PsiError: Vector containing error estimates for each component   *
! *             of the system of ODEs.                                 *
! *   PsiDim:   Length of PsiError.                                    *
! *                                                                    *
! * Output parameters:                                                 *
! *   Error:    Estimated error.                                       *
! *                                                                    *
! * V6.0 MB                                                            *
! * V7.0 GW addition of Data arrays                                    *
! *                                                                    *
! **********************************************************************

  Subroutine AbsBSError (Psi,PsiError,PsiDim,Error)

    Implicit None

    Integer, intent(in)                         :: PsiDim
    Integer                                     :: D
    Real(dop), intent(out)                      :: Error
    Complex(dop), dimension(PsiDim), intent(in) :: Psi,PsiError
    complex(dop)                                :: cdum

! dummy code to prevent unused variable warning on compilation
    cdum = psi(1)

! --- ESTIMATE THE ERROR ---

    Error = 0.0_dop
    Do D = 1,PsiDim
       Error = Max(Error,Dble(Abs(PsiError(D))))
    EndDo
         
    Return
  End Subroutine AbsBSError

! **********************************************************************
! *                                                                    *
! *                         SUBROUTINE POLYEXTRAPOL                    *
! *                                                                    *
! * Extrapolates the estimated values for psi to stepsize zero using a *
! * polynomial extrapolation.                                          *
! *                                                                    *
! * Input parameters:                                                  *
! *   Column:          Column in the extrapolation tableau to be       *
! *                    filled.                                         *
! *   Psi:             Vector to be filled in extrapolation tableau.   *
! *   PsiDim:          Length of Psi.                                  *
! *   Tableau:         Extrapolation tableau of size                   *
! *                    PsiDim*(IntOrder-1).                            *
! *   SquaredStepSize: Squared large step size.                        *
! *   IntOrder:        Maximum integration order.                      *
! *                                                                    *
! * Output parameters:                                                 *
! *   Psi:             Extrapolated vector.                            *
! *   DeltaPsi:        Error estimate of extrapolated vector.          *
! *   Tableau:         Updated extrapolation tableau.                  *
! *                                                                    *
! * Other parameters:                                                  *
! *   Difference:      Auxiliary vector of length PsiDim.              *
! *                                                                    *
! * V6.0 MB                                                            *
! *                                                                    *
! **********************************************************************

  Subroutine PolyExtrapol (Column,Psi,DeltaPsi,Difference,Tableau,&
                           PsiDim,SquaredStepSize,IntOrder)

    Implicit None

    Integer, Parameter                                      :: MaxOrder = 16
    Integer, intent(in)                                     :: Column,PsiDim,IntOrder
    Real(dop), intent(in)                                   :: SquaredStepSize
    Complex(dop), dimension(PsiDim), intent(out)            :: Difference,DeltaPsi
    Complex(dop), dimension(PsiDim), intent(inout)          :: Psi
    Complex(dop), dimension(PsiDim,IntOrder-1), intent(out) :: Tableau

    Integer                              :: D,P
    Real(dop)                            :: Factor,OldFactor,Help
    Real(dop), dimension(MaxOrder), save :: OldSquaredStepSize
    Complex(dop)                         :: Delta,InterimTableau
    

! --- INITIALIZE VARIABLES ---

    OldSquaredStepSize(Column) = SquaredStepSize
    Do D = 1,PsiDim
       DeltaPsi(D) = Psi(D)
    EndDo
      
! --- FILL IN TABLEAU ---

    If (Column .Eq. 1) Then

!    --- STORE FIRST ESTIMATE IN FIRST COLUMN ---

       Do D = 1,PsiDim
          Tableau(D,1) = Psi(D)
       EndDo
    Else
         
!    --- INITIALIZE COLUMN DIFFERENCE ---

       Do D = 1,PsiDim
          Difference(D) = Psi(D)
       EndDo

!    --- LOOP OVER EACH PREVIOUS COLUMN ---

       Do P = 1,Column-1
            
!       --- INITIALIZE VARIABLES ---

          Help = 1.0_dop/(OldSquaredStepSize(Column-P)-SquaredStepSize)
          Factor = SquaredStepSize*Help
          OldFactor = OldSquaredStepSize(Column-P)*Help

!       --- PROPAGATE TABLEAU ONE DIAGONAL FURTHER ---

          Do D = 1,PsiDim
             InterimTableau = Tableau(D,P)
             Tableau(D,P) = DeltaPsi(D)
             Delta = Difference(D)-InterimTableau
             DeltaPsi(D) = Factor*Delta
             Difference(D) = OldFactor*Delta
             Psi(D) = Psi(D)+DeltaPsi(D)
          EndDo
       EndDo

!    --- FILL FINAL COLUMN OF THE TABLEAU ---

       If (Column .LT. IntOrder) Then
          Do D = 1,PsiDim
             Tableau(D,Column) = DeltaPsi(D)
          EndDo
       EndIf
    EndIf

    Return
  End Subroutine PolyExtrapol

! **********************************************************************
! *                                                                    *
! *                        SUBROUTINE BSERRORMSG                       *
! *                                                                    *
! * Generates for a given error number returned by "BSStep" a          *
! * corresponding error message.                                       *
! *                                                                    *
! * Input parameters:                                                  *
! *   Error: Error code returned by BSStep.                            *
! *                                                                    *
! * Output parameters:                                                 *
! *   Msg:   Error message.                                            *
! *                                                                    *
! * V7.0 MB                                                            *
! *                                                                    *
! **********************************************************************

  Subroutine BSErrorMsg (Error,Msg)
        
    Implicit None

    Integer, intent(in)           :: Error
    Character(len=*), intent(out) :: Msg

! --- GENERATE ERROR MESSAGE ---

    If (Error .Eq. 1) Then
       Msg = 'Illegal integration order'
    ElseIf (Error .Eq. 2) Then
       Msg = 'Stepsize underflow'
    Else
       Msg = 'Unknown error occurred'
    EndIf

    Return
  End Subroutine BSErrorMsg
  
end module bslib
