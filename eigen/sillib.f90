! **********************************************************************
! ******** SIL libraries taken and adapted from the MCTDH code *********
! **********************************************************************

! **********************************************************************
! *                                                                    *
! *                 SHORT ITERATIVE LANCZOS (sillib.f90)               *
! *                                                                    *
! * Library module containing a short iterative Lanczos integrator.    *
! *                                                                    *
! * Contains:                                                          *
! *   SILStep:      The (real) Lanczos integration routine.            *
! *   SILErrorMsg:  Returns for a given error number a corresponding   *
! *                 error message.                                     *
! **********************************************************************
      module sillib

      implicit none

      private
      public :: silstep,silerrormsg

      save

      integer, parameter :: dop=selected_real_kind(8)

    contains


! **********************************************************************
! *                                                                    *
! *                 SUBROUTINE SILSTEP                                 *
! *                                                                    *
! * Integrates a system of complex linear first order differential     *
! * equations with constant and Hermitian Hamiltonian employing the    *
! * short iterative Lanczos method. The routine runs with both         *
! * variable step size and order. (First it increases the order to     *
! * achieve the desired accuracy. If this failes within the maximum    *
! * order, the stepsize is reduced.) SILStep makes only one single     *
! * integration step, so it has to be imbedded into a loop that calls  *
! * SILStep until the desired time interval is integrated. The ODE is  *
! * of the form i dPsi/dt = H|Psi> = Func(Time,Psi,DtPsi) =: DtPsi or  *
! * dPsi/dt = -i H|Psi> = Func(Time,Psi,DtPsi) =: DtPsi, depending on  *
! * the flag "StdForm". All computations are performed with double     *
! * precision.                                                         *
! *                                                                    *
! * Input parameters:                                                  *
! *   Psi:       The (complex) initial-value vector.                   *
! *   DtPsi:     Action of Hamiltonian on the initial-value vector,    *
! *              i.e. H|Psi> or -i H|Psi>, depending on "StdForm".     *
! *   PsiDim     Length of Psi and DtPsi vectors.                      *
! *   IntPeriod: Lenght of time interval to be integrated.             *
! *   IntOrder:  Maximum integration order.                            *
! *   TolError:  Maximum error that is tolerated.                      *
! *   Relax:     Flag for relaxation calculation. If true, Psi is      *
! *              relaxated, else propagated.                           *
! *   Restart:   Flag for restarting the integrator. If true, the      *
! *              Krylov space is built up before propagation, else the *
! *              old Krylov vectors are used.                          *
! *   StdForm:   If true, Func(Time,Psi,DtPsi) = -i H|Psi>, else       *
! *              Func(Time,Psi,DtPsi) = H|Psi>.                        *
! *   Steps:     Number of steps made so far (Steps is passed to       *
! *              "WriteStep").                                         *
! *   Krylov:    Matrix with minimum size PsiDim*(IntOrder-1) with     *
! *              columns containing the Krylov vectors H^2|psi> to     *
! *              H^TrueOrder|psi>. The Krylov vectors are needed on    *
! *              entry only when Restart is false.                     *
! *                                                                    *
! * Output parameters:                                                 *
! *   Psi:       Propagated Psi.                                       *
! *   DtPsi:     Normalised first Krylov vector, i. e. (H-<H>)|Psi>.   *
! *   Krylov:    Matrix with minimum size PsiDim*(IntOrder-1) with     *
! *              columns containing the Krylov vectors H^2|psi> to     *
! *              H^TrueOrder|psi>.                                     *
! *   Stepsize:  Time interval that actually has been integrated (can  *
! *              be lower than IntPeriod).                             *
! *   TrueOrder: The order that has actually been used (may be less    *
! *              than IntOrder).                                       *
! *   Steps:     Same as on entry.                                     *
! *   ErrorCode: Error code having the following meaning:              *
! *              0: everything was o. k.,                              *
! *              1: illegal integration order,                         *
! *              2: stepsize underflow,                                *
! *              3: diagonalization failed,                            *
! *              4: stepsize too large for restart.                    *
! *                                                                    *
! * External routines:                                                 *
! *   Func:      Computes the action of the Hamiltonian on Psi.        *
! *              Called as                                             *
! *              Func(AbsTime,Psi,DtPsiR)                              *
! *   WriteStep: Writes the stepsize and error to a file.              *
! *              Called as WriteStep(Steps,Order,Stepsize,Error).      *
! *                                                                    *
! * Some new Variables:                                                *
! *  rlaxnum : This is the argument of the keyword relaxation.         *
! *            If relaxation has no argument, rlaxnum=-1 .             *
! *  Relax   : logical set true for a relaxation run.                  *
! *                                                                    *
! * V6.0 MB                                                            *
! * V8.2 05/01 HDM  Addition of 'relaxation' to exited states.         *
! *                 Eigenvector-space now dynamically allocated.       *
! **********************************************************************

      Subroutine SILStep(Psi,DtPsi,PsiDim,IntPeriod,IntOrder,TolError,&
                         Relax,Restart,StdForm,Steps,Krylov,Stepsize,&
                         TrueOrder,ErrorCode,Time,Func,&
                         EigenVector,EigenVal,Diagonal,OffDg2,OffDiag)

!use wrintegrat, only: writestep

      Implicit None

      Real(dop), Parameter  :: RelativeMinStep = 1.0e-10_dop,Tiny = 1.0e-18_dop

      Logical(kind=4), intent(in)                               :: Relax,Restart,StdForm
      Integer, intent(inout)                                    :: Steps
      Integer, intent(out)                                      :: TrueOrder,ErrorCode
      Integer, intent(in)                                       :: PsiDim,IntOrder
      Real(dop), intent(out)                                    :: Stepsize
      Real(dop), intent(in)                                     :: IntPeriod,TolError,Time
      Complex(dop), dimension(PsiDim), intent(inout)            :: Psi,DtPsi
      Complex(dop), dimension(PsiDim,2:IntOrder), intent(inout) :: Krylov
      External                                                  :: Func

      Integer                                                      :: D,P,Q
      Integer, save                                                :: OldOrder
      Real(dop)                                                    :: Beta,Error,MinStepsize,beta1
      Real(dop), save                                              :: OldError,OldStepsize
      Real(dop), dimension(0:IntOrder,0:IntOrder+2), intent(inout) :: EigenVector
      Real(dop), dimension(0:IntOrder), intent(out)                :: Diagonal,EigenVal
      Real(dop), dimension(IntOrder), intent(out)                  :: OffDiag
      Real(dop), dimension(IntOrder+1), intent(out)                :: OffDg2
      Complex(dop)                                                 :: Alpha,CBeta,PreFactor,CInverse,sum


! --- CHECK INTEGRATION ORDER ---

      If ( IntOrder .LT. 2 ) Then
         ErrorCode = 1
         write(6,'(a,i4,a)') 'IntOrder =',IntOrder,' is too small!'
         Return
      EndIf

! --- INITIALIZE VARIABLES ---

      Stepsize = IntPeriod
      MinStepsize = Abs(RelativeMinStep*IntPeriod)
      ErrorCode = 0
      If (.not.StdForm) Then
         PreFactor = (1.0_dop,0.0_dop)
      ElseIf (Relax) Then
         PreFactor = (-1.0_dop,0.0_dop)
      Else
         PreFactor = (0.0_dop,-1.0_dop)
      EndIf


! --- RETURN RESULT IF DIFFERENTIAL EQUATION IS NOT A SYSTEM ---

      If (PsiDim .Eq. 1) Then
         TrueOrder = 1
         Error = 0.0_dop
!         Call WriteStep(Steps,TrueOrder,Stepsize,Error,Time)
         Psi(1) = Exp(StepSize*DtPsi(1)/Psi(1))*Psi(1)
         Return
      EndIf

! --- SKIP CALCULATION OF KRYLOV VECTORS IF DESIRED ---

      If (Restart) Then
         OldError=0.0_dop
         OldOrder=0
         OldStepsize=0.0_dop
      Else
         If (Abs(Stepsize) .GT. Abs(OldStepsize)+MinStepsize) Then
            ErrorCode = 4
            Return
         EndIf
         TrueOrder = OldOrder
         Error = Abs(OldError*(Stepsize/OldStepsize)**TrueOrder)
!         Call WriteStep(Steps,TrueOrder,Stepsize,Error,Time)
         Goto 200
      EndIf

! --- FURTHER INITIALIZE VARIABLES ---

      Alpha = (0.0_dop,0.0_dop)
      Beta = 0.0_dop
      TrueOrder = 1

! --- COMPUTE FIRST DIAGONAL ELEMENT ---

      Do D = 1,PsiDim
         Alpha = Alpha+DConjg(Psi(D))*DtPsi(D)
      EndDo
      Diagonal(0) = Alpha/PreFactor

! --- DETERMINE CORRESPONDING BASIS VECTOR AND OFF-DIAGONAL ELEMENT ---

      Do D = 1,PsiDim
         DtPsi(D) = DtPsi(D)-Alpha*Psi(D)
         Beta = Beta+Dble(DConjg(DtPsi(D))*DtPsi(D))
      EndDo

      Beta = Sqrt(Beta)
      OffDiag(1) = Beta
      beta1 = Beta
      CBeta = PreFactor*Beta

! --- NORMALIZE BASIS VECTOR ---

      If (Beta .LT. 1.0e-20_dop) Then
         Beta = 0.0_dop
         CInverse = (0.0_dop,0.0_dop)
      Else
         CInverse = 1.0_dop/CBeta
      EndIf
      Do D = 1,PsiDim
         DtPsi(D) = CInverse*DtPsi(D)
      EndDo

! --- COMPUTE ERROR ESTIMATE ---

      Error = Abs(Beta*Stepsize)

! --- WRITE STEPSIZE AND ERROR IF DESIRED ---

!      Call WriteStep(Steps,TrueOrder,Stepsize,Error,Time)

! --- BUILD UP THE KRYLOV SPACE ---
 
 100  Continue
        
! --- CHECK IF STEPSIZE IS TOO SMALL ---

      If (Abs(Stepsize) .LT. Abs(MinStepSize)) Then
         ErrorCode = 2
         Return
      EndIf

! --- RE-INITIALIZE VARIABLES ---

      Alpha = (0.0_dop,0.0_dop)
      Beta = 0.0_dop
      TrueOrder = TrueOrder+1


! --- EVALUATE FUNCTION WITH LAST BASIS VECTOR ---

      If (TrueOrder .Eq. 2) Then
         Call Func(Time,DtPsi,Krylov(1,2))
      Else
         Call Func(Time,Krylov(1,TrueOrder-1),Krylov(1,TrueOrder))
      EndIf

! --- COMPUTE DIAGONAL ELEMENT ---

! Note that the Krylov vectors number zero and one aren't stored in
! "Krylov" but in "Psi" and "DtPsi", respectively.

      If (TrueOrder .Eq. 2) Then
         Do D = 1,PsiDim
            Alpha = Alpha+DConjg(DtPsi(D))*Krylov(D,2)
         EndDo
      Else
         Do D = 1,PsiDim
            Alpha = Alpha+DConjg(Krylov(D,TrueOrder-1))&
                    *Krylov(D,TrueOrder)
         EndDo
      EndIf
      Diagonal(TrueOrder-1) = Alpha/PreFactor

! --- COMPUTE OFF-DIAGONAL ELEMENT AND BASIS VECTOR ---

      If (TrueOrder .Eq. 2) Then
         Do D = 1,PsiDim
            Krylov(D,2) = Krylov(D,2)-Alpha*DtPsi(D)-CBeta*Psi(D)
            Beta = Beta+Dble(DConjg(Krylov(D,2))*Krylov(D,2))
         EndDo
      ElseIf (TrueOrder .Eq. 3) Then
         Do D = 1,PsiDim
            Krylov(D,3) = Krylov(D,3)-Alpha*Krylov(D,2)-CBeta*DtPsi(D)
            Beta = Beta+Dble(DConjg(Krylov(D,3))*Krylov(D,3))
         EndDo
      Else
         Do D = 1,PsiDim
            Krylov(D,TrueOrder) = Krylov(D,TrueOrder)&
               -Alpha*Krylov(D,TrueOrder-1)-CBeta*Krylov(D,TrueOrder-2)
            Beta = Beta+DConjg(Krylov(D,TrueOrder))*Krylov(D,TrueOrder)
         EndDo
      EndIf

      Beta = Sqrt(Beta)
      OffDiag(TrueOrder) = Beta
      CBeta = PreFactor*Beta

! --- NORMALIZE BASIS VECTOR ---

      If (Beta .LT. 1.0e-20_dop) Then
         OffDiag(TrueOrder) = 0.0_dop
         Beta = 0.0_dop
         CInverse = (0.0_dop,0.0_dop)
      Else
         CInverse = 1.0_dop/CBeta
      EndIf
      
      Do D = 1,PsiDim
         Krylov(D,TrueOrder) = CInverse*Krylov(D,TrueOrder)
      EndDo

! --- COMPUTE ERROR ESTIMATE ---
      
      Error = Abs(Error*Beta*Stepsize/Dble(TrueOrder))

! --- WRITE STEPSIZE AND ERROR IF DESIRED ---

!      Call WriteStep(Steps,TrueOrder,Stepsize,Error,Time)

! --- CONTINUE IF NOT CONVERGED ---
      
      If((TrueOrder.LT.IntOrder) .And. (Error.GT.TolError)&
         .And. (Beta.gt.1.0e-20_dop) )  Goto 100
            
! --- ESTIMATE NEXT MATRIX ELEMENT ---

      Diagonal(TrueOrder) = Diagonal(TrueOrder-1)
!      Diagonal(TrueOrder) = Diagonal(TrueOrder-1) + beta

! --- DECREASE STEPSIZE IF LANCZOS DIDN'T CONVERGE ---
      
      If (Error .GT. TolError) Then
         Stepsize = Stepsize*(TolError/Error)**(1.0_dop/Dble(TrueOrder))
!         Call WriteStep(Steps,TrueOrder,Stepsize,TolError,Time)
      EndIf


! --- SAVE SOME VALUES IN CASE THAT NEXT CALL IS NON-RESTART CALL ---

      OldOrder = TrueOrder
      OldStepsize = Stepsize
      OldError = Error
    
! --- Save Diagonal and off diagonal values since the diagonalisation ---
! --- routine will overwrite those. ---

      call  cpvxd(Diagonal,EigenVal,TrueOrder+1)
      call  cpvxd(OffDiag,OffDg2,TrueOrder)

! --- DIAGONALIZE LANCZOS MATRIX --- CALL LAPACK ROUTINE --- 
! --- EigenVector(*,IntOrder+1) serves as work array 
      call DSTEQR('I',TrueOrder+1,EigenVal,OffDg2,EigenVector,&
                  IntOrder+1,EigenVector(1,IntOrder+1),ErrorCode)

      If (ErrorCode .NE. 0) Then
         ErrorCode = 3
         Return
      EndIf


 200  Continue                  ! Jump to here, if not RESTART.

! --- PROPAGATE WAVEFUNCTION ---

! Note again that the Krylov vectors number zero and one aren't stored
! in "Krylov" but in "Psi" and "DtPsi", respectively.

      Do P = 0,TrueOrder
         Sum = (0.0_dop,0.0_dop)
         Do Q = 0,TrueOrder
            Sum = Sum+Exp(PreFactor*EigenVal(Q)*Stepsize)*&
                 EigenVector(0,Q)*EigenVector(P,Q)
         EndDo
         If (P .Eq. 0) Then
            Do D = 1,PsiDim
               Psi(D) = Sum*Psi(D)
            EndDo
         ElseIf (P .Eq. 1) Then
            Do D = 1,PsiDim
               Psi(D) = Psi(D)+Sum*DtPsi(D)
            EndDo
         Else
            Do D = 1,PsiDim
               Psi(D) = Psi(D)+Sum*Krylov(D,P)
            EndDo
         EndIf
      EndDo
         
      Return
      End Subroutine SILStep


! **********************************************************************
! *                                                                    *
! *                       SUBROUTINE SILERRORMSG                       *
! *                                                                    *
! * Generates for a given error number returned by "SILStep/SILStep2/  *
! * SILStep3" a corresponding error message.                           *
! *                                                                    *
! * Input parameters:                                                  *
! *   Error: Error code returned by SILStep/SILStep2.                  *
! *                                                                    *
! * Output parameters:                                                 *
! *   Msg:   Error message.                                            *
! *                                                                    *
! * V7.0 MB                                                            *
! *                                                                    *
! **********************************************************************

      Subroutine silerrormsg (Error,Msg)

      Implicit None

      Integer, intent(in)          :: Error
      Character(len=*), intent(out) :: Msg

! --- GENERATE ERROR MESSAGE ---

      If (Error .Eq. 1) Then
         Msg = 'Illegal integration order'
      ElseIf (Error .Eq. 2) Then
         Msg = 'Stepsize underflow'
      ElseIf (Error .Eq. 3) Then
         Msg = 'Diagonalisation failed'
      ElseIf (Error .Eq. 4) Then
         Msg = 'Illegal initial stepsize'
      Else
         Msg = 'Unknown error occurred'
      EndIf

      Return
      End Subroutine silerrormsg

! **********************************************************************
! Library subroutine cpvxd
!
! copies a real vector to a different real vector
!     w(i) = v(i)
! **********************************************************************

      subroutine cpvxd (v,w,dim)

      integer                                :: i
      integer, intent(in)                    :: dim
      real(dop), dimension(dim), intent(in)  :: v
      real(dop), dimension(dim), intent(out) :: w

      do i = 1,dim
         w(i) = v(i)
      enddo

      return
      end subroutine

! **********************************************************************

      end module sillib
