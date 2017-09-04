! **********************************************************************
! ****** CSIL libraries taken and adapted from the Quantics code *******
! **********************************************************************

! **********************************************************************
! *                                                                    *
! *               COMPLEX SHORT ITERATIVE LANCZOS (csillib.f)          *
! *                                                                    *
! * Library module storing a complex short iterative Lanczos-Arnoldi   *
! * integrator.                                                        *
! *                                                                    *
! * Contains:                                                          *
! *   CSILStep:      The Lanczos-Arnoldi integration routine.          *
! *   CSILErrorMsg:  Returns for a given error number a corresponding  *
! *                  error message.                                    *
! *                                                                    *
! **********************************************************************
      module csillib

        implicit none

        private
        public :: csilstep,csilerrormsg

        save
        
        integer, parameter :: dop=selected_real_kind(8)

      contains
        
! **********************************************************************
! *                                                                    *
! *                SUBROUTINE CSILSTEP                                 *
! *                                                                    *
! * Integrates a system of complex linear first order differential     *
! * equations with time-independent Hamiltonian (not necessarily       *
! * hermitian) employing the short iterative Lanczos-Arnoldi method.   *
! * The routine runs with both variable step size and order. (First    *
! * the order is increased to achieve the desired accuracy. If this    *
! * failes within the maximum order, the stepsize is reduced.)         *
! * CSILStep makes only one single integration step, so it has to be   *
! * imbedded into a loop that calls CSILStep until the desired time    *
! * interval is integrated. The ODE has the form i dPsi/dt = H|Psi> =  *
! * Func(Time,Psi,DtPsi,...) =: DtPsi or dPsi/dt = -i H|Psi> =         *
! * Func(Time,Psi,DtPsi,...) =: DtPsi, depending on flag "StdForm".    *
! * All computations are performed with double precision.              *
! *                                                                    *
! * Input parameters:                                                  *
! *   Psi:       The (complex) initial-value vector.                   *
! *   DtPsi:     Action of Hamiltonian on the initial-value vector, i. *
! *              e. H|Psi> or -i H|Psi>, depending on "StdForm".       *
! *   PsiDim     Length of Psi and DtPsi vectors.                      *
! *   Noffd:     No. non-zero off-diagonal Hamiltonian matrix          *
! *              elements. This shouldn't really be here, but this     *
! *              number has to be passed to the matrix-vector          *
! *              multiplication routine.
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
! *   OldErrCri: If true, the old (i.e. standard) error criterion is   *
! *              used, otherwise an improved one is taken.             *
! *   Steps:     Number of steps made so far (Steps is passed to       *
! *              "WriteStep").                                         *
! *   Krylov:    Matrix with minimum size PsiDim*(IntOrder-1) with     *
! *              columns containing the Krylov vectors H^2|psi> to     *
! *              H^TrueOrder|psi>. The Krylov vectors are needed only  *
! *              when Restart is false.                                *
! *                                                                    *
! * Output parameters:                                                 *
! *   Psi:       Propagated Psi.                                       *
! *   DtPsi:     Normalised first Krylov vector, i. e. (H-<H>)|Psi>.   *
! *   Krylov:    Matrix with minimum size PsiDim*(IntOrder-1) the      *
! *              columns of which contain the Krylov vectors that      *
! *              correspond to the vectors H^2|psi> to                 *
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
! *              3: can't diagonalise matrix,                          *
! *              4: stepsize too large for restart,                    *
! *              5: can't solve linear equation.                       *
! *                                                                    *
! * Other parameters:                                                  *
! *   CData:     Complex*16 data needed in the Func routine.           *
! *   RData:     Real*8 data needed in the Func routine.               *
! *   IData:     Integer data needed in the Func routine.              *
! *   LData:     Logical data needed in the Func routine.              *
! *                                                                    *
! * External routines:                                                 *
! *   Func:      Computes the action of the Hamiltonian on Psi. Called *
! *              as Func(Time,Psi,DtPsi,CData,RData,IData,LData).      *
! *   eigchess:  Computes the eigenvalues and -vectors of a complex    *
! *              upper Hessenberg matrix. Called as:                   *
! *              eigchess(LeadingDimension,Dimension,Matrix,           *
! *              macheps,Eigenvalues,Eigenvectors,ErrorCode).          *
! *   LinEqz:    Solves the complex linear equation Ax = b. Called as: *
! *              LinEqz(LeadingDimension,Dimension,A,b,x,ErrorCode).   *
! *   WriteStep: Writes the stepsize and error to a file.              *
! *              Called as WriteStep(Steps,Order,Stepsize,Error).      *
! *                                                                    *
! * V6.0 MB                                                            *
! * V7.0 GW addition of CData,RData,IData,LData arrays                 *
! * V7.0 MB improved error criterion                                   *
! * V8.0 MB routine can handle H|psi> and -i H|psi>                    *
! *                                                                    *
! **********************************************************************


      Subroutine CSILStep(Psi,DtPsi,PsiDim,Noffd,IntPeriod,IntOrder,&
                          TolError,Relax,Restart,StdForm,OldErrCri,&
                          Steps,Stepsize,TrueOrder,ErrorCode,Time,&
                          macheps,Func,Hessenberg,EigVec,Krylov,AuxMat)

!        use wrintegrat, only: writestep
      use lineqmod, only: lineqz
      use eigchessmod
      
      Implicit None

      Real(dop), Parameter     :: RelativeMinStep = 1.0e-9_dop,Tiny = 1.0e-18_dop

      Integer, Parameter       :: MinOrder = 2,MaxOrder = 50

      Logical(kind=4), intent(in)                                :: Relax,Restart,StdForm,OldErrCri
      Integer, intent(inout)                                     :: Steps
      Integer, intent(out)                                       :: TrueOrder,ErrorCode
      Integer, intent(in)                                        :: PsiDim,IntOrder
      Integer*8, intent(in)                                      :: Noffd
      Real(dop), intent(in)                                      :: IntPeriod,TolError,Time,macheps
      Real(dop), intent(out)                                     :: Stepsize
      Complex(dop), dimension(PsiDim), intent(inout)             :: Psi,DtPsi
      Complex(dop), dimension(PsiDim,2:IntOrder), intent(inout)  :: Krylov
      Complex(dop), dimension(0:IntOrder,0:IntOrder), intent(out):: Hessenberg,EigVec,AuxMat
      External                                                   :: Func

      Integer                                   :: NextOrder,D,P,Q
      Integer, save                             :: OldOrder
      Real(dop), save                           :: OldError,OldStepsize
      Real(dop)                                 :: MinStepsize,Error,FormerError,SubDiag,&
                                                   Norm,Inverse
      Complex(dop)                              :: PreFactor,CTemp,CInverse
      Complex(dop), dimension(0:MaxOrder)       :: Aux1,Aux3,OldAux1
      Complex(dop), dimension(0:MaxOrder), save :: EigVal,Aux2


! --- CHECK INTEGRATION ORDER ---

      If ((IntOrder .LT. MinOrder) .Or. (IntOrder .GT. MaxOrder)) Then
         ErrorCode = 1
         Return
      EndIf

! --- INITIALIZE VARIABLES ---

      Stepsize = IntPeriod
      MinStepsize = Abs(RelativeMinStep*IntPeriod)
      ErrorCode = 0
      If (Relax) Then
         PreFactor = (-1.0_dop,0.0_dop)
      Else
         PreFactor = (0.0_dop,-1.0_dop)
      EndIf

! --- RETURN RESULT IF DIFFERENTIAL EQUATION IS NOT A SYSTEM ---

      If (PsiDim .Eq. 1) Then
         TrueOrder = 1
         Error = 0.0_dop
!         Call WriteStep(Steps,TrueOrder,Stepsize,Error,Time)
         If (StdForm) Then
            Psi(1) = Exp(StepSize*DtPsi(1)/Psi(1))*Psi(1)
         Else
            Psi(1) = Exp(PreFactor*StepSize*DtPsi(1)/Psi(1))*Psi(1)
         EndIf
         Return
      EndIf

! --- COMPUTE NORM OF INITIAL WAVEFUNCTION ---

      Norm = 0.0_dop
      Do D = 1,PsiDim
         Norm = Norm+dble(DConjg(Psi(D))*Psi(D))
      EndDo
      Norm = DSqrt(Norm)

! --- NORMALISE INITIAL WAVEFUNCTIONS ---

      Inverse = 1.0_dop/Norm
      If (Restart) Then
         Do D = 1,PsiDim
            Psi(D) = Inverse*Psi(D)
            DtPsi(D) = Inverse*DtPsi(D)
         EndDo
      Else
         Do D = 1,PsiDim
            Psi(D) = Inverse*Psi(D)
         EndDo
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

      FormerError = 1.0_dop
      TrueOrder = 0
      Do Q = 0,IntOrder
         OldAux1(Q) = (0.0_dop,0.0_dop)
         Do P = 0,IntOrder
            Hessenberg(P,Q) = (0.0_dop,0.0_dop)
         EndDo
      EndDo
      
! --- BUILD UP THE KRYLOV SPACE ---

 100  Continue
         
! --- CHECK WHETHER STEPSIZE IS TOO SMALL ---

      If (Abs(Stepsize) .LT. Abs(MinStepsize)) Then
         ErrorCode = 2
         Return
      EndIf

! --- LOOP OVER ALL LOWER KRYLOV VECTORS ---

      NextOrder = TrueOrder+1
      Do P = 0,TrueOrder

!    --- COMPUTE ELEMENT OF HESSENBERG MATRIX ---

! Note that the first two Krylov vectors are stored in "Psi" and
! "DtPsi", respectively, whereas the higher vectors reside in the
! corresponding columns of the "Krylov" matrix.

         CTemp = (0.0_dop,0.0_dop)
         If (P .Eq. 0) Then
            If (TrueOrder .Eq. 0) Then
               Do D = 1,PsiDim
                  CTemp = CTemp+DConjg(Psi(D))*DtPsi(D)
               EndDo
            Else
               Do D = 1,PsiDim
                  CTemp = CTemp+DConjg(Psi(D))*Krylov(D,NextOrder)
               EndDo
            EndIf
         ElseIf (P .Eq. 1) Then
            Do D = 1,PsiDim
               CTemp = CTemp+DConjg(DtPsi(D))*Krylov(D,NextOrder)
            EndDo
         Else
            Do D = 1,PsiDim
               CTemp = CTemp+DConjg(Krylov(D,P))*Krylov(D,NextOrder)
            EndDo
         EndIf
         Hessenberg(P,TrueOrder) = CTemp

!    --- ORTHOGONALISE CURRENT KRYLOV VECTOR ON LOWER ONES ---

         If (P .Eq. 0) Then
            If (TrueOrder .Eq. 0) Then
               Do D = 1,PsiDim
                  DtPsi(D) = DtPsi(D)-CTemp*Psi(D)
               EndDo
            Else
               Do D = 1,PsiDim
                  Krylov(D,NextOrder) = Krylov(D,NextOrder)-CTemp*Psi(D)
               EndDo
            EndIf
         ElseIf (P .Eq. 1) Then
            Do D = 1,PsiDim
               Krylov(D,NextOrder) = Krylov(D,NextOrder)-CTemp*DtPsi(D)
            EndDo
         Else
            Do D = 1,PsiDim
               Krylov(D,NextOrder) = Krylov(D,NextOrder)&
                                     -CTemp*Krylov(D,P)
            EndDo
         EndIf
      EndDo

! --- DETERMINE SUBDIAGONAL ELEMENT OF HESSENBERG MATRIX ---

      SubDiag = 0.0_dop
      If (TrueOrder .Eq. 0) Then
         Do D = 1,PsiDim
            SubDiag = SubDiag+Dble(DConjg(DtPsi(D))*DtPsi(D))
         EndDo
      Else
         Do D = 1,PsiDim
            SubDiag = SubDiag+Dble(&
                      DConjg(Krylov(D,NextOrder))*Krylov(D,NextOrder))
         EndDo
      EndIf
      SubDiag = DSqrt(SubDiag)
      If (StdForm) Then
         Hessenberg(NextOrder,TrueOrder) = PreFactor*DCmplx(SubDiag)
      Else
         Hessenberg(NextOrder,TrueOrder) = DCmplx(SubDiag)
      EndIf

! --- NORMALISE KRYLOV VECTOR ---

      If (SubDiag .LT. Tiny) Then
         SubDiag = 0.0_dop
         CInverse = (0.0_dop,0.0_dop)
      Else
         CInverse = (1.0_dop,0.0_dop)/Hessenberg(NextOrder,TrueOrder)
      EndIf
      If (TrueOrder .Eq. 0) Then
         Do D = 1,PsiDim
            DtPsi(D) = CInverse*DtPsi(D)
         EndDo
      Else
         Do D = 1,PsiDim
            Krylov(D,NextOrder) = CInverse*Krylov(D,NextOrder)
         EndDo
      EndIf
      TrueOrder = NextOrder

! --- ESTIMATE ERROR (OLD ERROR CRITERION)---

      FormerError = Abs(FormerError*Stepsize*SubDiag/Dble(TrueOrder))
      If (OldErrCri) Then

!    --- COMPUTE AND OUTPUT ERROR ---

         Error = FormerError
!         Call WriteStep(Steps,TrueOrder,Stepsize,Error,Time)

!    --- CONTINUE UNTIL CONVERGED ---

         If (((Error .GT. TolError) .And. (TrueOrder .LT. IntOrder))&
             .Or. (TrueOrder .LT. MinOrder)) Then

!       --- COMPUTE NEXT KRYLOV VECTOR ---

            If (TrueOrder .Eq. 1) Then
               Call Func(Time,PsiDim,Noffd,DtPsi,Krylov(1,2))
            Else
               Call Func(Time,PsiDim,Noffd,Krylov(1,TrueOrder),Krylov(1,TrueOrder+1))
            EndIf

            Goto 100
         EndIf
      EndIf

! --- ESTIMATE NEXT MATRIX ELEMENTS ---

      Hessenberg(TrueOrder-1,TrueOrder) =&
         Hessenberg(TrueOrder,TrueOrder-1)
      Hessenberg(TrueOrder,TrueOrder) =&
         Hessenberg(TrueOrder-1,TrueOrder-1)

! --- SAVE HESSENBERG MATRIX ---

      If (StdForm) Then
         CInverse = (1.0_dop,0.0_dop)/PreFactor
         Do Q = 0,TrueOrder
            Do P = 0,TrueOrder
               AuxMat(P,Q) = CInverse*Hessenberg(P,Q)
            EndDo
         EndDo
      Else
         Do Q = 0,TrueOrder
            Do P = 0,TrueOrder
               AuxMat(P,Q) = Hessenberg(P,Q)
            EndDo
         EndDo
      EndIf

! --- DIAGONALISE HESSENBERG MATRIX ---

      Call eigchess(IntOrder+1,TrueOrder+1,AuxMat,macheps,EigVal,&
                    EigVec,ErrorCode)
      If (ErrorCode .NE. 0) Then
         ErrorCode = 3
         Return
      EndIf

! --- SAVE EIGENVECTORS ---

      Do Q = 0,TrueOrder
         Do P = 0,TrueOrder
            AuxMat(P,Q) = EigVec(P,Q)
         EndDo
      EndDo

! --- SET-UP PSI-VECTOR IN KRYLOV SPACE ---

      Aux1(0) = (1.0_dop,0.0_dop)
      Do P = 1,TrueOrder
         Aux1(P) = (0.0_dop,0.0_dop)
      EndDo

! --- MULTIPLY INVERSE EIGENVECTOR MATRIX WITH PSI-VECTOR ---

      Call LinEqz(IntOrder+1,TrueOrder+1,AuxMat,Aux1,Aux2,ErrorCode)
      If (ErrorCode .NE. 0) Then
         ErrorCode = 5
         Return
      EndIf

! --- PROPAGATE (OR RELAX) THE PSI-VECTOR IN KRYLOV-SPACE ---

      Do P = 0,TrueOrder
         Aux3(P) = Exp(PreFactor*EigVal(P)*Stepsize)*Aux2(P)
      EndDo

! --- MULTIPLY EIGENVECTOR MATRIX WITH PSI-VECTOR ---

      Do P = 0,TrueOrder
         Aux1(P) = (0.0_dop,0.0_dop)
      EndDo
      Do Q = 0,TrueOrder
         Do P = 0,TrueOrder
            Aux1(P) = Aux1(P)+EigVec(P,Q)*Aux3(Q)
         EndDo
      EndDo

! --- ESTIMATE THE ERROR (NEW CRITERION)---

      If (.Not. OldErrCri) Then

!    --- COMPUTE AND OUTPUT ERROR ---

         Error = 0.0_dop
         If (TrueOrder .Eq. 1) Then
            Do P = 0,TrueOrder
               OldAux1(P) = Aux1(P)
            EndDo
         Else
            Do P = 0,TrueOrder
               Error = Error+Abs(OldAux1(P)-Aux1(P))**2
               OldAux1(P) = Aux1(P)
            EndDo
            Error = DSqrt(Error)
         EndIf
!         Call WriteStep(Steps,TrueOrder,Stepsize,Error,Time)

!    --- CONTINUE UNTIL CONVERGED ---

         If (((Error .GT. TolError) .And. (TrueOrder .LT. IntOrder))&
             .Or. (TrueOrder .LT. MinOrder)) Then

!       --- COMPUTE NEXT KRYLOV VECTOR ---

            If (TrueOrder .Eq. 1) Then
               Call Func(Time,PsiDim,Noffd,DtPsi,Krylov(1,2))
            Else
               Call Func(Time,PsiDim,Noffd,Krylov(1,TrueOrder),Krylov(1,TrueOrder+1))
            EndIf

            Goto 100
         EndIf
      EndIf

! --- SAVE SOME VALUES IN CASE THAT NEXT CALL IS NON-RESTART CALL ---

      OldOrder = TrueOrder
      OldStepsize = Stepsize
      OldError = Error

! --- DECREASE STEPSIZE IF LANCZOS DIDN'T CONVERGE ---
      
      If (Error .LE. TolError) Goto 300
      Stepsize = Min(Stepsize,&
         Stepsize*(TolError/FormerError)**(1.0_dop/Dble(TrueOrder)))
      OldStepsize = Stepsize
!      Call WriteStep(Steps,TrueOrder,Stepsize,TolError,Time)

! --- PROPAGATE (OR RELAX) THE PSI-VECTOR IN KRYLOV-SPACE ---

 200  Continue
      Do P = 0,TrueOrder
         Aux3(P) = Exp(PreFactor*EigVal(P)*Stepsize)*Aux2(P)
      EndDo

! --- MULTIPLY EIGENVECTOR MATRIX WITH PSI-VECTOR ---

      Do P = 0,TrueOrder
         Aux1(P) = (0.0_dop,0.0_dop)
      EndDo
      Do Q = 0,TrueOrder
         Do P = 0,TrueOrder
            Aux1(P) = Aux1(P)+EigVec(P,Q)*Aux3(Q)
         EndDo
      EndDo

! --- PROPAGATE PSI IN HILBERT SPACE ---

! Note again that the Krylov vectors number zero and one aren't stored
! in "Krylov" but in "Psi" and "DtPsi", respectively.

 300  Continue
      Do D = 1,PsiDim
         Psi(D) = Aux1(0)*Psi(D)+Aux1(1)*DtPsi(D)
      EndDo
      Do P = 2,TrueOrder
         Do D = 1,PsiDim
            Psi(D) = Psi(D)+Aux1(P)*Krylov(D,P)
         EndDo
      EndDo

! --- RE-ESTABLISH THE NORM OF PSI ---

      Do D = 1,PsiDim
         Psi(D) = Norm*Psi(D)
      EndDo

      Return
      End Subroutine CSILStep


! **********************************************************************
! *                                                                    *
! *                      SUBROUTINE CSILERRORMSG                       *
! *                                                                    *
! * Generates for a given error number returned by "CSILStep/          *
! * CSILStep2" a corresponding error message.                          *
! *                                                                    *
! * Input parameters:                                                  *
! *   Error: Error code returned by CSILStep/CSILStep2.                *
! *                                                                    *
! * Output parameters:                                                 *
! *   Msg:   Error message.                                            *
! *                                                                    *
! * V7.0 MB                                                            *
! *                                                                    *
! **********************************************************************

      Subroutine CSILErrorMsg (Error,Msg)

      Implicit None

      Integer, intent(in)           :: Error
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
      ElseIf (Error .Eq. 5) Then
         Msg = 'Cannot solve linear equation'
      Else
         Msg = 'Unknown error occurred'
      EndIf

      Return
      End Subroutine CSILErrorMsg


      end module csillib

