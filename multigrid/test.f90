 program test
   use accuracy
   use math
   integer(ik) :: imp, im, mult
   real(ark)   :: euler(3), mp, m, err
   complex(ark), allocatable :: mat(:,:)
   complex(ark) :: ref
   !
   call accuracyInitialize
   read(input,*) euler, mult
   allocate (mat(mult,mult))
   call MathYJMRotationMatrix(euler,mult,mat)
   do im=1,mult
     do imp=1,mult
       read (input,*) mp, m, ref
       err = abs(mat(imp,im)-ref)
       if (err>spacing(1e11_rk)) then
         write (out,"(1x,i4,1x,i4,3(2x,g22.14,1x,g22.14))") imp, im, mat(imp,im), ref, mat(imp,im)-ref
       end if
     end do
   end do
 end program test
!program test
!  use accuracy
!  use math
!  integer(ik) :: n
!  real(rk)    :: alp, bet, x, ref, pn, err
!  !
!  call accuracyInitialize
!  do
!    read *, n, alp, bet, x, ref
!    pn = MathJacobiPn(n,alp,bet,x)
!    err = abs(pn-ref)
!    if (err>300._rk*spacing(max(1._rk,abs(ref)))) then
!      write (out,"('n= ',i5,' alp= ',g22.14,' bet= ',g22.14,' x= ',g22.14,' pn= ',g22.14,' ref= ',g22.14,' err/1ulp= ',g22.14)") &
!             n, alp, bet, x, pn, ref, (pn-ref)/spacing(max(1._rk,abs(ref)))
!    end if
!  end do
!end program test
