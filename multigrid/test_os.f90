 program test
   use accuracy
   use math
   use gamess_internal
   use import_gamess
   !
   type (gam_operator_data) :: op
   real(rk)                 :: rho, t, gn(0:1000), gn_ref, rele
   integer(ik)              :: n
   !
   call accuracyInitialize
   op%op_name = '3C R/(R**2+A**2)'
   read_loop: do
     read(input,*,end=100) n, rho, t, op%imag_rc, gn_ref
     call os_basic_integral(op,rho,t,n,gn(0:n))
     rele = abs((gn(n)-gn_ref)/gn_ref)
     write (out,"(1x,i4,4(1x,g25.18),1x,f25.16,2(1x,g25.18))") n, rho, t, op%imag_rc, gn(n)-gn_ref, rele, gn(n), gn_ref
   end do read_loop
   100 continue
 end program test
