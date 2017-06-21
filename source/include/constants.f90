module constants
  integer, parameter :: d=selected_real_kind(8)
  integer, parameter :: lng=selected_int_kind(16)
  

  real(d), parameter :: rzero=0._d
  real(d), parameter :: rone=1._d
  real(d), parameter :: pi=3.14159265358979_d
  complex(d), parameter :: ci=(0._d,1._d)
  complex(d), parameter :: czero=(0._d,0._d)
  complex(d), parameter :: cone=(1._d,0._d)
  
  integer, dimension(64), parameter :: mtrow=(/ 1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,4,&
       3,2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,7,8,5,6,3,4,1,2,8,7,6,5,4,3,2,1 /)
end module constants
