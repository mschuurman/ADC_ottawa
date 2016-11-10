!   subroutine os_boys_table_real(n_max,t,fn)
!     integer(ik), intent(in) :: n_max       ! Max. desired order of the Boys' function
!     real(rk), intent(in)    :: t           ! Argument of the Boys' function
!     real(rk), intent(out)   :: fn(0:n_max) ! Table for the function values
!     !
      real(kind(t))           :: expt
      integer(ik)             :: n
      !
      !  Evaluation of the Boys' function in MathBoysF is accurate, but quite inefficient.
      !  If timing becomes an issue, table-based power-series expansion (as in OS) is
      !  probably the simplest option to improve the speed.
      !
      expt = exp(-t)
      fn(n_max) = MathBoysF(n_max,t)
      fn_boys_recursion: do n=n_max-1,0,-1
        fn(n) = (2*t*fn(n+1)+expt)/(2*n+1)
      end do fn_boys_recursion
!   end subroutine os_boys_table_real
