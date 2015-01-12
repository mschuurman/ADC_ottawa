!
!  Limited ARPACK interface for ci_tools.f90
!  The only purpose of this module is to avoid stupid type errors
!
module arpack
  use accuracy
  implicit none

  private
  public arpack_znaupd, arpack_zneupd

  interface arpack_znaupd
    subroutine znaupd(ido,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info)
      use accuracy
      implicit none
      integer           :: ido, info, ldv, lworkl, n, ncv, nev
      real(kind=drk)    :: tol
      character(len=1)  :: bmat
      character(len=2)  :: which
      integer           :: iparam(11), ipntr(14)
      complex(kind=drk) :: resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
      real(kind=drk)    :: rwork(ncv)
    end subroutine znaupd
  end interface arpack_znaupd

  interface arpack_zneupd
    subroutine zneupd(rvec,howmny,select,d,z,ldz,sigma,workev,bmat,n,which,nev,tol,resid,ncv,v,ldv,iparam, &
                      ipntr,workd,workl,lworkl,rwork,info)
      use accuracy
      implicit none
      character(len=1)  :: bmat, howmny
      character(len=2)  :: which
      logical           :: rvec
      integer           :: info, ldz, ldv, lworkl, n, ncv, nev
      complex(kind=drk) :: sigma
      real(kind=drk)    :: tol
      integer           :: iparam(11), ipntr(14)
      logical           :: select(ncv)
      real(kind=drk)    :: rwork(ncv)
      complex(kind=drk) :: d(nev), resid(n), v(ldv,ncv), z(ldz,nev), workd(3*n), workl(lworkl), workev(2*ncv)
    end subroutine zneupd
  end interface arpack_zneupd

end module arpack
