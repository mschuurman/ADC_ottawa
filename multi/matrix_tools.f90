!
!  Simple wrappers around matmul() and related routines to enable
!  OpenMP parallel execution. It looks like I can't trust gfortran
!  to do this for me :(
!
  module matrix_tools
    use accuracy
    use timer
    implicit none
    private
    public mt_matmul
    !
    integer(ik), parameter :: stripe     = 48  ! Size of a parallel stripe
    integer(ik), parameter :: min_blocks = 4   ! Smallest number of blocks to consider parallel execution
    !
    interface mt_matmul
      module procedure mt_matmul_crcr  ! complex(rk)  by complex(rk)
      module procedure mt_matmul_crrr  ! complex(rk)  by real(rk)
      module procedure mt_matmul_rrcr  ! real(rk)     by complex(rk)
      module procedure mt_matmul_rrrr  ! real(rk)     by real(rk)
!*qd  module procedure mt_matmul_cqcq  ! complex(xrk) by complex(xrk)
!*qd  module procedure mt_matmul_cqrq  ! complex(xrk) by real(xrk)
!*qd  module procedure mt_matmul_rqcq  ! reak(xrk)    by complex(xrk)
!*qd  module procedure mt_matmul_rqrq  ! real(xrk)    by real(xrk)
!*qd  module procedure mt_matmul_cqcr  ! complex(xrk) by complex(rk)
!*qd  module procedure mt_matmul_cqrr  ! complex(xrk) by real(rk)
!*qd  module procedure mt_matmul_rqcr  ! reak(xrk)    by complex(rk)
!*qd  module procedure mt_matmul_rqrr  ! real(xrk)    by real(rk)
!*qd  module procedure mt_matmul_crcq  ! complex(rk)  by complex(xrk)
!*qd  module procedure mt_matmul_crrq  ! complex(rk)  by real(xrk)
!*qd  module procedure mt_matmul_rrcq  ! reak(rk)     by complex(xrk)
!*qd  module procedure mt_matmul_rrrq  ! real(rk)     by real(xrk)
    end interface mt_matmul
    !
    contains
    !
    !  Calculate c=matmul(a,b), possibly using OpenMP parallel execution
    !  All functions reuse the same code; the only difference is in the
    !  arguments declaration
    !
    !  === standard by standard ===
    !
    function mt_matmul_crcr(a,b,serial) result(c)
      complex(rk), intent(in) :: a(:,:) !
      complex(rk), intent(in) :: b(:,:) !
      complex(rk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_crcr
    !
    function mt_matmul_crrr(a,b,serial) result(c)
      complex(rk), intent(in) :: a(:,:) !
      real   (rk), intent(in) :: b(:,:) !
      complex(rk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_crrr
    !
    function mt_matmul_rrcr(a,b,serial) result(c)
      real   (rk), intent(in) :: a(:,:) !
      complex(rk), intent(in) :: b(:,:) !
      complex(rk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_rrcr
    !
    function mt_matmul_rrrr(a,b,serial) result(c)
      real   (rk), intent(in) :: a(:,:) !
      real   (rk), intent(in) :: b(:,:) !
      real   (rk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_rrrr
    !
    !  === quad by quad ===
    !
    function mt_matmul_cqcq(a,b,serial) result(c)
      complex(xrk), intent(in) :: a(:,:) !
      complex(xrk), intent(in) :: b(:,:) !
      complex(xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_cqcq
    !
    function mt_matmul_cqrq(a,b,serial) result(c)
      complex(xrk), intent(in) :: a(:,:) !
      real   (xrk), intent(in) :: b(:,:) !
      complex(xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_cqrq
    !
    function mt_matmul_rqcq(a,b,serial) result(c)
      real   (xrk), intent(in) :: a(:,:) !
      complex(xrk), intent(in) :: b(:,:) !
      complex(xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_rqcq
    !
    function mt_matmul_rqrq(a,b,serial) result(c)
      real   (xrk), intent(in) :: a(:,:) !
      real   (xrk), intent(in) :: b(:,:) !
      real   (xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_rqrq
    !
    !  === quad by standard ===
    !
    function mt_matmul_cqcr(a,b,serial) result(c)
      complex(xrk), intent(in) :: a(:,:) !
      complex(rk),  intent(in) :: b(:,:) !
      complex(xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_cqcr
    !
    function mt_matmul_cqrr(a,b,serial) result(c)
      complex(xrk), intent(in) :: a(:,:) !
      real   (rk),  intent(in) :: b(:,:) !
      complex(xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_cqrr
    !
    function mt_matmul_rqcr(a,b,serial) result(c)
      real   (xrk), intent(in) :: a(:,:) !
      complex(rk),  intent(in) :: b(:,:) !
      complex(xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_rqcr
    !
    function mt_matmul_rqrr(a,b,serial) result(c)
      real   (xrk), intent(in) :: a(:,:) !
      real   (rk),  intent(in) :: b(:,:) !
      real   (xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_rqrr
    !
    !  === standard by quad ===
    !
    function mt_matmul_crcq(a,b,serial) result(c)
      complex(rk),  intent(in) :: a(:,:) !
      complex(xrk), intent(in) :: b(:,:) !
      complex(xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_crcq
    !
    function mt_matmul_crrq(a,b,serial) result(c)
      complex(rk),  intent(in) :: a(:,:) !
      real   (xrk), intent(in) :: b(:,:) !
      complex(xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_crrq
    !
    function mt_matmul_rrcq(a,b,serial) result(c)
      real   (rk),  intent(in) :: a(:,:) !
      complex(xrk), intent(in) :: b(:,:) !
      complex(xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_rrcq
    !
    function mt_matmul_rrrq(a,b,serial) result(c)
      real   (rk),  intent(in) :: a(:,:) !
      real   (xrk), intent(in) :: b(:,:) !
      real   (xrk)             :: c(size(a,dim=1),size(b,dim=2))
      !
      include 'matrix_tools_matmul_common.f90'
    end function mt_matmul_rrrq
  end module matrix_tools
