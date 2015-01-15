!
!  Mostly to be used for debugging - produces a reasonably tidy printout of
!  a matrix.
!
module printing
  use accuracy
  implicit none
  private
  public print_matrix

  interface print_matrix
     module procedure print_complex_matrix
     module procedure print_real_matrix
!*qd module procedure print_quad_complex_matrix
!*qd module procedure print_quad_matrix
  end interface ! print_matrix

  integer(ik), parameter :: page_width = 132_ik ! Max. width of the printout

  contains

  subroutine print_complex_matrix(m,w,f)
    complex(rk), intent(in)      :: m(:,:) ! Matrix to be printed
    integer(ik), intent(in)      :: w      ! Field width
    character(len=*), intent(in) :: f      ! Format to use for number conversion
    !
    integer(ik)        :: nr, nc ! number of rows and columns in the matrix
    integer(ik)        :: cpp    ! columns per page of the output
    integer(ik)        :: cs, ce ! starting and ending columns
    integer(ik)        :: cc, rc ! column and row being processed now
    character(len=200) :: fmt_head ! Format for the header
    character(len=200) :: fmt_line ! Format for the data line
    !
    nr = size(m,dim=1)
    nc = size(m,dim=2)
    !
    if (w<6) stop 'printing%print_complex_matrix - logic error 1'
    cpp = (page_width-6) / (3+2*w)
    if (cpp<1) stop 'printing%print_complex_matrix - logic error 2'
    !
    write (fmt_head,"('(6x,',i0,'(4x,i6,',i0,'x))')") cpp, 3+2*w-10
    write (fmt_line,"('(i6,',i0,'(2x,',a,',1x,',a,'))')") cpp, f, f
    !
    pages: do cs=1,nc,cpp
      ce = min(nc,cs+cpp-1)
      write (out,fmt=trim(fmt_head)) (cc,cc=cs,ce)
      rows: do rc=1,nr
        write (out,fmt=trim(fmt_line)) rc,(m(rc,cc),cc=cs,ce)
      end do rows
      write (out,"()")
    end do pages
  end subroutine print_complex_matrix

  subroutine print_quad_complex_matrix(m,w,f)
    complex(xrk), intent(in)     :: m(:,:) ! Matrix to be printed
    integer(ik), intent(in)      :: w      ! Field width
    character(len=*), intent(in) :: f      ! Format to use for number conversion
    !
    integer(ik)        :: nr, nc   ! number of rows and columns in the matrix
    integer(ik)        :: cpp      ! columns per page of the output
    integer(ik)        :: cs, ce   ! starting and ending columns
    integer(ik)        :: cc, rc   ! column and row being processed now
    character(len=200) :: fmt_head ! Format for the header
    character(len=200) :: fmt_line ! Format for the data line
    !
    nr = size(m,dim=1)
    nc = size(m,dim=2)
    !
    if (w<6) stop 'printing%print_quad_complex_matrix - logic error 1'
    cpp = (page_width-6) / (3+2*w)
    if (cpp<1) stop 'printing%print_quad_complex_matrix - logic error 2'
    !
    write (fmt_head,"('(6x,',i0,'(4x,i6,',i0,'x))')") cpp, 3+2*w-10
    write (fmt_line,"('(i6,',i0,'(2x,',a,',1x,',a,'))')") cpp, f, f
    !
    pages: do cs=1,nc,cpp
      ce = min(nc,cs+cpp-1)
      write (out,fmt=trim(fmt_head)) (cc,cc=cs,ce)
      rows: do rc=1,nr
        write (out,fmt=trim(fmt_line)) rc,(m(rc,cc),cc=cs,ce)
      end do rows
      write (out,"()")
    end do pages
  end subroutine print_quad_complex_matrix

  subroutine print_real_matrix(m,w,f)
    real(rk), intent(in)         :: m(:,:) ! Matrix to be printed
    integer(ik), intent(in)      :: w      ! Field width
    character(len=*), intent(in) :: f      ! Format to use for number conversion
    !
    integer(ik)        :: nr, nc ! number of rows and columns in the matrix
    integer(ik)        :: cpp    ! columns per page of the output
    integer(ik)        :: cs, ce ! starting and ending columns
    integer(ik)        :: cc, rc ! column and row being processed now
    character(len=200) :: fmt_head ! Format for the header
    character(len=200) :: fmt_line ! Format for the data line
    !
    nr = size(m,dim=1)
    nc = size(m,dim=2)
    !
    if (w<8) stop 'printing%print_real_matrix - logic error 1'
    cpp = (page_width-6) / (3+w)
    if (cpp<1) stop 'printing%print_real_matrix - logic error 2'
    !
    write (fmt_head,"('(6x,1x,',i0,'(3x,i6,',i0,'x))')") cpp, 1+w-9
    write (fmt_line,"('(i6,1x,',i0,'(1x,',a,'))')") cpp, f
    !
    pages: do cs=1,nc,cpp
      ce = min(nc,cs+cpp-1)
      write (out,fmt=trim(fmt_head)) (cc,cc=cs,ce)
      rows: do rc=1,nr
        write (out,fmt=trim(fmt_line)) rc,(m(rc,cc),cc=cs,ce)
      end do rows
      write (out,"()")
    end do pages
  end subroutine print_real_matrix

  subroutine print_quad_matrix(m,w,f)
    real(xrk), intent(in)        :: m(:,:) ! Matrix to be printed
    integer(ik), intent(in)      :: w      ! Field width
    character(len=*), intent(in) :: f      ! Format to use for number conversion
    !
    integer(ik)        :: nr, nc ! number of rows and columns in the matrix
    integer(ik)        :: cpp    ! columns per page of the output
    integer(ik)        :: cs, ce ! starting and ending columns
    integer(ik)        :: cc, rc ! column and row being processed now
    character(len=200) :: fmt_head ! Format for the header
    character(len=200) :: fmt_line ! Format for the data line
    !
    nr = size(m,dim=1)
    nc = size(m,dim=2)
    !
    if (w<8) stop 'printing%print_real_matrix - logic error 1'
    cpp = (page_width-6) / (3+w)
    if (cpp<1) stop 'printing%print_real_matrix - logic error 2'
    !
    write (fmt_head,"('(6x,1x,',i0,'(3x,i6,',i0,'x))')") cpp, 1+w-9
    write (fmt_line,"('(i6,1x,',i0,'(1x,',a,'))')") cpp, f
    !
    pages: do cs=1,nc,cpp
      ce = min(nc,cs+cpp-1)
      write (out,fmt=trim(fmt_head)) (cc,cc=cs,ce)
      rows: do rc=1,nr
        write (out,fmt=trim(fmt_line)) rc,(m(rc,cc),cc=cs,ce)
      end do rows
      write (out,"()")
    end do pages
  end subroutine print_quad_matrix

end module printing
