!
!  Simple OpenMP-parallel vector manipulation routines
!
module vector
  use accuracy
  use timer
  private
  public vector_copy, vector_add_to, vector_scale_and_add_to
  public vector_sum, vector_dot_product, vector_zero
  public vector_scale, vector_copy_and_scale, vector_merge_to
  public vector_minmax, vector_absmax, vector_clip_below
  public vector_max_change
  !
  contains
  !
  !  Externally visible subroutines
  !
  subroutine vector_copy(dst,src)
    real(rk), intent(out) :: dst(:)
    real(rk), intent(in)  :: src(:)
    !
    integer(ik) :: ipt, npt
    !
    npt = min(size(src),size(dst))
    !
    !$omp parallel do private (ipt)
    copy_vector: do ipt=1,npt
      dst(ipt) = src(ipt)
    end do copy_vector
    !$omp end parallel do
  end subroutine vector_copy
  !
  subroutine vector_add_to(dst,src)
    real(rk), intent(inout) :: dst(:)
    real(rk), intent(in)    :: src(:)
    !
    integer(ik) :: ipt, npt
    !
    npt = min(size(src),size(dst))
    !
    !$omp parallel do private (ipt)
    add_vectors: do ipt=1,npt
      dst(ipt) = dst(ipt) + src(ipt)
    end do add_vectors
    !$omp end parallel do
  end subroutine vector_add_to
  !
  subroutine vector_scale_and_add_to(dst,scale,src)
    real(rk), intent(inout) :: dst(:)
    real(rk), intent(in)    :: scale
    real(rk), intent(in)    :: src(:)
    !
    integer(ik) :: ipt, npt
    !
    npt = min(size(src),size(dst))
    !
    !$omp parallel do private (ipt)
    add_vectors: do ipt=1,npt
      dst(ipt) = dst(ipt) + scale * src(ipt)
    end do add_vectors
    !$omp end parallel do
  end subroutine vector_scale_and_add_to
  !
  subroutine vector_merge_to(dst,src,wgt_dst,wgt_src)
    real(rk), intent(inout) :: dst(:)
    real(rk), intent(in)    :: src(:)
    real(rk), intent(in)    :: wgt_src, wgt_dst
    !
    integer(ik) :: ipt, npt
    !
    npt = min(size(src),size(dst))
    !
    !$omp parallel do private (ipt)
    merge_vectors: do ipt=1,npt
      dst(ipt) = wgt_dst*dst(ipt) + wgt_src*src(ipt)
    end do merge_vectors
    !$omp end parallel do
  end subroutine vector_merge_to
  !
  subroutine vector_copy_and_scale(dst,src,a,b)
    real(rk), intent(out) :: dst(:)
    real(rk), intent(in)  :: src(:)
    real(rk), intent(in)  :: a, b ! dst = a + b*src
    !
    integer(ik) :: ipt, npt
    !
    npt = min(size(src),size(dst))
    !
    !$omp parallel do private (ipt)
    copy_vector: do ipt=1,npt
      dst(ipt) = a + b*src(ipt)
    end do copy_vector
    !$omp end parallel do
  end subroutine vector_copy_and_scale
  !
  subroutine vector_scale(dst,a,b)
    real(rk), intent(inout) :: dst(:)
    real(rk), intent(in)    :: a, b ! dst = a + b*dst
    !
    integer(ik) :: ipt, npt
    !
    npt = size(dst)
    !
    !$omp parallel do private (ipt)
    scale_vector: do ipt=1,npt
      dst(ipt) = a + b*dst(ipt)
    end do scale_vector
    !$omp end parallel do
  end subroutine vector_scale
  !
  subroutine vector_zero(dst)
    real(rk), intent(out) :: dst(:)
    !
    integer(ik) :: ipt, npt
    !
    npt = size(dst)
    !
    !$omp parallel do private (ipt)
    zero_vector: do ipt=1,npt
      dst(ipt) = 0
    end do zero_vector
    !$omp end parallel do
  end subroutine vector_zero
  !
  function vector_dot_product(src1,src2) result(s)
    real(rk), intent(in) :: src1(:), src2(:)
    real(rk)             :: s
    !
    integer(ik) :: ipt, npt
    !
    npt = min(size(src1),size(src2))
    s   = 0
    !
    !$omp parallel do private (ipt) reduction(+:s)
    sum_vector: do ipt=1,npt
      s = s + src1(ipt)*src2(ipt)
    end do sum_vector
    !$omp end parallel do
  end function vector_dot_product
  !
  function vector_sum(src) result(s)
    real(rk), intent(in) :: src(:)
    real(rk)             :: s
    !
    integer(ik) :: ipt, npt
    !
    npt = size(src)
    s   = 0
    !
    !$omp parallel do private (ipt) reduction(+:s)
    sum_vector: do ipt=1,npt
      s = s + src(ipt)
    end do sum_vector
    !$omp end parallel do
  end function vector_sum
  !
  function vector_max_change(src1,src2) result(d)
    real(rk), intent(in) :: src1(:), src2(:)
    real(rk)             :: d
    !
    integer(ik) :: ipt, npt
    !
    npt = min(size(src1),size(src2))
    d   = 0
    !
    !$omp parallel do private (ipt) reduction(max:d)
    diff_vector: do ipt=1,npt
      d = max(d,abs(src1(ipt)-src2(ipt)))
    end do diff_vector
    !$omp end parallel do
  end function vector_max_change
  !
  subroutine vector_minmax(src,min_src,max_src)
    real(rk), intent(in)  :: src(:)
    real(rk), intent(out) :: min_src, max_src
    !
    integer(ik) :: ipt, npt
    !
    npt     = size(src)
    min_src = src(1)
    max_src = src(1)
    !
    !$omp parallel do private(ipt) reduction(min:min_src) reduction(max:max_src)
    scan_vector: do ipt=2,npt
      min_src = min(min_src,src(ipt))
      max_src = max(max_src,src(ipt))
    end do scan_vector
    !$omp end parallel do
  end subroutine vector_minmax
  !
  function vector_absmax(src) result (am)
    real(rk), intent(in)  :: src(:)
    real(rk)              :: am
    !
    real(rk)    :: min_src, max_src
    !
    call vector_minmax(src,min_src,max_src)
    am = max(abs(min_src),abs(max_src))
  end function vector_absmax
  !
  subroutine vector_clip_below(dst,v_clip,n_clip)
    real(rk), intent(inout)  :: dst(:)
    real(rk), intent(in)     :: v_clip ! Max. allowed value
    integer(ik), intent(out) :: n_clip ! Number of values clipped
    !
    integer(ik) :: ipt, npt
    !
    npt    = size(dst)
    n_clip = 0
    !
    !$omp parallel do private (ipt) reduction(+:n_clip)
    clip_vector: do ipt=1,npt
      if (dst(ipt)<v_clip) then
        dst(ipt) = v_clip
        n_clip   = n_clip + 1
      end if
    end do clip_vector
    !$omp end parallel do
  end subroutine vector_clip_below
  !
end module vector
