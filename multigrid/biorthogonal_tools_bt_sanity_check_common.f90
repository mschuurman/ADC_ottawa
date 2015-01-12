!   subroutine bt_sanity_check_rk(evec,msg,sphalf,null_count)
!     complex(rk), intent(in)           :: evec(:,:,:)  ! Eigenvectors
!     character(len=*), intent(in)      :: msg          ! Name of the caller for error messages
!     real(rk), intent(in), optional    :: sphalf(:,:)  ! Metric tensor; if missing assume unity.
!     integer(ik), intent(in), optional :: null_count   ! Number of null orbitals expected
      !
      !  Keep code type-generic below this point!
      !
      real(kind(evec))                 :: eps
      complex(kind(evec)), allocatable :: ovt(:,:)
      complex(kind(evec))              :: ovr
      integer(ik)              :: imo, jmo, nev, lev, alloc, inull, nnull
      !
      call TimerStart('Bi-orthogonality check')
      lev = size(evec,dim=1)
      nev = size(evec,dim=2)
      if (nev>lev .or. 2/=size(evec,dim=3)) call stop('biorthogonal_tools%bt_sanity_check - bad input')
      !
      allocate (ovt(nev,nev),stat=alloc)
      if (alloc/=0) call stop('biorthogonal_tools%bt_sanity_check - no memory')
      !
      if (present(sphalf)) then
        ovt = mt_matmul(transpose(mt_matmul(sphalf,evec(:,:,1))),mt_matmul(sphalf,evec(:,:,2)))
      else 
        ovt = mt_matmul(transpose(evec(:,:,1)),evec(:,:,2))
      end if
      !
      eps = max(real(1e-14,kind(evec)), 1e4*spacing(maxval(abs(evec))))
      inull = 0
      check_diagonal: do imo=1,nev
        ovr = ovt(imo,imo)
        if (abs(ovr-1)<=eps) cycle check_diagonal
        if (abs(ovr)<=1e3*eps) then
          inull = inull + 1
          cycle check_diagonal
        end if
        write (out,"('In ',a,' norm of eigenvector ',i5,' is in error by',2(1x,g12.5))") trim(msg), imo, ovr-1
        if (abs(ovr-1)>1e3*eps) then
          write (out,"('This is too much. eps = ',g12.5)") eps
          call stop('biorthogonal_tools%bt_sanity_check - vectors are not normalized')
        end if
      end do check_diagonal
      nnull = 0
      if (present(null_count)) nnull = null_count
      if (inull/=nnull) then
        write (out,"(1x,i0,' null vectors were expected in ',a,', but ',i0,' were found.')") &
          nnull, trim(msg), inull
      end if
      !
      check_right: do jmo=1,nev
        check_left: do imo=1,nev
          if (imo==jmo) cycle check_left
          ovr = ovt(imo,jmo)
          !
          if (abs(ovr)>eps) then
            write (out,"('In ',a,' overlap of left and right eigenvectors ',i5,' and ',i5,' is ',2(1x,g12.5))") &
                   trim(msg), imo, jmo, ovr
            if (abs(ovr)>1e3_xrk*eps) then    
              write (out,"('This is too much. eps = ',g12.5)") eps
              call stop('biorthogonal_tools%bt_sanity_check - vectors are not biorthogonal')
            end if
          end if
        end do check_left
      end do check_right
      deallocate (ovt)
      call TimerStop('Bi-orthogonality check')
!   end subroutine bt_sanity_check_rk
