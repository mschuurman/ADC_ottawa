!   subroutine diis_extrapolate_rk(diis,iter,smat,rho,fmat,diis_fmat,diis_error,diis_eiej)
!     type(diis_state), intent(inout)  :: diis              ! DIIS state
!     integer(ik), intent(in)          :: iter              ! Current iteration; starts at 1
!     real(rk), intent(in)             :: smat(:,:)         ! Overlap matrix
!     complex(rk), intent(in)          :: rho (:,:)         ! Current density matrix
!     complex(rk), intent(inout)       :: fmat(:,:)         ! Input: Fock matrix corresponding to rho(:,:)
!                                                           ! Output: extrapolated Fock matrix
!     complex(rk), intent(inout)       :: diis_fmat (:,:,:) ! The appropriate-precision field in diis%
!     complex(rk), intent(inout)       :: diis_error(:,:,:) ! The appropriate-precision field in diis%
!     complex(rk), intent(inout)       :: diis_eiej (:,:)   ! The appropriate-precision field in diis%
      !
      !  The routine should be type-generic below this point
      !
      integer(ik)                      :: iv, alloc
      complex(kind(fmat)), allocatable :: diis_a(:,:)
      complex(kind(fmat)), allocatable :: diis_c(:,:)
      !
      call TimerStart('DIIS extrapolation')
      !
      !  (Re-)start DIIS iterations
      !
      if (iter==1) diis%nvec = 0
      if (diis%nvec==diis%max_nvec) diis%nvec = 0
      !
      !  Update DIIS Fock matrix and error vector tables
      !
      diis%nvec = diis%nvec + 1
      diis_fmat (:,:,diis%nvec) = fmat
      diis_error(:,:,diis%nvec) = mt_matmul(fmat,mt_matmul(transpose(rho),smat)) &
                                - mt_matmul(mt_matmul(smat,transpose(rho)),fmat)
      !
      !  Update error product table. We'll be going for minimizing the Euclidian
      !  norm of the error vectors.
      !
      !$omp parallel do default(none) shared(diis,diis_eiej,diis_error) private(iv)
      update_error_products: do iv=1,diis%nvec-1
        diis_eiej(iv,diis%nvec) = sum(conjg(diis_error(:,:,iv))*diis_error(:,:,diis%nvec))
        diis_eiej(diis%nvec,iv) = sum(conjg(diis_error(:,:,diis%nvec))*diis_error(:,:,iv))
      end do update_error_products
      !$omp end parallel do
      diis_eiej(diis%nvec,diis%nvec) = sum(conjg(diis_error(:,:,diis%nvec))*diis_error(:,:,diis%nvec))
      !
      if (verbose>0) then
        write (out,"(/'DIIS error norm =',2(1x,g16.8))") diis_eiej(diis%nvec,diis%nvec)
      end if
      !
      if (diis%nvec>1) then
        !
        !  We have more than one vector; go extrapolate
        !
        allocate (diis_a(diis%nvec+1,diis%nvec+1),diis_c(diis%nvec+1,1),stat=alloc)
        if (alloc/=0) then
          stop 'diis%diis_extrapolate - out of memory'
        end if
        !
        diis_a( :diis%nvec,:diis%nvec ) = diis_eiej(:diis%nvec,:diis%nvec)
        diis_a( :diis%nvec,diis%nvec+1) = 1
        diis_a(diis%nvec+1,:diis%nvec ) = 1
        diis_a(diis%nvec+1,diis%nvec+1) = 0
        diis_c( :diis%nvec,1) = 0
        diis_c(diis%nvec+1,1) = 1
        !
        call lapack_gelss(diis_a,diis_c)
        if (verbose>0) then
          write (out,"('DIIS extrapolation coefficients are:')")
          write (out,"((5(1x,g11.4,1x,g11.4,1x)))") diis_c(:diis%nvec,1)
          write (out,"()")
        end if
        !
        if (any(abs(diis_c(:diis%nvec,1))>diis%max_coeff)) then 
          !
          !  Failure; delete all vectors but the current one
          !
          diis_fmat (:,:,1) = diis_fmat (:,:,diis%nvec)
          diis_error(:,:,1) = diis_error(:,:,diis%nvec)
          diis_eiej (1,1)   = diis_eiej (diis%nvec,diis%nvec)
          diis%nvec = 1
          if (verbose>=0) then
            write (out,"('DIIS extrapolation failed, doing a simple iteration.')")
          end if
        else
          fmat = 0
          extrapolate_fmat: do iv=1,diis%nvec
            fmat = fmat + diis_c(iv,1) * diis_fmat(:,:,iv)
          end do extrapolate_fmat
        end if
        deallocate (diis_a,diis_c)
      end if
      !
      call TimerStop('DIIS extrapolation')
!   end subroutine diis_extrapolate_rk
