!   subroutine bt_eigenvectors_rk(nev,eval,evec,epsx)
!     integer(ik), intent(in)     :: nev        ! Number of eigenvalues
!     complex(rk), intent(inout) :: eval(:)     ! In:  Eigenvalues, in no particular order
!                                               ! Out: Eigenvalues, in order of increasing real part
!     complex(rk), intent(inout) :: evec(:,:,:) ! In/out: Left [(:,:,1)] and right [(:,:,2)] eigenvectors
!                                               ! On output, the eigenvectors are biorthogonal.
!     real(rk), intent(in)       :: epsx        ! The threshold for detecting degenerate eigenvalues
!                                               ! Zero or negative uses tighest sensible value given the accuracy
      !
      !  Maintain type-generic code below this point
      !
      integer(ik)      :: order(nev)
      integer(ik)      :: iev1, iev2, iev3
      real(kind(epsx)) :: eps
      !
      !  Lapack's geev returns solutions we are not quite happy with:
      !  the order of the eigenvalues is not to our liking; the left
      !  eigenvectors have been conjugated, and eigenvectors are 
      !  normalized individually, rather than as biorthogonal pairs.
      !
      !  Before we do anything else, we'll have to fix the order
      !  and conjugation.
      !
      call TimerStart('Bi-orthogonalize eigenvectors')
      call order_keys(real(eval,kind=kind(eval)),order)
      eval = eval(order)
      evec(:,:,1) = conjg(evec(:,order,1))
      evec(:,:,2) =       evec(:,order,2)
      !
      !  For some mysterious reasons, LAPACK's ZGEEV does not guarantee
      !  bi-orthogonal eigenvectors if degenerate eigenvalues are present
      !  (they nearly always are for us!). Fixing this up is a bit more
      !  work than I would have hoped for ...
      !
      eps  = eigenvalue_degeneracy_epsilon(epsx,eval,'biorthogonalize_eigenvectors')
      iev1 = 1
      scan_for_degeneracies: do while (iev1<=nev)
        iev3 = iev1
        iev2 = iev3 ! This assignment is not needed, since the loop below
                    ! will always execute at least once; however, it helps
                    ! to shut off a spurious warning from gfortran
        find_identical_eigenvalues: do while(iev3<=nev)
          if (abs(eval(iev3)-eval(iev1))>eps) exit find_identical_eigenvalues
          iev2 = iev3
          iev3 = iev3 + 1
        end do find_identical_eigenvalues
        call biorthogonalize_block_of_eigenvectors(evec(:,iev1:iev2,:))
        eval(iev1:iev2) = sum(eval(iev1:iev2))/(iev2-iev1+1)
        !
        iev1 = iev2 + 1
      end do scan_for_degeneracies
      !
      call bt_sanity_check(evec,'biorthogonalize_eigenvectors')
      call TimerStop('Bi-orthogonalize eigenvectors')
!   end subroutine bt_eigenvectors_rk
