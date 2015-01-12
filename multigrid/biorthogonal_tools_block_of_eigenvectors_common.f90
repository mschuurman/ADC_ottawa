!   subroutine biorthogonalize_block_of_eigenvectors_rk(mo_block)
!     complex(rk), intent(inout) :: mo_block(:,:,:)
      !
      !  Maintain type-generic code below this point
      !
      integer(ik)                          :: imo, jmo, imax
      integer(ik)                          :: nmos, alloc
      real(kind(mo_block))                 :: nrm
      complex(kind(mo_block))              :: ovr, scl
      complex(kind(mo_block)), allocatable :: sblock(:,:), cleft(:,:)
      !
      nmos = size(mo_block,dim=2)
      ! write (out,*) 'Doing block of ',nmos
      !
      !  Step 1: Gram-Schmidt orthogonalize right eigenvectors
      !
      right_vector: do imo=1,nmos
        !
        !  Orthogonalize
        !
        right_vector_ortho: do jmo=1,imo-1
          ovr = sum(conjg(mo_block(:,jmo,2))*mo_block(:,imo,2))
          ! write (out,*) ' overlap ',jmo,imo,' = ',ovr
          mo_block(:,imo,2) = mo_block(:,imo,2) - ovr * mo_block(:,jmo,2)
        end do right_vector_ortho
        !
        !  Normalize to Euclidian norm 1, making the largest component positive, real
        !
        nrm  = real(sum(conjg(mo_block(:,imo,2))*mo_block(:,imo,2)),kind=kind(mo_block))
        imax = maxloc(abs(mo_block(:,imo,2)),dim=1)
        ! write (out,*) ' initial length of ',imo,' was ',nrm
        ! write (out,*) ' largest element was ',mo_block(imax,imo,2),' at ',imax
        scl  = abs(mo_block(imax,imo,2)) / (mo_block(imax,imo,2) * sqrt(nrm))
        mo_block(:,imo,2) = scl*mo_block(:,imo,2)
      end do right_vector
      !
      !  Step 2: set up linear problem for left eigenvector rotation
      !
      allocate (sblock(nmos,nmos),cleft(nmos,nmos),stat=alloc)
      if (alloc/=0) then
        call stop('biorthogonal_tools%biorthogonalize_block_of_eigenvectors - out of memory')
      end if
      !
      cleft = 0
      overlap_right: do jmo=1,nmos
        overlap_left: do imo=1,nmos
          sblock(imo,jmo) = sum(mo_block(:,imo,1)*mo_block(:,jmo,2))
          ! write (out,*) 'left-right overlap of ', imo, jmo, sblock(imo,jmo)
        end do overlap_left
        cleft(jmo,jmo) = 1
      end do overlap_right
      sblock = transpose(sblock)
      call lapack_gesv(sblock,cleft)
      ! write (out,*) ' solutions for left norm are: '
      ! write (out,*) cleft
      !
      !  Rotate left eigenvectors
      !
      mo_block(:,:,1) = mt_matmul(mo_block(:,:,1),cleft)
      !
      deallocate (sblock,cleft)
!   end subroutine biorthogonalize_block_of_eigenvectors_rk
