!
!  Revision history: 2013 Feb 04 - Separated out routines specific to the real-space 
!                                  grid implementation
!
!
  module ecp_gamess
    use accuracy
    use timer
    use gamess_internal
    use ecp_convert_gamess
    use import_gamess
    use multigrid

    implicit none

    private
    public ecp_molecule, ecp_atom              ! Re-export of ecp_convert_gamess values
    public ecp_convert, ecp_destroy            ! Re-export of ecp_convert_gamess values
    public ecp_apply
    public ecp_build_cartesian_projector       ! Debugging entry point, use with extreme caution
    !
    !  Some global constants ...
    !
    integer(ik), parameter :: verbose     = 1
    !
    contains
    !
    !  ecp_apply assumes that ECP projectors do not overlap; if they do, you are in trouble.
    !  There is no checking; it is your responsibility to make sure ECPs make sense.
    !
    !  Actually, it is OK for the projectors to overlap: the src field is never modified,
    !  so that the ket part of the projectors is applied consistently. The dst field, on
    !  the other hand, is never sensed - so that the bra part of each projector is also
    !  fully consistent.
    !
    subroutine ecp_apply(src,dst,ecp,rot,save_memory)
      integer(ik), intent(in)           :: src         ! Wavefunction to apply projector to
      integer(ik), intent(in)           :: dst         ! Potential to update
      type(ecp_molecule), intent(inout) :: ecp         ! The ECP. It is not modified, but calling convetions
                                                       ! downstream require inout intent.
      real(rk), intent(in), optional    :: rot(:,:)    ! Rotation matrix; if present, will replace rotation matrix
                                                       ! within the ecp. Otherwise, the rotation matrix stored in
                                                       ! the ecp will be used
      logical, intent(in), optional     :: save_memory ! If true, Cartesian representations of the prejectors
                                                       ! will be created as needed, and destroyed immediately
                                                       ! afterwards.
      !
      integer(ik)              :: iprj, alloc
      complex(rk), allocatable :: wgt(:)
      logical                  :: be_slow
      !
      if (.not.ecp%active) return
      call TimerStart('Evaluate ECP')
      !
      be_slow = .false.
      if (present(save_memory)) be_slow = .not.save_memory
      !
      !  Since the projectors are in the form:
      !
      !    |dst> += |prj><prj|src>
      !
      !  the order in which they are applied does not make any difference.
      !
      apply_prjs: do iprj=1,ecp%necps
        if (present(rot)) then
          ecp%ecps(iprj)%projectors%rotmat = rot
        end if
        !
        !  Build Cartesian representation of the ECP if it is not present.
        !
        if (.not.allocated(ecp%ecps(iprj)%grid)) then
          call ecp_build_cartesian_projector(ecp%ecps(iprj))
        end if
        allocate (wgt(size(ecp%ecps(iprj)%vshift)),stat=alloc)
        if (alloc/=0) stop 'ecp_gamess%ecp_apply - out of memory'
        !
        !  Phase 1: Calculate overlaps between the ECP projectors and the wavefunction
        !
        call FieldECPProject(src,ecp%ecps(iprj)%projectors,ecp%ecps(iprj)%grid,ecp%ecps(iprj)%rmax,wgt)
        !
        !  Phase 2: Apply energy shifts
        !
        ! write (out,"(' ECP projector overlaps: '/(10g16.7))") wgt
        wgt = wgt * ecp%ecps(iprj)%vshift
        ! write (out,"(' ECP projector weights: '/(10g16.7))") wgt
        !
        !  Phase 3: add shifted projectros back to the wavefunction
        !
        call FieldECPApply(dst,ecp%ecps(iprj)%projectors,ecp%ecps(iprj)%grid,ecp%ecps(iprj)%rmax,wgt)
        deallocate (wgt)
        !
        !  If we are asked to conserve memory, release the ECP. This will make
        !  things _much_ slower! 
        !
        if (be_slow) deallocate(ecp%ecps(iprj)%grid)
      end do apply_prjs
      call TimerStop('Evaluate ECP')
    end subroutine ecp_apply
    !
    !  Transform projectors specified as a basis set expansion to projectors on grid
    !  This is a semi-debugging routine, do not call it!
    !
    subroutine ecp_build_cartesian_projector(ecp_at)
      type(ecp_atom), target, intent(inout) :: ecp_at ! Given basis-set representation of the ECP projectors,
                                                      ! construct the equivalent Cartesian form.
      !
      type(gam_structure), pointer :: prj               ! Projectors, represented by basis set expansion
      integer(ik)                  :: nprj              ! Number of projectors
      integer(ik)                  :: np   (3)          ! Extent of the full grid
      integer(ik)                  :: npp(2,3)          ! Extent of the projector sub-grid
      real(rk)                     :: dx   (3)          ! Grid spacing
      real(rk), allocatable        :: x(:), y(:), z(:)  ! Grid coordinates of the full grid
      integer(ik), allocatable     :: i123(:)           ! 1, 2, 3, ... - see gamess_load_orbitals if you care
      complex(rk), allocatable     :: p123(:,:,:,:)     ! Values of the projector; ditto
      real(rk), allocatable        :: pxyz(:,:,:,:)     ! Explicit XYZ coordinates for gamess_load_orbitals
      integer(ik)                  :: alloc
      integer(ik)                  :: ix, iy, iz, ip
      real(rk), allocatable        :: norm(:)
      !
      call TimerStart('Build Cartesian ECP projectors')
      !
      !  Get projector parameters
      !
      prj => ecp_at%projectors
      nprj = prj%nvectors
      call FieldECPGetExtent(prj,ecp_at%rmax,npp)
      !
      !  Allocate the serious stuff - the projector buffer
      !
      if (verbose>=0) then
        write (out,"('ECP radial extent for atom ',a,' is ',f10.3,' Bohr.')") &
               trim(prj%atoms(1)%name), ecp_at%rmax
        write (out,"('Cartesian ECP projectors for atom ',a,' require ',f10.3,' MBytes of RAM.')") &
               trim(prj%atoms(1)%name), product(npp(2,:)-npp(1,:)+0._rk)*nprj*rk_bytes/1024._rk**2
      end if
      allocate (ecp_at%grid(npp(1,1):npp(2,1),npp(1,2):npp(2,2),npp(1,3):npp(2,3),nprj),stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' nprj = ',i6,' npp = ',6i6)") alloc, nprj, npp
        stop 'ecp_gamess%build_cartesian_projector - allocate failure (1)'
      end if
      !
      !  Query grid parameters; we'll need grid spacings and coordinates
      !
      np = FieldGridNPoints(1_ik)
      dx = FieldGridSpacing(1_ik)
      !
      !  Allocate small stuff we'll be using as temp
      !
      allocate (x(np(1)),y(np(2)),z(np(3)),i123(nprj),norm(nprj), &
                p123(npp(1,1):npp(2,1),npp(1,2):npp(2,2),1,nprj), &
                pxyz(3,npp(1,1):npp(2,1),npp(1,2):npp(2,2),1), stat=alloc)
      if (alloc/=0) then
        write (out,"('Error ',i8,' in allocate. nprj = ',i6,' np = ',3i6,' npp = ',6i6)") &
               alloc, nprj, np, npp
        stop 'ecp_gamess%build_cartesian_projector - allocate failure (2)'
      end if
      call FieldGridCoordinates(1_ik,x,1_ik)
      call FieldGridCoordinates(2_ik,y,1_ik)
      call FieldGridCoordinates(3_ik,z,1_ik)
      !
      !  Fill orbital indices: The GAMESS import calling conventions are a bit strange ...
      !
      fill_123: do ix=1,size(i123)
        i123(ix) = ix
      end do fill_123
      !
      !  Fill the explicit grid point coordinates
      !
      fill_x: do iy=npp(1,2),npp(2,2)
        pxyz(1,:,iy,1) = x(npp(1,1):npp(2,1))
      end do fill_x
      fill_y: do ix=npp(1,1),npp(2,1)
        pxyz(2,ix,:,1) = y(npp(1,2):npp(2,2))
      end do fill_y
      !
      !  Do not worry about parallel execution: gamess_load_orbitals will already
      !  run in parallel where it makes sense.
      !
      pt_z: do iz=npp(1,3),npp(2,3)
        pxyz(3,:,:,1) = z(iz)
        call gamess_load_orbitals(mos=i123,dst=i123,grid=p123,structure=prj,coord=pxyz,dx=dx)
        !
        !  gamess_load_orbitals gives us complex result, even though all the projectors
        !  are actually real; keep the real part in out buffer, and drop the zeros
        !
        ecp_at%grid(:,:,iz:iz,:) = real(p123,kind=rk)
      end do pt_z
      !
      !  We have the raw projectors; although gamess_load_orbitals() should have taken
      !  care to oversample the volume elements where projectors vary rapidly, it is
      !  still prudent to check the norms of the projectors. Resetting them to one,
      !  on the other hand, is probably not a good idea.
      !
      !$omp parallel do default(none) private(ip) shared(nprj,norm,ecp_at,dx)
      renormalize_projectors: do ip=1,nprj
        norm(ip) = sqrt(product(dx)*sum(ecp_at%grid(:,:,:,ip)**2))
        ! Do not renormalize here.
        ! ecp_at%grid(:,:,:,ip) = ecp_at%grid(:,:,:,ip) * (1._rk/norm(ip))
      end do renormalize_projectors
      !$omp end parallel do
      !
      if (verbose>=1) then
        write (out,"()")
        write (out,"(t5,a5,2x,a12,2x,a12,2x,a12)") ' prj ', ' Shift, H ', ' <R>, Bohr ', '   Norm     ', &
                                                   '-----', '----------', '-----------', '------------'
        print_norms: do ip=1,nprj
          write (out,"(t5,i5,2x,f12.5,2x,f12.5,2x,f12.8)") ip, ecp_at%vshift(ip), ecp_at%rv(ip), norm(ip)
        end do print_norms
        write (out,"()") 
      end if
      !
      deallocate (x,y,z,i123,norm,p123,pxyz)
      call TimerStop('Build Cartesian ECP projectors')
    end subroutine ecp_build_cartesian_projector
    !
  end module ecp_gamess
