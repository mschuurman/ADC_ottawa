module kinetic
!
!  Calculation of kinetic error diffusion stencil, used for
!  preconditioning. The idea is to build a 3x3x3 orbital
!  shape, which will give Laplacian value of 1 at the centre, 
!  while minimizing the incidental Laplacians in the 5x5x5 
!  cube around it.
!
  use lapack
  implicit none
  private
  public kineticBuildStencil
!
  integer(ik), parameter :: verbose = 0
!
  contains

  subroutine kineticBuildStencil(dx,ps)
    real(rk), intent(in)  :: dx(3)              ! Grid step size along each direction
    real(rk), intent(out) :: ps(-1:1,-1:1,-1:1) ! Optimal wavefunction shape
    !
    integer(ik) :: var_cube(-3:3,-3:3,-3:3)     ! Indices of the wavefunction coefficients
                                                ! 0 means function is set to zero.
    integer(ik) :: form_loc(-2:2,-2:2,-2:2)     ! Location of the laplacian formula in
                                                ! the lap_code array
    real(rk)    :: lap_code(8,5**3)             ! Weights of the wavefunction coefficients
                                                ! in the laplacian, at all points of 5x5x5 
                                                ! cube.
    real(rk)    :: sx(3)                        ! Coordinate weights in laplacian
    real(rk)    :: error_mat(8,8)               ! Our error function is a general quadratic 
                                                ! form
    real(rk)    :: best_c(8)
    !
    sx = 1.0_rk / dx**2
    if (verbose>=2) write(out,"(/'Directional weights: ',3f12.5)") sx
    !
    ! For non-uniform grid spacings, writing out error diffusion
    ! equations by hand is a little too much work. Instead, we'll
    ! do some symbolic algebra, to build the equations on the fly.
    !
    call init_var_cube
    !
    ! Evaluate the Laplacian for the 5x5x5 surrounding cube
    !
    call laplacian
    !
    ! Calculate the error function, for the "diffused" part of the
    ! laplacian. This is given by a quadratic form in the eight
    ! independent wavefunction coefficients.
    !
    call error_function
    !
    ! Determine the best possible error diffusion mode
    !
    call best_mode
    !
    ! Expand coefficients to form the stencil
    !
    call expand_mode
    !
    if (verbose>=2) then
      write (out,"('Laplacian for unit wavefunction norm: ',f12.5)") &
             1.0_rk/sqrt(product(dx)*sum(ps**2))
    end if
    !
    ! That's that.
    !
    return
    !
    contains

      subroutine init_var_cube
        !
        !  Outer part of the cube must be zero
        !
        var_cube = 0
        !
        !  Our wavefunction must be of D2d symmetry around the 
        !  central point, to prevent artificial symmetry breaking
        !  Let's assign unique variables accordingly:
        !
        var_cube(     0,     0,     0) = 1   ! Centre
        var_cube(-1:1:2,     0,     0) = 2   ! +/-dx
        var_cube(     0,-1:1:2,     0) = 3   ! +/-dy
        var_cube(     0,     0,-1:1:2) = 4   ! +/-dz
        var_cube(-1:1:2,-1:1:2,     0) = 5   ! +/-dx, +/-dy
        var_cube(-1:1:2,     0,-1:1:2) = 6   ! +/-dx, +/-dz
        var_cube(     0,-1:1:2,-1:1:2) = 7   ! +/-dy, +/-dz
        var_cube(-1:1:2,-1:1:2,-1:1:2) = 8   ! +/-dx, +/-dx, +/-dz
      end subroutine init_var_cube

      subroutine laplacian
        integer(ik) :: ix, iy, iz, iform, ic, ip
        integer(ik) :: fc, fd(2,3)
        !
        if (verbose>=2) then
          write (out,"(/t15,'Laplacian formulas')")
          write (out,"(3a3,2x,8a10)") 'DX', 'DY', 'DZ', 'C', 'CX', 'CY', 'CZ', &
                 'CXY', 'CXZ', 'CYZ', 'CXYZ'
        end if
        !
        iform    = 0
        z: do iz=-2,2
          y: do iy=-2,2
            x: do ix=-2,2
              iform = iform + 1
              form_loc(ix,iy,iz) = iform
              !
              !  Get variables
              !
              fc      = var_cube(ix  ,iy  ,iz  )  ! Central coefficient
              fd(1,1) = var_cube(ix-1,iy  ,iz  )  ! Displacements along x
              fd(2,1) = var_cube(ix+1,iy  ,iz  )
              fd(1,2) = var_cube(ix  ,iy-1,iz  )  ! Displacements along y
              fd(2,2) = var_cube(ix  ,iy+1,iz  )
              fd(1,3) = var_cube(ix  ,iy  ,iz-1)  ! Displacements along z
              fd(2,3) = var_cube(ix  ,iy  ,iz+1)
              !
              !  Update formulas:  lap = (f(x+dx)+f(x-dx)-2*f(x))/dx**2 + ...
              !
              lap_code( :,iform) = 0
              if (fc/=0) lap_code(fc,iform) = -2.0_rk*sum(sx)
              dir: do ic=1,3
                do ip=1,2
                  if (fd(ip,ic)/=0) &
                    lap_code(fd(ip,ic),iform) = lap_code(fd(ip,ic),iform) + sx(ic)
                end do
              end do dir
              if (verbose>=2) then
                write (out,"(3i3,2x,8f10.5)") ix, iy, iz, lap_code(:,iform)
              end if
              !
            end do x
          end do y
        end do z
        if (iform/=size(lap_code,dim=2)) &
          stop 'kinetic%laplacian - logic'
      end subroutine laplacian

      subroutine error_function
        integer(ik) :: iform
        integer(ik) :: i, j
        !
        ! The loop below actually does the equivalent of DSYRK.
        ! Nowever, the cost is so minor that there is no point
        ! in calling the BLAS routine.
        !
        error_mat = 0
        form: do iform=1,size(lap_code,dim=2)
          if (iform==form_loc(0,0,0)) cycle form ! Central point - exclude
          do i=1,size(lap_code,dim=1)
            do j=1,size(lap_code,dim=1)
              error_mat(i,j) = error_mat(i,j) + lap_code(i,iform) * lap_code(j,iform)
            end do 
          end do
        end do form
        !
        if (verbose>=2) then
          write (out,"(/t15,'Error matrix:'/)")
          write (out,"((1x,8f12.2))") error_mat
        end if
        !
      end subroutine error_function
 
      subroutine best_mode
        real(rk)    :: error_eig(8)  ! Eigenvalues of the error function
        integer(ik) :: imode
        integer(ik) :: best          ! Index of the best error diffusion mode
        real(rk)    :: lap_c
        integer(ik) :: form_c        ! Expression for the cenral Laplacian
        integer(ik) :: ix, iy, iz
        !
        ! The best possible error-diffusing function corresponds to
        ! the smallest absolute eigenvalue of the error matrix. A
        ! small complication is that different modes will need 
        ! different amplitudes to give us laplacian of 1 at the 
        ! centre - so before making the choice, we'll have to adjust
        ! the normalization.
        !
        call lapack_syev(error_mat,error_eig)
        !
        if (verbose>=2) then
          write (out,"(/' Error eigenvalues: ',(t20,4e12.5))") error_eig
        end if
        !
        form_c = form_loc(0,0,0)
        renorm: do imode=1,8
          lap_c = 1.0_rk / dot_product(lap_code(:,form_c),error_mat(:,imode))
          error_mat(:,imode) = error_mat(:,imode) * lap_c
          error_eig(  imode) = error_eig(  imode) * lap_c**2
        end do renorm
        !
        if (verbose>=2) then
          write (out,"(/' Renormalized eigenvalues: ',(t20,4e12.5))") error_eig
        end if
        !
        best   = minloc(abs(error_eig),dim=1)
        best_c = error_mat(:,best)
        !
        if (verbose>=2) then
          write (out,"(/' Best error diffusion mode is ',i2,':')") best
          write (out,"((t12,4f14.7))") best_c
        end if
        !
        if (verbose>=2) then
          write (out,"(/t15,'Diffused Laplacian:')") 
          write (out,"(1x,3a3,2x,a14)") 'DX', 'DY', 'DZ', 'Laplacian'
          p_z: do iz=-2,2
            p_y: do iy=-2,2
              p_x: do ix=-2,2
                write (out,"(1x,3i3,2x,f14.7)") ix, iy, iz, &
                       dot_product(lap_code(:,form_loc(ix,iy,iz)),best_c)
              end do p_x
            end do p_y
          end do p_z
        end if
        !
      end subroutine best_mode

      subroutine expand_mode
        integer(ik) :: ix, iy, iz
        !
        forall (ix=-1:1)
          forall (iy=-1:1)
            forall (iz=-1:1)
              ps(ix,iy,iz) = best_c(var_cube(ix,iy,iz))
            end forall
          end forall
        end forall
        !
        if (verbose>=2) then
          write (out,"(/t15,'Optimal wavefunction for unit Laplacian:')") 
          write (out,"(1x,3a3,2x,a14)") 'DX', 'DY', 'DZ', 'Psi'
          p_z: do iz=-1,1
            p_y: do iy=-1,1
              p_x: do ix=-1,1
                write (out,"(1x,3i3,2x,f14.7)") ix, iy, iz, ps(ix,iy,iz)
              end do p_x
            end do p_y
          end do p_z
        end if
      end subroutine expand_mode

  end subroutine kineticBuildStencil

end module kinetic
