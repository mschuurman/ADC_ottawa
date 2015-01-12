module opendx
    !
    !  Generate OpenDX-formatted data file for plotting angular distributions.
    !  This routine has been cannibalized from spherical_interpolate.f90, and
    !  should be compatible with the all OpenDX scripts using that output.
    !
  use accuracy
  use timer
  implicit none
  private
  public odx_write_spherical, odx_write_cube
  !
  interface odx_write_spherical
    module procedure odx_write_sherical_real
  end interface odx_write_spherical
  !
  interface odx_write_cube
    module procedure odx_write_cube_real
  end interface odx_write_cube
  !
  integer(ik), parameter :: iu_output = 30  ! I/O unit to use to make OpenDX files.
  !
  contains
  !
  !  Generate OpenDX file for a uniform spherical-product grid
  !
  subroutine odx_write_sherical_real(file,data,comment)
    character(len=*), intent(in)           :: file      ! Name of the output file to write. Blank means
                                                        ! to use standard output.
    real(rk), intent(in)                   :: data(:,:) ! Data to write; theta is the first index
    character(len=*), intent(in), optional :: comment   ! A text string to add to the output
    !
    integer(ik) :: output
    integer(ik) :: itheta, iphi, icol
    integer(ik) :: ntheta, nphi
    real(rk)    :: step_theta, step_phi
    ! 
    call TimerStart('OpenDX Write Spherical Real')
    !
    ntheta = size(data,dim=1)
    nphi   = size(data,dim=2)
    !
    if (file/=' ') then
      output = iu_output
      open (output,form='formatted',status='replace',file=trim(file))
    else
      write (out,"('====== OpenDX file begins below this line ====')")
      output = out
    end if
    !
    step_phi   = twopi/nphi
    step_theta =    pi/ntheta
    !
    !  Write the OpenDX header. To make our life easier, we'll add an extra point for phi
    !  so that the grids fold on themselves after warping. This seems to be easier than 
    !  specifying connectivity explicitly
    !
    write (output,"('object 1 class array type float rank 0 items ',i0,' data 0')") ntheta * (nphi+1)
    write (output,"('attribute ""dep"" string ""positions""')")
    write (output,"('#')")
    write (output,"('object 2 class gridpositions counts ',i0,1x,i0)") ntheta, nphi+1
    write (output,"(' origin   ',2(1x,f14.8))") 0.5_rk * step_theta, 0.5_rk * step_phi
    write (output,"(' delta    ',2(1x,f14.8))") step_theta, 0._rk
    write (output,"(' delta    ',2(1x,f14.8))") 0._rk, step_phi
    write (output,"('attribute ""dep"" string ""positions""')")
    write (output,"('#')")
    write (output,"('object 3 class gridconnections counts ',i0,1x,i0)") ntheta, nphi+1
    write (output,"('attribute ""element type"" string ""quads""')")
    write (output,"('attribute ""dep"" string ""connections""')")
    write (output,"('attribute ""ref"" string ""positions""')")
    write (output,"('#')")
    if (present(comment)) then
      write (output,"('object 4 class string ""',a,'""')") trim(comment)
    end if
    write (output,"('object ""field0"" class field')")
    write (output,"('component ""data"" value 1')")
    write (output,"('component ""positions"" value 2')")
    write (output,"('component ""connections"" value 3')")
    if (present(comment)) then
      write (output,"('component ""comment"" value 4')")
    end if
    write (output,"('attribute ""name"" string ""field0""')")
    write (output,"('#')")
    write (output,"('end')")
    !
    !  Write the data. Data is expected in the "C" order (row-major)
    !
    icol = 0
    loop_theta: do itheta=1,ntheta
      loop_phi: do iphi=1,nphi
        call write_output(itheta,iphi)
      end do loop_phi
      ! Extra point at the end, duplicating iphi=1
      call write_output(itheta,1_ik)
    end do loop_theta
    !
    if (icol/=0) write (output,"()")
    !
    if (file/=' ') then
      close (output)
    else
      write (out,"('====== OpenDX file ends above this line ====')")
    end if
    !
    call TimerStop('OpenDX Write Spherical Real')
    !
    contains
      subroutine write_output(it,ip)
        integer(ik), intent(in) :: it, ip
        !
        write (output,"(g12.6,1x)",advance='no') data(it,ip)
        icol = icol + 1
        if (icol==6) then
          icol = 0
          write (output,"()")
        end if
     end subroutine write_output
  end subroutine odx_write_sherical_real
  !
  !  Generate OpenDX file for a uniformly spaced cube
  !
  subroutine odx_write_cube_real(file,data,origin,step,counts)
    character(len=*), intent(in) :: file      ! Name of the output file to write. Blank means
                                              ! to use standard output.
    real(rk), intent(in)         :: data(:,:) ! Data to write; Fortran index order. The first
                                              ! index is the field number; the second index is
                                              ! the grid point
    real(rk), intent(in)         :: origin(:) ! Grid origin
    real(rk), intent(in)         :: step(:,:) ! Ditection vectors for the grid axes
    integer(ik), intent(in)      :: counts(:) ! Number of grid points along each direction
    !
    integer(ik)       :: output     ! Unit for output
    integer(ik)       :: shape      ! Length of each data vector
    integer(ik)       :: items      ! Total number of data vectors; must match the grid size
    integer(ik)       :: spacedim   ! Space dimenstionality
    integer(ik)       :: griddim    ! Grid dimenstionality
    integer(ik)       :: idim, ipos
    integer(ik)       :: icol
    character(len=20) :: conn_type  ! Grid connection type
    ! 
    call TimerStart('OpenDX Write Cube Real')
    !
    !  Make sure all parameters we got are consistent
    !
    shape    = size(data,dim=1)
    items    = size(data,dim=2)
    spacedim = size(origin)
    griddim  = size(counts)
    !
    select case (griddim)
      case (2) ; conn_type = "quads"
      case (3) ; conn_type = "cubes"
      case default ; conn_type = "invalid"
    end select
    !
    if (product(counts)/=items) then
      write (out,"('Number of data items (',i0,') does not match grid size: ',5(1x,i0))") items, counts
      stop 'opendx%odx_write_cube_real - grid and data incompatible'
    end if
    if (size(step,dim=1)/=spacedim .or. size(step,dim=2)/=griddim) then
      stop 'opendx%odx_write_cube_real - inconsistent grid definition'
    end if
    !
    if (file/=' ') then
      output = iu_output
      open (output,form='formatted',status='replace',file=trim(file))
    else
      write (out,"('====== OpenDX file begins below this line ====')")
      output = out
    end if
    !
    !  Write the OpenDX header. 
    !  OpenDX native files use C array order, so that we need to reverse grid counts
    !
    write (output,"('object 1 class array type float rank 1 shape ',i0,' items ',i0,' data 0')") shape, items
    write (output,"('attribute ""dep"" string ""positions""')")
    write (output,"('#')")
    write (output,"('object 2 class gridpositions counts ',5(1x,i0))") counts(griddim:1:-1)
    write (output,"(' origin   ',5(1x,f14.8))") origin
    write_delta: do idim=griddim,1,-1
    write (output,"(' delta    ',5(1x,f14.8))") step(:,idim)
    end do write_delta
    write (output,"('attribute ""dep"" string ""positions""')")
    write (output,"('#')")
    write (output,"('object 3 class gridconnections counts ',5(1x,i0))") counts(griddim:1:-1)
    write (output,"('attribute ""element type"" string ""',a,'""')") trim(conn_type)
    write (output,"('attribute ""dep"" string ""connections""')")
    write (output,"('attribute ""ref"" string ""positions""')")
    write (output,"('#')")
    write (output,"('object ""field0"" class field')")
    write (output,"('component ""data"" value 1')")
    write (output,"('component ""positions"" value 2')")
    write (output,"('component ""connections"" value 3')")
    write (output,"('attribute ""name"" string ""field0""')")
    write (output,"('#')")
    write (output,"('end')")
    !
    !  Write the data.
    !
    icol = 0
    loop_samples: do ipos=1,items
      loop_elements: do idim=1,shape
        call write_output(data(idim,ipos))
      end do loop_elements
    end do loop_samples
    !
    if (icol/=0) write (output,"()")
    !
    if (file/=' ') then
      close (output)
    else
      write (out,"('====== OpenDX file ends above this line ====')")
    end if
    !
    call TimerStop('OpenDX Write Cube Real')
    !
    contains
      subroutine write_output(val)
        real(rk), intent(in) :: val
        !
        write (output,"(g12.6,1x)",advance='no') val
        icol = icol + 1
        if (icol==6) then
          icol = 0
          write (output,"()")
        end if
     end subroutine write_output
  end subroutine odx_write_cube_real
end module opendx
