!   subroutine gamess_2e_integrals_real(what,v,batch,a,b,c,d,op_param,accuracy)
!     character(len=*), intent(in)                      :: what        ! What to evaluate
!     real(rk), intent(out)                             :: v(:,:,:,:)  ! Buffer for the results
!                                                                      ! Integrals are in the charge-cloud order: that is,
!                                                                      ! first two indices correspond to the first variable;
!                                                                      ! the last two indices correspond to the second variable
!     integer(ik), intent(in)                           :: batch(:)    ! Batch indices
!     type(gam_structure), intent(in), target, optional :: a, b, c, d  ! Contexts to work on. In principle, we can handle
!                                                                      ! the situation with every index being on a separate
!                                                                      ! molecule - not that this is hard for an AO integral
!     real(rk), intent(in), optional                    :: op_param(:) ! Additional parameters for the operator
!     real(rk), intent(in), optional                    :: accuracy    ! Desired accuracy of the integrals; contributions
!                                                                      ! smaller than this can be neglected. 
!                                                                      ! The default is to be as accurate as possible.
!                                                                      ! Negative accuracy will restore the default behaviour
!     !
      type(gam_structure), pointer :: g_a, g_b, g_c, g_d            ! All missing contexts will be replaced by gam_def
      type(gam_operator_data)      :: op_data                       ! Operator data
      real(kind(v))                :: xyz(3,4)                      ! Coordinates of the basis functions within each batch
      integer(ik)                  :: sh_l (4)                      ! Angular momentum of each batch; all shells within the batch
                                                                    ! are supposed to have the same angular momentum.
      integer(ik)                  :: sh_ns(4)                      ! Number of shells in each batch
      integer(ik), allocatable     :: sh_np(:,:)                    ! Number of primitives in each shell [2nd index]
      real(kind(v)), allocatable   :: p_zet(:,:,:)                  ! Primitive exponents in each shell [2nd index] and centre [3rd index]
      real(kind(v)), allocatable   :: p_c  (:,:,:)                  ! Primitive contraction coefficients in each shell
      real(kind(v)), allocatable   :: p_c_max(:,:,:)                ! Largest absolute contraction coefficient over all instanced of a primitive
      logical, allocatable         :: p_master(:,:,:)               ! "Master" copy of the exponent; other contractions using the
                                                                    ! same exponent will be listed in p_z_same
      integer(ik), allocatable     :: p_z_same(:,:,:,:)             ! Next contraction with an identical exponent
                                                                    ! First index: 1 = shell index; 2 = contraction within shell
                                                                    ! 2nd, 3rd, and 4th indices are as the 1st etc in p_zet/p_c/p_master
      integer(ik)                  :: max_shells                    ! Largest possible batch size across all inputs
      integer(ik)                  :: max_contractions              ! Longest possible contraction across all inputs
      integer(ik)                  :: alloc
      character(len=40)            :: lwhat                         ! Local copy of "what"
      real(kind(v))                :: cutoff_2e                     ! Cutoff
      !
      call TimerStart('GAMESS 2e '//trim(what))
      !
      !  Sanity check on the batch() argument
      !
      if (size(batch)/=4) stop 'import_gamess%gamess_2e_integrals - bad batch argument'
      !
      !  Fill out structure pointers
      !
      g_a => gam_def ; if (present(a)) g_a => a 
      g_b => gam_def ; if (present(b)) g_b => b 
      g_c => gam_def ; if (present(c)) g_c => c 
      g_d => gam_def ; if (present(d)) g_d => d 
      !
      !  Allocate local buffers for the exponents and contraction coefficients within the batch
      !
      max_shells = max(g_a%max_shells,g_b%max_shells,g_c%max_shells,g_d%max_shells)
      max_contractions = max(g_a%max_contractions,g_b%max_contractions,g_c%max_contractions,g_d%max_contractions)
      !
      allocate (sh_np(max_shells,4), &
                p_zet     (max_contractions,max_shells,4), &
                p_master  (max_contractions,max_shells,4), &
                p_c       (max_contractions,max_shells,4), &
                p_c_max   (max_contractions,max_shells,4), &
                p_z_same(2,max_contractions,max_shells,4),stat=alloc)
      if (alloc/=0) then
        write (out,"('import_gamess%gamess_2e_integrals: Error ',i8,' allocating basis tables')") alloc
        stop 'import_gamess%gamess_2e_integrals - allocation'
      end if
      !
      !  Check batch indices and corresponding buffer sizes, and fill per-batch data.
      !
      call fetch_batch(1_ik,g_a)
      call fetch_batch(2_ik,g_b)
      call fetch_batch(3_ik,g_c)
      call fetch_batch(4_ik,g_d)
      !
      lwhat = what
      op_data%op_name = lwhat(4:)
      !
      !
      select case (what)
        case default
          write (out,"('import_gamess%gamess_2e_integrals: ',a,' are not implemented')") trim(what)
          stop 'import_gamess%gamess_2e_integrals - integral not implemented'
        case ('AO 4C DELTA','AO 4C ONE','AO 4C 1/R','AO 4C R','AO 4C R**2')
          !
          !  All parameter-free integrals in os_basic_integral, including:
          !
          !  'AO 4C DELTA'     - Delta-function of r12 (I know it's silly)
          !                    - It actually isn't - we can use this primitive to compute 
          !                      4-centre 1-electron overlap.
          !  'AO 4C ONE'       - Very complicated implementation of the product of two
          !                      overlap integrals
          !  'AO 4C 1/R'       - The "usual" 2-electron interaction
          !  'AO 4C R'         - Linear interaction
          !  'AO 4C R**2'      - Harmonic interaction
          !
        case ('AO 4C GAUSS','AO 4C GAUSS ONE','AO 4C GAUSS R')
          !
          !  All integrals taking one parameter in the interaction potential, including:
          !
          !  'AO 4C GAUSS'     - Gaussian interaction
          !  'AO 4C GAUSS ONE' - Ditto, but a more complicated implementation
          !  'AO 4C GAUSS R'   - Gaussian-damped linear interaction
          !
          if (.not.present(op_param)) then
            write (out,"('import_gamess%gamess_2e_integrals: exponent missing for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_2e_integrals: Missing operator parameter'
          end if
          if (size(op_param)/=1) then
            write (out,"('import_gamess%gamess_2e_integrals: wrong number of parameters for ',a)") trim(lwhat)
            stop 'import_gamess%gamess_2e_integrals: Wrong parameter count'
          end if
          op_data%omega = real(op_param(1),kind=kind(op_data%omega))
      end select
      !
      !  Done with preparations; calculate the integrals!
      !
      cutoff_2e = -1
      if (present(accuracy)) cutoff_2e = accuracy
      call os_2e_batch(op_data,v,xyz,sh_l,sh_ns,sh_np,p_zet,p_c,p_c_max,p_master,p_z_same,cutoff_2e)
      deallocate (sh_np,p_zet,p_master,p_c,p_c_max,p_z_same)
      !
      call TimerStop('GAMESS 2e '//trim(what))
      !
      contains
      subroutine fetch_batch(ind,gam)
        integer(ik), intent(in)         :: ind  ! Index in the integral
        type(gam_structure), intent(in) :: gam  ! Corresponding context
        !
        integer(ik)   :: ib           ! Batch number
        integer(ik)   :: atom         ! Atom index
        integer(ik)   :: bsize        ! Number of orbitals within this batch
        integer(ik)   :: shell        ! Shell index within atom's basis
        integer(ik)   :: is, is2, iso ! Sequential number of a shell within the batch
        integer(ik)   :: ip, ip2, ipo ! Sequential number of a primitive within a shell
        integer(ik)   :: sh_p1, sh_pn ! First and last contractions within the shell
        integer(ik)   :: nc           ! Number of contractions
        real(kind(v)) :: zref         ! Exponent
        !
        ib = batch(ind)
        if (ib<1 .or. ib>gam%nbatches) then
          write (out,"('import_gamess%gamess_2e_integrals: ',a,' given bad batch ',i10,' for for index ',i2)") &
                 trim(what), ib, ind
          stop 'import_gamess%gamess_2e_integrals: Bad batch index'
        end if
        atom  = gam%batches(1,ib)
        bsize = gam%batches(4,ib)
        if ( size(v,dim=ind)<bsize ) then
          write (out,"('import_gamess%gamess_2e_integrals: ',a,' given buffer of size ',i5,' for batch '" &
                   //",i10,' index ',i2,', but ',i10,' is needed.')") &
                 trim(what), size(v,dim=ind), ib, ind, bsize
          stop 'import_gamess%gamess_2e_integrals: Integral buffer too small'
        end if
        !
        !  Remember all the required information about this batch
        !
        xyz(:,ind) = real(gam%atoms(atom)%xyz,kind=kind(xyz)) / abohr ! We keep atom coordinates in Angstroms; this keeps biting me..
        sh_ns(ind) = gam%batches(2,ib)
        copy_shells: do is=1,sh_ns(ind)
          shell              = gam%batches(4+is,ib)
          sh_l (ind)         = gam%atoms(atom)%sh_l(shell)
          sh_p1              = gam%atoms(atom)%sh_p(shell)
          sh_pn              = gam%atoms(atom)%sh_p(shell+1)-1
          nc                 = sh_pn - sh_p1 + 1
          sh_np(is,ind)      = nc
          p_zet(1:nc,is,ind) = real(gam%atoms(atom)%p_zet(sh_p1:sh_pn),kind=kind(p_zet))
          p_c  (1:nc,is,ind) = real(gam%atoms(atom)%p_c  (sh_p1:sh_pn),kind=kind(p_c))
        end do copy_shells
        !
        !  Flag identical exponents for the primitives, and set up equivalence lists
        !
        p_master(  :,:,ind) = .true. ! Assume what all primitives are unique unless known otherwise
        p_z_same(:,:,:,ind) = 0
        p_c_max (  :,:,ind) = 0
        coalesce_shells: do is=1,sh_ns(ind)
          coalesce_primitives: do ip=1,sh_np(is,ind)
            if (.not.p_master(ip,is,ind)) cycle coalesce_primitives
            zref = p_zet(ip,is,ind)
            !
            !  We can assume what basis was constructed by a sane person, and
            !  a primitive appears only once within a shell. Hence, we can
            !  start scanning with the next shell
            !
            iso = is ; ipo = ip
            p_c_max(ip,is,ind) = abs(p_c(ip,is,ind))
            scan_shells: do is2=is+1,sh_ns(ind)
              scan_primitives: do ip2=1,sh_np(is2,ind)
                if (abs(p_zet(ip2,is2,ind)-zref)>10*spacing(zref)) cycle scan_primitives
                !
                !  This primitive was seen before; flag it, and add location of this
                !  primitive to the chain
                !
                p_master(ip2,is2,ind) = .false.
                p_z_same(:,ipo,iso,ind) = (/ is2, ip2 /)
                iso = is2 ; ipo = ip2
                p_c_max(ip,is,ind) = max(p_c_max(ip,is,ind),abs(p_c(ip2,is2,ind)))
              end do scan_primitives
            end do scan_shells
            !
          end do coalesce_primitives
        end do coalesce_shells
        !
        !  A bit of printing
        !
        if (verbose>=2) then
          write (out,"(/'Shell analysis, index = ',i1)") ind
          write (out,"(1x,a8,1x,a8,1x,a12,1x,a12,1x,a8,1x,a8,1x,a8,1x,a12)") &
                 ' Shell ', ' Primitive ', ' Exponent ', ' Contraction ', ' Master ', 'Same shl', 'Same prm', 'Max. contr.'
          print_shells: do is=1,sh_ns(ind)
            print_primitives: do ip=1,sh_np(is,ind)
              write (out,"(1x,i8,1x,i8,1x,f12.5,1x,f12.5,1x,l8,1x,i8,1x,i8,1x,f12.5)") &
                     is, ip, p_zet(ip,is,ind), p_c(ip,is,ind), p_master(ip,is,ind), p_z_same(:,ip,is,ind), p_c_max(ip,is,ind)
            end do print_primitives
            write (out,"()")
          end do print_shells
        end if
      end subroutine fetch_batch
!   end subroutine gamess_2e_integrals_real
