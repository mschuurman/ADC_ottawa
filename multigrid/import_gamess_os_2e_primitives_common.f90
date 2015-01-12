!   subroutine os_2e_primitives_real(what,xyz,zet,p2,nrec,vb,w_cut,zero)
!     type(gam_operator_data), intent(in) :: what        ! Operator parameters
!     real(ark), intent(in)               :: xyz(:,:)    ! Coordinates of the basis functions within each batch
!     real(rk), intent(in)                :: zet  (:)    ! Primitive exponents
!     integer(ik), intent(in)             :: p2   (:)    ! Ending position of the desired block of integrals
!     integer(ik), intent(in)             :: nrec        ! Maximum recursion order
!     real(rk), intent(out)               :: vb(0:p2(1),0:p2(2),0:p2(3),0:p2(4),0:nrec)
!     real(rk), intent(in)                :: w_cut       ! Absolute cutoff for the integrals
!     logical, intent(out)                :: zero        ! True if all integrals are smaller than w_cut
!     !
      real(kind(vb)) :: zeta, eta, rho, p(3), q(3), t, sab, scd        ! See Ahlrichs eqs. 3-6
      real(kind(vb)) :: s                  ! Integral being evaluated
      integer(ik)    :: i,   j,   k,   l   ! Integral indices
      integer(ik)    :: id1, jd1, kd1, ld1 ! "parent" indices - drop one power for the cordinate ic
      integer(ik)    :: id2, jd2, kd2, ld2 ! "parent" indices - drop two powers for the coordinate ic
      integer(ik)    :: il,  jl,  kl,  ll  ! Sum of all coordinate powers for given function
      integer(ik)    :: ic                 ! Coordinate to recurse along
      integer(ik)    :: m                  ! Recursion order
      real(kind(vb)) :: grow_factor        ! Estimate of the maximum amplification factor in integral recursion
      real(kind(vb)) :: gf_tmp
      !
      ! There is a lot of overlap between formulae below and os_common_primitives()
      ! 
      zeta = zet(1) + zet(2)                                        ! Ahlrichs eq. 3
      eta  = zet(3) + zet(4)                                        ! ... ditto
      rho  = zeta * eta / (zeta + eta)                              ! ... ditto
      p(:) = (zet(1)*xyz(:,1) + zet(2)*xyz(:,2))/zeta               ! Ahlrichs eq. 4
      q(:) = (zet(3)*xyz(:,3) + zet(4)*xyz(:,4))/eta                ! ... ditto
      t    = rho * sum((p-q)**2)                                    ! Ahlrichs eq. 5
      sab  = exp(-(zet(1)*zet(2)/zeta)*sum((xyz(:,1)-xyz(:,2))**2)) ! Ahlrichs eq. 6
      scd  = exp(-(zet(3)*zet(4)/eta )*sum((xyz(:,3)-xyz(:,4))**2)) ! ... ditto
      !
      ! write (out,"(' zeta = ',f25.15,' eta = ',f25.15,' rho = ',f25.15)") zeta, eta, rho
      ! write (out,"(' p = ',3f25.15)") p
      ! write (out,"(' q = ',3f25.15)") q
      ! write (out,"(' xyz(:,1) = ',3f25.15)") xyz(:,1)
      ! write (out,"(' xyz(:,2) = ',3f25.15)") xyz(:,2)
      ! write (out,"(' xyz(:,3) = ',3f25.15)") xyz(:,3)
      ! write (out,"(' xyz(:,4) = ',3f25.15)") xyz(:,4)
      ! write (out,"(' t = ',f25.15,' sab = ',f25.15,' scd = ',f25.15)") t, sab, scd
      !
      !  Prepare basic integrals for recursion; this takes care of the S-type
      !  integrals.
      !
      call os_basic_integral(what,rho,t,nrec,vb(0,0,0,0,0:nrec))
      vb(0,0,0,0,0:nrec) = vb(0,0,0,0,0:nrec) * sab * scd * (real(pi_xrk,kind(zeta))/(zeta+eta))**real(1.5_xrk,kind(zeta))
      !
      !  Ready to apply the cut-offs. We'll try to be extremely conservative.
      !
      grow_factor = 1
      ! (P-A) and (P-Q) prefactors in the recursion
      m      = sum(ang_nxyz(p2,:))
      gf_tmp = maxval(abs(p-q))
      do i=1,4
        gf_tmp = max(gf_tmp,maxval(abs(xyz(:,1) - p)))
      end do
      grow_factor = max(grow_factor,gf_tmp**m)
      grow_factor = max(grow_factor,max(rho,real(0.5_xrk,kind(grow_factor))*MathFactorial(m,zeta))/zeta)
      !
      ! write (out,*) '       w_cut = ', w_cut
      ! write (out,*) ' grow_factor = ', grow_factor
      ! write (out,*) '  max. basic = ', maxval(abs(vb(0,0,0,0,0:nrec)))
      if (maxval(abs(vb(0,0,0,0,0:nrec)))*grow_factor<=w_cut) then
        ! write (out,*) ' batch is zero'
        zero = .true.
        return
      end if
      !
      !  This integral is not obviously negligible, and will have to be evaluated.
      !
      zero = .false.
      !
      ! write (out,"(' basic integrals: '/(4f25.15))") vb(0,0,0,0,0:nrec)
      !
      !  Stage 1: keep the last three indices at zero; increment the first index
      !           Ahlrichs' equation 32. All other indices are fixed at 1 (S functions)
      !
      boot_i: do i=1,p2(1)
        call get_parent(i,ic,id1)
        il = sum(ang_nxyz(i,:))
        boot_i_order: do m=0,nrec-il
            s = (p(ic)-xyz(ic,1)) * vb(id1,0,0,0,m) 
            s = s - (rho/zeta) * (p(ic)-q(ic)) * vb(id1,0,0,0,m+1)
            id2 = drop_xyz(id1,ic)
            if (id2>=0) then
              s = s + (ang_nxyz(id1,ic)/(2*zeta)) * (vb(id2,0,0,0,m) - (rho/zeta) * vb(id2,0,0,0,m+1))
            end if
            vb(i,0,0,0,m) = s
        end do boot_i_order
      end do boot_i
      !
      !  Stage 2: keep the second and the fourth indices to zero; increment the third index
      !           This is Ahlrichs' equation 33.
      !
      boot_k: do k=1,p2(3)
        call get_parent(k,ic,kd1)
        kl  = sum(ang_nxyz(k,:))
        kd2 = drop_xyz(kd1,ic)
        boot_k_i: do i=0,p2(1)
          id1 = drop_xyz(i,ic)
          il  = sum(ang_nxyz(i,:))
          boot_k_order: do m=0,nrec-(il+kl)
            s = (q(ic)-xyz(ic,3)) * vb(i,0,kd1,0,m)
            s = s + (rho/eta) * (p(ic)-q(ic)) * vb(i,0,kd1,0,m+1)
            if (kd2>=0) then
              s = s + (ang_nxyz(kd1,ic)/(2*eta)) * (vb(i,0,kd2,0,m) - (rho/eta) * vb(i,0,kd2,0,m+1))
            end if
            if (id1>=0) then
              s = s + (ang_nxyz(i,ic)/(2*(zeta+eta))) * vb(id1,0,kd1,0,m+1)
            end if
            vb(i,0,k,0,m) = s
          end do boot_k_order
        end do boot_k_i
      end do boot_k
      !
      !  Stage 3: keep the fourth index at zero; increment the second index. This is a
      !           special case of Ahlrichs' equation 31. Note that eq. 31 as printed
      !           is missing a closing square bracket after the second "b-1" term.
      !
      boot_j: do j=1,p2(2)
        call get_parent(j,ic,jd1)
        jl  = sum(ang_nxyz(j,:))
        jd2 = drop_xyz(jd1,ic)
        boot_j_k: do k=0,p2(3)
          kd1 = drop_xyz(k,ic)
          kl  = sum(ang_nxyz(k,:))
          boot_j_i: do i=0,p2(1)
            id1 = drop_xyz(i,ic)
            il  = sum(ang_nxyz(i,:))
            boot_j_order: do m=0,nrec-(jl+kl+il)
              s = (p(ic)-xyz(ic,2)) * vb(i,jd1,k,0,m)
              s = s - (rho/zeta) * (p(ic)-q(ic)) * vb(i,jd1,k,0,m+1)
              if (jd2>=0) then
                s = s + (ang_nxyz(jd1,ic)/(2*zeta)) * (vb(i,jd2,k,0,m) - (rho/zeta) * vb(i,jd2,k,0,m+1))
              end if
              if (id1>=0) then
                s = s + (ang_nxyz(i,ic)/(2*zeta)) * (vb(id1,jd1,k,0,m) - (rho/zeta) * vb(id1,jd1,k,0,m+1))
              end if
              if (kd1>=0) then
                s = s + (ang_nxyz(k,ic)/(2*(zeta+eta))) * vb(i,jd1,kd1,0,m+1)
              end if
              vb(i,j,k,0,m) = s
            end do boot_j_order
          end do boot_j_i
        end do boot_j_k
      end do boot_j
      !
      !  Stage 4: Increment the fourth index. This is the general-case eq. 31
      !
      boot_l: do l=1,p2(4)
        call get_parent(l,ic,ld1)
        ll  = sum(ang_nxyz(l,:))
        ld2 = drop_xyz(ld1,ic)
        boot_l_k: do k=0,p2(3)
          kd1 = drop_xyz(k,ic)
          kl  = sum(ang_nxyz(k,:))
          boot_l_j: do j=0,p2(2)
            jd1 = drop_xyz(j,ic)
            jl  = sum(ang_nxyz(j,:))
            boot_l_i: do i=0,p2(1)
              id1 = drop_xyz(i,ic)
              il  = sum(ang_nxyz(i,:))
              boot_l_order: do m=0,nrec-(il+jl+kl+ll)
                s = (q(ic)-xyz(ic,4)) * vb(i,j,k,ld1,m)
                s = s - (rho/eta) * (q(ic)-p(ic)) * vb(i,j,k,ld1,m+1)
                if (ld2>=0) then
                  s = s + (ang_nxyz(ld1,ic)/(2*eta)) * (vb(i,j,k,ld2,m) - (rho/eta) * vb(i,j,k,ld2,m+1))
                end if
                if (kd1>=0) then
                  s = s + (ang_nxyz(k,ic)/(2*eta)) * (vb(i,j,kd1,ld1,m) - (rho/eta) * vb(i,j,kd1,ld1,m+1))
                end if
                if (id1>=0) then
                  s = s + (ang_nxyz(i,ic)/(2*(eta+zeta))) * vb(id1,j,k,ld1,m+1)
                end if
                if (jd1>=0) then
                  s = s + (ang_nxyz(j,ic)/(2*(eta+zeta))) * vb(i,jd1,k,ld1,m+1)
                end if
                vb(i,j,k,l,m) = s
              end do boot_l_order
            end do boot_l_i
          end do boot_l_j
        end do boot_l_k
      end do boot_l
      contains
      !
      !  Determine the index to recurse along and the parent integral
      !  This routine can be trivially made more efficient, but it probably does not matter (much)
      !
      subroutine get_parent(func,ic,id)
        integer(ik), intent(in)  :: func  ! Target index
        integer(ik), intent(out) :: ic    ! Coordinate to recurse along: 1, 2, or 3
        integer(ik), intent(out) :: id    ! Parent integral
        !
        ic = 1 ; id = drop_xyz(func,ic) ; if (id>=0) return
        ic = 2 ; id = drop_xyz(func,ic) ; if (id>=0) return
        ic = 3 ; id = drop_xyz(func,ic) ; if (id>=0) return
        stop 'import_gamess%os_2e_primitives - get_parent can''t get a parent!'
      end subroutine get_parent
!   end subroutine os_2e_primitives_real
