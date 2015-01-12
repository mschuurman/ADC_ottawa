      type(gam_structure), pointer :: l, r
      integer(ik)                  :: ir, ic, ic_p1, ic_pn
      integer(ik), parameter       :: cpp = 10               ! Columns per page
      logical                      :: usegfmt
      character(len=20)            :: sym, head
      character(len=20)            :: lab
      !
      call TimerStart('GAMESS print 1e')
      l => gam_def ; r => gam_def
      if (present(bra)) l => bra
      if (present(ket)) r => ket
      !
      !  Symmetry check
      !
      sym = 'HERMITIAN'
      if (present(symmetry)) sym = symmetry
      if (size(v,dim=1)==size(v,dim=2)) then
        write (out,"('Expected integral symmetry: ',a)") trim(sym)
        select case (sym)
          case default
            write (out,"('Symmetry code ',a,' is not recognized')") trim(sym)
            stop 'import_gamess%gamess_print_1e_integrals_real - bad argument'
          case ('NONE')
          case ('HERMITIAN')
            write (out,"('    Maximum deviation from symmetry = ',g12.5)") maxval(abs(v-transpose(v)))
            write (out,"('         Maximum deviation at index = ',2i6  )") maxloc(abs(v-transpose(v)))
          case ('ANTI-HERMITIAN')
            write (out,"('    Maximum deviation from symmetry = ',g12.5)") maxval(abs(v+transpose(v)))
            write (out,"('         Maximum deviation at index = ',2i6  )") maxloc(abs(v+transpose(v)))
        end select
      end if
            write (out,"('             Largest matrix element = ',g12.5)") maxval(v)
            write (out,"(' Largest matrix element is at index = ',2i6)") maxloc(v)
            write (out,"('            Smallest matrix element = ',g12.5)") minval(v)
            write (out,"('Smallest matrix element is at index = ',2i6)") minloc(v)
      !
      head = 'BOTH'
      if (present(heading)) head = heading
      usegfmt = any(abs(v)>=1e5) .or. all(abs(v)<=1e-3)
      print_pages: do ic_p1=1,r%nbasis,cpp
        ic_pn = min(r%nbasis,ic_p1+cpp-1)
        write (out,"(t17,10(1x,i12))") (ic,               ic=ic_p1,ic_pn)
        if (head=='BOTH') write (out,"(t23,10(1x,a12))") (r%bas_labels(ic), ic=ic_p1,ic_pn)
        print_rows: do ir=1,l%nbasis
          lab = ' '
          if (head/='NONE') lab = l%bas_labels(ir)
          if (usegfmt) then
            write (out,"(1x,i5,1x,a11,1x,10(1x,g12.5))") ir, lab, v(ir,ic_p1:ic_pn)
          else
            write (out,"(1x,i5,1x,a11,1x,10(1x,f12.5))") ir, lab, v(ir,ic_p1:ic_pn)
          endif
        end do print_rows
        write (out,"()")
      end do print_pages
      !
      call TimerStop('GAMESS print 1e')
