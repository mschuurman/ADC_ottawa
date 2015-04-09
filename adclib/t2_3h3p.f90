 real(d) function t2_3h3p_hc(apr,bpr,cpr,kpr,lpr,mpr)

 integer, intent(in) :: kpr,apr,bpr,lpr,cpr,mpr


!!! I only have to write the occupied-unoccupied part of the density matrix !!!
!!!!!!!   THE APR OR KPR HAVE ONLY ONE POSSIBLITY FOR THE SPIN
!!!!!!!   THE FACTOR OF TWO IS TAKEN INTO ACCOUNT IN THE CALLING
!!!!!!!   BY THE FUNCTIONS D2_1,D2_2,D2_3,D2_4.
!!!!!!!   THE FACTOR OF 2 TAKES INTO ACCOUHT THAT OR APR OR KPR
!!!!!!!   CAN TAKE ONE TIME ALPHA AND ONE TIME BETA VALUES AND
!!!!!!!   THE EXPRESSION DOES NOT DEPEND OF WHICH VALUE IT TAKES



    integer :: d,f,nsym1,nsym2,nsym3,nsym4,d1,f1,cnt
    real*8 :: DA,eabcklm,eterm,term





!!! part with summation over OCCUPIED index





      do f1=1,nOcc
        f=roccnum(f1)



!!! first term



             nsym1=MT(orbSym(apr),orbSym(bpr))
             nsym2=MT(orbSym(kpr),orbSym(f))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(cpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0



                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)



                eterm = e(kpr)+e(f)-e(apr)-e(bpr)


                DA = eabcklm*eterm

                term=term+vpqrs(cpr,lpr,f,mpr)*(+2.0*vpqrs(apr,kpr,bpr,f)-1.0*vpqrs(apr,f,bpr,kpr))
                term=term+vpqrs(cpr,mpr,f,lpr)*(-1.0*vpqrs(apr,kpr,bpr,f)+2.0*vpqrs(apr,f,bpr,kpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term


              end if



!!! second term (cpr and bpr) with respect to the first



             nsym1=MT(orbSym(apr),orbSym(cpr))
             nsym2=MT(orbSym(kpr),orbSym(f))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(bpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0



                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)



                eterm = e(kpr)+e(f)-e(apr)-e(cpr)


                DA = eabcklm*eterm

                term=term+vpqrs(bpr,lpr,f,mpr)*(+2.0*vpqrs(apr,kpr,cpr,f)-1.0*vpqrs(apr,f,cpr,kpr))
                term=term+vpqrs(bpr,mpr,f,lpr)*(-1.0*vpqrs(apr,kpr,cpr,f)+2.0*vpqrs(apr,f,cpr,kpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc-term


              end if



!!! third term (apr and bpr) with respect to the second



             nsym1=MT(orbSym(bpr),orbSym(cpr))
             nsym2=MT(orbSym(kpr),orbSym(f))
             nsym3=MT(orbSym(lpr),orbSym(mpr))
             nsym4=MT(orbSym(apr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0



                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)



                eterm = e(kpr)+e(f)-e(bpr)-e(cpr)


                DA = eabcklm*eterm

                term=term+vpqrs(apr,lpr,f,mpr)*(+2.0*vpqrs(bpr,kpr,cpr,f)-1.0*vpqrs(bpr,f,cpr,kpr))
                term=term+vpqrs(apr,mpr,f,lpr)*(-1.0*vpqrs(bpr,kpr,cpr,f)+2.0*vpqrs(bpr,f,cpr,kpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term



              end if




! fourth term (kpr and lpr) with respect to the first



             nsym1=MT(orbSym(apr),orbSym(bpr))
             nsym2=MT(orbSym(lpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(cpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0



                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)



                eterm = e(lpr)+e(f)-e(apr)-e(bpr)


                DA = eabcklm*eterm

                term=term+vpqrs(cpr,kpr,f,mpr)*(+2.0*vpqrs(apr,lpr,bpr,f)-1.0*vpqrs(apr,f,bpr,lpr))
                term=term+vpqrs(cpr,mpr,f,kpr)*(-1.0*vpqrs(apr,lpr,bpr,f)+2.0*vpqrs(apr,f,bpr,lpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc-term


              end if


! fifth term (bpr and cpr) with respect to the fourth




             nsym1=MT(orbSym(apr),orbSym(cpr))
             nsym2=MT(orbSym(lpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(bpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0



                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)



                eterm = e(lpr)+e(f)-e(apr)-e(cpr)


                DA = eabcklm*eterm

                term=term+vpqrs(bpr,kpr,f,mpr)*(+2.0*vpqrs(apr,lpr,cpr,f)-1.0*vpqrs(apr,f,cpr,lpr))
                term=term+vpqrs(bpr,mpr,f,kpr)*(-1.0*vpqrs(apr,lpr,cpr,f)+2.0*vpqrs(apr,f,cpr,lpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term


              end if


! sixth term (apr and bpr) with respect to the fifth





             nsym1=MT(orbSym(bpr),orbSym(cpr))
             nsym2=MT(orbSym(lpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(mpr))
             nsym4=MT(orbSym(apr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0



                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)



                eterm = e(lpr)+e(f)-e(bpr)-e(cpr)


                DA = eabcklm*eterm

                term=term+vpqrs(apr,kpr,f,mpr)*(+2.0*vpqrs(bpr,lpr,cpr,f)-1.0*vpqrs(bpr,f,cpr,lpr))
                term=term+vpqrs(apr,mpr,f,kpr)*(-1.0*vpqrs(bpr,lpr,cpr,f)+2.0*vpqrs(bpr,f,cpr,lpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc-term


              end if


! seventh term (lpr and mpr) with respect to the fourth




             nsym1=MT(orbSym(apr),orbSym(bpr))
             nsym2=MT(orbSym(mpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(lpr))
             nsym4=MT(orbSym(cpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0



                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)



                eterm = e(mpr)+e(f)-e(apr)-e(bpr)


                DA = eabcklm*eterm

                term=term+vpqrs(cpr,kpr,f,lpr)*(+2.0*vpqrs(apr,mpr,bpr,f)-1.0*vpqrs(apr,f,bpr,mpr))
                term=term+vpqrs(cpr,lpr,f,kpr)*(-1.0*vpqrs(apr,mpr,bpr,f)+2.0*vpqrs(apr,f,bpr,mpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term




              end if


! eighth term (lpr and mpr) with respect to the fifth





             nsym1=MT(orbSym(apr),orbSym(cpr))
             nsym2=MT(orbSym(mpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(lpr))
             nsym4=MT(orbSym(bpr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0



                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)



                eterm = e(mpr)+e(f)-e(apr)-e(cpr)


                DA = eabcklm*eterm

                term=term+vpqrs(bpr,kpr,f,lpr)*(+2.0*vpqrs(apr,mpr,cpr,f)-1.0*vpqrs(apr,f,cpr,mpr))
                term=term+vpqrs(bpr,lpr,f,kpr)*(-1.0*vpqrs(apr,mpr,cpr,f)+2.0*vpqrs(apr,f,cpr,mpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc-term



              end if


! nineth term (lpr and mpr) with respect to the sixth




             nsym1=MT(orbSym(bpr),orbSym(cpr))
             nsym2=MT(orbSym(mpr),orbSym(f))
             nsym3=MT(orbSym(kpr),orbSym(lpr))
             nsym4=MT(orbSym(apr),orbSym(f))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0



                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)



                eterm = e(mpr)+e(f)-e(bpr)-e(cpr)


                DA = eabcklm*eterm

                term=term+vpqrs(apr,kpr,f,lpr)*(+2.0*vpqrs(bpr,mpr,cpr,f)-1.0*vpqrs(bpr,f,cpr,mpr))
                term=term+vpqrs(apr,lpr,f,kpr)*(-1.0*vpqrs(bpr,mpr,cpr,f)+2.0*vpqrs(bpr,f,cpr,mpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term


             end if





      end do








!!! part with summation over  VIRTUAL  index







      do d1=nOcc+1,nBas
         d=roccnum(d1)



!!! first term





             nsym1=MT(orbSym(mpr),orbSym(d))
             nsym2=MT(orbSym(bpr),orbSym(cpr))
             nsym3=MT(orbSym(lpr),orbSym(kpr))
             nsym4=MT(orbSym(apr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0




                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)




                eterm = e(kpr)+e(lpr)-e(apr)-e(d)


                DA = eabcklm*eterm

                term=term+vpqrs(bpr,d,cpr,mpr)*(+2.0*vpqrs(apr,kpr,d,lpr)-1.0*vpqrs(apr,lpr,d,kpr))
                term=term+vpqrs(bpr,mpr,cpr,d)*(-1.0*vpqrs(apr,kpr,d,lpr)+2.0*vpqrs(apr,lpr,d,kpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term  

             end if





!!! second term (lpr and mpr) with respect to the first





             nsym1=MT(orbSym(lpr),orbSym(d))
             nsym2=MT(orbSym(bpr),orbSym(cpr))
             nsym3=MT(orbSym(mpr),orbSym(kpr))
             nsym4=MT(orbSym(apr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0




                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)




                eterm = e(kpr)+e(mpr)-e(apr)-e(d)


                DA = eabcklm*eterm

                term=term+vpqrs(bpr,d,cpr,lpr)*(+2.0*vpqrs(apr,kpr,d,mpr)-1.0*vpqrs(apr,mpr,d,kpr))
                term=term+vpqrs(bpr,lpr,cpr,d)*(-1.0*vpqrs(apr,kpr,d,mpr)+2.0*vpqrs(apr,mpr,d,kpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc-term  ! this part has minus sign

             end if



!!! third  term  (kpr and lpr) with respect to the second


             nsym1=MT(orbSym(kpr),orbSym(d))
             nsym2=MT(orbSym(bpr),orbSym(cpr))
             nsym3=MT(orbSym(mpr),orbSym(lpr))
             nsym4=MT(orbSym(apr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0




                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)




                eterm = e(lpr)+e(mpr)-e(apr)-e(d)


                DA = eabcklm*eterm

                term=term+vpqrs(bpr,d,cpr,kpr)*(+2.0*vpqrs(apr,lpr,d,mpr)-1.0*vpqrs(apr,mpr,d,lpr))
                term=term+vpqrs(bpr,kpr,cpr,d)*(-1.0*vpqrs(apr,lpr,d,mpr)+2.0*vpqrs(apr,mpr,d,lpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term  

             end if



!!! fourth term (apr and bpr) with respect to the first



             nsym1=MT(orbSym(mpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(cpr))
             nsym3=MT(orbSym(lpr),orbSym(kpr))
             nsym4=MT(orbSym(bpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0




                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)




                eterm = e(kpr)+e(lpr)-e(bpr)-e(d)


                DA = eabcklm*eterm

                term=term+vpqrs(apr,d,cpr,mpr)*(+2.0*vpqrs(bpr,kpr,d,lpr)-1.0*vpqrs(bpr,lpr,d,kpr))
                term=term+vpqrs(apr,mpr,cpr,d)*(-1.0*vpqrs(bpr,kpr,d,lpr)+2.0*vpqrs(bpr,lpr,d,kpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc-term   ! this part as a minus sign  

             end if



!!! fifth term (apr and bpr) with respect to the second



             nsym1=MT(orbSym(lpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(cpr))
             nsym3=MT(orbSym(mpr),orbSym(kpr))
             nsym4=MT(orbSym(bpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0




                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)




                eterm = e(kpr)+e(mpr)-e(bpr)-e(d)


                DA = eabcklm*eterm

                term=term+vpqrs(apr,d,cpr,lpr)*(+2.0*vpqrs(bpr,kpr,d,mpr)-1.0*vpqrs(bpr,mpr,d,kpr))
                term=term+vpqrs(apr,lpr,cpr,d)*(-1.0*vpqrs(bpr,kpr,d,mpr)+2.0*vpqrs(bpr,mpr,d,kpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term  

             end if


!!! sixth term (apr and bpr) with respect to the third



             nsym1=MT(orbSym(kpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(cpr))
             nsym3=MT(orbSym(mpr),orbSym(lpr))
             nsym4=MT(orbSym(bpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0




                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)




                eterm = e(lpr)+e(mpr)-e(bpr)-e(d)


                DA = eabcklm*eterm

                term=term+vpqrs(apr,d,cpr,kpr)*(+2.0*vpqrs(bpr,lpr,d,mpr)-1.0*vpqrs(bpr,mpr,d,lpr))
                term=term+vpqrs(apr,kpr,cpr,d)*(-1.0*vpqrs(bpr,lpr,d,mpr)+2.0*vpqrs(bpr,mpr,d,lpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc-term  ! this part has minus sign

             end if




!!! seventh term (bpr and cpr) with respect to the fourth






             nsym1=MT(orbSym(mpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(bpr))
             nsym3=MT(orbSym(lpr),orbSym(kpr))
             nsym4=MT(orbSym(cpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0




                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)




                eterm = e(kpr)+e(lpr)-e(cpr)-e(d)


                DA = eabcklm*eterm

                term=term+vpqrs(apr,d,bpr,mpr)*(+2.0*vpqrs(cpr,kpr,d,lpr)-1.0*vpqrs(cpr,lpr,d,kpr))
                term=term+vpqrs(apr,mpr,bpr,d)*(-1.0*vpqrs(cpr,kpr,d,lpr)+2.0*vpqrs(cpr,lpr,d,kpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term  

             end if



!!! eighth term (lpr and mpr) with respect to the seventh






             nsym1=MT(orbSym(lpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(bpr))
             nsym3=MT(orbSym(mpr),orbSym(kpr))
             nsym4=MT(orbSym(cpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0




                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)




                eterm = e(kpr)+e(mpr)-e(cpr)-e(d)


                DA = eabcklm*eterm

                term=term+vpqrs(apr,d,bpr,lpr)*(+2.0*vpqrs(cpr,kpr,d,mpr)-1.0*vpqrs(cpr,mpr,d,kpr))
                term=term+vpqrs(apr,lpr,bpr,d)*(-1.0*vpqrs(cpr,kpr,d,mpr)+2.0*vpqrs(cpr,mpr,d,kpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc-term   ! this part as a minus sign  

             end if



!!! nineth term  (kpr and lpr) with respect to the eighth




             nsym1=MT(orbSym(kpr),orbSym(d))
             nsym2=MT(orbSym(apr),orbSym(bpr))
             nsym3=MT(orbSym(mpr),orbSym(lpr))
             nsym4=MT(orbSym(cpr),orbSym(d))

             if((MT(nsym1,nsym2) .eq. 1) .and. (MT(nsym3,nsym4) .eq. 1)) then

                cnt=cnt+1
                term=0.0




                eabcklm = -e(apr)-e(bpr)-e(cpr)+e(kpr)+e(lpr)+e(mpr)




                eterm = e(lpr)+e(mpr)-e(cpr)-e(d)


                DA = eabcklm*eterm

                term=term+vpqrs(apr,d,bpr,kpr)*(+2.0*vpqrs(cpr,lpr,d,mpr)-1.0*vpqrs(cpr,mpr,d,lpr))
                term=term+vpqrs(apr,kpr,bpr,d)*(-1.0*vpqrs(cpr,lpr,d,mpr)+2.0*vpqrs(cpr,mpr,d,lpr))
                term=term/DA

                t2_3h3p_hc=t2_3h3p_hc+term  

             end if



      end do




  end function t2_3h3p_hc




