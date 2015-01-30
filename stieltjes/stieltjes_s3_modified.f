      subroutine stieltjes (num,e8,g8,printflag,overmax,filename,filint)
c Vitali Averbukh (2003). Comments to: vitali@tc.pci.uni-heidelberg.de
c This routine receives the sequence of the energy-gamma pairs and transforms it 
c into a shorter sequence of the Stieltjes energy-gamma pairs using the algorithm 
c of Muller-Plathe & Dierksen [Electronic Structure of Atoms, Molecules and Solids, 
c Proceedings of the II Escola Brasileira de Estruture Eletronica (World Scientific, 
c Singapore, 1990), p.1-29]. The tutorial is available on the net: 
c                              http://citeseer.nj.nec.com/465113.html
c Sequences of more than three points and up to NP points are allowed. 
c The recursive relations for the coefficients of the orthogonal polynomials
c are implemented in quadruple precision. The maximal polynomial order is determined 
c by the magnitude of the polynomial overlap approacing the polynomial norm.
c If only low-order (MAXORD < 5) approximation is available, the program produces the 
c interpolated Gamma(E) function. Otherwise, a sequence of Gamma(E) functions for the 
c approximation orders from 5 to MAXORD so that the convergence can be seen.
c The tridiagonal coefficient matrix is diagonalized by the REAL*16 version of the TQL2 
c routine (source available).
c The cumulative gamma is differentiated numerically.
c OVERMAX -  the maximal permitted norm-to-overlap ratio of the 
c            adjacent orthogonal polynomials
c Filenames added to Stieltjes routine so can save data from both
c length and velocity gauge. Filename must be given with quote marks
c either side

      implicit none

      integer num,nmax,npol
      integer np,maxord,ierr,imax,min,max,printflag
      integer k,n,i,j,ifail,converge,conv_max,iconv

c the maximal number of the energy points is NP
c the maximal order of polynomial in reduced moment problem is NMAX
      parameter (np=10000,nmax=1000)
      
      real*16 e_point(np),g_point(np) 
      real*16 e_new(nmax),g_new(nmax)
      real*16 bcoef(0:nmax),acoef(nmax),qpol(0:nmax,np)
      real*16 diag(nmax),offdiag(nmax),abvec(nmax,nmax)
      real*16 e_min,e_max,bprod,asum,qnorm,qoverlap,overmax

      real*8 e8(np),g8(np), e0(nmax),gamma(nmax)
      character*40 filename,filint !max length of filename string is 40 
c check the number of the input energy points 
      open(unit=2007,file=filename)
      open(unit=2008,file=filint)
      if (num.gt.np) then
         print*, '***WARNING*** Stieltjes: too many energy points'
         print11, 'NUM (',num,') > ',np
         return
      end if
 11   format(10x,a5,i3,a4,i3)
      if (num.le.3) then
         print*, '***ERROR*** Stieltjes: not enough energy points'
         print22, 'NUM = ',num
         return
      end if
 22   format(10x,a6,i3)

c transform the energy and gamma values into quadruple precision
      if (printflag.ne.0) write(6,*) 'Input energy points'
      do i=1,num
         e_point(i)=e8(i)
         g_point(i)=g8(i)
         if (printflag.ne.0) write(6,*) e8(i),g8(i),i
      end do
      if (printflag.ne.0) write(6,*)

c define the minimal energy and the maximal energies
      e_min=e_point(1)
      e_max=e_point(1)
      do i=2,num
         if( e_min.gt.e_point(i)) e_min=e_point(i)
         if( e_max.lt.e_point(i)) e_max=e_point(i)
      end do
      if (printflag.ne.0) print*, 'e_min =',e_min,'   e_max =',e_max

c initiate the recursive computation of the a,b coefficients and the orthogonal 
c polynomials according to (3.3.20-23) of Mueller-Plathe & Dierksen (1990)
       bcoef(0)=0.q0
       acoef(1)=0.q0
       do i=1,num
          bcoef(0)=bcoef(0)+g_point(i)
          acoef(1)=acoef(1)+g_point(i)/e_point(i)
       end do
       acoef(1)=acoef(1)/bcoef(0)

       do i=1,num
          qpol(0,i)=1.q0
          qpol(1,i)=1.q0/e_point(i)-acoef(1)
       end do

       bcoef(1)=0.q0
       acoef(2)=0.q0
       do i=1,num
          bcoef(1)=bcoef(1)+qpol(1,i)*g_point(i)/e_point(i)
          acoef(2)=acoef(2)+qpol(1,i)*g_point(i)/(e_point(i)**2)
       end do
       bcoef(1)= bcoef(1)/bcoef(0)
       acoef(2)=acoef(2)/(bcoef(0)*bcoef(1))-acoef(1)

c calculate the higher-order coefficients and polynomials recursively
c up to the (NUM-1)th order (total of NUM polynomials)
c if NUM > NMAX, we calculate only NMAX polynomials
       if (num.le.nmax) then
          npol=num
       else
          npol=nmax
       end if

       asum=acoef(1)
       do i=3,npol

          asum=asum+acoef(i-1)

          do j=1,num
             qpol(i-1,j)=(1.q0/e_point(j)-acoef(i-1))*qpol(i-2,j)
     &            -bcoef(i-2)*qpol(i-3,j)
          end do

          bprod=bcoef(0)
          do j=1,i-2
             bprod=bprod*bcoef(j)
          end do

          bcoef(i-1)=0.q0
          do j=1,num
             bcoef(i-1)=bcoef(i-1)+qpol(i-1,j)*g_point(j)
     &            /(e_point(j)**(i-1))
          end do
          bcoef(i-1)=bcoef(i-1)/bprod

          bprod=bprod*bcoef(i-1)

          acoef(i)=0.q0
          do j=1,num
             acoef(i)=acoef(i)+qpol(i-1,j)*g_point(j)/(e_point(j)**i)
          end do
          acoef(i)=acoef(i)/bprod-asum

c          print*, i,acoef(i),i-1,bcoef(i-1)

       end do

c calculate the NUM-th (NPOL-th) order polynomial just for the orthogonality check 
       do j=1,npol
          qpol(num,j)=(1.q0/e_point(j)-acoef(num))*qpol(num-1,j)
     &         -bcoef(num-1)*qpol(num-2,j)
       end do

c check the orthogonality of the polynomials to define the maximal approximation order 
c if the orthogonality is preserved for all orders, MAXORD is set to NPOL
       maxord=npol
       qnorm=bcoef(0)
       do i=1,npol
          qnorm=0.q0
          qoverlap=0.q0
          do j=1,num 
             qnorm=qnorm+qpol(i,j)**2*g_point(j)
             qoverlap=qoverlap+qpol(i,j)*qpol(i-1,j)*g_point(j)
          end do
          if (qabs(qoverlap).lt.1.q-50) qoverlap=1.q-50
c          print*, i,qoverlap,qnorm
          if (qnorm/qabs(qoverlap).le.overmax) then
c MAXORD=I-1 is appropriate since the polynomial failing 
c the orthogonality check should not be used
             maxord=i-1
             go to 10
          end if
       end do

 10    continue



c look how many Stieltjes orders are available
       if (maxord.lt.5) then
          min=maxord
          max=maxord
          print*, '***WARNING*** Stieltjes:' 
          print*, ' only very low-order approximation is available'
          print*, ' MAXORD=',maxord
       else
          min=5
          max=maxord
       end if


c perform the gamma calculation using the successive approximations 
c N=5,...,MAXORD (if MAXORD > 5)
c       write(2007,*) max
        write(2008,*) max
       do imax=min,max

c fill the coefficients matrix
          do i=1,imax
             diag(i)=acoef(i)
          end do
          do i=2,imax
             offdiag(i)=-qsqrt(bcoef(i-1))
          end do

c diagonalize the coefficients matrix
c initialize the arrays
          do i=1,nmax
             do j=1,nmax
                abvec(i,j)=0.q0
             end do
            abvec(i,i)=1.q0
          end do
          call tql2(nmax,imax,diag,offdiag,abvec,ierr)
          if (ierr.ne.0) then
             print*, '***WARNING*** Stieltjes:'
             print*, ' the eigenvalue no. ',ierr,' failed to converge'
          end if
c fill the Stieltjes energy and gamma arrays
c note that the eigenvalues are inverse energies and are given in ascending order 
          do i=1,imax
             e_new(i)=1.q0/diag(imax+1-i)
             g_new(i)=bcoef(0)*abvec(1,imax+1-i)**2
c             print*, i,e_new(i),g_new(i)
          end do

c calculate the gamma's by simple numerical differentiation at the middle 
c point of each [E_NEW(I),E_NEW(I+1)] interval
c          write(2007,*) imax
          write(2008,*) imax 
          write(2007,*) " "
          do i=1,imax-1
             e0(i)=(e_new(i)+e_new(i+1))/2.d0
             gamma(i)=(g_new(i+1)+g_new(i))/(2.d0*(e_new(i+1)-e_new(i)))
             write(2007,*) e0(i)*27.211396,gamma(i)
             write(2008,*) e0(i)*27.211396,gamma(i)
         end do

      end do
      close(2007)
      close(2008)

      return
      end

