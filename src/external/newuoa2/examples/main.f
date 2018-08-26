c
c  Test the algorithm using More's Dfo test problem library
c  in the paper "Benchmarking Derivative-Free Optimzation Algorithms"
c  Problem downloaded from www.mcs.anl.gov/~more/dfo
c
      program mjdp_h

      integer nmax, nptmax, nspace
      parameter (nmax=100,nptmax=2*nmax+1,
     &     nspace=(nptmax+11)*(nptmax+nmax)+nmax*(3*nmax+11)/2)
      double precision w(nspace)
      integer version  
C
C     THE ARGUMENT N DENOTES THE NUMBER OF UNKNOWNS 
C     THE ARGUMENT X HAS THE STARTING GUESS FOR THE UNKNOWNS
C
C     SET THINGS UP AND CALL NEWOA() FUNCTION
C
      INTEGER N, IPRINT, MAXFUN, NPT, MV, nprob, xs
      DOUBLE PRECISION X(nmax), e_noise, factor
      DOUBLE PRECISION RHOBEG, RHOEND 
      common /test/ nprob, mv
      common /noise/e_noise
      integer nread, nwrite
      parameter (nread = 1, nwrite = 2 )

      integer mmax
      parameter (mmax = 400 )
      double precision v_err(mmax), f, f1
c
      IPRINT = 2
c   Maxfun was set later to be 100*(n+1)
c      MAXFUN = 100000
      RHOBEG = 1.d0
      RHOEND = 1.0D-8
c
c  version = 1, Original Mjdp
c  version = 2, Modified Mjdp
c
      version = 2
c
c  noise level
c
      e_noise = 0.d-3

      open (nread,file='dfo.dat',status='old')
      open (nwrite,file='dfo.info')
 
      do while (1 .eq. 1) 
c
c  read problem data nprob, n, mv, xs
c  nprob: Problem #,        n: dimension, 
c  mv: number of equations, xs: 10^{xs}*standard starting point
c 
      read (nread,*) nprob, n, mv, xs 
c      print *,  nprob, n, mv, xs

      if (nprob.eq.0) then
         close(nread)
         close(nwrite)
         stop
      endif
C Set Maxfun = 100(n+1) 
C i.e. the solver can use at most 200 (simplex) gradients
      MAXFUN = 100*(n+1)
C
C dfoxs provides the starting point X
C
      if (xs.ne.0) then
         factor = 1.d1**xs
      else
         factor = 1.d0
      endif

      call dfoxs(n,x,nprob,factor)

      NPT = 2*N+1
      
      PRINT 20, nprob, N, MV, NPT, xs
   20 FORMAT (//4X,"Results with NTEST=",I5, ' N =',I5,
     &       ' MV=',I5,' and NPT =',I5, ' xs=', I2)

      if (version.eq.1) then

C        CALL NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)

      elseif (version.eq.2) then
c
c  Check if the problem is least square problem or not.
c       
c        call dfovec(n, mv, x, v_err)
c        call f_value(mv,v_err,f)
c        call CALFUN (N,X,f1)
c        print *, f, f1
c        if (dabs(f-f1).gt.1.d-15) then
c          print *, "It is not a least square problem,"
c          print *, "Can not use the modified version right now."
c          print *, "Have to use use original Mjdp. Set version=1."
c          stop
c        endif

        CALL NEWUOA_H(N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,MV)

      endif
C
C Calculate the final function value again 
C
        call dfovec(n, mv, x, v_err)
        f = 0.d0
        do i=1,mv
          f = f + v_err(i)**2
        enddo
        print *, "Final function value, f_final=", f
 
      enddo

      END

      subroutine dfoxs(n,x,nprob,factor)
      integer n, nprob
      double precision factor
      double precision x(n)
c     **********
c
c     Subroutine dfoxs
c
c     This subroutine specifies the standard starting points for the
c     functions defined by subroutine dfovec. The subroutine returns
c     in x a multiple (factor) of the standard starting point. 
c
c     The subroutine statement is
c
c       subroutine dfoxs(n,x,nprob,factor)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an output array of length n which contains the standard
c         starting point for problem nprob multiplied by factor.
c
c       nprob is a positive integer input variable which defines the
c         number of the problem. nprob must not exceed 22.
c
c       factor is an input variable which specifies the multiple of
c         the standard starting point.
c
c     Argonne National Laboratory. 
c     Jorge More' and Stefan Wild. September 2007.
c
c     **********
      integer i, j
      double precision sum, temp

c     Selection of initial point.

      if (nprob .le. 3) then

c        Linear function - full rank or rank 1.

         do j = 1, n
            x(j) = 1.0d0
         end do
      
      else if (nprob .eq. 4) then

c        Rosenbrock function.

         x(1) = -1.2d0
         x(2) = 1.0d0

      else if (nprob .eq. 5) then

c        Helical valley function.

         x(1) = -1.0d0
         x(2) = 0.0d0
         x(3) = 0.0d0

      else if (nprob .eq. 6) then

c        Powell singular function.

         x(1) = 3.0d0
         x(2) = -1.0d0
         x(3) = 0.0d0
         x(4) = 1.0d0

      else if (nprob .eq. 7) then

c        Freudenstein and Roth function.

         x(1) = 0.5d0
         x(2) = -2.0d0

      else if (nprob .eq. 8) then

c        Bard function.

         x(1) = 1.0d0
         x(2) = 1.0d0
         x(3) = 1.0d0

      else if (nprob .eq. 9) then

c        Kowalik and Osborne function.

         x(1) = 0.25d0
         x(2) = 0.39d0
         x(3) = 0.415d0
         x(4) = 0.39d0

      else if (nprob .eq. 10) then

c        Meyer function.

        x(1) = 0.02d0
        x(2) = 4000.0d0
        x(3) = 250.0d0

      else if (nprob .eq. 11) then

c        Watson function.

         do j = 1, n
            x(j) = 0.5d0
         end do

      else if (nprob .eq. 12) then

c        Box 3-dimensional function.

         x(1) = 0.0d0
         x(2) = 10.0d0
         x(3) = 20.0d0

      else if (nprob .eq. 13) then
 
c        Jennrich and Sampson function.

         x(1) = 0.3d0
         x(2) = 0.4d0

      else if (nprob .eq. 14) then

c        Brown and Dennis function.

         x(1) = 25.0d0
         x(2) = 5.0d0
         x(3) = -5.0d0
         x(4) = -1.0d0

      else if (nprob .eq. 15) then

c        Chebyquad function.

         do j = 1, n
            x(j) = j/dble(n+1)
         end do

      else if (nprob .eq. 16) then

c        Brown almost-linear function.

         do j = 1, n
            x(j) = 0.5d0
         end do

      else if (nprob .eq. 17) then

c        Osborne 1 function.

         x(1) = 0.5d0
         x(2) = 1.5d0
         x(3) = 1.0d0
         x(4) = 0.01d0
         x(5) = 0.02d0

      else if (nprob .eq. 18) then

c        Osborne 2 function.

         x(1) = 1.3d0
         x(2) = 0.65d0
         x(3) = 0.65d0
         x(4) = 0.7d0
         x(5) = 0.6d0
         x(6) = 3.0d0
         x(7) = 5.0d0
         x(8) = 7.0d0
         x(9) = 2.0d0
         x(10) = 4.5d0
         x(11) = 5.5d0

      else if (nprob .eq. 19) then

c        Bdqrtic function.

         do j = 1, n
            x(j) = 1.0d0
         end do
         
      else if (nprob .eq. 20) then

c        Cube function.

         do j = 1, n
            x(j) = 0.5d0
         end do
                  
      else if (nprob .eq. 21) then

c        Mancino function.

         do i = 1, n
            sum = 0.0d0 
            do j = 1, n
               temp = sqrt(dble(i)/dble(j))
               sum = sum + 
     *               temp*((sin(log(temp)))**5+(cos(log(temp)))**5)
            end do
            x(i) = -8.710996D-4*((i-50)**3 + sum)
         end do
         
      else if (nprob .eq. 22) then

c        Heart8 function.

         x(1) = -0.3d0
         x(2) = -0.39d0
         x(3) = 0.3d0
         x(4) = -0.344d0
         x(5) = -1.2d0
         x(6) = 2.69d0
         x(7) = 1.59d0
         x(8) = -1.5d0
                           
      else
         write (*,*) "Parameter nprob > 22 in subroutine dfoxs"
      end if

c     Compute multiple of initial point.

      do j = 1, n
         x(j) = factor*x(j)
      end do

      return

      end
