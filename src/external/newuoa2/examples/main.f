c
c  Test the algorithm using More's Dfo test problem library
c  in the paper "Benchmarking Derivative-Free Optimzation Algorithms"
c  Problem downloaded from www.mcs.anl.gov/~more/dfo
c
      program mjdp_h

      use, intrinsic :: iso_fortran_env
      use newuoa2_real64

      implicit none

      integer, parameter :: PREC = real64

      integer nmax, nptmax, nspace
      parameter (nmax=100,nptmax=2*nmax+1,
     &     nspace=(nptmax+11)*(nptmax+nmax)+nmax*(3*nmax+11)/2)
      real (PREC) w(nspace)
      integer version  
C
C     THE ARGUMENT N DENOTES THE NUMBER OF UNKNOWNS 
C     THE ARGUMENT X HAS THE STARTING GUESS FOR THE UNKNOWNS
C
C     SET THINGS UP AND CALL NEWOA() FUNCTION
C
      INTEGER N, IPRINT, MAXFUN, NPT, MV, nprob, xs, i
      real (PREC) :: X(nmax), e_noise, factor
      real (PREC) :: RHOBEG, RHOEND
      common /test/ nprob, mv
      common /noise/e_noise
      integer nread, nwrite
      parameter (nread = 1, nwrite = 2 )

      integer mmax
      parameter (mmax = 400 )
      real (PREC) :: v_err(mmax), f, f1
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
         return
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

        CALL newuoa2 (FOBJ, N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,MV)

      endif
C
C Calculate the final function value again 
C
        call fobj (n, mv, x, v_err)
        f = 0.d0
        do i=1,mv
          f = f + v_err(i)**2
        enddo
        print *, "Final function value, f_final=", f
 
      enddo

      contains

      subroutine fobj (n, m, x, fx)
c       Objective functionw wrapper: selects the appropriate objective to be
c       called based on value of NPROB.
        integer, intent(in) :: n, m
        real (PREC), intent(in), dimension(*) :: x
        real (PREC), intent(out), dimension(*) :: fx

c       Note: value of NPROB is taken from outer scope
c       Note: DFOVEC function flips the order of N, M, expects X(N), FX(M)
        call dfovec (m, n, x, fx, nprob)

      end subroutine

      END

