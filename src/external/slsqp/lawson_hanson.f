C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
C
C  CONSTRUCTION AND/OR APPLICATION OF A SINGLE
C  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B
C
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
c                     Subroutine Arguments
c
C     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
c            Householder transformation, or Algorithm H2 to apply a
c            previously constructed transformation.
C     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT.
C     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO
C            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
c            vector.  IUE is the storage increment between elements.
c            On exit when MODE = 1, U() and UP contain quantities
c            defining the vector U of the Householder transformation.
c            on entry with MODE = 2, U() and UP should contain
c            quantities previously computed with MODE = 1.  These will
c            not be modified during the entry with MODE = 2.
C     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
c            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
c            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
c            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
C     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
C     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().
C     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0
C            NO OPERATIONS WILL BE DONE ON C().
C     ------------------------------------------------------------------
      SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)
C     ------------------------------------------------------------------
      integer I, I2, I3, I4, ICE, ICV, INCR, IUE, J
      integer L1, LPIVOT, M, MODE, NCV
      double precision B, C(*), CL, CLINV, ONE, SM
c     double precision U(IUE,M)
      double precision U(IUE,*)
      double precision UP
      parameter(ONE = 1.0d0)
C     ------------------------------------------------------------------
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
      CL=abs(U(1,LPIVOT))
      IF (MODE.EQ.2) GO TO 60
C                            ****** CONSTRUCT THE TRANSFORMATION. ******
          DO 10 J=L1,M
   10     CL=MAX(abs(U(1,J)),CL)
      IF (CL) 130,130,20
   20 CLINV=ONE/CL
      SM=(U(1,LPIVOT)*CLINV)**2
          DO 30 J=L1,M
   30     SM=SM+(U(1,J)*CLINV)**2
      CL=CL*SQRT(SM)
      IF (U(1,LPIVOT)) 50,50,40
   40 CL=-CL
   50 UP=U(1,LPIVOT)-CL
      U(1,LPIVOT)=CL
      GO TO 70
C            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C
   60 IF (CL) 130,130,70
   70 IF (NCV.LE.0) RETURN
      B= UP*U(1,LPIVOT)
C                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C
      IF (B) 80,130,130
   80 B=ONE/B
      I2=1-ICV+ICE*(LPIVOT-1)
      INCR=ICE*(L1-LPIVOT)
          DO 120 J=1,NCV
          I2=I2+ICV
          I3=I2+INCR
          I4=I3
          SM=C(I2)*UP
              DO 90 I=L1,M
              SM=SM+C(I3)*U(1,I)
   90         I3=I3+ICE
          IF (SM) 100,120,100
  100     SM=SM*B
          C(I2)=C(I2)+SM*UP
              DO 110 I=L1,M
              C(I4)=C(I4)+SM*U(1,I)
  110         I4=I4+ICE
  120     CONTINUE
  130 RETURN
      END

      SUBROUTINE HFTI (A,MDA,M,N,B,MDB,NB,TAU,KRANK,RNORM,H,G,IP)
c
C  SOLVE LEAST SQUARES PROBLEM USING ALGORITHM, HFTI.
c  Householder Forward Triangulation with column Interchanges.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
      integer I, II, IP1, J, JB, JJ, K, KP1, KRANK
      integer L, LDIAG, LMAX, M, MDA, MDB, N, NB
c     integer IP(N)
c     double precision A(MDA,N),B(MDB,NB),H(N),G(N),RNORM(NB)
      integer IP(*)
      double precision A(MDA,*),B(MDB, *),H(*),G(*),RNORM( *)
      double precision FACTOR, HMAX, SM, TAU, TMP, ZERO
      parameter(FACTOR = 0.001d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
C
      K=0
      LDIAG=min(M,N)
      IF (LDIAG.LE.0) GO TO 270
          DO 80 J=1,LDIAG
          IF (J.EQ.1) GO TO 20
C
C     UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
C    ..
          LMAX=J
              DO 10 L=J,N
              H(L)=H(L)-A(J-1,L)**2
              IF (H(L).GT.H(LMAX)) LMAX=L
   10         CONTINUE
          IF(DIFF(HMAX+FACTOR*H(LMAX),HMAX)) 20,20,50
C
C     COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
C    ..
   20     LMAX=J
              DO 40 L=J,N
              H(L)=0.
                  DO 30 I=J,M
   30             H(L)=H(L)+A(I,L)**2
              IF (H(L).GT.H(LMAX)) LMAX=L
   40         CONTINUE
          HMAX=H(LMAX)
C    ..
C     LMAX HAS BEEN DETERMINED
C
C     DO COLUMN INTERCHANGES IF NEEDED.
C    ..
   50     CONTINUE
          IP(J)=LMAX
          IF (IP(J).EQ.J) GO TO 70
              DO 60 I=1,M
              TMP=A(I,J)
              A(I,J)=A(I,LMAX)
   60         A(I,LMAX)=TMP
          H(LMAX)=H(J)
C
C     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A AND B.
C    ..
   70     CALL H12 (1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,N-J)
   80     CALL H12 (2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
C
C     DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU.
C    ..
          DO 90 J=1,LDIAG
          IF (ABS(A(J,J)).LE.TAU) GO TO 100
   90     CONTINUE
      K=LDIAG
      GO TO 110
  100 K=J-1
  110 KP1=K+1
C
C     COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
C
      IF (NB.LE.0) GO TO 140
          DO 130 JB=1,NB
          TMP=ZERO
          IF (KP1.GT.M) GO TO 130
              DO 120 I=KP1,M
  120         TMP=TMP+B(I,JB)**2
  130     RNORM(JB)=SQRT(TMP)
  140 CONTINUE
C                                           SPECIAL FOR PSEUDORANK = 0
      IF (K.GT.0) GO TO 160
      IF (NB.LE.0) GO TO 270
          DO 150 JB=1,NB
              DO 150 I=1,N
  150         B(I,JB)=ZERO
      GO TO 270
C
C     IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
C     DECOMPOSITION OF FIRST K ROWS.
C    ..
  160 IF (K.EQ.N) GO TO 180
          DO 170 II=1,K
          I=KP1-II
  170     CALL H12 (1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  180 CONTINUE
C
C
      IF (NB.LE.0) GO TO 270
          DO 260 JB=1,NB
C
C     SOLVE THE K BY K TRIANGULAR SYSTEM.
C    ..
              DO 210 L=1,K
              SM=ZERO
              I=KP1-L
              IF (I.EQ.K) GO TO 200
              IP1=I+1
                  DO 190 J=IP1,K
  190             SM=SM+A(I,J)*B(J,JB)
  200         continue
  210         B(I,JB)=(B(I,JB)-SM)/A(I,I)
C
C     COMPLETE COMPUTATION OF SOLUTION VECTOR.
C    ..
          IF (K.EQ.N) GO TO 240
              DO 220 J=KP1,N
  220         B(J,JB)=ZERO
              DO 230 I=1,K
  230         CALL H12 (2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,MDB,1)
C
C      RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
C      COLUMN INTERCHANGES.
C    ..
  240         DO 250 JJ=1,LDIAG
              J=LDIAG+1-JJ
              IF (IP(J).EQ.J) GO TO 250
              L=IP(J)
              TMP=B(L,JB)
              B(L,JB)=B(J,JB)
              B(J,JB)=TMP
  250         CONTINUE
  260     CONTINUE
C    ..
C     THE SOLUTION VECTORS, X, ARE NOW
C     IN THE FIRST  N  ROWS OF THE ARRAY B(,).
C
  270 KRANK=K
      RETURN
      END


      double precision FUNCTION DIFF(X,Y)
c
c  Function used in tests that depend on machine precision.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 7, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
      double precision X, Y
      DIFF=X-Y
      RETURN
      END

C     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE)
C
C  Algorithm NNLS: NONNEGATIVE LEAST SQUARES
C
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 15, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
c
C     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
C     N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM
C
C                      A * X = B  SUBJECT TO X .GE. 0
C     ------------------------------------------------------------------
c                     Subroutine Arguments
c
C     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE
C                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N
C                     MATRIX, A.           ON EXIT A() CONTAINS
C                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN
C                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY
C                     THIS SUBROUTINE.
C     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON-
C             TAINS Q*B.
C     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL
C             CONTAIN THE SOLUTION VECTOR.
C     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE
C             RESIDUAL VECTOR.
C     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN
C             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.
C             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z
C     ZZ()     AN M-ARRAY OF WORKING SPACE.
C     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
C                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS
C                 P AND Z AS FOLLOWS..
C
C                 INDEX(1)   THRU INDEX(NSETP) = SET P.
C                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.
C                 IZ1 = NSETP + 1 = NPP1
C                 IZ2 = N
C     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING
C             MEANINGS.
C             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
C             2     THE DIMENSIONS OF THE PROBLEM ARE BAD.
C                   EITHER M .LE. 0 OR N .LE. 0.
C             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS.
C
C     ------------------------------------------------------------------
      SUBROUTINE NNLS (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE)
C     ------------------------------------------------------------------
      integer I, II, IP, ITER, ITMAX, IZ, IZ1, IZ2, IZMAX, J, JJ, JZ, L
      integer M, MDA, MODE,N, NPP1, NSETP, RTNKEY
c     integer INDEX(N)
c     double precision A(MDA,N), B(M), W(N), X(N), ZZ(M)
      integer INDEX(*)
      double precision A(MDA,*), B(*), W(*), X(*), ZZ(*)
      double precision ALPHA, ASAVE, CC, FACTOR, RNORM
      double precision SM, SS, T, TEMP, TWO, UNORM, UP, WMAX
      double precision DUMMY(1)
      double precision ZERO, ZTEST
      parameter(FACTOR = 0.01d0)
      parameter(TWO = 2.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      MODE=1
      IF (M .le. 0 .or. N .le. 0) then
         MODE=2
         RETURN
      endif
      ITER=0
      ITMAX=3*N
C
C                    INITIALIZE THE ARRAYS INDEX() AND X().
C
          DO 20 I=1,N
          X(I)=ZERO
   20     INDEX(I)=I
C
      IZ2=N
      IZ1=1
      NSETP=0
      NPP1=1
C                             ******  MAIN LOOP BEGINS HERE  ******
   30 CONTINUE
C                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
C                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.
C
      IF (IZ1 .GT.IZ2.OR.NSETP.GE.M) GO TO 350
C
C         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
C
      DO 50 IZ=IZ1,IZ2
         J=INDEX(IZ)
         SM=ZERO
         DO 40 L=NPP1,M
   40        SM=SM+A(L,J)*B(L)
         W(J)=SM
   50 continue
C                                   FIND LARGEST POSITIVE W(J).
   60 continue
      WMAX=ZERO
      DO 70 IZ=IZ1,IZ2
         J=INDEX(IZ)
         IF (W(J) .gt. WMAX) then
            WMAX=W(J)
            IZMAX=IZ
         endif
   70 CONTINUE
C
C             IF WMAX .LE. 0. GO TO TERMINATION.
C             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
C
      IF (WMAX .le. ZERO) go to 350
      IZ=IZMAX
      J=INDEX(IZ)
C
C     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.
C     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID
C     NEAR LINEAR DEPENDENCE.
C
      ASAVE=A(NPP1,J)
      CALL H12 (1,NPP1,NPP1+1,M,A(1,J),1,UP,DUMMY,1,1,0)
      UNORM=ZERO
      IF (NSETP .ne. 0) then
          DO 90 L=1,NSETP
   90       UNORM=UNORM+A(L,J)**2
      endif
      UNORM=sqrt(UNORM)
      IF (DIFF(UNORM+ABS(A(NPP1,J))*FACTOR,UNORM) .gt. ZERO) then
C
C        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
C        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).
C
         DO 120 L=1,M
  120        ZZ(L)=B(L)
         CALL H12 (2,NPP1,NPP1+1,M,A(1,J),1,UP,ZZ,1,1,1)
         ZTEST=ZZ(NPP1)/A(NPP1,J)
C
C                                     SEE IF ZTEST IS POSITIVE
C
         IF (ZTEST .gt. ZERO) go to 140
      endif
C
C     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.
C     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
C     COEFFS AGAIN.
C
      A(NPP1,J)=ASAVE
      W(J)=ZERO
      GO TO 60
C
C     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
C     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER
C     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN
C     COL J,  SET W(J)=0.
C
  140 continue
      DO 150 L=1,M
  150    B(L)=ZZ(L)
C
      INDEX(IZ)=INDEX(IZ1)
      INDEX(IZ1)=J
      IZ1=IZ1+1
      NSETP=NPP1
      NPP1=NPP1+1
C
      IF (IZ1 .le. IZ2) then
         DO 160 JZ=IZ1,IZ2
            JJ=INDEX(JZ)
            CALL H12 (2,NSETP,NPP1,M,A(1,J),1,UP,A(1,JJ),1,MDA,1)
  160    continue
      endif
C
      IF (NSETP .ne. M) then
         DO 180 L=NPP1,M
  180       A(L,J)=ZERO
      endif
C
      W(J)=ZERO
C                                SOLVE THE TRIANGULAR SYSTEM.
C                                STORE THE SOLUTION TEMPORARILY IN ZZ().
      RTNKEY = 1
      GO TO 400
  200 CONTINUE
C
C                       ******  SECONDARY LOOP BEGINS HERE ******
C
C                          ITERATION COUNTER.
C
  210 continue
      ITER=ITER+1
      IF (ITER .gt. ITMAX) then
         MODE=3
         write (*,'(/a)') ' NNLS quitting on iteration count.'
         GO TO 350
      endif
C
C                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.
C                                  IF NOT COMPUTE ALPHA.
C
      ALPHA=TWO
      DO 240 IP=1,NSETP
         L=INDEX(IP)
         IF (ZZ(IP) .le. ZERO) then
            T=-X(L)/(ZZ(IP)-X(L))
            IF (ALPHA .gt. T) then
               ALPHA=T
               JJ=IP
            endif
         endif
  240 CONTINUE
C
C          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL
C          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.
C
      IF (ALPHA.EQ.TWO) GO TO 330
C
C          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO
C          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.
C
      DO 250 IP=1,NSETP
         L=INDEX(IP)
         X(L)=X(L)+ALPHA*(ZZ(IP)-X(L))
  250 continue
C
C        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I
C        FROM SET P TO SET Z.
C
      I=INDEX(JJ)
  260 continue
      X(I)=ZERO
C
      IF (JJ .ne. NSETP) then
         JJ=JJ+1
         DO 280 J=JJ,NSETP
            II=INDEX(J)
            INDEX(J-1)=II
            CALL G1 (A(J-1,II),A(J,II),CC,SS,A(J-1,II))
            A(J,II)=ZERO
            DO 270 L=1,N
               IF (L.NE.II) then
c
c                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))
c
                  TEMP = A(J-1,L)
                  A(J-1,L) = CC*TEMP + SS*A(J,L)
                  A(J,L)   =-SS*TEMP + CC*A(J,L)
               endif
  270       CONTINUE
c
c                 Apply procedure G2 (CC,SS,B(J-1),B(J))
c
            TEMP = B(J-1)
            B(J-1) = CC*TEMP + SS*B(J)
            B(J)   =-SS*TEMP + CC*B(J)
  280    continue
      endif
c
      NPP1=NSETP
      NSETP=NSETP-1
      IZ1=IZ1-1
      INDEX(IZ1)=I
C
C        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
C        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
C        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY
C        THAT ARE NONPOSITIVE WILL BE SET TO ZERO
C        AND MOVED FROM SET P TO SET Z.
C
      DO 300 JJ=1,NSETP
         I=INDEX(JJ)
         IF (X(I) .le. ZERO) go to 260
  300 CONTINUE
C
C         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
C
      DO 310 I=1,M
  310     ZZ(I)=B(I)
      RTNKEY = 2
      GO TO 400
  320 CONTINUE
      GO TO 210
C                      ******  END OF SECONDARY LOOP  ******
C
  330 continue
      DO 340 IP=1,NSETP
          I=INDEX(IP)
  340     X(I)=ZZ(IP)
C        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.
      GO TO 30
C
C                        ******  END OF MAIN LOOP  ******
C
C                        COME TO HERE FOR TERMINATION.
C                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.
C
  350 continue
      SM=ZERO
      IF (NPP1 .le. M) then
         DO 360 I=NPP1,M
  360       SM=SM+B(I)**2
      else
         DO 380 J=1,N
  380       W(J)=ZERO
      endif
      RNORM=sqrt(SM)
      RETURN
C
C     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE
C     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().
C
  400 continue
      DO 430 L=1,NSETP
         IP=NSETP+1-L
         IF (L .ne. 1) then
            DO 410 II=1,IP
               ZZ(II)=ZZ(II)-A(II,JJ)*ZZ(IP+1)
  410       continue
         endif
         JJ=INDEX(IP)
         ZZ(IP)=ZZ(IP)/A(IP,JJ)
  430 continue
      go to (200, 320), RTNKEY
      END

      SUBROUTINE G1 (A,B,CTERM,STERM,SIG)
c
C     COMPUTE ORTHOGONAL ROTATION MATRIX..
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
C     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))
C                        (-S,C)         (-S,C)(B)   (   0          )
C     COMPUTE SIG = SQRT(A**2+B**2)
C        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT
C        SIG MAY BE IN THE SAME LOCATION AS A OR B .
C     ------------------------------------------------------------------
      double precision A, B, CTERM, ONE, SIG, STERM, XR, YR, ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      if (abs(A) .gt. abs(B)) then
         XR=B/A
         YR=sqrt(ONE+XR**2)
         CTERM=sign(ONE/YR,A)
         STERM=CTERM*XR
         SIG=abs(A)*YR
         RETURN
      endif

      if (B .ne. ZERO) then
         XR=A/B
         YR=sqrt(ONE+XR**2)
         STERM=sign(ONE/YR,B)
         CTERM=STERM*XR
         SIG=abs(B)*YR
         RETURN
      endif

      SIG=ZERO
      CTERM=ZERO
      STERM=ONE
      RETURN
      END

      SUBROUTINE LDP (G,MDG,M,N,H,X,XNORM,W,INDEX,MODE)
c
C  Algorithm LDP: LEAST DISTANCE PROGRAMMING
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1974 MAR 1, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
      integer I, IW, IWDUAL, IY, IZ, J, JF, M, MDG, MODE, N, NP1
c     integer INDEX(M)
c     double precision G(MDG,N), H(M), X(N), W(*)
      integer INDEX(*)
      double precision G(MDG,*), H(*), X(*), W(*)
      double precision FAC, ONE, RNORM, XNORM, ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      IF (N.LE.0) GO TO 120
          DO 10 J=1,N
   10     X(J)=ZERO
      XNORM=ZERO
      IF (M.LE.0) GO TO 110
C
C     THE DECLARED DIMENSION OF W() MUST BE AT LEAST (N+1)*(M+2)+2*M.
C
C      FIRST (N+1)*M LOCS OF W()   =  MATRIX E FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR F FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR Z FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR Y FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR WDUAL FOR PROBLEM NNLS.
C     COPY G**T INTO FIRST N ROWS AND M COLUMNS OF E.
C     COPY H**T INTO ROW N+1 OF E.
C
      IW=0
          DO 30 J=1,M
              DO 20 I=1,N
              IW=IW+1
   20         W(IW)=G(J,I)
          IW=IW+1
   30     W(IW)=H(J)
      JF=IW+1
C                                STORE N ZEROS FOLLOWED BY A ONE INTO F.
          DO 40 I=1,N
          IW=IW+1
   40     W(IW)=ZERO
      W(IW+1)=ONE
C
      NP1=N+1
      IZ=IW+2
      IY=IZ+NP1
      IWDUAL=IY+M
C
      CALL NNLS (W,NP1,NP1,M,W(JF),W(IY),RNORM,W(IWDUAL),W(IZ),INDEX,
     *           MODE)
C                      USE THE FOLLOWING RETURN IF UNSUCCESSFUL IN NNLS.
      IF (MODE.NE.1) RETURN
      IF (RNORM) 130,130,50
   50 FAC=ONE
      IW=IY-1
          DO 60 I=1,M
          IW=IW+1
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
   60     FAC=FAC-H(I)*W(IW)
C
      IF (DIFF(ONE+FAC,ONE)) 130,130,70
   70 FAC=ONE/FAC
          DO 90 J=1,N
          IW=IY-1
              DO 80 I=1,M
              IW=IW+1
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
   80         X(J)=X(J)+G(I,J)*W(IW)
   90     X(J)=X(J)*FAC
          DO 100 J=1,N
  100     XNORM=XNORM+X(J)**2
      XNORM=sqrt(XNORM)
C                             SUCCESSFUL RETURN.
  110 MODE=1
      RETURN
C                             ERROR RETURN.       N .LE. 0.
  120 MODE=2
      RETURN
C                             RETURNING WITH CONSTRAINTS NOT COMPATIBLE.
  130 MODE=4
      RETURN
      END

      SUBROUTINE BNDACC (G,MDG,NB,IP,IR,MT,JT)
c
C  SEQUENTIAL ALGORITHM FOR BANDED LEAST SQUARES PROBLEM..
C  ACCUMULATION PHASE.      FOR SOLUTION PHASE USE BNDSOL.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
C     THE CALLING PROGRAM MUST SET IR=1 AND IP=1 BEFORE THE FIRST CALL
C     TO BNDACC FOR A NEW CASE.
C
C     THE SECOND SUBSCRIPT OF G( ) MUST BE DIMENSIONED AT LEAST
C     NB+1 IN THE CALLING PROGRAM.
c     ------------------------------------------------------------------
      integer I, J, IE, IG, IG1, IG2, IP, IR, JG, JT, K, KH, L, LP1
      integer MDG, MH, MT, MU, NB, NBP1
c     double precision G(MDG,NB+1)
      double precision G(MDG,*)
      double precision RHO, ZERO
      parameter(ZERO = 0.0d0)
c     ------------------------------------------------------------------
C
C              ALG. STEPS 1-4 ARE PERFORMED EXTERNAL TO THIS SUBROUTINE.
C
      NBP1=NB+1
      IF (MT.LE.0) RETURN
C                                             ALG. STEP 5
      IF (JT.EQ.IP) GO TO 70
C                                             ALG. STEPS 6-7
      IF (JT.LE.IR) GO TO 30
C                                             ALG. STEPS 8-9
      DO 10 I=1,MT
        IG1=JT+MT-I
        IG2=IR+MT-I
        DO 10 J=1,NBP1
   10   G(IG1,J)=G(IG2,J)
C                                             ALG. STEP 10
      IE=JT-IR
      DO 20 I=1,IE
        IG=IR+I-1
        DO 20 J=1,NBP1
   20   G(IG,J)=ZERO
C                                             ALG. STEP 11
      IR=JT
C                                             ALG. STEP 12
   30 MU=min(NB-1,IR-IP-1)
      IF (MU.EQ.0) GO TO 60
C                                             ALG. STEP 13
      DO 50 L=1,MU
C                                             ALG. STEP 14
        K=min(L,JT-IP)
C                                             ALG. STEP 15
        LP1=L+1
        IG=IP+L
        DO 40 I=LP1,NB
          JG=I-K
   40     G(IG,JG)=G(IG,I)
C                                             ALG. STEP 16
        DO 50 I=1,K
        JG=NBP1-I
   50   G(IG,JG)=ZERO
C                                             ALG. STEP 17
   60 IP=JT
C                                             ALG. STEPS 18-19
   70 MH=IR+MT-IP
      KH=min(NBP1,MH)
C                                             ALG. STEP 20
      DO 80 I=1,KH
        CALL H12 (1,I,max(I+1,IR-IP+1),MH,G(IP,I),1,RHO,
     *            G(IP,I+1),1,MDG,NBP1-I)
   80 continue
C                                             ALG. STEP 21
      IR=IP+KH
C                                             ALG. STEP 22
      IF (KH.LT.NBP1) GO TO 100
C                                             ALG. STEP 23
      DO 90 I=1,NB
   90   G(IR-1,I)=ZERO
C                                             ALG. STEP 24
  100 CONTINUE
C                                             ALG. STEP 25
      RETURN
      END


      SUBROUTINE BNDSOL (MODE,G,MDG,NB,IP,IR,X,N,RNORM)
c
C  SEQUENTIAL SOLUTION OF A BANDED LEAST SQUARES PROBLEM..
C  SOLUTION PHASE.   FOR THE ACCUMULATION PHASE USE BNDACC.
c
c  The original version of this code was developed by
c  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
c  1973 JUN 12, and published in the book
c  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
c  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
C     MODE = 1     SOLVE R*X=Y   WHERE R AND Y ARE IN THE G( ) ARRAY
C                  AND X WILL BE STORED IN THE X( ) ARRAY.
C            2     SOLVE (R**T)*X=Y   WHERE R IS IN G( ),
C                  Y IS INITIALLY IN X( ), AND X REPLACES Y IN X( ),
C            3     SOLVE R*X=Y   WHERE R IS IN G( ).
C                  Y IS INITIALLY IN X( ), AND X REPLACES Y IN X( ).
C
C     THE SECOND SUBSCRIPT OF G( ) MUST BE DIMENSIONED AT LEAST
C     NB+1 IN THE CALLING PROGRAM.
      integer I, I1, I2, IE, II, IP, IR, IX, J, JG, L
      integer MDG, MODE, N, NB, NP1, IRM1
      double precision G(MDG,*), RNORM, RSQ, S, X(N), ZERO
      parameter(ZERO = 0.0d0)
C
      RNORM=ZERO
      GO TO (10,90,50), MODE
C                                   ********************* MODE = 1
C                                   ALG. STEP 26
   10      DO 20 J=1,N
   20      X(J)=G(J,NB+1)
      RSQ=ZERO
      NP1=N+1
      IRM1=IR-1
      IF (NP1.GT.IRM1) GO TO 40
           DO 30 J=NP1,IRM1
   30      RSQ=RSQ+G(J,NB+1)**2
      RNORM=SQRT(RSQ)
   40 CONTINUE
C                                   ********************* MODE = 3
C                                   ALG. STEP 27
   50      DO 80 II=1,N
           I=N+1-II
C                                   ALG. STEP 28
           S=ZERO
           L=max(0,I-IP)
C                                   ALG. STEP 29
           IF (I.EQ.N) GO TO 70
C                                   ALG. STEP 30
           IE=min(N+1-I,NB)
                DO 60 J=2,IE
                JG=J+L
                IX=I-1+J
   60           S=S+G(I,JG)*X(IX)
C                                   ALG. STEP 31
   70      continue
           IF (G(I,L+1) .eq. ZERO) go to 130
   80      X(I)=(X(I)-S)/G(I,L+1)
C                                   ALG. STEP 32
      RETURN
C                                   ********************* MODE = 2
   90      DO 120 J=1,N
           S=ZERO
           IF (J.EQ.1) GO TO 110
           I1=max(1,J-NB+1)
           I2=J-1
                DO 100 I=I1,I2
                L=J-I+1+max(0,I-IP)
  100           S=S+X(I)*G(I,L)
  110      L=max(0,J-IP)
           IF (G(J,L+1) .eq. ZERO) go to 130
  120      X(J)=(X(J)-S)/G(J,L+1)
      RETURN
C
  130 write (*,'(/a/a,4i6)')' ZERO DIAGONAL TERM IN BNDSOL.',
     *      ' MODE,I,J,L = ',MODE,I,J,L
      STOP
      END
