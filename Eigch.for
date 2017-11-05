C-- EIGCH - ÑàAÉ KOMèã MATP -------------------------------------------
         SUBROUTINE EIGCH (A,N,JOBN,D,Z,IZ,WK,IER)
         INTEGER  N,JOBN,IZ,IER
         REAL     A(*),D(N),Z(*),WK(*)
         INTEGER  JER,K,I,NE,NTAU,NA,NI,NI2,IM1,J,IIZ,NZ,IIZ1,
     1            IJOB,JR,IR,IJ,JI,NP1,
     2            JZ,JZI,L,M,II,IL,KK,LZ,MZ,LK,KZ
         REAL     ANORM,ASUM,PI,SUMZ,SUMR,SUMI,S,TEN,RDELP,
     1            ZERO,ONE,THOUS,AN,SIGNA
         DATA     RDELP/0.9537E-06/
         DATA     ZERO,ONE/0.0,1.0/,TEN/10.0/,THOUS/1000.0/
         IER = 0
         JER = 0
         IF (JOBN.LT.10) GOTO 15
         JR = N + N - 2
         IJ = 2
         K  = 2
         DO 10 J=1,N
            DO 5 I=1,J
               A(K-1) = A(IJ-1)
               A(K) = -A(IJ)
               K = K+2
               IJ = IJ + 2
 5          CONTINUE
            IJ = IJ + JR
            JR = JR - 2
 10      CONTINUE
 15      IJOB = MOD(JOBN,10)
         IF (IJOB.GE.0.AND.IJOB.LE.3) GOTO 20
         IER = 66
         IJOB = 1
         GOTO 25
 20      IF (IJOB.EQ.0) GOTO 45
 25      IF (IZ.GE.N) GOTO 30
         IER = 67
         IJOB = 0
 30      K = 2
         DO 40 I=1,N
            IF (A(K).EQ.ZERO) GOTO 35
            A(K) = ZERO
            IER = 68
 35         K = K+I+I+2
 40      CONTINUE
         IF (IJOB.EQ.3) GOTO 110
 45      NE = 1
         NTAU = NE+N
         NA = NTAU+N+N
         NI = (N*(N+1))/2
         NI2 = NI+NI
         IF (IJOB.NE.2) GOTO 55
         K = NA
         DO 50 I=1,NI2
            WK(K) = A(I)
            K = K+1
 50      CONTINUE
 55      IF (NI.LT.2) GOTO 70
         IM1 = 1
         DO 65 I=2,NI
            K = IM1+I
            PI = A(K)
            DO 60 J=1,IM1
               A(K) = A(K-1)
               K = K-1
 60         CONTINUE
            A(I) = PI
            IM1 = I
 65      CONTINUE
 70      CALL EHOUSH (A(1),A(NI+1),N,D,WK(NE),WK(NTAU))
         IIZ = 1
         IF (IJOB.NE.0) IIZ = IZ+IZ
         IF (IIZ.EQ.1) GOTO 85
         NZ = (IZ+IZ)*N
         DO 75 I=1,NZ
            Z(I) = ZERO
 75      CONTINUE
         K = 1
         IIZ1 = IIZ+1
         DO 80 I=1,N
            Z(K) = ONE
            K = K+IIZ1
 80      CONTINUE
 85      CALL EQRT2S (D,WK(NE),N,Z(1),IIZ,JER)
         IF (IJOB.EQ.0) GOTO 9000
         CALL EHBCKH (A(1),A(NI+1),WK(NTAU),N,Z(1),Z(IZ+1),IIZ)
         JZ = 0
         DO 100 J=1,N
            JZI = JZ+IZ
            DO 90 I=1,N
               K = JZI+I
               WK(I) = Z(K)
 90         CONTINUE
            K = JZ+N
            L = K+N-1
            M = N
            DO 95 I=1,N
               Z(L) = Z(K)
               Z(L+1) = WK(M)
               K = K-1
               L = L-2
               M = M-1
 95         CONTINUE
            JZ = JZ+IZ+IZ
 100     CONTINUE
         IF (IJOB.NE.2) GOTO 9000
         K = NA
         DO 105 I=1,NI2
            A(I) = WK(K)
            K = K+1
 105     CONTINUE
         WK(1) = THOUS
         IF (JER.NE.0) GOTO 9000
 110     ANORM = ZERO
         II = 1
         DO 120 I=1,N
            ASUM = ZERO
            IL = II
            KK = 2
            DO 115 L=1,N
               ASUM = ASUM+CABS(CMPLX(A(IL),A(IL+1)))
               IF (L.GE.I) KK = L+L
               IL = IL+KK
 115        CONTINUE
            ANORM = AMAX1(ANORM,ASUM)
            II = II+I+I
 120     CONTINUE
         IF (ANORM.EQ.ZERO) ANORM = ONE
         PI = ZERO
         DO 135 I=1,N
            II = 1
            S = ZERO
            SUMZ = ZERO
            LZ = (IZ+IZ)*(I-1)+1
            LZ = IZ*(I-1)*2+1
            MZ = LZ
            DO 130 L=1,N
               LK = II
               KK = 2
               SUMZ = SUMZ+CABS(CMPLX(Z(LZ),Z(LZ+1)))
               SUMR = -D(I)*Z(LZ)
               SUMI = -D(I)*Z(LZ+1)
               KZ = MZ
               DO 125 K=1,N
                  SIGNA = ONE
                  IF (K.GT.L) SIGNA = -ONE
                  SUMR = SUMR+A(LK)*Z(KZ)-SIGNA*A(LK+1)*Z(KZ+1)
                  SUMI = SUMI+A(LK)*Z(KZ+1)+SIGNA*A(LK+1)*Z(KZ)
                  IF (K.GE.L) KK = K+K
                  LK = LK+KK
                  KZ = KZ+2
 125           CONTINUE
               S = S+CABS(CMPLX(SUMR,SUMI))
               LZ = LZ+2
               II = II+L+L
 130        CONTINUE
            IF (SUMZ.EQ.ZERO) GOTO 135
            PI = AMAX1(PI,S/SUMZ)
 135     CONTINUE
         AN = N
         PI = PI/(ANORM*TEN*AN*RDELP)
         WK(1) = PI
         IF (JOBN.LT.10) GOTO 9000
         NP1 = N+1
         IJ = (N-1) * NP1
         IJ = IJ + IJ + 2
         K = N * NP1
         DO 145 JR=1,N
            J = N+1-JR
            DO 140 IR=1,J
               A(IJ-1) = A(K-1)
               A(IJ) = -A(K)
               K = K-2
               IJ = IJ - 2
 140        CONTINUE
            IJ = IJ -JR - JR
 145     CONTINUE
         JR = N + N
         II = 2
         JI = 2
         DO 155 I=1,N
            IJ = II
            DO 150 J=1,I
               A(IJ-1) = A(JI-1)
               A(IJ) = -A(JI)
               JI = JI+2
               IJ = IJ+JR
 150        CONTINUE
            JI = JI + JR - I - I
            II = II + 2
 155     CONTINUE
 9000    CONTINUE
         IF (IER.NE.0) CALL UERTST (IER,6HEIGCH )
         IF (JER.EQ.0) GOTO 9005
         IER = JER
         CALL UERTST (IER,6HEIGCH )
 9005    RETURN
         END
C----------------------------------------------------------------------
         SUBROUTINE EHBCKH (AR,AI,TAU,N,ZR,ZI,IZ)
         INTEGER  N,IZ
         REAL     AR(1),AI(1),TAU(2,1),ZR(IZ,1),ZI(IZ,1)
         INTEGER  J,K,NR,L,NRM1,INX1,INX2,K1
         REAL     DELTA,ZERO,ALPHA1,ALPHA2
         DATA     ZERO/0.0/
         DO 5 J=1,N
            DO 5 K=1,N
               ZI(J,K)=-ZR(J,K)*TAU(2,J)
               ZR(J,K)= ZR(J,K)*TAU(1,J)
 5       CONTINUE
         IF (N.LE.2) GOTO 30
         DO 25 L=3,N
            NR = N-L+2
            NRM1 = NR-1
            INX1 = (NR*(NRM1))/2+NR
            INX2 = INX1-1
            IF (AI(INX1) .EQ. ZERO) GOTO 25
            DELTA = AI(INX1) * SQRT(AR(INX2)**2+AI(INX2)**2)
            DO 20 J=1,N
               ALPHA1 = ZERO
               ALPHA2 = ZERO
               DO 10 K=NR,N
                  K1 = (K*(K-1))/2+NRM1
                  ALPHA1=ALPHA1+AR(K1)*ZR(K,J)+AI(K1)*ZI(K,J)
                  ALPHA2=ALPHA2-AI(K1)*ZR(K,J)+AR(K1)*ZI(K,J)
 10            CONTINUE
               ALPHA1 = ALPHA1/DELTA
               ALPHA2 = ALPHA2/DELTA
               DO 15 K=NR,N
                  K1 = (K*(K-1))/2+NRM1
                  ZR(K,J)=ZR(K,J)-AR(K1)*ALPHA1+AI(K1)*ALPHA2
                  ZI(K,J)=ZI(K,J)-AR(K1)*ALPHA2-AI(K1)*ALPHA1
 15            CONTINUE
 20         CONTINUE
 25      CONTINUE
 30      RETURN
         END
C----------------------------------------------------------------------
         SUBROUTINE EHOUSH (AR,AI,N,D,E,TAU)
         INTEGER  N
         REAL     AR(1),AI(1),D(1),E(1),TAU(2,1)
         INTEGER  NM1,NN,I,NR,NRM1,L,INDX,J,JJ,INX1,INX2,JP1,KK,IX,IM1
         REAL     RHO,TOLER,ZERO,ONE,T1,T2,TESTBB,VR,ROOT,DELTA,
     *            RATIC,RDELP,Q1,Q2,X1,X2,TT1,TT2,BB
         DATA     ZERO/0.0/,ONE/1.0/
         DATA     RDELP/0.9537E-06/
         NM1 = N-1
         TOLER=ZERO
         NN=(N*(N+1))/2
         DO 5 I=1,NN
            T1=ABS(AR(I))
            T2=ABS(AI(I))
            IF(T2.GT.T1) T1=T2
            IF(T1.GT.TOLER) TOLER=T1
 5       CONTINUE
         TESTBB=RDELP*TOLER
         IF(N.LE.2) GOTO 65
         DO 60 NR=2,NM1
            NRM1=NR-1
            VR=ZERO
            TAU(1,NR)=ZERO
            TAU(2,NR)=ZERO
            TAU(2,1)=ZERO
            DO 10 L=NR,N
               INDX=(L*(L-1))/2+NRM1
               VR=AR(INDX)**2+AI(INDX)**2+VR
 10         CONTINUE
            INDX=(NR*NRM1)/2+NRM1
            IF((TESTBB)**2 .GE. VR) GOTO 60
            ROOT=CABS(CMPLX(AR(INDX),AI(INDX)))*SQRT(VR)
            IF(ROOT.NE.ZERO) GOTO 15
            AR(INDX)=SQRT(VR)
            DELTA=VR
            TAU(1,1)=-AR(INDX)
            GOTO 20
 15         DELTA=VR+ROOT
            RATIO=VR/ROOT
            TAU(1,1)=-RATIO*AR(INDX)
            TAU(2,1)= RATIO*AI(INDX)
            AR(INDX)=(RATIO+ONE)*AR(INDX)
            AI(INDX)=(RATIO+ONE)*AI(INDX)
 20         DO 35 J=NR,N
               JJ=(J*(J-1))/2
               INDX=JJ+NRM1
               TAU(1,J)=AR(INDX)/DELTA
               TAU(2,J)=AI(INDX)/DELTA
               D(J)=ZERO
               E(J)=ZERO
               DO 25 L=NR,J
                  INX1=(L*(L-1))/2+NRM1
                  INX2=JJ+L
                  D(J)=D(J)+AR(INX2)*AR(INX1)-AI(INX2)*AI(INX1)
                  E(J)=E(J)+AR(INX2)*AI(INX1)+AI(INX2)*AR(INX1)
 25            CONTINUE
               JP1=J+1
               IF(JP1.GT.N) GOTO 40
               DO 30 L=JP1,N
                  KK=(L*(L-1))/2
                  INX1=KK+NRM1
                  INX2=KK+J
                  D(J)=D(J)+AR(INX2)*AR(INX1)+AI(INX2)*AI(INX1)
                  E(J)=E(J)+AR(INX2)*AI(INX1)-AI(INX2)*AR(INX1)
 30            CONTINUE
 35         CONTINUE
 40         RHO=ZERO
            DO 45 L=NR,N
               RHO=RHO+D(L)*TAU(1,L)+E(L)*TAU(2,L)
 45         CONTINUE
            IX=(NRM1*(NR-2))/2
            DO 55 I=NR,N
               IX=IX+I-1
               INX2=IX+NRM1
               DO 50 J=NR,I
                  INX1=IX+J
                  X1=TAU(1,I)*D(J)+TAU(2,I)*E(J)
                  X2=TAU(2,I)*D(J)-TAU(1,I)*E(J)
                  Q1=D(I)-RHO*AR(INX2)
                  Q2=E(I)-RHO*AI(INX2)
                  T1=Q1*TAU(1,J)+Q2*TAU(2,J)
                  T2=Q2*TAU(1,J)-Q1*TAU(2,J)
                  AR(INX1)=AR(INX1)-X1-T1
                  AI(INX1)=AI(INX1)-X2-T2
 50            CONTINUE
 55         CONTINUE
            TAU(1,NR)=TAU(1,1)
            TAU(2,NR)=TAU(2,1)
 60      CONTINUE
 65      INDX=0
         DO 70 I=1,N
            INDX=INDX+I
            D(I)=AR(INDX)
 70      CONTINUE
         TAU(1,1)=ONE
         TAU(2,1)=ZERO
         E(1)=ZERO
         IF(N.EQ.1) GOTO 85
         INDX=(N*NM1)/2+NM1
         TAU(1,N)= AR(INDX)
         TAU(2,N)=-AI(INDX)
         INDX=1
         DO 80 I=2,N
            INDX=INDX+I
            IM1=I-1
            BB= SQRT(TAU(1,I)**2+TAU(2,I)**2)
            E(I)=BB
            AI(INDX)=BB
            IF(TESTBB.LT.BB) GOTO 75
            TAU(1,I)=ONE
            TAU(2,I)=ZERO
            BB=ONE
 75         TT1=TAU(1,I)*TAU(1,IM1)-TAU(2,I)*TAU(2,IM1)
            TT2=TAU(1,I)*TAU(2,IM1)+TAU(2,I)*TAU(1,IM1)
            TAU(1,I)=TT1/BB
            TAU(2,I)=TT2/BB
 80      CONTINUE
 85      RETURN
         END
C----------------------------------------------------------------------
         SUBROUTINE EQRT2S (D,E,N,Z,IZ,IER)
         DIMENSION  D(1),E(1),Z(IZ,1)
         DATA       RDELP/0.9537E-06/
         DATA       ZERO,ONE/0.0,1.0/
         IER = 0
         IF (N .EQ. 1) GOTO 9005
         DO 5 I=2,N
            E(I-1) = E(I)
 5       CONTINUE
         E(N) = ZERO
         B = ZERO
         F = ZERO
         DO 60 L=1,N
            J = 0
            H = RDELP*(ABS(D(L))+ABS(E(L)))
            IF (B.LT.H)  B = H
            DO 10 M=L,N
               K = M
               IF (ABS(E(K)) .LE. B) GOTO 15
 10         CONTINUE
 15         M = K
            IF (M .EQ. L ) GOTO 55
 20         IF (J .EQ. 30) GOTO 85
            J  = J+1
            L1 = L+1
            G  = D(L)
            P  = (D(L1)-G)/(E(L)+E(L))
            R  = SQRT(P*P+ONE)
            D(L) = E(L)/(P+SIGN(R,P))
            H  = G - D(L)
            DO 25 I = L1,N
               D(I) = D(I) - H
 25         CONTINUE
            F = F+H
            P = D(M)
            C = ONE
            S = ZERO
            MM1 = M-1
            MM1PL = MM1+L
            IF (L.GT.MM1) GOTO 50
            DO 45 II=L,MM1
               I = MM1PL-II
               G = C*E(I)
               H = C*P
               IF (ABS(P).LT.ABS(E(I))) GOTO 30
               C = E(I)/P
               R = SQRT(C*C+ONE)
               E(I+1) = S*P*R
               S = C/R
               C = ONE/R
               GOTO 35
 30            C = P/E(I)
               R = SQRT(C*C+ONE)
               E(I+1) = S*E(I)*R
               S = ONE/R
               C = C*S
 35            P = C*D(I)-S*G
               D(I+1) = H+S*(C*G+S*D(I))
               IF (IZ.LT.N) GOTO 45
               DO 40 K=1,N
                  H = Z(K,I+1)
                  Z(K,I+1) = S*Z(K,I)+C*H
                  Z(K,I) = C*Z(K,I)-S*H
 40            CONTINUE
 45         CONTINUE
 50         E(L) = S*P
            D(L) = C*P
            IF (ABS(E(L)) .GT. B) GOTO 20
 55         D(L) = D(L)+F
 60      CONTINUE
         DO 80 I=1,N
            K = I
            P = D(I)
            IP1 = I+1
            IF (IP1.GT.N) GOTO 70
            DO 65 J=IP1,N
               IF (D(J).GE.P) GOTO 65
               K = J
               P = D(J)
 65         CONTINUE
 70         IF (K.EQ.I) GOTO 80
            D(K) = D(I)
            D(I) = P
            IF (IZ.LT.N) GOTO 80
            DO 75 J=1,N
               P = Z(J,I)
               Z(J,I) = Z(J,K)
               Z(J,K) = P
 75         CONTINUE
 80      CONTINUE
         GOTO 9005
 85      IER = 128+L
 9000    CONTINUE
         CALL UERTST(IER,6HEQRT2S)
 9005    RETURN
         END
C----------------------------------------------------------------------
         SUBROUTINE UERTST (IER,NAME)
         INTEGER    IER,NAME(3)
         INTEGER    NAMSET(3),NAMEQ(3)
         DATA       NAMSET/2HUE,2HRS,2HET/
         DATA       NAMEQ/2H  ,2H  ,2H  /
         DATA       LEVEL/4/,IEQDF/0/,IEQ/1H=/
         IF (IER.GT.999) GOTO 25
         IF (IER.LT.-32) GOTO 55
         IF (IER.LE.128) GOTO 5
         IF (LEVEL.LT.1) GOTO 30
C        CALL UGETIO(1,NIN,IOUNIT)
         IF (IEQDF.EQ.1) PRINT 35,IER,NAMEQ,IEQ,NAME
         IF (IEQDF.EQ.0) PRINT 35,IER,NAME
         GOTO 30
 5       IF (IER.LE.64) GOTO 10
         IF (LEVEL.LT.2) GOTO 30
C        CALL UGETIO(1,NIN,IOUNIT)
         IF (IEQDF.EQ.1) PRINT 40, IER,NAMEQ,IEQ,NAME
         IF (IEQDF.EQ.0) PRINT 40, IER,NAME
         GOTO 30
 10      IF (IER.LE.32) GOTO 15
         IF (LEVEL.LT.3) GOTO 30
C        CALL UGETIO(1,NIN,IOUNIT)
         IF (IEQDF.EQ.1) PRINT 45, IER,NAMEQ,IEQ,NAME
         IF (IEQDF.EQ.0) PRINT 45, IER,NAME
         GOTO 30
 15      CONTINUE
         DO 20 I=1,3
            IF (NAME(I).NE.NAMSET(I)) GOTO 25
 20      CONTINUE
         LEVOLD = LEVEL
         LEVEL = IER
         IER = LEVOLD
         IF (LEVEL.LT.0) LEVEL = 4
         IF (LEVEL.GT.4) LEVEL = 4
         GOTO 30
 25      CONTINUE
         IF (LEVEL.LT.4) GOTO 30
C        CALL UGETIO(1,NIN,IOUNIT)
         IF (IEQDF.EQ.1) PRINT 50, IER,NAMEQ,IEQ,NAME
         IF (IEQDF.EQ.0) PRINT 50, IER,NAME
 30      IEQDF = 0
         RETURN
 35      FORMAT(19H *** TERMINAL ERROR,10X,7H(IER = ,I3,
     *          20H) FROM IMSL ROUTINE ,3A2,A1,3A2)
 40      FORMAT(36H *** WARNING WITH FIX ERROR  (IER = ,I3,
     *          20H) FROM IMSL ROUTINE ,3A2,A1,3A2)
 45      FORMAT(18H *** WARNING ERROR,11X,7H(IER = ,I3,
     *          20H) FROM IMSL ROUTINE ,3A2,A1,3A2)
 50      FORMAT(20H *** UNDEFINED ERROR,9X,7H(IER = ,I5,
     *          20H) FROM IMSL ROUTINE ,3A2,A1,3A2)
C
 55      IEQDF = 1
         DO 60 I=1,3
 60         NAMEQ(I) = NAME(I)
 65      RETURN
         END
C----------------------------------------------------------------------
         SUBROUTINE UGETIO (IOPT,NIN,NOUT)
         INTEGER            IOPT,NIN,NOUT
         INTEGER            NIND,NOUTD
         DATA               NIND/9/,NOUTD/12/
         IF (IOPT.EQ.3) GOTO 10
         IF (IOPT.EQ.2) GOTO 5
         IF (IOPT.NE.1) GOTO 9005
         NIN = NIND
         NOUT = NOUTD
         GOTO 9005
 5       NIND =NIN
         GOTO 9005
 10      NOUTD = NOUT
 9005    RETURN
         END
