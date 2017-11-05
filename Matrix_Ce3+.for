CCC ++++++ Ce ++++++++++++++++++++++++++++++++++++++++++++++++26.03.10++
       REAL*8     XX
       DIMENSION  XX(9),X(9),Q(51)

       REAL*4 MD21(21),ED21(6),BD21(6,6),XM,ZM,SM,XiX,XiY,XiZ,T,GJ
       REAL*4 XiX_Scaled, XiZ_Scaled
       REAL*4 MK21(42),EK(6),ZK(338),CR(6,6),CJ(6,6)
       REAL*4 C,Cel,Clat,Td,delta_Ground_1st,delta_Ground_2nd, EKM
       REAL*4 FI, TETA, AL, H, H2, AL3
       REAL*4 WWW, KondoTemp, cutoffPar, arg_SF
       INTEGER mGr,MExc
       DATA        PI /3.1415926/       


       OPEN(5, FILE='Input.dat', STATUS = 'OLD')
       OPEN(6, FILE='Output.dat', STATUS = 'OLD')
       OPEN(7, FILE='Energy.dat', STATUS = 'OLD')


        GJ = 6./7.
        AL = -1.999

C--------Field orientation-----------------------------------------------
		FI = 0.
		TETA = 0.

C--------FLEX------------------------------------------------------------
		NX = 9
		
C-----Debye temperature-------------------------------------------------
		Td = 250.

C-----Kondo temperature-------------------------------------------------
		KondoTemp = 8.5

C-----Cut-off parameter A-------------------------------------------------
		cutoffPar = 3.

C-----Degeneracy of ground mGr and first excited MExc states-------------------------------------------------
		mGr = 2
        MExc = 4


		
C-----------------------------------------------------------------------		
		T = 2.
		T2 = 320.
		
		H = 0.00
		H2 = 0.
C------------------AL1-коэфициент для расчета магнитной восприимчивости
C--------ед. измерения коэф. (см^3)*К (в подпрогр. воспр. без AL1 измер. в 1/K)
				AL1= 0.4578E-24
C--------для расчета маг. воспр. в mol/emu надо AL1 умножить на число Авогадро
				AL1=AL1*6.0248E23

C--------ед. измерения коэф. emu/G (воспр. из намагниченности)
	AL3 = 9.274E-24
C--------для расчета маг. воспр. в emu/mol надо AL3 умножить на число Авогадро
	AL3=AL3*6.0248E23
				
				
            READ(5,1000) (XX(J),J=1,NX)


 4       DO 5 J5=1,NX
 5          X(J5)= XX(J5)
 
 
         WRITE(6,86) AL, GJ, FI, TETA
         WRITE(6,6) ( X(J), J=1,NX )

 9     FORMAT (T3,'NT$B (TM)  BTC  12.06.2014    AL=',F9.4,'  GJ=',F9.4)

 6       FORMAT(/T3,'X(J)',10F12.5/)
 
 
      DO 73 WHILE (T.LT.T2)
	  
	     CALL XXXCE1( X,XX)
	     CALL FMD ( X,MD21)
         
	     CALL ZEEMAN(AL,H,MD21,FI,TETA)	     
	     CALL NDIAD (MD21,ED21,BD21,6,0)

         CALL NPZVD (H,ED21,BD21,6)

         CALL FM1(T,ED21,BD21, GJ, XDM,ZDM )
         CALL FM(FI,TETA,XDM,ZDM,SM)

c         Xizz = AL3 * SM / H
         
c         CALL HC (Td,T,ED21,C,Cel,Clat )
         CALL FMS( AL1,T,ED21,BD21,6,XiX,XiZ )

c         CALL FK6( AL,FI,TETA,H,MD21,MK21 )
cC           CALL PMK9 (MK9,9)
C....................  ДИAГ KOMПЛ MATP (1-COБ ZN И COБ V)......
c         CALL NDIAK (MK21,EK,ZK,6, 1)
c         CALL DMSV6 ( ZK,CR,CJ,6)
c         CALL NPSZVK (H,EK,ZK,6)
c         CALL FM6 (T,EK,CR,CJ, GJ, XM,YM,ZM )
c         CALL FMK ( FI,TETA,XM,YM,ZM,SM )
c         CALL FKMS ( AL1,T,EK,CR,CJ,6,XiX,XiY,XiZ )
c        CALL FMXI ( FI,TETA,XiX,XiY,XiZ,Xi )
cc         CALL NEKM (T,EK,EKM ) 

c   Gaps in meV
          
         delta_Ground_1st=ABS(ED21(6)-ED21(4)) * 0.124
         delta_Ground_2nd=ABS(ED21(6)-ED21(2)) * 0.124
         
         CALL SF(T,KondoTemp, cutoffPar, mGr, MExc,
     *   delta_Ground_1st,arg_SF)


         IF((1./arg_SF).EQ.0.E0) THEN
             EM=-EXP(-1.E0)
             WWW = 1./WAPR(EM, 0, NERROR, 0, 0)
         ELSE
             WWW = 1./WAPR((1./arg_SF), 0, NERROR, 1, 0)
         ENDIF
         
c         WWW = 1./BISECT(1./arg_SF, 0, NERROR, 0)
         
c         WRITE(*,*) delta_Ground_1st, delta_Ground_2nd
c         WRITE(*,*) BD21(2,6),BD21(6,6)
         
c         WRITE(*,*) T, C/T, Cel/T, Clat/T

c         Xi = (2.*XiX+XiZ)/3.

         XiZ_Scaled = XiZ * (1. - 2.*WWW)

         WRITE(6,*) T, 1./XiZ
c         WRITE(*,*) T, ((1.0/Xi) - 0.50846)
c         WRITE(*,*) T, 1./(XiZ*COS(TETA*PI/180.)+XiX*SIN(TETA*PI/180.))
c         WRITE(*,*) T, Cel

         T = T + 1.
         
73    CONTINUE


C.....................................................................
 1055    write(6,105)
 105     FORMAT (  5('  '/) )
 107     FORMAT (  1('  '/) )
 10      FORMAT (I2,T10,F5.2)
 1000    FORMAT (9G10.5)
 64     FORMAT ( 1F7.2, 6G15.7,3G12.5)
 65      FORMAT ( 2F7.1, 6G13.5, G10.2, 2G13.5 )
 86     FORMAT (/' NCe$Z (MAGNETIZATION AND ZEEMAN)    '
     *    /' AL=',F9.4,'	GJ=',F9.4,/' FI=',F5.0,' TETA=',F5.0)
      CLOSE(5)
      CLOSE(6)
      CLOSE(7)      

         END


C- XXXCE - ПPEOБPAЗOB. X(J)  ( Ce ) ------cm^-1------------ 25.03.10 ---
         SUBROUTINE XXXCE (X,XX)
         REAL*4     X(*),G1,G2,G3
         REAL*8     XX(*)
               G1= -2./35.
               G2=  2./315.
               G3=  0.
            X( 1)= XX( 1) * 1./2.     *     2.* G1
            X( 2)= XX( 2) * 1./8.     *    60.* G2
            X( 3)= XX( 3) * 1./16.    *     0.* G3
            X( 4)= XX( 4) * 1./0.8165 *     1.* G1
            X( 5)= XX( 5) * 1./1.2649 *     3.* G2
            X( 6)= XX( 6) * 1./1.5614 *     0.* G3
            X( 7)= XX( 7) * 1./0.9562 *    12.* G2
            X( 8)= XX( 8) * 1./1.4254 *     0.* G3
            X( 9)= XX( 9) * 1./1.0527 *     0.* G3
         RETURN
         END
C- XXXCE1 - ПPEOБPAЗOB. X(J)  ( Ce ) ------cm^-1----------- 25.03.10 ---
C- with Lambda_qk ratios of CF parameters
         SUBROUTINE XXXCE1 (X,XX)
         REAL*4     X(*)
         REAL*8     XX(*)
               G1= -2./35.
               G2=  2./315.
               G3=  0.
            X( 1)= XX( 1) *     2.* G1
            X( 2)= XX( 2) *    60.* G2
            X( 3)= XX( 3) *     0.* G3
            X( 4)= XX( 4) *     1.* G1
            X( 5)= XX( 5) *     3.* G2
            X( 6)= XX( 6) *     0.* G3
            X( 7)= XX( 7) *    12.* G2
            X( 8)= XX( 8) *     0.* G3
            X( 9)= XX( 9) *     0.* G3
         RETURN
         END
C- XXXCE3 - ПPEOБPAЗOB. X(J)  ( Ce ) ------cm^-1------------ 25.03.10 ---
         SUBROUTINE XXXCE3 (X,XX)
         REAL*4     X(*),G1,G2,G3
         REAL*8     XX(*)
               G1= -2./35.
               G2=  2./315.
               G3=  0.
            X( 1)= XX( 1) * 2.     *     2.* G1
            X( 2)= XX( 2) * 8.     *    60.* G2
            X( 3)= XX( 3) * 16.    *     0.* G3
            X( 4)= XX( 4) * 0.8165 *     1.* G1
            X( 5)= XX( 5) * 1.2649 *     3.* G2
            X( 6)= XX( 6) * 1.5614 *     0.* G3
            X( 7)= XX( 7) * 0.9562 *    12.* G2
            X( 8)= XX( 8) * 1.4254 *     0.* G3
            X( 9)= XX( 9) * 1.0527 *     0.* G3
         RETURN
         END
C
C- XXXCE2 - ПPEOБPAЗOB. X(J)  ( Ce ) ---Kelvin------------- 25.03.10 ---
         SUBROUTINE XXXCE2 (X,XX)
         REAL*4     X(*),G1,G2,G3
         REAL*8     XX(*)
               G1= -2./35.
               G2=  2./315.
               G3=  0.
            X( 1)= XX( 1) * 1./2.     *     2.* G1/1.439
            X( 2)= XX( 2) * 1./8.     *    60.* G2/1.439
            X( 3)= XX( 3) * 1./16.    *     0.* G3/1.439
            X( 4)= XX( 4) * 1./0.8165 *     1.* G1/1.439
            X( 5)= XX( 5) * 1./1.2649 *     3.* G2/1.439
            X( 6)= XX( 6) * 1./1.5614 *     0.* G3/1.439
            X( 7)= XX( 7) * 1./0.9562 *    12.* G2/1.439
            X( 8)= XX( 8) * 1./1.4254 *     0.* G3/1.439
            X( 9)= XX( 9) * 1./1.0527 *     0.* G3/1.439
         RETURN
         END
         
C- FMD ----------- ФOPM БAЗOB ДEЙCTB MATP MD
         SUBROUTINE FMD ( X,M )
         REAL*4     X(*),M(21)
         A1 =    5.* X(1)+X(2)
         A2 = -X(1) - 3. * X(2)
         A3 = SQRT(10.) * X(4) + 3. * SQRT(10.) * X(5)
         A4 = -4. * X(1) + 2. * X(2)
         A5 =  3. * SQRT(2.) * X(4) - 5. * SQRT(2.) * X(5)
         A6 = SQRT(5.) * X(7)
         
         DO 1 J=1, 21
1	M(J)=0
              M( 1) = A1
              M( 3) = A2
              M( 4) = A3
              M( 6) = A4
              M( 8) = A5
              M(10) = A4
              M(11) = A6
              M(13) = A5
              M(15) = A2
              M(17) = A6
              M(19) = A3
              M(21) = A1

         RETURN
         END

C- NPZVD - PRINT COÁÑÒ. ÇÍÀ×. È ÑÎÁÑÒ. BEKT. ÄEÉCTB. MATP. N<=8     07.06.89
         SUBROUTINE NPZVD (H,E,B,N)
         DIMENSION  E(N), B(N,N)
	     REAL*4 E
         REAL*4 H
         
         WRITE (7, 11) H, (E(J), J=1,N)
         WRITE (7, 14)
         DO 10 J=1,N
 10         WRITE (7, 12) (B(J,I),I=1,N)
 11      FORMAT( 1F7.1,136(F13.4))
 12      FORMAT( 1X,136F9.4)
 14      FORMAT( 1X,' ')
 15      RETURN
         END
 
C-FM -------------------------------------------------------------------
       SUBROUTINE FM ( FI,TETA,XM,ZM,SM )
       REAL*4   FI,TETA,XM,YM,ZM,SM
       DATA        PI /3.1415926/
       
       FIR = FI*PI/180.
       TER = TETA*PI/180.
        SF = SIN(FIR)
        CF = COS(FIR)
        ST = SIN(TER)
        CT = COS(TER)

        XM = XM*ST*CF
        ZM = ZM*CT
        SM = XM+ZM 
       RETURN
       END
                
C+ FM1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE FM1 (T,ED,BD, GJ, XDM,ZDM )
       REAL*4  T,ED(6),BD(6,6)
         CALL  FXDM (T,ED,BD,XDM )
         CALL  FZDM (T,ED,BD,ZDM )
          XDM = XDM * GJ
          ZDM = ZDM * GJ
       RETURN
       END
C+ FXDM +++++++++++++++++++++++++++++++++++++++++++++++++++ 29.03.10 +++
       SUBROUTINE FXDM ( T,ED,B,XDM )
       real*4  T,ED(6),B(6,6), XDM
       DATA  S / -1.439 /
             XDM = 0.
             ZD  = 0.
       DO 1 N=1,6
          IF (N.NE.6) GOTO 2
              EED = 1.
              GOTO 3
 2        EP = S*(ED(N)-ED(6))/T
          IF (EP.LE.-140.) GOTO 1
          EED= EXP(EP)
 3        ZD = ZD+EED
          XDM = XDM + EED * 2*(
     *   0.5*SQRT(5.)*(B(1,N)*B(2,N)+B(5,N)*B(6,N))
     *   +  SQRT(2.)*(B(2,N)*B(3,N)+B(4,N)*B(5,N))
     *   +  1.5 * B(3,N)*B(4,N))
 1     CONTINUE
       XDM = XDM/ZD
       RETURN
       END
       

C+ FZDM ++++++++++++++++++++++++++++++++++++++++++++++++++ 26.03.10 +++
       SUBROUTINE FZDM ( T,ED,B,ZDM )
       REAL*4  T,ED(6),B(6,6),ZDM
       DATA  S / -1.439 /
             ZDM = 0.
             ZD  = 0.
        DO 1 N=1,6
           IF (N.NE.6) GOTO 2
               EED = 1.
               GOTO 3
 2         EP = S*(ED(N)-ED(6))/T
           IF (EP.LE.-140.) GOTO 1
              EED= EXP(EP)
 3         ZD = ZD+EED
           ZDM= ZDM + EED * (
     *          2.5*(B(1,N)**2-B(6,N)**2)+1.5*(B(2,N)**2-B(5,N)**2)
     *        + 0.5*(B(3,N)**2-B(4,N)**2) )
 1      CONTINUE
        ZDM = ZDM/ZD
       RETURN
       END
C+ HC ++++++++++++++++++++++++++++++++++++++++++++++++++++ 26.03.10 +++
       SUBROUTINE HC ( Td,T,ED,C,Cel,Clat )
       REAL*4  T,ED(6),C,Cel,Clat,SUM1,SUM2
 
         S=-1.439
         ZZD=1.

          DO 33 J=1,5
              EP = S*(ED(J)-ED(6))/T
              IF (EP.LE.-140.) GOTO 33
              ZZD=ZZD+EXP(EP)
c	WRITE(*,*) ZZD, exp(ep)
33        CONTINUE 
        
         SUM2=0
         SUM1=0
          DO 777 I=1, 5
              EP=S*((ED(I)-ED(6))/T)
              IF (EP.LE.-140.) GOTO 777
              SUM1=SUM1+(ED(I)-ED(6))*EXP(EP)
              SUM2=SUM2+((ED(I)-ED(6))**2)*EXP(EP)

 777      CONTINUE

C Specific heat in (cm^-1)^2/(K^2 f.u.)

       Cel=-SUM1**2/((ZZD**2)*(T**2))+SUM2/((T**2)*ZZD)

C Specific heat in J/(K f.u.) to multiply by 10^-23

       Cel=2.859*Cel

C Lattice (phonon) specific heat in J/(K f.u.) to multiply by 10^-23

       Clat=322.6*((T/Td)**3)
      
       C=Cel+Clat

       RETURN
       END
C+ FMS РАСЧЕТ НАЧАЛЬНОЙ МАГНИТНОЙ ВОСПРИИМЧИВОСТИ+++++++++++29.03.10++++
       SUBROUTINE FMS ( AL1,T,ED,BD,NN,XX,XZ )
         REAL*4  T,ED(NN),BD(NN,NN),XX,XZ,Z,EEZ,EZ,XN,ZN,SXM,SZM,AL1
         REAL*4   XM,ZM
         DATA  S / -1.439 /
           XX  = 0.
           XZ  = 0.
           Z   = 0.
            SXM = 0.
            SZM = 0.
	       XN=0.
	       ZN=0.
	       XM=0.
	       ZM=0.


	   DO 1 N=1, NN

            IF (N.NE.NN)  GOTO 22
               EEZ = 1.
               GOTO 33
 22         EZ = S*(ED(N)-ED(NN))/T

            IF (EZ.LE.-140.) GOTO 1
               EEZ= EXP(EZ)
 33            Z  = Z + EEZ

            CALL FXN ( BD,NN,N,XN )
            CALL FZN ( BD,NN,N,ZN )

            SXM = 0.
            SZM = 0.

            DO 2 M=1,NN

            IF (M.EQ.N .OR. ABS(ED(M)-ED(N)).LT.0.1)  GOTO 2
            
                  CALL FXM ( BD,NN,M,N,XM )
                  CALL FZM ( BD,NN,M,N,ZM )

               SXM=SXM - XM/(S*(ED(M)-ED(N)))
               SZM=SZM - ZM/(S*(ED(M)-ED(N)))

 2          CONTINUE

            XX = XX + EEZ*(XN/T + 2.*SXM)

            XZ = XZ + EEZ*(ZN/T + 2.*SZM)

 1       CONTINUE

         XX = (XX / Z) * AL1
         XZ = (XZ / Z) * AL1

       RETURN
       END
C+ FXN +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE FXN ( B,NN,N,XN )
       REAL*4  B(NN,NN),XN
      XN=(
     *      SQRT(5.)*(B(1,N)*B(2,N)+B(5,N)*B(6,N))
     *  +2.*SQRT(2.)*(B(2,N)*B(3,N)+B(4,N)*B(5,N))
     *  +  3. *                     B(3,N)*B(4,N))**2


       RETURN
       END
C+ FXM +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE FXM ( B,NN,M,N,XM )
       REAL*4  B(NN,NN),XM
      XM=(
     *   0.5* SQRT(5.)*(B(1,N)*B(2,M)+B(2,N)*B(1,M)
     *                 +B(5,N)*B(6,M)+B(6,N)*B(5,M))
     * + SQRT(2.)     *(B(2,N)*B(3,M)+B(3,N)*B(2,M)
     *                + B(4,N)*B(5,M)+B(5,N)*B(4,M))
     * + 1.5 *         (B(3,N)*B(4,M)+B(4,N)*B(3,M)))**2
 
      RETURN
      END
C+ FZN +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FZN ( B,NN,N,ZN )
      REAL*4  B(NN,NN),ZN
      ZN=(
     *          2.5*(B(1,N)**2-B(6,N)**2)
     *        + 1.5*(B(2,N)**2-B(5,N)**2)
     *        + 0.5*(B(3,N)**2-B(4,N)**2))**2

      RETURN
      END
C+ FZM +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FZM ( B,NN,M,N,ZM )
      REAL*4  B(NN,NN),ZM
                
      ZM=(
     *       2.5*(B(1,N)*B(1,M)-B(6,N)*B(6,M))
     *     + 1.5*(B(2,N)*B(2,M)-B(5,N)*B(5,M))
     *     + 0.5*(B(3,N)*B(3,M)-B(4,N)*B(4,M)))**2


      RETURN
      END


C- Zeeman -----------------------------------------29.03.10------------ 
	SUBROUTINE ZEEMAN(ALP,H,M,FI,TETA)
      real*4 M(21),ALP,AL,H,FI,TETA
       DATA        PI /3.1415926/

       FIR = FI*PI/180.
       TER = TETA*PI/180.
        SF = SIN(FIR)
        CF = COS(FIR)
        ST = SIN(TER)
        CT = COS(TER)

      AL = ALP*H/50


              M( 1) = M( 1) + 2.5 *AL *CT
              M( 3) = M( 3) + 1.5 *AL *CT
              M( 6) = M( 6) + 0.5 *AL *CT
              M(10) = M(10) - 0.5 *AL *CT
              M(15) = M(15) - 1.5 *AL *CT
              M(21) = M(21) - 2.5 *AL *CT
              
              M( 2) = M( 2) + 0.5 * SQRT(5.) *AL *ST *CF
              M( 5) = M( 5) +       SQRT(2.) *AL *ST *CF
              M( 9) = M( 9) + 1.5            *AL *ST *CF              
              M(14) = M(14) +       SQRT(2.) *AL *ST *CF
              M(20) = M(20) + 0.5 * SQRT(5.) *AL *ST *CF              
      RETURN
      END
C======================================================================      
C- NPSZVK - PRINT COБ ЗHAЧ И COБ BEKT KOMПЛ MATP. N<=17    28.12.89
         SUBROUTINE NPSZVK ( H,E,Z,N )
         REAL*4   E(N), Z(*),H
         WRITE (7, 100) H, ( E(J), J=1,N)

         WRITE (7, 111)
         L = 0
         N2 = N*2
         DO 1 J=1,N2
            JK =  N2*(N-1)+J
            WRITE (7, 110) ( Z(K),K=J,JK,N2 )
            IF ( L .EQ. 0 ) GOTO 1
               WRITE (7, 111)
               L = -1
  1      L = L + 1
 111     FORMAT( 1X,' ')
 100     FORMAT(1F7.1, 9(G13.5))
 110     FORMAT( 1X,9(F7.2))
         RETURN
         END              
      
C- DMSV6 - PAЗДEЛ ДEЙCTB И MHИM ЧACT COБ BEKT Z HA CR,CJ ------
         SUBROUTINE DMSV6 ( Z,CR,CJ,N)
         REAL*4     Z(*),CR(N, N),CJ(N, N)
         K=1
         DO 1 L=1,N
            DO 1 M=1,N
               CR(M,L) = Z(K)
               CJ(M,L) = Z(K+1)
 1             K = K + 2
         RETURN
         END      
      

C- FK6 - ФOPM ДEЙCT,KOMП MATP MK6 ИЗ M6+ДOБABK
       SUBROUTINE  FK6 ( ALP,FI,TETA,H,MD6,MK6 )
       REAL*4      MD6(21), MK6(42),FI,TETA,H,AL
       DATA        PI / 3.14159265359/

       FIR = FI*PI/180.
       TER = TETA*PI/180.
        SF = SIN(FIR)
        CF = COS(FIR)
        ST = SIN(TER)
        CT = COS(TER)

         AL = ALP*H/50.
         
         DO 2 J=1,21
            L=J*2-1
            MK6(L)=MD6(J)
  2         MK6(L+1)=0.
C............  ДOБABKИ K ДИAГ ЭЛEM-M MK9 ( TOЛЬKO ДEЙCTBИT.)

              MK6( 1) = MK6( 1) + 2.5 *AL *CT
              MK6( 5) = MK6( 5) + 1.5 *AL *CT
              MK6(11) = MK6(11) + 0.5 *AL *CT
              MK6(19) = MK6(19) - 0.5 *AL *CT
              MK6(29) = MK6(29) - 1.5 *AL *CT
              MK6(41) = MK6(41) - 2.5 *AL *CT

C...........  ДOБABKИ K ПPЯMOУГ ЭЛEM-M  MK17
             MK6( 3) = MK6( 3) + 0.5 * SQRT(5.) *AL *ST *CF
             MK6( 9) = MK6( 9) +       SQRT(2.) *AL *ST *CF
             MK6(17) = MK6(17) + 1.5            *AL *ST *CF              
             MK6(27) = MK6(27) +       SQRT(2.) *AL *ST *CF
             MK6(39) = MK6(39) + 0.5 * SQRT(5.) *AL *ST *CF                    
              
             MK6( 4) = MK6( 4) + 0.5 * SQRT(5.) *AL *ST *SF
             MK6(10) = MK6(10) +       SQRT(2.) *AL *ST *SF
             MK6(18) = MK6(18) + 1.5            *AL *ST *SF              
             MK6(28) = MK6(28) +       SQRT(2.) *AL *ST *SF
             MK6(40) = MK6(40) + 0.5 * SQRT(5.) *AL *ST *SF         
             

       RETURN
       END
C-FMK_______РАСЧЕТ ПРОЕКЦИИ НАМАГНИЧЕННОСТИ НА НАПРАВЛЕНИЕ ПОЛЯ
       SUBROUTINE FMK ( FI,TETA,XM,YM,ZM,SM )
       REAL*4   FI,TETA,XM,YM,ZM,SM
       DATA        PI / 3.14159265359/

       FIR = FI*PI/180.
       TER = TETA*PI/180.
        SF = SIN(FIR)
        CF = COS(FIR)
        ST = SIN(TER)
        CT = COS(TER)

        XM = XM*ST*CF
        YM = YM*ST*SF
        ZM = ZM*CT
        SM = XM+YM+ZM
       RETURN
       END
C- FM6 ----------------------------------------------------------------
       SUBROUTINE FM6 (T,EK6,CR6,CJ6, GJ, XKM,YKM,ZKM )
       REAL*4  T,EK6(6),CR6(6,6),CJ6(6,6)
         CALL  FXKM (T,EK6,CR6,CJ6,XKM )
         CALL  FYKM (T,EK6,CR6,CJ6,YKM )
         CALL  FZKM (T,EK6,CR6,CJ6,ZKM )
          XKM = XKM * GJ
          YKM = YKM * GJ
          ZKM = ZKM * GJ
       RETURN
       END

C- FXKM ---------------------------------------------------------------
       SUBROUTINE FXKM (T,EK,CR,CJ,XKM )
       REAL*4  T,EK(6),CR(6,6),CJ(6,6),XKM
       DATA  S / -1.439 /
             XKM = 0.
             ZK  = 0.
        DO 1 N=1,6
           IF (N.NE.1) GOTO 4
              EEK = 1.
              GOTO 5
 4         EP = S*(EK(N)-EK(1))/T
           IF (EP.LE.-140.) GOTO 1
              EEK= EXP(EP)
 5            ZK = ZK+EEK
            XKM= XKM + EEK*2.0*(
     *   0.5*SQRT(5.)*(CR(1,N)*CR(2,N)+ CJ(1,N)*CJ(2,N)
     *                +CR(5,N)*CR(6,N)+ CJ(5,N)*CJ(6,N))
     *   +  SQRT(2.)* (CR(2,N)*CR(3,N)+ CJ(2,N)*CJ(3,N)
     *                +CR(4,N)*CR(5,N)+ CJ(4,N)*CJ(5,N))
     *   +  1.5 *     (CR(3,N)*CR(4,N)+ CJ(3,N)*CJ(4,N)))

 1      CONTINUE
        XKM = XKM/ZK
      RETURN
      END
C- FYKM ---------------------------------------------------------------
      SUBROUTINE FYKM ( T,EK,CR,CJ,YKM )
      REAL*4  T,EK(6),CR(6,6),CJ(6,6),YKM
        DATA  S / -1.439 /
              YKM = 0.
              ZK  = 0.
        DO 1 N=1,6
           IF (N.NE.1) GOTO 4
              EEK = 1.
              GOTO 5
 4         EP = S*(EK(N)-EK(1))/T
           IF (EP.LE.-140.) GOTO 1
              EEK= EXP(EP)
 5            ZK = ZK+EEK
           YKM= YKM + EEK*2.0* (
     *   0.5*SQRT(5.)*(CR(1,N)*CJ(2,N)- CJ(1,N)*CR(2,N)
     *                +CR(5,N)*CJ(6,N)- CJ(5,N)*CR(6,N))
     *   +  SQRT(2.)* (CR(2,N)*CJ(3,N)- CJ(2,N)*CR(3,N)
     *                +CR(4,N)*CJ(5,N)- CJ(4,N)*CR(5,N))
     *   +  1.5 *     (CR(3,N)*CJ(4,N)- CJ(3,N)*CR(4,N)))

 1      CONTINUE
        YKM = YKM/ZK
      RETURN
      END
C- FZKM ---------------------------------------------------------------
      SUBROUTINE FZKM (T,EK,CR,CJ,ZKM )
      REAL*4  T,EK(6),CR(6,6),CJ(6,6),ZKM
        DATA  S / -1.439 /
              ZKM = 0.
              ZK  = 0.
         DO 1 N=1,6
            IF (N.NE.1) GOTO 4
               EEK = 1.
               GOTO 5
 4          EP = S*(EK(N)-EK(1))/T
            IF (EP.LE.-140.) GOTO 1
               EEK= EXP(EP)
 5             ZK = ZK+EEK
            ZKM= ZKM + EEK*(
     *          2.5*(CR(1,N)**2+CJ(1,N)**2-CR(6,N)**2-CJ(6,N)**2)
     *        + 1.5*(CR(2,N)**2+CJ(2,N)**2-CR(5,N)**2-CJ(5,N)**2)
     *        + 0.5*(CR(3,N)**2+CJ(3,N)**2-CR(4,N)**2-CJ(4,N)**2))
          
 1       CONTINUE
         ZKM = ZKM/ZK
      RETURN
      END
 
C+ FKMS РАСЧЕТ НАЧАЛЬНОЙ МАГНИТНОЙ ВОСПРИИМЧИВОСТИ+++++++++++06.04.10++++
       SUBROUTINE FKMS ( AL1,T,EK,CR,CJ,NN,XX,XY,XZ )
         REAL*4  T,EK(NN),CR(NN,NN),CJ(NN,NN),Z,EEZ,EZ,XN,ZN,SXM,SZM,AL1
         REAL*4   XM,ZM,XX,XY,XZ
         DATA  S / -1.439 /
           XX  = 0.
           XY  = 0.
           XZ  = 0.
           Z   = 0.
            SXM = 0.
            SYM = 0.
            SZM = 0.
	       XN=0.
           YN=0.
	       ZN=0.
	       XM=0.
           YM=0.
	       ZM=0.


	   DO 1 N=1, NN

            IF (N.NE.1)  GOTO 22
               EEZ = 1.
               GOTO 33
 22         EZ = S*(EK(N)-EK(1))/T

            IF (EZ.LE.-140.) GOTO 1
               EEZ= EXP(EZ)
 33            Z  = Z + EEZ

            CALL FKXN ( CR,CJ,NN,N,XN )
            CALL FKYN ( CR,CJ,NN,N,YN )
            CALL FKZN ( CR,CJ,NN,N,ZN )

            SXM = 0.
            SYM = 0.
            SZM = 0.

            DO 2 M=1,NN

            IF (M.EQ.N .OR. ABS(EK(M)-EK(N)).LT.0.1)  GOTO 2
            
                  CALL FKXM ( CR,CJ,NN,M,N,XM )
                  CALL FKYM ( CR,CJ,NN,M,N,YM )
                  CALL FKZM ( CR,CJ,NN,M,N,ZM )

               SXM=SXM - XM/(S*(EK(M)-EK(N)))
               SYM=SYM - YM/(S*(EK(M)-EK(N)))
               SZM=SZM - ZM/(S*(EK(M)-EK(N)))

 2          CONTINUE

            XX = XX + EEZ*(XN/T + 2.*SXM)
            XY = XY + EEZ*(YN/T + 2.*SYM)
            XZ = XZ + EEZ*(ZN/T + 2.*SZM)

 1       CONTINUE

         XX = (XX / Z) * AL1
         XY = (XY / Z) * AL1
         XZ = (XZ / Z) * AL1

       RETURN
       END
C+ FKXN +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE FKXN ( CR,CJ,NN,N,XN )
       REAL*4  CR(NN,NN),CJ(NN,NN),XN
      XN=(
     *      SQRT(5.)*(CR(1,N)*CR(2,N)+ CJ(1,N)*CJ(2,N)
     *               +CR(5,N)*CR(6,N)+ CJ(5,N)*CJ(6,N))
     *  +2.*SQRT(2.)*(CR(2,N)*CR(3,N)+ CJ(2,N)*CJ(3,N)
     *               +CR(4,N)*CR(5,N)+ CJ(4,N)*CJ(5,N))
     *   +  3. *     (CR(3,N)*CR(4,N)+ CJ(3,N)*CJ(4,N)))**2


       RETURN
       END
C+ FKXM +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE FKXM ( CR,CJ,NN,M,N,XM )
       REAL*4  CR(NN,NN),CJ(NN,NN),XM
      XM=(
     *   0.5* SQRT(5.)*(CR(1,N)*CR(2,M)+CR(2,N)*CR(1,M)
     * +                CJ(1,N)*CJ(2,M)+CJ(2,N)*CJ(1,M)
     * +                CR(5,N)*CR(6,M)+CR(6,N)*CR(5,M)
     * +                CJ(5,N)*CJ(6,M)+CJ(6,N)*CJ(5,M))
     * + SQRT(2.)   *  (CR(2,N)*CR(3,M)+CR(3,N)*CR(2,M)
     * +                CJ(2,N)*CJ(3,M)+CJ(3,N)*CJ(2,M)
     * +                CR(4,N)*CR(5,M)+CR(5,N)*CR(4,M)
     * +                CJ(4,N)*CJ(5,M)+CJ(5,N)*CJ(4,M))
     * + 1.5 *         (CR(3,N)*CR(4,M)+CR(4,N)*CR(3,M)
     * +                CJ(3,N)*CJ(4,M)+CJ(4,N)*CJ(3,M)))**2

 
      RETURN
      END
      
C+ FKYN +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE FKYN ( CR,CJ,NN,N,YN )
       REAL*4  CR(NN,NN),CJ(NN,NN),YN
      YN=(
     *      SQRT(5.)*(CR(1,N)*CJ(2,N)- CJ(1,N)*CR(2,N)
     *               +CR(5,N)*CJ(6,N)- CJ(5,N)*CR(6,N))
     *  +2.*SQRT(2.)*(CR(2,N)*CJ(3,N)- CJ(2,N)*CR(3,N)
     *               +CR(4,N)*CJ(5,N)- CJ(4,N)*CR(5,N))
     *   +  3. *     (CR(3,N)*CJ(4,N)- CJ(3,N)*CR(4,N)))**2


       RETURN
       END
C+ FKYM +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE FKYM ( CR,CJ,NN,M,N,YM )
       REAL*4  CR(NN,NN),CJ(NN,NN),YM
      YM=(
     *   0.5* SQRT(5.)*(CR(1,N)*CJ(2,M)+CR(2,N)*CJ(1,M)
     *                 -CJ(1,N)*CR(2,M)-CJ(2,N)*CR(1,M)
     * +                CR(5,N)*CJ(6,M)+CR(6,N)*CJ(5,M)
     *                 -CJ(5,N)*CR(6,M)-CJ(6,N)*CR(5,M))
     * + SQRT(2.)   *  (CR(2,N)*CJ(3,M)+CR(3,N)*CJ(2,M)
     *                 -CJ(2,N)*CR(3,M)-CJ(3,N)*CR(2,M)
     * +                CR(4,N)*CJ(5,M)+CR(5,N)*CJ(4,M)
     *                 -CJ(4,N)*CR(5,M)-CJ(5,N)*CR(4,M))
     * + 1.5 *         (CR(3,N)*CJ(4,M)+CR(4,N)*CJ(3,M)
     *                 -CJ(3,N)*CR(4,M)-CJ(4,N)*CR(3,M)))**2

 
      RETURN
      END
      
      
            
      
C+ FKZN +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FKZN ( CR,CJ,NN,N,ZN )
      REAL*4  CR(NN,NN),CJ(NN,NN),ZN
      ZN=(
     *          2.5*(CR(1,N)**2+CJ(1,N)**2-CR(6,N)**2-CJ(6,N)**2)
     *        + 1.5*(CR(2,N)**2+CJ(2,N)**2-CR(5,N)**2-CJ(5,N)**2)
     *        + 0.5*(CR(3,N)**2+CJ(3,N)**2-CR(4,N)**2-CJ(4,N)**2))**2


      RETURN
      END
C+ FKZM +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE FKZM ( CR,CJ,NN,M,N,ZM )
      REAL*4  CR(NN,NN),CJ(NN,NN),ZM
                
      ZM=(
     *       2.5*(CR(1,N)*CR(1,M)+CJ(1,N)*CJ(1,M)
     * -          CR(6,N)*CR(6,M)-CJ(6,N)*CJ(6,M))
     *     + 1.5*(CR(2,N)*CR(2,M)+CJ(2,N)*CJ(2,M)
     * -          CR(5,N)*CR(5,M)-CJ(5,N)*CJ(5,M))
     *     + 0.5*(CR(3,N)*CR(3,M)+CJ(3,N)*CJ(3,M)
     * -          CR(4,N)*CR(4,M)-CJ(4,N)*CJ(4,M)))**2

      RETURN
      END
C-FMXI_______РАСЧЕТ ПРОЕКЦИИ НАМАГНИЧЕННОСТИ НА НАПРАВЛЕНИЕ ПОЛЯ
       SUBROUTINE FMXI ( FI,TETA,XiX,XiY,XiZ,Xi )
       REAL*4   FI,TETA,XiX,XiY,XiZ,Xi
       DATA        PI / 3.14159265359/

       FIR = FI*PI/180.
       TER = TETA*PI/180.
        SF = SIN(FIR)
        CF = COS(FIR)
        ST = SIN(TER)
        CT = COS(TER)

        XiX = XiX*ST*CF
        XiY = XiY*ST*SF
        XiZ = XiZ*CT
        Xi = XiX+XiY+XiZ
       RETURN
       END
C- Entropy---------------------------------------------------22.06.2010-
C -- to be done ------
      SUBROUTINE NEKM (T,EK,EKM )
      REAL*4  T,EK(6),EKM,EEK,EP,ZK
        DATA  S / -1.439 /

       EEK = 0.
       EP = 0.
       EKM = 0.
       ZK  = 0.

        DO 1 N=1,6
           IF (N.NE.1) GOTO 4
              EEK = 1.
              GOTO 5
4          EP = S*(EK(N)-EK(1))/T
           IF (EP.LE.-140.) GOTO 1
              EEK= EXP(EP)
5             ZK = ZK+EEK
1       CONTINUE
 
        DO 2 J = 1, 6
2          SE = SE + S*EK(J)

        EKM = log(ZK)
        
        RETURN
        END
        
C
C     __________________________________________________________________
C
C     Approximating the W function
C     ____________________________
C
      REAL FUNCTION WAPR(X,NB,NERROR,L,M)
C
C     WAPR - output
C     X - argument of W(X)
C     NB is the branch of the W function needed:
C        NB = 0 - upper branch
C        NB <> 0 - lower branch
C
C     NERROR is the output error flag:
C        NERROR = 0 -> routine completed successfully
C        NERROR = 1 -> X is out of range
C
C     Range: -exp(-1) <= X for the upper branch of the W function
C            -exp(-1) < X < 0 for the lower branch of the W function
C
C     L - determines how WAPR is to treat the argument X
C        L = 1 -> X is the offset from -exp(-1), so compute
C                 W(X-exp(-1))
C        L <> 1 -> X is the desired X, so compute W(X)
C
C     M - print messages from WAPR?
C         M = 1 -> Yes
C         M <> 1 -> No
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     NN is the output device number
C
C     NBITS is the number of bits (less 1) in the mantissa of the
C        floating point number number representation of your machine.
C        It is used to determine the level of accuracy to which the W
C        function should be calculated.
C
C        Most machines use a 24-bit matissa for single precision and
C        53-56 bits for double precision. The IEEE standard is 53
C        bits. The Fujitsu VP2200 uses 56 bits. Long word length
C        machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
C        single precision.
C
      IMPLICIT REAL (A-H,O-Z)
      SAVE
      PARAMETER(NN=6)
C
C     The following COMMON statement is needed only when testing this
C     function using BISECT, otherwise it can be removed.
C
      COMMON/WAPCOM/NBITS
      DATA INIT,NITER/0,1/
C     DATA NBITS/23/
C
      

      IF(INIT.EQ.0) THEN
         INIT=1
C
C        Code to calculate NBITS for the host machine. NBITS is the
C        mantissa length less one. This value is chosen to allow for
C        rounding error in the final bit. This need only be run once on
C        any particular machine. It can then be included in the above
C        DATA statement.
C
         DO 10 I=1,2000
            B=2.E0**(-I)
            V=1.E0+B
            IF(V.EQ.1.E0)THEN
               NBITS=I-1
               J=-ALOG10(B)
               IF(M.EQ.1) WRITE(NN,40)NBITS,J
               GO TO 20
            ENDIF
10       CONTINUE
C
C        Remove to here after NBITS has been calculated once
C
C        The case of large NBITS
C
20       IF(NBITS.GE.56) NITER=2
C
C        Various mathematical constants
C
         EM=-EXP(-1.E0)
         EM9=-EXP(-9.E0)
         C13=1.E0/3.E0
         C23=2.E0*C13
         EM2=2.E0/EM
         E12=-EM2
         TB=.5E0**NBITS
         TB2=SQRT(TB)
         X0=TB**(1.E0/6.E0)*.5E0
         X1=(1.E0-17.E0*TB**(2.E0/7.E0))*EM
         AN22=3.6E2/83.E0
         AN11=135./83.E0
         AN3=8.E0/3.E0
         AN4=135.E0/83.E0
         AN5=166.E0/39.E0
         AN6=3167.E0/3549.E0
         S2=SQRT(2.E0)
         S21=2.E0*S2-3.E0
         S22=4.E0-3.E0*S2
         S23=S2-2.E0
      ENDIF
      NERROR=0
      IF(L.EQ.1) THEN
         DELX=X
         IF(DELX.LT.0.E0) THEN
            IF(M.EQ.1) WRITE(NN,50)
            NERROR=1
            RETURN
         ENDIF
         XX=X+EM
         IF(E12*DELX.LT.TB**2.AND.M.EQ.1) WRITE(NN,60)DELX
      ELSE
         IF(X.LT.EM) THEN
            NERROR=1
            RETURN
         ELSE IF(X.EQ.EM) THEN
            WAPR=-1.E0
            RETURN
         ENDIF
         XX=X
         DELX=XX-EM
         IF(DELX.LT.TB2.AND.M.EQ.1) WRITE(NN,70)XX
      ENDIF
      IF(NB.EQ.0) THEN
C
C        Calculations for Wp
C
         IF(ABS(XX).LE.X0) THEN
            WAPR=XX/(1.E0+XX/(1.E0+XX/(2.E0+XX/(.6E0+.34E0*XX))))
            RETURN
         ELSE IF(XX.LE.X1) THEN
            RETA=SQRT(E12*DELX)
            WAPR=RETA/(1.E0+RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*
     &           AN6+AN5))+AN3)))-1.E0
            RETURN
         ELSE IF(XX.LE.2.E1) THEN
            RETA=S2*SQRT(1.E0-XX/EM)
            AN2=4.612634277343749E0*SQRT(SQRT(RETA+
     &          1.09556884765625E0))
            WAPR=RETA/(1.E0+RETA/(3.E0+(S21*AN2+S22)*RETA/
     &           (S23*(AN2+RETA))))-1.E0
         ELSE
            ZL=ALOG(XX)
            WAPR=ALOG(XX/ALOG(XX/ZL**EXP(-1.124491989777808E0/
     &          (.4225028202459761E0+ZL))))
         ENDIF
      ELSE
C
C        Calculations for Wm
C
         IF(XX.GE.0.E0) THEN
            NERROR=1
            RETURN
         ELSE IF(XX.LE.X1) THEN
            RETA=SQRT(E12*DELX)
            WAPR=RETA/(RETA/(3.E0+RETA/(RETA/(AN4+RETA/(RETA*
     &           AN6-AN5))-AN3))-1.E0)-1.E0
            RETURN
         ELSE IF(XX.LE.EM9) THEN
            ZL=ALOG(-XX)
            T=-1.E0-ZL
            TS=SQRT(T)
            WAPR=ZL-(2.E0*TS)/(S2+(C13-T/(2.7E2+
     &           TS*127.0471381349219E0))*TS)
         ELSE
            ZL=ALOG(-XX)
            ETA=2.E0-EM2*XX
            WAPR=ALOG(XX/ALOG(-XX/((1.E0-.5043921323068457E0*
     &           (ZL+1.E0))*(SQRT(ETA)+ETA/3.E0)+1.E0)))
         ENDIF
      ENDIF
      DO 30 I=1,NITER
         ZN=ALOG(XX/WAPR)-WAPR
         TEMP=1.E0+WAPR
         TEMP2=TEMP+C23*ZN
         TEMP2=2.E0*TEMP*TEMP2
         WAPR=WAPR*(1.E0+(ZN/TEMP)*(TEMP2-ZN)/(TEMP2-2.E0*ZN))
30    CONTINUE
      RETURN
40    FORMAT(/,' NBITS is',I4,'.',/,' Expect',I4,
     &       ' significant digits from WAPR.')
50    FORMAT(/,' Warning: the offset x is negative (it must be > 0)')
60    FORMAT(' Warning: For this offset (',E16.8,'),',/,
     &       ' W is negligibly different from -1')
70    FORMAT(' Warning: x (= ',E16.8,') is close to -exp(-1).',/,
     &       ' Enter x as an offset to -exp(-1) for greater accuracy')
      END
C
C     END of WAPR
C     __________________________________________________________________
C
C************************************************************************

C - SF ----------- Scaling Factor ----------------
         SUBROUTINE SF(T,TKondo,A,mGround,MExcited,delta,arg)
         REAL*4 T,TKondo,A, delta, term1, term2, arg, kB
         INTEGER mGround,MExcited
         PARAMETER (kB = 8.61733E-2)
                  
         term1 = (TKondo/(A * T))**mGround
         term2 = ((kB * TKondo + delta)/(A * kB * T + delta))**MExcited
         
         arg = term1 * term2
         
         RETURN
         END

C     _________________________________________________________________
C
C     Solution using bisection
C     ________________________
C
      REAL FUNCTION BISECT(XX,NB,NER,L)
C
C     XX is the argument, W(XX)
C     NB is the branch of the W function needed:
C        NB = 0 - upper branch
C        NB <> 0 - lower branch
C
C     NER is the error condition for the bisection routine
C        NER = 0 - routine completed successfully
C        NER = 1 - routine did not converge (indicates NBITS should be
C                  reduced)
C
C     L - indicates if XX is offset from -exp(-1)
C        L = 1, Yes, it is offset
C        L <> 1, No, true XX is entered
C
C     The parameter TOL, which determines the accuracy of the bisection
C     method, is calculated using NBITS (assuming the final bit is lost
C     due to rounding error).
C
C     N0 is the maximum number of iterations used in the bisection
C     method.
C
C     For XX close to 0 for Wp, the exponential approximation is used.
C     The approximation is exact to O(XX**8) so, depending on the value
C     of NBITS, the range of application of this formula varies. Outside
C     this range, the usual bisection method is used.
C
      IMPLICIT REAL (A-H,O-Z)
      SAVE
      COMMON/WAPCOM/NBITS
      PARAMETER(N0=500)
      NER=0
      IF(L.EQ.1) THEN
         X=XX-EXP(-1.E0)
      ELSE
         X=XX
      ENDIF
      IF(NB.EQ.0) THEN
         IF(ABS(X).LT.1.E0/(2.E0**NBITS)**(1.E0/7.E0)) THEN
            BISECT=X*EXP(-X*EXP(-X*EXP(-X*EXP(-X*EXP(-X*
     &             EXP(-X))))))
            RETURN
         ELSE
            U=CRUDE(X,NB)+1.E-3
            TOL=ABS(U)/2.E0**NBITS
            D=MAX(U-2.E-3,-1.E0)
            DO 10 I=1,N0
               R=.5E0*(U-D)
               BISECT=D+R
               IF(X.LT.EXP(1.E0))THEN
C
C                 Find root using w*exp(w)-x to avoid ln(0) error.
C
                  F=BISECT*EXP(BISECT)-X
                  FD=D*EXP(D)-X
               ELSE
C
C                 Find root using ln(w/x)+w to avoid overflow error.
C
                  F=ALOG(BISECT/X)+BISECT
                  FD=ALOG(D/X)+D
               ENDIF
               IF(F.EQ.0.E0.OR.ABS(R).LE.TOL)RETURN
               IF(FD*F.GT.0.E0)THEN
                  D=BISECT
               ELSE
                  U=BISECT
               ENDIF
10          CONTINUE
         ENDIF
      ELSE
         D=CRUDE(X,NB)-1.E-3
         U=MIN(D+2.E-3,-1.E0)
         TOL=ABS(U)/2.E0**NBITS
         DO 20 I=1,N0
            R=.5E0*(U-D)
            BISECT=D+R
            F=BISECT*EXP(BISECT)-X
            IF(F.EQ.0.E0.OR.ABS(R).LE.TOL) RETURN
            FD=D*EXP(D)-X
            IF(FD*F.GT.0.E0)THEN
               D=BISECT
            ELSE
               U=BISECT
            ENDIF
20       CONTINUE
      ENDIF
      NER=1
      RETURN
      END
C
C     __________________________________________________________________
C
C     Crude approximations for the W function (used by BISECT)
C     ________________________________________________________
C
      REAL FUNCTION CRUDE(XX,NB)
      IMPLICIT REAL (A-H,O-Z)
      COMMON/WAPCOM/NBITS
      DATA INIT/0/
C
      IF(INIT.EQ.0) THEN
         INIT=1
C
C        Various mathematical constants
C
         EM=-EXP(-1.E0)
         EM9=-EXP(-9.E0)
         C13=1.E0/3.E0
         EM2=2.E0/EM
         S2=SQRT(2.E0)
         S21=2.E0*S2-3.E0
         S22=4.E0-3.E0*S2
         S23=S2-2.E0
      ENDIF
      IF(NB.EQ.0) THEN
C
C        Calculations for crude Wp
C
         IF(XX.LE.2.E1) THEN
            RETA=S2*SQRT(1.E0-XX/EM)
            AN2=4.612634277343749E0*SQRT(SQRT(RETA+
     &          1.09556884765625E0))
            CRUDE=RETA/(1.E0+RETA/(3.E0+(S21*AN2+S22)*RETA/
     &            (S23*(AN2+RETA))))-1.E0
         ELSE
            ZL=ALOG(XX)
            CRUDE=ALOG(XX/ALOG(XX/ZL**EXP(-1.124491989777808E0/
     &            (.4225028202459761E0+ZL))))
         ENDIF
      ELSE
C
C        Calculations for crude Wm
C
         IF(XX.LE.EM9) THEN
            ZL=ALOG(-XX)
            T=-1.E0-ZL
            TS=SQRT(T)
            CRUDE=ZL-(2.E0*TS)/(S2+(C13-T/(2.7E2+
     &            TS*127.0471381349219E0))*TS)
         ELSE
            ZL=ALOG(-XX)
            ETA=2.E0-EM2*XX
            CRUDE=ALOG(XX/ALOG(-XX/((1.E0-.5043921323068457E0*
     &            (ZL+1.E0))*(SQRT(ETA)+ETA/3.E0)+1.E0)))
         ENDIF
      ENDIF
      RETURN
      END
C
