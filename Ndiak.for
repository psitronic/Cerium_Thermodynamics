C- NDIAK - ��A� KOM�� MATP  ------------------------------ 14.11.87 --
C        M - �CX KOM�� MATP-BEKTOP �O CTO���AM     1R
C        E - CO� �HA�            ( N ��CE� )       1I
C        Z - CO�CTB BEKTOP� ( 2N X N ��CE� )       2R  3R
C        N - �OP��OK MATP                          2I  3I
C        JOBN= 1 - CO� �HA�EH � CO�CTB BEKTOP�     4R  5R  6R
C        JOBN= 0 - TO��KO CO� �HA�EH               4I  5I  6I
C---------------------------------------------------------------------
         SUBROUTINE NDIAK ( M, E, Z, N ,JOBN )
         REAL*4     M(*),E(N),Z(*),WK(45)
         IZ = N
C................... CO�CTB �HA�EH � CO�CTB BEKTOP�
         CALL EIGCH ( M,N,JOBN,E,Z,IZ,WK,IER )
         RETURN
         END
