C- NDIAD -- ��A� �E�CTB MATP N-�O �OP ---------------------------------
C        M - �CX �E�CTB MATP. �TO BEKTOP �� ��EM-B �O CTO���AM
C        E - ��A�-E ��EMEHT�=CO� �HA� (N �T�K)     1 2 4 7 11
C        B - CO�CTB BEKTOP�                          3 5 8 12
C        N - �OP��OK MATP                              6 9 13
C        LSV=0 - CO� �HA�EH � CO�CTB BEKTOP�            10 14
C        LSV=1 - TO��KO CO� �HA�EH                         15
C----------------------------------------------------------------------
         SUBROUTINE NDIAD ( M,E,B,N,LSV)
         INTEGER N,K,J
         REAL*4     M(*),E(N),B(N,N)
C........CO�CTB �HA�EH � CO�CTB BEKTOP�

         CALL  EIGEN ( M, B, N, LSV )
         K=0
         DO 1 J=1,N
            K=K+J
 1          E(J) = M(K)
         RETURN
         END
