C- NDIAD -- ÄÈAÃ ÄEÉCTB MATP N-ÃO ÏOP ---------------------------------
C        M - ÈCX ÄEÉCTB MATP. ÝTO BEKTOP ÈÇ ÝËEM-B ÏO CTOËÁÖAM
C        E - ÄÈAÃ-E ÝËEMEHTÛ=COÁ ÇHA× (N ØTÓK)     1 2 4 7 11
C        B - COÁCTB BEKTOPÛ                          3 5 8 12
C        N - ÏOPßÄOK MATP                              6 9 13
C        LSV=0 - COÁ ÇHA×EH È COÁCTB BEKTOPÛ            10 14
C        LSV=1 - TOËÜKO COÁ ÇHA×EH                         15
C----------------------------------------------------------------------
         SUBROUTINE NDIAD ( M,E,B,N,LSV)
         INTEGER N,K,J
         REAL*4     M(*),E(N),B(N,N)
C........COÁCTB ÇHA×EH È COÁCTB BEKTOPÛ

         CALL  EIGEN ( M, B, N, LSV )
         K=0
         DO 1 J=1,N
            K=K+J
 1          E(J) = M(K)
         RETURN
         END
