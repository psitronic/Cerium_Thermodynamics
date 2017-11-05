C- NDIAK - ÑàAÉ KOMèã MATP  ------------------------------ 14.11.87 --
C        M - àCX KOMèã MATP-BEKTOP èO CTOãÅñAM     1R
C        E - COÅ áHAó            ( N óàCEã )       1I
C        Z - COÅCTB BEKTOPõ ( 2N X N óàCEã )       2R  3R
C        N - èOPüÑOK MATP                          2I  3I
C        JOBN= 1 - COÅ áHAóEH à COÅCTB BEKTOPõ     4R  5R  6R
C        JOBN= 0 - TOãúKO COÅ áHAóEH               4I  5I  6I
C---------------------------------------------------------------------
         SUBROUTINE NDIAK ( M, E, Z, N ,JOBN )
         REAL*4     M(*),E(N),Z(*),WK(45)
         IZ = N
C................... COÅCTB áHAóEH à COÅCTB BEKTOPõ
         CALL EIGCH ( M,N,JOBN,E,Z,IZ,WK,IER )
         RETURN
         END
