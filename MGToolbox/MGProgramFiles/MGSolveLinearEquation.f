C********************************************************************************
C** PURPOSE:  The matrix equation A*X = B is solved for X, where A is an       **
C**           n by n matrix, and X and B are n by 1 matrices.                  **
C**                                                                            **
C** INPUT:    N: The integer number n (number of equations).                   **
C**           A: N by N array of double precision numbers.                     **
C**           B: N by 1 array of double precision numbers.                     **
C**                                                                            **
C** OUTPUT:   Returns through X an n by 1 array of double precision numbers.   **
C**                                                                            **
C** Portions copyright (c) 1988-2021 by Paul Mitiguy.  Use permitted under the **
C**     3-clause BSD license: https://opensource.org/licenses/BSD-3-Clause.    **
C**     This copyright notice must appear in all copies and distributions.     **
C**     This particular copyright notice applies only to the SOLVE subroutine. **
C********************************************************************************
        SUBROUTINE SOLVE(N, A, B, X)
        IMPLICIT DOUBLE PRECISION (A - Z)
        INTEGER N,IPS(100),I,J,K,IP,KP,KP1,NM1,IDXPIV,IP1,IM1,NP1,IBACK
        DIMENSION A(N,N),SCALES(100),B(N),X(N)

C**     Make sure SCALES array is sufficiently dimensioned
        IF(N .GE. 100) GOTO 530

C*************** Beginning of LU decomposition of A ********************
        ZERO = 0.0D0
        DO 5 I=1,N
           IPS(I) = I
           ROWNRM = 0.0D0
           DO 20 J=1,N
              ROWNRM = DMAX1(ROWNRM,DABS(A(I,J)))
   20      CONTINUE
           IF(ROWNRM.EQ.ZERO) GOTO 500
           SCALES(I) = 1.0D0 / ROWNRM
    5   CONTINUE
        NM1 = N-1
        DO 17 K=1,NM1
           BIG = 0.0D0
           DO 11 I=K,N
              IP = IPS(I)
              SIZE = DABS(A(IP,K))*SCALES(IP)
              IF(SIZE .LE. BIG) GO TO 11
              BIG = SIZE
              IDXPIV = I
   11      CONTINUE
           IF(BIG .EQ. ZERO) GOTO 520
           IF(IDXPIV .EQ. K) GO TO 15
           J = IPS(K)
           IPS(K) = IPS(IDXPIV)
           IPS(IDXPIV) = J
   15      KP = IPS(K)
           PIVOT = A(KP,K)
           KP1 = K+1
           DO 16 I=KP1,N
              IP = IPS(I)
              EM = A(IP,K)/PIVOT
              A(IP,K) = EM
              DO 16 J = KP1,N
                 A(IP,J) = A(IP,J) - EM*A(KP,J)
   16      CONTINUE
   17   CONTINUE
        IF(A(IPS(N),N) .EQ. ZERO) GOTO 520

C**     Note: The LU decomposition of A is returned in A
C***************** Beginning of back substitution **********************
        NP1 = N+1
        X(1) = B(IPS(1))
        DO 2 I=2,N
           IP = IPS(I)
           IM1 = I-1
           SUM = 0.0D0
           DO 1 J=1,IM1
              SUM = SUM + A(IP,J)*X(J)
    1      CONTINUE
           X(I) = B(IP) - SUM
    2   CONTINUE
        X(N) = X(N)/A(IPS(N),N)
        DO 4 IBACK=2,N
           I = NP1-IBACK
           IP = IPS(I)
           IP1 = I+1
           SUM = 0.0D0
           DO 3 J=IP1,N
              SUM = SUM + A(IP,J)*X(J)
    3      CONTINUE
    4   X(I) = (X(I)-SUM)/A(IP,I)
        RETURN

  500  WRITE(*,600) I
       STOP
  520  WRITE(*,620) K
       STOP
  530  WRITE(*,630) N
       STOP
  600  FORMAT(/1X,'Error singular matrix: All elements in row ',I3,
     & '  of COEF are zeros.'/)
  620  FORMAT(/1X,'Error singular matrix: Zero pivot in column ',I3,
     & '  during LU-decomposition of COEF matrix.'/)
  630  FORMAT(/1X,'Error: Increase the SCALES array size in SOLVE',
     & ' to more than ', I3,/)
       END

