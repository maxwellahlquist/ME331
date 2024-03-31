C********************************************************************************
C** PURPOSE:  Solve a set of first order ordinary differential equations of    **
C**           the form dy(i)/dt = F(t,y(1), ..., y(numY) (i = 1, ..., numY).   **
C**                                                                            **
C** INPUT:                                                                     **
C**     eqns: Subroutine that evaluates dy(i)/dt (i = 1, ..., numY), the       **
C**           first derivatives of y(1), ..., y(numY) with respect to t.       **
C**                                                                            **
C**     numY: The number of differential equations to be solved.               **
C**                                                                            **
C**        y: One-dimensional array whose elements are y(1), ..., y(numY).     **
C**                                                                            **
C**        t: Independent variable.                                            **
C**                                                                            **
C**    tStep: Maximum integration stepsize.                                    **
C**                                                                            **
C** absError: Allowable absolute error in y(i)  (i=1, ..., numY).              **
C**                                                                            **
C** relError: Allowable relative error in y(i)  (i=1, ..., numY).              **
C**                                                                            **
C**      com: When com = 2, the Kutta-Merson algorithm (L. Fox, Numerical      **
C**           Solutions of Ordinary and Partial Differential Equations,        **
C**           Palo Alto: Addison-Wesley, 1962, pp. 24-25) is employed to       **
C**           perform the numerical solution of the differential equations.    **
C**           Accordingly, dy(i)/dt (i = 1, ..., numY) are evaluated at        **
C**           every integration boundary, including those at tInitial,         **
C**           tFinal, and ones created when tStep is halved to satisfy the     **
C**           requirements imposed by absError and relError.  Integration      **
C**           is self-starting at each boundary, and the occurrence, at        **
C**           boundaries, of discontinuities in derivatives does not lead      **
C**           to failure of the integration process.                           **
C**                                                                            **
C**           When com = 1, a modified form of the Kutta-Merson algorithm      **
C**           is employed.  It is nearly 20% faster than the one used when     **
C**           com = 2 because no recalculation of derivatives at inte-         **
C**           gration boundaries between tInitial and tFinal takes place.      **
C**           Integration is self-starting at tInitial and tFinal only.        **
C**           Integration may fail if any of dy(i)/dt (i = 1, ..., numY)       **
C**           is discontinuous between tInitial and tFinal.                    **
C**                                                                            **
C**           When com = 0, the function eqns is called and dy(i)/dt           **
C**           (i = 1, ..., numY) are evaluated, but no integration             **
C**           is performed.                                                    **
C**                                                                            **
C** OUTPUT:   The value of  t + tStep  is returned in t, and the values of     **
C**           y(i) at  t + tStep  are returned in y.                           **
C**                                                                            **
C** Portions copyright (c) 1988-2021 by Paul Mitiguy.  Use permitted under the **
C**     3-clause BSD license: https://opensource.org/licenses/BSD-3-Clause.    **
C**     This copyright notice must appear in all copies and distributions.     **
C**     This particular copyright notice applies only to the KUTTA subroutine. **
C********************************************************************************
      SUBROUTINE KUTTA (eqns,numY,Y,T,tStep,absError,relError,com,*)
      EXTERNAL         eqns
      INTEGER          numY, com, NUMCUTS, I
      LOGICAL          STEPDBL, ENTRY
      DOUBLE PRECISION Y(numY), F0, F1, F2, Y1, Y2
      DOUBLE PRECISION T, tStep, absError, relError, ERROR, TEST
      DOUBLE PRECISION tFinal, TT, HC, H, H2, H3, H6, H8
      COMMON/CKUTTA/   F0(100),F1(100),F2(100),Y1(100),Y2(100)
      DATA             HC, NUMCUTS / 0.0D0, 20 /

C**   Make sure F0, F1, F2, Y1, Y2 arrays are sufficiently dimensioned
      if( numY .GE. 100 ) THEN
        WRITE(*,2020)
        RETURN 1
      ENDIF

C**   If com=0, call eqns subroutine and return.
      IF( com .EQ. 0) THEN
        CALL eqns(T, Y, F0, 1)
        RETURN
      ENDIF

C**   Check for initial entry and adjust current value of stepsize.
      IF(numY .EQ. 0) THEN
        HC = tStep
        RETURN
      ENDIF
      IF(tStep .EQ. 0) RETURN 1
      IF(HC*tStep .LT. 0) HC = -HC
      IF(HC .EQ. 0)       HC = tStep

C**   Set local variables
      H = HC
      TT = T + H
      tFinal = T + tStep
      T  = tFinal
      ENTRY = .TRUE.

C**   Check round-off problems.
100   IF( TT+H .EQ. TT ) THEN
        T = TT
        WRITE(*,2010) H, T
        CALL eqns(T, Y, F0, 0)
        RETURN 1
      ENDIF
C**   Main Kutta-Merson step
      H2 = H * 0.5D0
      H3 = H / 3.0D0
      H6 = H / 6.0D0
      H8 = H * 0.125D0
      IF( com .EQ. 2 .OR. ENTRY )  CALL eqns(TT-H, Y, F0, 1)
      ENTRY = .FALSE.
      DO 110  I=1,numY
110     Y1(I) = Y(I) + H3*F0(I)
      CALL eqns(TT-2.0*H3, Y1, F1, 0)
      DO 120  I=1,numY
120     Y1(I) = Y(I) + H6*(F0(I) + F1(I))
      CALL eqns(TT-2.0*H3, Y1, F1, 0)
      DO 130  I=1,numY
130     Y1(I) = Y(I) + H8*(F0(I) + 3.0D0*F1(I) )
      CALL eqns(TT-H2,     Y1, F2, 0)
      DO 140  I=1,numY
140     Y1(I) = Y(I) + H2*(F0(I) - 3.0D0*F1(I)+ 4.0D0*F2(I) )
      CALL eqns(TT,        Y1, F1, 0)
      DO 150  I=1,numY
150     Y2(I) = Y(I) + H6*(F0(I) +  4.0D0*F2(I) + F1(I) )
C**   Assume that step needs to be doubled.  Check error criterion
      STEPDBL = .TRUE.
      DO 160 I=1,numY
        ERROR = DABS(Y1(I) - Y2(I)) * 0.2D0
        TEST  = DABS(Y1(I)) * relError
        IF(ERROR .GE. TEST .AND. ERROR .GE. absError) THEN
          HC = H2
          H  = HC
          TT = TT - H2
          NUMCUTS = NUMCUTS - 1
          IF(NUMCUTS .GE. 0) GO TO 100
          T = TT - H
          WRITE(*,2000) T
          CALL eqns(T, Y, F0, 0)
          RETURN 1
        ENDIF
      IF(STEPDBL .AND. 64.0D0*ERROR .GT. TEST
     &           .AND. 64.0D0*ERROR .GT. absError) STEPDBL=.FALSE.
160   CONTINUE
      DO 170  I = 1,numY
170     Y(I) = Y2(I)
C**   Double the STEPSIZE, maybe.
      IF( STEPDBL .AND. DABS(H+H) .LE. DABS(tStep) .AND.
     &     DABS(TT+H+H) .LE. DABS(tFinal) )  THEN
        HC = H + H
        H  = HC
        NUMCUTS = NUMCUTS + 1
      ENDIF
      IF( TT .EQ. tFinal ) THEN
        CALL eqns(tFinal, Y, F0, 2)
        RETURN
      ENDIF
      TT = TT + H
      IF( (H .GT. 0 .AND. TT .GT. tFinal-0.1D0*H) .OR.
     &    (H .LT. 0 .AND. TT .LT. tFinal-0.1D0*H)  )  THEN
        H  = tFinal - (TT-H)
        TT = tFinal
      ENDIF
      IF( com .EQ. 1 ) THEN
        DO 180  I = 1,numY
180       F0(I) = F1(I)
      ENDIF
      GOTO 100

2000  FORMAT(/1X,'THE STEPSIZE HAS BEEN HALVED TOO MANY TIMES; T = ',
     &1PD12.4,/1X,'ERROR: NUMERICAL INTEGRATION FAILED TO CONVERGE.',//)
2010  FORMAT(/1X,'THE STEPSIZE OF ',1PD22.14,' IS TOO SMALL RELATIVE ',
     &'TO THE TERMINAL TIME OF',/1PD22.14,'.  INTEGRATION HALTED BECA',
     &'USE OF NUMERICAL ROUND-OFF.',/,'THE STEPSIZE MAY HAVE BEEN CUT ',
     &'TOO MANY TIMES.'//)
2020  FORMAT(/1X,'ERROR: INCREASE THE ARRAY SIZE IN KUTTA.'//)
      END

