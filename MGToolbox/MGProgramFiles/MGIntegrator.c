/* -------------------------------------------------------------------------- **
** PURPOSE:  Solve a set of first order ordinary differential equations of    **
**           the form dy(i)/dt = f( t, y(1), ..., y(n) )  (i = 1 ... n)       **
**           with a Runga-Kutta-Merson algorithm that numerically integrates  **
**           the differential equations from tStart to tStart + hEntry.       **
**           (L. Fox, Numerical Solutions of Ordinary and Partial Differ-     **
**           ential Equations, Palo Alto: Addison-Wesley, 1962, pp. 24-25)    **
**                                                                            **
** INPUT:                                                                     **
**     eqns  Function that evaluates dy(i)/dt (i = 1 ... n), the first        **
**           derivatives of y(1) ... y(n) with respect to t.                  **
**                                                                            **
**     numY  The number n of differential equations to be solved.             **
**        y  One-dimensional array whose elements are y(1) ... y(n).          **
**                                                                            **
**   tStart  Value of independent variable t on entry to this function.       **
**                                                                            **
**   hEntry  Suggested numerical integration stepsize for this function.      **
**           If *hEntry == 0, the function eqns is called and dy(i)/dt        **
**           (i = 1 ... n) are evaluated, but no integration is performed.    **
**           On return, *hEntry is the actual stepsize the integrator used.   **
**                                                                            **
**    hNext  Returns estimated stepsize for the next call to this function.   **
**                                                                            **
** absError  Allowable absolute error in y(i)  (i = 1 ... n).                 **
** relError  Allowable relative error in y(i)  (i = 1 ... n).                 **
**                                                                            **
** OUTPUT:   Returns an error message if integration failed, otherwise NULL.  **
**       y   The values of y(i) (i=1 ... n) at t + *hEntry are returned in y. **
**  hEntry   The actual stepsize used by integrator is returned in hEntry.    **
**   hNext   The integrator's estimate (based on error analysis of current    **
**           step) is returned in hNext, or 0 if integration failed.          **
**                                                                            **
** Copyright (c) 2009-2021 by Paul Mitiguy.  Use permitted under the 3-clause **
**           BSD license: https://opensource.org/licenses/BSD-3-Clause.       **
**           This copyright notice must appear in all copies & distributions. **
**           This copyright notice applies to the functions MGIntegrator,     **
**           MGIntegrateOneStep, MGIntegrateForwardOrBackward                 **
** -------------------------------------------------------------------------- */
static char*  MGIntegrator( char*(*eqns)(double, double*, double*, char),
                            int numY, double y[], double tStart, double* hEntry, double* hNext,
                            double smallestAllowableStepsize, double absError, double relError )
{
   double f0[myNumberOfODES], f1[myNumberOfODES], f2[myNumberOfODES];
   double y1[myNumberOfODES], y2[myNumberOfODES];
   double h = hEntry ? *hEntry : 0;

   /* Always calculate derivatives at tStart (integration boundary here).   */
   /* If h == 0, just call eqns at tStart and return.                       */
   char* errorMessage = (*eqns)( tStart, y, f0, 1 );
   if( !errorMessage  &&  h == 0 )  return NULL;

   /* Avoid endless loop due to adding tiny h to tStart (precision problem) */
   while( !errorMessage  &&  tStart + 0.1*h != tStart )
   {
      int    i, indexFailedVariable = -1;        /* Variable that failed.   */
      double errorRatioMax = 0;                  /* Testing error criterion */
      double h2 = h * 0.5;                                 /* Half    of h  */
      double h3 = h / 3.0;                                 /* Third   of h  */
      double h6 = h / 6.0;                                 /* Sixth   of h  */
      double h8 = h * 0.125;                               /* Eighth  of h  */
      for( i=0;  i<numY;  i++ )  y1[i] = y[i] + h3*f0[i];
      if( (errorMessage = (*eqns)( tStart+h3, y1, f1, 0 )) != NULL ) break;
      for( i=0;  i<numY;  i++ )  y1[i] = y[i] + h6*(f0[i] + f1[i]);
      if( (errorMessage = (*eqns)( tStart+h3, y1, f1, 0 )) != NULL ) break;
      for( i=0;  i<numY;  i++ )  y1[i] = y[i] + h8*(f0[i] + 3*f1[i]);
      if( (errorMessage = (*eqns)( tStart+h2, y1, f2, 0 )) != NULL ) break;
      for( i=0;  i<numY;  i++ )  y1[i] = y[i] + h2*(f0[i] - 3*f1[i] + 4*f2[i]);
      if( (errorMessage = (*eqns)( tStart+h,  y1, f1, 0 )) != NULL)  break;
      for( i=0;  i<numY;  i++ )  y2[i] = y[i] + h6*(f0[i] + 4*f2[i] + f1[i]);

      /* Both y1[i] and y2[i] provide estimates to the new value of y.      */
      /* Decide if the relative and absolute error tolerances are met.      */
      /* If they are not, reduce the stepsize and restart the integration.  */
      for( i = 0;  i < numY;  i++ )              /* Check all variables     */
      {
         double errorInEstimate = 0.2 * fabs( y2[i] - y1[i] );
         double relTest  = fabs( y2[i] ) * relError;
         double largerOfAbsErrOrRelTest = absError > relTest ? absError : relTest;
         double errorRatio = errorInEstimate / largerOfAbsErrOrRelTest;
         if( errorRatio > errorRatioMax ) { errorRatioMax = errorRatio;  indexFailedVariable = i; }
      }

      /* If integration succeeded, update values of y and t before return.  */
      /* Return actual stepsize and estimate for next integration stepsize. */
      if( errorRatioMax < 1 )
      {
         for( i = 0;  i < numY;  i++ )  y[i] = y2[i];
         *hEntry = h;   *hNext = (errorRatioMax < 1.0/64.0) ? 2*h : h;
         return NULL;
      }

      /* Otherwise, errorRatioMax >= 1, so absError or relTest failed.      */
      /* In other words, errorInEstimate >= largerOfAbsErrOrRelTest.        */
      /* Try to halve stepsize and restart integration.                     */
      if( fabs(h=h2) <= fabs(smallestAllowableStepsize) )
      {
         static char stepsizeCutMessage[128];
         sprintf(errorMessage = stepsizeCutMessage,
                 "Error: Numerical integration stepsize cut too many times at t = %17.9E (variable %d).",
                 tStart, indexFailedVariable);
      }
   }

   /* Check if loop terminated due to numerical round-off.                  */
   if( !errorMessage )
      errorMessage = "Error: Numerical round off makes stepsize h too small relative to tStart, so tStart+h = tStart."
                   "\nIntegration stepsize may have been cut too many times.";

   /* Print error message that numerical integrator failed.    */
   /* If h != 0, call eqns to fill for error display.          */
   printf( "\n Error: Numerical integration failed at t = %17.9E.\n\n", tStart );
   if( h != 0 ) (*eqns)( tStart, y, f0, 1 );
   return errorMessage;
}


/* -------------------------------------------------------------------------- **
** PURPOSE:  Controls the independent variable, e.g., t while solving a set   **
**           of first order ordinary differential equations of the form       **
**           dy(i)/dt = f( t, y(1), ..., y(n) )  (i = 1, ..., n).             **
**           Discontinuities in functions and t should be handled here.       **
**                                                                            **
** INPUT:    eqns, numY, y, absError, relError are described above.           **
**           If numY = 0, the static value of myPreviousStepsize is reset to  **
**           tStep.  This is useful when integrating several sets of ODEs.    **
**                                                                            **
**        t  The current value of the independent variable.                   **
** tStepMax  Maximum integration stepsize for this integrator step.  On entry **
**           if tStepMax = 0, eqns is called and dy(i)/dt (i=1 ... n) are     **
**           evaluated, but no integration performed.                         **
**                                                                            **
** OUTPUT:   Returns an error message if integration failed, otherwise NULL.  **
**        t  The value of  t + tStep  is returned in t.                       **
**        y  The values of y(i) (i=1 ... n) at t + tStep are returned in y.   **
** -------------------------------------------------------------------------- */
static char*  MGIntegrateOneStep( char*(*eqns)(double, double*, double*, char),
                                  int numY, double y[], double* t, double tStepMax,
                                  double* stepsizeSuggested, double smallestAllowableStepsize,
                                  double absError, double relError )
{
   double hAccumulated = 0;                      /* How far to tStepMax.    */
   double h = *stepsizeSuggested;                /* Current stepsize.       */
   int    isStepFinished = (tStepMax == 0);

   /* Make as many little integration steps as necessary to get to end.     */
   /* Each integration step starts at *t and ends at *t + h.                */
   while( !isStepFinished )
   {
      /* Numerically integrate y[i] and maybe get a smaller value of h.     */
      /* Set hNext to integrator's estimate of next integration step-size.  */
      double hBeforeCall = h, hNext;             /* Suggested stepsize.     */
      char* errorMessage = MGIntegrator( eqns, numY, y, *t, &h, &hNext, smallestAllowableStepsize, absError, relError );
      if( errorMessage ) return errorMessage;    /* Integration failed      */

      /* Any time or function discontinuities should be handled here.       */
      /* Increment value of t using stepsize returned by the integrator.    */
      *t += h;                                   /* Prepare for next step.  */

      /* Change the stepsize if not taking the last step.                   */
      /* Reduce the stepsize if the integrator reduced the stepsize.        */
      /* Increase the stepsize if the integrator suggests a larger one.     */
      /* Update stepsizeSuggested for next call to integrator.              */
      /* Ensure stepsizeSuggested is not larger than tStepMax.              */
      isStepFinished = fabs( (hAccumulated+=h) + smallestAllowableStepsize) > fabs(tStepMax);
      if( !isStepFinished  ||  fabs(h) < fabs(hBeforeCall)  ||  fabs(hNext) > fabs(*stepsizeSuggested) )
      {
         if( fabs(hNext) > fabs(tStepMax) ) hNext = tStepMax;
         *stepsizeSuggested = h = hNext;
      }

      /* Integrator may suggest a stepsize that is bigger than allowed.     */
      /* If next jump is past or very close to end of tStepMax, adjust h.   */
      if( fabs(hAccumulated + 1.1*h) > fabs(tStepMax) )  h = tStepMax - hAccumulated;
   }

   /* Integration finished.  Calculate derivatives of y at final value of t */
   return MGIntegrator( eqns, numY, y, *t, NULL, NULL, 0, 0, 0 );
}


/* ------------------------------------------------------------------------ */
static char*  MGIntegrateForwardOrBackward( int numVariables, double varArrayToIntegrate[], double outputToFill[], double tInitial, double tFinal, double tStepMax, double absError, double relError, int printIntScreen, int printIntFile )
{
   double stepsizeSuggested = tStepMax,   smallestAllowableStepsize = 1.0E-7 * tStepMax;
   double t = tInitial,   tFinalMinusEpsilon = tFinal - smallestAllowableStepsize;
   int    isFirstCall = 1,  isIntegrateForward = (tFinal > tInitial),  integrateDirection = isIntegrateForward ? 1 : -1;
   int    isPrintToScreen,  isPrintToFile,  printCounterScreen = 0,   printCounterFile = 0;

   /* Ensure valid parameters and sign(tStepMax) moves integration in proper direction. */
   char* errorMessage = (numVariables <= 0  ||  !varArrayToIntegrate  ||  absError <= 0  || relError < 0  ||  integrateDirection * tStepMax < 0) ?
                        "Error: Invalid argument to MGIntegrateForwardOrBackward" : NULL;

   /* Initialize integrator with call at t = tInitial, thereafter integrate */
   int isIntegrationFinished = 0;
   while( !isIntegrationFinished  &&  !errorMessage )
   {
      /* Near the end of numerical integration, perhaps take a partial step (decrease tStepMax). */
      if( (isIntegrateForward && t+tStepMax > tFinal)  ||  (!isIntegrateForward && t+tStepMax < tFinal) )
      {
         tStepMax = tFinal - t;
         if( integrateDirection * stepsizeSuggested > fabs(tStepMax) )  stepsizeSuggested = tStepMax;
      }
      /* Due to round-off or other, maybe set stepsizeSuggested = tStepMax. */
      else if( fabs(stepsizeSuggested) >= 0.99 * fabs(tStepMax) ) stepsizeSuggested = tStepMax;

      errorMessage = MGIntegrateOneStep( MGeqns, numVariables, varArrayToIntegrate, &t, (isFirstCall ? 0 : tStepMax), &stepsizeSuggested, smallestAllowableStepsize, absError, relError );
      isIntegrationFinished = errorMessage  ||  (isIntegrateForward && t > tFinalMinusEpsilon)  ||  (!isIntegrateForward && t < tFinalMinusEpsilon);
      if( (isPrintToScreen = printIntScreen  &&  (isIntegrationFinished || --printCounterScreen <= 0))  != 0 ) printCounterScreen = printIntScreen;
      if( (isPrintToFile   = printIntFile    &&  (isIntegrationFinished || --printCounterFile   <= 0))  != 0 ) printCounterFile   = printIntFile;
      OutputToScreenOrFile( t, outputToFill, isPrintToScreen, isPrintToFile );
      isFirstCall = 0;
   }
   return errorMessage;
}

