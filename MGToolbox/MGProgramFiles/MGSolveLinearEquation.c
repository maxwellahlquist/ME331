/* -------------------------------------------------------------------------- **
** PURPOSE:  The matrix equation A*X = B is solved for X, where A is an       **
**           n by n matrix, and X and B are n by 1 matrices.                  **
**                                                                            **
** INPUT:    numberOfEqns: The integer n (number of equations).               **
**           A: An array of n pointers to rows of doubles.                    **
**              In effect, A is an n by n array of double precision numbers.  **
**           B: An n by 1 array of double precision numbers.                  **
**              The elements of B are unchanged on return.                    **
**                                                                            **
** OUTPUT:   On return, X contains the solution, an n by 1 array of doubles.  **
**           Returns an error message if no solution could be found.          **
**                                                                            **
** Copyright (c) 2009-2021 by Paul Mitiguy.  Use permitted under the 3-clause **
**      BSD license: https://opensource.org/licenses/BSD-3-Clause.            **
**      This copyright notice must appear in all copies & distributions.      **
**      This copyright notice applies to the function MGSolveLinearEquation.  **
** -------------------------------------------------------------------------- */
static char*  MGSolveLinearEquation( int numberOfEqns, double* A[], double B[], double X[] )
{
   static char errorMessage[96];
   double      rowScaleFactor[ myNumberOfCoupledLinearEqns ];
   int         i, j, k;

   /* Ensure there is sufficient space for row-scaling factors */
   if( numberOfEqns > myNumberOfCoupledLinearEqns ) return "Error in MGSolveLinearEquation: numberOfEqns > myNumberOfCoupledLinearEqnsIncrease.";

   /* Begin decomposition */
   for( i = 0;  i < numberOfEqns;  i++ )
   {
      /* Find the element in each row with the maximum absolute value */
      double *Ai = A[i];                         /* Short-cut to this row   */
      double rowMax = 0.0;                       /* Max absolute column     */
      for( j = 0;  j < numberOfEqns;  j++ )      /* Check for zero row      */
        if( rowMax < fabs(Ai[j]) )  rowMax = fabs(Ai[j]);

      /* Issue error if row of zeros are found */
      if( rowMax == 0.0 ) { sprintf( errorMessage, "Error in MGSolveLinearEquation: All elements in row %d of coefficient matrix are zero.", i+1 );  return errorMessage; }

      /* Keep track of row scaling factor */
      rowScaleFactor[i] = 1.0 / rowMax;

      /* Keep the B matrix unchanged by copying its elements into x */
      X[i] = B[i];
   }

   /* Change ordering of rows */
   for( k = 0;  k < numberOfEqns;  k++ )
   {
      double largestPivot = 0.0;                 /* Largest relative pivot  */
      int    swapi = k;
      for( i = k;  i < numberOfEqns;  i++ )      /* Check remaining rows    */
      {
         double relSizeThisColumn = fabs(A[i][k]) * rowScaleFactor[i];
         if( relSizeThisColumn > largestPivot ) { swapi = i;  largestPivot = relSizeThisColumn; }
      }

      /* Issue warning if zero pivot is encountered */
      if( largestPivot == 0.0 )  { sprintf( errorMessage, "Error in MGSolveLinearEquation: Zero pivot in column %d during LU-decomposition of COEF matrix.", k+1 );  return errorMessage; }

      /* Maybe switch rows of A and x by changing row pointers and x values */
      if( swapi != k )
      {
         double *swapa = A[k],  swapx = X[k],  swapScale = rowScaleFactor[k];
         A[k] = A[swapi];  A[swapi] = swapa;     /* Change row pointers of A*/
         X[k] = X[swapi];  X[swapi] = swapx;     /* Change row values of X  */
         rowScaleFactor[k] = rowScaleFactor[swapi];  rowScaleFactor[swapi] = swapScale;
      }

      /* Calculate value of pivot and change lower rows */
      if( k < numberOfEqns - 1 )
      {
         double pivotReciprocal = 1.0 / A[k][k];
         for( i = k+1;  i < numberOfEqns;  i++ )
         {
            double *Ak = A[k],  *Ai = A[i];         /* Shortcuts to rows       */
            double ratio = Ai[k] * pivotReciprocal; /* Multiplicative factor   */
            X[i] -= ratio * X[k];                   /* Modify elements of X    */
            for( j = k+1;  j < numberOfEqns;  j++ )  Ai[j] -= ratio * Ak[j];
         }
      }
   }

   /* Solve upperTriangularMatrix * X = Z  */
   for( i = numberOfEqns-1;  i >= 0;  i-- )
   {
      double sum = 0.0, *Ai = A[i];
      for( j = i+1;  j < numberOfEqns;  j++ )  sum += Ai[j] * X[j];
      X[i] = (X[i] - sum) / Ai[i];
   }
   return NULL;
}

