function [Output] = CodePointNearCircle
%===========================================================================
% File: CodePointNearCircle.m created Apr 10 2021 by MotionGenesis 5.9.
% Portions copyright (c) 2009-2021 Motion Genesis LLC. are licensed under
% the 3-clause BSD license.  https://opensource.org/licenses/BSD-3-Clause.
% This copyright notice must appear in all copies and distributions.
% MotionGenesis Professional Licensee: Motion Genesis LLC.
%===========================================================================


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
x                               =  4;                      % m                   Guess
y                               =  2;                      % m                   Guess

absError                        =  1.0E-09;                %                     Absolute Error
shouldPrintToScreen             =  1;                      % NoUnits             0 or 1
shouldPrintToFile               =  1;                      % NoUnits             0 or 1
%-------------------------------+--------------------------+-------------------+-----------------


VAR(1) = x;
VAR(2) = y;
SolutionToAlgebraicEquations = SolveSetOfNonlinearAlgebraicEquations( transpose(VAR), absError );
ErrorInConvergenceOfSolution = CalculateFunctionEvaluatedAtX( SolutionToAlgebraicEquations );
MaxErrorInConvergence = max( abs(ErrorInConvergenceOfSolution) );
Output = CalculateOutput( MaxErrorInConvergence );
OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files.


%===========================================================================
function RHS = CalculateFunctionEvaluatedAtX( VAR )
%===========================================================================
% Update variables with new values
x = VAR(1);
y = VAR(2);
RHS(1) = 9 - x^2 - y^2;

end


%===========================================================================
function GradientAtX = CalculateGradientAtX( VAR )
%===========================================================================
RHS = CalculateFunctionEvaluatedAtX( VAR );
COEF = zeros( 1, 2 );
COEF(1,1) = 2*x;
COEF(1,2) = 2*y;

GradientAtX = pinv(COEF) * transpose(RHS);
end



%===========================================================================
function Output = CalculateOutput( MaxErrorInConvergence )
%===========================================================================
Output = zeros( 1, 3 );
Output(1) = MaxErrorInConvergence;
Output(2) = x;
Output(3) = y;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      fclose( FileIdentifier(1) );
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the file CodePointNearCircle.1\n' );
   end
   if( shouldPrintToScreen ),
      fprintf( 1, '\n Output returned to the calling function:\n' );
      fprintf( 1, ' ResidualSum (UNITS)    x (m)    y (m) ... \n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%  ResidualSum         x              y\n' );
      fprintf( 1,                '%%    (UNITS)          (m)            (m)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier(1) = fopen('CodePointNearCircle.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file CodePointNearCircle.1'); end
      fprintf(FileIdentifier(1), '%% FILE: CodePointNearCircle.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%  ResidualSum         x              y\n' );
      fprintf(FileIdentifier(1), '%%    (UNITS)          (m)            (m)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:3) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:3) );  end
end


%===========================================================================
function WriteNumericalData( fileIdentifier, Output )
%===========================================================================
numberOfOutputQuantities = length( Output );
if( numberOfOutputQuantities > 0 ),
   for( i = 1 : numberOfOutputQuantities ),
      fprintf( fileIdentifier, ' %- 14.6E', Output(i) );
   end
   fprintf( fileIdentifier, '\n' );
end
end



%===========================================================================
function SolutionToAlgebraicEquations = SolveSetOfNonlinearAlgebraicEquations( X, AbsoluteError )
%===========================================================================
alphaBest = 0.1;
for( numberOfGradientCalculations = 0:  1:  100),
   magErrorInSetOfFunctions = CalculateErrorInSetOfFunctions( X );
   if( magErrorInSetOfFunctions < AbsoluteError ), break; end
   dX = CalculateGradientAtX( X );
   if( max(abs(dX)) == 0.0 ), warning('Nonlinear solver failed (dX == 0)'); break;  end
   alphaBest = CalculateBestValueOfAlpha( X, alphaBest, dX, magErrorInSetOfFunctions );
   if( alphaBest == 0.0 ), warning('Nonlinear solver failed (alphaBest == 0)'); break;  end
   X = X + alphaBest*dX;
end
if( numberOfGradientCalculations >= 100 ), warning('Nonlinear solver failed (too many gradient calculations)'); end
SolutionToAlgebraicEquations = X;
end



%===========================================================================
function alphaBest = CalculateBestValueOfAlpha( X, alphaLastIteration, dX, magErrorInSetOfFunctionsAtAlpha0 )
%===========================================================================
if( alphaLastIteration > 0.5 ), alphaLastIteration = 0.5;  end
alpha = [ 0;  0.8*alphaLastIteration;  1.2*alphaLastIteration;  2.0*alphaLastIteration ];
error(4) = CalculateErrorInSetOfFunctions( X + alpha(4) * dX );
error(3) = CalculateErrorInSetOfFunctions( X + alpha(3) * dX );
error(2) = CalculateErrorInSetOfFunctions( X + alpha(2) * dX );
error(1) = magErrorInSetOfFunctionsAtAlpha0;
alphaBest = 0.0;
for( numberOfTries = 0:  1:  20),
   if( error(1) < error(2)  ||  error(2) < error(3) ),
      alpha(4) = alpha(3);  alpha(3) = alpha(2);  alpha(2) = 0.5*(alpha(1)+alpha(2));
      error(4) = error(3);  error(3) = error(2);  error(2) = CalculateErrorInSetOfFunctions( X+alpha(2)*dX );
   elseif( error(3) < error(4) ),
      space12 = alpha(2) - alpha(1);
      space23 = alpha(3) - alpha(2);
      space34 = alpha(4) - alpha(3);
      if( space12 >= 0.99*space23  &&  space12 >= 0.99*space34 ),
         alpha(4) = alpha(3);  alpha(3) = alpha(2);  alpha(2) = 0.5*(alpha(1)+alpha(2));
         error(4) = error(3);  error(3) = error(2);  error(2) = CalculateErrorInSetOfFunctions( X+alpha(2)*dX );
      elseif( space34 >= 0.99*space23 ),
         alpha(4) = 0.5*(alpha(3)+alpha(4));
         error(4) = CalculateErrorInSetOfFunctions( X+alpha(4)*dX );
      else
         alpha(4) = alpha(3);  alpha(3) = 0.5*(alpha(2)+alpha(3));
         error(4) = error(3);  error(3) = CalculateErrorInSetOfFunctions( X+alpha(3)*dX );
      end
   elseif( error(4) < error(1) ),
      alphaBest = alpha(4);   break;
   end
end
end



%===========================================================================
function ErrorInSetOfFunctions = CalculateErrorInSetOfFunctions( X )
%===========================================================================
ErrorInSetOfFunctions = norm( CalculateFunctionEvaluatedAtX(X) );
end


%============================================
end    % End of function CodePointNearCircle
%============================================
