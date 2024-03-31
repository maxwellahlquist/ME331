function [SolutionToAlgebraicEquations,Output] = CodeTest
%===========================================================================
% File: CodeTest.m created Jan 01 2020 by MotionGenesis 5.9.
% Portions copyright (c) 2009-2020 Motion Genesis LLC. are licensed under
% the 3-clause BSD license.  https://opensource.org/licenses/BSD-3-Clause.
% This copyright notice must appear in all copies and distributions.
% MotionGenesis Professional Licensee: Motion Genesis LLC.
%===========================================================================
x=0; y=0;
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;

%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
R                               =  3;                      % meters              Constant
shouldPrintToScreen             =  1;                      % NoUnits             0 or 1
shouldPrintToFile               =  1;                      % NoUnits             0 or 1
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions

% Begin for-loop(s)
for( t = 0:  30 * DEGtoRAD:  360 * DEGtoRAD ),                        % Converted from degrees

SolutionToAlgebraicEquations = DoCalculations;
Output = CalculateOutput;
OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile );

end   % End for-loop

OutputToScreenOrFile( [], 0, 0 );   % Close output files.


%===========================================================================
function SolutionToAlgebraicEquations = DoCalculations
%===========================================================================
COEF = zeros( 2, 2 );
COEF(1,1) = 1;
COEF(1,2) = 1;
COEF(2,1) = 2;
COEF(2,2) = 3;
RHS = zeros( 1, 2 );
RHS(1) = R;
RHS(2) = 4*sin(t);
SolutionToAlgebraicEquations = COEF \ transpose(RHS);

% Update variables after uncoupling equations
x = SolutionToAlgebraicEquations(1);
y = SolutionToAlgebraicEquations(2);


end


%===========================================================================
function Output = CalculateOutput
%===========================================================================
Output = zeros( 1, 3 );
Output(1) = t*RADtoDEG;
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
      fprintf( 1, '\n Output is in the file CodeTest.1\n' );
      fprintf( 1, '\n Note: To automate plotting, issue the command OutputPlot in MotionGenesis.\n' );
      fprintf( 1, '\n To load and plot columns 1 and 2 with a solid line and columns 1 and 3 with a dashed line, enter:\n' );
      fprintf( 1, '    someName = load( ''CodeTest.1'' );\n' );
      fprintf( 1, '    plot( someName(:,1), someName(:,2), ''-'', someName(:,1), someName(:,3), ''--'' )\n' );
   end
   if( shouldPrintToScreen ),
      fprintf( 1, '\n Output returned to the calling function:\n' );
      fprintf( 1, ' x (UNITS)    y (UNITS) ... \n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t              x              y\n' );
      fprintf( 1,                '%%   (degrees)      (meters)       (meters)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier(1) = fopen('CodeTest.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file CodeTest.1'); end
      fprintf(FileIdentifier(1), '%% FILE: CodeTest.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t              x              y\n' );
      fprintf(FileIdentifier(1), '%%   (degrees)      (meters)       (meters)\n\n' );
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


%=================================
end    % End of function CodeTest
%=================================
