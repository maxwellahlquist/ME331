% File: gyrostat.r
% Copyright (c) 1992-2021 by Paul Mitiguy.  All rights reserved.
% Licensed distribution by Motion Genesis LLC (with permision).
%----------------------------------------------------------------------------
IF( Strcmpi(#2#,CYLINDER) == false  &&  Strcmpi(#2#,SPHERE) == false )
  { ECHO("\k\aError: 2nd Argument to GYROSTAT command should be CYLINDER or SPHERE\n") }
ELSEIF( Strcmpi(#1#,FRSTAR) == 0 )
{
IF( Strcmpi(#2#,CYLINDER) == 0 )
  {
-(#5#)*(DOT( GetPartialVelocities(W_#4#_#3#>), ALF_#3#_#NEWTONIANFRAME#>) +&
        DOT( GetPartialVelocities(W_#4#_#NEWTONIANFRAME#>), ALF_#4#_#3#>) +&
        DOT( GetPartialVelocities(W_#3#_#NEWTONIANFRAME#>), CROSS(W_#3#_#NEWTONIANFRAME#>, W_#4#_#3#>)) )
  }
ELSE
  {
(#5#)*( DOT( GetPartialVelocities(W_#3#_#NEWTONIANFRAME#>), ALF_#3#_#NEWTONIANFRAME#>) - &
        DOT( GetPartialVelocities(W_#4#_#NEWTONIANFRAME#>), ALF_#4#_#NEWTONIANFRAME#>) )
  }
}
ELSEIF( Strcmpi(#1#,KE) == 0 )
{
IF( Strcmpi(#2#,CYLINDER) == 0 )
  {0.5*(#5#)*MAGSQ(W_#4#_#3#>) + (#5#)*DOT(W_#3#_#NEWTONIANFRAME#>,W_#4#_#3#>)}
ELSE
  {0.5*(#5#)*( MAGSQ(W_#4#_#NEWTONIANFRAME#>) - MAGSQ(W_#3#_#NEWTONIANFRAME#>) )}
}
ELSEIF( Strcmpi(#1#,ANGMOM) == 0 )  { (#5#) * W_#4#_#3#> }
ELSE { ECHO("\k\aError: 1st Argument to GYROSTAT command should be FRSTAR, KE or ANGMOM\n") }
