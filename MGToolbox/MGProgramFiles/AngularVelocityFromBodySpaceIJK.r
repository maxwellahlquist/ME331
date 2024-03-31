% File: AngularVelocityFromBodySpaceIJK.r
% Portions copyright (c) 2009-2021 Motion Genesis LLC.  All rights reserved.
%----------------------------------------------------------------------------
% Syntax: AngularVelocityFromBodySpaceIJK( B, BODY123,   q1, q2, q3 )
% Syntax: AngularVelocityFromBodySpaceIJK( B, SPACE123,  q1, q2, q3 )
%----------------------------------------------------------------------------
IF(#NUM_ARGS# != 5) {echo("\a\k\nError: Expecting 5 arguments for this command.\n")}
ELSEIF( strcmpi(#2#,BODY123) == 0 )
  {
  Vector(#1#,  (Dt(#3#)*cos(#4#)*cos(#5#)+Dt(#4#)*sin(#5#)),  &
              (-Dt(#3#)*cos(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)),  &
              ( Dt(#3#)*sin(#4#)+Dt(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY231) == 0 )
  {
  Vector(#1#,  (Dt(#3#)*sin(#4#)+Dt(#5#)),                    &
               (Dt(#3#)*cos(#4#)*cos(#5#)+Dt(#4#)*sin(#5#)),  &
              (-Dt(#3#)*cos(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY312) == 0 )
  {
  Vector(#1#,  (-Dt(#3#)*cos(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)), &
                (Dt(#3#)*sin(#4#)+Dt(#5#)),                   &
                (Dt(#3#)*cos(#4#)*cos(#5#)+Dt(#4#)*sin(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY132) == 0 )
  {
  Vector(#1#,  (Dt(#3#)*cos(#4#)*cos(#5#)-Dt(#4#)*sin(#5#)),  &
              (-Dt(#3#)*sin(#4#)+Dt(#5#)),                    &
               (Dt(#3#)*cos(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY213) == 0 )
  {
  Vector(#1#,  (Dt(#3#)*cos(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)),  &
               (Dt(#3#)*cos(#4#)*cos(#5#)-Dt(#4#)*sin(#5#)),  &
              (-Dt(#3#)*sin(#4#)+Dt(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY321) == 0 )
  {
  Vector(#1#, (-Dt(#3#)*sin(#4#)+Dt(#5#)),                    &
               (Dt(#3#)*cos(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)),  &
               (Dt(#3#)*cos(#4#)*cos(#5#)-Dt(#4#)*sin(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY121) == 0 )
  {
  Vector(#1#,  (Dt(#3#)*cos(#4#)+Dt(#5#)),                    &
               (Dt(#3#)*sin(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)),  &
               (Dt(#3#)*sin(#4#)*cos(#5#)-Dt(#4#)*sin(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY131) == 0 )
  {
  Vector(#1#,  (Dt(#3#)*cos(#4#)+Dt(#5#)),                    &
              (-Dt(#3#)*sin(#4#)*cos(#5#)+Dt(#4#)*sin(#5#)),  &
               (Dt(#3#)*sin(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY212) == 0 )
  {
  Vector(#1#,  (Dt(#3#)*sin(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)),  &
               (Dt(#3#)*cos(#4#)+Dt(#5#)),                    &
              (-Dt(#3#)*sin(#4#)*cos(#5#)+Dt(#4#)*sin(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY232) == 0 )
  {
  Vector(#1#,  (Dt(#3#)*sin(#4#)*cos(#5#)-Dt(#4#)*sin(#5#)),  &
               (Dt(#3#)*cos(#4#)+Dt(#5#)),                    &
               (Dt(#3#)*sin(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY313) == 0 )
  {
  Vector(#1#,  (Dt(#3#)*sin(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)),  &
               (Dt(#3#)*sin(#4#)*cos(#5#)-Dt(#4#)*sin(#5#)),  &
               (Dt(#3#)*cos(#4#)+Dt(#5#)) )
  }
ELSEIF( strcmpi(#2#,BODY323) == 0 )
  {
  Vector(#1#, (-Dt(#3#)*sin(#4#)*cos(#5#)+Dt(#4#)*sin(#5#)),  &
               (Dt(#3#)*sin(#4#)*sin(#5#)+Dt(#4#)*cos(#5#)),  &
               (Dt(#3#)*cos(#4#)+Dt(#5#)) )
  }
ELSEIF( strcmpi(#2#,SPACE123) == 0 )
  {
  Vector(#1#,  (Dt(#3#)-Dt(#5#)*sin(#4#)),                    &
               (Dt(#4#)*cos(#3#)+Dt(#5#)*sin(#3#)*cos(#4#)),  &
              (-Dt(#4#)*sin(#3#)+Dt(#5#)*cos(#3#)*cos(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE231) == 0 )
  {
  Vector(#1#, (-Dt(#4#)*sin(#3#)+Dt(#5#)*cos(#3#)*cos(#4#)),  &
               (Dt(#3#)-Dt(#5#)*sin(#4#)),                    &
               (Dt(#4#)*cos(#3#)+Dt(#5#)*sin(#3#)*cos(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE312) == 0 )
  {
  Vector(#1#,  (Dt(#4#)*cos(#3#)+Dt(#5#)*sin(#3#)*cos(#4#)),  &
              (-Dt(#4#)*sin(#3#)+Dt(#5#)*cos(#3#)*cos(#4#)),  &
               (Dt(#3#)-Dt(#5#)*sin(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE132) == 0 )
  {
  Vector(#1#,  (Dt(#3#)+Dt(#5#)*sin(#4#)),                    &
               (Dt(#4#)*sin(#3#)+Dt(#5#)*cos(#3#)*cos(#4#)),  &
               (Dt(#4#)*cos(#3#)-Dt(#5#)*sin(#3#)*cos(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE213) == 0 )
  {
  Vector(#1#,  (Dt(#4#)*cos(#3#)-Dt(#5#)*sin(#3#)*cos(#4#)),  &
               (Dt(#3#)+Dt(#5#)*sin(#4#)),                    &
               (Dt(#4#)*sin(#3#)+Dt(#5#)*cos(#3#)*cos(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE321) == 0 )
  {
  Vector(#1#,  (Dt(#4#)*sin(#3#)+Dt(#5#)*cos(#3#)*cos(#4#)),  &
               (Dt(#4#)*cos(#3#)-Dt(#5#)*sin(#3#)*cos(#4#)),  &
               (Dt(#3#)+Dt(#5#)*sin(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE121) == 0 )
  {
  Vector(#1#,  (Dt(#3#)+Dt(#5#)*cos(#4#)),                    &
               (Dt(#4#)*cos(#3#)+Dt(#5#)*sin(#3#)*sin(#4#)),  &
              (-Dt(#4#)*sin(#3#)+Dt(#5#)*cos(#3#)*sin(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE131) == 0 )
  {
  Vector(#1#,  (Dt(#3#)+Dt(#5#)*cos(#4#)),                    &
               (Dt(#4#)*sin(#3#)-Dt(#5#)*cos(#3#)*sin(#4#)),  &
               (Dt(#4#)*cos(#3#)+Dt(#5#)*sin(#3#)*sin(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE212) == 0 )
  {
  Vector(#1#,  (Dt(#4#)*cos(#3#)+Dt(#5#)*sin(#3#)*sin(#4#)),  &
               (Dt(#3#)+Dt(#5#)*cos(#4#)),                    &
               (Dt(#4#)*sin(#3#)-Dt(#5#)*cos(#3#)*sin(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE232) == 0 )
  {
  Vector(#1#, (-Dt(#4#)*sin(#3#)+Dt(#5#)*cos(#3#)*sin(#4#)),  &
               (Dt(#3#)+Dt(#5#)*cos(#4#)),                    &
               (Dt(#4#)*cos(#3#)+Dt(#5#)*sin(#3#)*sin(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE313) == 0 )
  {
  Vector(#1#,  (Dt(#4#)*cos(#3#)+Dt(#5#)*sin(#3#)*sin(#4#)),  &
              (-Dt(#4#)*sin(#3#)+Dt(#5#)*cos(#3#)*sin(#4#)),  &
               (Dt(#3#)+Dt(#5#)*cos(#4#)) )
  }
ELSEIF( strcmpi(#2#,SPACE323) == 0 )
  {
  Vector(#1#,  (Dt(#4#)*sin(#3#)-Dt(#5#)*cos(#3#)*sin(#4#)),  &
               (Dt(#4#)*cos(#3#)+Dt(#5#)*sin(#3#)*sin(#4#)),  &
               (Dt(#3#)+Dt(#5#)*cos(#4#)) )
  }
ELSE {echo("\a\kError: #2# is invalid as the third argument of AngularVelocityFromBodySpaceIJK.")}

