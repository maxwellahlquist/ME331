% SPLINE   Fits a function through two points with various techniques.
% Portions copyright (c) 2009-2021 Motion Genesis LLC.  All rights reserved.
%----------------------------------------------------------------------------
%    Argument #1: Spline type: Step, Bell, Transition, Pulse, 1,3,5,7,9
%    Argument #2: Independent variable
%    Argument #3: Initial value of independent variable
%    Argument #4: Final value of independent variable
%    Argument #5: Initial value of function
%    Argument #6: Final   value of function (maximum value for Pulse, Bell)
%
%  For splines of order 3, 5, 7, 9
%    Argument #7: Initial value of function's first derivative
%    Argument #8: Final   value of function's first derivative
%
%  For splines of order 5, 7, 9
%    Argument #9:  Initial value of function's second derivative
%    Argument #10: final   value of function's second derivative
%
%  For splines of order 7, 9
%    Argument #11: Initial value of function's third derivative
%    Argument #12: Final   value of function's third derivative
%
%  For splines of order 9
%    Argument #13: Initial value of function's fourth derivative
%    Argument #14: Final   value of function's fourth derivative
%----------------------------------------------------------------------------
%  TRANSITION: Constructs a function y(t) such that
%             y(t=t0) = y0      y(t=tF) = yF
%            y'(t=t0) = 0      y'(t=tF) = 0
%           y''(t=t0) = 0     y''(t=tF) = 0
%    Type:
%          y = spline(transition, t, t0, tF, y0, yF)
%
%  Which produces:
%      y(t) = y0 + (yF-y0)*dT - (yF-y0)/(2*pi) * sin(2*pi*dT)
%      where dT = (t-t0)/(tF-t0)
%
%  EXAMPLE:
%        Y = spline(transition, X, 0, 1, 0, 10)
%     -> Y = 10*X - 1.5915*sin(6.2832*X)
%----------------------------------------------------------------------------
%  PULSE: Constructs a function y(t) such that
%             y(t=t0) =   y(y=tF) = y0        y(t=tF/2) = yMax
%            y'(t=t0) =  y'(y=tF) = 0        y'(t=tF/2) = 0
%           y''(t=t0) = y''(y=tF) = 0
%    Type:
%          y = spline(pulse, t, t0, tF, y0, yMax)
%
%    Which returns:
%        y = y0 + (yMax-y0)/(1-5*pi/16)*( pi*(-dT+2*dT^3-dT^4)+sin(pi*dT) )
%      where dT = (t-t0)/(tF-t0)
%
%  EXAMPLE:
%       R = spline(pulse, X, 0, 1, 0, 10)
%    -> R = 3442.4*X^3 + 547.88*sin(3.1416*X) - 1721.2*X - 1721.2*X^4
%----------------------------------------------------------------------------
%  1: Constructs a first-order polynomial (line) y(t) such that
%             y(t=to) = y0     y(t=tf) = yF
%    Type:
%          y = spline(1, t, t0, tF, y0, yF)
%
%  EXAMPLE:
%       y = spline(1, t, 0, 10, 4, 24)
%    -> y = 4 + 2*t
%----------------------------------------------------------------------------
%  THIRD-ORDER POLYNOMIAL:
%    To construct a third-order polynomial y(t) such that
%              y(t=to) = a     y(t=tf) = b
%             y'(t=to) = c    y'(t=tf) = d
%    Type:
%          y = spline(3,t,to,tf,a,b,c,d)
%
%  EXAMPLE:
%       Y = spline(3,T,0,100,100,0,1,-1)
%    -> Y = 100 + T + 0.0002*T^3 - 0.04*T^2
%----------------------------------------------------------------------------
%  FIFTH-ORDER POLYNOMIAL:
%    To construct a fifth-order polynomial y(t) such that
%              y(t=to) = a     y(t=tf) = b
%             y'(t=to) = c    y'(t=tf) = d
%            y''(t=to) = e   y''(t=tf) = f
%    Type:
%          y = spline(5,t,to,tf,a,b,c,d,e,f)
%
%  EXAMPLE:
%       Y = spline(5,x,0,1,0,100,0,1,0,0)
%    -> Y = 597*X^5 + 996*X^3 - 1493*X^4
%----------------------------------------------------------------------------
IF( #NUM_ARGS# < 6 ) { echo("\a\kError: The command SPLINE has too few arguments.") }
ELSEIF( Strcmpi(#3#,#4#) == 0 )   % To == Tf
  {
  echo("\a\kError: The initial and final values of an independent variable\ecannot be equal to each other.\n")
  }
ELSEIF( Strcmpi(#1#, STEP ) == 0 )  % Paul Mitiguy + Greg Stevens Step Function
  {
  AUTOZ( .5*(#6#+(#5#)) + 1/pi*(#6#-(#5#))*atan(31.820515953773853*( 2*(#2#) - (#4#+(#3#))) / (#4#-(#3#)) )  )
  }
ELSEIF( Strcmpi(#1#, BELL ) == 0 )  % Paul Mitiguy + Greg Stevens Bell Function
  {
  AUTOZ( #5# + (#6#-(#5#))*EXP(-( 2.1459660262893472*(2*(#2#) - (#4#+(#3#)) ) / (#4#-(#3#)) )^2 )   )
  }
ELSEIF( Strcmpi(#1#, TRANSITION ) == 0 ) % Kane's transition function
  {
% Y(T) = Yo + (Yf-Yo)*DT - (Yf-Yo)/(2*pi)*sin(2*pi*DT)
% where DT =(T-To)/(Tf-To)
  AUTOZ( #5# + ((#6#)-(#5#)) * ((#2#)-(#3#)) / ((#4#)-(#3#))&
         +((#5#)-(#6#))*(0.5/pi)*sin(2*pi*((#2#)-(#3#))/((#4#)-(#3#)) )   )
  }
ELSEIF( Strcmpi(#1#, PULSE ) == 0 ) % Kane's pulse function
  {
% Y(T) = Yo + (Ymax-Yo)/(1-5*pi/16)*(pi*(-DT+2*DT^3-DT^4)+sin(pi*DT))
% where DT =(t-to)/(tf-to)
  AUTOZ(  #5# + ((#6#)-(#5#))/(1-5*pi/16) * &
  (pi*(-((#2#-(#3#))/(#4#-(#3#)))+2*((#2#-(#3#))/(#4#-(#3#)))^3-((#2#-(#3#))/(#4#-(#3#)))^4)+sin(pi*((#2#-(#3#))/(#4#-(#3#)))))   )
  }
ELSEIF( Strcmpi(#1#,1) == 0 || Strcmpi(#1#, LINE ) == 0 )  % Keith 1st order Polynomial
  {
  AUTOZ(  (#6#-(#5#))/(#4#-(#3#))*(#2#)+((#5#)*(#4#)-(#6#)*(#3#))/((#4#)-(#3#))  )
  }
ELSEIF( Strcmpi(#1#,3) == 0 || Strcmpi(#1#, CUBIC ) == 0 )  % Keith 3rd order Polynomial
  {
  IF(#NUM_ARGS# < 8) {echo("\a\kError: The command SPLINE has too few arguments.")}
  ELSE
    {
    ELEMENT(  [1, (#2#), (#2#)^2, (#2#)^3] * &
    AUTOZ(inverse(minors,&
        [1,  (#3#),     (#3#)^2,     (#3#)^3;&
         1,  (#4#),     (#4#)^2,     (#4#)^3;&
         0,    1,      2*(#3#),    3*(#3#)^2;&
         0,    1,      2*(#4#),    3*(#4#)^2]) * &
    [#5#; #6#; #7#; #8#]),  1)
    }
  }
ELSEIF( Strcmpi(#1#,5) == 0 )  % Keith 5th order Polynomial
  {
  IF(#NUM_ARGS# < 10) {echo("\a\kError: The command SPLINE has too few arguments.")}
  ELSE
    {
    ELEMENT( [1, (#2#), (#2#)^2, (#2#)^3, (#2#)^4, (#2#)^5] * &
    AUTOZ(inverse(minors,&
       [1,  (#3#),  (#3#)^2,   (#3#)^3,    (#3#)^4,     (#3#)^5;&
        1,  (#4#),  (#4#)^2,   (#4#)^3,    (#4#)^4,     (#4#)^5;&
        0,    1,    2*(#3#),  3*(#3#)^2,   4*(#3#)^3,   5*(#3#)^4;&
        0,    1,    2*(#4#),  3*(#4#)^2,   4*(#4#)^3,   5*(#4#)^4;&
        0,    0,      2,       6*(#3#),   12*(#3#)^2,  20*(#3#)^3;&
        0,    0,      2,       6*(#4#),   12*(#4#)^2,  20*(#4#)^3])*&
    [#5#; #6#; #7#; #8#; #9#; #10#]), 1)
    }
  }
ELSEIF( Strcmpi(#1#,7) == 0 )  % Keith 7th order Polynomial
  {
  IF(#NUM_ARGS# < 12) {echo("\a\kError: The command SPLINE has too few arguments.")}
  ELSE
    {
   ELEMENT( [1, (#2#), (#2#)^2, (#2#)^3, (#2#)^4, (#2#)^5, (#2#)^6, (#2#)^7]*&
   AUTOZ(inverse([&
   1,(#3#),  (#3#)^2,   (#3#)^3,   (#3#)^4,   (#3#)^5,     (#3#)^6,    (#3#)^7;&
   1,(#4#),  (#4#)^2,   (#4#)^3,   (#4#)^4,   (#4#)^5,     (#4#)^6,    (#4#)^7;&
   0, 1,   2*(#3#),   3*(#3#)^2, 4*(#3#)^3, 5*(#3#)^4,   6*(#3#)^5,  7*(#3#)^6;&
   0, 1,   2*(#4#),   3*(#4#)^2, 4*(#4#)^3, 5*(#4#)^4,   6*(#4#)^5,  7*(#4#)^6;&
   0, 0,      2,      6*(#3#),  12*(#3#)^2,20*(#3#)^3,  30*(#3#)^4, 42*(#3#)^5;&
   0, 0,      2,      6*(#4#),  12*(#4#)^2,20*(#4#)^3,  30*(#4#)^4, 42*(#4#)^5;&
   0, 0,      0,      6,        24*(#3#)  ,60*(#3#)^2, 120*(#3#)^3,210*(#3#)^4;&
   0, 0,      0,      6,        24*(#4#)  ,60*(#4#)^2, 120*(#4#)^3,210*(#4#)^4])*&
   [#5#; #6#; #7#; #8#; #9#; #10#; #11#; #12#]), 1)
    }
  }
ELSEIF( Strcmpi(#1#,9) == 0 )  % Keith 9th order Polynomial
  {
  IF(#NUM_ARGS# < 14) {echo("\a\kError: The command SPLINE has too few arguments.")}
  ELSE
    {
ELEMENT( &
 [1,(#2#),(#2#)^2,(#2#)^3,(#2#)^4,(#2#)^5,(#2#)^6,(#2#)^7,(#2#)^8,(#2#)^9] *&
 AUTOZ(inverse([&
 1,(#3#),(#3#)^2,(#3#)^3,(#3#)^4,(#3#)^5,(#3#)^6,(#3#)^7,(#3#)^8, (#3#)^9;&
 1,(#4#),(#4#)^2,(#4#)^3,(#4#)^4,(#4#)^5,(#4#)^6,(#4#)^7,(#4#)^8, (#4#)^9;&
 0,1,2*(#3#),3*(#3#)^2,4*(#3#)^3,5*(#3#)^4,6*(#3#)^5,7*(#3#)^6,8*(#3#)^7,9*(#3#)^8;&
 0,1,2*(#4#),3*(#4#)^2,4*(#4#)^3,5*(#4#)^4,6*(#4#)^5,7*(#4#)^6,8*(#4#)^7,9*(#4#)^8;&
 0,0,2,6*(#3#),12*(#3#)^2,20*(#3#)^3,30*(#3#)^4,42*(#3#)^5,56*(#3#)^6,72*(#3#)^7;&
 0,0,2,6*(#4#),12*(#4#)^2,20*(#4#)^3,30*(#4#)^4,42*(#4#)^5,56*(#4#)^6,72*(#4#)^7;&
 0,0,0,6,24*(#3#),60*(#3#)^2,120*(#3#)^3,210*(#3#)^4,336*(#3#)^5,504*(#3#)^6;&
 0,0,0,6,24*(#4#),60*(#4#)^2,120*(#4#)^3,210*(#4#)^4,336*(#4#)^5,504*(#4#)^6;&
 0,0,0,0,24,120*(#3#),360*(#3#)^2,840*(#3#)^3,1680*(#3#)^4,3024*(#3#)^5;&
 0,0,0,0,24,120*(#4#),360*(#4#)^2,840*(#4#)^3,1680*(#4#)^4,3024*(#4#)^5])*&
 [#5#; #6#; #7#; #8#; #9#; #10#; #11#; #12#; #13#; #14#]), 1)
    }
  }
ELSE{ ECHO("\a\kError: The first argument of SPLINE cannot be #1#; it must be\eone of 1, 3, 5, 7, 9, STEP, BELL, PULSE, or TRANSITION.\n") }
% end of file
