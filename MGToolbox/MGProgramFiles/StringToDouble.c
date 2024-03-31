static double StringToDouble( char* s )
{
   char* endNumber;   double x = strtod(s, &endNumber);
   while( isspace(*endNumber) ) endNumber++;
   if(*endNumber || (x == HUGE_VAL || x == -HUGE_VAL)) { x = HUGE_VAL;  printf("\nError: Unable to convert to valid double precision number from string %s", s); }
   return x;
}
