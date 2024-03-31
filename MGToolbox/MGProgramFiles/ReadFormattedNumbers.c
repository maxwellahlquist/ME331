static char*  ReadFormattedString( FILE* Fptr, double* returnNumber )
{
   static char line[256], *errorMessage = NULL;
   static long lineNumber = 0;   lineNumber++;
   if( !fgets(line, 256, Fptr) ) sprintf( errorMessage = line, "Error: End of input file found while reading line %ld", lineNumber );
   else if( returnNumber  &&  (strlen(line) < 60  ||  (*returnNumber = StringToDouble(line+59)) == HUGE_VAL) )
      { printf("\n\n%s\n", line );  sprintf( errorMessage = line, "Error in input file: Expecting number starting at 60th character of line %ld", lineNumber); }
   return errorMessage;
}

/* ------------------------------------------------------------------------ */
static char*  ReadFormattedNumbers( FILE* Fptr, double* next, ... )
{
   char* errorMessage = NULL;
   va_list args;                       /* Variable argument list  */
   for( va_start(args, next);  next && !errorMessage;  next = va_arg(args, double*) ) errorMessage = ReadFormattedString(Fptr, next);
   va_end( args );                     /* Help function make normal return  */

   /* If no error message, also read to end of line (read newline char).    */
   return errorMessage ? errorMessage : ReadFormattedString(Fptr, NULL);
}

