/* parse-tempo2-par.cpp -- Parse a tempo2 parameter file to a readable file by
 * the configfile class.

   Rutger van Haasteren 24 August 2007 haasteren@strw.leidenuniv.nl

   Copyright (C) 2005-2007 Rutger van Haasteren.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.  */


#include <string.h>
#include <stdio.h>


/* Read one line from standard input, */
/* copying it to line array (but no more than max chars). */
/* Does not place terminating \n in line array. */
/* Returns line length, or 0 for empty line, or EOF for end-of-file. */
int getline(FILE* pFile, char line[], int max)
{
  int nch = 0;
  int c;
  max = max - 1;			/* leave room for '\0' */

  while((c = fgetc(pFile)) != EOF)
  {
    if(c == '\n')
      break;

    if(nch < max)
    {
      line[nch] = c;
      nch = nch + 1;
    }
  }

  if(c == EOF && nch == 0)
    return EOF;

  line[nch] = '\0';
  return nch;
}


bool getdataline(FILE* pFile, char line[], int max) {
  int nch;
  nch = getline(pFile, line, max);
  while( nch != EOF ){
    // Do some checks
    if(nch > 2 && line[0] != '#') //***********//
      return true;

    nch = getline(pFile, line, max);
  }
  return false;
} // getdataline


// Parse a tempo2 parameter file to a readable file by the configfile class.
int main(int argc, char **argv) {
  FILE *pFileIn, *pFileOut;
  char *strFileIn, *strFileOut, buf[80];
  char strBufLeft[80], strBufRight[80];

  int nTemp=0, nFiles;
  bool bExponentiate=false;

  for(int i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      // process the switches
      for(int j=1; j<strlen(argv[i]); j++) {
        switch(argv[i][j]) {
        case 'e':
          bExponentiate=true;
          break;
        default:
          break;
        } // switch
      } // for j
    } else {
      nTemp++;
      if(nTemp==1) strFileIn = argv[i];
      if(nTemp==2) strFileOut = argv[i];
    } // if argv
  } // for i

  if(nTemp < 2) {
    printf("Usage is: parse-tempo2-par [-benu] inputfile.txt outputfile.txt\n");
    printf("\n");
    printf("               Switches: -X  None as of yet\n");
//    printf("               Switches: -b  Bound the exponent to be 2.0 < e < 2.8\n");
//    printf("                         -e  Exponentiate the result, exp(min(z)-z)\n");
//    printf("\n");
    return 1;
  } else {
    if(nTemp == 2) {
      nFiles = 1;
    } else {
      nFiles = nTemp-1;
    } // if nTemp
  } // if nTemp

  // Read the files:
  try {
    if(! (pFileIn = fopen(strFileIn, "r+")) ) throw 1;
    if(! (pFileOut = fopen(strFileOut, "w+")) ) throw 2;


    // Read the input file, and write the outputfile:
    while(getline(pFileIn, buf, 80) > 0) {
      if(buf[0] == '#') {
	fprintf(pFileOut, "%s\n", buf);
      } else {
	sscanf(buf, "%s %s", strBufLeft, strBufRight);
	fprintf(pFileOut, "%s = %s\n", strBufLeft, strBufRight);
      } // if buf
    } // while

    fclose(pFileOut);
    fclose(pFileIn);
  } catch(int nError) {
    printf("Error number: %i \n", nError);
  } // try

  return 0;
}
