/* normalize.cpp -- normalize a 2-column gp data file to have unit area

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


// Last update:   2006-10-18


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

// Normalize a 2-column gp data file to have unit area
int main(int argc, char **argv) {
  FILE *pFileIn, *pFileOut;
  char *strFileIn, *strFileOut, buf[80];
  double *pdBufferX, *pdBufferY, *pdTemp;
  double dArea=0;
  int nCount=0, nBufSize=0, i;

  if(argc > 2) {
    strFileIn=argv[1];
    strFileOut=argv[2];
  } else {
    printf("Usage is: normalize inputfile.txt outputfile.txt\n");
    return 1;
  }

  // Read the files:
  try {
    if(! (pFileIn = fopen(strFileIn, "r+")) ) throw 1;
    if(! (pFileOut = fopen(strFileOut, "w+")) ) throw 2;

    if(getline(pFileIn, buf, 80) <= 0) throw 3;
    while(buf[0] == '#') {
      fprintf(pFileOut, "%s\n", buf);
//      printf("%s\n", buf);
      if(getline(pFileIn, buf, 80) <= 0) throw 4;
    } // while buf

    pdBufferX = new double[11];
    pdBufferY = new double[11];
    nBufSize = 10;
    sscanf(buf, "%lf %lf", pdBufferX, pdBufferY);
//    printf("%e %e\n", pdBufferX[0], pdBufferY[0]);
    nCount++;

    for(i=1; ; i++) {
      if(i > nBufSize) {
        pdTemp = new double[nBufSize + 11];
        for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferX[j];
        delete[] pdBufferX; pdBufferX = pdTemp;
        pdTemp = new double[nBufSize + 11];
        for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferY[j];
        delete[] pdBufferY; pdBufferY = pdTemp;
        nBufSize += 10;
      } // if i

      if(fscanf(pFileIn, "%lf %lf", &pdBufferX[i], &pdBufferY[i]) <= 0) break;
//      printf("%e %e\n", pdBufferX[i], pdBufferY[i]);
      nCount++;
    } // for i

    for(i=1; i<nCount; i++) {
      dArea += pdBufferY[i] * (pdBufferX[i] - pdBufferX[i-1]);
    } // for i

    for(i=0; i<nCount; i++) {
      pdBufferY[i] /= dArea;;
    } // for i

    for(int i=0; i<nCount; i++) {
      fprintf(pFileOut, "%e %e\n", pdBufferX[i], pdBufferY[i]);
    } // for i

    delete[] pdBufferX;
    delete[] pdBufferY;
    fclose(pFileOut);
    fclose(pFileIn);
  } catch(int nError) {
    printf("Error number: %i \n", nError);
  } // try

  return 0;
}



