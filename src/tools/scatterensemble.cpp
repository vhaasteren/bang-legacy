/* scatterensemble.cpp -- Make a scatter-plot of the ML points for use with
   gnuplot

   Rutger van Haasteren 24 August 2007 haasteren@strw.leidenuniv.nl

   Copyright (C) 2005-2007 Rutger van Haasteren.

   This program is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation; either version 3, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License along with
   this program; if not, write to the Free Software Foundation, Inc., 51
   Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.  */


#include <string.h>
#include <stdio.h>


void d_swap(double &x, double &y) {
  double dTemp;
  dTemp = x;
  x = y;
  y = dTemp;
} // d_swap()


void QuickSort(double *pdBufferX, double *pdBufferY, double *pdBufferZ, double *pdBufferComp, int left, int right) {
//  left is the lower index, right is the upper index
//  of the region of array a that is to be sorted
    int i=left, j=right;
    double pivot = pdBufferComp[(left+right)/2];

    //  partition
    while(i <= j) {    
        while (pdBufferComp[i]<pivot) i++; 
        while (pdBufferComp[j]>pivot) j--;
        if (i<=j)
        {
            d_swap(pdBufferX[i], pdBufferX[j]);
            d_swap(pdBufferY[i], pdBufferY[j]);
            d_swap(pdBufferZ[i], pdBufferZ[j]);
            i++; j--;
        }
    }

    //  recursion
    if (left<j) QuickSort(pdBufferX, pdBufferY, pdBufferZ, pdBufferComp, left, j);
    if (i<right) QuickSort(pdBufferX, pdBufferY, pdBufferZ, pdBufferComp, i, right);
} // QuickSort()




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


// Make a scatter-plot
int main(int argc, char **argv) {
  int nDirs=0, nTemp=0, k, l, nParameters;
  FILE *pFile;
  char strBuf[80];
  double *pdAmp, *pdExp;
  double *pdBuf;


  for(int i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      // process the switches
      for(int j=1; j<strlen(argv[i]); j++) {
        switch(argv[i][j]) {
        case 'x':
	  // Place switches here
          break;
        default:
          break;
        } // switch
      } // for j
    } else {
      nDirs++;
    } // if argv
  } // for i

  if(nDirs < 1) {
    printf("Usage is: scatterensemble [-switches] inputdir1 [inputdir2 ...]\n");
    printf("\n");
    printf("               Switches: -x  No switches yes");
    printf("\n");
    return 1;
  } else {
    // Do something here
  } // if nTemp

  pdAmp = new double[nDirs];
  pdExp = new double[nDirs];


  // First read all the ML points
  for(int i=0; i<nDirs; i++) {
    try {
      strcpy(strBuf, argv[i+1]);
      if(strBuf[strlen(strBuf)-1] != '/') strcat(strBuf, "/");
      strcat(strBuf, "optimumdata.dat");
      
      if(! (pFile = fopen(strBuf, "rb+")) ) throw 1;

      if(! fread(&k , sizeof(int) , 1 , pFile) ) throw 2;
      if(! fread(&l , sizeof(int) , 1 , pFile) ) throw 3;

      nParameters = k*2 + 2;
      pdBuf = new double[nParameters];

      if(! fread(pdBuf , sizeof(double) , nParameters , pFile) ) throw 4;
      if(fclose(pFile) ) throw 0;

      pdAmp[i] = pdBuf[0];
      pdExp[i] = pdBuf[1];

      printf("File %s: [%f, %f]\n", strBuf, pdBuf[0], pdBuf[1]);
      delete[] pdBuf;
    } catch(int nError) {
      printf("Error number: %i, in file: %s\n", nError, strBuf);
      pdAmp[i] = 1.0;
      pdExp[i] = 2.330;
//      return 0;
    } // try
  } // for i


  strcpy(strBuf, "workdata");
  if(strBuf[strlen(strBuf)-1] != '/') strcat(strBuf, "/");
  strcat(strBuf, "plotensemble.txt");

  // Now write the ML points in a plot
  try {
    if(! (pFile = fopen(strBuf, "w+")) ) throw 5;

    fprintf(pFile, "# Ensemble plot of ML estimation points\n");
    fprintf(pFile, "# Generated by scatterensemble - Rutger van Haasteren\n");

    for(int i=0; i<nDirs; i++) {
      fprintf(pFile, "%.18e %.18e\n", pdAmp[i], pdExp[i]);
    } // for i

    if(fclose(pFile) ) throw 6;
  } catch(int nError) {
    printf("Error number: %i\n", nError);
    return 0;
  } // try

  delete[] pdAmp;
  delete[] pdExp;
  return 0;
} // main