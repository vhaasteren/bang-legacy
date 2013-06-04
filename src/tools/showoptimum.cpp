/* showoptimum.cpp -- show the parameters of a fisherdata.dat file

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
#include <math.h>


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


// Read the file
int main(int argc, char **argv) {
  FILE *pFile;
  int k, n;
  double d;

  try {
    if(! (pFile = fopen("fisherdata.dat", "rb+")) ) throw 1;

    if(! fread(&k , sizeof(int) , 1 , pFile) ) throw 2;
    printf("k: %i\n", k);
    if(! fread(&n , sizeof(int) , 1 , pFile) ) throw 3;
    printf("l: %i\n", n);

    if(! fread(&d , sizeof(double) , 1 , pFile) ) throw 4;
    printf("A: %f\n", d);
    if(! fread(&d , sizeof(double) , 1 , pFile) ) throw 5;
    printf("g: %f\n", d);

    printf("PN: ");
    for(int i=0; i<2*k; i++) {
      if(! fread(&d , sizeof(double) , 1 , pFile) ) throw 6;
      printf("%f  ", d);
    } // for i
    printf("\n");

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    printf("Error: %i\n", nError);
    return 1;
  } // try

  return 0;
}




