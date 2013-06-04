/* multiplyplots.cpp -- Multiply several plots

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


void QuickSort(double *pdBufferX, double *pdBufferY, double *pdBufferComp, int left, int right) {
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
            i++; j--;
        }
    }

    //  recursion
    if (left<j) QuickSort(pdBufferX, pdBufferY, pdBufferComp, left, j);
    if (i<right) QuickSort(pdBufferX, pdBufferY, pdBufferComp, i, right);
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


// Normalize a 3-column gp data file to have unit area in 2nd column
int main(int argc, char **argv) {
  FILE **ppFileIn, *pFileOut;
  char **pstrFileIn, *strFileOut, buf[80];

  double *pdBufferX, *pdBufferY, *pdTemp;
  double *pdResultX, *pdResultY;
  double dArea=0, dPrevious;
  int nCount=0, nBufSize=0, i, nTemp=0, nFiles, nFileCount=0, nStrideStart;
  bool bExponentiate=false, bBoundExp=false, bNormalize=false, bUnitize=false;

  for(int i=1; i<argc; i++) {
    if(argv[i][0] != '-') nFileCount++;
  } // for i

  if(nFileCount < 2) {
    printf("Usage is: multiplymcmcplots [-benu] inputfile1.txt [inputfile2.txt ...] outputfile.txt\n");
    printf("\n");
    printf("               Switches: -b  Bound the exponent to be 2.0 < e < 2.8\n");
    printf("                         -e  Exponentiate the result, exp(min(z)-z)\n");
    printf("                         -n  Normalize the reslult, z/area\n");
    printf("                         -u  Unitize the result, z/max(z)\n");
    printf("\n");
    return 1;
  } // if nFileCount 

  nFileCount--;
  pstrFileIn = new char *[nFileCount];
  ppFileIn = new FILE *[nFileCount];

  for(int i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      // process the switches
      for(int j=1; j<strlen(argv[i]); j++) {
        switch(argv[i][j]) {
        case 'b':
          bBoundExp=true;
          break;
        case 'e':
          bExponentiate=true;
          break;
        case 'n':
          bNormalize=true;
          break;
        case 'u':
          bUnitize=true;
          break;
        default:
          break;
        } // switch
      } // for j
    } else {
      // Set Filenames
      if(nTemp != nFileCount)
  	pstrFileIn[nTemp] = argv[i];
      else
	strFileOut = argv[i];
      nTemp++;
    } // if argv
  } // for i

  // Read the files:
  try {
    if(! (pFileOut = fopen(strFileOut, "w+")) ) throw 2;
    for(int a=0; a<nFileCount; a++)
      if(! (ppFileIn[a] = fopen(pstrFileIn[a], "r+")) ) throw 1;
    if(nBufSize == 0) {
      pdBufferX = new double[11];
      pdBufferY = new double[11];
      nBufSize = 10;
    } // if nBufSize

    // Loop through all the files
    for(int a=0; a<nFileCount; a++) {
      if(getline(ppFileIn[a], buf, 80) <= 0) throw 3;
      // Scan the first comments
      while(buf[0] == '#') {
	if(a==0) fprintf(pFileOut, "%s\n", buf);
	if(getline(ppFileIn[a], buf, 80) <= 0) throw 4;
      } // while buf

      // The first real line of the file
      sscanf(buf, "%lf %lf", &pdBufferX[nCount], &pdBufferY[nCount]); //, &pdBufferZ[nCount]);
      nCount++;

      // Loop through all datalines
      for(i=1; ; i++) {
	if(nCount > nBufSize) {
	  pdTemp = new double[nBufSize + 11];
	  for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferX[j];
	  delete[] pdBufferX; pdBufferX = pdTemp;
	  pdTemp = new double[nBufSize + 11];
	  for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferY[j];
	  delete[] pdBufferY; pdBufferY = pdTemp;
	  nBufSize += 10;
	} // if i

	if(getdataline(ppFileIn[a], buf, 80) == false) break;
	sscanf(buf, "%lf %lf", &pdBufferX[nCount], &pdBufferY[nCount]); //, &pdBufferZ[nCount]);
  //      if(fscanf(pFileIn1, "%lf %lf %lf", &pdBufferX[i], &pdBufferY[i], &pdBufferZ[i]) <= 0) break;
  //      printf("%e %e %e\n", pdBufferX[i], pdBufferY[i], pdBufferZ[i]);
	nCount++;
      } // for i

      // Now sort the data
      QuickSort(pdBufferX, pdBufferY, pdBufferX, 0, nCount-1);
#if 0
      dPrevious=pdBufferX[0];
      nStrideStart=0;
      for(i=0; i<nCount; i++) {
	if(pdBufferX[i] > dPrevious || i == (nCount-1)) {
	  QuickSort(pdBufferX, pdBufferY, pdBufferY, nStrideStart, i-1);
	  nStrideStart = i;
	  dPrevious = pdBufferX[i];
	} // pdBufferX
      } // for i
#endif

      // if this is the first file, allocate memory for the multiplied buffer
      if(a == 0) {
  	pdResultX = new double[nCount];
  	pdResultY = new double[nCount];
	for(i=0; i<nCount; i++) {
	  pdResultX[i] = pdBufferX[i];
	  pdResultY[i] = 1;
	} // for i
      } // if a

      // Multiply the result by the current file
      for(i=0; i<nCount; i++)
	pdResultY[i] *= pdBufferY[i];

      // Reset the counter
      if(a != nFileCount - 1) nCount = 0;
    } // for a


    double dMin=pdResultY[0];
    double dMax=pdResultY[0];
    double dArea=0;
    for(int i=1; i<nCount; i++) {
      if(pdResultY[i] < dMin) dMin = pdResultY[i];
      if(pdResultY[i] > dMax) dMax = pdResultY[i];
      dArea += pdResultY[i];
    } // for i

    printf("dMin: %lf,  dMax: %lf dArea: %lf\n", dMin, dMax, dArea);

    for(int i=0; i<nCount; i++) {
      if(bExponentiate) {
        pdResultY[i] = exp(dMin-pdResultY[i]);
      } else if(bUnitize) {
        pdResultY[i] = pdResultY[i]/dMax;
      } else if(bNormalize) {
        pdResultY[i] = pdResultY[i]/dArea;
      }
    } // for i
// 0.011109, 0.043937, 0.13534, 0.32465, 0.60653, 0.88250

    // Now write the data:
    dPrevious=pdResultX[0];
    for(int i=0; i<nCount; i++) {
      if(pdResultX[i] > dPrevious) {
//        fprintf(pFileOut, "\n");
        dPrevious = pdResultX[i];
      } // pdResultX
      fprintf(pFileOut, "%.18e %.18e\n", pdResultX[i], pdResultY[i]);
    } // for i

    // Calculate the Area again to be able to set confidence intervals
    dArea = 0;
    for(int i=0; i<nCount; i++) {
      dArea += pdResultY[i];
    } // for i
    QuickSort(pdResultX, pdResultY, pdResultY, 0, nCount-1);
    double dCum=0, dSigma1=0, dSigma2=0, dSigma3=0;
    int nStep=0;
    for(int i=nCount-1; i>=0; i--) {
      dCum += pdResultY[i];
      switch(nStep) {
	case 0:
	  if(dCum > 0.68268949 * dArea) {
	    nStep++;
	    dSigma1 = pdResultY[i];
	  } // if
	  break;
	case 1:
	  if(dCum > 0.95449974 * dArea) {
	    nStep++;
	    dSigma2 = pdResultY[i];
	  } // if
	  break;
	case 2:
	  if(dCum > 0.99730024 * dArea) {
	    nStep++;
	    dSigma3 = pdResultY[i];
	  } // if
	  break;
	default:
	  break;
      } // switch
    } // for i
    printf("dSigma1: %f,  dSigma2: %f,  dSigma3: %f\n", dSigma1, dSigma2, dSigma3);

    fclose(pFileOut);
    for(int a=0; a<nFileCount; a++)
      fclose(ppFileIn[a]);
  } catch(int nError) {
    printf("Error number: %i \n", nError);
  } // try

  delete[] ppFileIn;
  delete[] pdBufferX;
  delete[] pdBufferY;
  delete[] pdResultX;
  delete[] pdResultY;
  return 0;
}
