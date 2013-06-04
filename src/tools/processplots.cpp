/* processplots.cpp -- Process plotting data of the pta program for use with gnuplot

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


// Normalize a 3-column gp data file to have unit area in 2nd column
int main(int argc, char **argv) {
  FILE *pFileIn1, *pFileIn2, *pFileOut;
  char *strFileIn1, *strFileIn2, *strFileOut, buf[80];

  double *pdBufferX, *pdBufferY, *pdBufferZ, *pdTemp;
  double dArea=0, dMean=0, dSigma=0.0;
  int nCount=0, nBufSize=0, i, nTemp=0, nFiles, nMCount=0;
  bool bExponentiate=false, bBoundExp=false, bNormalize=false,
       bUnitize=false, bRepeat=false;

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
        case 'r':
          bRepeat=true;
          break;
        default:
          break;
        } // switch
      } // for j
    } else {
      // Set computer name
      nTemp++;
      if(nTemp==1) strFileIn1 = argv[i];
      if(nTemp==2) strFileIn2 = argv[i];
      if(nTemp==3) strFileOut = argv[i];
    } // if argv
  } // for i

  if(nTemp < 2) {
    printf("Usage is: processplots [-benu] inputfile1.txt [inputfile2.txt] outputfile.txt\n");
    printf("\n");
    printf("               Switches: -b  Bound the exponent to be 2.0 < e < 2.8\n");
    printf("                         -e  Exponentiate the result, exp(min(z)-z)\n");
    printf("                         -n  Normalize the reslult, z/area\n");
    printf("                         -u  Unitize the result, z/max(z)\n");
    printf("\n");
    return 1;
  } else {
    if(nTemp == 2) {
      strFileOut = strFileIn2;
      nFiles = 1;
    } else {
      nFiles = 2;
    } // if nTemp
  } // if nTemp

  // Read the files:
  try {
    if(! (pFileIn1 = fopen(strFileIn1, "r+")) ) throw 1;
    if(nFiles==2) { if(! (pFileIn2 = fopen(strFileIn2, "r+")) ) throw 1; }
    if(! (pFileOut = fopen(strFileOut, "w+")) ) throw 2;


    // First read file 1
    if(getline(pFileIn1, buf, 80) <= 0) throw 3;
    while(buf[0] == '#') {//
      fprintf(pFileOut, "%s\n", buf);
      if(getline(pFileIn1, buf, 80) <= 0) throw 4;
    } // while buf

    pdBufferX = new double[11];
    pdBufferY = new double[11];
    pdBufferZ = new double[11];
    nBufSize = 10;
    sscanf(buf, "%lf %lf %lf", pdBufferX, pdBufferY, pdBufferZ);
    nCount++;

    for(i=1; ; i++) {
      if(i > nBufSize) {
        pdTemp = new double[nBufSize + 11];
        for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferX[j];
        delete[] pdBufferX; pdBufferX = pdTemp;
        pdTemp = new double[nBufSize + 11];
        for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferY[j];
        delete[] pdBufferY; pdBufferY = pdTemp;
        pdTemp = new double[nBufSize + 11];
        for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferZ[j];
        delete[] pdBufferZ; pdBufferZ = pdTemp;
        nBufSize += 10;
      } // if i

      if(getdataline(pFileIn1, buf, 80) == false) break;
      sscanf(buf, "%lf %lf %lf", &pdBufferX[i], &pdBufferY[i], &pdBufferZ[i]);
//      if(fscanf(pFileIn1, "%lf %lf %lf", &pdBufferX[i], &pdBufferY[i], &pdBufferZ[i]) <= 0) break;
//      printf("%e %e %e\n", pdBufferX[i], pdBufferY[i], pdBufferZ[i]);
      nCount++;
    } // for i


    // Now read file 2
    if(nFiles == 2) {
      if(getdataline(pFileIn2, buf, 80) == false) throw 3;
      sscanf(buf, "%lf %lf %lf", &pdBufferX[nCount], &pdBufferY[nCount], &pdBufferZ[nCount]);
      nCount++;

      for(i=nCount; ; i++) {
        if(nCount > nBufSize) {
          pdTemp = new double[nBufSize + 11];
          for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferX[j];
          delete[] pdBufferX; pdBufferX = pdTemp;
          pdTemp = new double[nBufSize + 11];
          for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferY[j];
          delete[] pdBufferY; pdBufferY = pdTemp;
          pdTemp = new double[nBufSize + 11];
          for(int j=0; j<=nBufSize; j++) pdTemp[j] = pdBufferZ[j];
          delete[] pdBufferZ; pdBufferZ = pdTemp;
          nBufSize += 10;
        } // if i

        if(getdataline(pFileIn2, buf, 80) == false) break;
        sscanf(buf, "%lf %lf %lf", &pdBufferX[nCount], &pdBufferY[nCount], &pdBufferZ[nCount]);

//        printf("%e %e %e\n", pdBufferX[nCount], pdBufferY[nCount], pdBufferZ[nCount]);
        nCount++;
      } // for i
    } // if nFiles

    // Now sort the data
    QuickSort(pdBufferX, pdBufferY, pdBufferZ, pdBufferX, 0, nCount-1);
    for(int i=1; i<nCount; i++) if(pdBufferX[i] < pdBufferX[i-1]) { printf("Niet gesorteerd: %i\n",i); break;}

    // Also sort the individual lines
    double dPrevious=pdBufferX[0];
    int nStrideStart=0;
    for(int i=0; i<nCount; i++) {
      if(pdBufferX[i] > dPrevious || i == (nCount-1)) {
        QuickSort(pdBufferX, pdBufferY, pdBufferZ, pdBufferY, nStrideStart, i-1);
        nStrideStart = i;
        dPrevious = pdBufferX[i];
      } // pdBufferX
    } // for i

    double dMin=pdBufferZ[0];
    double dMax=pdBufferZ[0];
    double dArea=0;
    for(int i=1; i<nCount; i++) {
      if( (pdBufferY[i] < 2.8 && pdBufferY[i] > 2.0) || (! bBoundExp) ) {
        if(pdBufferZ[i] < dMin) dMin = pdBufferZ[i];
        if(pdBufferZ[i] > dMax) dMax = pdBufferZ[i];
        dArea += pdBufferZ[i];
      } // if pdBufferY
    } // for i

    printf("dMin: %lf,  dMax: %lf dArea: %lf\n", dMin, dMax, dArea);

    for(int i=0; i<nCount; i++) {
      if(bExponentiate) {
        pdBufferZ[i] = exp(dMin-pdBufferZ[i]);
      } else if(bUnitize) {
        pdBufferZ[i] = pdBufferZ[i]/dMax;
      } else if(bNormalize) {
        pdBufferZ[i] = pdBufferZ[i]/dArea;
      }
    } // for i
// 0.011109, 0.043937, 0.13534, 0.32465, 0.60653, 0.88250

    // Now write the data:
    dPrevious=pdBufferX[0];
    dMean=0;
    dSigma=0;
    nMCount=0;
    for(int i=0; i<nCount; i++) {
      if(pdBufferX[i] > dPrevious) {
        fprintf(pFileOut, "\n");
        dPrevious = pdBufferX[i];
      } // pdBufferX
      if(bRepeat) {
	dMean+=pdBufferZ[i];
	dSigma+=pdBufferZ[i]*pdBufferZ[i];
	nMCount++;
	if(pdBufferY[i] != pdBufferY[ (i==nCount-1) ? 0 : i+1] ) {
	  fprintf(pFileOut, "%.18e %.18e %.18e %.18e\n", pdBufferX[i], pdBufferY[i],
	      dMean/nMCount, sqrt((dSigma/nMCount - dMean*dMean/(nMCount*nMCount))));
	  nMCount = 0;
	  dMean = 0;
	  dSigma = 0;
	} // if nMCount
      } else {
  	fprintf(pFileOut, "%.18e %.18e %.18e\n", pdBufferX[i], pdBufferY[i], pdBufferZ[i]);
      } // if bRepeat
    } // for i

    // Calculate the Area again to be able to set confidence intervals
    dArea = 0;
    for(int i=0; i<nCount; i++) {
      dArea += pdBufferZ[i];
    } // for i
    QuickSort(pdBufferX, pdBufferY, pdBufferZ, pdBufferZ, 0, nCount-1);
    double dCum=0, dSigma1=0, dSigma2=0, dSigma3=0;
    int nStep=0;
    for(int i=nCount-1; i>=0; i--) {
      dCum += pdBufferZ[i];
      switch(nStep) {
	case 0:
	  if(dCum > 0.68268949 * dArea) {
	    nStep++;
	    dSigma1 = pdBufferZ[i];
	  } // if
	  break;
	case 1:
	  if(dCum > 0.95449974 * dArea) {
	    nStep++;
	    dSigma2 = pdBufferZ[i];
	  } // if
	  break;
	case 2:
	  if(dCum > 0.99730024 * dArea) {
	    nStep++;
	    dSigma3 = pdBufferZ[i];
	  } // if
	  break;
	default:
	  break;
      } // switch
    } // for i
    printf("dSigma1: %f,  dSigma2: %f,  dSigma3: %f\n", dSigma1, dSigma2, dSigma3);

    fclose(pFileOut);
    fclose(pFileIn1);
    if(nFiles == 2) fclose(pFileIn2);
  } catch(int nError) {
    printf("Error number: %i \n", nError);
  } // try

  delete[] pdBufferX;
  delete[] pdBufferY;
  delete[] pdBufferZ;
  return 0;
}
