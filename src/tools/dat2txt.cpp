/* dat2txt.cpp -- Convert a data file of a mcmc run of the pta analysis
                  program to a readable txt file.

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

int main(int argc, char **argv) {
  FILE *pFileIn, *pFileOut;
  char *strFileIn, *strFileOut;
  int nMCMCSteps, nPulsarCount;

  double *pdBuffer1, *pdBuffer2, dBuffer;

  if(argc > 2) {
    strFileIn=argv[1];
    strFileOut=argv[2];
  } else {
    printf("Usage is: dat2txt file1.dat file2.txt\n");
    return 1;
  }

  // Read the files:
  try {
    if(! (pFileIn = fopen(strFileIn, "rb+")) ) throw 1;
    if(! (pFileOut = fopen(strFileOut, "w+")) ) throw 2;

    if(! fread(&(nMCMCSteps) , sizeof(int) , 1 , pFileIn) ) throw 3;
    if(! fread(&(nPulsarCount) , sizeof(int) , 1 , pFileIn) ) throw 4;

    fprintf(pFileOut, "# Output of the data of the mcmc calculation\n");
    fprintf(pFileOut, "nMCMCSteps = %i\n", nMCMCSteps);
    fprintf(pFileOut, "nPulsarCount = %i\n\n", nPulsarCount);

    pdBuffer1 = new double[nPulsarCount];
    pdBuffer2 = new double[nPulsarCount];
    for(int i=0; i<nMCMCSteps; i++) {
      fprintf(pFileOut, "[MC%i]\n", i);
      if(! fread(&dBuffer , sizeof(double) , 1 , pFileIn) ) throw 5;
      fprintf(pFileOut, "dAmp = %e\n", dBuffer);
      if(! fread(&dBuffer , sizeof(double) , 1 , pFileIn) ) throw 6;
      fprintf(pFileOut, "dExp = %e\n", dBuffer);
      if(! fread(pdBuffer1 , sizeof(double) , nPulsarCount , pFileIn) ) throw 7;
      if(! fread(pdBuffer2 , sizeof(double) , nPulsarCount , pFileIn) ) throw 8;
      for(int j=0; j<nPulsarCount; j++) {
        fprintf(pFileOut, "vdAmpi[%i] = %e\n", j, pdBuffer1[j]);
      } // for j
      for(int j=0; j<nPulsarCount; j++) {
        fprintf(pFileOut, "vdExpi[%i] = %e\n", j, pdBuffer2[j]);
      } // for j
      if(! fread(&dBuffer , sizeof(double) , 1 , pFileIn) ) throw 14;
      fprintf(pFileOut, "dLogLik = %e\n\n", dBuffer);
    } // for i
    delete pdBuffer1;
    delete pdBuffer2;
  } catch(int nError) {
    printf("Error number: %i \n", nError);
  } // try
  return 0;
}



