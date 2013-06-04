/* combinemcmcdatafiles.cpp -- Combine mcmc data files of the pta
                               analysis program

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

#include <stddef.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <string.h>

const char *strVersion="0.40";


int main(int argc, char **argv) {
  FILE **ppFile, *pFileOutput;
  int *pnMCMCSteps;
  int *pnPPPCount;
  int nFile, i, j, nFileLength;
  int nTotalMCMCSteps=0, nPPPulsarCount=0;
  double *pdBuffer, dBuffer;
  char strFileVersion[16];

  if(argc > 1) {
    ppFile = new FILE*[argc-1];
    pnMCMCSteps = new int[argc-1];
    pnPPPCount = new int[argc-1];
  } else {
    printf("Usage is: combine file1.dat file2.dat file3.dat ...\n");
    return 1;
  }

  // Read the files:
  try {
    for(nFile=1; nFile<argc; nFile++) {
      i = nFile - 1;
      if(! (ppFile[i] = fopen(argv[nFile], "rb+")) ) throw 1;

      // Find the file length
      if(fseek(ppFile[i], 0, SEEK_END) ) throw 2;
      nFileLength = ftell(ppFile[i]);
      if(fseek(ppFile[i], 0, SEEK_SET) ) throw 3;

      // Read version and amount of parameters
      if(! fread(strFileVersion , sizeof(char) , 16 , ppFile[i]) ) throw 2;
      if(! fread(&(pnPPPCount[i]) , sizeof(int) , 1 , ppFile[i]) ) throw 3;

      if(i > 0) {
	if(nPPPulsarCount != pnPPPCount[i]) throw 4;
      } else {
	nPPPulsarCount = pnPPPCount[i];
      } // if i

      pnMCMCSteps[i] = (nFileLength-sizeof(char)*16 - sizeof(int)) / (sizeof(double)*(nPPPulsarCount+1));
    if(pnMCMCSteps[i] <= 0) throw 6;
      if(strcmp(strFileVersion, strVersion) != 0) {
	printf("Version/Fileversion of file %s: %s, %s\n", argv[nFile], strVersion, strFileVersion);
      } // if strcmp


      // Keep track of counting and sanity check of the files
      nTotalMCMCSteps+= pnMCMCSteps[i];
    } // for nFile

    // Now we know enough to actually build one big file: total.dat
    pdBuffer = new double[nPPPulsarCount];
    if(! (pFileOutput = fopen("mcmcdata.total.dat", "wb+")) ) throw 5;
    if(! fwrite(strVersion, sizeof(char), 16, pFileOutput) ) throw 6;
    if(! fwrite(&nPPPulsarCount , sizeof(int) , 1 , pFileOutput) ) throw 7;

    for(nFile=1; nFile<argc; nFile++) {
      i = nFile -1;

      for(j=0; j<pnMCMCSteps[i]; j++) {
	// Parameters
	if(! fread(pdBuffer, sizeof(double), nPPPulsarCount, ppFile[i])) throw 8;
	if(! fwrite(pdBuffer, sizeof(double), nPPPulsarCount, pFileOutput)) throw 9;

	// dLogLik
        if(! fread(&dBuffer , sizeof(double) , 1 , ppFile[i]) ) throw 14;
        if(! fwrite(&dBuffer , sizeof(double) , 1 , pFileOutput) ) throw 15;
      } // for j
      if(fclose(ppFile[i]) ) throw 0;
    } // for nFile
      if(fclose(pFileOutput) ) throw 0;

  } catch(int nError) {
    printf("Error number: %i \n", nError);
  } // try

  delete[] ppFile;
  delete[] pnMCMCSteps;
  delete[] pnPPPCount;
  delete[] pdBuffer;
  return 0;
}



