/* processensembleresults1d.cpp -- Process plotting data of the pta program

   Rutger van Haasteren 12 August 2011 haasteren@strw.leidenuniv.nl

   Copyright (C) 2005-2011 Rutger van Haasteren.

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

// Also use cubic spline interpolation from the GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_cdf.h>


void d_swap(double &x, double &y) {
  double dTemp;
  dTemp = x;
  x = y;
  y = dTemp;
} // d_swap()


void QuickSort(double *pdBuffer, int left, int right) {
//  left is the lower index, right is the upper index
//  of the region of array a that is to be sorted
    int i=left, j=right;
    double pivot = pdBuffer[(left+right)/2];

    //  partition
    while(i <= j) {    
        while (pdBuffer[i]<pivot) i++; 
        while (pdBuffer[j]>pivot) j--;
        if (i<=j)
        {
            d_swap(pdBuffer[i], pdBuffer[j]);
            i++; j--;
        }
    }

    //  recursion
    if (left<j) QuickSort(pdBuffer, left, j);
    if (i<right) QuickSort(pdBuffer, i, right);
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
  FILE *pFileIn, *pFileOut;
  char *strFileIn, *strFileOut, *strFileOutML,
       *strFileOutMLCorr, *strFileOutGP, buf[80];
  double *pdPlotX, *pdPlotY, *pdPlotYML, *pdPlotYMLCorr;
  double *pdPlotXInt, *pdPlotYInt, *pdPlotYMLInt, *pdPlotYMLCorrInt;
  double *pdPercX, **ppdPercY, **ppdPercYML, **ppdPercYMLCorr;
  double *pdSigmaL, *pdSigmaH, *pdPdfSort, *pdPdf;
  double *pdWidth, *pdWidthML, *pdWidthMLCorr;
  double *pdTrueX;
  int nDataSets, nParameters, nPlotPoints, nIntPoints;
  int nTemp=0, nTrueX, nIndex;
  double dTrueLL, dTrueLLML, dTrueLLMLCorr, dArea, dAreaML,
	 dAreaMLCorr, dCumArea, dCumAreaML, dCumAreaMLCorr, dLevel;
  gsl_interp_accel *pAccel;
  gsl_spline *pSpline;

  bool bExponentiate=false, bBoundExp=false, bNoTrueX=false,
       bUnitize=false, bRepeat=false, bNoML=false;

  for(int i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      // process the switches
      for(int j=1; j<strlen(argv[i]); j++) {
        switch(argv[i][j]) {
        case 'm':
          bNoML=true;
          break;
        case 'e':
          bExponentiate=true;
          break;
        case 'n':
          bNoTrueX=true;
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
      if(nTemp==1) strFileIn = argv[i];
      if(nTemp==2) strFileOut = argv[i];
      if(nTemp==3) strFileOutML = argv[i];
      if(nTemp==4) strFileOutMLCorr = argv[i];
      if(nTemp==5) strFileOutGP = argv[i];
    } // if argv
  } // for i

  if(nTemp < 3 || (nTemp < 5 && bNoML == false)) {
    printf("Usage is: processensembleresults1d [-benu] inputfile.dat outputfile.txt outputfileML.txt  outputfileMLCorr.txt gnuplotfile.gp\n");
    printf("\n");
    printf("                         -n  No true X in the data\n");
    printf("                         -m  No ML in the data\n");
    return 1;
  } // if nTemp

  if(nTemp == 3) {
    strFileOutGP = strFileOutML;
  } // if nTemp

  // Read the files:
  try {
    if(! (pFileIn = fopen(strFileIn, "rb")) ) throw 1;

    if(! fread(&nDataSets, sizeof(int), 1, pFileIn) ) throw 2;
    if(! fread(&nParameters, sizeof(int), 1, pFileIn) ) throw 3;
    if(! fread(&nPlotPoints, sizeof(int), 1, pFileIn) ) throw 4;

    nIntPoints = 10*nDataSets;
    pdPlotX = new double[nPlotPoints];
    pdPlotY = new double[nPlotPoints];
    pdPlotYML = new double[nPlotPoints];
    pdPlotYMLCorr = new double[nPlotPoints];
    pdPlotXInt = new double[nIntPoints];
    pdPlotYInt = new double[nIntPoints];
    pdPlotYMLInt = new double[nIntPoints];
    pdPlotYMLCorrInt = new double[nIntPoints];
    pdTrueX = new double[nParameters];

    pdSigmaL = new double[nDataSets+1];
    pdSigmaH = new double[nDataSets+1];
    pdPdfSort = new double[nDataSets+1];
    pdPdf = new double[nDataSets+1];
    pdPercX = new double[nDataSets+1];
    ppdPercY = new double*[nParameters];
    ppdPercYML = new double*[nParameters];
    ppdPercYMLCorr = new double*[nParameters];
    pdWidth = new double[nParameters];
    pdWidthML = new double[nParameters];
    pdWidthMLCorr = new double[nParameters];

    for(int p=0; p<nParameters; p++) {
      ppdPercY[p] = new double[nDataSets+1];
      ppdPercYML[p] = new double[nDataSets+1];
      ppdPercYMLCorr[p] = new double[nDataSets+1];

      pdWidth[p] = 0;
      pdWidthML[p] = 0;
      pdWidthMLCorr[p] = 0;
    } // for p

    if(! bNoTrueX) {
      if(! fread(pdTrueX, sizeof(double), nParameters, pFileIn) ) throw 4;
    } else {
      printf("Assuming true values aren't there. Setting zero...\n");
      for(int p=0; p<nParameters; p++) {
	pdTrueX[p] = 0;
      } // for p
    } // if bNoTrueX

    // Allocate the accelerator and the spline object
    pAccel = gsl_interp_accel_alloc();
    pSpline = gsl_spline_alloc(gsl_interp_cspline, nPlotPoints);

    // Set the x-axis percentages, and the cdf sigma percentages
    for(int i=0; i<nDataSets+1; i++) {
      // Calculate the binomial pdf
      for(int j=0; j<nDataSets+1; j++) {
	if(j==0) {
	  pdPdf[j] = gsl_cdf_binomial_P(j, double(i)/double(nDataSets), nDataSets);
	} else {
	  pdPdf[j] = gsl_cdf_binomial_P(j, double(i)/double(nDataSets), nDataSets) - gsl_cdf_binomial_P(j-1, double(i)/double(nDataSets), nDataSets);
	} // if j
	pdPdfSort[j] = pdPdf[j];
      } // for j
      QuickSort(pdPdfSort, 0, nDataSets);
      // Find the 1-sigma line. 1, 2, 3 sigma:
      // 0.68268949   0.95449974   0.99730024
      dCumArea = 0; dArea = 0;
      for(int j=nDataSets; j>=0; j--) {
	dCumArea += pdPdfSort[j];
//	if(dCumArea >= 0.68268949) {
	if(dCumArea >= 0.95449974) {
	  dArea = pdPdfSort[j];
	  break;
	} // if dCumArea
      } // for j

      pdSigmaL[i] = 0;
      for(int j=0; j<nDataSets+1; j++) {
	if(pdPdf[j] > dArea) {
	  pdSigmaL[i] = j * 100.0 / nDataSets;
	  break;
	} // if pdPdfSort
      } // for j

      pdSigmaH[i] = 0;
      for(int j=nDataSets; j>=0; j--) {
	if(pdPdf[j] > dArea) {
	  pdSigmaH[i] = j * 100.0 / nDataSets;
	  break;
	} // if pdPdfSort
      } // for j
      if(i==nDataSets) {
	pdSigmaL[i] = 100;
	pdSigmaH[i] = 100;
      } // if i

      pdPercX[i] = i * 100.0 / nDataSets;
    } // for k

    // Set all the y-axis percentages to zero
    for(int p=0; p<nParameters; p++) {
      for(int i=0; i<nDataSets; i++) {
	ppdPercY[p][i] = 0;
	ppdPercYML[p][i] = 0;
	ppdPercYMLCorr[p][i] = 0;
      } // for k
    } // for p

    for(int i=0; i<nDataSets; i++) {
      for(int p=0; p<nParameters; p++) {
	if(! fread(pdPlotX, sizeof(double), nPlotPoints, pFileIn) ) throw 5;
	if(! fread(pdPlotY, sizeof(double), nPlotPoints, pFileIn) ) throw 6;
	if(! bNoML) {
	  if(! fread(pdPlotYML, sizeof(double), nPlotPoints, pFileIn) ) throw 7;
	  if(! fread(pdPlotYMLCorr, sizeof(double), nPlotPoints, pFileIn) ) throw 7;
	} // if bNoML

	// Initialise the spline interpolator for pdPlotY
	gsl_spline_init(pSpline, pdPlotX, pdPlotY, nPlotPoints);

	// Interpolate the pdPlotY curve on the following values
	for(int k=0; k<nIntPoints; k++) {
	  pdPlotXInt[k] = pdPlotX[0] + k*(pdPlotX[nPlotPoints-1]-pdPlotX[0])/(nIntPoints-1);
	  if(k==nIntPoints-1) pdPlotXInt[k] = pdPlotX[nPlotPoints-1];
	  pdPlotYInt[k] = gsl_spline_eval(pSpline, pdPlotXInt[k], pAccel);
	} // for k


	if(! bNoML) {
	  // Initialise the spline interpolator for pdPlotYML
	  gsl_spline_init(pSpline, pdPlotX, pdPlotYML, nPlotPoints);

	  // Interpolate the pdPlotYML curve on the following values
	  for(int k=0; k<nIntPoints; k++) {
	    pdPlotXInt[k] = pdPlotX[0] + k*(pdPlotX[nPlotPoints-1]-pdPlotX[0])/(nIntPoints-1);
	    if(k==nIntPoints-1) pdPlotXInt[k] = pdPlotX[nPlotPoints-1];
	    pdPlotYMLInt[k] = gsl_spline_eval(pSpline, pdPlotXInt[k], pAccel);
	  } // for k


	  // Initialise the spline interpolator for pdPlotYMLCorr
	  gsl_spline_init(pSpline, pdPlotX, pdPlotYMLCorr, nPlotPoints);

	  // Interpolate the pdPlotYMLCorr curve on the following values
	  for(int k=0; k<nIntPoints; k++) {
	    pdPlotXInt[k] = pdPlotX[0] + k*(pdPlotX[nPlotPoints-1]-pdPlotX[0])/(nIntPoints-1);
	    if(k==nIntPoints-1) pdPlotXInt[k] = pdPlotX[nPlotPoints-1];
	    pdPlotYMLCorrInt[k] = gsl_spline_eval(pSpline, pdPlotXInt[k], pAccel);
	  } // for k

#if 0
	  if(i==0 && p == 0) {
	    double dMax, dMaxML, dMaxMLCorr;
	    dMax = pdPlotYInt[0];
	    dMaxML = pdPlotYMLInt[0];
	    dMaxMLCorr = pdPlotYMLCorrInt[0];
	    for(int k=0; k<nIntPoints; k++) {
	      if(dMax < pdPlotYInt[k]) dMax = pdPlotYInt[k];
	      if(dMaxML < pdPlotYMLInt[k]) dMaxML = pdPlotYMLInt[k];
	      if(dMaxMLCorr < pdPlotYMLCorrInt[k]) dMaxMLCorr = pdPlotYMLCorrInt[k];
	    } // for k

	    if(! (pFileOut = fopen("test.txt", "w+")) ) throw 8;
	    for(int k=0; k<nIntPoints; k++) {
	      fprintf(pFileOut, "%e   %e   %e   %e\n",
		  pdPlotXInt[k], pdPlotYInt[k]/dMax,
		  pdPlotYMLInt[k]/dMaxML, pdPlotYMLCorrInt[k]/dMaxMLCorr);
	    } // for k
	    fclose(pFileOut);
	  } // if i
#endif
	} // if bNoMLCorr

	// Search for the entry that contains the true value x=0, and calculate
	// the area under the curves
	nTrueX = -1; dArea = 0; dAreaML = 0; dAreaMLCorr = 0;
	for(int k=0; k<nIntPoints; k++) {
	  dArea += pdPlotYInt[k];
	  if(! bNoML) {
	    dAreaML += pdPlotYMLInt[k];
	    dAreaMLCorr += pdPlotYMLCorrInt[k];
	  } // if bNoML
	  if(k < nIntPoints-1) {
	    if(pdPlotXInt[k] <= pdTrueX[p] && pdPlotXInt[k+1] > pdTrueX[p]) {
	      nTrueX = k;
	    } // if pdPlotXInt
	  } // if k
	} // for k

	if(nTrueX == -1) {
	  // This should be very very rare (only with wrong model or so)
	  if(pdPlotXInt[0] > pdTrueX[p]) {
	    nTrueX = 0;
	  } else {
	    nTrueX = nIntPoints - 1;
	  } // if pdPlotXInt
	} // if nTrueX
	dTrueLL = pdPlotYInt[nTrueX];
	if(! bNoML) {
	  dTrueLLML = pdPlotYMLInt[nTrueX];
	  dTrueLLMLCorr = pdPlotYMLCorrInt[nTrueX];
	} // if bNoML

	// Sort the two distributions small to big, and figure out the inner percentage
	dCumArea = 0; dCumAreaML = 0; dCumAreaMLCorr = 0;
	QuickSort(pdPlotYInt, 0, nIntPoints-1);
	QuickSort(pdPlotYMLInt, 0, nIntPoints-1);
	QuickSort(pdPlotYMLCorrInt, 0, nIntPoints-1);
	nIndex = nIntPoints - 1;
	while(nIndex >= 0 && pdPlotYInt[nIndex] >= dTrueLL) {
	  dCumArea += pdPlotYInt[nIndex];
	  nIndex--;
	} // while pdPlotYInt
	if(! bNoML) {
	  nIndex = nIntPoints - 1;
	  while(nIndex >= 0 && pdPlotYMLInt[nIndex] >= dTrueLLML) {
	    dCumAreaML += pdPlotYMLInt[nIndex];
	    nIndex--;
	  } // while pdPlotYML

	  nIndex = nIntPoints - 1;
	  while(nIndex >= 0 && pdPlotYMLCorrInt[nIndex] >= dTrueLLMLCorr) {
	    dCumAreaMLCorr += pdPlotYMLCorrInt[nIndex];
	    nIndex--;
	  } // while pdPlotYMLCorrInt
	} // if bNoML

	// Fill the percentage plot
	for(int ii=0; ii<nDataSets+1; ii++) {
	  if(dCumArea*100.0 <= dArea*pdPercX[ii]) {
	    ppdPercY[p][ii] += 1.0;
	  } // if dCumArea
	  if(! bNoML) {
	    if(dCumAreaML*100.0 <= dAreaML*pdPercX[ii]) {
	      ppdPercYML[p][ii] += 1.0;
	    } // if dCumArea
	    if(dCumAreaMLCorr*100.0 <= dAreaMLCorr*pdPercX[ii]) {
	      ppdPercYMLCorr[p][ii] += 1.0;
	    } // if dCumArea
	  } // if bNoML
	} // for ii

	if(! bNoML) {
	  // Also find the width of the plots
	  nIndex = nIntPoints - 1;
	  dCumArea = 0;
	  while(nIndex >= 0 && dCumArea <= 0.68268949 * dArea) {
	    dCumArea += pdPlotYInt[nIndex];
	    nIndex--;
	  } // while pdPlotYML
	  pdWidth[p] += nIntPoints - nIndex - 1;

	  nIndex = nIntPoints - 1;
	  dCumArea = 0;
	  while(nIndex >= 0 && dCumArea <= 0.68268949 * dAreaML) {
	    dCumArea += pdPlotYMLInt[nIndex];
	    nIndex--;
	  } // while pdPlotYML
	  pdWidthML[p] += nIntPoints - nIndex - 1;

	  nIndex = nIntPoints - 1;
	  dCumArea = 0;
	  while(nIndex >= 0 && dCumArea <= 0.68268949 * dAreaMLCorr) {
	    dCumArea += pdPlotYMLCorrInt[nIndex];
	    nIndex--;
	  } // while pdPlotYML
	  pdWidthMLCorr[p] += nIntPoints - nIndex - 1;
	} // if bNoML
      } // for p
    } // for i
    fclose(pFileIn);


    if(! (pFileOut = fopen(strFileOut, "w+")) ) throw 8;
    for(int i=0; i<nDataSets+1; i++) {
      fprintf(pFileOut, "%7.3f   %7.3f   %7.3f", pdPercX[i], pdSigmaL[i], pdSigmaH[i]);
      for(int p=0; p<nParameters; p++) {
	fprintf(pFileOut, "   %7.3f", ppdPercY[p][i]*100.0 / nDataSets);
      } // for p
      fprintf(pFileOut, "\n");
    } // for i
    fclose(pFileOut);

    if(! bNoML) {
      // The ML plot
      if(! (pFileOut = fopen(strFileOutML, "w+")) ) throw 8;
      for(int i=0; i<nDataSets+1; i++) {
	fprintf(pFileOut, "%7.3f   %7.3f   %7.3f", pdPercX[i], pdSigmaL[i], pdSigmaH[i]);
	for(int p=0; p<nParameters; p++) {
	  fprintf(pFileOut, "   %7.3f", ppdPercYML[p][i]*100.0 / nDataSets);
	} // for p
	fprintf(pFileOut, "\n");
      } // for i
      fclose(pFileOut);

      // The MLCorr plot
      if(! (pFileOut = fopen(strFileOutMLCorr, "w+")) ) throw 8;
      for(int i=0; i<nDataSets+1; i++) {
	fprintf(pFileOut, "%7.3f   %7.3f   %7.3f", pdPercX[i], pdSigmaL[i], pdSigmaH[i]);
	for(int p=0; p<nParameters; p++) {
	  fprintf(pFileOut, "   %7.3f", ppdPercYMLCorr[p][i]*100.0 / nDataSets);
	} // for p
	fprintf(pFileOut, "\n");
      } // for i
      fclose(pFileOut);

      // The width comparison table
      fprintf(stderr, "Width of the parameters:\n");
      fprintf(stderr, "                  MCMC      ML     MLCorr\n");
      for(int p=0; p<nParameters; p++) {
	fprintf(stderr, "Parameter [%2i]: %6.1f   %6.1f   %6.1f\n", p,
	    pdWidth[p]*40/(nIntPoints*nDataSets),
	    pdWidthML[p]*40/(nIntPoints*nDataSets),
	    pdWidthMLCorr[p]*40/(nIntPoints*nDataSets));
      } // for p
    } // if bNoML

    // Check if we have a file with parameters description available
    pFileIn = NULL;
    if(! bNoML) {
      if(! (pFileIn = fopen("parameters.txt", "r"))) {
	pFileIn = NULL;
      } // if pFileIn
    } // if bNoML

    // For convenience, also write a gnuplot file
    if(! (pFileOut = fopen(strFileOutGP, "w+")) ) throw 8;
    fprintf(pFileOut, "set terminal postscript eps enhanced monochrome\n");
    fprintf(pFileOut, "set output \"edfbayes.eps\"\n\n");

    fprintf(pFileOut, "set title \"Empirical distribution functions for Bayesian analysis\"\n");
    fprintf(pFileOut, "set xlabel \"a\"\n");
    fprintf(pFileOut, "set ylabel \"F_{i,1000}(a)\"\n");
    fprintf(pFileOut, "unset zlabel\n\n");

    fprintf(pFileOut, "set key at 0.18, 0.95\n\n");

    fprintf(pFileOut, "set xrange [0:1]\n");
    fprintf(pFileOut, "set yrange [0:1]\n");

    fprintf(pFileOut, "\n");

    fprintf(pFileOut, "set style line  1 lt  1 lw   5\n");
    fprintf(pFileOut, "set style line  2 lt  2 lw   2\n");
    fprintf(pFileOut, "set style line  3 lt  5 lw   2\n");
    fprintf(pFileOut, "set style line  4 lt 10 lw   2\n");
    fprintf(pFileOut, "set style line  5 lt  3 lw 0.5\n");
    fprintf(pFileOut, "set style line  6 lt  6 lw 0.2\n");
    fprintf(pFileOut, "set style line  7 lt  7 lw 0.2\n");
    fprintf(pFileOut, "set style line  8 lt  8 lw 0.2\n");
    fprintf(pFileOut, "set style line  9 lt  9 lw 0.2\n");
    fprintf(pFileOut, "set style line 10 lt  4 lw 0.2\n");
    fprintf(pFileOut, "set style line 11 lt 11 lw 0.2\n");
    fprintf(pFileOut, "set style line 12 lt 12 lw 0.2\n");

    fprintf(pFileOut, "\n");

    fprintf(pFileOut, "plot x+0.052 ti \"K-S b_{+}\" w l ls 1, \\\n", strFileOut);
    fprintf(pFileOut, "     x-0.052 ti \"K-S b_{+}\" w l ls 1, \\\n", strFileOut);



//    fprintf(pFileOut, "     \"%s\" u ($1/100):($5/100) ti \"2 {/Symbol s}\" w l lw 3, \\\n", strFileOut);
    for(int p=0; p<nParameters; p++) {
      if(pFileIn && getdataline(pFileIn, buf, 80)) {
	fprintf(pFileOut, "     \"%s\" u ($1/100):($%i/100) ti \"%s\" w l ls %i", strFileOut,
	    p+4, buf, (p%8)+5);
      } else {
	fprintf(pFileOut, "     \"%s\" u ($1/100):($%i/100) ti \"par %i\" w l ls %i", strFileOut,
	    p+4, p+1, (p%8)+5);
      } // if pFileIn

      if(p != nParameters - 1) {
	fprintf(pFileOut, ", \\\n");
      } else {
	fprintf(pFileOut, "\n");
      } // if p
    } // for p
    fprintf(pFileOut, "\n");
    fclose(pFileOut);
    if(pFileIn) fclose(pFileIn);

    // Free the accelerator and the spline object
    gsl_spline_free(pSpline);
    gsl_interp_accel_free(pAccel);

    // Temporary presentation stuff
    if(0) {
      if(! (pFileOut = fopen("temp.txt", "w+")) ) throw 8;
      for(int k=0; k<nIntPoints; k++) {
	pdPlotXInt[k]= k * 10.0 / nIntPoints;
	pdPlotYInt[k] = 1.5 * exp( -0.6*(3-pdPlotXInt[k])*(3-pdPlotXInt[k])) + exp( -0.8*(6-pdPlotXInt[k])*(6-pdPlotXInt[k]));
	fprintf(pFileOut, "%e    %e\n", pdPlotXInt[k], pdPlotYInt[k]);
      } // for k
      fclose(pFileOut);
    }

    for(int p=0; p<nParameters; p++) {
      delete[] ppdPercY[p];
      delete[] ppdPercYML[p];
      delete[] ppdPercYMLCorr[p];
    } // for p
    delete[] pdPdfSort;
    delete[] pdPdf;
    delete[] pdSigmaL;
    delete[] pdSigmaH;
    delete[] pdPercX;
    delete[] ppdPercY;
    delete[] ppdPercYML;
    delete[] ppdPercYMLCorr;
    delete[] pdPlotX;
    delete[] pdPlotY;
    delete[] pdPlotYML;
    delete[] pdPlotYMLCorr;
    delete[] pdPlotXInt;
    delete[] pdPlotYInt;
    delete[] pdPlotYMLInt;
    delete[] pdPlotYMLCorrInt;
    delete[] pdTrueX;
    delete[] pdWidth;
    delete[] pdWidthML;
    delete[] pdWidthMLCorr;
  } catch(int nError) {
    printf("Error number: %i \n", nError);
  } // try

  return 0;
}
