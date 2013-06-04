/* analyze.cpp -- Analysis cpp-file for PTA bayesian analysis program
 *
 * Rutger van Haasteren 15 August 2007 haasteren@strw.leidenuniv.nl
 *
 * Copyright (C) 2006-2010 Rutger van Haasteren.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 * */


#include "config.h"

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_cdf.h>
// #include "moremath.h"
#include "linal.h"
#include "linalfunc.h"
#include "banfunc.h"
#include "corefunctions.h"
#include "analyze.h"
#include "filefunctions.h"



using namespace std;

// TODO: Apply the FITPARAMETER constants to these functions


/* Performs a MCMC integration with a Gaussian proposal distribution. The
 * final chain is written to a datafile: strFileName.
 *
 * Todo: during the burnin-time, the random walkers do walk, but the proposal
 * distribution width is not adjusted. The process of setting the width of
 * the distribution to yield an optimal acceptance rate must be automized. */
void MCMCGauss(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName) {
// dWidthExp dWidthAmp
// Generating the data was:
// vdTemp = oData.mdEVC*(Sqrt(oData.vdLambdaC*2)&&InvErf(vdRand*2-1));
//
// Generate new point with: current + delta = Current + Sqrt(2)*vdSigma&&InvErf(vdRand*2-1);

  int nMCMCSteps, nPlotPoints, nStep, nStepIndex;	// Total steps, number of plotting points on a graph, and step number
  int nAcceptedSteps;
  bool bBurnIn=true;				// Are we in the burn-in period?
  int nCheckSteps=0;				// Amount of steps, accepted or not accepted, since last accept-check
  int nIndex, nTemp;
  int nParCount;
  time_t tStartBurnin, tStart, tTime;		// Time stamps
  double dAcceptRatio;				// Likelihood-ratio
  char strBuf[100];
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);	// Random number generator
  gsl_rng_set(rng, (unsigned int)time(NULL));

  MCMCDataPoint *pDat;				// All data aqcuired
  CNumber ndRand1(0), ndRand2(0);

  nPlotPoints = oConstants.nPlotPoints;
  nMCMCSteps = oConstants.nMCMCSteps;
  pDat = new MCMCDataPoint[MAX_MCMC_BUFFER];

//  InitProgressBar("Generating MCMC points...");
//  DrawProgressBar(0);
  fprintf(stderr, "Generating MCMC points...\n");

  tStart = clock();

  // First set some values for the first random selection point
  ndRand1.Randomize();
  StartParametersFromSources(pDat[0].oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromStartParameters(pDat[0].oParameters, oConstants.poSources, oConstants.nSources);
  pDat[0].dLogLik = LogLikelihoodTimesPrior(oData, pDat[0].oParameters, oConstants);
  if(oConstants.bCalcTMPars) {
    pDat[0].mdCXiInv = oData.mdCXiInv;
    pDat[0].vdChi = oData.vdChi;
  } // if bCalcTMPars

#if 0  // Print the likelihood value
  printf("dLogLik = %e\n", pDat[0].dLogLik);
#endif

  // Write the header of the MCMC parameter file
  WriteMCMCDataFileStart(oConstants.nParameters, oConstants.nMarParameters, strFileName);

  // Now repeat this process many times to get a chain
  nStep=1;
  nStepIndex=1;
  nAcceptedSteps=1;
  nCheckSteps=1;
  while(nStep < nMCMCSteps) {
    pDat[nStepIndex].oParameters = pDat[0].oParameters;
    // Generate new parameters
    for(int s=0; s<oConstants.nSources; s++) {
      for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
	nIndex = oConstants.poSources[s].nFirstParIndex + p;
	if(SourceWorks(oConstants, s) && (oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p))) {
	  for(;;) { // Keep checking for a good value
	    pDat[nStepIndex].oParameters.pdPar[nIndex] = pDat[nStepIndex-1].oParameters.pdPar[nIndex] +
	      pDat[nStepIndex-1].oParameters.pdParWidthMCMC[nIndex] * 
	      oConstants.dGlobalMCMCWidthFactor *
	      gsl_ran_ugaussian(rng);
	    if(oConstants.poSources[s].oSourceType.nWrapTag & FITPARAMETER(p)) {
	      // Check whether we need to wrap this parameter
	      if(pDat[nStepIndex].oParameters.pdPar[nIndex] < pDat[nStepIndex].oParameters.pdParMinBound[nIndex]) {
		pDat[nStepIndex].oParameters.pdPar[nIndex] = pDat[nStepIndex].oParameters.pdParMaxBound[nIndex] - (pDat[nStepIndex].oParameters.pdParMinBound[nIndex] - pDat[nStepIndex].oParameters.pdPar[nIndex]);
	      } // if pdPar < pdParMinBound

	      if(pDat[nStepIndex].oParameters.pdPar[nIndex] > pDat[nStepIndex].oParameters.pdParMaxBound[nIndex]) {
		pDat[nStepIndex].oParameters.pdPar[nIndex] = pDat[nStepIndex].oParameters.pdParMinBound[nIndex] + (pDat[nStepIndex].oParameters.pdPar[nIndex] - pDat[nStepIndex].oParameters.pdParMaxBound[nIndex]);
	      } // if pdPar < pdParMinBound
	    } // if nWrapTag

	    if(pDat[nStepIndex].oParameters.pdPar[nIndex] >= pDat[nStepIndex].oParameters.pdParMinBound[nIndex] &&
		pDat[nStepIndex].oParameters.pdPar[nIndex] <= pDat[nStepIndex].oParameters.pdParMaxBound[nIndex])
	      break;
	  } // for
	} // if nFitTag
      } // for p
    } // for s

    try{
      pDat[nStepIndex].dLogLik = LogLikelihoodTimesPrior(oData, pDat[nStepIndex].oParameters, oConstants);
      if(oConstants.bCalcTMPars) {
	pDat[nStepIndex].mdCXiInv = oData.mdCXiInv;
	pDat[nStepIndex].vdChi = oData.vdChi;
      } // if bCalcTMPars
#if 0  // Print the likelihood value
      printf("\ndLogLik = %e\n", pDat[nStepIndex].dLogLik);
#endif
    } catch(ELinearError err) {  // Error handling
      // TODO: Put something useful here
      pDat[nStepIndex].dLogLik = 10E200;
    } // try

    dAcceptRatio = exp(pDat[nStepIndex-1].dLogLik - pDat[nStepIndex].dLogLik);

    if(gsl_rng_uniform(rng) < dAcceptRatio) {
      nAcceptedSteps++;
    } else {
      // Reject the point
      for(int s=0; s<oConstants.nSources; s++) {
	for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
	  nIndex = oConstants.poSources[s].nFirstParIndex + p;
	  if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	      pDat[nStepIndex].oParameters.pdPar[nIndex] = pDat[nStepIndex-1].oParameters.pdPar[nIndex];
	  } // if nFitTag
	} // for p
      } // for s
      pDat[nStepIndex].dLogLik = pDat[nStepIndex-1].dLogLik;
    } // if gsl_rng_uniform

    // We except this point, but first check whether we are still
    // in the burn-in period.
    if(bBurnIn && nStep >= oConstants.nBurnInSteps) {
      pDat[0].oParameters = pDat[nStepIndex].oParameters;
      pDat[0].dLogLik = pDat[nStepIndex].dLogLik;
      bBurnIn = false;
      nStep = 0; nAcceptedSteps = 0;;
      tStartBurnin = tStart;
      tStart = clock();
    } // bBurnIn

    nStep++; nStepIndex++;

    if(nStepIndex == MAX_MCMC_BUFFER) {
      // We need to write this stuff to a file
      if(! bBurnIn)
	WriteMCMCDataFileAppend(MAX_MCMC_BUFFER, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName, oConstants.bCalcTMPars);
      pDat[0] = pDat[nStepIndex-1];
      nStepIndex = 1;
    } // if nStepIndex

#if 0
    if(bBurnIn) {
      if(nCheckSteps % 1000 == 0) { // This is the accept-check in the burn-in period
        // TODO: Adjust the sample distribution width accordingly wrt the acceptance rate
        nAcceptedSteps = 0; nCheckSteps = 0;
      } // nCheckSteps
    } // bBurnIn
    nCheckSteps++;
#endif

//    if(nStep % 20 == 0) {
      tTime = clock();
//      DrawProgressBar(int(100.0*nStep/nMCMCSteps));
// Show some nice progress indicator:
      if(nStepIndex >= 2) {
	sprintf(strBuf, "%3i, %2i: ", int(100.0*nStep/nMCMCSteps), int(nAcceptedSteps*100.0/nStep));

	if(bBurnIn)
	  sprintf(strBuf, "B-%s", strBuf);
	else
	  sprintf(strBuf, "S-%s", strBuf);

	nParCount = 0;
	for(int s=0; s<oConstants.nSources && nParCount < 5; s++) {
	  for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters && nParCount < 5; p++) {
	    if(SourceWorks(oConstants, s) && (oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p))) {
	      nParCount++;
	      sprintf(strBuf, "%s%6.2e | ", strBuf, pDat[nStepIndex-2].oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+p]);
	    } // if oConstants
	  } // for p
	} // for s

      fprintf(stderr, "%s\r", strBuf);
#if 0  // Check what the values of _all_ parameters are
      printf("\n");
	for(int s=0; s<oConstants.nSources; s++) {
	  for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
	    printf("%6.2e | ", pDat[nStepIndex-2].oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+p]);
//	    printf("%6.2e | ", pDat[nStepIndex-2].oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+p]);
	  }
	}
      printf("\n");
#endif
    } // if nStep
  } // while nStep

  tTime = clock();
//  FinishProgressBar();
  fprintf(stderr, "\nComputing time: %lf sec.\n", double(tTime - tStartBurnin)/CLOCKS_PER_SEC); 

  // Write the data to a datafile
//  WriteMCMCData(nMCMCSteps, oConstants.nParameters, pDat, strFileName);
  WriteMCMCDataFileAppend(nStepIndex-1, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName, oConstants.bCalcTMPars);

  // This should call for the destructors...
  delete[] pDat;
  gsl_rng_free(rng);
} // MCMCGauss


/* The state of the cloud of samples
 *
 * The first index is the number of the sample in the cloud. The second index is
 * the number of the parameter.
 *
 * ppdX[i][k], sample i, parameter k
 * pdLL[i], loglikelihood of sample i
 * */
struct SSampleState {
  double **ppdX;
  double **ppdCXiInv;
  double **ppdChi;;
  double *pdLL;
};


/* This function returns a scalar, drawn from a 1/sqrt(z) distribution, between
 * the boundaries 1/sqrt(a) and sqrt(a)
 * */
double DrawEnsembleScalar(gsl_rng *rng, double sqrta) {
  return gsl_pow_2(
      		(gsl_rng_uniform(rng)*(sqrta*sqrta-1.0)/sqrta)
		+ 1.0/sqrta
	);
} //  DrawEnsembleScalar


/* This function generates a new set of samples in parameter space, all drawn
 * using the ensemble sampling algorithm
 * */
void GenerateNewCloud(SSampleState &oOldState, SSampleState &oNewState, int nEnsembleSamples, int nParameterSpaceDimensions, int *pnParameterIndex, int pnAccepted[], gsl_rng *rng, SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dZ, dLogAcceptRatio;
  int nOtherSample;
  bool bInsideDomain=true;
  SParametersType oTemp;

  oTemp = oParameters;

  // For each sample, generate a new one
  for(int k=0; k<nEnsembleSamples; k++) {
    // Draw a sample from the 1/sqrt(z) distribution (this is where the
    // acceptance ratio is set) (IMPORTANT)
//    dZ = DrawEnsembleScalar(rng, 2);
    dZ = DrawEnsembleScalar(rng, oConstants.dEnsembleScalar);

    // Select one of the nEnsembleSamples-1 other samples
    nOtherSample = int(gsl_rng_uniform(rng)*(nEnsembleSamples-1));
    if(nOtherSample >= k) nOtherSample++;

    // Calculate the position of the new sample
    for(int p=0; p<nParameterSpaceDimensions; p++) {
      oNewState.ppdX[k][p] = oOldState.ppdX[nOtherSample][p] +
	dZ*(oOldState.ppdX[k][p] - oOldState.ppdX[nOtherSample][p]);
      oTemp.pdPar[pnParameterIndex[p]] = oNewState.ppdX[k][p];
    } // for p

    // Figure out if this value is inside our prior domain
    bInsideDomain=true;
    for(int p=0; p<nParameterSpaceDimensions; p++) {
      if(oNewState.ppdX[k][p] < oParameters.pdParMinBound[pnParameterIndex[p]] ||
	  oNewState.ppdX[k][p] > oParameters.pdParMaxBound[pnParameterIndex[p]]) {
	bInsideDomain=false;
      } // if ppdX
    } // for p

    if(bInsideDomain) {
      // We have a new sample, now calculate the Likelihood & acceptance ratio
      oNewState.pdLL[k] = -LogLikelihoodTimesPrior(oData, oTemp, oConstants);
      dLogAcceptRatio = (nParameterSpaceDimensions-1)*log(dZ) + oNewState.pdLL[k] - oOldState.pdLL[k];
//      fprintf(stderr, "n-1 = %i  log(Z) = %e  Z = %e  newll = %e  oldll = %e  dlogaccr = %e  accr = %e\n",
//	  nParameterSpaceDimensions-1, log(dZ), dZ, oNewState.pdLL[k], oOldState.pdLL[k], dLogAcceptRatio, exp(dLogAcceptRatio));
    } // if bInsideDomain

    // Now accept or reject the sample
    if(bInsideDomain && (log(gsl_rng_uniform(rng)) < dLogAcceptRatio) ) {
      // Accept
      pnAccepted[k]++;

      if(oConstants.bCalcTMPars) {
	for(int m=0; m<oConstants.nMarParameters; m++) {
	  oNewState.ppdChi[k][m] = oData.vdChi.m_pdData[m];
	} // for m
	for(int m=0; m<oConstants.nMarParameters*oConstants.nMarParameters; m++) {
	  oNewState.ppdCXiInv[k][m] = oData.mdCXiInv.m_pdData[m];
	} // for m
      } // if bCalcTMPars
    } else {
      // Reject
      for(int p=0; p<nParameterSpaceDimensions; p++) {
	oNewState.ppdX[k][p] = oOldState.ppdX[k][p];
      } // for p
      oNewState.pdLL[k] = oOldState.pdLL[k];
      if(oConstants.bCalcTMPars) {
	for(int m=0; m<oConstants.nMarParameters; m++) {
	  oNewState.ppdChi[k][m] = oOldState.ppdChi[k][m];
	} // for m
	for(int m=0; m<oConstants.nMarParameters*oConstants.nMarParameters; m++) {
	  oNewState.ppdCXiInv[k][m] = oOldState.ppdCXiInv[k][m];
	} // for m
      } // if bCalcTMPars
    } // if rng
  } // for k

  return;
} // GenerateNewCloud



/*  This function samples from the likelihood function using an ensemble
 *  sampler.*/
void EnsembleSampler(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName) {
  int nCloudSteps, nCloudSize, nParameterSpaceDimensions, nBurnSteps;
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  SSampleState *poState;
  int *pnAccepted;
  int *pnParameterIndex;
  int nStepIndex, nPrevStepIndex;
  int nMaxBuffer;
  int nMarParameters;
  int nAllCloudAccepted;
  SParametersType oTemp;
  MCMCDataPoint *pDat;

  gsl_rng_set(rng, (unsigned int)time(NULL));

  // Figure out how big the parameter space is
  nParameterSpaceDimensions = NumberOfVaryingParameters(oConstants);

  // Figure out the dimensions of the calculation
  nCloudSize = CLOUD_SIZE_PER_DIM * nParameterSpaceDimensions;
  nCloudSteps = oConstants.nMCMCSteps / nCloudSize;
  nBurnSteps = oConstants.nBurnInSteps / nCloudSize;
  nMaxBuffer = int(MAX_MCMC_BUFFER / nCloudSize) + 2;
//  nMaxBuffer = 2;
  nMarParameters = oConstants.nMarParameters;

  // Produce the starting point for the parameters
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromStartParameters(oParameters, oConstants.poSources, oConstants.nSources);
  oTemp = oParameters;

  // Allocate memory
//  poState = new SSampleState[nCloudSteps];
  poState = new SSampleState[nMaxBuffer];
  pnAccepted = new int[nCloudSize];
  pnParameterIndex = new int[nParameterSpaceDimensions];
//  for(int i=0; i<nCloudSteps; i++) {
  for(int i=0; i<nMaxBuffer; i++) {
    poState[i].pdLL = new double[nCloudSize];
    poState[i].ppdX = new double*[nCloudSize];
    if(oConstants.bCalcTMPars) {
      poState[i].ppdCXiInv = new double*[nCloudSize];
      poState[i].ppdChi = new double*[nCloudSize];
    } // if bCalcTMPars
    for(int k=0; k<nCloudSize; k++) {
      poState[i].ppdX[k] = new double[nParameterSpaceDimensions];
      if(oConstants.bCalcTMPars) {
	poState[i].ppdCXiInv[k] = new double[nMarParameters*nMarParameters];
	poState[i].ppdChi[k] = new double[nMarParameters];
      } // if bCalcTMPars
    } // for k
  } // for i

  pDat = new MCMCDataPoint[nMaxBuffer * nCloudSize];

  // Set the array of indices that allows quick reference in parameter structs
  SetParameterIndexArray(oConstants, pnParameterIndex);

  // Initialise the first samples of the cloud
  for(int k=0; k<nCloudSize; k++) {
    pnAccepted[k] = 0;
    for(int p=0; p<nParameterSpaceDimensions; p++) {
      poState[0].ppdX[k][p] = GenerateValidProposal(
	  oParameters.pdPar[pnParameterIndex[p]],
	  oParameters.pdParWidthMCMC[pnParameterIndex[p]],
	  oParameters.pdParMinBound[pnParameterIndex[p]],
	  oParameters.pdParMaxBound[pnParameterIndex[p]],
	  rng, false);
      oTemp.pdPar[pnParameterIndex[p]] = poState[0].ppdX[k][p];
    } // for p

    poState[0].pdLL[k] = -LogLikelihoodTimesPrior(oData, oTemp, oConstants);
    if(oConstants.bCalcTMPars) {
      for(int m=0; m<nMarParameters; m++) {
	poState[0].ppdChi[k][m] = oData.vdChi.m_pdData[m];
      } /* for m */
      for(int m=0; m<nMarParameters*nMarParameters; m++) {
	poState[0].ppdCXiInv[k][m] = oData.mdCXiInv.m_pdData[m];
      } /* for m */
    } // if bCalcTMPars
  } // for k

  // Initialise the pDat structure
  for(int i=0; i<nMaxBuffer*nCloudSize; i++) {
    pDat[i].dLogLik = -poState[0].pdLL[nCloudSize-1];
    pDat[i].oParameters = oParameters;
    if(oConstants.bCalcTMPars) {
      pDat[i].mdCXiInv = oData.mdCXiInv;
      pDat[i].vdChi = oData.vdChi;
    } // if bCalcTMPars
  } // for i

  // Run for some burn-in period
  nStepIndex = 1; nPrevStepIndex = 0;
  for(int i=1; i<nBurnSteps; i++) {
    GenerateNewCloud(poState[nPrevStepIndex], poState[nStepIndex], nCloudSize, nParameterSpaceDimensions, pnParameterIndex, pnAccepted, rng, oData, oParameters, oConstants);
    nStepIndex++; nPrevStepIndex = nStepIndex - 1;

    // Print the acceptance ratio
    fprintf(stderr, "Burn-in %2i%%:  ", int(i * 100.0/nBurnSteps));
    nAllCloudAccepted = 0;
    for(int k=0; k<nCloudSize; ++k) {
      nAllCloudAccepted += pnAccepted[k];
    } // for k
    fprintf(stderr, "a.r. %2i ", int(nAllCloudAccepted*100.0/(i*nCloudSize)));
    fprintf(stderr, "  ||  ");
    for(int p=0; (p<nParameterSpaceDimensions && p<4); p++) {
      if(p != 0) fprintf(stderr, "| ");
      fprintf(stderr, "%6.2e ", poState[nPrevStepIndex].ppdX[0][p]);
    } // for p
    fprintf(stderr, "\r");

    if(nStepIndex == nMaxBuffer) {
      nStepIndex = 0;
    } // if nStepIndex
  } // for i
  fprintf(stderr, "\n");

  for(int k=0; k<nCloudSize; k++) {
    pnAccepted[k] = 0;
    for(int p=0; p<nParameterSpaceDimensions; p++) {
      poState[0].ppdX[k][p] = poState[nPrevStepIndex].ppdX[k][p];
    } // for p
    poState[0].pdLL[k] = poState[nPrevStepIndex].pdLL[k];
  } // for k

  // Start the file where we save the chain
  WriteMCMCDataFileStart(oConstants.nParameters, oConstants.nMarParameters, strFileName);

  // Run the actual chain
  nStepIndex = 1; nPrevStepIndex = 0;
  for(int i=1; i<nCloudSteps; i++) {
    GenerateNewCloud(poState[nPrevStepIndex], poState[nStepIndex], nCloudSize, nParameterSpaceDimensions, pnParameterIndex, pnAccepted, rng, oData, oParameters, oConstants);
    nStepIndex++; nPrevStepIndex = nStepIndex - 1;

    fprintf(stderr, "Sampling %2i%%:  ", int(i * 100.0/nCloudSteps));
    nAllCloudAccepted = 0;
    for(int k=0; k<nCloudSize; ++k) {
      nAllCloudAccepted += pnAccepted[k];
    } // for k
    fprintf(stderr, "a.r. %2i ", int(nAllCloudAccepted*100.0/(i*nCloudSize)));
    fprintf(stderr, "  ||  ");
    for(int p=0; (p<nParameterSpaceDimensions && p<4); p++) {
      if(p != 0) fprintf(stderr, "| ");
      fprintf(stderr, "%6.2e ", poState[nPrevStepIndex].ppdX[0][p]);
    } // for p
    fprintf(stderr, "\r");

    if(nStepIndex == nMaxBuffer) {
      // Write the buffer to disk
      // First copy all the necessary info to the MCMCData struct
      for(int j=0; j<nMaxBuffer; j++) {
	for(int k=0; k<nCloudSize; k++) {
	  pDat[j*nCloudSize+k].oParameters = oParameters;
	  pDat[j*nCloudSize+k].dLogLik = -poState[j].pdLL[k];
	  for(int p=0; p<nParameterSpaceDimensions; p++) {
	    pDat[j*nCloudSize+k].oParameters.pdPar[pnParameterIndex[p]] =
	      poState[j].ppdX[k][p];
	  } // for p
	  if(oConstants.bCalcTMPars) {
	    for(int m=0; m<nMarParameters; m++) {
	      pDat[j*nCloudSize+k].vdChi.m_pdData[m] = poState[j].ppdChi[k][m];
	    } // for m
	    for(int m=0; m<nMarParameters*nMarParameters; m++) {
	      pDat[j*nCloudSize+k].mdCXiInv.m_pdData[m] =  poState[j].ppdCXiInv[k][m];
	    } // for m
	  } // if bCalcTMPars
	} // for k
      } // for j
      WriteMCMCDataFileAppend(nMaxBuffer*nCloudSize, oConstants.nParameters, nMarParameters, pDat, strFileName, oConstants.bCalcTMPars);
      nStepIndex = 0;
    } // if nStepIndex
  } // for i
  fprintf(stderr, "\n");

  if(nStepIndex > 1) {
    for(int j=0; j<nStepIndex-1; j++) {
      for(int k=0; k<nCloudSize; k++) {
	pDat[j*nCloudSize+k].oParameters = oParameters;
	pDat[j*nCloudSize+k].dLogLik = -poState[j].pdLL[k];
	for(int p=0; p<nParameterSpaceDimensions; p++) {
	  pDat[j*nCloudSize+k].oParameters.pdPar[pnParameterIndex[p]] =
	    poState[j].ppdX[k][p];
	} // for p
	if(oConstants.bCalcTMPars) {
	  for(int m=0; m<nMarParameters; m++) {
	    pDat[j*nCloudSize+k].vdChi.m_pdData[m] =  poState[j].ppdChi[k][m];
	  } // for m
	  for(int m=0; m<nMarParameters*nMarParameters; m++) {
	    pDat[j*nCloudSize+k].mdCXiInv.m_pdData[m] =  poState[j].ppdCXiInv[k][m];
	  } // for m
	} // if bCalcTMPars
      } // for k
    } // for j

    WriteMCMCDataFileAppend((nStepIndex-1)*nCloudSize, oConstants.nParameters, nMarParameters, pDat, strFileName, oConstants.bCalcTMPars);
  } // if nStepIndex


  // Free all memory allocated in this function
  for(int i=0; i<nMaxBuffer; i++) {
    for(int j=0; j<nCloudSize; j++) {
      delete[] poState[i].ppdX[j];
      if(oConstants.bCalcTMPars) {
	delete[] poState[i].ppdCXiInv[j];
	delete[] poState[i].ppdChi[j];
      } // if bCalcTMPars
    } // for j
    delete[] poState[i].ppdX;
    if(oConstants.bCalcTMPars) {
      delete[] poState[i].ppdCXiInv;
      delete[] poState[i].ppdChi;
    } // if bCalcTMPars
    delete[] poState[i].pdLL;
  } // for i
  delete []pnParameterIndex;
  delete[] pnAccepted;
  delete[] poState;
  delete[] pDat;
  gsl_rng_free(rng);
} // EnsembleSampler


/* This function reads an ordinary mcmc datafile, and then evaluates the
 * likelihood for the whole ensemble
 * */
void MCMCEnsembleImportanceSample(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nSkipSamples) {
  int nMCMCSteps, nStepIndex;
  SParametersType oParameters;
  char strEnsembleFile[160];

  CMCMCDataPoint *pDatP;
  MCMCDataEnsemble *pDatE;

  nMCMCSteps = oConstants.nMCMCSteps;

  pDatE = new MCMCDataEnsemble[MAX_MCMC_BUFFER];
  // pDatP is automatically allocated by ReadMCMCDataFile
  // pDatP = new MCMCDataPoint[nMCMCSteps];
  if(! ReadMCMCDataFile(nMCMCSteps, oConstants.nParameters, oConstants.nMarParameters, pDatP, strFileName) ) {
    fprintf(stderr, "ERROR: Cannot read %s\n", strFileName);
    return;
  } // if ReadMCMCDataFile

  // Fill the parameters struct
  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

  if(nSkipSamples == 0) {
    strcpy(strEnsembleFile, "mcmcensembledata.dat");
    if(! WriteMCMCEnsembleDataFileStart(oData.nDataSets, oConstants.n, oConstants.nParameters, oConstants.nMarParameters, oData.mdDataSets.m_pdData, strEnsembleFile)) {
      fprintf(stderr, "ERROR: Cannot create \"%s\"\n", strEnsembleFile);
      fprintf(stderr, "       Check whether it exists already\n");
      delete[] pDatP;
      delete[] pDatE;
      return;
    } else {
      fprintf(stderr, "Creating \"%s\" for data\n", strEnsembleFile);
    } // if WriteMCMCEnsembleDataFileStart
  } else {
    strcpy(strEnsembleFile, "mcmcensemble-resume.dat");
    fprintf(stderr, "Appending to \"%s\"\n", strEnsembleFile);
  } // if nSkipSamples

  nStepIndex=0;
  // Now fill all Ensemble structs, and write them to file
  for(int i=nSkipSamples; i<nMCMCSteps; i++) {
//    fprintf(stderr, "\nnParameters: %i\n", oConstants.nParameters);
    for(int p=0; p<oConstants.nParameters; p++) {
      oParameters.pdPar[p] = pDatP[i].pdPar[p];
    } // for p

#if 1 // THIS IS ONLY FOR MLDR ALGORITHM!!! CHANGE FOR IMPORTANCE RESAMPLING
    oParameters.pdPar[14] = 0;
    pDatP[i].pdPar[14] = 0;
#endif

    // Calculate all the likelihoods
    LogLikelihoodEnsemble(oData, oParameters, oConstants);

#if 0  // Just some checks

    for(int p=14; p<oConstants.nParameters; ++p) {
      fprintf(stderr, "%2i: %e\n", p, oParameters.pdPar[p]);
    } // for p
    for(int d=2; d<5 /*oData.nDataSets*/; ++d) {
      for(int j=0; j<oConstants.n; ++j) {
	oData.vdData[i] = double(oData.mdDataSets[i][d]);
//	fprintf(stderr, "%e    %e    %e\n", double(oData.vdData[j]), double(oData.mdDataSets[j][2]),  double(oData.vdData[j]) - double(oData.mdDataSets[j][2]));
      } // for i
      fprintf(stderr, "Kernel: %e   Set %2i: %e   control: %e   w prior: %e\n",
	  pDatP[i].dLogLik,
	  d,
	  double(oData.vdLL[d]),
	  LogLikelihood(oData, oParameters, oConstants),
	  LogLikelihoodTimesPrior(oData, oParameters, oConstants));
    } // for d
#endif

    for(int p=0; p<oConstants.nParameters; p++) {
      pDatE[nStepIndex].oParameters.pdPar[p] = pDatP[i].pdPar[p];
    } // for p
    pDatE[nStepIndex].mdCXiInv = oData.mdCXiInv;
    pDatE[nStepIndex].mdChi = oData.mdChi;
    pDatE[nStepIndex].vdLogLik = oData.vdLL;
    pDatE[nStepIndex].dLogLik = pDatP[i].dLogLik;

    nStepIndex++;
    if(i % 100 == 0)
      fprintf(stderr, "Ensemble: %2i%%\r", int(i*100.0/nMCMCSteps));

    if(nStepIndex == MAX_MCMC_BUFFER) {
      WriteMCMCEnsembleDataFileAppend(MAX_MCMC_BUFFER, oData.nDataSets, oConstants.nParameters, oConstants.nMarParameters, pDatE, strEnsembleFile);
      nStepIndex = 0;
    } // if nStepIndex
  } // for i

  if(nStepIndex > 1) {
    WriteMCMCEnsembleDataFileAppend(nStepIndex-1, oData.nDataSets, oConstants.nParameters, oConstants.nMarParameters, pDatE, strEnsembleFile);
  } // if nStepIndex

  // This should call for the destructors...
  delete[] pDatP;
  delete[] pDatE;
} // MCMCEnsembleImportanceSample


/* This function performs the actual MCMC integration and returns the result in
 * the vdX, vdY and vdYErr vectors
 * */
void Calculate1DMCMCIntegrationOnTheFly(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nPlotParameterNumber, CVector &vdX, CVector &vdY) {
  int nMCMCSteps, nPlotPoints, nRandIndex;
  int nFileParameters, nFileTMParameters;
  unsigned long nFileLength;
  FILE *pFile;
  char strMsg[160], strFileVersion[16];
  CNumber ndRand1;
  CMCMCDataPoint oDat;
  CVector vdPlot, vdFunctionValues, vdTemp, vdTemp2, vdMeanValues, vdMeanSquareValues, vdSigma, vdRawValues;
  int nIndex;
  double dPrevious;

  nPlotPoints = oConstants.nSmallPlotPoints;
  nMCMCSteps = oConstants.nMCMCSteps;

  strcpy(strMsg, "Making MCMC plot from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "' on the fly...\n");

  fprintf(stderr, strMsg);
  InitProgressBar("Integrating");
  DrawProgressBar(0);

  try {
    // Open the file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(nFileParameters != oConstants.nParameters) {
      printf("WARNING: nFileParameters != nParameters\n");
      printf("nFileParameters: %i,  nParameters: %i\n", nFileParameters, oConstants.nParameters);
    } // if nFileParameters
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;
    if(nFileTMParameters != oConstants.nMarParameters) {
      printf("WARNING: nFileTMParameters != TMnParameters\n");
      printf("nFileTMParameters: %i,  nTMParameters: %i\n", nFileTMParameters, oConstants.nMarParameters);
    } // if nFileTMParameters

    // Calculate nMCMCSteps
    if(oConstants.bCalcTMPars) {
      nMCMCSteps = (nFileLength-sizeof(char)*16 - 2*sizeof(int)) /
	(
	  sizeof(double)*(nFileParameters+1+nFileTMParameters*(nFileTMParameters+1))
	);
      if(nMCMCSteps <= 0) throw 7;
    } else {
      nMCMCSteps = (nFileLength-sizeof(char)*16 - 2*sizeof(int)) /
	(
	  sizeof(double)*(nFileParameters+1)
	);
    } // if bCalcTMPars

    // Leave the file open 'till the end
    oDat.Initialize(nFileParameters, nFileTMParameters);


    // Prepare for integrating: initialize the variables
    vdPlot.Initialize(nPlotPoints + 1);
    vdRawValues.Initialize(nMCMCSteps);
    vdFunctionValues = vdPlot;
    StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

    for(int i=0; i<=nPlotPoints; i++) {
      vdPlot[i] = oParameters.pdParMinBound[nPlotParameterNumber] +
	i*((oParameters.pdParMaxBound[nPlotParameterNumber] -
	      oParameters.pdParMinBound[nPlotParameterNumber])/nPlotPoints);
      vdFunctionValues[i] = 0;
    } // for i

    vdTemp = vdFunctionValues;
    vdSigma = vdFunctionValues;
    vdMeanValues = vdFunctionValues;
    vdMeanSquareValues = vdFunctionValues;

    // Read the first element
    if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
    if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 9;
    if(oConstants.bCalcTMPars) {
      if(! fread(oDat.pdChi, sizeof(double), nFileTMParameters, pFile) ) throw 10;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 11;
    } // if oConstants.bCalcTMPars

    nIndex = int((nPlotPoints-1)*(oDat.pdPar[nPlotParameterNumber] - oParameters.pdParMinBound[nPlotParameterNumber])/(oParameters.pdParMaxBound[nPlotParameterNumber]- oParameters.pdParMinBound[nPlotParameterNumber]));
    vdRawValues[0] = oDat.pdPar[nPlotParameterNumber];

    // Now integrate the values! Work through all computer data //
    if(nIndex >= 0 && nIndex <= nPlotPoints) vdFunctionValues[nIndex] += 1;// else {printf("%i ", nIndex); }
    dPrevious = oDat.dLogLik;
    for(int nStep=1; nStep<nMCMCSteps; nStep++) {
      // Read the next element
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 9;
    if(oConstants.bCalcTMPars) {
	if(! fread(oDat.pdChi, sizeof(double), nFileTMParameters, pFile) ) throw 10;
	if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 11;
    } // if bCalcTMPars

      nIndex = int((nPlotPoints-1)*(oDat.pdPar[nPlotParameterNumber] - oParameters.pdParMinBound[nPlotParameterNumber])/(oParameters.pdParMaxBound[nPlotParameterNumber]- oParameters.pdParMinBound[nPlotParameterNumber]));
      vdRawValues[nStep] = oDat.pdPar[nPlotParameterNumber];

      if(nIndex >= 0 && nIndex <= nPlotPoints) vdFunctionValues[nIndex] += 1;
      dPrevious = oDat.dLogLik;

      if(nStep % 10 == 0) DrawProgressBar(int(100.0*nStep/nMCMCSteps));
    } // for nStep

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    fclose(pFile);
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
  } // try
  FinishProgressBar();

  vdTemp = vdFunctionValues/(Max(vdFunctionValues));
  WritePlot("mcmcplot-1d-raw.txt", vdRawValues);

  vdX = vdPlot;
  vdY = vdTemp;
} // Calculate1DMCMCIntegrationOnTheFly



/* This function performs the actual MCMC integration and returns the result in
 * the vdX, vdY and vdYErr vectors
 * */
void Calculate1DMCMCIntegration(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nPlotParameterNumber, CVector &vdX, CVector &vdY, CVector &vdYErr, bool bCalcErr) {
  int nMCMCSteps, nPlotPoints, nRandIndex;
  CNumber ndRand1;
  CMCMCDataPoint *pDat;

  nPlotPoints = oConstants.nSmallPlotPoints;
  nMCMCSteps = oConstants.nMCMCSteps;

//  pDat = new MCMCDataPoint[nMCMCSteps];
  if(! ReadMCMCDataFile(nMCMCSteps, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName, oConstants.bCalcTMPars) ) {
    PrintFailed();
    return;
  }

  InitProgressBar("Making MCMC plot...");
  DrawProgressBar(0);

  // Prepare for integrating: initialize the variables
  CVector vdPlot, vdFunctionValues, vdTemp, vdTemp2, vdMeanValues, vdMeanSquareValues, vdSigma;
  int nIndex;
  double dPrevious;
  vdPlot.Initialize(nPlotPoints + 1);
  vdFunctionValues = vdPlot;
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

  for(int i=0; i<=nPlotPoints; i++) {
    vdPlot[i] = oParameters.pdParMinBound[nPlotParameterNumber] +
      i*((oParameters.pdParMaxBound[nPlotParameterNumber] -
	    oParameters.pdParMinBound[nPlotParameterNumber])/nPlotPoints);
    vdFunctionValues[i] = 0;
  } // for i

  vdTemp = vdFunctionValues;
  vdSigma = vdFunctionValues;
  vdMeanValues = vdFunctionValues;
  vdMeanSquareValues = vdFunctionValues;

  nIndex = int((nPlotPoints-1)*(pDat[0].pdPar[nPlotParameterNumber] - oParameters.pdParMinBound[nPlotParameterNumber])/(oParameters.pdParMaxBound[nPlotParameterNumber]- oParameters.pdParMinBound[nPlotParameterNumber]));

  // Now integrate the values! Work through all computer data //
  if(nIndex >= 0 && nIndex <= nPlotPoints) vdFunctionValues[nIndex] += 1;// else {printf("%i ", nIndex); }
  dPrevious = pDat[0].dLogLik;
  for(int nStep=1; nStep<nMCMCSteps; nStep++) {
    nIndex = int((nPlotPoints-1)*(pDat[nStep].pdPar[nPlotParameterNumber] - oParameters.pdParMinBound[nPlotParameterNumber])/(oParameters.pdParMaxBound[nPlotParameterNumber]- oParameters.pdParMinBound[nPlotParameterNumber]));

    if(nIndex >= 0 && nIndex <= nPlotPoints) vdFunctionValues[nIndex] += 1;
    dPrevious = pDat[nStep].dLogLik;
  } // for nStep


  if(bCalcErr) {
    // Now do the same thing, many times again, to estimate the bootstrap-errors!
    for(int i=0; i<oConstants.nBootstrapAttempts; i++) {
      ndRand1.Randomize();
      nRandIndex = int( double(ndRand1)*(nMCMCSteps-2)+1.5 );
      nIndex = int((nPlotPoints-1)*(pDat[nRandIndex].pdPar[nPlotParameterNumber] - oParameters.pdParMinBound[nPlotParameterNumber])/(oParameters.pdParMaxBound[nPlotParameterNumber]- oParameters.pdParMinBound[nPlotParameterNumber]));

      // Now integrate the values! Work through all computer data //
      if(nIndex >= 0 && nIndex <= nPlotPoints) vdTemp[nIndex] += 1;
      for(int nStep=1; nStep<nMCMCSteps; nStep++) {
	ndRand1.Randomize();
	nRandIndex = int( double(ndRand1)*(nMCMCSteps-2)+1.5 );

	// Add a new one
	nIndex = int((nPlotPoints-1)*(pDat[nRandIndex].pdPar[nPlotParameterNumber] - oParameters.pdParMinBound[nPlotParameterNumber])/(oParameters.pdParMaxBound[nPlotParameterNumber]- oParameters.pdParMinBound[nPlotParameterNumber]));

	if(nIndex >= 0 && nIndex <= nPlotPoints) vdTemp[nIndex] += 1;
      } // for i
      for(int j=0; j<nPlotPoints; j++) {
	vdMeanValues[j] += double(vdTemp[j]);
	vdMeanSquareValues[j] += double(vdTemp[j]) * double(vdTemp[j]);
	vdTemp[j] = 0;
      } // for j
      if(i % 500 == 0)
	DrawProgressBar(int(100.0*i/oConstants.nBootstrapAttempts));
    } // int i

    for(int i=0; i<nPlotPoints; i++) {
      vdSigma[i] = sqrt((double(vdMeanSquareValues[i])/oConstants.nBootstrapAttempts) - (double(vdMeanValues[i])/oConstants.nBootstrapAttempts)*(double(vdMeanValues[i])/oConstants.nBootstrapAttempts));
    } // for i
  } // if bCalcErr

  FinishProgressBar();

  vdTemp = vdFunctionValues/(Max(vdFunctionValues));
  if(bCalcErr)
    vdTemp2 = vdSigma/(Max(vdFunctionValues));

  vdX = vdPlot;
  vdY = vdTemp;
  vdYErr = vdTemp2;

  // This should call for the destructors...
  delete[] pDat;
} // Calculate1DMCMCIntegration

/* Integrates the MCMC data over all but one parameter. This parameter is
 * usually the GWB amplitude.
 *
 * */
void IntegrateMCMCData(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nPlotParameterNumber, bool bCalcErr) {
  CVector vdX, vdY, vdYErr;
  char strBuf[160];

  Calculate1DMCMCIntegration(oData, oParameters, oConstants, strFileName, nPlotParameterNumber, vdX, vdY, vdYErr, bCalcErr);

  // Writa data to file so it can be used in with gnuplot
  strcpy(strBuf, oConstants.strDataDir);
  if(strBuf[strlen(strBuf)-1] != '/') strcat(strBuf, "/");
  strcat(strBuf, "plotmcmcdata-1d.txt");

  if(bCalcErr)
    WritePlot(strBuf, vdX, vdY, vdYErr);
  else
    WritePlot(strBuf, vdX, vdY);
} // IntegrateMCMCData


/* This function performs the actual MCMC integration and returns the result in
 * the vdX, vdY and mdZ objects
 * */
void Calculate3DMCMCIntegrationOnTheFly(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ) {
  int nMCMCSteps, nPlotPoints, nRandIndex;
  int nSourceNumber, nParameterNumber;
  int nFileParameters, nFileTMParameters;
  unsigned long nFileLength;
  FILE *pFile;
  char strMsg[160], strFileVersion[16];
  CMCMCDataPoint oDat;
  CNumber ndRand1;
  CVector vdPlotX, vdPlotY, vdFunctionValues, vdTemp, vdTemp2, vdTemp3, vdMeanValues, vdMeanSquareValues, vdSigma, vdRawValuesX, vdRawValuesY;
  CMatrix mdFunctionValues;
  int nIndexX, nIndexY;
  double dPrevious;

  nPlotPoints = int(sqrt(double(oConstants.nPlotPoints)));
  nMCMCSteps = oConstants.nMCMCSteps;

  strcpy(strMsg, "Making MCMC plot from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "' on the fly...\n");

  fprintf(stderr, "Using integration numbers: %i,%i\n", nParameter1, nParameter2);
  fprintf(stderr, strMsg);
  InitProgressBar("Integrating");
  DrawProgressBar(0);

  try {
    // Open the file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(nFileParameters != oConstants.nParameters) {
      printf("WARNING: nFileParameters != nParameters\n");
      printf("nFileParameters: %i,  nParameters: %i\n", nFileParameters, oConstants.nParameters);
    } // if nFileParameters
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;
    if(nFileTMParameters != oConstants.nMarParameters) {
      printf("WARNING: nFileTMParameters != TMnParameters\n");
      printf("nFileTMParameters: %i,  nTMParameters: %i\n", nFileTMParameters, oConstants.nMarParameters);
    } // if nFileTMParameters

    // Calculate nMCMCSteps
    if(oConstants.bCalcTMPars) {
      nMCMCSteps = (nFileLength-sizeof(char)*16 - 2*sizeof(int)) /
	(
	  sizeof(double)*(nFileParameters+1+nFileTMParameters*(nFileTMParameters+1))
	);
      if(nMCMCSteps <= 0) throw 7;
    } else {
      nMCMCSteps = (nFileLength-sizeof(char)*16 - 2*sizeof(int)) /
	(
	  sizeof(double)*(nFileParameters+1)
	);
    } // if bCalcTMPars

    // Leave the file open 'till the end
    oDat.Initialize(nFileParameters, nFileTMParameters);
    vdRawValuesX.Initialize(nMCMCSteps);
    vdRawValuesY.Initialize(nMCMCSteps);

    // Prepare for integrating: initialize the variables
    vdPlotX.Initialize(nPlotPoints + 1);
    vdPlotY.Initialize(nPlotPoints + 1);
    mdFunctionValues.Initialize(nPlotPoints + 1, nPlotPoints + 1);
    vdFunctionValues = vdPlotX;

    ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

    for(int i=0; i<=nPlotPoints; i++) {
      vdPlotX[i] = oParameters.pdParMinBound[nParameter1] +
	i*((oParameters.pdParMaxBound[nParameter1] -
	      oParameters.pdParMinBound[nParameter1])/nPlotPoints);
      vdPlotY[i] = oParameters.pdParMinBound[nParameter2] +
	i*((oParameters.pdParMaxBound[nParameter2] -
	      oParameters.pdParMinBound[nParameter2])/nPlotPoints);
      vdFunctionValues[i] = 0;
      for(int j=0; j<= nPlotPoints; j++)
	mdFunctionValues[i][j] = 0;
    } // for i
    vdTemp = vdFunctionValues;
    vdSigma = vdFunctionValues;
    vdMeanValues = vdFunctionValues;
    vdMeanSquareValues = vdFunctionValues;

    // Read the first element
    if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
    if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 9;
    if(oConstants.bCalcTMPars) {
      if(! fread(oDat.pdChi, sizeof(double), nFileTMParameters, pFile) ) throw 10;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 11;
    } // if bCalcTMPars

    nIndexX = int((nPlotPoints-1)*(oDat.pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
    nIndexY = int((nPlotPoints-1)*(oDat.pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));
    vdRawValuesX[0] = oDat.pdPar[nParameter1];
    vdRawValuesY[0] = oDat.pdPar[nParameter2];

    if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
      vdFunctionValues[nIndexX] += 1;// else {printf("%i ", nIndex);
      mdFunctionValues[nIndexX][nIndexY] += 1;
    } // if nIndex
    dPrevious = oDat.dLogLik;

    // Now integrate the values! Work through all computer data //
    for(int nStep=1; nStep<nMCMCSteps; nStep++) {
      // Read the next element
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 9;
      if(oConstants.bCalcTMPars) {
	if(! fread(oDat.pdChi, sizeof(double), nFileTMParameters, pFile) ) throw 10;
	if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 11;
      } // if bCalcTMPars

      nIndexX = int((nPlotPoints-1)*(oDat.pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
      nIndexY = int((nPlotPoints-1)*(oDat.pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));
      vdRawValuesX[nStep] = oDat.pdPar[nParameter1];
      vdRawValuesY[nStep] = oDat.pdPar[nParameter2];

      if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
	vdFunctionValues[nIndexX] += 1;
	mdFunctionValues[nIndexX][nIndexY] += 1;
      } // if nIndex

      dPrevious = oDat.dLogLik;

      DrawProgressBar(int(100.0*nStep/nMCMCSteps));
    } // for nStep

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    fclose(pFile);
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
  } // try
  FinishProgressBar();

  vdX = vdPlotX;
  vdY = vdPlotY;
  mdZ = mdFunctionValues;

  WritePlot("mcmcplot-2d-raw.txt", vdRawValuesX, vdRawValuesY);

  PrintSuccess();
} // Calculate3DMCMCIntegrationOnTheFly



/* This function performs the actual MCMC integration and returns the result in
 * the vdX, vdY and mdZ objects
 * */
void Calculate3DMCMCIntegration(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ) {
  int nMCMCSteps, nPlotPoints, nRandIndex;
  int nSourceNumber, nParameterNumber;
//  int nParameter1, nParameter2;
  CNumber ndRand1;
  CMCMCDataPoint *pDat;

  nPlotPoints = int(sqrt(double(oConstants.nPlotPoints)));
  nMCMCSteps = oConstants.nMCMCSteps;

//  pDat = new MCMCDataPoint[nMCMCSteps];
  if(! ReadMCMCDataFile(nMCMCSteps, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName, oConstants.bCalcTMPars) ) {
    PrintFailed();
    return;
  }

  printf("Using integration numbers: %i,%i\n", nParameter1, nParameter2);
  PrintStatus("Integrating datapoints and generating plot...");

  // Prepare for integrating: initialize the variables
  CVector vdPlotX, vdPlotY, vdFunctionValues, vdTemp, vdTemp2, vdTemp3, vdMeanValues, vdMeanSquareValues, vdSigma;
  CMatrix mdFunctionValues;
  int nIndexX, nIndexY;
  double dPrevious;
  vdPlotX.Initialize(nPlotPoints + 1);
  vdPlotY.Initialize(nPlotPoints + 1);
  mdFunctionValues.Initialize(nPlotPoints + 1, nPlotPoints + 1);
  vdFunctionValues = vdPlotX;

  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

  for(int i=0; i<=nPlotPoints; i++) {
    vdPlotX[i] = oParameters.pdParMinBound[nParameter1] +
      i*((oParameters.pdParMaxBound[nParameter1] -
	    oParameters.pdParMinBound[nParameter1])/nPlotPoints);
    vdPlotY[i] = oParameters.pdParMinBound[nParameter2] +
      i*((oParameters.pdParMaxBound[nParameter2] -
	    oParameters.pdParMinBound[nParameter2])/nPlotPoints);
    vdFunctionValues[i] = 0;
    for(int j=0; j<= nPlotPoints; j++)
      mdFunctionValues[i][j] = 0;
  } // for i
  vdTemp = vdFunctionValues;
  vdSigma = vdFunctionValues;
  vdMeanValues = vdFunctionValues;
  vdMeanSquareValues = vdFunctionValues;

  nIndexX = int((nPlotPoints-1)*(pDat[0].pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
  nIndexY = int((nPlotPoints-1)*(pDat[0].pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

  if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
    vdFunctionValues[nIndexX] += 1;// else {printf("%i ", nIndex);
    mdFunctionValues[nIndexX][nIndexY] += 1;
  } // if nIndex
  dPrevious = pDat[0].dLogLik;

  // Now integrate the values! Work through all computer data //
  for(int nStep=1; nStep<nMCMCSteps; nStep++) {
    nIndexX = int((nPlotPoints-1)*(pDat[nStep].pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
    nIndexY = int((nPlotPoints-1)*(pDat[nStep].pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

    if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
      vdFunctionValues[nIndexX] += 1;
      mdFunctionValues[nIndexX][nIndexY] += 1;
    } // if nIndex

      dPrevious = pDat[nStep].dLogLik;
  } // for nStep

  vdX = vdPlotX;
  vdY = vdPlotY;
  mdZ = mdFunctionValues;

  PrintSuccess();

  // This should call for the destructors...
  delete[] pDat;
} // Calculate3DMCMCIntegration



/* This function performs the actual MCMC integration for one of the ensemble
 * sets and returns the result in the vdX, vdY and mdZ objects
 * */
void Calculate3DMCMCEnsembleIntegration(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, int nDataSet, CVector &vdX, CVector &vdY, CMatrix &mdZ) {
  int nMCMCSteps, nPlotPoints, nRandIndex;
  int nSourceNumber, nParameterNumber;
  double dMaxLL;
  CMCMCDataEnsemble *pDat;
  SParametersType oParameters;

  nPlotPoints = int(sqrt(double(oConstants.nPlotPoints)));
  nMCMCSteps = oConstants.nMCMCSteps;

  if(! ReadMCMCEnsembleDataFile(nMCMCSteps, oData.nDataSets, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName, oData.mdDataSets)) {
    return;
  } // if ReadMCMCEnsembleDataFile


  printf("Using integration numbers: %i,%i\n", nParameter1, nParameter2);
  PrintStatus("Integrating ensemble and generating plot...");

  // Prepare for integrating: initialize the variables
  CVector vdPlotX, vdPlotY, vdFunctionValues, vdTemp, vdTemp2, vdTemp3, vdMeanValues, vdMeanSquareValues, vdSigma;
  CMatrix mdFunctionValues;
  int nIndexX, nIndexY;
  vdPlotX.Initialize(nPlotPoints + 1);
  vdPlotY.Initialize(nPlotPoints + 1);
  mdFunctionValues.Initialize(nPlotPoints + 1, nPlotPoints + 1);
  vdFunctionValues = vdPlotX;

  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

  for(int i=0; i<=nPlotPoints; i++) {
    vdPlotX[i] = oParameters.pdParMinBound[nParameter1] +
      i*((oParameters.pdParMaxBound[nParameter1] -
	    oParameters.pdParMinBound[nParameter1])/nPlotPoints);
    vdPlotY[i] = oParameters.pdParMinBound[nParameter2] +
      i*((oParameters.pdParMaxBound[nParameter2] -
	    oParameters.pdParMinBound[nParameter2])/nPlotPoints);
    vdFunctionValues[i] = 0;
    for(int j=0; j<= nPlotPoints; j++)
      mdFunctionValues[i][j] = 0;
  } // for i
  vdTemp = vdFunctionValues;
  vdSigma = vdFunctionValues;
  vdMeanValues = vdFunctionValues;
  vdMeanSquareValues = vdFunctionValues;

  // Look for the maximum LL value first
  dMaxLL = pDat[0].dLogLik - pDat[0].pdLogLik[nDataSet];
  for(int nStep=1; nStep<nMCMCSteps; nStep++) {
    if(dMaxLL < pDat[nStep].dLogLik - pDat[nStep].pdLogLik[nDataSet]) {
      dMaxLL = pDat[nStep].dLogLik - pDat[nStep].pdLogLik[nDataSet];
    } // if dMaxLL
  } // for nStep

  nIndexX = int((nPlotPoints-1)*(pDat[0].pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
  nIndexY = int((nPlotPoints-1)*(pDat[0].pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

  if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
    vdFunctionValues[nIndexX] += exp(pDat[0].dLogLik - pDat[0].pdLogLik[nDataSet] - dMaxLL);
    mdFunctionValues[nIndexX][nIndexY] += exp(pDat[0].dLogLik - pDat[0].pdLogLik[nDataSet] - dMaxLL);
  } // if nIndex

  // Now integrate the values! Work through all computer data //
  for(int nStep=1; nStep<nMCMCSteps; nStep++) {
    nIndexX = int((nPlotPoints-1)*(pDat[nStep].pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
    nIndexY = int((nPlotPoints-1)*(pDat[nStep].pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

    if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
      vdFunctionValues[nIndexX] += exp(pDat[nStep].dLogLik - pDat[nStep].pdLogLik[nDataSet] - dMaxLL);
      mdFunctionValues[nIndexX][nIndexY] += exp(pDat[nStep].dLogLik - pDat[nStep].pdLogLik[nDataSet] - dMaxLL);
    } // if nIndex

  } // for nStep

  vdX = vdPlotX;
  vdY = vdPlotY;
  mdZ = mdFunctionValues;

  PrintSuccess();

  if(1) {
    int nIndex=0;

    PrintStatus("Changing the data of the pulsar to selected set...");
    for(int a=0; a<oConstants.k; a++) {
      for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	oConstants.poPulsars[a].pdResiduals[i] = double(oData.mdDataSets[i][nDataSet]);
	oData.vdData[nIndex] = oConstants.poPulsars[a].pdResiduals[i];
	nIndex++;
      } // for i
    } // for a
    PrintSuccess();
  }

  // This should call for the destructors...
  delete[] pDat;
} // Calculate3DMCMCEnsembleIntegration



/* This function performs the actual MCMC integration for one of the ensemble
 * sets and returns the result in the vdX, vdY and mdZ objects
 * */
void Calculate3DMCMCEnsembleIntegrationOnTheFly(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, int nDataSet, CVector &vdX, CVector &vdY, CMatrix &mdZ) {
  int nMCMCSteps, nPlotPoints, nRandIndex;
  int nSourceNumber, nParameterNumber;
  double dMaxLL;
  CMCMCDataEnsemble oDat;
  SParametersType oParameters;
  FILE *pFile;
  int nFileParameters, nFileTMParameters, nFileDataSets, nFileObservations;
  int nDataSetsLength;
  unsigned long nFileLength, nFilePos;

  char strMsg[160], strFileVersion[16];
  strcpy(strMsg, "Integrating data from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "' on the fly...");

  nPlotPoints = int(sqrt(double(oConstants.nPlotPoints)));
  nMCMCSteps = oConstants.nMCMCSteps;

  try {
//    if(! ReadMCMCEnsembleDataFile(nMCMCSteps, oData.nDataSets, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName, oData.mdDataSets)) {
//      return;
//    } // if ReadMCMCEnsembleDataFile
    // Open the file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileDataSets, sizeof(int), 1, pFile) ) throw 5;
    oData.nDataSets = nFileDataSets;
    if(! fread(&nFileObservations, sizeof(int), 1, pFile) ) throw 5;
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(nFileParameters != oConstants.nParameters) {
      printf("WARNING: nFileParameters != nParameters\n");
      printf("nFileParameters: %i,  nParameters: %i\n", nFileParameters, oConstants.nParameters);
    } // if nFileParameters
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;
    if(nFileTMParameters != oConstants.nMarParameters) {
      printf("WARNING: nFileTMParameters != nTMParameters\n");
      printf("nFileTMParameters: %i,  nTMParameters: %i\n", nFileTMParameters, oConstants.nMarParameters);
    } // if nFileTMParameters
    nDataSetsLength = nFileObservations*nFileDataSets;
    oData.mdDataSets.Initialize(nFileObservations, nFileDataSets);
    if(! fread(oData.mdDataSets.m_pdData, sizeof(double), nDataSetsLength, pFile) ) throw 6;

    // Calculate nMCMCSteps
    nMCMCSteps = (nFileLength-sizeof(char)*16 - 4*sizeof(int) - nDataSetsLength*sizeof(double)) / (
        sizeof(double)*(
	  nFileParameters + nFileDataSets + 1 +
	  nFileTMParameters*nFileDataSets +
	  nFileTMParameters*nFileTMParameters));

    if(nMCMCSteps <= 0) throw 7;

    // Leave the file open 'till the end
    oDat.Initialize(nFileDataSets, nFileParameters, nFileTMParameters);


    printf("Using integration numbers: %i,%i\n", nParameter1, nParameter2);
    PrintStatus(strMsg);

    // Prepare for integrating: initialize the variables
    CVector vdPlotX, vdPlotY, vdFunctionValues, vdTemp, vdTemp2, vdTemp3, vdMeanValues, vdMeanSquareValues, vdSigma;
    CMatrix mdFunctionValues;
    int nIndexX, nIndexY;
    vdPlotX.Initialize(nPlotPoints + 1);
    vdPlotY.Initialize(nPlotPoints + 1);
    mdFunctionValues.Initialize(nPlotPoints + 1, nPlotPoints + 1);
    vdFunctionValues = vdPlotX;

    ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

    for(int i=0; i<=nPlotPoints; i++) {
      vdPlotX[i] = oParameters.pdParMinBound[nParameter1] +
	i*((oParameters.pdParMaxBound[nParameter1] -
	      oParameters.pdParMinBound[nParameter1])/nPlotPoints);
      vdPlotY[i] = oParameters.pdParMinBound[nParameter2] +
	i*((oParameters.pdParMaxBound[nParameter2] -
	      oParameters.pdParMinBound[nParameter2])/nPlotPoints);
      vdFunctionValues[i] = 0;
      for(int j=0; j<= nPlotPoints; j++)
	mdFunctionValues[i][j] = 0;
    } // for i
    vdTemp = vdFunctionValues;
    vdSigma = vdFunctionValues;
    vdMeanValues = vdFunctionValues;
    vdMeanSquareValues = vdFunctionValues;

    // Look for the maximum LL value in the file. Save the current file
    // position.
    nFilePos = ftell(pFile);
    if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
    if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
    if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
    if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
    if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;
    dMaxLL = oDat.dLogLik - oDat.pdLogLik[nDataSet];
    for(int nStep=1; nStep<nMCMCSteps; nStep++) {
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;
      if(dMaxLL < oDat.dLogLik - oDat.pdLogLik[nDataSet]) {
	dMaxLL = oDat.dLogLik - oDat.pdLogLik[nDataSet];
      } // if dMaxLL
    } // for nStep
    // Set the file position back to where we were
    if(fseek(pFile, nFilePos, SEEK_SET)) throw 3;

    // Read the first element again
    if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
    if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
    if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
    if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
    if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;

    nIndexX = int((nPlotPoints-1)*(oDat.pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
    nIndexY = int((nPlotPoints-1)*(oDat.pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

    if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
      vdFunctionValues[nIndexX] += exp(oDat.dLogLik - oDat.pdLogLik[nDataSet] - dMaxLL);
      mdFunctionValues[nIndexX][nIndexY] += exp(oDat.dLogLik - oDat.pdLogLik[nDataSet] - dMaxLL);
    } // if nIndex

    // Now integrate the values! Work through all computer data //
    for(int nStep=1; nStep<nMCMCSteps; nStep++) {
      // Read the next element
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;
      nIndexX = int((nPlotPoints-1)*(oDat.pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
      nIndexY = int((nPlotPoints-1)*(oDat.pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

      if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
	vdFunctionValues[nIndexX] += exp(oDat.dLogLik - oDat.pdLogLik[nDataSet] - dMaxLL);
	mdFunctionValues[nIndexX][nIndexY] += exp(oDat.dLogLik - oDat.pdLogLik[nDataSet] - dMaxLL);
      } // if nIndex

    } // for nStep

    vdX = vdPlotX;
    vdY = vdPlotY;
    mdZ = mdFunctionValues;

    fclose(pFile);
    PrintSuccess();

    if(1) {
      int nIndex=0;

      PrintStatus("Changing the data of the pulsar to selected set...");
      for(int a=0; a<oConstants.k; a++) {
	for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	  oConstants.poPulsars[a].pdResiduals[i] = double(oData.mdDataSets[i][nDataSet]);
	  oData.vdData[nIndex] = oConstants.poPulsars[a].pdResiduals[i];
	  nIndex++;
	} // for i
      } // for a
      PrintSuccess();
    } // Only for testing purposes
  } catch(int nError) {
    fclose(pFile);
    PrintFailed();
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
  } // try
} // Calculate3DMCMCEnsembleIntegrationOnTheFly


/* This function performs the actual MCMC integration for all of the ensemble
 * sets and returns the average result in the vdX, vdY and mdZ objects
 *
 * This is as it should be done in the MLDR combined algorithm. The per-set
 * loglikelihood does not matter: all datasets are given equal weight.
 * */
void Calculate3DMCMCEnsembleMLDRIntegrationOnTheFly(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ) {
  int nMCMCSteps, nPlotPoints, nRandIndex;
  int nSourceNumber, nParameterNumber;
  double dMaxLL;
  CMCMCDataEnsemble oDat;
  SParametersType oParameters;
  FILE *pFile;
  int nFileParameters, nFileTMParameters, nFileDataSets, nFileObservations;
  int nDataSetsLength;
  unsigned long nFileLength, nFilePos;

  char strMsg[160], strFileVersion[16];
  strcpy(strMsg, "Integrating data from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "' on the fly (MLDR)...");

  nPlotPoints = int(sqrt(double(oConstants.nPlotPoints)));
  nMCMCSteps = oConstants.nMCMCSteps;

  try {
//    if(! ReadMCMCEnsembleDataFile(nMCMCSteps, oData.nDataSets, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName, oData.mdDataSets)) {
//      return;
//    } // if ReadMCMCEnsembleDataFile
    // Open the file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileDataSets, sizeof(int), 1, pFile) ) throw 5;
    oData.nDataSets = nFileDataSets;
    if(! fread(&nFileObservations, sizeof(int), 1, pFile) ) throw 5;
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(nFileParameters != oConstants.nParameters) {
      printf("WARNING: nFileParameters != nParameters\n");
      printf("nFileParameters: %i,  nParameters: %i\n", nFileParameters, oConstants.nParameters);
    } // if nFileParameters
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;
    if(nFileTMParameters != oConstants.nMarParameters) {
      printf("WARNING: nFileTMParameters != nTMParameters\n");
      printf("nFileTMParameters: %i,  nTMParameters: %i\n", nFileTMParameters, oConstants.nMarParameters);
    } // if nFileTMParameters
    nDataSetsLength = nFileObservations*nFileDataSets;
    oData.mdDataSets.Initialize(nFileObservations, nFileDataSets);
    if(! fread(oData.mdDataSets.m_pdData, sizeof(double), nDataSetsLength, pFile) ) throw 6;

    // Calculate nMCMCSteps
    nMCMCSteps = (nFileLength-sizeof(char)*16 - 4*sizeof(int) - nDataSetsLength*sizeof(double)) / (
        sizeof(double)*(
	  nFileParameters + nFileDataSets + 1 +
	  nFileTMParameters*nFileDataSets +
	  nFileTMParameters*nFileTMParameters));

    if(nMCMCSteps <= 0) throw 7;

    // Leave the file open 'till the end
    oDat.Initialize(nFileDataSets, nFileParameters, nFileTMParameters);


    printf("Using integration numbers: %i,%i\n", nParameter1, nParameter2);
    PrintStatus(strMsg);

    // Prepare for integrating: initialize the variables
    CVector vdPlotX, vdPlotY, vdFunctionValues, vdTemp, vdTemp2, vdTemp3, vdMeanValues, vdMeanSquareValues, vdSigma;
    CMatrix mdFunctionValues;
    int nIndexX, nIndexY;
    vdPlotX.Initialize(nPlotPoints + 1);
    vdPlotY.Initialize(nPlotPoints + 1);
    mdFunctionValues.Initialize(nPlotPoints + 1, nPlotPoints + 1);
    vdFunctionValues = vdPlotX;

    ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

    for(int i=0; i<=nPlotPoints; i++) {
      vdPlotX[i] = oParameters.pdParMinBound[nParameter1] +
	i*((oParameters.pdParMaxBound[nParameter1] -
	      oParameters.pdParMinBound[nParameter1])/nPlotPoints);
      vdPlotY[i] = oParameters.pdParMinBound[nParameter2] +
	i*((oParameters.pdParMaxBound[nParameter2] -
	      oParameters.pdParMinBound[nParameter2])/nPlotPoints);
      vdFunctionValues[i] = 0;
      for(int j=0; j<= nPlotPoints; j++)
	mdFunctionValues[i][j] = 0;
    } // for i
    vdTemp = vdFunctionValues;
    vdSigma = vdFunctionValues;
    vdMeanValues = vdFunctionValues;
    vdMeanSquareValues = vdFunctionValues;

    // Look for the maximum LL value in the file. Save the current file
    // position.
    nFilePos = ftell(pFile);
    
    if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
    if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
    if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
    if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
    if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;
    dMaxLL = oDat.dLogLik - oDat.pdLogLik[2];
    for(int nSet=2; nSet<nFileDataSets; ++nSet) {
      if(dMaxLL < oDat.dLogLik - oDat.pdLogLik[nSet]) {
	dMaxLL = oDat.dLogLik - oDat.pdLogLik[nSet];
      } // if dMaxLL
    } // for nSet
    for(int nStep=1; nStep<nMCMCSteps; nStep++) {
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;
      for(int nSet=2; nSet<nFileDataSets; ++nSet) {
	if(dMaxLL < oDat.dLogLik - oDat.pdLogLik[nSet]) {
	  dMaxLL = oDat.dLogLik - oDat.pdLogLik[nSet];
	} // if dMaxLL
      } // for nSet
    } // for nStep
    // Set the file position back to where we were
    if(fseek(pFile, nFilePos, SEEK_SET)) throw 3;

    // Read the first element again
    if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
    if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
    if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
    if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
    if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;

    nIndexX = int((nPlotPoints-1)*(oDat.pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
    nIndexY = int((nPlotPoints-1)*(oDat.pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

    for(int nSet=2; nSet<nFileDataSets; ++nSet) {
      if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
	vdFunctionValues[nIndexX] += exp(oDat.dLogLik - oDat.pdLogLik[nSet] - dMaxLL);
	mdFunctionValues[nIndexX][nIndexY] += exp(oDat.dLogLik - oDat.pdLogLik[nSet] - dMaxLL);
      } // if nIndex
    } // for nSet

    // Now integrate the values! Work through all computer data //
    for(int nStep=1; nStep<nMCMCSteps; nStep++) {
      // Read the next element
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;
      nIndexX = int((nPlotPoints-1)*(oDat.pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
      nIndexY = int((nPlotPoints-1)*(oDat.pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

      for(int nSet=2; nSet<nFileDataSets; ++nSet) {
	if(nIndexX >= 0 && nIndexX <= nPlotPoints && nIndexY >= 0 && nIndexY <= nPlotPoints) {
#if 0
	  fprintf(stderr, "(%i, %i) = exp(%e - %e - %e) = exp(%e) = %e\n",
	      nStep, nSet,
	      oDat.dLogLik, oDat.pdLogLik[nSet], dMaxLL,
	      oDat.dLogLik - oDat.pdLogLik[nSet] - dMaxLL, exp(oDat.dLogLik - oDat.pdLogLik[nSet] - dMaxLL));
#endif
	  vdFunctionValues[nIndexX] += exp(oDat.dLogLik - oDat.pdLogLik[nSet] - dMaxLL);
	  mdFunctionValues[nIndexX][nIndexY] += exp(oDat.dLogLik - oDat.pdLogLik[nSet] - dMaxLL);
	} // if nIndex
      } // for nSet

    } // for nStep

    vdX = vdPlotX;
    vdY = vdPlotY;
    mdZ = mdFunctionValues;

    fclose(pFile);
    PrintSuccess();

  } catch(int nError) {
    fclose(pFile);
    PrintFailed();
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
  } // try
} // Calculate3DMCMCEnsembleMLDRIntegrationOnTheFly






/* This function checks whether or not we have to display marginalisation
 * parameters, or stochastic parameters. For either, it finds the min/max
 * bounds, and returns the indices for in mdBlock or in pdPar
 * */
void FindParameterBounds(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameter1, int nParameter2, bool &bMarPar1, bool &bMarPar2, int &nP1, int &nP2, double &dMinBound1, double &dMinBound2, double &dMaxBound1, double &dMaxBound2, CMatrix &mdCXi) {
  int nStochasticIndex=0, nMarIndex=0, nCIndex=0;

  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.eID == SID_Deterministic &&
	  (oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p))) {
	if(nParameter1 == oConstants.poSources[s].nFirstParIndex + p){
	  nP1 = nMarIndex;
	  bMarPar1 = true;
	  dMinBound1 = -15 * sqrt(double(mdCXi[nP1][nP1]));
	  dMaxBound1 = 15 * sqrt(double(mdCXi[nP1][nP1]));
	}
	if(nParameter2 == oConstants.poSources[s].nFirstParIndex + p){
	  nP2 = nMarIndex;
	  bMarPar2 = true;
	  dMinBound2 = -15 * sqrt(double(mdCXi[nP2][nP2]));
	  dMaxBound2 = 15 * sqrt(double(mdCXi[nP2][nP2]));
	}
	nMarIndex++;
      } else if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	if(nParameter1 == oConstants.poSources[s].nFirstParIndex + p){
	  nP1 = nParameter1;
	  bMarPar1 = false;
	  dMinBound1 = oParameters.pdParMinBound[nP1];
	  dMaxBound1 = oParameters.pdParMaxBound[nP1];
	}
	if(nParameter2 == oConstants.poSources[s].nFirstParIndex + p){
	  nP2 = nParameter2;
	  bMarPar2 = false;
	  dMinBound2 = oParameters.pdParMinBound[nP2];
	  dMaxBound2 = oParameters.pdParMaxBound[nP2];
	}
	nStochasticIndex++;
      }// if FITPARAMETER
    } // for p
  } // for  s
} // FindParameterBounds



/* This function performs the actual MCMC integration and returns the result in
 * the vdX, vdY and mdZ objects
 *
 * This version can do Marginalisation parameters as well
 * */
void Calculate3DMCMCIntegrationIncMar(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, int &nPR1, int &nPR2, CVector &vdX, CVector &vdY, CMatrix &mdZ, CMatrix &mdZML, CMatrix &mdZT2) {
  int nMCMCSteps, nPlotPoints, nRandIndex;
  int nSourceNumber, nParameterNumber;
  double dEvidence;
  CMCMCDataPoint *pDat;
  int nMax;
  double dLLMax;

  nPlotPoints = int(sqrt(double(oConstants.nPlotPoints)));
  nMCMCSteps = oConstants.nMCMCSteps;

//  pDat = new MCMCDataPoint[nMCMCSteps];
  if(! ReadMCMCDataFile(nMCMCSteps, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName) ) {
    PrintFailed();
    return;
  }

  // Prepare for integrating: initialize the variables
  CVector vdPlotX, vdPlotY, vdFunctionValues, vdTemp, vdTemp2, vdTemp3, vdMeanValues, vdMeanSquareValues, vdSigma;
  CMatrix mdFunctionValues, mdFunctionValuesML, mdFunctionValuesT2;
  int nIndexX, nIndexY;
  vdPlotX.Initialize(nPlotPoints + 1);
  vdPlotY.Initialize(nPlotPoints + 1);
  mdFunctionValues.Initialize(nPlotPoints + 1, nPlotPoints + 1);
  vdFunctionValues = vdPlotX;

  bool bMarPar1=false, bMarPar2=false;				// Marginalised p y/n
  double dMinBound1, dMinBound2, dMaxBound1, dMaxBound2;	// Bounds
  int nP1, nP2;							// True index

  CMatrix mdLG, mdLGTrans, mdTemp, mdCXi, mdTemp2, mdE, mdBlockTrans;
  CVector vdXi, vdChi;
  CNumber ndNum;
  double dLogDetLCL, dTemp;
  int nIndex1, nIndex2;

  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

  // Find the MCMC point with the highest likelihood value
  nMax = 0; dLLMax = pDat[0].dLogLik;
  for(int nStep=1; nStep<nMCMCSteps; nStep++) {
    if(dLLMax < pDat[nStep].dLogLik) {
      dLLMax = pDat[nStep].dLogLik;
      nMax = nStep;
    } // if dLLMax
  } // for nStep

  mdCXi.Initialize(oConstants.nMarParameters, oConstants.nMarParameters);

  for(int i=0; i<oConstants.nMarParameters*oConstants.nMarParameters; i++) {
    mdCXi.m_pdData[i] = pDat[nMax].pdCXiInv[i];
  } // for i

  // Figure out whether or not we have marginalised parameters to display
  FindParameterBounds(oData, oParameters, oConstants, nParameter1, nParameter2, bMarPar1, bMarPar2, nP1, nP2, dMinBound1, dMinBound2, dMaxBound1, dMaxBound2, mdCXi);

  printf("Using integration numbers: %i,%i = %i,%i\n", nParameter1, nParameter2, nP1, nP2);
  PrintStatus("Integrating datapoints and generating plot...");

  for(int i=0; i<=nPlotPoints; i++) {
    vdPlotX[i] = dMinBound1 + i*((dMaxBound1 - dMinBound1)/nPlotPoints);
    vdPlotY[i] = dMinBound2 + i*((dMaxBound2 - dMinBound2)/nPlotPoints);
    vdFunctionValues[i] = 0;
    for(int j=0; j<= nPlotPoints; j++)
      mdFunctionValues[i][j] = 0;
  } // for i

  vdSigma = vdFunctionValues;
  vdMeanValues = vdFunctionValues;
  vdMeanSquareValues = vdFunctionValues;

  vdChi.Initialize(oConstants.nMarParameters);
  if(bMarPar1==true && bMarPar2==false) {
    mdLG.Initialize(oConstants.nMarParameters, 1);
    vdXi.Initialize(1);
    for(int i=0; i<oConstants.nMarParameters; i++) {
      mdLG[i][0] = (i == nP1 ? 1 : 0);
    } // for i
    mdLGTrans = mdLG[LO_TRANSPOSE];
  }
  if(bMarPar1==false && bMarPar2==true) {
    mdLG.Initialize(oConstants.nMarParameters, 1);
    vdXi.Initialize(1);
    for(int i=0; i<oConstants.nMarParameters; i++) {
      mdLG[i][0] = (i == nP2 ? 1 : 0);
    } // for i
    mdLGTrans = mdLG[LO_TRANSPOSE];
  }
  if(bMarPar1==true && bMarPar2==true) {
    mdLG.Initialize(oConstants.nMarParameters, 2);
    vdXi.Initialize(2);
    for(int i=0; i<oConstants.nMarParameters; i++) {
      mdLG[i][0] = (i == nP1 ? 1 : 0);
      mdLG[i][1] = (i == nP2 ? 1 : 0);
    } // for i
    mdLGTrans = mdLG[LO_TRANSPOSE];
  }

  for(int nStep=0; nStep<nMCMCSteps; nStep++) {
//  for(int nStep=0; nStep<100; nStep++) {
    // 4 possibilities
    if(bMarPar1==false && bMarPar2==false) {
      // Both stochastic parameters
      nIndexX = int((nPlotPoints-1)*(pDat[nStep].pdPar[nP1] - dMinBound1) /
	  (dMaxBound1 - dMinBound1));
      nIndexY = int((nPlotPoints-1)*(pDat[nStep].pdPar[nP2] - dMinBound2) /
	  (dMaxBound2 - dMinBound2));
      if(nIndexX >= 0 && nIndexX <= nPlotPoints &&
	 nIndexY >= 0 && nIndexY <= nPlotPoints ) {
	mdFunctionValues[nIndexX][nIndexY] += 1;
      } // if nIndex
    } // if bMarPar
    if(bMarPar1==true && bMarPar2==false) {
      for(int i=0; i<oConstants.nMarParameters*oConstants.nMarParameters; i++) {
	mdCXi.m_pdData[i] = pDat[nStep].pdCXiInv[i];
      } // for i
      for(int i=0; i<oConstants.nMarParameters; i++) {
	vdChi.m_pdData[i] = pDat[nStep].pdChi[i];
      } // for i
      mdTemp = mdLGTrans * (mdCXi * mdLG);
      mdTemp.InvertChol(&dLogDetLCL);

      nIndexY = int((nPlotPoints-1)*(pDat[nStep].pdPar[nP2] - dMinBound2) /
	  (dMaxBound2 - dMinBound2));

      if(nIndexY >= 0 && nIndexY <= nPlotPoints) {
	for(int i=0; i<=nPlotPoints; i++) {
	  vdXi[0] = double(vdPlotX[i]) - double(vdChi[nP1]);
	  vdTemp = mdTemp * vdXi;
	  dTemp = 0;
	  for(int k=0; k<vdXi.m_pnDimSize[0]; k++)
	    dTemp += double(vdXi[k])*double(vdTemp[k]);
	  mdFunctionValues[i][nIndexY] += (1.0 / sqrt(2*M_PI)) * exp(-0.5*
	      (dLogDetLCL + dTemp));
	} // for i
      } // if nIndex
    } // if bMarPar
    if(bMarPar1==false && bMarPar2==true) {
      for(int i=0; i<oConstants.nMarParameters*oConstants.nMarParameters; i++) {
	mdCXi.m_pdData[i] = pDat[nStep].pdCXiInv[i];
      } // for i
      for(int i=0; i<oConstants.nMarParameters; i++) {
	vdChi.m_pdData[i] = pDat[nStep].pdChi[i];
      } // for i
      mdTemp = mdLGTrans * (mdCXi * mdLG);
      mdTemp.InvertChol(&dLogDetLCL);

      nIndexX = int((nPlotPoints-1)*(pDat[nStep].pdPar[nP1] - dMinBound1) /
	  (dMaxBound1 - dMinBound1));

      if(nIndexX >= 0 && nIndexX <= nPlotPoints) {
	for(int j=0; j<=nPlotPoints; j++) {
	  vdXi[0] = double(vdPlotY[j]) - double(vdChi[nP2]);
	  vdTemp = mdTemp * vdXi;
	  dTemp = 0;
	  for(int k=0; k<vdXi.m_pnDimSize[0]; k++)
	    dTemp += double(vdXi[k])*double(vdTemp[k]);
	  mdFunctionValues[nIndexX][j] += (1.0 / sqrt(2*M_PI)) * exp(-0.5*
	      (dLogDetLCL + dTemp));
	} // for j
      } // if nIndex
    } // if bMarPar
    if(bMarPar1==true && bMarPar2==true) {
      for(int i=0; i<oConstants.nMarParameters*oConstants.nMarParameters; i++) {
	mdCXi.m_pdData[i] = pDat[nStep].pdCXiInv[i];
      } // for i
      for(int i=0; i<oConstants.nMarParameters; i++) {
	vdChi.m_pdData[i] = pDat[nStep].pdChi[i];
      } // for i
      mdTemp = mdLGTrans * (mdCXi * mdLG);
      mdTemp.InvertChol(&dLogDetLCL);

      for(int i=0; i<=nPlotPoints; i++) {
	for(int j=0; j<=nPlotPoints; j++) {
	  vdXi[0] = double(vdPlotX[i]) - double(vdChi[nP1]);
	  vdXi[1] = double(vdPlotY[j]) - double(vdChi[nP2]);
	  vdTemp = mdTemp * vdXi;
	  dTemp = 0;
	  for(int k=0; k<vdXi.m_pnDimSize[0]; k++)
	    dTemp += double(vdXi[k])*double(vdTemp[k]);
	  mdFunctionValues[i][j] += (1.0 / sqrt(4*M_PI*M_PI)) * exp(-0.5*
	      (dLogDetLCL + dTemp));
	} // for j
      } // for i
    } // if bMarPar
  } // for nStep

  if(bMarPar1)
    nPR1 = nP1;
  else
    nPR1 = -1;
  if(bMarPar2)
    nPR2 = nP2;
  else
    nPR2 = -1;

  vdX = vdPlotX;
  vdY = vdPlotY;

  mdZ = mdFunctionValues;

  // For 2 marginalised parameters, also calculate the ML plot and the T2 plot
  if(bMarPar1==true && bMarPar2==true) {
    mdFunctionValuesML.Initialize(nPlotPoints + 1, nPlotPoints + 1);
    mdFunctionValuesT2.Initialize(nPlotPoints + 1, nPlotPoints + 1);

    for(int i=0; i<oConstants.nMarParameters*oConstants.nMarParameters; i++) {
      mdCXi.m_pdData[i] = pDat[nMax].pdCXiInv[i];
    } // for i

    for(int i=0; i<oConstants.nMarParameters; i++) {
      vdChi.m_pdData[i] = pDat[nMax].pdChi[i];
    } // for i

    mdTemp = mdLGTrans * (mdCXi * mdLG);
    mdTemp.InvertChol(&dLogDetLCL);

    for(int i=0; i<=nPlotPoints; i++) {
      for(int j=0; j<=nPlotPoints; j++) {
	vdXi[0] = double(vdPlotX[i]) - double(vdChi[nP1]);
	vdXi[1] = double(vdPlotY[j]) - double(vdChi[nP2]);
	vdTemp = mdTemp * vdXi;
	dTemp = 0;
	for(int k=0; k<vdXi.m_pnDimSize[0]; k++)
	  dTemp += double(vdXi[k])*double(vdTemp[k]);
	mdFunctionValuesML[i][j] += (1.0 / sqrt(4*M_PI*M_PI)) * exp(-0.5*
	    (dLogDetLCL + dTemp));
      } // for j
    } // for i
    mdZML = mdFunctionValuesML;

    // The tempo2 errors plot (no efac included)
    mdE.Initialize(oConstants.n, oConstants.n);
    nIndex1 = 0; nIndex2 = 0;
    for(int a=0; a<oConstants.k; a++) {
      for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	for(int b=0; b<oConstants.k; b++) {
	  for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
	    if(nIndex1 == nIndex2) {
	      mdE[nIndex1][nIndex2] = 1.0 / (oConstants.poPulsars[a].pdDeltaResiduals[i]*oConstants.poPulsars[a].pdDeltaResiduals[i]);
	    } else {
	      mdE[nIndex1][nIndex2] = 0;
	    }// if nIndex
	    nIndex2++;
	  } // for j
	} // for b
	nIndex2 = 0;
	nIndex1++;
      } // for i
    } // for p
    mdBlockTrans = oData.mdBlock[LO_TRANSPOSE];

    mdCXi = mdBlockTrans * (mdE * oData.mdBlock);
    mdCXi.InvertChol();

    // Calculate vdChi
    vdChi = mdCXi * (mdBlockTrans * (mdE * oData.vdData));

    mdTemp = mdLGTrans * (mdCXi * mdLG);
    mdTemp.InvertChol(&dLogDetLCL);

    for(int i=0; i<=nPlotPoints; i++) {
      for(int j=0; j<=nPlotPoints; j++) {
	vdXi[0] = double(vdPlotX[i]) - double(vdChi[nP1]);
	vdXi[1] = double(vdPlotY[j]) - double(vdChi[nP2]);
	vdTemp = mdTemp * vdXi;
	dTemp = 0;
	for(int k=0; k<vdXi.m_pnDimSize[0]; k++)
	  dTemp += double(vdXi[k])*double(vdTemp[k]);
	mdFunctionValuesT2[i][j] += (1.0 / sqrt(4*M_PI*M_PI)) * exp(-0.5*
	    (dLogDetLCL + dTemp));
      } // for j
    } // for i
    mdZT2 = mdFunctionValuesT2;
  } // if bMarPar

  PrintSuccess();

  // This should call for the destructors...
  delete[] pDat;
} // Calculate3DMCMCIntegrationIncMar


/* For every dataset in the ensemble, the ML sample is found, and the value of
 * the two given parameters is stored in the vectors vdX and vdY
 */
void Calculate2DEnsembleMLPoints(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY) {
  int nMCMCSteps;
  CMCMCDataEnsemble oDat;
  SParametersType oParameters;
  FILE *pFile;
  int nFileParameters, nFileTMParameters, nFileDataSets, nFileObservations;
  int nDataSetsLength;
  double *pdMaxLL=NULL;
  unsigned long nFileLength, nFilePos;
  char strMsg[160], strFileVersion[16];

  strcpy(strMsg, "Making 1D plots from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "' on the fly...");

  nMCMCSteps = oConstants.nMCMCSteps;

  try {
    // Open the ensemble data file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileDataSets, sizeof(int), 1, pFile) ) throw 5;
    oData.nDataSets = nFileDataSets;
    if(! fread(&nFileObservations, sizeof(int), 1, pFile) ) throw 5;
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(nFileParameters != oConstants.nParameters) {
      printf("WARNING: nFileParameters != nParameters\n");
      printf("nFileParameters: %i,  nParameters: %i\n", nFileParameters, oConstants.nParameters);
    } // if nFileParameters
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;
    if(nFileTMParameters != oConstants.nMarParameters) {
      printf("WARNING: nFileTMParameters != nTMParameters\n");
      printf("nFileTMParameters: %i,  nTMParameters: %i\n", nFileTMParameters, oConstants.nMarParameters);
    } // if nFileTMParameters
    nDataSetsLength = nFileObservations*nFileDataSets;
    oData.mdDataSets.Initialize(nFileObservations, nFileDataSets);
    if(! fread(oData.mdDataSets.m_pdData, sizeof(double), nDataSetsLength, pFile) ) throw 6;

    // Calculate nMCMCSteps
    nMCMCSteps = (nFileLength-sizeof(char)*16 - 4*sizeof(int) - nDataSetsLength*sizeof(double)) / (
        sizeof(double)*(
	  nFileParameters + nFileDataSets + 1 +
	  nFileTMParameters*nFileDataSets +
	  nFileTMParameters*nFileTMParameters));

    if(nMCMCSteps <= 0) throw 7;
    // Leave the file open 'till the end

    // This structure contains a single MCMC ensemble step
    oDat.Initialize(nFileDataSets, nFileParameters, nFileTMParameters);

    // Remember the index of the max LL value
    pdMaxLL = new double[nFileDataSets];

    vdX.Initialize(nFileDataSets);
    vdY.Initialize(nFileDataSets);


    // Figure the parameter plotting boundaries (but we actually want to
    // calculate these)
    ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

    // Walk through the whole chain to figure out the plotting boundaries
    // Save the current file position.
    nFilePos = ftell(pFile);

    fprintf(stderr, "First run through mcmc\r");
    for(int nStep=0; nStep<nMCMCSteps; nStep++) {
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;

      for(int i=0; i<nFileDataSets; i++) {
	if(nStep == 0 || pdMaxLL[i] > oDat.pdLogLik[i]) {
	  pdMaxLL[i] = oDat.pdLogLik[i];
	  vdX[i] = oDat.pdPar[nParameter1];
	  vdY[i] = oDat.pdPar[nParameter2];
	} // if nStep
      } // for i
      fprintf(stderr, "Run through mcmc: %5.2f%%\r", nStep*100.0/nMCMCSteps);
    } // for nStep
    fprintf(stderr, "Run through mcmc: done    \n");

    fclose(pFile);
  } catch(int nError) {
    if(nError != 1) fclose(pFile);
    PrintFailed();
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
  } // try

  if(pdMaxLL) delete[] pdMaxLL;
} // Calculate2DEnsembleMLPoints


/* For every marginalisation parameter p, and for every dataset in the ensemble
 * i, this function does the following. It marginalises over every parameter but
 * p, and it produces the marginalised posterior. The plot is saved in the
 * binary file strFileName.
 * */
void Calculate1DEnsembleIntegrationMarParameters(SDataType &oData, SConstantsType &oConstants, const char *strFileName) {
  int nMCMCSteps, nPlotPoints;
  CMCMCDataEnsemble oDat;
  SParametersType oParameters;
  FILE *pFile, *pFileOut;
  int nFileParameters, nFileTMParameters, nFileDataSets, nFileObservations;
  int nDataSetsLength;
  int nIndex;
  int *pnMaxLLIndex=NULL;
  double *pdMaxLL=NULL;
  double *pdMaxLLDiff=NULL;
  double **ppdDev=NULL;
  double **ppdDevSquare=NULL;
  double *pdRms=NULL;
  double dTemp, dMean, dMeanDevSquare;
  unsigned long nFileLength, nFilePos;
  char strMsg[160], strFileVersion[16];
  CMatrix mdMin, mdMax;
  CMatrix mdPlotX, mdPlotY, mdPlotYML, mdPlotYMLCorr;

  strcpy(strMsg, "Making 1D plots from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "' on the fly...");

  nPlotPoints = oConstants.nSmallPlotPoints;
  nMCMCSteps = oConstants.nMCMCSteps;

  try {
    // Open the ensemble data file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileDataSets, sizeof(int), 1, pFile) ) throw 5;
    oData.nDataSets = nFileDataSets;
    if(! fread(&nFileObservations, sizeof(int), 1, pFile) ) throw 5;
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(nFileParameters != oConstants.nParameters) {
      printf("WARNING: nFileParameters != nParameters\n");
      printf("nFileParameters: %i,  nParameters: %i\n", nFileParameters, oConstants.nParameters);
    } // if nFileParameters
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;
    if(nFileTMParameters != oConstants.nMarParameters) {
      printf("WARNING: nFileTMParameters != nTMParameters\n");
      printf("nFileTMParameters: %i,  nTMParameters: %i\n", nFileTMParameters, oConstants.nMarParameters);
    } // if nFileTMParameters
    nDataSetsLength = nFileObservations*nFileDataSets;
    oData.mdDataSets.Initialize(nFileObservations, nFileDataSets);
    if(! fread(oData.mdDataSets.m_pdData, sizeof(double), nDataSetsLength, pFile) ) throw 6;

    // Calculate nMCMCSteps
    nMCMCSteps = (nFileLength-sizeof(char)*16 - 4*sizeof(int) - nDataSetsLength*sizeof(double)) / (
        sizeof(double)*(
	  nFileParameters + nFileDataSets + 1 +
	  nFileTMParameters*nFileDataSets +
	  nFileTMParameters*nFileTMParameters));

    if(nMCMCSteps <= 0) throw 7;
    // Leave the file open 'till the end

    // This structure contains a single MCMC ensemble step
    oDat.Initialize(nFileDataSets, nFileParameters, nFileTMParameters);

    // Remember the index of the max LL value
    pnMaxLLIndex = new int[nFileDataSets];
    pdMaxLL = new double[nFileDataSets];
    pdMaxLLDiff = new double[nFileDataSets];
    ppdDev = new double*[nFileTMParameters];
    ppdDevSquare = new double*[nFileTMParameters];
    pdRms = new double[nFileTMParameters];

    for(int p=0; p<nFileTMParameters; p++) {
      ppdDev[p] = new double[nFileDataSets];
      ppdDevSquare[p] = new double[nFileDataSets];
      for(int i=0; i<nFileDataSets; i++) {
	ppdDev[p][i] = 0;
	ppdDevSquare[p][i] = 0;
      } // for i
      pdRms[p] = 0;
    } // for p

    // Figure the parameter plotting boundaries (but we actually want to
    // calculate these
    ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    mdMin.Initialize(nFileDataSets, nFileTMParameters);
    mdMax.Initialize(nFileDataSets, nFileTMParameters);

    // Walk through the whole chain to figure out the plotting boundaries
    // Save the current file position.
    nFilePos = ftell(pFile);

    fprintf(stderr, "First run through mcmc\r");
    for(int nStep=0; nStep<nMCMCSteps; nStep++) {
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;

      for(int i=0; i<nFileDataSets; i++) {
	for(int p=0; p<nFileTMParameters; p++) {
	  if(nStep == 0) {
	    mdMin[i][p] = oDat.pdChi[p+i*nFileTMParameters] - 4*sqrt(oDat.pdCXiInv[p+p*nFileTMParameters]);
	    mdMax[i][p] = oDat.pdChi[p+i*nFileTMParameters] + 4*sqrt(oDat.pdCXiInv[p+p*nFileTMParameters]);
	  } else {
	    if(mdMin[i][p] > oDat.pdChi[p+i*nFileTMParameters] - 4*sqrt(oDat.pdCXiInv[p+p*nFileTMParameters]))
	      mdMin[i][p] = oDat.pdChi[p+i*nFileTMParameters] - 4*sqrt(oDat.pdCXiInv[p+p*nFileTMParameters]);
	    if(mdMax[i][p] < oDat.pdChi[p+i*nFileTMParameters] + 4*sqrt(oDat.pdCXiInv[p+p*nFileTMParameters]))
	      mdMax[i][p] = oDat.pdChi[p+i*nFileTMParameters] + 4*sqrt(oDat.pdCXiInv[p+p*nFileTMParameters]);
	  } // if nStep
	} // for p

	if(nStep == 0 || pdMaxLL[i] > oDat.pdLogLik[i]) {
	  pnMaxLLIndex[i] = nStep;
	  pdMaxLL[i] = oDat.pdLogLik[i];
	  for(int p=0; p<nFileTMParameters; p++) {
	    ppdDev[p][i] = oDat.pdChi[p+i*nFileTMParameters];
	    ppdDevSquare[p][i] = oDat.pdChi[p+i*nFileTMParameters] * oDat.pdChi[p+i*nFileTMParameters];
	  } // for p
	} // if nStep
	if(nStep == 0 || pdMaxLLDiff[i] < oDat.dLogLik - oDat.pdLogLik[i]) {
	  pdMaxLLDiff[i] = oDat.dLogLik - oDat.pdLogLik[i];
	} // if nStep
      } // for i
      fprintf(stderr, "First run through mcmc: %5.2f%%\r", nStep*100.0/nMCMCSteps);
    } // for nStep
    fprintf(stderr, "First run through mcmc: done    \n");
    // Set the file position back to where we were
    if(fseek(pFile, nFilePos, SEEK_SET)) throw 3;

    // Set the x-axis values of our plots
    mdPlotX.Initialize(nPlotPoints+1, nFileDataSets*nFileTMParameters);
    mdPlotY.Initialize(nPlotPoints+1, nFileDataSets*nFileTMParameters);
    mdPlotYML.Initialize(nPlotPoints+1, nFileDataSets*nFileTMParameters);
    mdPlotYMLCorr.Initialize(nPlotPoints+1, nFileDataSets*nFileTMParameters);
    for(int i=0; i<nFileDataSets; i++) {
      for(int p=0; p<nFileTMParameters; p++) {
	for(int k=0; k<nPlotPoints+1; k++) {
	  mdPlotX[k][i*nFileTMParameters+p] = double(mdMin[i][p]) + k*((double(mdMax[i][p]) - double(mdMin[i][p]))/nPlotPoints);
	  mdPlotY[k][i*nFileTMParameters+p] = 0;
	  mdPlotYML[k][i*nFileTMParameters+p] = 0;
	  mdPlotYMLCorr[k][i*nFileTMParameters+p] = 0;
	} // for k
      } // for p
    } // for i

    // Calculate the rms of the parameter estimates
    for(int p=0; p<nFileTMParameters; p++) {
      dMean=0; dMeanDevSquare=0;
      for(int i=0; i<nFileDataSets; i++) {
	dMean += ppdDev[p][i];
	dMeanDevSquare += ppdDevSquare[p][i];
      } // for i
      dMean /= nFileDataSets;
      dMeanDevSquare /= nFileDataSets;
      pdRms[p] = sqrt(dMeanDevSquare - dMean*dMean);
    } // for p


    fprintf(stderr, "Second run through mcmc\r");
    // With the x-values available, we can start calculating the y-values
    for(int nStep=0; nStep<nMCMCSteps; nStep++) {
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;

      for(int i=0; i<nFileDataSets; i++) {
	for(int p=0; p<nFileTMParameters; p++) {
	  for(int k=0; k<nPlotPoints+1; k++) {
	    mdPlotY[k][i*nFileTMParameters+p] +=
	      exp(oDat.dLogLik - oDat.pdLogLik[i] - pdMaxLLDiff[i]) *
	      (1.0 / sqrt(2 * M_PI * oDat.pdCXiInv[p+p*nFileTMParameters])) *
	      exp( -0.5 * 
		  (double(mdPlotX[k][i*nFileTMParameters+p]) - oDat.pdChi[p+i*nFileTMParameters]) *
		  (double(mdPlotX[k][i*nFileTMParameters+p]) - oDat.pdChi[p+i*nFileTMParameters]) / oDat.pdCXiInv[p+p*nFileTMParameters]);
	  } // for k
	} // for p

	// If necessary, also calculate the ML plot
	if(pnMaxLLIndex[i] == nStep) {
	  for(int p=0; p<nFileTMParameters; p++) {
	    for(int k=0; k<nPlotPoints+1; k++) {
	      mdPlotYML[k][i*nFileTMParameters+p] =
		(1.0 / sqrt(2 * M_PI * oDat.pdCXiInv[p+p*nFileTMParameters])) *
		exp( -0.5 * 
		    (double(mdPlotX[k][i*nFileTMParameters+p]) - oDat.pdChi[p+i*nFileTMParameters]) *
		    (double(mdPlotX[k][i*nFileTMParameters+p]) - oDat.pdChi[p+i*nFileTMParameters]) / oDat.pdCXiInv[p+p*nFileTMParameters]);
	      mdPlotYMLCorr[k][i*nFileTMParameters+p] =
		(1.0 / (pdRms[p]*sqrt(2 * M_PI))) *
		exp( -0.5 * 
		    (double(mdPlotX[k][i*nFileTMParameters+p]) - oDat.pdChi[p+i*nFileTMParameters]) *
		    (double(mdPlotX[k][i*nFileTMParameters+p]) - oDat.pdChi[p+i*nFileTMParameters]) / (pdRms[p]*pdRms[p]));
	    } // for k
	  } // for p
	} // if pnMaxLLIndex
      } // for i
      fprintf(stderr, "Second run through mcmc: %5.2f%%\r", nStep*100.0/nMCMCSteps);
    } // for nStep
    fprintf(stderr, "Second run through mcmc: done    \n");
    fclose(pFile);
  } catch(int nError) {
    if(nError != 1) fclose(pFile);
    PrintFailed();
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
  } // try

  try {
    /*
    fprintf(stderr, "\n");
    for(int k=0; k<=nPlotPoints; k++) {
      fprintf(stderr, "%20.17e    %20.17e    %20.17e\n", double(mdPlotX[k][0]), double(mdPlotY[k][0]), double(mdPlotYML[k][0]));
    } // for k
    fprintf(stderr, "\n");
    */

    PrintStatus("Saving plots to disk...");
    nIndex=0;
    nPlotPoints++;
    dTemp = 0;

    if(! (pFileOut = fopen("ensemblemarplots.dat", "wb")) ) throw 1;

    if(! fwrite(&nFileDataSets , sizeof(int) , 1 , pFileOut) ) throw 2;
    if(! fwrite(&nFileTMParameters , sizeof(int) , 1 , pFileOut) ) throw 3;
    if(! fwrite(&nPlotPoints , sizeof(int) , 1 , pFileOut) ) throw 4;

    // Also write the true values
    for(int p=0; p<nFileTMParameters; p++) {
      if(! fwrite(&dTemp , sizeof(double) , 1 , pFileOut) ) throw 4;
    } // for p

    for(int i=0; i<nFileDataSets; i++) {
      for(int p=0; p<nFileTMParameters; p++) {
	if(! fwrite(mdPlotX.m_pdData + nIndex*nPlotPoints, sizeof(double) , nPlotPoints , pFileOut) ) throw 5;
	if(! fwrite(mdPlotY.m_pdData + nIndex*nPlotPoints, sizeof(double) , nPlotPoints , pFileOut) ) throw 6;
	if(! fwrite(mdPlotYML.m_pdData + nIndex*nPlotPoints, sizeof(double) , nPlotPoints , pFileOut) ) throw 7;
	if(! fwrite(mdPlotYMLCorr.m_pdData + nIndex*nPlotPoints, sizeof(double) , nPlotPoints , pFileOut) ) throw 8;

	nIndex++;
      } // for p
    } // for i

    fclose(pFileOut);
    PrintSuccess();
  } catch(int nError) {
    if(nError != 1) fclose(pFileOut);
    PrintFailed();
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
  } // try

  if(pnMaxLLIndex) delete[] pnMaxLLIndex;
  if(pdMaxLL) delete[] pdMaxLL;
  if(pdMaxLLDiff) delete[] pdMaxLLDiff;
  for(int p=0; p<nFileTMParameters; p++) {
    delete[] ppdDev[p];
    delete[] ppdDevSquare[p];
  } // for p
  delete[] ppdDev;
  delete[] ppdDevSquare;
  if(pdRms) delete[] pdRms;
} // Calculate1DEnsembleIntegrationMarParameters


/* For every MCMC parameter p, and for every dataset in the ensemble i, this
 * function does the following. It marginalises over every parameter but p, and
 * it produces the marginalised posterior. The plot is saved in the binary file
 * strFileName.
 * */
void Calculate1DEnsembleIntegrationMCMCParameters(SDataType &oData, SConstantsType &oConstants, const char *strFileName) {
  int nMCMCSteps, nPlotPoints;
  CMCMCDataEnsemble oDat;
  SParametersType oParameters;
  FILE *pFile, *pFileOut;
  int nFileParameters, nFileTMParameters, nFileDataSets, nFileObservations;
  int nDataSetsLength;
  int nIndex;
  int *pnMaxLLIndex=NULL;
  int *pnParameterIndex=NULL;
  int nMCMCParameters=0;
  double *pdMaxLL=NULL;
  double *pdMaxLLDiff=NULL;
  unsigned long nFileLength, nFilePos;
  char strMsg[160], strFileVersion[16];
  CVector vdMin, vdMax;
  CMatrix mdPlotX, mdPlotY;

  strcpy(strMsg, "Making 1D plots from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "' on the fly...");

  nPlotPoints = oConstants.nSmallPlotPoints;
  nMCMCSteps = oConstants.nMCMCSteps;

  // Get the varying parameters, and indices from parameters.conf, not from the
  // datafile. That is too involved, and we're rewriting all code anyways
  nMCMCParameters = NumberOfVaryingParameters(oConstants);
  pnParameterIndex = new int[nMCMCParameters];

  // Set the array of indices that allows quick reference in parameter structs
  SetParameterIndexArray(oConstants, pnParameterIndex);

  try {
    // Open the ensemble data file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileDataSets, sizeof(int), 1, pFile) ) throw 5;
    oData.nDataSets = nFileDataSets;
    if(! fread(&nFileObservations, sizeof(int), 1, pFile) ) throw 5;
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(nFileParameters != oConstants.nParameters) {
      printf("WARNING: nFileParameters != nParameters\n");
      printf("nFileParameters: %i,  nParameters: %i\n", nFileParameters, oConstants.nParameters);
    } // if nFileParameters
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;
    if(nFileTMParameters != oConstants.nMarParameters) {
      printf("WARNING: nFileTMParameters != nTMParameters\n");
      printf("nFileTMParameters: %i,  nTMParameters: %i\n", nFileTMParameters, oConstants.nMarParameters);
    } // if nFileTMParameters
    nDataSetsLength = nFileObservations*nFileDataSets;
    oData.mdDataSets.Initialize(nFileObservations, nFileDataSets);
    if(! fread(oData.mdDataSets.m_pdData, sizeof(double), nDataSetsLength, pFile) ) throw 6;

    // Calculate nMCMCSteps
    nMCMCSteps = (nFileLength-sizeof(char)*16 - 4*sizeof(int) - nDataSetsLength*sizeof(double)) / (
        sizeof(double)*(
	  nFileParameters + nFileDataSets + 1 +
	  nFileTMParameters*nFileDataSets +
	  nFileTMParameters*nFileTMParameters));

    if(nMCMCSteps <= 0) throw 7;
    // Leave the file open 'till the end

    // This structure contains a single MCMC ensemble step
    oDat.Initialize(nFileDataSets, nFileParameters, nFileTMParameters);

    // Remember the index of the max LL value
    pnMaxLLIndex = new int[nFileDataSets];
    pdMaxLL = new double[nFileDataSets];
    pdMaxLLDiff = new double[nFileDataSets];

    // Figure the parameter plotting boundaries (but we actually want to
    // calculate these
    ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    vdMin.Initialize(nMCMCParameters);
    vdMax.Initialize(nMCMCParameters);

    // Set the min and max for all parameters (No, just the varying ones)
    for(int p=0; p<nMCMCParameters; p++) {
      vdMin[p] = oParameters.pdParMinBound[pnParameterIndex[p]];
      vdMax[p] = oParameters.pdParMaxBound[pnParameterIndex[p]];
    } // for p

    // Set the x-axis values of our plots
    mdPlotX.Initialize(nPlotPoints+1, nFileDataSets*nMCMCParameters);
    mdPlotY.Initialize(nPlotPoints+1, nFileDataSets*nMCMCParameters);
    for(int i=0; i<nFileDataSets; i++) {
      for(int p=0; p<nMCMCParameters; p++) {
	for(int k=0; k<nPlotPoints+1; k++) {
	  mdPlotX[k][i*nMCMCParameters+p] = double(vdMin[p]) + k*((double(vdMax[p]) - double(vdMin[p]))/nPlotPoints);
	  mdPlotY[k][i*nMCMCParameters+p] = 0;
	} // for k
      } // for p
    } // for i

    // Walk through the whole chain to figure out the offset of the exponent
    nFilePos = ftell(pFile);

    fprintf(stderr, "First run through mcmc\r");
    for(int nStep=0; nStep<nMCMCSteps; nStep++) {
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;

      for(int i=0; i<nFileDataSets; i++) {
	if(nStep == 0 || pdMaxLLDiff[i] < oDat.dLogLik - oDat.pdLogLik[i]) {
	  pdMaxLLDiff[i] = oDat.dLogLik - oDat.pdLogLik[i];
	} // if nStep
      } // for i
      fprintf(stderr, "First run through mcmc: %5.2f%%\r", nStep*100.0/nMCMCSteps);
    } // for nStep
    fprintf(stderr, "First run through mcmc: done    \n");
    // Set the file position back to where we were
    if(fseek(pFile, nFilePos, SEEK_SET)) throw 3;



    fprintf(stderr, "Second run through mcmc\r");
    // With the x-values available, we can start calculating the y-values
    for(int nStep=0; nStep<nMCMCSteps; nStep++) {
      if(! fread(oDat.pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(oDat.pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(oDat.dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(oDat.pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(oDat.pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;

      for(int i=0; i<nFileDataSets; i++) {
	for(int p=0; p<nMCMCParameters; p++) {
	  nIndex = int(nPlotPoints*
	      (oDat.pdPar[pnParameterIndex[p]] - double(vdMin[p])) /
	      (double(vdMax[p]) - double(vdMin[p])));
	  if(nIndex>=0 && nIndex<=nPlotPoints) {
	    mdPlotY[nIndex][i*nMCMCParameters+p] += exp(oDat.dLogLik - oDat.pdLogLik[i] - pdMaxLLDiff[i]);
	  } // if nIndex
	} // for p

      } // for i
      fprintf(stderr, "Second run through mcmc: %5.2f%%\r", nStep*100.0/nMCMCSteps);
    } // for nStep
    fprintf(stderr, "Second run through mcmc: done    \n");
    fclose(pFile);
  } catch(int nError) {
    if(nError != 1) fclose(pFile);
    PrintFailed();
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
  } // try

  try {
    PrintStatus("Saving plots to disk...");
    nIndex=0;
    nPlotPoints++;

    if(! (pFileOut = fopen("ensemblemcmcplots.dat", "wb")) ) throw 1;

    if(! fwrite(&nFileDataSets , sizeof(int) , 1 , pFileOut) ) throw 2;
    if(! fwrite(&nMCMCParameters , sizeof(int) , 1 , pFileOut) ) throw 3;
    if(! fwrite(&nPlotPoints , sizeof(int) , 1 , pFileOut) ) throw 4;

    // Also write the true values
    for(int p=0; p<nMCMCParameters; p++) {
      if(! fwrite(&(oParameters.pdPar[pnParameterIndex[p]]) , sizeof(double) , 1 , pFileOut) ) throw 4;
    } // for p

    for(int i=0; i<nFileDataSets; i++) {
      for(int p=0; p<nMCMCParameters; p++) {
	if(! fwrite(mdPlotX.m_pdData + nIndex*nPlotPoints, sizeof(double) , nPlotPoints , pFileOut) ) throw 5;
	if(! fwrite(mdPlotY.m_pdData + nIndex*nPlotPoints, sizeof(double) , nPlotPoints , pFileOut) ) throw 6;

	nIndex++;
      } // for p
    } // for i

    fclose(pFileOut);
    PrintSuccess();
  } catch(int nError) {
    if(nError != 1) fclose(pFileOut);
    PrintFailed();
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
  } // try

  if(pnParameterIndex) delete[] pnParameterIndex;
  if(pnMaxLLIndex) delete[] pnMaxLLIndex;
  if(pdMaxLL) delete[] pdMaxLL;
  if(pdMaxLLDiff) delete[] pdMaxLLDiff;
} // Calculate1DEnsembleIntegrationMCMCParameters



/* This function performs the actual MCMC integration and returns the result in
 * the vdX, vdY and mdZ objects
 * */
void Calculate3DMCMCIntegrationSkymap(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ) {
  int nMCMCSteps, nPlotPointsX, nPlotPointsY, nRandIndex;
  int nSourceNumber, nParameterNumber;
  double dEvidence;
//  int nParameter1, nParameter2;
  CMCMCDataPoint *pDat;

  nPlotPointsX = 90-1;
  nPlotPointsY = 45-1;
  nPlotPointsX = 30-1;
  nPlotPointsY = 15-1;
  nMCMCSteps = oConstants.nMCMCSteps;

//  pDat = new MCMCDataPoint[nMCMCSteps];
  if(! ReadMCMCDataFile(nMCMCSteps, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName) ) {
    PrintFailed();
    return;
  }

  printf("Using integration numbers: %i,%i\n", nParameter1, nParameter2);
  PrintStatus("Integrating datapoints and generating plot...");

  // Prepare for integrating: initialize the variables
  CVector vdPlotX, vdPlotY, vdFunctionValues, vdTemp, vdTemp2, vdTemp3, vdMeanValues, vdMeanSquareValues, vdSigma;
  CMatrix mdFunctionValues;
  int nIndexX, nIndexY;
  vdPlotX.Initialize(nPlotPointsX + 1);
  vdPlotY.Initialize(nPlotPointsY + 1);
  mdFunctionValues.Initialize(nPlotPointsX + 1, nPlotPointsY + 1);
  vdFunctionValues = vdPlotX;

  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

  for(int i=0; i<=nPlotPointsX; i++) {
    vdPlotX[i] = oParameters.pdParMinBound[nParameter1] +
      i*((oParameters.pdParMaxBound[nParameter1] -
	    oParameters.pdParMinBound[nParameter1])/nPlotPointsX);
  } // for i

  for(int i=0; i<=nPlotPointsY; i++) {
    vdPlotY[i] = oParameters.pdParMinBound[nParameter2] +
      i*((oParameters.pdParMaxBound[nParameter2] -
	    oParameters.pdParMinBound[nParameter2])/nPlotPointsY);
  } // for i

  for(int i=0; i<=nPlotPointsX; i++) {
    vdFunctionValues[i] = 0;
    for(int j=0; j<= nPlotPointsY; j++)
      mdFunctionValues[i][j] = 0;
  } // for i
  vdTemp = vdFunctionValues;
  vdSigma = vdFunctionValues;
  vdMeanValues = vdFunctionValues;
  vdMeanSquareValues = vdFunctionValues;

  nIndexX = int((nPlotPointsX-1)*(pDat[0].pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
  nIndexY = int((nPlotPointsY-1)*(pDat[0].pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

  if(nIndexX >= 0 && nIndexX <= nPlotPointsX && nIndexY >= 0 && nIndexY <= nPlotPointsY) {
    vdFunctionValues[nIndexX] += 1;// else {printf("%i ", nIndex);
    mdFunctionValues[nIndexX][nIndexY] += 1;
  } // if nIndex

  // Now integrate the values! Work through all computer data //
  for(int nStep=1; nStep<nMCMCSteps; nStep++) {
    nIndexX = int((nPlotPointsX-1)*(pDat[nStep].pdPar[nParameter1] - oParameters.pdParMinBound[nParameter1])/(oParameters.pdParMaxBound[nParameter1]- oParameters.pdParMinBound[nParameter1]));
    nIndexY = int((nPlotPointsY-1)*(pDat[nStep].pdPar[nParameter2] - oParameters.pdParMinBound[nParameter2])/(oParameters.pdParMaxBound[nParameter2]- oParameters.pdParMinBound[nParameter2]));

    if(nIndexX >= 0 && nIndexX <= nPlotPointsX && nIndexY >= 0 && nIndexY <= nPlotPointsY) {
      vdFunctionValues[nIndexX] += 1;
      mdFunctionValues[nIndexX][nIndexY] += 1;
    } // if nIndex

  } // for nStep

  vdX = vdPlotX;
  vdY = vdPlotY;
  mdZ = mdFunctionValues;

  PrintSuccess();

  // This should call for the destructors...
  delete[] pDat;
} // Calculate3DMCMCIntegrationSkymap


/* Integrate the MCMC data over all parameters, except 2 particular
 * parameters: -1 keeps the GWB parameters
 *              n keeps the n'th pulsar noise paramters
 *
 * Last updated: 2007-08-23 */
//void Integrate3DMCMCData(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nPulsarNumber, int nUseParameters) {
void Integrate3DMCMCData(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2) {
  CVector vdX, vdY, vdZ;
  CMatrix mdZ;
  CVector vdTemp, vdTemp2;
  CMatrix mdTemp;
  char strBuf[160];
  int nPlotPoints;
  int i, j;

  Calculate3DMCMCIntegration(oData, oParameters, oConstants, strFileName, nParameter1, nParameter2, vdX, vdY, mdZ);

  // Writa data to file so it can be used in with gnuplot
  strcpy(strBuf, oConstants.strDataDir);
  if(strBuf[strlen(strBuf)-1] != '/') strcat(strBuf, "/");
  strcat(strBuf, "plotmcmcdata.txt");

  nPlotPoints = vdX.m_pnDimSize[0]-1;
  vdTemp = vdX;
  vdTemp2 = vdY;
  mdTemp = mdZ;

  vdX.Initialize( (nPlotPoints + 1)*(nPlotPoints + 1) );
  vdY.Initialize( (nPlotPoints + 1)*(nPlotPoints + 1) );
  vdZ.Initialize( (nPlotPoints + 1)*(nPlotPoints + 1) );

  for(i=0; i<=nPlotPoints; i++) {
    for(j=0; j<=nPlotPoints; j++) {
      vdX[i+(nPlotPoints+1)*j] = double(vdTemp[i]);
      vdY[i+(nPlotPoints+1)*j] = double(vdTemp2[j]);
      vdZ[i+(nPlotPoints+1)*j] = double(mdTemp[i][j]);
    } // for j
  } // for i

  WritePlot(strBuf, vdX, vdY, vdZ);
} // Integrate3DMCMCData


/* This function checks how many points are in the ellipsoid region defined by
 * mdC^{-1} = mdInvC. It returns the enlargement factor dEnlargement with which
 * to enlarge the covariance matrix such that 30% of the points are in the
 * region.
 * */
double CalcEnlargement(int nDimensions, int nMCMCSteps, double **ppdX, double *pdZ, double *pdMu, CMatrix &mdInvC, bool bPrintResults) {
  const double dFraction = 0.2;
  double dReturnValue, *pdDistance;
  int nIndex;
  CNumber ndDistance;
  CVector vdTheta;

  vdTheta.Initialize(nDimensions);
  pdDistance = new double[nMCMCSteps];

  for(int i=0; i<nMCMCSteps; i++) {
    for(int d=0; d<nDimensions; d++) {
      vdTheta[d] = ppdX[i][d] - pdMu[d];
    } // for d
    ndDistance = vdTheta * mdInvC * vdTheta;
    pdDistance[i] = double(ndDistance);
  } // for i

  // Sort this list
  QuickSort(pdDistance, 0, nMCMCSteps-1);

  nIndex = int(nMCMCSteps*dFraction);
  dReturnValue = 0.5*(pdDistance[nIndex] + pdDistance[nIndex+1]);
  delete[] pdDistance;

  return dReturnValue;
} // CalcEnlargement

/* This function checks how many points are in the ellipsoid region defined by
 * mdC^{-1} = mdInvC. It returns that number;
 * */
int PointsInEllipsoid(int nDimensions, int nMCMCSteps, double **ppdX, double *pdZ, double *pdMu, CMatrix &mdInvC, double dEnlargement) {
  int nNumberOfPoints;
  CNumber ndDistance;
  CVector vdTheta;

  vdTheta.Initialize(nDimensions);
  nNumberOfPoints = 0;
  for(int i=0; i<nMCMCSteps; i++) {
    for(int d=0; d<nDimensions; d++) {
      vdTheta[d] = ppdX[i][d] - pdMu[d];
    } // for d
    ndDistance = vdTheta * mdInvC * vdTheta;
    if(double(ndDistance) < dEnlargement)
      nNumberOfPoints++;
  } // for i
  return nNumberOfPoints;
} // PointsInEllipsoid

void CalcCovariance(int nDimensions, int nMCMCSteps, double **ppdX, double *pdZ, double *pdMu, CMatrix &mdInvC, double &dEnlargement, double &dLogDetC, bool bPrintResults) {
  int nCovariancePoints;
  double dCentreMassWeight, dZMu, dZMax;
  // Determine the best fitting Gaussian covariance matrix for this set of
  // datapoints. First sort the whole set according to Z-value (likelihood)
  //
  // First use the entire dataset. We are using a Gaussian anyways.
  // Select the top 20% for determining the covariance matrix and the centre
  nCovariancePoints = 0.2 * nMCMCSteps;

  dZMu = 0;
  dZMax = pdZ[nMCMCSteps-1];
  for(int i=nMCMCSteps - nCovariancePoints; i<nMCMCSteps; i++)
    dZMu += pdZ[i];
  dZMu /= nCovariancePoints;

  dCentreMassWeight = 0.0;
  for(int i=0; i<nDimensions; i++) {
    pdMu[i] = 0.0;
    for(int j=0; j<nDimensions; j++) {
      mdInvC[i][j] = 0.0;
    } // for j
  } // for i

  for(int i=nMCMCSteps - nCovariancePoints; i<nMCMCSteps; i++) {
    dCentreMassWeight += 1; //exp(pdZ[i] - dZMax);
    for(int a=0; a<nDimensions; a++) {
      pdMu[a] += ppdX[i][a] * 1; //exp(pdZ[i]-dZMax);
    } // for a
  } // for i

  for(int a=0; a<nDimensions; a++)
    pdMu[a] /= dCentreMassWeight;

  for(int i=nMCMCSteps - nCovariancePoints; i<nMCMCSteps; i++) {
    for(int b=0; b<nDimensions; b++) {
      for(int a=0; a<nDimensions; a++) {
	mdInvC[a][b] += (ppdX[i][a] - pdMu[a])*(ppdX[i][b] - pdMu[b]);
      } // for b
    } // for a
  } // for i
  mdInvC /= nCovariancePoints;

  // Regularise the covariance matrix here. This is probably not stable enough
  // for volume calculations
  try {
    mdInvC.InvertChol(&dLogDetC);
  } catch (ELinearError err) {  // Error handling
    switch(err) {
    case ELENotDefined:
      printf("Not Defined\n");
      break;
    case ELEWrongClassType:
      printf("Wrong Class Type\n");
      break;
    case ELEBadIndex:
      printf("Wrong Index number\n");
      break;
    case ELEDimensionMisMatch:
      printf("Matrix Dimensions Do Not Match\n");
      break;
    case ELELapack:
      printf("Lapack error!\n");
      break;
    default:
      printf("Default error!\n");
      break;
    } // switch
  } // try

  // Calculate the enlargement factor of the covariance matrix
  dEnlargement = CalcEnlargement(nDimensions, nMCMCSteps, ppdX, pdZ, pdMu, mdInvC, bPrintResults);
} // CalcCovariance


/* Calculate the integral value, given the MCMC chain
 * */
double Integral(int nDimensions, int nMCMCSteps, double **ppdX, double *pdZ, double *pdMu, double dEnlargement, CMatrix &mdInvC, double dLogDetC, bool bPrintResults) {
  double dPriorVolume, dTotalWeight, dSigmaWeight,
	 dCoefficient, dLogCoefficient;
  int nNumberOfPoints;
  double dVolumeEllipsoid, dLogVolumeEllipsoid;
  double dRadiusSquared, dFactorial;
  double dIntegral;
  double dZMu, dZMin;
  CVector vdTheta;
  CNumber ndDistance;

  // Calculate the prior volume
  dPriorVolume = 0;
  dTotalWeight = 0;
  dSigmaWeight = 0;
  dCoefficient = 0;

  // We know the volume of the (approx) 1-sigma volume
  dFactorial = pow(M_PI, double(nDimensions)/2.0) / gsl_sf_gamma(1+double(nDimensions)/2.0);
  dVolumeEllipsoid = pow(dEnlargement, double(nDimensions)/2.0) * exp(0.5*dLogDetC) * dFactorial;
  dLogVolumeEllipsoid = log(dEnlargement)*double(nDimensions)/2.0 + 0.5*dLogDetC + log(dFactorial);

  nNumberOfPoints = PointsInEllipsoid(nDimensions, nMCMCSteps, ppdX, pdZ, pdMu, mdInvC, dEnlargement);
  if(nNumberOfPoints >= 0 && nNumberOfPoints <= nMCMCSteps) {
    dZMin = pdZ[nMCMCSteps - nNumberOfPoints];
  } else {
    dZMin = pdZ[0];
  } // if nNumberOfPoints

  nNumberOfPoints = 0;
  vdTheta.Initialize(nDimensions);
  for(int i=0; i<nMCMCSteps; i++) {
    dTotalWeight += exp(dZMin - pdZ[i]);
    for(int d=0; d<nDimensions; d++)
      vdTheta[d] = ppdX[i][d] - pdMu[d];
    ndDistance = vdTheta * mdInvC * vdTheta;
    dRadiusSquared = double(ndDistance);
    if(dRadiusSquared < dEnlargement) {
      // We are in the region
      dSigmaWeight += exp(dZMin - pdZ[i]);
      nNumberOfPoints++;
    } // if function
  } // for i

  dCoefficient = dVolumeEllipsoid / dSigmaWeight;
  dLogCoefficient = dLogVolumeEllipsoid - log(dSigmaWeight);
  dPriorVolume = dCoefficient * dTotalWeight;

  // Now perform the integral
  dIntegral = log(nMCMCSteps) + dLogCoefficient + dZMin;

  if(bPrintResults) {
    printf("Prior volume:           %15.13e\n", dPriorVolume);
    printf("Coef. Est. (CE) volume: %15.13e\n", dVolumeEllipsoid);
    printf("dSigmaWeight:           %15.13e\n", dSigmaWeight);
    printf("Coefficient:            %15.13e\n", dCoefficient);
    printf("Points in CE volume:    %i / %i\n", nNumberOfPoints, nMCMCSteps);
  } // if bPrintResults

  return dIntegral;
} // Integral

/* This function calculates the Bayesian evidence from a MCMC chain. This
 * function does not yet check for the right model: even the amount of fitted
 * parameters in the mcmc file can differ from parameters.conf.
 * */
double MCMCEvidence(SDataType &oData, SConstantsType &oConstants, CMCMCDataPoint *pDat, int nMCMCSteps, int nParameters) {
  int nFitParameters = oConstants.nFitParameters;
  int nElementIndex, nBootstraps = 10;
  int nDimensions = nFitParameters;
  int *pnParameterIndices, nIndexPar, nIndexFit;
  double *pdZ, *pdBootstrapZ;
  double *pdPlotX, *pdPlotY;
  double **ppdX, **ppdBootstrapX;
  double dLogEvidence;
  int nNumberOfPoints;
  double dLogEvidenceRMS, dLogEvidenceMean, *pdLogEvidence;
  double *pdMu, dEnlargement;
  CMatrix mdC;
  double dLogDetC;

  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, (unsigned int)time(NULL));

  // The covariance matrix
  mdC.Initialize(nDimensions, nDimensions);

  // Allocate some memory
  pdMu = new double[nDimensions];

  pdLogEvidence = new double[nBootstraps];
  ppdX = new double *[nMCMCSteps];
  ppdBootstrapX = new double *[nMCMCSteps];
  pdZ = new double[nMCMCSteps];
  pdBootstrapZ = new double[nMCMCSteps];
  pnParameterIndices = new int[nDimensions];
  for(int i=0; i<nMCMCSteps; i++) {
    ppdX[i] = new double[nDimensions];
    ppdBootstrapX[i] = new double[nDimensions];
  } // for i

  // We should only work with the fitted parameters
  nIndexFit = 0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      nIndexPar = oConstants.poSources[s].nFirstParIndex + p;
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	pnParameterIndices[nIndexFit] = nIndexPar;
	nIndexFit++;
      } // if FITPARAMETER
    } // for p
  } // for s

  // Fill all values
  for(int i=0; i<nMCMCSteps; i++) {
    for(int d=0; d<nDimensions; d++) {
      ppdX[i][d] = pDat[i].pdPar[pnParameterIndices[d]];
      pdZ[i] = -pDat[i].dLogLik;
    } // for d
  } // for i

  QuickSort(ppdX, pdZ, nDimensions, 0, nMCMCSteps-1);

  CalcCovariance(nDimensions, nMCMCSteps, ppdX, pdZ, pdMu, mdC, dEnlargement, dLogDetC, false);
  nNumberOfPoints = PointsInEllipsoid(nDimensions, nMCMCSteps, ppdX, pdZ, pdMu, mdC, dEnlargement);
  dLogEvidence = Integral(nDimensions, nMCMCSteps, ppdX, pdZ, pdMu, dEnlargement, mdC, dLogDetC, false);
  printf("dLogEvidence = %18.16e\n", dLogEvidence);

  for(int i=0; i<nBootstraps; i++) {
    for(int j=0; j<nMCMCSteps; j++) {
      nElementIndex = gsl_rng_uniform_int(rng, nMCMCSteps);
      for(int k=0; k<nDimensions; k++) {
	ppdBootstrapX[j][k] = ppdX[nElementIndex][k];
      } // for k
      pdBootstrapZ[j] = pdZ[nElementIndex];
    } // for j

    QuickSort(ppdBootstrapX, pdBootstrapZ, nDimensions, 0, nMCMCSteps-1);
    CalcCovariance(nDimensions, nMCMCSteps, ppdBootstrapX, pdBootstrapZ, pdMu, mdC, dEnlargement, dLogDetC, false);
    pdLogEvidence[i] = Integral(nDimensions, nMCMCSteps, ppdBootstrapX, pdBootstrapZ, pdMu, dEnlargement, mdC, dLogDetC, false);
  } // for i

  // Now calculate the bootstrap error
  for(int i=0; i<nBootstraps; i++) {
    dLogEvidenceMean += pdLogEvidence[i];
  } // for i
  dLogEvidenceMean /= nBootstraps;
  for(int i=0; i<nBootstraps; i++) {
    dLogEvidenceRMS += (pdLogEvidence[i] - dLogEvidenceMean) * (pdLogEvidence[i] - dLogEvidenceMean);
  } // for i
  dLogEvidenceRMS = sqrt(dLogEvidenceRMS )/ nBootstraps;
  printf("dEvidence = exp(%18.16e +/- %18.16e%%)\n", dLogEvidence, dLogEvidenceRMS);

  // Free some memory
  delete[] pdMu;

  delete[] pdLogEvidence;
  delete[] pdZ;
  delete[] pdBootstrapZ;
  delete[] pnParameterIndices;
  for(int i=0; i<nMCMCSteps; i++) {
    delete[] ppdX[i];
    delete[] ppdBootstrapX[i];
  } // for i
  delete[] ppdX;
  delete[] ppdBootstrapX;

  return dLogEvidence;
} // MCMCEvidence

/* This function performs the actual MCMC integration and returns the result in
 * the vdX, vdY and mdZ objects
 * */
void CalculateMCMCEvidence(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName) {
  int nMCMCSteps, nPlotPoints, nRandIndex;
  int nSourceNumber, nParameterNumber;
  double dEvidence;

  CMCMCDataPoint *pDat;

  nPlotPoints = int(sqrt(double(oConstants.nPlotPoints)));
  nMCMCSteps = oConstants.nMCMCSteps;

//  pDat = new MCMCDataPoint[nMCMCSteps];
  if(! ReadMCMCDataFile(nMCMCSteps, oConstants.nParameters, oConstants.nMarParameters, pDat, strFileName) ) {
    PrintFailed();
    return;
  }

  dEvidence = MCMCEvidence(oData, oConstants, pDat, nMCMCSteps, oConstants.nParameters);

  // This should call for the destructors...
  delete[] pDat;
} // CalculateMCMCEvidence

/* Create clock corrected residuals, and write those to file
 */
void CreateClockCorrectedResiduals(SDataType &oData, SConstantsType &oConstants) {
#if 0
  SParametersType oParameters;
  CVector vdNewData, vdDetResiduals, vdClockSignal, vdClockSignalNoQSD, vdTemp;
  char strFileName[80];
  double dLogDetC;
  double dMinTOA, dMaxTOA, dT;
  int nClockTOAs=100, nClockSource=-1;
  CMatrix mdG, mdGTrans, mdBlockTemp, mdU, mdV,
	  mdGtot, mdGtotTrans, mdC,
	  mdGtCGt, mdM, mdMTrans, mdTemp,
	  mdSigma, mdTempLeft, mdSigmaInv,
	  mdMClock, mdGClock, mdGClockTrans;
  CVector vdS, vdResTotComp;

  for(int s=0; s<oConstants.nSources; ++s) {
    if(oConstants.poSources[s].eCorrelation == CID_Uniform) {
      nClockSource = s;
      break;
    } /* if eID */
  } /* for s */
  if(nClockSource == -1) {
    printf("No valid clock source found in parameters file\n");
    return;
  } /* if nClockSource */


  /* Set the parameters to the values in parameters.conf */
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

  /* Calculate the non-det residuals */
  vdDetResiduals = oData.vdData;
  ResidualsFromDetSources(oConstants, oParameters, oData, &vdDetResiduals);
  vdNewData = oData.vdData - vdDetResiduals;

  /* Construct the covariance matrix, with optimal clock parameters (taken from
   * "parameters.conf". Will do from MCMC chain later
   */
  SetCoherenceMatrix(oData, oParameters, oConstants);
  mdC = oData.mdC;

  /* Construct the G-matrix, or obtain it from the data reduction pipeline if it
   * has been set.
   */
  if(oConstants.bUseReducedBasis && oData.mdGReduce.Defined() && false) {
    printf("Obtaining G matrix from memory\n");
    /* Obtain mdGReduce from memory
     */
    mdG = oData.mdGReduce;
    mdGTrans = oData.mdGReduce[LO_TRANSPOSE];
  } else {
    printf("Calculating G matrix\n");
    /* Create mdGReduce
     *
     * Start with SVD of the design matrix
     */
    mdBlockTemp = oData.mdBlock;
    mdBlockTemp.SVD(mdU, mdV, vdS);

    /* Produce the timing model data reduction matrix
     */
    mdG.Initialize(oConstants.n, oConstants.n-oConstants.nMarParameters);
    for(int i=0; i<oConstants.n; ++i) {
      for(int j=oConstants.nMarParameters; j<oConstants.n; ++j) {
	mdG[i][j-oConstants.nMarParameters] = double(mdU[i][j]);
      } // for j
    } // for i
    mdGTrans = mdG[LO_TRANSPOSE];
  } /* if bUseReducedBasis */

  /* Figure out the first and the last absolute baseline TOAs across all
   * pulsars, and construct the observation dates of the clock corrections
   * (do about 100 of 'm)
   */
  dMinTOA = oConstants.poPulsars[0].pdTOA[0];
  dMaxTOA = dMinTOA;
  for(int a=0; a<oConstants.k; ++a) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; ++i) {
      if(dMinTOA > oConstants.poPulsars[a].pdTOA[i]) {
	dMinTOA = oConstants.poPulsars[a].pdTOA[i];
      } /* if dMinTOA */
      if(dMaxTOA < oConstants.poPulsars[a].pdTOA[i]) {
	dMaxTOA = oConstants.poPulsars[a].pdTOA[i];
      } /* if dMinTOA */
    } /* for i */
  } /* for a */
//  dT = dMaxTOA - dMinTOA;
//  dMinTOA -= 0.50 * dT;
//  dMaxTOA += 0.50 * dT;

  /* Add an extra pulsar to the array. This pulsar will have no timing-residual
   * sources, except the clock correlations with all other pulsars. These timing
   * residuals will then be estimated
   */
  oConstants.k++;
  oConstants.poPulsars[oConstants.k-1].nObservations = nClockTOAs;
  oConstants.poPulsars[oConstants.k-1].nTempo2Parameters = 0;
  oConstants.poPulsars[oConstants.k-1].bInterpolatorsSet = false;
  oConstants.n += nClockTOAs;
  oData.mdC.Initialize(oConstants.n, oConstants.n);
  oData.mdInvC.Initialize(oConstants.n, oConstants.n);
  for(int i=0; i<nClockTOAs; ++i) {
    oConstants.poPulsars[oConstants.k-1].pdTOA[i] = dMinTOA + i * (dMaxTOA - dMinTOA) / nClockTOAs;
  } /* for i */
  vdResTotComp.Initialize(mdG.m_pnDimSize[1]+nClockTOAs-3);
  vdTemp = mdGTrans * vdNewData;
  vdNewData = mdG * vdTemp;
  for(int i=0; i<mdG.m_pnDimSize[1]+nClockTOAs-3; ++i) {
    if(i < (mdG.m_pnDimSize[1])) {
      vdResTotComp[i] = double(vdTemp[i]);
    } else {
      vdResTotComp[i] = 0;
    } /* if i */
  } /* for i */

  /* Construct the clock QSD removal matrix
   */
  mdMClock.Initialize(nClockTOAs, 3);
  for(int i=0; i<nClockTOAs; ++i) {
    for(int j=0; j<3; ++j) {
      mdMClock[i][j] = pow(oConstants.poPulsars[oConstants.k-1].pdTOA[i], j);
    } /* for j */
  } /* for i */
  mdMClock.SVD(mdU, mdV, vdS);
  mdGClock.Initialize(nClockTOAs, nClockTOAs-3);
  for(int i=0; i<nClockTOAs; ++i) {
    for(int j=3; j<nClockTOAs; ++j) {
      mdGClock[i][j-3] = double(mdU[i][j]);
    } /* for j */
  } /* for i */
  mdGClockTrans = mdGClock[LO_TRANSPOSE];

  /* Set the clock source to also work on this new pulsar
   */
  oConstants.poSources[nClockSource].pbScope[oConstants.k-1] = true;

  /* Construct the G-matrix (Gt), including the new clock correction TOAs
   */
  mdGtot.Initialize(mdG.m_pnDimSize[0]+mdGClock.m_pnDimSize[0],
      mdG.m_pnDimSize[1]+mdGClock.m_pnDimSize[1]);
  for(int i=0; i<mdGtot.m_pnDimSize[0]; ++i) {
    for(int j=0; j<mdGtot.m_pnDimSize[1]; ++j) {
      if(i < mdG.m_pnDimSize[0] && j < mdG.m_pnDimSize[1]) {
	mdGtot[i][j] = double(mdG[i][j]);
      } else if(i < mdG.m_pnDimSize[0] || j < mdG.m_pnDimSize[1]) {
	mdGtot[i][j] = 0.0;
      } else {
	mdGtot[i][j] = double(mdGClock[i-mdG.m_pnDimSize[0]][j-mdG.m_pnDimSize[1]]);
      }/* if i && j */
    } /* for j */
  } /* for i */
  mdGtotTrans = mdGtot[LO_TRANSPOSE];

  /* Construct the larger covariance matrix for only the clock source
   */
  SetCoherenceMatrixClockSources(oData, oParameters, oConstants);

  /* Replace the data-data covariance with the one from the matrix we calculated
   * earlier. That contains all other sources.
   */
  for(int i=0; i<mdC.m_pnDimSize[0]; ++i) {
    for(int j=0; j<mdC.m_pnDimSize[1]; ++j) {
      oData.mdC[i][j] = double(mdC[i][j]);
    } /* for j */
  } /* for i */

  /* Set the auxiliary matrices required for the calculation of the clock
   * corrections
   */
  mdM.Initialize(mdGtot.m_pnDimSize[1], nClockTOAs);
  for(int i=0; i<mdGtot.m_pnDimSize[1]; ++i) {
    for(int j=0; j<nClockTOAs; ++j) {
      if(i < mdG.m_pnDimSize[1]) {
	mdM[i][j] = 0;
      } else {
	mdM[i][j] = double(mdGClock[j][i-mdG.m_pnDimSize[1]]);
      } /* if i == j */
    } /* for j */
  } /* for i */
  mdMTrans = mdM[LO_TRANSPOSE];
//  PrintMatrix(mdM);

//  PrintVector(vdResTotComp);

  /* Calculate the ML clock corrections, and the clock correction covariance
   * matrix
   */
  mdTemp = oData.mdC * mdGtot;
  mdGtCGt = mdGtotTrans * mdTemp;
//  PrintMatrix(mdGtCGt);
  mdGtCGt.InvertSVD();
//  mdGtCGt.Invert();
//  mdGtCGt.InvertChol();
  mdTempLeft = mdMTrans * mdGtCGt;
  mdSigma = mdTempLeft * mdM;
  mdSigmaInv = mdSigma;
//  PrintMatrix(mdSigma);
  mdSigmaInv.InvertSVD();
//  mdSigmaInv.InvertChol();
//  mdSigmaInv.Invert();
  mdTemp = mdSigmaInv * mdTempLeft;
  vdClockSignal = mdGClock * (mdGClockTrans * (mdTemp * vdResTotComp));
  vdClockSignalNoQSD = mdTemp * vdResTotComp;

  mdTemp = mdSigmaInv * mdGClock;
  mdSigmaInv = mdTemp * mdGClockTrans;
  mdTemp = mdGClockTrans * mdSigmaInv;
  mdSigmaInv = mdGClock * mdTemp;

  /* Write the files to disk */
  FILE *pFile;
  try {
    sprintf(strFileName, "./clockcorr.txt");
    if(! (pFile = fopen(strFileName, "w+")) ) throw 1;
    fprintf(pFile, "# TOA  clockcorr  deltaclockcorr\n");
    for(int i=0; i<nClockTOAs; ++i) {
      fprintf(pFile, "%e  %e  %e  %e\n",
	  2008 + ((54500 + oConstants.poPulsars[oConstants.k-1].pdTOA[i] / (3600*24))-54466)/365.25,
	  double(vdClockSignal[i]),
	  sqrt(double(mdSigmaInv[i][i])),
	  double(mdSigmaInv[i][i]));
    } /* for i */
    fclose(pFile);
  } catch(int nError) {
    fprintf(stderr, "Error number: %i\n", nError);
  } /* try */

  /* Remove the extra pulsar from the array again */
  oConstants.poSources[nClockSource].pbScope[oConstants.k-1] = false;
  oConstants.k--;
  oConstants.n -= nClockTOAs;
  oData.mdC.Initialize(oConstants.n, oConstants.n);
  oData.mdInvC.Initialize(oConstants.n, oConstants.n);
#endif
} // CreateClockCorrectedResiduals

/* Writes all the MCMC datapoints to the file strFileName.
 *
 * Last update: 30-01-2008
 * */
bool WriteMCMCData(int nMCMCSteps, int nParameters, MCMCDataPoint *pDat, const char *strFileName) {
  FILE *pFile;

  char strMsg[160], strFileVersion[16], strVersion[16];
  strcpy(strMsg, "Writing data to '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  strcpy(strVersion, VERSION);

  try {
    if(! (pFile = fopen(strFileName, "wb+")) ) throw 1;

    if(! fwrite(strVersion , sizeof(char) , 16 , pFile) ) throw 2;
    if(! fwrite(&nMCMCSteps , sizeof(int) , 1 , pFile) ) throw 2;

    if(! fwrite(&nParameters , sizeof(int) , 1 , pFile) ) throw 2;
    for(int i=0; i<nMCMCSteps; i++) {
      if(! fwrite(pDat[i].oParameters.pdPar, sizeof(double), nParameters, pFile)) throw 4;
      if(! fwrite(&(pDat[i].dLogLik) , sizeof(double) , 1 , pFile) ) throw 8;
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return false;
  } // try

  PrintSuccess();
  return true;
} // WriteMCMCData


/* Writes all the MCMC datapoints to the file strFileName. This function writes
 * the data of all the MCMC chains in one file
 *
 * Last update: 2-11-2008
 * */
bool WriteMCMCData(int nProcesses, int nMCMCSteps, int nParameters, MCMCDataPoint **ppDat, char *strFileName) {
  FILE *pFile;

  char strMsg[160], strFileVersion[16], strVersion[16];
  int nSteps=nMCMCSteps*nProcesses;

  strcpy(strVersion, VERSION);

  try {
    if(! (pFile = fopen(strFileName, "wb+")) ) throw 1;

    if(! fwrite(strVersion , sizeof(char) , 16 , pFile) ) throw 2;
    if(! fwrite(&nSteps , sizeof(int) , 1 , pFile) ) throw 2;

    if(! fwrite(&nParameters , sizeof(int) , 1 , pFile) ) throw 2;
    for(int i=0; i<nMCMCSteps; i++) {
      for(int j=0; j<nProcesses; j++) {
	if(! fwrite(ppDat[j][i].oParameters.pdPar, sizeof(double), nParameters, pFile)) throw 4;
	if(! fwrite(&(ppDat[j][i].dLogLik) , sizeof(double) , 1 , pFile) ) throw 8;
      } // for j
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return false;
  } // try

  return true;
} // WriteMCMCData


/* Reads all the MCMC datapoints from the datafile strFileName.
 *
 * Note: deprecated
 * */
bool ReadMCMCData(int &nMCMCSteps, int nParameters, MCMCDataPoint *&pDat, const char *strFileName) {
  FILE *pFile;
  int nTemp;

  char strMsg[160], strFileVersion[16];
  strcpy(strMsg, "Reading data from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  try {
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 2;
    if(! fread(&nMCMCSteps , sizeof(int) , 1 , pFile) ) throw 2;
    if(! fread(&nTemp, sizeof(int), 1, pFile) ) throw 3;
    if(nTemp != nParameters) printf("WARNING: nTemp != nParameters\n");

    pDat = new MCMCDataPoint[nMCMCSteps];
    for(int i=0; i<nMCMCSteps; i++) {
      if(! fread(pDat[i].oParameters.pdPar, sizeof(double), nTemp ,pFile) ) throw 4;
      if(! fread(&(pDat[i].dLogLik) , sizeof(double) , 1 , pFile) ) throw 8;
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
//    printf("\n\nError numbers: %i\n\n", nError);
    return false;
  } // try

  PrintSuccess();
  return true;
} // ReadMCMCData



/* This function starts a file for MCMC parameter data. The file is created, or
 * truncated if it exists, and info about the ban version and parameters is
 * written.
 *
 * Last update: 4-11-2008
 * */
bool WriteMCMCDataFileStart(int nParameters, int nTMParameters, const char *strFileName) {
  FILE *pFile;

  char strVersion[16];

  strcpy(strVersion, VERSION);

  try {
    if(! (pFile = fopen(strFileName, "wb+")) ) throw 1;

    if(! fwrite(strVersion , sizeof(char) , 16 , pFile) ) throw 2;

    if(! fwrite(&nParameters , sizeof(int) , 1 , pFile) ) throw 3;

    if(! fwrite(&nTMParameters , sizeof(int) , 1 , pFile) ) throw 4;

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return false;
  } // try

  return true;
} // WriteMCMCDataFileStart



/* This function appends data to a MCMC parameter file.
 *
 * Last update: 4-11-2008
 * */
bool WriteMCMCDataFileAppend(int nMCMCSteps, int nParameters, int nTMParameters, MCMCDataPoint *pDat, const char *strFileName, bool bCalcTMPars) {
  FILE *pFile;
  double *pdBuf=NULL;

  if(! bCalcTMPars) {
    pdBuf = new double[nTMParameters*nTMParameters];
  } // if bCalcTMPars

  try {
    if(! (pFile = fopen(strFileName, "a")) ) throw 1;

    for(int i=0; i<nMCMCSteps; i++) {
      if(! fwrite(pDat[i].oParameters.pdPar, sizeof(double), nParameters, pFile)) throw 2;
      if(! fwrite(&(pDat[i].dLogLik) , sizeof(double) , 1 , pFile) ) throw 3;
      if(bCalcTMPars) {
	if(! fwrite(pDat[i].vdChi.m_pdData, sizeof(double), nTMParameters, pFile)) throw 4;
	if(! fwrite(pDat[i].mdCXiInv.m_pdData, sizeof(double), nTMParameters*nTMParameters, pFile)) throw 5;
      } else {
//	if(! fwrite(pdBuf, sizeof(double), nTMParameters, pFile)) throw 4;
//	if(! fwrite(pdBuf, sizeof(double), nTMParameters*nTMParameters, pFile)) throw 5;
      } // if bCalcTMPars
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return false;
  } // try

  if(pdBuf) {
    delete[] pdBuf;
  } // if bCalcTMPars

  return true;
} // WriteMCMCDataFileAppend


/* This function reads the MCMC parameter data from a file
 *
 * Last update: 4-11-2008
 * */
bool ReadMCMCDataFile(int &nMCMCSteps, int nParameters, int nTMParameters, CMCMCDataPoint *&pDat, const char *strFileName, bool bCalcTMPars) {
  FILE *pFile;
  int nFileParameters, nFileTMParameters;
  unsigned long nFileLength;
  double *pdBuf=NULL;

  char strMsg[160], strFileVersion[16];
  strcpy(strMsg, "Reading data from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  try {
    // Open the file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(nFileParameters != nParameters) {
      printf("WARNING: nFileParameters != nParameters\n");
      printf("nFileParameters: %i,  nParameters: %i\n", nFileParameters, nParameters);
    } // if nFileParameters
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;
    if(nFileTMParameters != nTMParameters) {
      printf("WARNING: nFileTMParameters != TMnParameters\n");
      printf("nFileTMParameters: %i,  nTMParameters: %i\n", nFileTMParameters, nTMParameters);
    } // if nFileTMParameters
    if(! bCalcTMPars) {
      pdBuf = new double[nFileTMParameters * nFileTMParameters];
    } // if bCalcTMPars

    // Calculate nMCMCSteps
    if(bCalcTMPars) {
      nMCMCSteps = (nFileLength-sizeof(char)*16 - 2*sizeof(int)) /
	(
	  sizeof(double)*(nFileParameters+1+nFileTMParameters*(nFileTMParameters+1))
	);
    } else {
      nMCMCSteps = (nFileLength-sizeof(char)*16 - 2*sizeof(int)) /
	(
	  sizeof(double)*(nFileParameters+1)
	);
    } // if bCalcTMPars
    if(nMCMCSteps <= 0) throw 7;

//    printf("\n\nnFileLength: %i,  nMCMCSteps: %i\n", nFileLength, nMCMCSteps);
//    printf("n = sizeof(MCMCDataPoint) = %i\n", sizeof(MCMCDataPoint));
//    printf("n*nMCMCSteps = %i\n", sizeof(MCMCDataPoint)*nMCMCSteps);
//    printf("Trying: 'pDat = new MCMCDataPoint[nMCMCSteps];'\n");
    pDat = new CMCMCDataPoint[nMCMCSteps];
    for(int i=0; i<nMCMCSteps; i++) {
      pDat[i].Initialize(nFileParameters, nFileTMParameters, bCalcTMPars);
//    printf("Succes!\n");
    } // for i

    for(int i=0; i<nMCMCSteps; i++) {
      if(! fread(pDat[i].pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(&(pDat[i].dLogLik), sizeof(double), 1, pFile) ) throw 9;
      if(bCalcTMPars) {
	if(! fread(pDat[i].pdChi, sizeof(double), nFileTMParameters, pFile) ) throw 10;
	if(! fread(pDat[i].pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 11;
      } else {
//	if(! fread(pdBuf, sizeof(double), nFileTMParameters, pFile) ) throw 10;
//	if(! fread(pdBuf, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 11;
      } // if bCalcTMPars
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    fclose(pFile);
    PrintFailed();
    printf("\n\nError number: %i\n\n", nError);
    return false;
  } // try

  if(pdBuf) {
    delete[] pdBuf;
  } // if bCalcTMPars

  PrintSuccess();
  return true;
} // ReadMCMCDataFile


/* This function starts a file for MCMC ensemble parameter data. The file is
 * created, or truncated if it exists, and info about the ban version and
 * parameters is written.
 *
 * Last update: 4-08-2011
 * */
bool WriteMCMCEnsembleDataFileStart(int nDataSets, int nObservations, int nParameters, int nTMParameters, double *pdData, const char *strFileName) {
  FILE *pFile;

  char strVersion[16];

  strcpy(strVersion, VERSION);

  try {
    if((pFile = fopen(strFileName, "r")) ) {
      fclose(pFile);
      throw 0;
    } // if pFile
    if(! (pFile = fopen(strFileName, "wb")) ) throw 1;

    if(! fwrite(strVersion , sizeof(char) , 16 , pFile) ) throw 2;

    if(! fwrite(&nDataSets , sizeof(int) , 1 , pFile) ) throw 3;

    if(! fwrite(&nObservations , sizeof(int) , 1 , pFile) ) throw 4;

    if(! fwrite(&nParameters , sizeof(int) , 1 , pFile) ) throw 5;

    if(! fwrite(&nTMParameters , sizeof(int) , 1 , pFile) ) throw 6;

    if(! fwrite(pdData , sizeof(double) , nDataSets*nObservations , pFile) ) throw 7;

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
//    fprintf(stderr, "\nError number: %i\n\n", nError);
    return false;
  } // try

  return true;
} // WriteMCMCEnsembleDataFileStart


/* This function appends data to a MCMC ensemble parameter file.
 *
 * Last update: 4-08-2011
 * */
bool WriteMCMCEnsembleDataFileAppend(int nMCMCSteps, int nDataSets, int nParameters, int nTMParameters, MCMCDataEnsemble *pDat, const char *strFileName) {
  FILE *pFile;

  try {
    if(! (pFile = fopen(strFileName, "a")) ) throw 1;

    for(int i=0; i<nMCMCSteps; i++) {
      if(! fwrite(pDat[i].oParameters.pdPar, sizeof(double), nParameters, pFile)) throw 2;
      if(! fwrite(pDat[i].vdLogLik.m_pdData , sizeof(double) , nDataSets , pFile) ) throw 3;
      if(! fwrite(&(pDat[i].dLogLik) , sizeof(double) , 1 , pFile) ) throw 3;
      if(! fwrite(pDat[i].mdChi.m_pdData, sizeof(double), nTMParameters*nDataSets, pFile)) throw 4;
      if(! fwrite(pDat[i].mdCXiInv.m_pdData, sizeof(double), nTMParameters*nTMParameters, pFile)) throw 5;
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return false;
  } // try

  return true;
} // WriteMCMCEnsembleDataFileAppend


/* This function reads the MCMC ensemble parameter data from a file
 *
 * Last update: 4-08-2011
 * */
bool ReadMCMCEnsembleDataFile(int &nMCMCSteps, int &nDataSets, int nParameters, int nTMParameters, CMCMCDataEnsemble *&pDat, const char *strFileName, CMatrix &mdDataSets) {
  FILE *pFile;
  int nFileParameters, nFileTMParameters, nFileDataSets, nFileObservations;
  int nDataSetsLength;
  unsigned long nFileLength;

  char strMsg[160], strFileVersion[16];
  strcpy(strMsg, "Reading data from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  try {
    // Open the file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileDataSets, sizeof(int), 1, pFile) ) throw 5;
    nDataSets = nFileDataSets;
//    if(nFileParameters != nParameters) {
//      printf("WARNING: nFileDataSets != nDataSets\n");
//      printf("nFileDataSets: %i,  nDataSets: %i\n", nFileDataSets, nDataSets);
//    } // if nFileDataSets
    if(! fread(&nFileObservations, sizeof(int), 1, pFile) ) throw 5;
//    if(nFileObservations != nObservations) {
//      printf("WARNING: nFileObservations != nObservations\n");
//      printf("nFileObservations: %i,  nObservations: %i\n", nFileObservations, nObservations);
//    } // if nFileObservations
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(nFileParameters != nParameters) {
      printf("WARNING: nFileParameters != nParameters\n");
      printf("nFileParameters: %i,  nParameters: %i\n", nFileParameters, nParameters);
    } // if nFileParameters
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;
    if(nFileTMParameters != nTMParameters) {
      printf("WARNING: nFileTMParameters != TMnParameters\n");
      printf("nFileTMParameters: %i,  nTMParameters: %i\n", nFileTMParameters, nTMParameters);
    } // if nFileTMParameters

    nDataSetsLength = nFileObservations*nFileDataSets;
    mdDataSets.Initialize(nFileObservations, nFileDataSets);
    if(! fread(mdDataSets.m_pdData, sizeof(double), nDataSetsLength, pFile) ) throw 6;

    // Calculate nMCMCSteps
    nMCMCSteps = (nFileLength-sizeof(char)*16 - 4*sizeof(int) - nDataSetsLength*sizeof(double)) / (
        sizeof(double)*(
	  nFileParameters + nFileDataSets + 1 +
	  nFileTMParameters*nFileDataSets +
	  nFileTMParameters*nFileTMParameters));

    if(nMCMCSteps <= 0) throw 7;

    pDat = new CMCMCDataEnsemble[nMCMCSteps];
    for(int i=0; i<nMCMCSteps; i++) {
      pDat[i].Initialize(nFileDataSets, nFileParameters, nFileTMParameters);
    } // for i

    for(int i=0; i<nMCMCSteps; i++) {
      if(! fread(pDat[i].pdPar, sizeof(double), nFileParameters, pFile) ) throw 8;
      if(! fread(pDat[i].pdLogLik, sizeof(double), nFileDataSets, pFile) ) throw 9;
      if(! fread(&(pDat[i].dLogLik), sizeof(double), 1, pFile) ) throw 10;
      if(! fread(pDat[i].pdChi, sizeof(double), nFileDataSets*nFileTMParameters, pFile) ) throw 11;
      if(! fread(pDat[i].pdCXiInv, sizeof(double), nFileTMParameters*nFileTMParameters, pFile) ) throw 12;
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    fclose(pFile);
    PrintFailed();
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
    return false;
  } // try

  PrintSuccess();
  return true;
} // ReadMCMCEnsembleDataFile


/* This function reads the datasets from an MCMC ensemble parameter data
 *
 * Last update: 8-08-2011
 * */
bool ReadMCMCEnsembleDataFileSets(int &nDataSets, int &nMCMCSteps, const char *strFileName, CMatrix &mdDataSets) {
  FILE *pFile;
  int nFileParameters, nFileTMParameters, nFileDataSets, nFileObservations;
  int nDataSetsLength;
  unsigned long nFileLength;

  char strMsg[160], strFileVersion[16];
  strcpy(strMsg, "Reading data from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  try {
    // Open the file
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    // Get the length of the file (in bytes)
    if(fseek(pFile, 0, SEEK_END) ) throw 2;
    nFileLength = ftell(pFile);
    if(fseek(pFile, 0, SEEK_SET) ) throw 3;

    // Read the header
    if(! fread(strFileVersion, sizeof(char), 16, pFile) ) throw 4;
    if(! fread(&nFileDataSets, sizeof(int), 1, pFile) ) throw 5;
    nDataSets = nFileDataSets;

    if(! fread(&nFileObservations, sizeof(int), 1, pFile) ) throw 5;
    if(! fread(&nFileParameters, sizeof(int), 1, pFile) ) throw 5;
    if(! fread(&nFileTMParameters, sizeof(int), 1, pFile) ) throw 6;

    nDataSetsLength = nFileObservations*nFileDataSets;
    mdDataSets.Initialize(nFileObservations, nFileDataSets);
    if(! fread(mdDataSets.m_pdData, sizeof(double), nDataSetsLength, pFile) ) throw 6;

    // Calculate nMCMCSteps
    nMCMCSteps = (nFileLength-sizeof(char)*16 - 4*sizeof(int) - nDataSetsLength*sizeof(double)) / (
        sizeof(double)*(
	  nFileParameters + nFileDataSets + 1 +
	  nFileTMParameters*nFileDataSets +
	  nFileTMParameters*nFileTMParameters));

    if(nMCMCSteps <= 0) throw 7;

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    fclose(pFile);
    PrintFailed();
    fprintf(stderr, "\n\nError number: %i\n\n", nError);
    return false;
  } // try

  PrintSuccess();
  return true;
} // ReadMCMCEnsembleDataFileSets




/* Calculate the rms of the residuals of a pulsar
 * */
double CalcResidualsRms(SDataType &oData, SConstantsType &oConstants, int nCalcRmsPulsar) {
  CMatrix mdTemp, mdTemp2;
  CVector vdReducedResiduals, vdQSD, vdR;

  // Now vdR is the actual correlation, vdZeta is the theoretical correlation
  // vdTheta is the angle between the pairs of pulsars

  double dRBar=0, dRSigma=0;
  int nIndex=0;

  try {
    // Calculate the best-fit quadratic spindown parameters (or gen parameters)
    mdTemp2 = oData.mdGenBlock[LO_TRANSPOSE];
    mdTemp = mdTemp2 * oData.mdGenBlock;
    mdTemp.InvertChol();
    vdQSD = mdTemp * (oData.mdGenBlock[LO_TRANSPOSE] * oData.vdData);

    // Subtract linear and quadratic terms from the residuals
    vdReducedResiduals = oData.vdData - oData.mdGenBlock * vdQSD;

    vdR.Initialize(oConstants.poPulsars[nCalcRmsPulsar].nObservations);

    for(int a=0; a<nCalcRmsPulsar; a++)
      nIndex += oConstants.poPulsars[a].nObservations;

    for(int i=0; i<oConstants.poPulsars[nCalcRmsPulsar].nObservations; i++) {
      vdR[i] = double(vdReducedResiduals[nIndex]);
      nIndex++;
    } // for i

    for(int i=0; i<oConstants.poPulsars[nCalcRmsPulsar].nObservations; i++)
      dRBar += double(vdR[i]);
    dRBar /= oConstants.poPulsars[nCalcRmsPulsar].nObservations;

    for(int i=0; i<oConstants.poPulsars[nCalcRmsPulsar].nObservations; i++)
      dRSigma += (double(vdR[i]) - dRBar)*(double(vdR[i]) - dRBar);

    dRSigma = sqrt(dRSigma / oConstants.poPulsars[nCalcRmsPulsar].nObservations);
  } catch(ELinearError err) {
    dRSigma = 0;
  } // try

  return dRSigma;
} // CalcResidualsRms

/* Calculate the chi-squared value of the residuals pulsar nCalcChisqPulsar
 * */
double CalcResidualsChisq(SDataType &oData, SConstantsType &oConstants, int nCalcChisqPulsar) {
  CMatrix mdTemp, mdTemp2, mdInvC;
  CVector vdReducedResiduals, vdQSD;
  int nIndex;
  double dChisq;

  try {
    // Calculate the inverse of mdC
    mdInvC.Initialize(oData.vdData.m_pnDimSize[0], oData.vdData.m_pnDimSize[0]);
    for(int i=0; i<oData.vdData.m_pnDimSize[0]; i++) {
      for(int j=0; j<oData.vdData.m_pnDimSize[0]; j++) {
	mdInvC[i][j] = 0;
	if(i==j)
	  mdInvC[i][j] = 1 / gsl_pow_2(double(oData.vdDataErr[j]));
      } // for j
    } // for i

    // Calculate the best-fit quadratic spindown parameters (or gen parameters)
    mdTemp2 = oData.mdBlock[LO_TRANSPOSE];
    mdTemp = mdTemp2 * mdInvC * oData.mdBlock;
    mdTemp.InvertChol();
    vdQSD = mdTemp * (oData.mdBlock[LO_TRANSPOSE] * (mdInvC * oData.vdData));

    // Subtract linear and quadratic terms from the residuals
    vdReducedResiduals = oData.vdData - oData.mdBlock * vdQSD;

    // Jump to the observations of this particular pulsar
    nIndex = 0;
    dChisq = 0;
    for(int a=0; a<nCalcChisqPulsar; a++)
      nIndex += oConstants.poPulsars[a].nObservations;

    // Calculate the Chi-squared value
    for(int i=0; i<oConstants.poPulsars[nCalcChisqPulsar].nObservations; i++) {
      dChisq += gsl_pow_2(double(vdReducedResiduals[nIndex]) / double(oData.vdDataErr[nIndex]));
      //printf("Residual / deltaResidual [0] = %e / %e\n", double(vdReducedResiduals[nIndex]), double(oData.vdDataErr[nIndex]));
      nIndex++;
    } // for i
  } catch(ELinearError err) {
    dChisq = 0.0001;
  } // try

  return dChisq;
} // CalcResidualsChisq


/* This function calculates the points of the likelihood function of one single
 * parameter. The values are returned in the vectors vdX, vdY
 * */
void Calculate1DPlot(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameterNumber, CVector &vdX, CVector &vdY) {
    CVector vdPlot, vdFunctionValues, vdTemp, vdThetaMin, vdThetaMax, vdTheta;
    int nPlotPoints;
    int k, l, nParameters, nCount;
    double dEvidence, dMin;

    nParameters = oConstants.nParameters;
    vdThetaMin.Initialize(nParameters);
    vdThetaMax.Initialize(nParameters);
    vdTheta.Initialize(nParameters);

    // Initialize the Width-of-the-Parameters-Vector
    StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
    for(int i=0; i<nParameters; i++) {
      vdThetaMin[i] = oParameters.pdParMinBound[i];
      vdThetaMax[i] = oParameters.pdParMaxBound[i];
    } // for i

    nPlotPoints = oConstants.nSmallPlotPoints;	// Number of points used in plot

    k = oConstants.k;
    l = oConstants.l;

    vdPlot.Initialize(nPlotPoints + 1);
    for(int i=0; i<=nPlotPoints; i++) {
      vdPlot[i] = double(vdThetaMin[nParameterNumber]) + i*((double(vdThetaMax[nParameterNumber]) - double(vdThetaMin[nParameterNumber]))/nPlotPoints);
    } // for i
    vdFunctionValues = vdPlot;

    InitProgressBar("Generating 1D plot...");
    DrawProgressBar(0);

// From here, initialize the values of the parameters:
    ParametersToVector(oParameters, vdTheta);

//    fprintf(stderr, "%3i/%3i:               \r", 0, nPlotPoints);
    for(int i=0; i<=nPlotPoints; i++) {
      vdTheta[nParameterNumber] = double(vdPlot[i]);
      VectorToParameters(oParameters, vdTheta);

      vdFunctionValues[i] = LogLikelihoodTimesPrior(oData, oParameters, oConstants);
//      vdFunctionValues[i] = LogLikelihoodSets(oData, oParameters, oConstants);

//      fprintf(stderr, "\n%3i/%3i: (%e, %e)", i, nPlotPoints, double(vdPlot[i]), double(vdFunctionValues[i]));
      DrawProgressBar(int(100.0 * i / pow(nPlotPoints+1, 1)));
    } // for i

    FinishProgressBar();

    WritePlot("logplotdata-1d.txt", vdPlot, vdFunctionValues);

    // Rescale the function for plotting
    dMin = Min(vdFunctionValues);
    vdTemp = Exp_1(vdFunctionValues - dMin);

    // Calculate the integral
    dEvidence = 0;
    for(int i=1; i<nPlotPoints; i++) {
      dEvidence += 0.5 * exp(dMin - double(vdFunctionValues[i])) *
	(double(vdPlot[i+1]) - double(vdPlot[i-1]));
    } // for i

    vdX = vdPlot;
    vdY = vdTemp;

    printf("dEvidence = exp(%e)\n", log(dEvidence) - dMin);
} // Calculate1DPlot


/* Makes a plot of the likelihood function. Usually of the GWB amplitude.
 *
 * */
void MakePlot(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameterNumber) {
  char strFileName[160];
  CVector vdX, vdY;
  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcat(strFileName, "plotdata.txt");

  Calculate1DPlot(oData, oParameters, oConstants, nParameterNumber, vdX, vdY);

  // Writa data to file so it can be used in with gnuplot
  WritePlot(strFileName, vdX, vdY);
//    WritePlot(strFileName, vdPlot, vdFunctionValues);
//    WritePlot("../data/workdata/logplot.txt", vdPlot, vdFunctionValues);
} // MakePlot



/* This function calculates the values of the 3D plot.
 *
 * TODO: Sourcenumber and PlotParameterNumber aren't used anymore. Remove them
 * from the code
 * */
void Calculate3DPlot(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ) {
  CVector vdPlotX, vdPlotY, vdPlotZ, vdTemp, vdTheta, vdThetaMin, vdThetaMax;
  int nPlotPoints;
  int k, l, nIndex1, nIndex2;
  int nPulsarNumber = -1, nSourceNumber, nPlotParameterNumber;
  double dTemp;

  nPlotPoints = int(sqrt(double(oConstants.nSmallPlotPoints)));	// Number of points used in plot

  k = oConstants.k;
  l = oConstants.l;

  vdPlotX.Initialize( (nPlotPoints + 1)*(nPlotPoints + 1) );
  vdPlotY.Initialize( (nPlotPoints + 1)*(nPlotPoints + 1) );
  vdPlotZ.Initialize( (nPlotPoints + 1)*(nPlotPoints + 1) );
  
  vdX.Initialize(nPlotPoints+1);
  vdY.Initialize(nPlotPoints+1);
  mdZ.Initialize(nPlotPoints+1, nPlotPoints+1);

  vdTheta.Initialize(oConstants.nParameters);
  vdThetaMin = vdTheta;
  vdThetaMax = vdTheta;
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  for(int i=0; i<oConstants.nParameters; i++) {
    vdThetaMin[i] = oParameters.pdParMinBound[i];
    vdThetaMax[i] = oParameters.pdParMaxBound[i];
  } // for i

#if 0
  if(nPulsarNumber == -1) {
    for(int s=0; s<oConstants.nSources; s++) {
      if(oConstants.poSources[s].eCorrelation == CID_GR) {
	// This is a gravitational wave background
	nSourceNumber = s;
	break;
      } // if eCorrelation
    } // for s
  } else {
    nSourceNumber = nPulsarNumber;
  } // if nPulsarNumber
  nPlotParameterNumber = oConstants.poSources[nSourceNumber].nFirstParIndex;
#endif

  InitProgressBar("Generating 3D plot...");
  DrawProgressBar(0);
  for(int i=0; i<=nPlotPoints; i++) {
      for(int j=0; j<=nPlotPoints; j++) {
	vdPlotX[j+(nPlotPoints+1)*i] = oParameters.pdParMinBound[nParameter1] +
	  i*((oParameters.pdParMaxBound[nParameter1] -
		oParameters.pdParMinBound[nParameter1])/nPlotPoints);
	vdPlotY[j+(nPlotPoints+1)*i] = oParameters.pdParMinBound[nParameter2] +
	  j*((oParameters.pdParMaxBound[nParameter2] -
		oParameters.pdParMinBound[nParameter2])/nPlotPoints);
      } // for j
  } // for i


  for(int i=0; i<(nPlotPoints + 1)*(nPlotPoints + 1); i++) {
    oParameters.pdPar[nParameter1] = double(vdPlotX[i]);
    oParameters.pdPar[nParameter2] = double(vdPlotY[i]);

    try {
//      nIndex1 = (i - nIndex1)/(nPlotPoints+1);
      nIndex1 = (i )/(nPlotPoints+1);
      nIndex2 = i % (nPlotPoints+1);
      vdPlotZ[i] = LogLikelihoodTimesPrior(oData, oParameters, oConstants);
//      fprintf(stderr, "\n %e\n", double(vdPlotZ[i]));
    } catch (ELinearError err) {  // Error handling
      switch(err) {
      case ELENotDefined:
	printf("Not Defined\n");
	break;
      case ELEWrongClassType:
	printf("Wrong Class Type\n");
	break;
      case ELEBadIndex:
	printf("Wrong Index number\n");
	break;
      case ELEDimensionMisMatch:
	printf("Matrix Dimensions Do Not Match\n");
	break;
      case ELELapack:
	printf("Lapack error!\n");
	break;
      default:
	printf("Default error!\n");
	break;
      } // switch
      printf("\nValues of the parameters of calculation:\n");
      printf("parameter 1: %e\n", oParameters.pdPar[nParameter1]);
      printf("parameter 2: %e\n", oParameters.pdPar[nParameter2]);

      vdPlotZ[i] = vdPlotZ[i-1];
    } // try

    DrawProgressBar(int(100.0*i/pow(nPlotPoints+1,2)));
  } // for i
  FinishProgressBar();

  // Rescale the function for plotting
//  dTemp = Min(vdPlotZ);
//  vdTemp = Exp_1(vdPlotZ - dTemp);


  for(int i=0; i<(nPlotPoints + 1)*(nPlotPoints + 1); i++) {
//    nIndex1 = (i - nIndex1)/(nPlotPoints+1);
    nIndex1 = (i )/(nPlotPoints+1);
    nIndex2 = i % (nPlotPoints+1);

    vdX[nIndex1] = double(vdPlotX[i]);
    vdY[nIndex2] = double(vdPlotY[i]);
//    mdZ[nIndex1][nIndex2] = double(vdTemp[i]);
    mdZ[nIndex1][nIndex2] = double(vdPlotZ[i]);
//    mdZ[nIndex1][nIndex2] = dTemp - double(vdPlotZ[i]);
  } // for i
} // Calculate3DPlot


// Make a 3D plot of the likelihood function (w.r.t. a pulsar-noise parameter)
void Make3DPlot(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameter1, int nParameter2) {
  CVector vdX, vdY, vdZ;
  CMatrix mdZ;
  CVector vdTemp, vdTemp2;
  CMatrix mdTemp;
  int nPlotPoints, i, j;
  char strFileName[160];

  Calculate3DPlot(oData, oParameters, oConstants, nParameter1, nParameter2, vdX, vdY, mdZ);

  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcat(strFileName, "plot3d.txt");

  nPlotPoints = vdX.m_pnDimSize[0]-1;
  vdTemp = vdX;
  vdTemp2 = vdY;
  mdTemp = mdZ;

  vdX.Initialize( (nPlotPoints + 1)*(nPlotPoints + 1) );
  vdY.Initialize( (nPlotPoints + 1)*(nPlotPoints + 1) );
  vdZ.Initialize( (nPlotPoints + 1)*(nPlotPoints + 1) );

  for(i=0; i<=nPlotPoints; i++) {
    for(j=0; j<=nPlotPoints; j++) {
      vdX[i+(nPlotPoints+1)*j] = double(vdTemp[i]);
      vdY[i+(nPlotPoints+1)*j] = double(vdTemp2[j]);
      vdZ[i+(nPlotPoints+1)*j] = double(mdTemp[i][j]);
    } // for j
  } // for i

  // Writa data to file so it can be used in with gnuplot
  WritePlot(strFileName, vdX, vdY, vdZ);
} // Make3DPlot


//#define DO_PROFILING
/* The likelihood function for reduced bases
 */
double LogLikelihoodReducedBasis(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  int nParIndex, nAddIndex;
  CMatrix mdTemp, mdTemp2, mdGT, mdRedTempTot;
  double dLogDetGCG, dReturnValue;
  SParametersType oParametersTemp;
  int nP;

#ifdef DO_PROFILING
  static int nEvals = 0;

  static time_t tA1=0, tA2=0, tA3=0, tA4=0, tB=0;
  time_t tStartA1, tEndA1,
	 tStartA2, tEndA2,
	 tStartA3, tEndA3,
	 tStartA4, tEndA4,
	 tStartB, tEndB;
#endif

  // It is unclear at this point to what extend this block is actually necessary
  if(oData.bRedLikFirstRun) {
    /* It is the first time we run this function with this configuration, so we
     * must initialise a couple of values
     */

    /* Allocate memory for all GCG matrices */
    oData.mdGCG.Initialize(oData.mdGReduce.m_pnDimSize[1], oData.mdGReduce.m_pnDimSize[1]);
    for(int a=0; a<oConstants.k; ++a) {
      oConstants.poPulsars[a].mdGCG.Initialize(
	  oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1],
	  oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1]);
    } // for a

    /* Now loop over all sources, and allocate GCG combinations */
    for(int s=0; s<oConstants.nSources; ++s) {
      // Is this even required?
      if(oConstants.poSources[s].oIntAcc.bSet) {
          if(oConstants.poSources[s].oIntAcc.nPulsar < 0) {
              // It's for _all_ pulsars
              oConstants.poSources[s].oIntAcc.mdGCiG.Initialize(
                      oData.mdGReduce.m_pnDimSize[1],
                      oData.mdGReduce.m_pnDimSize[1]);
          } else {
              oConstants.poSources[s].oIntAcc.mdGCiG.Initialize(
                      oConstants.poPulsars[oConstants.poSources[s].oIntAcc.nPulsar].mdGReduce.m_pnDimSize[1],
                      oConstants.poPulsars[oConstants.poSources[s].oIntAcc.nPulsar].mdGReduce.m_pnDimSize[1]);
          } // if nPulsar
      } // if oIntAcc
    } // for s
    oData.bRedLikFirstRun = false;
  } // if bRedLikFirstRun


#ifdef DO_PROFILING
    tStartA1 = clock();
    tEndA1 = tStartA1;
#endif

  /* Set all the GCG elements equal to zero
   */
  for(int a=0; a<oConstants.k; ++a) {
    for(int i=0; i<oConstants.poPulsars[a].mdGCG.m_pnDimSize[0]; ++i) {
      for(int j=0; j<oConstants.poPulsars[a].mdGCG.m_pnDimSize[1]; ++j) {
	oConstants.poPulsars[a].mdGCG[i][j] = 0;
      } // for j
    } // for i
  } // for a
  for(int i=0; i<oData.mdGCG.m_pnDimSize[0]; ++i) {
    for(int j=0; j<oData.mdGCG.m_pnDimSize[1]; ++j) {
      oData.mdGCG[i][j] = 0;
    } // for j
  } // for i

#ifdef DO_PROFILING
      tEndA1 = clock();
      tStartA2 = tEndA1;
#endif

  /* Now loop over all sources, and create GCG combinations */
  for(int s=0; s<oConstants.nSources; ++s) {
    // If this source is an efac/equad
    if(oConstants.poSources[s].oAmpAcc.bSet) {
      // The amplitude acceleration thingy is set
      nParIndex=0;
      if(oConstants.poSources[s].oAmpAcc.nPulsar < 0) {
        for(int i=0; i<oData.mdGCG.m_pnDimSize[0]; ++i) {
          for(int j=0; j<oData.mdGCG.m_pnDimSize[1]; ++j) {
            oData.mdGCG[i][j] += double(oConstants.poSources[s].oAmpAcc.mdGCaG[i][j]) *
              exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+nParIndex]) *
              exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+nParIndex]);
          } // for j
        } // for i
      } else {
        for(int i=0; i<oConstants.poPulsars[oConstants.poSources[s].oAmpAcc.nPulsar].mdGCG.m_pnDimSize[0]; ++i) {
          for(int j=0; j<oConstants.poPulsars[oConstants.poSources[s].oAmpAcc.nPulsar].mdGCG.m_pnDimSize[0]; ++j) {
            oConstants.poPulsars[oConstants.poSources[s].oAmpAcc.nPulsar].mdGCG[i][j] +=
              double(oConstants.poSources[s].oAmpAcc.mdGCaG[i][j]) *
              exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+nParIndex]) *
              exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+nParIndex]);
          } // for j
        } // for i
      } // if nPulsar
//      PrintMatrix(oConstants.poPulsars[oConstants.poSources[s].oAmpAcc.nPulsar].mdGCG);
    } else if(oConstants.poSources[s].oIntAcc.bSet) {
      // The interpolation acceleration is set. Interpolate the signal's
      // covariance matrix first.
      CalcPlIntMatrix(
	  oConstants.poSources[s].oIntAcc.ppIntAccel,
	  oConstants.poSources[s].oIntAcc.ppIntSpline,
	  oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+1],
	  oConstants.poSources[s].oIntAcc.mdGCiG);
      if(oConstants.poSources[s].oIntAcc.nPulsar < 0) {
        for(int i=0; i<oData.mdGCG.m_pnDimSize[0]; ++i) {
          for(int j=0; j<oData.mdGCG.m_pnDimSize[1]; ++j) {
            oData.mdGCG[i][j] += double(oConstants.poSources[s].oIntAcc.mdGCiG[i][j]) *
              exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex]) *
              exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex]);
          } // for j
        } // for i
      } else {
        for(int i=0; i<oConstants.poPulsars[oConstants.poSources[s].oIntAcc.nPulsar].mdGCG.m_pnDimSize[0]; ++i) {
          for(int j=0; j<oConstants.poPulsars[oConstants.poSources[s].oIntAcc.nPulsar].mdGCG.m_pnDimSize[0]; ++j) {
            oConstants.poPulsars[oConstants.poSources[s].oIntAcc.nPulsar].mdGCG[i][j] +=
              double(oConstants.poSources[s].oIntAcc.mdGCiG[i][j]) *
              exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex]) *
              exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex]);
          } // for j
        } // for i
      } // if nPulsar
    } else {
#if 0
      if(! mdRedTempTot.Defined()) {
	mdRedTempTot.Initialize(oConstants.n, oConstants.n);
	for(int i=0; i<oConstants.n; ++i) {
	  for(int j=0; j<oConstants.n; ++j) {
	    mdRedTempTot[i][j] = 0;
	  } // for j
	} // for i
      } // if mdRedTempTot

	StartParametersFromSources(oParametersTemp, oConstants.poSources, oConstants.nSources);
	ParametersFromStartParameters(oParametersTemp, oConstants.poSources, oConstants.nSources);
	mdTemp.Initialize(oConstants.n, oConstants.n);
	SetCoherenceMatrixOneSource(oData, oParametersTemp, oConstants, mdTemp, s, -1);
	for(int i=0; i<oConstants.n; ++i) {
	  for(int j=0; j<oConstants.n; ++j) {
	    mdRedTempTot[i][j] += double(mdTemp[i][j]);
	  } // for j
	} // for i
#endif

      // This also happens for the linear sources
      // There is no acceleration for this source. Do it the old-fashioned way:
      // brute force
//      printf("WARNING: brute-force reduction not implemented yet. Use interpolation!\n");
    } // if oAmpAcc, oIntAcc
  } // for s

#if 0
  if(mdRedTempTot.Defined()) {
	mdTemp2 = oData.mdGReduceT * mdRedTempTot;
	mdTemp = mdTemp2 * oData.mdGReduce;
	for(int i=0; i<oData.mdGCG.m_pnDimSize[0]; ++i) {
	  for(int j=0; j<oData.mdGCG.m_pnDimSize[1]; ++j) {
	    oData.mdGCG[i][j] += double(mdTemp[i][j]);
	  } // for j
	} // for i
  } // if mdRedtempTot
#endif

#ifdef DO_PROFILING
      tEndA2 = clock();
      tStartA3 = tEndA2;
#endif

  /* Add the pulsar GCG to the total GCG's */
  nAddIndex = 0;
  for(int a=0; a<oConstants.k; ++a) {
    for(int i=0; i<oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1]; ++i) {
      for(int j=0; j<oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1]; ++j) {
	oData.mdGCG[nAddIndex+i][nAddIndex+j] += double(oConstants.poPulsars[a].mdGCG[i][j]);
      } /* for j */
    } /* for i */
    nAddIndex += oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1];
  } /* for a */

#ifdef DO_PROFILING
      tEndA3 = clock();
      tStartA4 = tEndA3;
#endif

  /* The LL is now that of a 'normal' RGP: */
  mdTemp = oData.mdGCG.InverseChol(&dLogDetGCG);

#ifdef DO_PROFILING
      tEndA4 = clock();
      tStartB = tEndA4;
#endif

  dReturnValue = 0.5 * oData.mdGCG.m_pnDimSize[0]*log(2*M_PI);
  dReturnValue += 0.5*( oData.vdDataReduced * mdTemp * oData.vdDataReduced );
  dReturnValue += 0.5*(dLogDetGCG);

#ifdef DO_PROFILING
  tEndB = clock();

  nEvals++;
  if(nEvals > 1) {
    tA1 += tEndA1 - tStartA1;
    tA2 += tEndA2 - tStartA2;
    tA3 += tEndA3 - tStartA3;
    tA4 += tEndA4 - tStartA4;
    tB += tEndB - tStartB;

    if(nEvals % 30 == 0) {
      FILE *pFile;
      try {
	if(! (pFile = fopen("profile.txt", "w+")) ) throw 1;

	fprintf(pFile, "%i   %.3f   %.3f   %.3f   %.3f  %.3f\n",
	    nEvals-1,
	    double(tA1) / CLOCKS_PER_SEC,
	    double(tA2) / CLOCKS_PER_SEC,
	    double(tA3) / CLOCKS_PER_SEC,
	    double(tA4) / CLOCKS_PER_SEC,
	    double(tB) / CLOCKS_PER_SEC);

	fprintf(pFile, "Dimensions mdGCG: (%i x %i)\n", oData.mdGCG.m_pnDimSize[0], oData.mdGCG.m_pnDimSize[1]);

	if(fclose(pFile) ) throw 0;
      } catch(int nError) {
	printf("Error number: %i\n", nError);
      } // try
    } // if nEvals
  } /* if nEvals */
#endif

  return dReturnValue;
} /* LogLikelihoodReducedBasis */



#if 0

// #define DO_PROFILING
/* The likelihood function for reduced bases
 */
double LogLikelihoodReducedBasisOld(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  CMatrix mdTemp, mdTemp2;
  double dLogDetGCG, dReturnValue, dTemp;
  int nAddIndex;
  int nIndex1, nIndex2;
//  static bool bFirstRun=true;
//  static int nGWBPar=-1, nAllPar=-1, nClockPar=-1;
#ifdef DO_PROFILING
  static int nEvals = 0;

  static time_t tA1=0, tA2=0, tA3=0, tA30=0, tA4=0, tB=0, tC=0, tF=0, tT=0;
  time_t tStartA1, tEndA1,
	 tStartA2, tEndA2,
	 tStartA30, tEndA30,
	 tStartA3, tEndA3,
	 tStartA4, tEndA4,
	 tStartB, tEndB,
	 tStartC, tEndC,
	 tStartF, tEndF,
	 tStartT, tEndT;
#endif

  /* Set the GWB and the "All" sources */
  if(oData.bRedLikFirstRun) {
    for(int s=0; s<oConstants.nSources; ++s) {
      if(oConstants.poSources[s].eCorrelation == CID_GR) {
	oData.nRedLikGWBPar = oConstants.poSources[s].nFirstParIndex;
	break;
      } /* if eCorrelation */
    } /* for s */
    for(int s=0; s<oConstants.nSources; ++s) {
      if(oConstants.poSources[s].eCorrelation == CID_All) {
	oData.nRedLikAllPar = oConstants.poSources[s].nFirstParIndex;
	break;
      } /* if eCorrelation */
    } /* for s */
    for(int s=0; s<oConstants.nSources; ++s) {
      if(oConstants.poSources[s].eCorrelation == CID_Uniform) {
	oData.nRedLikClockPar = oConstants.poSources[s].nFirstParIndex;
	break;
      } /* if eCorrelation */
    } /* for s */
  } /* if bRedLikFirstRun */

#ifdef DO_PROFILING
    tStartA1 = clock();
    tStartA2 = tStartA1;
    tEndA1 = tStartA1;
    tEndA2 = tStartA1;
    tStartT = tStartA1;
#endif

  /* Set all the individual pulsar covariance matrices */
  for(int a=0; a<oConstants.k; ++a) {
    if(oConstants.poPulsars[a].bInterpolatorsSet) {
      /* We can greatly accelerate this by using the interpolators */
      if(oConstants.poPulsars[a].nParEfacIndex < 0) {
	fprintf(stderr, "ERROR: nParEfacIndex not set. This should not happen!!\n");
	return 0;
      } /* if nParEfacIndex */

      if(oConstants.poPulsars[a].nParPlIndex >= 0) {
	CalcPlIntMatrix(
	    oConstants.poPulsars[a].ppIntAccel,
	    oConstants.poPulsars[a].ppIntSpline,
	    oParameters.pdPar[oConstants.poPulsars[a].nParPlIndex+1],
	    oConstants.poPulsars[a].mdGCplG);
      } else if(oData.nRedLikAllPar >= 0) {
	/* If the noise is set the same for all pulsars */
	CalcPlIntMatrix(
	    oConstants.poPulsars[a].ppIntAccel,
	    oConstants.poPulsars[a].ppIntSpline,
	    oParameters.pdPar[oData.nRedLikAllPar+1],
	    oConstants.poPulsars[a].mdGCplG);
      } /* if nSourcePl */

      for(int i=0; i<oConstants.poPulsars[a].mdGCG.m_pnDimSize[0]; ++i) {
	for(int j=0; j<oConstants.poPulsars[a].mdGCG.m_pnDimSize[1]; ++j) {
	  oConstants.poPulsars[a].mdGCG[i][j] =
	    double(oConstants.poPulsars[a].mdGCefG[i][j]) *
	    oParameters.pdPar[oConstants.poPulsars[a].nParEfacIndex] *
	    oParameters.pdPar[oConstants.poPulsars[a].nParEfacIndex];

	  if(oConstants.poPulsars[a].nParEquadIndex >= 0) {
	    oConstants.poPulsars[a].mdGCG[i][j] +=
	      double(oConstants.poPulsars[a].mdGCeqG[i][j]) *
	      oParameters.pdPar[oConstants.poPulsars[a].nParEquadIndex] *
	      oParameters.pdPar[oConstants.poPulsars[a].nParEquadIndex];
	  } /* if nSourceEquad */

	  if(oConstants.poPulsars[a].nParPlIndex >= 0) {
	    oConstants.poPulsars[a].mdGCG[i][j] +=
	      double(oConstants.poPulsars[a].mdGCplG[i][j]) *
	      oParameters.pdPar[oConstants.poPulsars[a].nParPlIndex] *
	      oParameters.pdPar[oConstants.poPulsars[a].nParPlIndex];
	  } else if(oData.nRedLikAllPar >= 0) {
	    oConstants.poPulsars[a].mdGCG[i][j] +=
	      double(oConstants.poPulsars[a].mdGCplG[i][j]) *
	      oParameters.pdPar[oData.nRedLikAllPar] *
	      oParameters.pdPar[oData.nRedLikAllPar];
	  } /* if nSourcePl */
	} /* for j */
      } /* for i */
#ifdef DO_PROFILING
      tEndA2 = clock();
#endif
    } else {
      /* Set the covariance matrix */
      SetCoherenceMatrixSinglePulsarNoCorrSig(oData, oParameters, oConstants, a);
#ifdef DO_PROFILING
      tEndA1 = clock();
      tStartA2 = tEndA1;
#endif

      /* Calculate the GCG products */
      mdTemp = oConstants.poPulsars[a].mdC * oConstants.poPulsars[a].mdGReduce;
      oConstants.poPulsars[a].mdGCG = oConstants.poPulsars[a].mdGReduceT * mdTemp;
#ifdef DO_PROFILING
      tEndA2 = clock();
#endif
    } /* if bInterpolatorsSet */
  } /* for a */

  /* Now calculate the full GWB / clock product */
  if(oData.bRedLikFirstRun) {
    if(oData.bGWBInterpolatorsSet) {
      oData.mdGCGCorr.Initialize(oData.mdGReduce.m_pnDimSize[1], oData.mdGReduce.m_pnDimSize[1]);
    } else if(oData.nRedLikGWBPar >= 0 && oConstants.bCorAmpAccel) {
      dTemp = oParameters.pdPar[oData.nRedLikGWBPar];
      oParameters.pdPar[oData.nRedLikGWBPar] = 1.0;
      SetCoherenceMatrixCorrSig(oData, oParameters, oConstants);
      mdTemp = oData.mdCCor * oData.mdGReduce;
      oData.mdGCGCorr = oData.mdGReduceT * mdTemp;
      oParameters.pdPar[oData.nRedLikGWBPar] = dTemp;
    } else {
      oData.mdGCGCorr.Initialize(oData.mdGReduce.m_pnDimSize[1], oData.mdGReduce.m_pnDimSize[1]);
    }/* if oData.nRedLikGWBPar */

    oData.mdGCG.Initialize(oData.mdGReduce.m_pnDimSize[1], oData.mdGReduce.m_pnDimSize[1]);
    oData.bRedLikFirstRun = false;
  } /* if bRedLikFirstRun */

#ifdef DO_PROFILING
  tStartA3 = clock();
  tStartA30 = tStartA3;
  tStartA4 = tStartA3;
  tEndA3 = tStartA3;
  tEndA30 = tStartA3;
  tEndA4 = tStartA3;
#endif
  if(oConstants.bCorAmpAccel) {
    if(oData.bGWBInterpolatorsSet) {
      CalcPlIntMatrix(
	  oData.ppGWBIntAccel,
	  oData.ppGWBIntSpline,
	  oParameters.pdPar[oData.nRedLikGWBPar+1],
	  oData.mdGCGCorr);
#ifdef DO_PROFILING
      tEndA30 = clock();
      tStartA3 = tEndA30;
#endif
    } /* if bGWBInterpolatorsSet */

    /* Add the GWB signal to the mdGCG */
    for(int i=0; i<oData.mdGReduce.m_pnDimSize[1]; ++i) {
      for(int j=0; j<oData.mdGReduce.m_pnDimSize[1]; ++j) {
	if(oData.nRedLikGWBPar >= 0) {
	  oData.mdGCG[i][j] = double(oData.mdGCGCorr[i][j]) * (oParameters.pdPar[oData.nRedLikGWBPar] * oParameters.pdPar[oData.nRedLikGWBPar]);
	} else {
	  oData.mdGCG[i][j] = 0;
	}/* if oData.nRedLikGWBPar */
      } /* for j */
    } /* for i */

    /* Add the clock signal to the mdGCG */
    if(oData.bClockInterpolatorsSet && oData.nRedLikClockPar >= 0) {
      CalcPlIntMatrix(
	  oData.ppClockIntAccel,
	  oData.ppClockIntSpline,
	  oParameters.pdPar[oData.nRedLikClockPar+1],
	  oData.mdGCGCorr);

      for(int i=0; i<oData.mdGReduce.m_pnDimSize[1]; ++i) {
	for(int j=0; j<oData.mdGReduce.m_pnDimSize[1]; ++j) {
	  oData.mdGCG[i][j] = double(oData.mdGCGCorr[i][j]) * (oParameters.pdPar[oData.nRedLikClockPar] * oParameters.pdPar[oData.nRedLikClockPar]);
	} /* for j */
      } /* for i */
    } /* if bClockInterpolatorsSet */
#ifdef DO_PROFILING
    tEndA3 = clock();
#endif
  } else {
    SetCoherenceMatrixCorrSig(oData, oParameters, oConstants);
#ifdef DO_PROFILING
    tEndA3 = clock();
    tStartA4 = clock();
#endif
    mdTemp = oData.mdCCor * oData.mdGReduce;
    oData.mdGCG = oData.mdGReduceT * mdTemp;
#ifdef DO_PROFILING
    tEndA4 = clock();
#endif
  } /* if bCorAmpAccel */

#ifdef DO_PROFILING
  tStartB = clock();
#endif
  /* Add the pulsar covariances */
  nAddIndex = 0;
  for(int a=0; a<oConstants.k; ++a) {
    for(int i=0; i<oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1]; ++i) {
      for(int j=0; j<oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1]; ++j) {
	oData.mdGCG[nAddIndex+i][nAddIndex+j] += double(oConstants.poPulsars[a].mdGCG[i][j]);
      } /* for j */
    } /* for i */
    nAddIndex += oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1];
  } /* for a */

#ifdef DO_PROFILING
  tEndB = clock();
  tStartC = tEndB;
#endif

  /* The LL is now that of a 'normal' RGP: */
  mdTemp = oData.mdGCG.InverseChol(&dLogDetGCG);

#ifdef DO_PROFILING
  tEndC = clock();
  tStartF = tEndC;
#endif

  dReturnValue = 0.5 * oData.mdGCG.m_pnDimSize[0]*log(2*M_PI);
  dReturnValue += 0.5*( oData.vdDataReduced * mdTemp * oData.vdDataReduced );
  dReturnValue += 0.5*(dLogDetGCG);


#ifdef DO_PROFILING
  tEndF = clock();
  tEndT = tEndF;

  nEvals++;
  if(nEvals > 1) {
    tA1 += tEndA1 - tStartA1;
    tA2 += tEndA2 - tStartA2;
    tA30 += tEndA30 - tStartA30;
    tA3 += tEndA3 - tStartA3;
    tA4 += tEndA4 - tStartA4;
    tB += tEndB - tStartB;
    tC += tEndC - tStartC;
    tF += tEndF - tStartF;
    tT += tEndT - tStartT;

    /*
    fprintf(stderr, "\n%i   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f\n",
	nEvals-1,
	double(tA1) / CLOCKS_PER_SEC,
	double(tA2) / CLOCKS_PER_SEC,
	double(tA3) / CLOCKS_PER_SEC,
	double(tA4) / CLOCKS_PER_SEC,
	double(tB) / CLOCKS_PER_SEC,
	double(tC) / CLOCKS_PER_SEC,
	double(tF) / CLOCKS_PER_SEC,
	double(tT) / CLOCKS_PER_SEC);
	*/

    FILE *pFile;
    try {
      if(! (pFile = fopen("profile.txt", "w+")) ) throw 1;

      fprintf(pFile, "%i   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f\n",
	  nEvals-1,
	  double(tA1) / CLOCKS_PER_SEC,
	  double(tA2) / CLOCKS_PER_SEC,
	  double(tA30) / CLOCKS_PER_SEC,
	  double(tA3) / CLOCKS_PER_SEC,
	  double(tA4) / CLOCKS_PER_SEC,
	  double(tB) / CLOCKS_PER_SEC,
	  double(tC) / CLOCKS_PER_SEC,
	  double(tF) / CLOCKS_PER_SEC,
	  double(tT) / CLOCKS_PER_SEC);

      if(fclose(pFile) ) throw 0;
    } catch(int nError) {
      printf("Error number: %i\n", nError);
    } // try
  } /* if nEvals */
#endif

  return dReturnValue;
} /* LogLikelihoodReducedBasisOld */
#endif

/* Calculate the LogLikelihood value. Outline:
**
** - x^T(invC - invC*M*inv(M^T*invC*M)*M^T*invC )x
**
** Last updated: 2007-08-23
*/
double LogLikelihood(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dLogDetC=0, dLogDetMCM = 0, dReturnValue=0;
  CMatrix mdTemp;
  CVector vdNewData, vdDetResiduals;
#ifdef DO_PROFILING
  static int nEvals = 0;

  static time_t tA0=0, tC=0, tD=0, tE=0, tF=0, tT=0;
  time_t tStartA0, tEndA0,
	 tStartC, tEndC,
	 tStartD, tEndD,
	 tStartE, tEndE,
	 tStartF, tEndF,
	 tStartT, tEndT;
#endif
  // oConstants.bUseReducedBasis = false;

  if(! oConstants.bUseReducedBasis) {

#ifdef DO_PROFILING
    tStartA0 = clock();
    tStartT = tStartA0;
#endif
    // First set the new coherence matrix
    SetCoherenceMatrix(oData, oParameters, oConstants);
#ifdef DO_PROFILING
    tEndA0 = clock();

    tStartC = tEndA0;
#endif

    // Calculate the inverse and the determinant of the matrix C
    oData.mdInvC = oData.mdC.InverseChol(&dLogDetC);
  //  oData.mdInvC = oData.mdC.Inverse();

#ifdef DO_PROFILING
    tEndC = clock();
    tStartD = tEndC;
#endif
    // Build the marginalization matrices
    mdTemp = oData.mdBlockTrans * oData.mdInvC * oData.mdBlock;

#ifdef DO_PROFILING
    tEndD = clock();
    tStartE = tEndD;
#endif

    mdTemp.InvertChol(&dLogDetMCM);

#ifdef DO_PROFILING
    tEndE = clock();
    tStartF = tEndD;
#endif

    // Subtract the deterministic sources from the residuals
    vdDetResiduals = oData.vdData;
    ResidualsFromDetSources(oConstants, oParameters, oData, &vdDetResiduals);
    vdNewData = oData.vdData - vdDetResiduals;

    // Save this for the MCMC
    oData.mdCXiInv = mdTemp;
    oData.vdChi = mdTemp * (oData.mdBlockTrans * (oData.mdInvC * vdNewData));

    // Calculate the LogLikelihood
    dReturnValue += 0.5 * (oConstants.n - oConstants.nMarParameters)*log(2*M_PI);
    dReturnValue += 0.5*( vdNewData * oData.mdInvC * vdNewData );
    dReturnValue -= 0.5*(vdNewData*(oData.mdInvC*(oData.mdBlock*(mdTemp*(oData.mdBlockTrans*(oData.mdInvC*vdNewData))))));
    dReturnValue += 0.5*(dLogDetC + dLogDetMCM);

#ifdef DO_PROFILING
    tEndF = clock();
    tEndT = tEndF;

    nEvals++;
    if(nEvals > 1) {
      tA0 += tEndA0 - tStartA0;
      tC += tEndC - tStartC;
      tD += tEndD - tStartD;
      tE += tEndE - tStartE;
      tF += tEndF - tStartF;
      tT += tEndT - tStartT;

      /*
      fprintf(stderr, "\n%i   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f\n",
	  nEvals-1,
	  double(tA0) / CLOCKS_PER_SEC,
	  double(tC) / CLOCKS_PER_SEC,
	  double(tD) / CLOCKS_PER_SEC,
	  double(tE) / CLOCKS_PER_SEC,
	  double(tF) / CLOCKS_PER_SEC,
	  double(tT) / CLOCKS_PER_SEC);
	  */
      FILE *pFile;
      try {
	if(! (pFile = fopen("profile.txt", "w+")) ) throw 1;

	fprintf(pFile, "%i   %.3f   %.3f   %.3f   %.3f   %.3f   %.3f\n",
	    nEvals-1,
	    double(tA0) / CLOCKS_PER_SEC,
	    double(tC) / CLOCKS_PER_SEC,
	    double(tD) / CLOCKS_PER_SEC,
	    double(tE) / CLOCKS_PER_SEC,
	    double(tF) / CLOCKS_PER_SEC,
	    double(tT) / CLOCKS_PER_SEC);

	if(fclose(pFile) ) throw 0;
      } catch(int nError) {
	printf("Error number: %i\n", nError);
      } // try

    } /* if nEvals */
#endif

  } else {
    dReturnValue = LogLikelihoodReducedBasis(oData, oParameters, oConstants);
  } /* if bUseReducedBasis */

  return dReturnValue;
} // LogLikelihood


/* Calculate the LogLikelihood value for an ensemble of datasets. Outline:
**
** - x^T(invC - invC*M*inv(M^T*invC*M)*M^T*invC )x
**
** Last updated: 2011-08-04
*/
void LogLikelihoodEnsemble(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dLogDetC=0, dLogDetMCM = 0;
  CMatrix mdTemp, mdTemp2;
  CMatrix mdNewData;;
  CVector vdDetResiduals, vdNewData, vdTemp1, vdTemp2;
  CNumber ndNum;
#if 0
  time_t tTimeNow, tTimePrevious;
  double dTimeSeconds;
  tTimePrevious = clock();
#endif

  if(! oData.mdDataSets.Defined() ) {
    fprintf(stderr, "WARNING: mdDataSets not defined!\n");
    return;
  } // if Defined

  if(! (oData.mdChi.Defined() && oData.vdLL.Defined() && oData.mdBlockTrans.Defined())) {
    oData.mdBlockTrans = oData.mdBlock[LO_TRANSPOSE];
    oData.mdChi.Initialize(oConstants.nMarParameters, oData.nDataSets);
    oData.vdLL.Initialize(oData.nDataSets);
  } // if Defined

  // First set the new coherence matrix
  SetCoherenceMatrix(oData, oParameters, oConstants);

  // Calculate the inverse and the determinant of the matrix C
  oData.mdInvC = oData.mdC.InverseChol(&dLogDetC);

  // Build the marginalization matrices
  mdTemp2 = oData.mdInvC * oData.mdBlock;
  mdTemp = oData.mdBlockTrans * mdTemp2;

  // Calculate the inverse and the determinant of the matrix M^{T}C^{-1}M
  mdTemp.InvertChol(&dLogDetMCM);

  // Subtract the deterministic sources from the residuals
  vdDetResiduals.Initialize(oConstants.n);
  mdNewData.Initialize(oConstants.n, oData.nDataSets);
  vdNewData.Initialize(oConstants.n);
  ResidualsFromDetSources(oConstants, oParameters, oData, &vdDetResiduals);
  for(int i=0; i<oConstants.n; i++) {
    for(int j=0; j<oData.nDataSets; j++) {
      mdNewData[i][j] = double(oData.mdDataSets[i][j]) - double(vdDetResiduals[i]);
    } // for j
  } // for i

  // Save this for the MCMC
  oData.mdCXiInv = mdTemp;
  for(int j=0; j<oData.nDataSets; j++) {
    for(int i=0; i<oConstants.n; i++) {
      vdNewData[i] = double(mdNewData[i][j]);
    } // for i
    // Save this for the MCMC
    vdTemp1 = oData.mdInvC * vdNewData;
    vdTemp2 = oData.mdBlockTrans * vdTemp1;
    oData.vdChi = mdTemp * vdTemp2;
//    oData.vdChi = mdTemp * (oData.mdBlockTrans * (oData.mdInvC * vdNewData));
    for(int i=0; i<oConstants.nMarParameters; i++) {
      oData.mdChi[i][j] = double(oData.vdChi[i]);
    } // for i

    // The LogLikelihood
    oData.vdLL[j] = 0.5 * (oConstants.n - oConstants.nMarParameters)*log(2*M_PI);
    vdTemp1 = oData.mdInvC * vdNewData;

    for(int i=0; i<oConstants.n; i++) {
      oData.vdLL[j] += 0.5 * double(vdNewData[i]) * double(vdTemp1[i]);
    } // for j
//    ndNum = vdNewData * vdTemp1;
//    oData.vdLL[j] += 0.5 * double(ndNum);

    vdTemp2 = oData.mdBlockTrans * vdTemp1;
    vdTemp1 = mdTemp * vdTemp2;
    vdTemp2 = oData.mdBlock * vdTemp1;
    vdTemp1 = oData.mdInvC * vdTemp2;
//    ndNum = vdNewData * vdTemp1;
//    oData.vdLL[j] -= 0.5 * double(ndNum);
    for(int i=0; i<oConstants.n; i++) {
      oData.vdLL[j] -= 0.5 * double(vdNewData[i]) * double(vdTemp1[i]);
    } // for j
//    oData.vdLL[j] -= 0.5*(vdNewData*(oData.mdInvC*(oData.mdBlock*(mdTemp*(oData.mdBlockTrans*(oData.mdInvC*vdNewData))))));
    oData.vdLL[j] += 0.5*(dLogDetC + dLogDetMCM);
  } // for j

#if 0
  tTimeNow = clock();
  dTimeSeconds = double(tTimeNow - tTimePrevious)/CLOCKS_PER_SEC;
  printf("Time elapsed: %f seconds  (size %i x %i)\n", dTimeSeconds, oData.mdInvC.m_pnDimSize[0], oData.mdInvC.m_pnDimSize[1]);
#endif
} // LogLikelihoodEnsemble



#ifdef HAVE_MPI
/* Performs a MCMC integration with a Gaussian proposal distribution. The
 * final chain is written to a datafile: strFileName.
 *
 * Todo: during the burnin-time, the random walkers do walk, but the proposal
 * distribution width is not adjusted. The process of setting the width of
 * the distribution to yield an optimal acceptance rate must be automized.
 *
 * This version of MCMCGauss serves as a server in the client/server MPI model
 *
 * TODO: This function does not send QUIT msgs to the clients anymore. Remove
 * them.
 * */
void MCMCGaussServer(MPI_Datatype *pmpiParameters, SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameterNumber) {
  // MPI variables
  SMPIParametersType *poMPIParameters;
  MPI_Request *pmpiEvalRequest, *pmpiMsgRequest;
  MPI_Status *pmpiEvalStatus, *pmpiMsgStatus;
  int *pnEvalFlag, *pnMsgFlag, *pnEvalError, *pnMsgError;
  bool *pbProcessHasJob, *pbProcessTerminated, *pbProcessFirstEval;
  bool bContinue=true;
  int nIndex=0, nIndexFit=0;
  int nMsg, nTemp, nParCount, nUpdateScreen=0;
  int nKey=-1;
  int nProcesses, nProcessRank, nNameLength;
  char strProcessorName[MPI_MAX_PROCESSOR_NAME], strMsg[160], strFileName[160];

  // MCMC variables
  double dAcceptRatio=0;
  int nMCMCSteps, nPlotPoints, *pnStep, *pnStepIndex;	// Total steps, number of plotting points on a graph, and step number
  MCMCDataPoint **ppDat;			// All data aqcuired so far
  bool *pbBurnIn;				// Are we in the burn-in period?
  int *pnCheckSteps;				// Amount of steps, accepted or not accepted, since last accept-check
  int *pnAcceptedSteps;				// Amount of accepted steps since last accept-check
  double *pdAcceptRatio;			// Likelihood-ratio

  char strBuf[100];
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);	// Random number generator
  gsl_rng_set(rng, (unsigned int)time(NULL));

  CVector vdThetaMin, vdThetaMax, vdThetaWidth, vdTheta;

  vdTheta.Initialize(oConstants.nFitParameters);
  vdThetaMin.Initialize(oConstants.nFitParameters);
  vdThetaMax.Initialize(oConstants.nFitParameters);
  vdThetaWidth.Initialize(oConstants.nFitParameters);

  // Initialize MPI environment
  MPI_Comm_size(MPI_COMM_WORLD,&nProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &nProcessRank);
  MPI_Get_processor_name(strProcessorName, &nNameLength);

  sprintf(strMsg, "Starting MPI Server on %s, %i of %i",
      strProcessorName, nProcessRank, nProcesses);
  PrintStatus(strMsg);

  // For each process we need to keep track of some data
  pbProcessHasJob = new bool[nProcesses-1];
  pbProcessTerminated = new bool[nProcesses-1];
  pbProcessFirstEval = new bool[nProcesses-1];
  poMPIParameters = new SMPIParametersType[nProcesses-1];
  pmpiEvalRequest = new MPI_Request[nProcesses-1];
  pmpiMsgRequest = new MPI_Request[nProcesses-1];
  pmpiEvalStatus = new MPI_Status[nProcesses-1];
  pmpiMsgStatus = new MPI_Status[nProcesses-1];
  pnEvalFlag = new int[nProcesses-1];
  pnMsgFlag = new int[nProcesses-1];
  pnEvalError = new int[nProcesses-1];
  pnMsgError = new int[nProcesses-1];

  // Initialize the pointers other
  nPlotPoints = oConstants.nPlotPoints;
  nMCMCSteps = oConstants.nMCMCSteps;
  ppDat = new MCMCDataPoint*[nProcesses-1];
  for(int i=0; i<nProcesses-1; i++)
    ppDat[i] = new MCMCDataPoint[MAX_MCMC_BUFFER];
  pbBurnIn = new bool[nProcesses-1];
  pnCheckSteps = new int[nProcesses-1];
  pnAcceptedSteps = new int[nProcesses-1];
  pnStep = new int[nProcesses-1];
  pnStepIndex = new int[nProcesses-1];
  pdAcceptRatio = new double[nProcesses-1];

  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcat(strFileName, "mcmcdata.mpitotal.dat");

  PrintSuccess();

  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  VectorFromFitParameters(oParameters, oConstants, vdTheta);
  MinVectorFromFitParameters(oParameters, oConstants, vdThetaMin);
  MaxVectorFromFitParameters(oParameters, oConstants, vdThetaMax);
  WidthMCMCVectorFromFitParameters(oParameters, oConstants, vdThetaWidth);

  // Start the MCMC Chains
  PrintUpdate("Clients are generating MCMC points...");

  // Initialize the chains
  for(int i=0; i<nProcesses-1; i++) {
    StartParametersFromSources(ppDat[i][0].oParameters, oConstants.poSources, oConstants.nSources);
    ParametersFromStartParameters(ppDat[i][0].oParameters, oConstants.poSources, oConstants.nSources);
    if(i==0) {
      ppDat[0][0].dLogLik = LogLikelihoodTimesPrior(oData, ppDat[0][0].oParameters, oConstants);
    } else {
      ppDat[i][0].dLogLik = ppDat[0][0].dLogLik;
    } // if i
  } // for i

  // Initialize all MPI-server parameters
  for(int i=0; i<oConstants.nFitParameters; i++)
    poMPIParameters[0].pdPar[i] = double(vdTheta[i]);
  poMPIParameters[0].dLogLik = ppDat[0][0].dLogLik;
  poMPIParameters[0].nStatus = MPI_BAN_EVAL_SUCCESS;
  //MPI_Send( a, 100, MPI_DOUBLE, 1, 17, MPI_COMM_WORLD );
  for(int i=0; i<nProcesses-1; i++) {
    if(i != 0) poMPIParameters[i] = poMPIParameters[0];
    pnEvalFlag[i] = 0;
    pnMsgFlag[i] = 0;
    pnEvalError[i] = 0;
    pnMsgError[i] = 0;
    pbProcessHasJob[i] = false;
    pbProcessFirstEval[i] = true;
    pbProcessTerminated[i] = false;
    pbBurnIn[i] = true;
    pnCheckSteps[i] = 1;
    pnAcceptedSteps[i] = 1;
    pnStep[i] = 1;
    pnStepIndex[i] = 1;
    pdAcceptRatio[i] = 0;
  } // for i

  // Write the header of the MCMC parameter file
  WriteMCMCDataFileStart(oConstants.nParameters, strFileName);

  // Continue the conversation with the clients until we are done
  while(bContinue) {
    // First check whether we need to give orders to clients
    for(int i=0; i<nProcesses-1; i++) {
      if(pbProcessHasJob[i] == false && pbProcessTerminated[i] == false && pnStep[i] < nMCMCSteps) {
	// We need to give this process a new job. First initialize the
	// parameters we should send
	pbProcessHasJob[i] = true;

	ppDat[i][pnStepIndex[i]].oParameters = ppDat[0][0].oParameters;
	nIndex=0;
	nIndexFit=0;
	// Generate new parameters for the next proposal
	for(int s=0; s<oConstants.nSources; s++) {
	  for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
	    nIndex = oConstants.poSources[s].nFirstParIndex + p;
	    if(SourceWorks(oConstants, s) && (oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p))) {
	      for(;;) { // Keep checking for a good value
		poMPIParameters[i].pdPar[nIndexFit] = ppDat[i][pnStepIndex[i]-1].oParameters.pdPar[nIndex] + double(vdThetaWidth[nIndexFit])*oConstants.dGlobalMCMCWidthFactor*gsl_ran_ugaussian(rng);
		ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] = poMPIParameters[i].pdPar[nIndexFit];

		if(oConstants.poSources[s].oSourceType.nWrapTag & FITPARAMETER(p)) {
		  // Check whether we need to wrap this parameter
		  if(ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] < ppDat[i][pnStepIndex[i]].oParameters.pdParMinBound[nIndex]) {
		    ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] = ppDat[i][pnStepIndex[i]].oParameters.pdParMaxBound[nIndex] - (ppDat[i][pnStepIndex[i]].oParameters.pdParMinBound[nIndex] - ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex]);
		  } // if pdPar < pdParMinBound

		  if(ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] > ppDat[i][pnStepIndex[i]].oParameters.pdParMaxBound[nIndex]) {
		    ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] = ppDat[i][pnStepIndex[i]].oParameters.pdParMinBound[nIndex] + (ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] - ppDat[i][pnStepIndex[i]].oParameters.pdParMaxBound[nIndex]);
		  } // if pdPar < pdParMinBound
		} // if nWrapTag


		if(ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] >= double(vdThetaMin[nIndexFit]) && ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] <= double(vdThetaMax[nIndexFit]))
		  break;
	      } // for
	      nIndexFit++;
	    } // if nFitTag
	  } // for p
	} // for s

	poMPIParameters[i].dLogLik = 0;
	poMPIParameters[i].nStatus = MPI_BAN_EVAL_SUCCESS;

	MPI_Send(poMPIParameters+i, 1, *pmpiParameters, i+1, MPI_BAN_TAG_EVAL, MPI_COMM_WORLD);

	// Directly also post a non-blocking receive for this value
	pnEvalError[i] = MPI_Irecv(poMPIParameters+i, 1, *pmpiParameters, i+1, MPI_BAN_TAG_EVAL, MPI_COMM_WORLD, pmpiEvalRequest+i);
      } // if pbProcessHasJob
    } // for i

    // Now check whether some clients are done with their evaluation
    for(int i=0; i<nProcesses-1; i++) {
      pnEvalError[i] = MPI_Test(pmpiEvalRequest+i, pnEvalFlag+i, pmpiEvalStatus+i);
      // Check whether a client is done with its calculations
      if(pnEvalFlag[i] && pbProcessHasJob[i]) {
	pbProcessHasJob[i] = false;

	// Retrieve the likelihood value
	ppDat[i][pnStepIndex[i]].dLogLik = poMPIParameters[i].dLogLik;

	if(! pbProcessFirstEval[i]) {
	  // Now accept or reject this point
	  dAcceptRatio = exp(ppDat[i][pnStepIndex[i]-1].dLogLik - ppDat[i][pnStepIndex[i]].dLogLik);
	  if(gsl_rng_uniform(rng) < dAcceptRatio) {
	    pnAcceptedSteps[i]++;
	  }  else {
	    // Reject the point
	    for(int s=0; s<oConstants.nSources; s++) {
	      for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
		nIndex = oConstants.poSources[s].nFirstParIndex + p;
		if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
		    ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] = ppDat[i][pnStepIndex[i]-1].oParameters.pdPar[nIndex];
		} // if nFitTag
	      } // for p
	    } // for s
	    ppDat[i][pnStepIndex[i]].dLogLik = ppDat[i][pnStepIndex[i]-1].dLogLik;
	  } // if gsl_rng_uniform(rng)

	  // We except this point, but first check whether we are still
	  // in the burn-in period.
	  if(pbBurnIn[i] && pnStep[i] >= oConstants.nBurnInSteps) {
	    ppDat[i][0].oParameters = ppDat[i][pnStepIndex[i]].oParameters;
	    ppDat[i][0].dLogLik = ppDat[i][pnStepIndex[i]].dLogLik;
	    pbBurnIn[i] = false;
	    pnStep[i] = 0; pnAcceptedSteps[i]=0;
	  } // bBurnIn

	  pnStep[i]++; pnStepIndex[i]++;


	  if(pnStepIndex[i] == MAX_MCMC_BUFFER) {
	    // We need to write this stuff to a file
	    if(! pbBurnIn[i])
	      WriteMCMCDataFileAppend(MAX_MCMC_BUFFER, oConstants.nParameters, ppDat[i], strFileName);
	    ppDat[i][0] = ppDat[i][pnStepIndex[i]-1];
	    pnStepIndex[i] = 1;
	  } // if pnStepIndex

	  // If we are in the burn-in period, check whether we continue with that
	  if(pbBurnIn[i]) {
	    if(pnCheckSteps[i] % 1000 == 0) { // This is the accept-check in the burn-in period
	      // Adjust the sample distribution width accordingly wrt the acceptance rate
	      pnAcceptedSteps[i] = 0; pnCheckSteps[i] = 0;
	    } // nCheckSteps
	  } // bBurnIn
	  pnCheckSteps[i]++;
	} else {
	  pbProcessFirstEval[i] = false;
	    // Reject the point
	    for(int s=0; s<oConstants.nSources; s++) {
	      for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
		nIndex = oConstants.poSources[s].nFirstParIndex + p;
		if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
		    ppDat[i][pnStepIndex[i]].oParameters.pdPar[nIndex] = ppDat[i][0].oParameters.pdPar[nIndex];
		} // if nFitTag
	      } // for p
	    } // for s
	    ppDat[i][pnStepIndex[i]].dLogLik = ppDat[i][0].dLogLik;
	} // if ! pbProcessFirstEval
      } // if nEvalFlag
    } // for i

    nTemp=0;
    while(pbProcessTerminated[nTemp]) nTemp++;
    if(nTemp == nProcesses-1) nTemp = 0;

    if(nUpdateScreen != pnStepIndex[nTemp]) {
      nUpdateScreen = pnStepIndex[nTemp];
      nParCount=0;
      if(pnStepIndex[nTemp] >= 2) {
	sprintf(strBuf, "C-%i:%3i%% %2i: ", nTemp, int(100.0*pnStep[nTemp]/nMCMCSteps), int(pnAcceptedSteps[nTemp]*100.0/pnStep[nTemp]));
	for(int s=0; s<oConstants.nSources && nParCount < 6; s++) {
	  for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters && nParCount < 6; p++) {
	    if(SourceWorks(oConstants, s) && (oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p))) {
	      nParCount++;
	      sprintf(strBuf, "%s%6.2e | ", strBuf, ppDat[nTemp][pnStepIndex[nTemp]-2].oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+p]);
	    } // if oConstants
	  } // for p
	} //for s
	PrintUpdate(strBuf);
      } // if pnStep
    } // if nUpdateScreen

    if(n_min(pnStep, nProcesses-1) < oConstants.nMCMCSteps)
      bContinue = true;
    else
      bContinue = false;

    // Check whether we should send some quit messages
    for(int i=0; i<nProcesses-1; i++) {
      if(pnStep[i] >= oConstants.nMCMCSteps) {
	if(! pbProcessTerminated[i]) {
	  pbProcessTerminated[i] = true;
	  nMsg = MPI_BAN_MSG_QUIT;
//	  MPI_Send(&nMsg, 1, MPI_INT, i+1, MPI_BAN_TAG_MSG, MPI_COMM_WORLD);

	  // Also write the MCMC data
	  WriteMCMCDataFileAppend(pnStepIndex[i]-1, oConstants.nParameters, ppDat[i], strFileName);
	  ppDat[i][0] = ppDat[i][pnStepIndex[i]-1];
	} // if pbProcessTerminated
      } // if pnStep[i]
    } // for i

    // Check whether the user has requested to stop
    // 27 = <ESC>, 113 = 'q', 81 = 'Q'
    // Doesn't work in combination with mpi (too bad)
    // nKey = KbHit();
    if(nKey == 27 || nKey == 113 || nKey == 81)
      bContinue = false;
  } // while bContinue

  // Send the kill to all the clients (just to be sure...)
  for(int i=0; i<nProcesses-1; i++) {
    if(! pbProcessTerminated[i]) {
      nMsg = MPI_BAN_MSG_QUIT;
//      MPI_Send(&nMsg, 1, MPI_INT, i+1, MPI_BAN_TAG_MSG, MPI_COMM_WORLD);
    } // if pbProcessTerminated
  } // for i

  if(nKey == 27 || nKey == 113 || nKey == 81) {
    printf("User has requested to stop the evaluation.\n");
    printf("A quit message has been sent to all the clients.\n");
  } // if nKey

  // Write all the MCMC data to file (Already done above)
//  WriteMCMCDataFileAppend(pnStepIndex[i]-1, oConstants.nParameters, ppDat[i], strFileName);
//  WriteMCMCDataFile(nProcesses-1, oConstants.nMCMCSteps, oConstants.nParameters, ppDat, strFileName);

  // Clean up memory
  delete[] pbProcessHasJob;
  delete[] pbProcessTerminated;
  delete[] pbProcessFirstEval;
  delete[] poMPIParameters;
  delete[] pmpiEvalRequest;
  delete[] pmpiMsgRequest;
  delete[] pmpiEvalStatus;
  delete[] pmpiMsgStatus;
  delete[] pnEvalFlag;
  delete[] pnMsgFlag;
  delete[] pnEvalError;
  delete[] pnMsgError;

  for(int i=0; i<nProcesses-1; i++)
    delete[] (ppDat[i]);
  delete[] ppDat;
  delete[] pbBurnIn;
  delete[] pnCheckSteps;
  delete[] pnAcceptedSteps;
  delete[] pnStep;
  delete[] pnStepIndex;
  delete[] pdAcceptRatio;
} // MCMCGaussServer



/* This function serves as a server that commands the clients to make a simple
 * 1D plot of the likelihoodfunction
 *
 * */
void MakePlotServer(MPI_Datatype *pmpiParameters, SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameterNumber) {
  CVector vdPlot, vdFunctionValues, vdTemp, vdThetaMin, vdThetaMax, vdTheta;
  int nPlotPoints;
  int k, l, nParameters, nCount;
  int nProcesses, nProcessRank, nNameLength;
  char strProcessorName[MPI_MAX_PROCESSOR_NAME], strMsg[160];

  SMPIParametersType *poMPIParameters;
  MPI_Request *pmpiEvalRequest, *pmpiMsgRequest;
  MPI_Status *pmpiEvalStatus, *pmpiMsgStatus;
  int *pnEvalFlag, *pnMsgFlag, *pnEvalError, *pnMsgError;
  bool *pbProcessHasJob, *pbProcessTerminated;
  bool bContinue=true;
  int nIndex=0;
  int nMsg;
  int nKey=-1;

  // Initialize MPI environment
  MPI_Comm_size(MPI_COMM_WORLD,&nProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &nProcessRank);
  MPI_Get_processor_name(strProcessorName, &nNameLength);

  sprintf(strMsg, "Starting MPI Server on %s, %i of %i",
      strProcessorName, nProcessRank, nProcesses);
  PrintStatus(strMsg);

  // For each process we need to keep track of some data
  pbProcessHasJob = new bool[nProcesses-1];
  pbProcessTerminated = new bool[nProcesses-1];
  poMPIParameters = new SMPIParametersType[nProcesses-1];
  pmpiEvalRequest = new MPI_Request[nProcesses-1];
  pmpiMsgRequest = new MPI_Request[nProcesses-1];
  pmpiEvalStatus = new MPI_Status[nProcesses-1];
  pmpiMsgStatus = new MPI_Status[nProcesses-1];
  pnEvalFlag = new int[nProcesses-1];
  pnMsgFlag = new int[nProcesses-1];
  pnEvalError = new int[nProcesses-1];
  pnMsgError = new int[nProcesses-1];
  for(int i=0; i<nProcesses-1; i++) {
    pbProcessHasJob[i] = false;
    pbProcessTerminated[i] = false;
  } // for i


  PrintSuccess();

  InitProgressBar("Making 1D plot...");
  DrawProgressBar(0);

  nParameters = oConstants.nFitParameters;
  vdThetaMin.Initialize(oConstants.nFitParameters);
  vdThetaMax.Initialize(oConstants.nFitParameters);
  vdTheta.Initialize(oConstants.nFitParameters);

  // Initialize the Width-of-the-Parameters-Vector
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  VectorFromFitParameters(oParameters, oConstants, vdTheta);
  MinVectorFromFitParameters(oParameters, oConstants, vdThetaMin);
  MaxVectorFromFitParameters(oParameters, oConstants, vdThetaMax);

  nPlotPoints = oConstants.nSmallPlotPoints;	// Number of points used in plot

  k = oConstants.k;
  l = oConstants.l;

  vdPlot.Initialize(nPlotPoints + 1);
  for(int i=0; i<=nPlotPoints; i++) {
    vdPlot[i] = double(vdThetaMin[nParameterNumber]) + i*((double(vdThetaMax[nParameterNumber]) - double(vdThetaMin[nParameterNumber]))/nPlotPoints);
  } // for i
  vdFunctionValues = vdPlot;

// From here, initialize the values of the parameters:
  VectorFromFitParameters(oParameters, oConstants, vdTheta);

  // Continue the conversation with the clients until we are done
  nIndex=0;
  while(bContinue) {
    // First check whether we need to give orders to clients
    for(int i=0; i<nProcesses-1; i++) {
      if(pbProcessHasJob[i] == false && pbProcessTerminated[i] == false && nIndex <= nPlotPoints) {
	// We need to give this process a new job. First initialize the
	// parameters we should send
	pbProcessHasJob[i] = true;
	vdTheta[nParameterNumber] = double(vdPlot[nIndex++]);
	for(int j=0; j<vdTheta.m_pnDimSize[0]; j++)
	  poMPIParameters[i].pdPar[j] = double(vdTheta[j]);
	poMPIParameters[i].dLogLik = 0;
	poMPIParameters[i].nStatus = MPI_BAN_EVAL_SUCCESS;
	MPI_Send(poMPIParameters+i, 1, *pmpiParameters, i+1, MPI_BAN_TAG_EVAL, MPI_COMM_WORLD);
	DrawProgressBar(int(100.0 * (nIndex-1) / pow(nPlotPoints+1, 1)));

	// Directly also post a non-blocking receive for this value
	pnEvalError[i] = MPI_Irecv(poMPIParameters+i, 1, *pmpiParameters, i+1, MPI_BAN_TAG_EVAL, MPI_COMM_WORLD, pmpiEvalRequest+i);
      } // if pbProcessHasJob
    } // for i

    // Now check whether some clients are done with their evaluation
    for(int i=0; i<nProcesses-1; i++) {
      pnEvalError[i] = MPI_Test(pmpiEvalRequest+i, pnEvalFlag+i, pmpiEvalStatus+i);
      // Check whether a client is done with its calculations
      if(pnEvalFlag[i]) {
	pbProcessHasJob[i] = false;
	for(int j=0; j<=nPlotPoints; j++) {
	  if(double(vdPlot[j]) == double(poMPIParameters[i].pdPar[nParameterNumber])) {
	    vdFunctionValues[j] = poMPIParameters[i].dLogLik;
	    break;
	  } // if vdPlot
	} // for j
      } // if nEvalFlag
    } // for i

    // Check whether we should send some quit messages
    if(nIndex > nPlotPoints) {
      for(int i=0; i<nProcesses-1; i++) {
	if(pbProcessTerminated[i] == false && pbProcessHasJob[i] == false) {
	  pbProcessTerminated[i] = true;
	  nMsg = MPI_BAN_MSG_QUIT;
//	  MPI_Send(&nMsg, 1, MPI_INT, i+1, MPI_BAN_TAG_MSG, MPI_COMM_WORLD);
	} // if pbProcessTerminated
      } // if pnStep[i]
    } // for i

    // Check whether we are done with the calculations
    if(nIndex > nPlotPoints) {
      bContinue = false;
      for(int i=0; i<nProcesses-1; i++) {
	if(pbProcessHasJob[i])
	  bContinue = true;
      } // for i
    } // if nPlotPoints

    // Check whether the user has requested to stop
    // 27 = <ESC>, 113 = 'q', 81 = 'Q'
//    nKey = KbHit();
    if(nKey == 27 || nKey == 113 || nKey == 81)
      bContinue = false;
  } // while bContinue

  // Send the kill to all the clients
  for(int i=0; i<nProcesses-1; i++) {
    nMsg = MPI_BAN_MSG_QUIT;
//    MPI_Send(&nMsg, 1, MPI_INT, i+1, MPI_BAN_TAG_MSG, MPI_COMM_WORLD);
  } // for i
  FinishProgressBar();

  if(nKey == 27 || nKey == 113 || nKey == 81) {
    printf("User has requested to stop the evaluation.\n");
    printf("A quit message has been sent to all the clients.\n");
  } // if nKey

  // Rescale the function for plotting
  vdTemp = Exp_1(vdFunctionValues - Min(vdFunctionValues));

  char strFileName[160];
  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcat(strFileName, "plotdata.txt");

  // Writa data to file so it can be used in with gnuplot
  WritePlot(strFileName, vdPlot, vdTemp);
//  WritePlot(strFileName, vdPlot, vdFunctionValues);

  // Clean up memory
  delete[] pbProcessHasJob;
  delete[] poMPIParameters;
  delete[] pmpiEvalRequest;
  delete[] pmpiMsgRequest;
  delete[] pmpiEvalStatus;
  delete[] pmpiMsgStatus;
  delete[] pnEvalFlag;
  delete[] pnMsgFlag;
  delete[] pnEvalError;
  delete[] pnMsgError;
} // MakePlotServer



/* This function serves as a server that commands the clients to make a simple
 * 3D plot of the likelihoodfunction
 *
 * TODO: Implement this further. Now this is a copy of MakePlotServer
 * */
void Make3DPlotServer(MPI_Datatype *pmpiParameters, SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameterNumber) {
  CVector vdPlot, vdFunctionValues, vdTemp, vdThetaMin, vdThetaMax, vdTheta;
  int nPlotPoints;
  int k, l, nParameters, nCount;
  int nProcesses, nProcessRank, nNameLength;
  char strProcessorName[MPI_MAX_PROCESSOR_NAME], strMsg[160];

  SMPIParametersType *poMPIParameters;
  MPI_Request *pmpiEvalRequest, *pmpiMsgRequest;
  MPI_Status *pmpiEvalStatus, *pmpiMsgStatus;
  int *pnEvalFlag, *pnMsgFlag, *pnEvalError, *pnMsgError;
  bool *pbProcessHasJob, *pbProcessTerminated;
  bool bContinue=true;
  int nIndex=0;
  int nMsg;
  int nKey=-1;

  // Initialize MPI environment
  MPI_Comm_size(MPI_COMM_WORLD,&nProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &nProcessRank);
  MPI_Get_processor_name(strProcessorName, &nNameLength);

  sprintf(strMsg, "Starting MPI Server on %s, %i of %i",
      strProcessorName, nProcessRank, nProcesses);
  PrintStatus(strMsg);

  // For each process we need to keep track of some data
  pbProcessHasJob = new bool[nProcesses-1];
  pbProcessTerminated = new bool[nProcesses-1];
  poMPIParameters = new SMPIParametersType[nProcesses-1];
  pmpiEvalRequest = new MPI_Request[nProcesses-1];
  pmpiMsgRequest = new MPI_Request[nProcesses-1];
  pmpiEvalStatus = new MPI_Status[nProcesses-1];
  pmpiMsgStatus = new MPI_Status[nProcesses-1];
  pnEvalFlag = new int[nProcesses-1];
  pnMsgFlag = new int[nProcesses-1];
  pnEvalError = new int[nProcesses-1];
  pnMsgError = new int[nProcesses-1];
  for(int i=0; i<nProcesses-1; i++) {
    pbProcessHasJob[i] = false;
    pbProcessTerminated[i] = false;
  } // for i


  PrintSuccess();

  InitProgressBar("Making 1D plot...");
  DrawProgressBar(0);

  nParameters = oConstants.nFitParameters;
  vdThetaMin.Initialize(oConstants.nFitParameters);
  vdThetaMax.Initialize(oConstants.nFitParameters);
  vdTheta.Initialize(oConstants.nFitParameters);

  // Initialize the Width-of-the-Parameters-Vector
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  VectorFromFitParameters(oParameters, oConstants, vdTheta);
  MinVectorFromFitParameters(oParameters, oConstants, vdThetaMin);
  MaxVectorFromFitParameters(oParameters, oConstants, vdThetaMax);

  nPlotPoints = oConstants.nSmallPlotPoints;	// Number of points used in plot

  k = oConstants.k;
  l = oConstants.l;

  vdPlot.Initialize(nPlotPoints + 1);
  for(int i=0; i<=nPlotPoints; i++) {
    vdPlot[i] = double(vdThetaMin[nParameterNumber]) + i*((double(vdThetaMax[nParameterNumber]) - double(vdThetaMin[nParameterNumber]))/nPlotPoints);
  } // for i
  vdFunctionValues = vdPlot;

// From here, initialize the values of the parameters:
  VectorFromFitParameters(oParameters, oConstants, vdTheta);

  // Continue the conversation with the clients until we are done
  nIndex=0;
  while(bContinue) {
    // First check whether we need to give orders to clients
    for(int i=0; i<nProcesses-1; i++) {
      if(pbProcessHasJob[i] == false && pbProcessTerminated[i] == false && nIndex <= nPlotPoints) {
	// We need to give this process a new job. First initialize the
	// parameters we should send
	pbProcessHasJob[i] = true;
	vdTheta[nParameterNumber] = double(vdPlot[nIndex++]);
	for(int j=0; j<vdTheta.m_pnDimSize[0]; j++)
	  poMPIParameters[i].pdPar[j] = double(vdTheta[j]);
	poMPIParameters[i].dLogLik = 0;
	poMPIParameters[i].nStatus = MPI_BAN_EVAL_SUCCESS;
	MPI_Send(poMPIParameters+i, 1, *pmpiParameters, i+1, MPI_BAN_TAG_EVAL, MPI_COMM_WORLD);
	DrawProgressBar(int(100.0 * (nIndex-1) / pow(nPlotPoints+1, 1)));

	// Directly also post a non-blocking receive for this value
	pnEvalError[i] = MPI_Irecv(poMPIParameters+i, 1, *pmpiParameters, i+1, MPI_BAN_TAG_EVAL, MPI_COMM_WORLD, pmpiEvalRequest+i);
      } // if pbProcessHasJob
    } // for i

    // Now check whether some clients are done with their evaluation
    for(int i=0; i<nProcesses-1; i++) {
      pnEvalError[i] = MPI_Test(pmpiEvalRequest+i, pnEvalFlag+i, pmpiEvalStatus+i);
      // Check whether a client is done with its calculations
      if(pnEvalFlag[i]) {
	pbProcessHasJob[i] = false;
	for(int j=0; j<=nPlotPoints; j++) {
	  if(double(vdPlot[j]) == double(poMPIParameters[i].pdPar[nParameterNumber])) {
	    vdFunctionValues[j] = poMPIParameters[i].dLogLik;
	    break;
	  } // if vdPlot
	} // for j
      } // if nEvalFlag
    } // for i

    // Check whether we should send some quit messages
    if(nIndex > nPlotPoints) {
      for(int i=0; i<nProcesses-1; i++) {
	if(pbProcessTerminated[i] == false && pbProcessHasJob[i] == false) {
	  pbProcessTerminated[i] = true;
	  nMsg = MPI_BAN_MSG_QUIT;
//	  MPI_Send(&nMsg, 1, MPI_INT, i+1, MPI_BAN_TAG_MSG, MPI_COMM_WORLD);
	} // if pbProcessTerminated
      } // if pnStep[i]
    } // for i

    // Check whether we are done with the calculations
    if(nIndex > nPlotPoints) {
      bContinue = false;
      for(int i=0; i<nProcesses-1; i++) {
	if(pbProcessHasJob[i])
	  bContinue = true;
      } // for i
    } // if nPlotPoints

    // Check whether the user has requested to stop
    // 27 = <ESC>, 113 = 'q', 81 = 'Q'
//    nKey = KbHit();
    if(nKey == 27 || nKey == 113 || nKey == 81)
      bContinue = false;
  } // while bContinue

  // Send the kill to all the clients
  for(int i=0; i<nProcesses-1; i++) {
    nMsg = MPI_BAN_MSG_QUIT;
//    MPI_Send(&nMsg, 1, MPI_INT, i+1, MPI_BAN_TAG_MSG, MPI_COMM_WORLD);
  } // for i
  FinishProgressBar();

  if(nKey == 27 || nKey == 113 || nKey == 81) {
    printf("User has requested to stop the evaluation.\n");
    printf("A quit message has been sent to all the clients.\n");
  } // if nKey

  // Rescale the function for plotting
  vdTemp = Exp_1(vdFunctionValues - Min(vdFunctionValues));

  char strFileName[160];
  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcat(strFileName, "plotdata.txt");

  // Writa data to file so it can be used in with gnuplot
  WritePlot(strFileName, vdPlot, vdTemp);
//  WritePlot(strFileName, vdPlot, vdFunctionValues);

  // Clean up memory
  delete[] pbProcessHasJob;
  delete[] poMPIParameters;
  delete[] pmpiEvalRequest;
  delete[] pmpiMsgRequest;
  delete[] pmpiEvalStatus;
  delete[] pmpiMsgStatus;
  delete[] pnEvalFlag;
  delete[] pnMsgFlag;
  delete[] pnEvalError;
  delete[] pnMsgError;
} // MakePlotServer



/* This function serves as the client in a client/server model. Useful if the
 * program is run on a cluster with MPI
 *
 * The client receives messages from the server from which it extracts
 * information on what to do. In this case the messages will consist of a set of
 * parameters; what point of the likelihood function is next to be evaluated.
 *
 * TODO: Add error-handling. Very important in multi-thread/processor
 *       applications
 * */
void MPIAnalyzeClient(MPI_Datatype *pmpiParameters, SDataType &oData, SConstantsType &oConstants) {
  int nIndex, s, p, nEvalError, nMsgError, nMsg, nEvalFlag, nMsgFlag;
  bool bContinue = true;
  bool bEvalPosted = false;
  bool bMsgPosted = false;
  SParametersType oParameters;
  SMPIParametersType oMPIParameters;
  MPI_Request mpiEvalRequest, mpiMsgRequest;
  MPI_Status mpiEvalStatus, mpiMsgStatus;

  // Initialize the parameters
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromStartParameters(oParameters, oConstants.poSources, oConstants.nSources);

  // Set the new parameters from the pdFitPar
  nIndex = 0;
  for(s=0; s<oConstants.nSources; s++) {
    for(p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	oMPIParameters.pdPar[nIndex++] = oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+p];
      } // if nFitTag
    } // for p
  } // for s

  SetBlockMatrix(oData, oConstants);
  SetGeometricPart(oData, oConstants);
  oData.mdC.Initialize(oConstants.n, oConstants.n);

  // Now initialize the Loglikelihood function
  LogLikelihoodTimesPrior(oMPIParameters.pdPar, &oData, &oParameters, &oConstants);

  // Model: if not done yet, first post non-blocking receives for both types of
  // messages. After that, parse them
  while(bContinue) {
    // Check for Msg post
    if(! bMsgPosted) {
      nMsgError = MPI_Irecv(&nMsg, 1, MPI_INT, MPI_ANY_SOURCE, MPI_BAN_TAG_MSG, MPI_COMM_WORLD, &mpiMsgRequest);
      bMsgPosted = true;
    } // if bEvalPosted

    // Check for Evaluation post
    if(! bEvalPosted) {
      nEvalError = MPI_Irecv(&oMPIParameters, 1, *pmpiParameters, MPI_ANY_SOURCE, MPI_BAN_TAG_EVAL, MPI_COMM_WORLD, &mpiEvalRequest);
      bEvalPosted = true;
    } // if bEvalPosted

    // Test for Evaluation receive
    nEvalError = MPI_Test(&mpiEvalRequest, &nEvalFlag, &mpiEvalStatus );
    if(nEvalFlag) {
      // We have received a new evaluation order
      bEvalPosted = false;

      // Now execute the order
      oMPIParameters.dLogLik = LogLikelihoodTimesPrior(oMPIParameters.pdPar);

      // Send the result back to the server
      MPI_Send(&oMPIParameters, 1, *pmpiParameters, 0, MPI_BAN_TAG_EVAL, MPI_COMM_WORLD);
    } // if nEvalFlag

    // Test for Msg receive
    nMsgError = MPI_Test(&mpiMsgRequest, &nMsgFlag, &mpiMsgStatus );
    if(nMsgFlag) {
      // We have received a message. Post a new receive
      bMsgPosted = false;
      switch(nMsg) {
	case MPI_BAN_MSG_START:
	  break;
	case MPI_BAN_MSG_QUIT:
	  bContinue = false;
	  break;
	default:
	  break;
      } // switch(nMsg)
    } // if nMsgFlag
  } // while bContinue
return ;
} // MPIAnalyzeClient


/* This function should be called after doing useful analysis work. It stops and
 * terminates all the clients.
 *
 * TODO: Add error-handling. Very important in multi-thread/processor
 *       applications
 * */
void MPIStopClients() {
  int nProcesses, nProcessRank, nNameLength, nMsg;
  char strProcessorName[MPI_MAX_PROCESSOR_NAME];

  // Initialize MPI environment
  MPI_Comm_size(MPI_COMM_WORLD,&nProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &nProcessRank);
  MPI_Get_processor_name(strProcessorName, &nNameLength);

  // Send the kill to all the clients
  for(int i=0; i<nProcesses-1; i++) {
    nMsg = MPI_BAN_MSG_QUIT;
    MPI_Send(&nMsg, 1, MPI_INT, i+1, MPI_BAN_TAG_MSG, MPI_COMM_WORLD);
  } // for i
} // MPIStopClients

/* Calculate the LogLikelihood value as in LogLikelihood above, but now it
 * is specially crafted for use in a client-server model.
 *
 * If one of the default pointers is NULL, we are expected to calculate the
 * loglikelihoodfunction using previous set configuration
 *
 * Outline:
 *
 * - x^T(invC - invC*M*inv(M^T*invC*M)*M^T*invC )x
 *
 * Last updated: 2008-10-30
 */
double LogLikelihood(double *pdFitPar, SDataType *poData, SParametersType *poParameters, SConstantsType *poConstants) {
  double dLogDetC=0, dLogDetMCM = 0, dReturnValue=0;
  static SDataType *poLocalData=NULL;
  static SParametersType *poLocalParameters=NULL;
  static SConstantsType *poLocalConstants=NULL;
  static bool bInitialized=false;
  int s, p, nIndex;
  CMatrix mdTemp;
  CVector vdNewData, vdDetResiduals;

#if 0
  time_t tTimeNow, tTimePrevious;
  double dTimeSeconds;
  tTimePrevious = clock();
#endif

  if(poConstants != NULL) {
    poLocalConstants = poConstants;
    poLocalParameters = poParameters;
    poLocalData = poData;
    bInitialized = true;
  } else if(! bInitialized) {
    return 0;
  } // poConstants

  // Set the new parameters from the pdFitPar
  nIndex = 0;
  for(s=0; s<poLocalConstants->nSources; s++) {
    for(p=0; p<poLocalConstants->poSources[s].oSourceType.nParameters; p++) {
      if(poLocalConstants->poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	poLocalParameters->pdPar[poLocalConstants->poSources[s].nFirstParIndex+p] =
	  pdFitPar[nIndex++];
      } // if nFitTag
    } // for p
  } // for s

  // First set the new coherence matrix
  SetCoherenceMatrix(*poLocalData, *poLocalParameters, *poLocalConstants);

  // Calculate the inverse and the determinant of the matrix C
  poLocalData->mdInvC = poLocalData->mdC.InverseChol(&dLogDetC);

  // Build the marginalization matrices
  mdTemp = oData.mdBlockTrans * poLocalData->mdInvC * poLocalData->mdBlock;

  // Calculate the inverse and the determinant of the matrix M^{T}C^{-1}M
  mdTemp.InvertChol(&dLogDetMCM);

  // Subtract the deterministic sources from the residuals
  vdDetResiduals = poLocalData->vdData;
  ResidualsFromDetSources(*poLocalConstants, *poLocalParameters, *poLocalData, &vdDetResiduals);
  vdNewData = poLocalData->vdData - vdDetResiduals;


  // Calculate the LogLikelihood
  dReturnValue += 0.5 * (poLocalConstants->n - poLocalConstants->nMarParameters)*log(2*M_PI);
  dReturnValue += 0.5*( vdNewData * poLocalData->mdInvC * vdNewData );
  dReturnValue -= 0.5*(vdNewData*(poLocalData->mdInvC*(poLocalData->mdBlock*(mdTemp*(oData.mdBlockTrans*(poLocalData->mdInvC*vdNewData))))));
  dReturnValue += 0.5*(dLogDetC + dLogDetMCM);

#if 0
  tTimeNow = clock();
  dTimeSeconds = double(tTimeNow - tTimePrevious)/CLOCKS_PER_SEC;
  printf("Time elapsed: %f seconds\n", dTimeSeconds);
#endif

  return dReturnValue;
} // LogLikelihood

#endif // HAVE_MPI



/* Calculate the LogLikelihood value. Outline:
**
** - x^T(invC - invC*M*inv(M^T*invC*M)*M^T*invC )x
**
** Do not use yet. No QSD removal!!!
**
*/
double LogLikelihoodSets(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dLogDetC=0, dReturnValue=0;

  SetCoherenceMatrix(oData, oParameters, oConstants);

  oData.mdInvC = oData.mdC.InverseChol(&dLogDetC);

  dReturnValue += 0.5*Trace( oData.mdInvC * oData.mdX );

  dReturnValue += 0.5*oConstants.m*dLogDetC;

  return dReturnValue;
} // LogLikelihoodSets



#if 0
/* Calculate the LogLikelihood value.
 * Assume white noise, identical for all pulsars
 *
 * Last updated: 2008-01-31
 * */
double LogLikelihood(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  static CMatrix mdD, mdInvD, mdA, mdTransA, mdAM, mdTransAM;
  static CVector vdAX;
  CLinearObject *cloTemp;
  SParametersType oParametersTemp;
  CMatrix mdTemp;
  CMatrix mdInvDCur;
  double dLogDetC=0, dLogDetMCM = 0, dReturnValue=0;

  // Check whether we need to diagonalize everything
  if(! mdD.Defined() ) {
    // Initialize the D matrices
    mdD.Initialize(oConstants.n, oConstants.n);
    for(int i=0; i<oConstants.n; i++) {
      for(int j=0; j<oConstants.n; j++) {
	mdD[i][j] = 0;
      } // for j
    } // for i
    mdInvD = mdD;I need to adjust my code to read unevenly sampled data


    // Set the parameters to only include GWB
    oParametersTemp = oParameters;
    for(int i=0; i<oConstants.k; i++) {
      oParametersTemp.mdPar[i][0] = 0;
    } // for i

    // Now calculate the GWB coherence matrix
    SetCoherenceMatrix(oData, oParametersTemp, oConstants);

    // Unitize the coherence matrix (GWB amplitude = 1)
    oData.mdC /= oParametersTemp.dAmp * oParametersTemp.dAmp;

    // Calculate the diagonalization components
    cloTemp = oData.mdC.Eigen("V");
    oData.vdLambdaC = *((CVector *)&cloTemp[0]);
    oData.mdEVC = *((CMatrix *)&cloTemp[1]);
    delete[] cloTemp;

    // Calculate the static parameters
    mdA = oData.mdEVC;
    mdTransA = oData.mdEVC[LO_TRANSPOSE];
    mdAM = mdTransA * oData.mdBlock;
    mdTransAM = mdAM[LO_TRANSPOSE];
    vdAX = mdTransA * oData.vdData;

    // Calculate the new diagonal matrices
    for(int i=0; i<oConstants.n; i++) {
      mdD[i][i] = double(oData.vdLambdaC[i]);
      mdInvD[i][i] = 1.0/double(oData.vdLambdaC[i]);
    } // for i


  } // if mdD.Defined

  mdInvDCur = mdD;
  for(int i=0; i<oConstants.n; i++) {
    mdInvDCur[i][i] = 1.0/(double(mdD[i][i])*gsl_pow_2(oParameters.dAmp)
	+ gsl_pow_2(double(oParameters.mdPar[0][0])) );
    dLogDetC += log( fabs( double(mdD[i][i])*gsl_pow_2(oParameters.dAmp)
	+ gsl_pow_2(double(oParameters.mdPar[0][0])) ) );
  } // for i

  mdTemp = mdTransAM * mdInvDCur * mdAM;
  mdTemp.InvertChol(&dLogDetMCM);

  dReturnValue += 0.5*( vdAX * mdInvDCur * vdAX );
//  dReturnValue += 0.5*(dLogDetC);

  dReturnValue -= 0.5*( vdAX * (mdInvDCur * (mdAM * (mdTemp * (mdTransAM *
	      ( mdInvDCur * vdAX ) ) ) ) ) );
  dReturnValue += 0.5*(dLogDetC + dLogDetMCM);

  return dReturnValue;
} // LogLikelihood
#endif



/* This test function calculates the power spectrum of pulsar timing residuals
 * using a least-squares spectral analysis (Lomb-Scargle). The result is written
 * to a datafile
 * */
void PowerSpectrum(SDataType &oData, SConstantsType &oConstants, int nPulsar) {
  CMatrix mdTemp;
  CVector vdReducedResiduals, vdQSD;
  CVector vdPulsarResiduals, vdPRTimes;
  CVector vdFreq, vdLSSA;
  int nIndex=0;
  FILE *pFile;

  // Calculate the best-fit quadratic spindown parameters
  mdTemp = oData.mdGenBlock[LO_TRANSPOSE] * oData.mdGenBlock;
  mdTemp.InvertChol();
  vdQSD = mdTemp * (oData.mdGenBlock[LO_TRANSPOSE] * oData.vdData);

  // Subtract linear and quadratic terms from the residuals
  vdReducedResiduals = oData.vdData - oData.mdGenBlock * vdQSD;

  for(int a=0; a<oConstants.k; a++) {
    if(nPulsar!=a) {
      nIndex += oConstants.poPulsars[a].nObservations;
    } else {
      vdPulsarResiduals.Initialize(oConstants.poPulsars[a].nObservations);
      vdPRTimes = vdPulsarResiduals;
      for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	vdPRTimes[i] = oConstants.poPulsars[a].pdTOA[i];
	vdPulsarResiduals[i] = double(vdReducedResiduals[nIndex]);
	nIndex++;
      } // for i
    } // if a
  } // for a

  PrintStatus("Calculating power spectrum...");
  //LSSA(vdPRTimes, vdPulsarResiduals, vdFreq, vdLSSA, 4);
  PrintSuccess();

  try {
    if(! (pFile = fopen("ls.txt", "w+")) ) throw 1;

    for(int i=0; i<vdPRTimes.m_pnDimSize[0]; i++) {
      fprintf(pFile, "%e  %e\n", double(vdPRTimes[i]), double(vdPulsarResiduals[i]));
    } // for i

    if(fclose(pFile) ) throw 0;

    if(! (pFile = fopen("lssa.txt", "w+")) ) throw 1;

    for(int i=0; i<vdFreq.m_pnDimSize[0]; i++) {
      fprintf(pFile, "%e  %e\n", double(vdFreq[i]), double(vdLSSA[i]));
    } // for i

    if(fclose(pFile) ) throw 0;

  } catch(int nError) {
    printf("Error number: %i\n", nError); return ;
  } // try


  return ;
} // PowerSpectrum

/*
 * This function reads a reduced basis matrix from a data file
 */
void WriteReduceMatrix(const char strFileName[], CMatrix &mdMat) {
  FILE *pFile;

  try {
    if(! (pFile = fopen(strFileName, "wb+")) ) throw 1;

    if(! fwrite(&mdMat.m_pnDimSize[0], sizeof(int), 1, pFile) ) throw 2;
    if(! fwrite(&mdMat.m_pnDimSize[1], sizeof(int), 1, pFile) ) throw 3;

    if(! fwrite(mdMat.m_pdData, sizeof(double), mdMat.m_pnDimSize[0]*mdMat.m_pnDimSize[1], pFile) ) throw 3;

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    fprintf(stderr, "Error: Could not write to file '%s'", strFileName);
    return;
  } // try
} // WriteReduceMatrix

/*
 * This function reads a reduced basis matrix from a data file
 */
void ReadReduceMatrix(const char strFileName[], CMatrix &mdMat) {
  FILE *pFile;
  int a, b;

  try {
    if(! (pFile = fopen(strFileName, "rb")) ) throw 1;

    if(! fread(&a, sizeof(int), 1, pFile) ) throw 2;
    if(! fread(&b, sizeof(int), 1, pFile) ) throw 3;

    mdMat.Initialize(a, b);

    if(! fread(mdMat.m_pdData, sizeof(double), mdMat.m_pnDimSize[0]*mdMat.m_pnDimSize[1], pFile) ) throw 3;

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    fprintf(stderr, "Error: Could not read from file '%s'", strFileName);
    return;
  } // try
} // ReadReduceMatrix


/* This function produces a file that can be used to interpolate a power-law
 * signal in a a data reduction scheme. Will be used instead of a data reduction
 * process.
 *
 * Can only work for one pulsar at a time (or should).
 */
void CreatePlInterpolation(SDataType &oData, SConstantsType &oConstants) {
#if 0
  int nSourceNumber;
  const int nPulsar=0;
  FILE *pFile;
  double dGamma;
  const char strFileName[]="interpol.int";
  CVector vdParameters, vdParametersNull;
  CMatrix mdTemp;
  vdParameters.Initialize(3);

  if(oConstants.k > 1) {
    fprintf(stderr, "ERROR: -calcint should only be used for one pulsar\n");
    return;
  } /* if oConstants.k */

  /* Look up the source number */
  nSourceNumber = -1;
  for(int s=0; s<oConstants.nSources; s++) {
    if(SourceWorksOnPulsar(oConstants, s, nPulsar) &&
	oConstants.poSources[s].eCorrelation == CID_SinglePulsar) {
      if(oConstants.poSources[s].oSourceType.eID == SID_PowerLaw) {
	nSourceNumber = s;
	break;
      } /* if SID_PowerLaw */
    } /* if SourceWorksOnPulsar */
  } /* for s */

  if(nSourceNumber == -1) {
    fprintf(stderr, "ERROR: No valid source found for interpolation\n");
    return;
  } // if nSourceNumber

  try {
    if(! (pFile = fopen(strFileName, "wb+")) ) throw 1;

    InitProgressBar("Creating elements of H^{T}C^{PN}H...");
    DrawProgressBar(0);

    /* First write the data reduction matrix */
    if(! fwrite(&oConstants.poPulsars[nPulsar].mdGReduce.m_pnDimSize[0], sizeof(int), 1, pFile) ) throw 2;
    if(! fwrite(&oConstants.poPulsars[nPulsar].mdGReduce.m_pnDimSize[1], sizeof(int), 1, pFile) ) throw 3;
    if(! fwrite(&oConstants.nInterpolationBins, sizeof(int), 1, pFile) ) throw 4;

    if(! fwrite(oConstants.poPulsars[nPulsar].mdGReduce.m_pdData, sizeof(double), oConstants.poPulsars[nPulsar].mdGReduce.m_pnDimSize[0]*oConstants.poPulsars[nPulsar].mdGReduce.m_pnDimSize[1], pFile) ) throw 5;

    for(int nIndex=0; nIndex<oConstants.nInterpolationBins; ++nIndex) {
      dGamma = -0.999 + (4.999 + 0.999) * nIndex / (oConstants.nInterpolationBins-1);

      /* Set the parameters vector */
      vdParameters[0] = 1.0;
      vdParameters[1] = dGamma;
      vdParameters[2] = oConstants.poSources[nSourceNumber].pdPar[2];

      /* Initialise the spectrum integral */
      SpectrumIntegral(oConstants.poSources[nSourceNumber].oSourceType.eID, 0, 0, vdParameters);

      /* Fill the mdC matrix */
      for(int i=0; i<oConstants.poPulsars[nPulsar].nObservations; ++i) {
	for(int j=0; j<oConstants.poPulsars[nPulsar].nObservations; ++j) {
	  oConstants.poPulsars[nPulsar].mdC[i][j] = SpectrumIntegral(
	      oConstants.poSources[nSourceNumber].oSourceType.eID,
	      oConstants.poPulsars[nPulsar].pdTOA[i],
	      oConstants.poPulsars[nPulsar].pdTOA[j],
	      vdParametersNull);
	} /* for j */
      } /* for i */

      /* Calculate mdGCG */
      mdTemp = oConstants.poPulsars[nPulsar].mdC * oConstants.poPulsars[nPulsar].mdGReduce;
      oConstants.poPulsars[nPulsar].mdGCplG = oConstants.poPulsars[nPulsar].mdGReduceT * mdTemp;

      /* Write both dGamma and mdGCplG to disk */
      if(! fwrite(&dGamma, sizeof(double), 1, pFile) ) throw 5;
      if(! fwrite(oConstants.poPulsars[nPulsar].mdGCplG.m_pdData, sizeof(double), oConstants.poPulsars[nPulsar].mdGCplG.m_pnDimSize[0]*oConstants.poPulsars[nPulsar].mdGCplG.m_pnDimSize[1], pFile) ) throw 5;

      DrawProgressBar(int(100.0*nIndex/int(oConstants.nInterpolationBins)));
    } /* for nIndex */
    fclose(pFile);
    FinishProgressBar();
  } catch(int nError) {
    FinishProgressBar();
    fprintf(stderr, "Error: Could not write to file '%s'", strFileName);
  } /* try */
#endif
} // CreatePlInterpolation

/* This function creates two files. One with elements of HCH from the
 * interpolation, one with the true elements.
 */
void CreatePlInterpolTest(SDataType &oData, SConstantsType &oConstants) {
#if 0
  int nSourceNumber;
  const int nPulsar=0;
  double dGamma;
  CVector vdParameters, vdParametersNull;
  CMatrix mdTemp;

  vdParameters.Initialize(3);
  const int nExamples=5;
  CVector vdTrue[nExamples], vdInter[nExamples], vdGamma;

  vdGamma.Initialize(oConstants.nInterpolationBins);
  for(int i=0; i<nExamples; ++i) {
    vdTrue[i].Initialize(oConstants.nInterpolationBins);
    vdInter[i].Initialize(oConstants.nInterpolationBins);
  } /* for i */

  if(oConstants.k > 1) {
    fprintf(stderr, "ERROR: -testint should only be used for one pulsar\n");
    return;
  } /* if oConstants.k */

  if(! oConstants.poPulsars[nPulsar].bInterpolatorsSet) {
    fprintf(stderr, "ERROR: interpolators not set\n");
    return;
  } /* if bInterpolatorsSet */

  /* Look up the source number */
  nSourceNumber = -1;
  for(int s=0; s<oConstants.nSources; s++) {
    if(SourceWorksOnPulsar(oConstants, s, nPulsar) &&
	oConstants.poSources[s].eCorrelation == CID_SinglePulsar) {
      if(oConstants.poSources[s].oSourceType.eID == SID_PowerLaw) {
	nSourceNumber = s;
	break;
      } /* if SID_PowerLaw */
    } /* if SourceWorksOnPulsar */
  } /* for s */

  if(nSourceNumber == -1) {
    fprintf(stderr, "ERROR: No valid source found for interpolation\n");
    return;
  } // if nSourceNumber

    for(int nIndex=0; nIndex<oConstants.nInterpolationBins; ++nIndex) {
      dGamma = -0.999 + (4.99 + 0.999) * nIndex / (oConstants.nInterpolationBins-1);
      vdGamma[nIndex] = dGamma;

    /* Set the parameters vector */
    vdParameters[0] = 1.0;
    vdParameters[1] = dGamma;
    vdParameters[2] = oConstants.poSources[nSourceNumber].pdPar[2];

    /* Initialise the spectrum integral */
    SpectrumIntegral(oConstants.poSources[nSourceNumber].oSourceType.eID, 0, 0, vdParameters);

    /* Fill the mdC matrix */
    for(int i=0; i<oConstants.poPulsars[nPulsar].nObservations; ++i) {
      for(int j=0; j<oConstants.poPulsars[nPulsar].nObservations; ++j) {
	oConstants.poPulsars[nPulsar].mdC[i][j] = SpectrumIntegral(
	    oConstants.poSources[nSourceNumber].oSourceType.eID,
	    oConstants.poPulsars[nPulsar].pdTOA[i],
	    oConstants.poPulsars[nPulsar].pdTOA[j],
	    vdParametersNull);
      } /* for j */
    } /* for i */

    /* Calculate mdGCG */
    mdTemp = oConstants.poPulsars[nPulsar].mdC * oConstants.poPulsars[nPulsar].mdGReduce;
    oConstants.poPulsars[nPulsar].mdGCplG = oConstants.poPulsars[nPulsar].mdGReduceT * mdTemp;

    /* This little nOffsetIndex doesn'tdo anything now, but it can be used to
     * check elements other than on the first line. Just let the for-loop go
     * over a few rounds.
     */
    int nOffsetIndex=0;
    for(int i=0; i<0; ++i) {
      nOffsetIndex += oConstants.poPulsars[nPulsar].mdGCplG.m_pnDimSize[1] - i;
    } // for i

    for(int i=0; i<nExamples; ++i) {
	vdTrue[i][nIndex] = double(oConstants.poPulsars[nPulsar].mdGCplG[3][i+3]);
	vdInter[i][nIndex] = gsl_spline_eval(
	    oConstants.poPulsars[nPulsar].ppIntSpline[nOffsetIndex+i],
	    dGamma,
	    oConstants.poPulsars[nPulsar].ppIntAccel[nOffsetIndex+i]);
    } /* for i */
  } /* for nIndex */

  WritePlot("truegamma.txt", vdGamma, vdTrue[0], vdTrue[1], vdTrue[2], vdTrue[3], vdTrue[4]);
  WritePlot("interpolgamma.txt", vdGamma, vdInter[0], vdInter[1], vdInter[2], vdInter[3], vdInter[4]);
#endif
} /* CreatePlInterpolTest */


/* This function produces a file that can be used to interpolate a GWB power-law
 * signal in a data reduction scheme. Will be used instead of a data reduction
 * process.
 *
 * This one does it all for the GWB
 */
void CreateGWBInterpolation(SDataType &oData, SConstantsType &oConstants) {
#if 0
  int nSourceNumber;
  FILE *pFile;
  double dGamma;
  const char strFileName[]="interpol.int";
  int nIndex1, nIndex2;
  CVector vdParameters, vdParametersNull;
  CMatrix mdTemp, mdGCG;
  vdParameters.Initialize(3);

  /* Look up the source number */
  nSourceNumber = -1;
  for(int s=0; s<oConstants.nSources; s++) {
    if(oConstants.poSources[s].eCorrelation == CID_GR && 
	oConstants.poSources[s].oSourceType.eID == SID_PowerLaw) {
      nSourceNumber = s;
      break;
    } /* if SID_PowerLaw */
  } /* for s */

  if(nSourceNumber == -1) {
    fprintf(stderr, "ERROR: No valid source found for interpolation\n");
    return;
  } // if nSourceNumber

  if(! oData.mdC.Defined()) {
    oData.mdC.Initialize(oData.mdGReduce.m_pnDimSize[0], oData.mdGReduce.m_pnDimSize[0]);
  } /* if Defined */

  try {
    if(! (pFile = fopen(strFileName, "wb+")) ) throw 1;

    InitProgressBar("Creating elements of H^{T}C^{GW}H...");
    DrawProgressBar(0);


    /* First write the data reduction matrix */
    if(! fwrite(&oData.mdGReduce.m_pnDimSize[0], sizeof(int), 1, pFile) ) throw 2;
    if(! fwrite(&oData.mdGReduce.m_pnDimSize[1], sizeof(int), 1, pFile) ) throw 2;

    if(! fwrite(&oConstants.nInterpolationBins, sizeof(int), 1, pFile) ) throw 4;

    if(! fwrite(oData.mdGReduce.m_pdData, sizeof(double), oData.mdGReduce.m_pnDimSize[0]*oData.mdGReduce.m_pnDimSize[1], pFile) ) throw 5;

    for(int nIndex=0; nIndex<oConstants.nInterpolationBins; ++nIndex) {
      dGamma = -0.999 + (4.999 + 0.999) * nIndex / (oConstants.nInterpolationBins-1);

      /* Set the parameters vector */
      vdParameters[0] = 1.0;
      vdParameters[1] = dGamma;
      vdParameters[2] = oConstants.poSources[nSourceNumber].pdPar[2];

      /* Initialise the spectrum integral */
      SpectrumIntegral(oConstants.poSources[nSourceNumber].oSourceType.eID, 0, 0, vdParameters);

      nIndex1 = nIndex2 = 0;
      for(int a=0; a<oConstants.k; a++) {
	for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	  for(int b=0; b<oConstants.k; b++) {
	    for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
		  oData.mdC[nIndex1][nIndex2] = 
		    double(oData.mdGeometric[a][b]) * SpectrumIntegral(
			    oConstants.poSources[nSourceNumber].oSourceType.eID,
			    oConstants.poPulsars[a].pdTOA[i],
			    oConstants.poPulsars[b].pdTOA[j],
			    vdParametersNull); 
	      nIndex2++;
	    } /* for j */
	  } /* for b */
	  nIndex1++;
	  nIndex2 = 0;
	} /* for i */
      } /* for a */

      /* Calculate mdGCG */
      mdTemp = oData.mdC * oData.mdGReduce;
      mdGCG = oData.mdGReduceT * mdTemp;

      /* Write both dGamma and mdGCplG to disk */
      if(! fwrite(&dGamma, sizeof(double), 1, pFile) ) throw 5;
      if(! fwrite(mdGCG.m_pdData, sizeof(double), mdGCG.m_pnDimSize[0]*mdGCG.m_pnDimSize[1], pFile) ) throw 5;

      DrawProgressBar(int(100.0*nIndex/int(oConstants.nInterpolationBins)));
    } /* for nIndex */
    fclose(pFile);
    FinishProgressBar();
  } catch(int nError) {
    FinishProgressBar();
    fprintf(stderr, "Error: Could not write to file '%s'", strFileName);
  } /* try */
#endif
} // CreateGWBInterpolation

/* This function does the matrix multiplication GCG in blocks, which should
 * speed up the process by a factor of nBlocks
 */
void MatMultGCG(CMatrix &mdG, CMatrix &mdC, CMatrix &mdResult, int nBlocks, int *pnD1, int *pnD2) {
  CMatrix *pmdG, *pmdGT, mdCs, mdGCGs, mdTemp;
  int nAddIndexFull=0, nAddIndexRed=0, nAddA=0, nAddB=0, nAddRedA=0, nAddRedB=0;

  /* We first construct the small block matrices G and GT, which are then used
   * in the double summation over all blocks of C
   */
  pmdG = new CMatrix[nBlocks];
  pmdGT = new CMatrix[nBlocks];
  for(int a=0; a<nBlocks; ++a) {
    pmdG[a].Initialize(pnD1[a], pnD2[a]);
    for(int i=0; i<pnD1[a]; ++i) {
      for(int j=0; j<pnD2[a]; ++j) {
	pmdG[a][i][j] = double(mdG[nAddIndexFull+i][nAddIndexRed+j]);
      } // for j
    } // for i
    pmdGT[a] = pmdG[a][LO_TRANSPOSE];
    nAddIndexFull += pnD1[a];
    nAddIndexRed += pnD2[a];
  } // for a

  mdResult.Initialize(mdG.m_pnDimSize[1], mdG.m_pnDimSize[1]);

  /* Now loop over all the a^2 (non-square) blocks of mdC */
  for(int a=0; a<nBlocks; ++a) {
    for(int b=0; b<nBlocks; ++b) {
      // Allocate memory for the non-square full block of C
      mdCs.Initialize(pnD1[a], pnD1[b]);
      for(int i=0; i<pnD1[a]; ++i) {
	for(int j=0; j<pnD1[b]; ++j) {
	  mdCs[i][j] = double(mdC[nAddA+i][nAddB+j]);
	} // for j
      } // for i
      mdTemp = pmdGT[a] * mdCs;
      mdGCGs = mdTemp * pmdG[b];
      for(int i=0; i<pnD2[a]; ++i) {
	for(int j=0; j<pnD2[b]; ++j) {
	  mdResult[nAddRedA+i][nAddRedB+j] = double(mdGCGs[i][j]);
	} // for j
      } // for i

      nAddB += pnD1[b];
      nAddRedB += pnD2[b];
    } // for b
    nAddA += pnD1[a];
    nAddB = 0;
    nAddRedA += pnD2[a];
    nAddRedB = 0;
  } // for a

  /* Delete the allocated memory for the block matrices */
  delete[] pmdG;
  delete[] pmdGT;
} // MatMultGCG

/* This function calculates the covariance matrices of nInterpolationBins bins
 * used to initialise the GSL cubic spline interpolation. Values are returned
 * through arguments
 *
 * nPulsar: -1 = produce covariance matrix for whole array, > 0 is produce
 * covariance matrix for one specific pulsar
 * nSource: the number of the source for which to produce these results
 * nInterpolationBins: how many bins we need for our interpolation estimation
 * pdAllData: an array, nUniques * nInterpolationBins long (return data)
 * pdGamma: an array with the gamma interpolation values (return data)
 *
 * nUniques is the l*(l+1) / 2, with l the number of observations of nPulsar
 */
void CalcSourceInterpolationBins(SDataType &oData, SConstantsType &oConstants, int nPulsar, int nSource, int nInterpolationBins, double *pdAllData, double *pdGamma) {
  CVector vdX;
  int nUniques, nRedObs, nP, nU, *pnObs, *pnRedObs;
  CMatrix mdC, mdG, mdGT, mdGCG, mdTemp;
  SParametersType oParameters;

  // First figure out the number of unique elements of the covariance matrix.
  if(nPulsar < 0) { // The whole array counts
    mdG = oData.mdGReduce;

    pnObs = new int[oConstants.k];
    pnRedObs = new int[oConstants.k];
    for(int a=0; a<oConstants.k; ++a) {
      pnObs[a] = oConstants.poPulsars[a].mdGReduce.m_pnDimSize[0];
      pnRedObs[a] = oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1];
    } // for a
  } else { // Only one pulsar
    mdG = oConstants.poPulsars[nPulsar].mdGReduce;
  } // if nPulsar
  nRedObs = mdG.m_pnDimSize[1];
  mdGT = mdG[LO_TRANSPOSE];
  nUniques = nRedObs*(nRedObs+1)/2; // How many unique elements are in the b*b matrix

//  printf("%i   %i   %i\n", nPulsar, nRedObs, nUniques);

  // Wow, it is so unclear what the difference is between these two functions :S
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromStartParameters(oParameters, oConstants.poSources, oConstants.nSources);
  nP = oConstants.poSources[nSource].nFirstParIndex;
  oParameters.pdPar[nP] = 0.0;  // Unit amplitude (exp(0) = 1.0)

  InitProgressBar("Interpolation of GCG");
  DrawProgressBar(0);
  for(int nIndex=0; nIndex<nInterpolationBins; ++nIndex) {
    oParameters.pdPar[nP+1] = pdGamma[nIndex];
    SetCoherenceMatrixOneSource(oData, oParameters, oConstants, mdC, nSource, nPulsar);

    if(nPulsar < 0) {
      MatMultGCG(mdG, mdC, mdGCG, oConstants.k, pnObs, pnRedObs);
    } else {
      // Calculate GCG
      mdTemp = mdC * mdG;
      mdGCG = mdGT * mdTemp;
    } // if nPulsar

    /* Set the uniques */
    nU = 0;
    for(int a=0; a<mdG.m_pnDimSize[1]; ++a) {
      for(int b=a; b<mdG.m_pnDimSize[1]; ++b) {
	pdAllData[nU + nUniques*nIndex] = double(mdGCG[a][b]);
	nU++;
      } // for a
    } // for a
    DrawProgressBar(int(100.0*nIndex/int(nInterpolationBins)));
  } // for nIndex
  FinishProgressBar();

  if(nPulsar < 0) {
    delete[] pnObs;
    delete[] pnRedObs;
  } // if nPulsar
} // CalcSourceInterpolationBins



/* This function sets the GSL interpolation based on a series
 * (nInterpolationBins) of values of \gamma of the covariance matrix.
 */
void SetSourceInterpolation(SDataType &oData, SConstantsType &oConstants, int nPulsar, int nSource, int nInterpolationBins, double *pdAllData, double *pdGamma) {
  // Old code: a = observations, b = number of generalised residuals
//  oConstants.poPulsars[nPulsar].mdGReduce.Initialize(a, b);
//  vdGamma.Initialize(nInterpolationBins);
  CVector vdX;
  int nUniques;
  int nRedObs;

  // First figure out the number of unique elements of the covariance matrix.
  if(nPulsar < 0) { // The whole array counts
    nRedObs = oData.mdGReduce.m_pnDimSize[1];
  } else { // Only one pulsar
    nRedObs = oConstants.poPulsars[nPulsar].mdGReduce.m_pnDimSize[1];
  } // if nPulsar
  /* How many unique elements are in the b*b matrix */
  nUniques = nRedObs*(nRedObs+1)/2;

  /* For each unique index, initialise an interpolator */
  oConstants.poSources[nSource].oIntAcc.nPulsar = nPulsar;
  oConstants.poSources[nSource].oIntAcc.nUniques = nUniques;
  oConstants.poSources[nSource].oIntAcc.ppIntAccel = new gsl_interp_accel *[nUniques];
  oConstants.poSources[nSource].oIntAcc.ppIntSpline = new gsl_spline *[nUniques];

  vdX.Initialize(nInterpolationBins);
  for(int i=0; i<nUniques; ++i) {
    /* Set up the 'dataset' */
    for(int nIndex=0; nIndex<nInterpolationBins; ++nIndex) {
      vdX[nIndex] = pdAllData[i + nUniques*nIndex];
    } /* for nIndex */

    // Allocate the accelerator and the spline object
    oConstants.poSources[nSource].oIntAcc.ppIntAccel[i] = gsl_interp_accel_alloc();
    oConstants.poSources[nSource].oIntAcc.ppIntSpline[i] = gsl_spline_alloc(
	gsl_interp_cspline, nInterpolationBins);

    /* Initialise the spline interpolator */
    gsl_spline_init(oConstants.poSources[nSource].oIntAcc.ppIntSpline[i],
	pdGamma, vdX.m_pdData, nInterpolationBins);
  } /* for i */

  oConstants.poSources[nSource].oIntAcc.bSet = true;
} // SetSourceInterpolation

/* This function prepares the acceleration of an amplitude-only source, like
 * efac/equad
 */
void SetSourceAmpAccelerator(SDataType &oData, SConstantsType &oConstants, int nPulsar, int nSource) {
  SParametersType oParameters;
  CMatrix mdTemp, mdTemp2;

  oConstants.poSources[nSource].oAmpAcc.nPulsar = nPulsar;

  // Set all the parameter boundary values (start, width, bounds, etc)
  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);

  // Set _all_ the parameter values to their start values
  ParametersFromStartParameters(oParameters, oConstants.poSources, oConstants.nSources);

  // Unit amplitude (exp(0.0) = 1.0)
  oParameters.pdPar[oConstants.poSources[nSource].nFirstParIndex] = 0.0;
  SetCoherenceMatrixOneSource(oData, oParameters, oConstants, mdTemp, nSource, oConstants.poSources[nSource].oAmpAcc.nPulsar);

  if(oConstants.poSources[nSource].oAmpAcc.nPulsar < 0) {
    mdTemp2 = oData.mdGReduceT * mdTemp;
    oConstants.poSources[nSource].oAmpAcc.mdGCaG = mdTemp2 * oData.mdGReduce;
  } else {
    mdTemp2 = oConstants.poPulsars[oConstants.poSources[nSource].oAmpAcc.nPulsar].mdGReduceT * mdTemp;
    oConstants.poSources[nSource].oAmpAcc.mdGCaG = mdTemp2 * oConstants.poPulsars[oConstants.poSources[nSource].oAmpAcc.nPulsar].mdGReduce;
  } // if nPulsar

  oConstants.poSources[nSource].oAmpAcc.bSet = true;
} // SetSourceAmpAccelerator


/* This function reads an interpolation file, and produces the gsl interpolators
 */
void ReadPlInterpolation(SDataType &oData, SConstantsType &oConstants, int nPulsar, const char strFileName[]) {
#if 0
  int a, b;
  FILE *pFile;
  int nInterpolationBins, nUniques, nUniqueIndex;
  double dGamma;
  CVector vdGamma, vdX;
  CMatrix mdGCplG, mdTot;
  gsl_interp_accel *pAc;
  gsl_spline *pSp;

  try {
    if(! (pFile = fopen(strFileName, "rb")) ) throw 1;

    if(! fread(&a, sizeof(int), 1, pFile) ) throw 2;
    if(! fread(&b, sizeof(int), 1, pFile) ) throw 3;

    if(! fread(&nInterpolationBins, sizeof(int), 1, pFile) ) throw 4;

    /* How many unique elements are in the b*b matrix */
    nUniques = (b * (b+1)) / 2;

    /* Allocate memory for the bins */
    oConstants.poPulsars[nPulsar].mdGReduce.Initialize(a, b);
    vdGamma.Initialize(nInterpolationBins);
    mdGCplG.Initialize(b, b);
    mdTot.Initialize(nUniques, nInterpolationBins);

    /* Read the data compression matrix */
    if(! fread(oConstants.poPulsars[nPulsar].mdGReduce.m_pdData, sizeof(double), a*b, pFile) ) throw 5;

    for(int nIndex=0; nIndex<nInterpolationBins; ++nIndex) {

      if(! fread(&dGamma, sizeof(double), 1, pFile) ) throw 5;
      if(! fread(mdGCplG.m_pdData, sizeof(double), b*b, pFile) ) throw 5;
      vdGamma[nIndex] = dGamma;

      /* Fill the collection matrix */
      nUniqueIndex = 0;
      for(int i=0; i<b; ++i) {
	for(int j=i; j<b; ++j) {
	  mdTot[nUniqueIndex][nIndex] = double(mdGCplG[i][j]);
	  nUniqueIndex++;
	} /* for j */
      } /* for i */
    } /* for nIndex */
    fclose(pFile);

    /* For each unique index, initialise an interpolator */
    oConstants.poPulsars[nPulsar].ppIntAccel = new gsl_interp_accel *[nUniques];
    oConstants.poPulsars[nPulsar].ppIntSpline = new gsl_spline *[nUniques];
    vdX.Initialize(nInterpolationBins);
    for(int i=0; i<nUniques; ++i) {
      /* Set up the 'dataset' */
      for(int nIndex=0; nIndex<nInterpolationBins; ++nIndex) {
	vdX[nIndex] = double(mdTot[i][nIndex]);
      } /* for nIndex */

      // Allocate the accelerator and the spline object
      oConstants.poPulsars[nPulsar].ppIntAccel[i] = gsl_interp_accel_alloc();
      oConstants.poPulsars[nPulsar].ppIntSpline[i] = gsl_spline_alloc(
	  gsl_interp_cspline, nInterpolationBins);
      pAc = oConstants.poPulsars[nPulsar].ppIntAccel[i];
      pSp = oConstants.poPulsars[nPulsar].ppIntSpline[i];

      /* Initialise the spline interpolator */
      gsl_spline_init(pSp, vdGamma.m_pdData, vdX.m_pdData, nInterpolationBins);
    } /* for i */

    oConstants.poPulsars[nPulsar].nUniques = nUniques;
    oConstants.poPulsars[nPulsar].bInterpolatorsSet = true;
  } catch(int nError) {
    fprintf(stderr, "Error: Could not read from file '%s'", strFileName);
  } /* try */
#endif
} // ReadPlInterpolation

/* This function reads an interpolation file, and produces the gsl interpolators
 */
void ReadGWBInterpolation(SDataType &oData, SConstantsType &oConstants, const char strFileName[]) {
#if 0
  int a, b;
  FILE *pFile;
  int nInterpolationBins, nUniques, nUniqueIndex;
  double dGamma;
  CVector vdGamma, vdX;
  CMatrix mdGCplG, mdTot;
  gsl_interp_accel *pAc;
  gsl_spline *pSp;

  try {
    if(! (pFile = fopen(strFileName, "rb")) ) throw 1;

    if(! fread(&a, sizeof(int), 1, pFile) ) throw 2;
    if(! fread(&b, sizeof(int), 1, pFile) ) throw 3;

    if(! fread(&nInterpolationBins, sizeof(int), 1, pFile) ) throw 4;

    /* How many unique elements are in the b*b matrix */
    nUniques = (b * (b+1)) / 2;

    /* Allocate memory for the bins */
    oData.mdGReduce.Initialize(a, b);
    vdGamma.Initialize(nInterpolationBins);
    mdGCplG.Initialize(b, b);
    mdTot.Initialize(nUniques, nInterpolationBins);

    /* Read the data compression matrix */
    if(! fread(oData.mdGReduce.m_pdData, sizeof(double), a*b, pFile) ) throw 5;

    for(int nIndex=0; nIndex<nInterpolationBins; ++nIndex) {

      if(! fread(&dGamma, sizeof(double), 1, pFile) ) throw 5;
      if(! fread(mdGCplG.m_pdData, sizeof(double), b*b, pFile) ) throw 5;
      vdGamma[nIndex] = dGamma;

      /* Fill the collection matrix */
      nUniqueIndex = 0;
      for(int i=0; i<b; ++i) {
	for(int j=i; j<b; ++j) {
	  mdTot[nUniqueIndex][nIndex] = double(mdGCplG[i][j]);
	  nUniqueIndex++;
	} /* for j */
      } /* for i */
    } /* for nIndex */
    fclose(pFile);

    /* For each unique index, initialise an interpolator */
    oData.ppGWBIntAccel = new gsl_interp_accel *[nUniques];
    oData.ppGWBIntSpline = new gsl_spline *[nUniques];
    vdX.Initialize(nInterpolationBins);
    for(int i=0; i<nUniques; ++i) {
      /* Set up the 'dataset' */
      for(int nIndex=0; nIndex<nInterpolationBins; ++nIndex) {
	vdX[nIndex] = double(mdTot[i][nIndex]);
      } /* for nIndex */

      // Allocate the accelerator and the spline object
      oData.ppGWBIntAccel[i] = gsl_interp_accel_alloc();
      oData.ppGWBIntSpline[i] = gsl_spline_alloc(
	  gsl_interp_cspline, nInterpolationBins);
      pAc = oData.ppGWBIntAccel[i];
      pSp = oData.ppGWBIntSpline[i];

      /* Initialise the spline interpolator */
      gsl_spline_init(pSp, vdGamma.m_pdData, vdX.m_pdData, nInterpolationBins);
    } /* for i */

    oData.nGWBUniques = nUniques;
    oData.bGWBInterpolatorsSet = true;
  } catch(int nError) {
    fprintf(stderr, "Error: Could not read from file '%s'", strFileName);
  } /* try */
#endif
} // ReadGWBInterpolation

/* This function reads an interpolation file, and produces the gsl interpolators
 */
void ReadClockInterpolation(SDataType &oData, SConstantsType &oConstants, const char strFileName[]) {
#if 0
  int a, b;
  FILE *pFile;
  int nInterpolationBins, nUniques, nUniqueIndex;
  double dGamma;
  CVector vdGamma, vdX;
  CMatrix mdGCplG, mdTot;
  gsl_interp_accel *pAc;
  gsl_spline *pSp;

  try {
    if(! (pFile = fopen(strFileName, "rb")) ) throw 1;

    if(! fread(&a, sizeof(int), 1, pFile) ) throw 2;
    if(! fread(&b, sizeof(int), 1, pFile) ) throw 3;

    if(! fread(&nInterpolationBins, sizeof(int), 1, pFile) ) throw 4;

    /* How many unique elements are in the b*b matrix */
    nUniques = (b * (b+1)) / 2;

    /* Allocate memory for the bins */
    oData.mdGReduce.Initialize(a, b);
    vdGamma.Initialize(nInterpolationBins);
    mdGCplG.Initialize(b, b);
    mdTot.Initialize(nUniques, nInterpolationBins);

    /* Read the data compression matrix */
    if(! fread(oData.mdGReduce.m_pdData, sizeof(double), a*b, pFile) ) throw 5;

    for(int nIndex=0; nIndex<nInterpolationBins; ++nIndex) {

      if(! fread(&dGamma, sizeof(double), 1, pFile) ) throw 5;
      if(! fread(mdGCplG.m_pdData, sizeof(double), b*b, pFile) ) throw 5;
      vdGamma[nIndex] = dGamma;

      /* Fill the collection matrix */
      nUniqueIndex = 0;
      for(int i=0; i<b; ++i) {
	for(int j=i; j<b; ++j) {
	  mdTot[nUniqueIndex][nIndex] = double(mdGCplG[i][j]);
	  nUniqueIndex++;
	} /* for j */
      } /* for i */
    } /* for nIndex */
    fclose(pFile);

    /* For each unique index, initialise an interpolator */
    oData.ppClockIntAccel = new gsl_interp_accel *[nUniques];
    oData.ppClockIntSpline = new gsl_spline *[nUniques];
    vdX.Initialize(nInterpolationBins);
    for(int i=0; i<nUniques; ++i) {
      /* Set up the 'dataset' */
      for(int nIndex=0; nIndex<nInterpolationBins; ++nIndex) {
	vdX[nIndex] = double(mdTot[i][nIndex]);
      } /* for nIndex */

      // Allocate the accelerator and the spline object
      oData.ppClockIntAccel[i] = gsl_interp_accel_alloc();
      oData.ppClockIntSpline[i] = gsl_spline_alloc(
	  gsl_interp_cspline, nInterpolationBins);
      pAc = oData.ppClockIntAccel[i];
      pSp = oData.ppClockIntSpline[i];

      /* Initialise the spline interpolator */
      gsl_spline_init(pSp, vdGamma.m_pdData, vdX.m_pdData, nInterpolationBins);
    } /* for i */

    oData.nClockUniques = nUniques;
    oData.bClockInterpolatorsSet = true;
  } catch(int nError) {
    fprintf(stderr, "Error: Could not read from file '%s'", strFileName);
  } /* try */
#endif
} // ReadClockInterpolation

/* This matrix interpolates the GCplG matrix.
 */
void CalcPlIntMatrix(gsl_interp_accel **ppAc, gsl_spline **ppSp, double dGamma, CMatrix &mdGCplG) {
  int nIndex=0;
  for(int i=0; i<mdGCplG.m_pnDimSize[0]; ++i) {
    for(int j=i; j<mdGCplG.m_pnDimSize[1]; ++j) {
      mdGCplG[i][j] = gsl_spline_eval(
	  ppSp[nIndex],
	  dGamma,
	  ppAc[nIndex]);
      if(i != j) {
	mdGCplG[j][i] = double(mdGCplG[i][j]);
      } /* if i */
      nIndex++;
    } /* for j */
  } /* for i */
} // CalcPlIntMatrix

/* For now use the start parameters as the data reduction values
 *
 * nReducePoints is the number of datapoints we reduce our dataset to. This
 * overrides the fidelity setting
 *
 * The result is returned in mdG
 */
void ReduceData(SDataType &oData, SConstantsType &oConstants, CMatrix &mdResult, double dFidelity, int nReducePoints) {
  SParametersType oParameters;
  CMatrix mdR, mdRInv, mdU, mdV, mdG, mdGTrans, mdS, mdA, mdTemp,
	  mdEV, mdLo, mdUp, mdSqrt, mdH;
  CVector vdS, vdLambda;
  CMatrix mdBlockTemp;
  CLinearObject *cloTemp;
  double dFrac, dTotal, dSum;
  bool bUseCholesky = false;
  int nIndex=0;

//  dFrac = oConstants.dReducedFidelity;
  dFrac = dFidelity;

  StartParametersFromSources(oParameters, oConstants.poSources, oConstants.nSources);
  ParametersFromStartParameters(oParameters, oConstants.poSources, oConstants.nSources);

  /* Check the singular value decomposition of the design matrix
   */
  mdBlockTemp = oData.mdBlock;
  mdBlockTemp.SVD(mdU, mdV, vdS);

  /* Produce the timing model data reduction matrix
   */
  mdG.Initialize(oConstants.n, oConstants.n-oConstants.nMarParameters);
  for(int i=0; i<oConstants.n; ++i) {
    for(int j=oConstants.nMarParameters; j<oConstants.n; ++j) {
      mdG[i][j-oConstants.nMarParameters] = double(mdU[i][j]);
    } // for j
  } // for i
  mdGTrans = mdG[LO_TRANSPOSE];

  /* Produce the noise matrix
   */
  SetCoherenceMatrixNonReduceSources(oData, oParameters, oConstants);
  mdR = oData.mdC;
  mdTemp = mdR * mdG;
  mdR = mdGTrans * mdTemp;

//  PrintMatrix(oData.mdC);
//  PrintMatrix(mdR);

  /* Produce the signal matrix
   */
  SetCoherenceMatrixReduceSources(oData, oParameters, oConstants);
  mdTemp = oData.mdC * mdG;
  mdS = mdGTrans * mdTemp;

  /*
   * Use an SVD per default. A Cholesky decomposition can also be used, but
   * seems less numerically stable...
   */
  if(bUseCholesky) {
    try {
      mdRInv = mdR.InverseChol();

      mdLo = mdRInv.Cholesky();
      mdUp = mdLo[LO_TRANSPOSE];
      mdTemp = mdLo * mdS;
      mdA = mdTemp * mdUp;

      mdBlockTemp = mdA;
      mdBlockTemp.SVD(mdU, mdV, vdS);

      mdSqrt = mdLo;
    } catch(ELinearError err) {
      bUseCholesky = false;
    } // try
  } // if bUseCholesky

  if(! bUseCholesky) {
    /* Cholesky didn't work on the noise matrix. Use the SVD
     */
    mdTemp = mdR;
    mdTemp.SVD(mdU, mdV, vdS);

    mdSqrt.Initialize(vdS.m_pnDimSize[0], vdS.m_pnDimSize[0]);
    for(int i=0; i<vdS.m_pnDimSize[0]; ++i) {
      for(int j=0; j<vdS.m_pnDimSize[0]; ++j) {
	if(i == j && double(vdS[i]) > 0) {
	  mdSqrt[i][j] = 1.0 / sqrt(double(vdS[i]));
	} else {
	  mdSqrt[i][j] = 0;
	} /* if i == j */
      } /* for j */
    } /* for i */

    mdA = mdU[LO_TRANSPOSE];
    mdTemp = mdSqrt * mdA;
    mdSqrt = mdV * mdTemp;

    mdTemp = mdSqrt * mdS;
    mdA = mdTemp * mdSqrt;

    mdBlockTemp = mdA;
    mdBlockTemp.SVD(mdU, mdV, vdS);

//    fprintf(stderr, "\nmdA.SVD:\n");
//    for(int i=0; i<vdS.m_pnDimSize[0]; ++i) {
//      fprintf(stderr, "%e  ", double(vdS[i])*double(vdS[i]) / ((1 + double(vdS[i]))*(1+double(vdS[i]))));
//    } /* for i */
//    fprintf(stderr, "\n");
//    PrintVector(vdS);
 } // if bUseCholesky

  /* Select the upper 99% of the signal
   */
  dTotal = 0;
  for(int i=0; i<mdA.m_pnDimSize[0]; ++i) {
//    dTotal += double(mdA[i][i]);
    dTotal += double(vdS[i]) * double(vdS[i]) / ((1 + double(vdS[i]))*(1 + double(vdS[i])));
  } // for i

  dSum = 0; nIndex = vdS.m_pnDimSize[0];
  for(int i=0; i<vdS.m_pnDimSize[0]; ++i) {
//    dSum += double(vdS[i]);
    dSum += double(vdS[i]) * double(vdS[i]) / ((1 + double(vdS[i]))*(1 + double(vdS[i])));

    nIndex = i;
    if(nReducePoints > 0 && nReducePoints == i+1) {
      break;
    } else if(nReducePoints <= 0 && dSum >= dFrac*dTotal) {
      break;
    } // if nReducePoints
  } // for i

#if 0
  /* Save the fidelity */
  {
    CVector vdCum, vdLC, vdX, vdLX, vdLSA;
    vdCum = vdS;
    vdLC = vdS;
    vdX = vdS;
    vdLX = vdS;
    vdLSA = vdS;
    vdCum[0] =  double(vdS[0]) * double(vdS[0]) / ((1 + double(vdS[0]))*(1 + double(vdS[0])));
    vdLSA[0] = double(vdS[0]) * double(vdS[0]);
    for(int i=1; i<vdS.m_pnDimSize[0]; ++i) {
      vdCum[i] = double(vdCum[i-1]) + double(vdS[i]) * double(vdS[i]) / ((1 + double(vdS[i]))*(1 + double(vdS[i])));
      vdLSA[i] = double(vdLSA[i-1]) +  double(vdS[i]) * double(vdS[i]);
    } // for i
    vdCum /= dTotal;
    vdLSA /= double(vdLSA[vdLSA.m_pnDimSize[0]-1]);
    for(int i=0; i<vdS.m_pnDimSize[0]; ++i) {
      vdLX[i] = log(i+1);
      vdX[i] = i+1;
      vdLC[i] = log(double(vdCum[i]));
    } // for i
    WritePlot("fidelity.txt", vdX, vdCum, vdLSA, vdLX, vdLC);
  }
#endif

  fprintf(stderr, "Compressed basis: %i->%i with %4.1f%% fidelity\n",
      vdS.m_pnDimSize[0], nIndex+1, 100.0*dSum/dTotal);
//  fprintf(stderr, "dTotal = %e   dSum = %e\n", dTotal, dSum);

  /* Calculate the reduced basis in the timing model reduced basis
   */
  mdTemp = mdSqrt * mdU;
  mdH.Initialize(mdU.m_pnDimSize[0], nIndex+1);
  for(int i=0; i<mdH.m_pnDimSize[0]; ++i) {
    for(int j=0; j<nIndex+1; ++j) {
      mdH[i][j] = double(mdTemp[i][j]);
    } // for j
  } // for i

  /* Form the total reduced basis mdG
   */
  oData.mdGReduce = mdG * mdH;

  /* Return the result
   */
  mdResult = oData.mdGReduce;
} // ReduceData


/* This function returns the value of log(Prior)
 *
 * Input
 *   d: The value of the parameter (the one that integrates to 1)
 *   dMean:      The center of the Lorentzian function
 *   dWidth:     The width of the Lorentzian function (how well determined the
 *               variable is by construction
 *
 * Output:
 *   return value: The value of the function
 * */
double LogPrior(SParametersType &oParameters, SConstantsType &oConstants) {
  double dReturnValue=0;
  int nTemp, nIndexFit;

  // Add a prior to the loglikelihood value
  nIndexFit = 0;
  for(int s=0; s<oConstants.nSources; s++) {
    if(SourceWorks(oConstants, s)) {
      for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
	nTemp = oConstants.poSources[s].nFirstParIndex + p;

	if(oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p)) {
	  // We are analytically marginalising over this parameter
	  // Figure out a way to incorporate the correct prior here...
	} else if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	  // We are fitting for this parameter (in MCMC)
	  if(oParameters.pdParWidthPrior[nTemp] > 0) {
	    // We are using a Lorentzian prior
	    dReturnValue += log(LorentzianPrior(oParameters.pdPar[nTemp],
		  oParameters.pdParMeanPrior[nTemp], oParameters.pdParWidthPrior[nTemp]));
	  } else if(oParameters.pdParWidthPrior[nTemp] < 0) {
	    dReturnValue -= log( log( oParameters.pdParMaxBound[nTemp] ) - log(oParameters.pdParMinBound[nTemp])) + log(oParameters.pdPar[nTemp]);
	  } else {
	    // Use the min and max bound as the prior volume.
	    dReturnValue -= log(oParameters.pdParMaxBound[nTemp] - oParameters.pdParMinBound[nTemp]);
	  } // if pdParWidthPrior
	  nIndexFit++;
	} // if nMarTag
      } // for p
    } // if SourceWorks
  } // for s
//  return 0;  // return 0 FOR MLDR check. Remove or THIS MAY CAUSE YOU PROBLEMS LATER
  return dReturnValue;
} // LogPrior


/* This function returns the value of a Lorentzian function. Useful when adding
 * a prior while doing the mcmc.
 *
 * Input
 *   dParameter: The value of the parameter (the one that integrates to 1)
 *   dMean:      The center of the Lorentzian function
 *   dWidth:     The width of the Lorentzian function (how well determined the
 *               variable is by construction
 *
 * Output:
 *   return value: The value of the function
 * */
double LorentzianPrior(double dParameter, double dMean, double dWidth) {
  return dWidth / (M_PI * (dWidth*dWidth + (dParameter-dMean)*(dParameter-dMean)));
} // LorentzianPrior

#ifdef HAVE_MPI
double LogLikelihoodTimesPrior(double *pdFitPar, SDataType *poData, SParametersType *poParameters, SConstantsType *poConstants) {
  double dReturnValue;
  int nTemp;

  static bool bInitialized=false;
  static SParametersType *poLocalParameters=NULL;
  static SConstantsType *poLocalConstants=NULL;
  static SDataType *poLocalData=NULL;

  if(poConstants != NULL) {
    poLocalConstants = poConstants;
    poLocalParameters = poParameters;
    poLocalData = poData;
    bInitialized = true;
  } else if(! bInitialized) {
    return 0;
  } // poConstants

  dReturnValue = LogLikelihood(pdFitPar, poLocalData, poLocalParameters, poLocalConstants);
  dReturnValue -= LogPrior(pdFitPar, poLocalParameters, poLocalConstants);

  return dReturnValue;
} // LogLikelihoodTimesPrior
#endif // HAVE_MPI

double LogLikelihoodTimesPrior(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dReturnValue;
  dReturnValue = LogLikelihood(oData, oParameters, oConstants);
  dReturnValue -= LogPrior(oParameters, oConstants);

  return dReturnValue;
} // LogLikelihoodTimesPrior
