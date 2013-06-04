/* banginterface.cpp -- interface to the bang functions in C
 *
 * vim: tabstop=2:softtabstop=2:shiftwidth=2:expandtab
 *
 * Rutger van Haasteren 09 December 2012 vhaasteren@gmail.com
 *
 * Copyright (C) 2011-2012 Rutger van Haasteren.
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

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.  */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_cdf.h>

#include "bangui.h"
#include "linal.h"
#include "linalfunc.h"
#include "corefunctions.h"
#include "banfunc.h"
#include "correlation.h"
#include "analyze.h"

// This is taken from tempo2.h. Use here because you are too lazy to change shit
#define MAX_FILELEN 500

BanWindow *g_pMainWindow=NULL;
int g_nDimSize=0;
int *g_pnParameterIndex=NULL;

/* The usage of this interface is as follows. Call all in this order or things
 * go wrong.
 *
 * required
 * InitWrapper
 *     - nPsr: number of pulsars
 *     - pnPsrObs: array with number of observations per pulsar
 *     - pdPhi: array with right ascension of each pulsar
 *     - pdTheta: array with declination (pi/2 - decj) of each pulsar
 *     - pdTOAs: array with TOAs of all pulsars combined
 *     - pdPreRes: array with prefit residuals of all pulsars combined
 *     - pdErr: array with error bars of all pulsars combined
 *     - pdDesignMatrix: array with the full design matrix
 *     - pnTMPars: array with the number of timing model parameters per pulsar
 *     - pdTMPVal: array with the prefit value of the timing model parameters
 *     - pstrPulsarNames: array with strings, with the names of all pulsars
 *     - pstrFlags: array with strings, with the '-instr' flag names for all
 *       observations for all pulsars
 *     - pstrTMPDescription: array with strings, with descriptions of all timing
 *       model parameters for all pulsars
 *
 * required
 * AddLinearSource
 *     - nPsr: number of the pulsar for which to add the linear source
 *     - nTMP: number of timing model parameters for this source
 * Note: This is sort of incompatible with the HDF5 design tree of information,
 *       but it is just how the code works now. It is not completely obvious to
 *       what extend this information is useful, but just call this function and
 *       it will all be fine.
 *
 * optional
 * SetPulsarCompression
 *     - nPulsar: number of the pulsar for which to set the compression matrix
 *     - nRedObs: number of generalised residuals for this pulsar
 *     - pdCompressionBasis: the compression matrix/basis that defines the
 *       compression
 *
 * Previous two done k times (k number of pulsars)
 *
 * required
 * ProcessNewBasis (no arguments)
 *
 *
 * required
 * AddEfacSource (Add error bars to the covariance matrix, including efac)
 *     - strFlag: the flag value this source works on
 *     - dMin: minimum of the domain of this parameter
 *     - dMax: maximum of the domain of this parameter
 *     - dStart: starting point in parameter space (or fixed value)
 *     - dStepSize: MCMC proposal width, and scale parameter
 *     - nVary: true if included in the (MCMC) search, false otherwise
 *
 * SetSourceAmpAccelerator (required for efac, and equad)
 *     - nPulsar: pulsar for which this source is valid (-1 for all, like GWB)
 *     - nSource: internal number of the source for which to set the interpolator
 *
 * AddEquadSource (white noise)
 *     - strFlag: the flag value this source works on
 *     - dMin: minimum of the domain of this parameter
 *     - dMax: maximum of the domain of this parameter
 *     - dStart: starting point in parameter space (or fixed value)
 *     - dStepSize: MCMC proposal width, and scale parameter
 *     - nVary: true if included in the (MCMC) search, false otherwise
 * 
 * AddRedNoiseSource (power-law red noise)
 *     - strFlag: the flag value this source works on
 *     - nPsr: number of the pulsar this source works on (used if strFlag=="")
 *     - pdMin: array of the minimum (3 elements)
 *     - pdMax: array of the maximum ''
 *     - pdStart: array of the start ''
 *     - pdStepSize: proposal width ''
 *     - pnVary: array ''
 * 
 * AddDMVSource (power-law dispersion measure variations)
 *     - strFlag: the flag value this source works on
 *     - nPsr: number of the pulsar this source works on (used if strFlag=="")
 *     - pdMin: array of the minimum (3 elements)
 *     - pdMax: array of the maximum ''
 *     - pdStart: array of the start ''
 *     - pdStepSize: proposal width ''
 *     - pnVary: array ''
 *
 * SetSourceInterpolator (for power-law signals, also for GWB)
 *     - nPulsar: pulsar for which this source is valid (-1 for all, like for GWB)
 *     - nSource: internal number of the source for which to set the interpolator
 *     - nInterpolationBins: number of bins for the cubic spline interpolators
 *     - pdData: the raw interpolation data, read from the hdf5 file
 *     - pdGamma: the values of the spectral index in the interpolation bins
 *
 * AddGWBSource (GWB, power-law, works on all observations)
 *     - pdMin: array of the minimum (3 elements)
 *     - pdMax: array of the maximum ''
 *     - pdStart: array of the start ''
 *     - pdStepSize: proposal width ''
 *     - pnVary: array ''
 *
 *
 * required
 * CompleteInitWrapper (finish setting up the wrapper)
 *
 *
 * ** Now the wrapper is all set up. Obtain info through: **
 * GetParameterDimensions
 *    returns the number of dimensions in the parameter space
 *
 * GetBoundaryConditions
 *    - pdMin: returns the minimum values of the varying parameters
 *    - pdMax: returns the maximum values of the varying parameters
 *    - pdStart: returns the start value of the varying parameters
 *    - pdWidth: returns the stepsize of the varying parameters
 *
 * WrapLL (returns the value of the log-likelihood)
 *    - pdPar: array with the parameters
 *
 * WrapHCLL
 *    - pdPar: same as above, but now with these parameters on a hypercube
 *
 * PrintMdoel (prints the raw model representation of the wrapper)
 *
 *
 * ** That's all folks **
 *
 * For niceness, clear memory with:
 *
 * FinishWrapper
 */


void SetPulsarData(int nPsr, int *pnPsrObs, double *pdPhi, double *pdTheta, double *pdTOAs, double *pdPreRes, double *pdErr, double *pdFreq, double *pdDesignMatrix, int *pnTMPars, double *pdTMPVal, char **pstrPulsarNames, char **pstrFlags, char **pstrTMPDescription);
extern "C" void ProcessNewBasis();
extern "C" void FinishWrapper();
extern "C" void InitWrapper(int nPsr, int *pnPsrObs, double *pdPhi, double *pdTheta, double *pdTOAs, double *pdPreRes, double *pdErr, double *pdFreq, double *pdDesignMatrix, int *pnTMPars, double *pdTMPVal, char **pstrPulsarNames, char **pstrFlags, char **pstrTMPDescription);

extern "C" void SetPulsarCompression(int nPulsar, int nRedObs, double *pdCompressionBasis);
extern "C" int CompressionMatrix(double *pdCompressedBasis, double dFidelity=0.99, int nReducePoints=-1);
extern "C" void CovarianceMatrix(double *pdCovariance, double *pdPar);
extern "C" void CalcInterpolData(int nPulsar, int nSource, int nInterpolationBins, double *pdData, double *pdGamma);
extern "C" void SetSourceInterpolator(int nPulsar, int nSource, int nInterpolationBins, double *pdData, double *pdGamma);
extern "C" void SetSourceAmpAccelerator(int nPulsar, int nSource);

extern "C" void AddEfacSource(char strFlag[], int nPsr, double dMin, double dMax, double dStart, double dStepSize, int nVary, bool bRed);
extern "C" void AddEquadSource(char strFlag[], int nPsr, double dMin, double dMax, double dStart, double dStepSize, int nVary, bool bRed);
extern "C" void AddRedNoiseSource(char strFlag[], int nPsr, double *pdMin, double *pdMax, double *pdStart, double *pdStepSize, long *pnVary, bool bRed);
extern "C" void AddDMVSource(char strFlag[], int nPsr, double *pdMin, double *pdMax, double *pdStart, double *pdStepSize, long *pnVary, bool bRed);
extern "C" void AddGWBSource(double *pdMin, double *pdMax, double *pdStart, double *pdStepSize, long *pnVary, bool bRed);
extern "C" void AddDipoleGWBSource(int nMode, double *pdMin, double *pdMax, double *pdStart, double *pdStepSize, long *pnVary, bool bRed);
extern "C" void CompleteInitWrapper();

extern "C" void PrintModel();

extern "C" void SetUseReducedBasis(bool bUseReducedBasis=true);

extern "C" int GetParameterDimensions();
extern "C" void GetBoundaryConditions(double *pdMin, double *pdMax, double *pdStart, double *pdWidth);
extern "C" void GetAmplitudeIndices(double *pdIndices);
extern "C" double WrapLL(double *pdPar);
extern "C" double WrapHCLL(double *pdPar);



/* Initialise the standard data structures.
 *
 * Units:
 * pdPhi [rad]
 * pdTheta [rad] (Note: this is pi/2 - decj
 * pdTOAs [sec]
 * pdPreRes [sec]
 * pdErr [sec]
 * pdFreq [MHz]
 * pdDesignMatrix [varies]
 * pdTMPVal [varies]
 *
 */
extern "C" void InitWrapper(int nPsr, int *pnPsrObs, double *pdPhi, double *pdTheta, double *pdTOAs, double *pdPreRes, double *pdErr, double *pdFreq, double *pdDesignMatrix, int *pnTMPars, double *pdTMPVal, char **pstrPulsarNames, char **pstrFlags, char **pstrTMPDescription) {
  // Clean up possible previous data
  FinishWrapper(); 

  // Initialise the global variable containing the analysis
  // This now no longer reads "ReadGlobalConstants"
  g_pMainWindow = new BanWindow(0, NULL, true, NULL);

  // Initialise the basic data structures, except the actual model (sources)
  SetPulsarData(nPsr, pnPsrObs, pdPhi, pdTheta, pdTOAs, pdPreRes, pdErr, pdFreq, pdDesignMatrix, pnTMPars, pdTMPVal, pstrPulsarNames, pstrFlags, pstrTMPDescription);
} // InitWrapper


/* Set all the standard information that is obtained from tempo2
 *
 * Units:
 * pdPhi [rad]
 * pdTheta [rad] (Note: this is pi/2 - decj
 * pdTOAs [sec]
 * pdPreRes [sec]
 * pdErr [sec]
 * pdFreq [MHz]
 * pdDesignMatrix [varies]
 * pdTMPVal [varies]
 *
 */
void SetPulsarData(int nPsr, int *pnPsrObs, double *pdPhi, double *pdTheta, double *pdTOAs, double *pdPreRes, double *pdErr, double *pdFreq, double *pdDesignMatrix, int *pnTMPars, double *pdTMPVal, char **pstrPulsarNames, char **pstrFlags, char **pstrTMPDescription) {
  int nIndex=0, nIndexTM=0;
  g_pMainWindow->m_oConstants.k = 0;
  g_pMainWindow->m_oConstants.n = 0;
  g_pMainWindow->m_oConstants.l = 0;
  g_pMainWindow->m_oConstants.m = 0;
  g_pMainWindow->m_oConstants.nMarParameters = 0;
  g_pMainWindow->m_oConstants.nParameters = 0;
  g_pMainWindow->m_oConstants.nFitParameters = 0;

  // Reserve memory for the pulsars
  // TODO: very sloppy to do this here. Is owned by g_pMainWindow, so should be
  // done by that class. But for now this works
  g_pMainWindow->AllocatePulsars(nPsr, pnPsrObs);

  // Set the total number of observations and timing model parameters
  for(int a=0; a<nPsr; ++a) {
    g_pMainWindow->m_oConstants.n += pnPsrObs[a];
    g_pMainWindow->m_oConstants.nMarParameters += pnTMPars[a];
    strcpy(g_pMainWindow->m_oConstants.poPulsars[a].strPulsarName, pstrPulsarNames[a]);

    for(int p=0; p<pnTMPars[a]; ++p) {
      strcpy(g_pMainWindow->m_oConstants.poPulsars[a].pstrTempo2Descriptions[p], pstrTMPDescription[nIndexTM]);
      nIndexTM++;
    } // for p
  } // for a
  nIndexTM = 0;

  // Allocate memory for the data
  g_pMainWindow->m_oData.vdData.Initialize(g_pMainWindow->m_oConstants.n);
  g_pMainWindow->m_oData.vdDataErr.Initialize(g_pMainWindow->m_oConstants.n);

  nIndexTM = 0; nIndex = 0;
  // Set the number of observations per pulsar, and copy the data
  for(int a=0; a<nPsr; ++a) {
    g_pMainWindow->m_oConstants.poPulsars[a].nObservations = pnPsrObs[a];

    g_pMainWindow->m_oConstants.poPulsars[a].dPhi = pdPhi[a];
    g_pMainWindow->m_oConstants.poPulsars[a].dTheta = pdTheta[a];

    // Set the design matrices. (We'll assume we always have at least one
    // TMParameter)
    g_pMainWindow->m_oConstants.poPulsars[a].nTempo2Parameters = pnTMPars[a];
    g_pMainWindow->m_oConstants.poPulsars[a].mdTempo2ParameterDerivative.Initialize(pnPsrObs[a], pnTMPars[a]);

    // Set the actual data
    for(int i=0; i<pnPsrObs[a]; ++i) {
      g_pMainWindow->m_oConstants.poPulsars[a].pdTOA[i] = pdTOAs[nIndex];
      g_pMainWindow->m_oConstants.poPulsars[a].pdFreq[i] = pdFreq[nIndex];
      g_pMainWindow->m_oConstants.poPulsars[a].pdResiduals[i] = pdPreRes[nIndex];
      g_pMainWindow->m_oConstants.poPulsars[a].pdDeltaResiduals[i] = pdErr[nIndex];
      g_pMainWindow->m_oData.vdData[nIndex] = pdPreRes[nIndex];
      g_pMainWindow->m_oData.vdDataErr[nIndex] = pdErr[nIndex];
      strcpy(g_pMainWindow->m_oConstants.poPulsars[a].pstrFlags[i], pstrFlags[nIndex]);

      if(strlen(pstrFlags[nIndex]) > 0) {
	g_pMainWindow->m_oConstants.poPulsars[a].pbFlagSet[i] = true;
      } else {
	g_pMainWindow->m_oConstants.poPulsars[a].pbFlagSet[i] = false;
      } // if strlen

      for(int t=0; t<pnTMPars[a]; ++t) {
	g_pMainWindow->m_oConstants.poPulsars[a].mdTempo2ParameterDerivative[i][t] = pdDesignMatrix[nIndexTM+t + g_pMainWindow->m_oConstants.nMarParameters * nIndex];
      } // for t

      nIndex++;
    } // for i
    nIndexTM += pnTMPars[a];
  } // for a

  // Set the design matrix M (internally), and create the matrix of Hellings &
  // Downs coefficients
  g_pMainWindow->SetMatrices(false);

  // Set all the reduced basis and ABC-method stuff to not-initialised
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    /*
    g_pMainWindow->m_oConstants.poPulsars[a].bInterpolatorsSet = false;
    g_pMainWindow->m_oConstants.poPulsars[a].nParEfacIndex = -1;
    g_pMainWindow->m_oConstants.poPulsars[a].nParEquadIndex = -1;
    g_pMainWindow->m_oConstants.poPulsars[a].nParPlIndex = -1;
    g_pMainWindow->m_oConstants.poPulsars[a].nUniques = 0;
    g_pMainWindow->m_oConstants.poPulsars[a].ppIntAccel = NULL;
    g_pMainWindow->m_oConstants.poPulsars[a].ppIntSpline = NULL;
    */

    // Define the G-matrices (this is not even the SVD-one, as in van Haasteren
    // et al. (2012). But not used. Use van Haasteren 2009 likelihood
    g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce.Initialize(g_pMainWindow->m_oConstants.poPulsars[a].nObservations, g_pMainWindow->m_oConstants.poPulsars[a].nObservations);
    for(int i=0; i<g_pMainWindow->m_oConstants.poPulsars[a].nObservations; ++i) {
      for(int j=0; j<g_pMainWindow->m_oConstants.poPulsars[a].nObservations; ++j) {
	if(i==j) {
	  g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce[i][j] = 1.0;
	} else {
	  g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce[i][j] = 0.0;
	} /* if i==j */
      } /* for j */
    } /* for i */
  } // for a

  // Process this non-reduced basis
  ProcessNewBasis();
} // SetPulsarData

/* This function sets the compression matrix */
extern "C" void SetPulsarCompression(int nPulsar, int nRedObs, double *pdCompressionBasis) {
  g_pMainWindow->m_oConstants.poPulsars[nPulsar].mdGReduce.Initialize(
      g_pMainWindow->m_oConstants.poPulsars[nPulsar].nObservations,
      nRedObs);
  for(int i=0; i<g_pMainWindow->m_oConstants.poPulsars[nPulsar].nObservations; ++i) {
    for(int j=0; j<nRedObs; ++j) {
      g_pMainWindow->m_oConstants.poPulsars[nPulsar].mdGReduce[i][j] =
	pdCompressionBasis[j+i*nRedObs];
    } // for j
  } // for i

  g_pMainWindow->m_oConstants.bUseReducedBasis = true;
} // SetPulsarCompression

/* This function sets all the derived data compression matrices after the bases
 * have been defined. Also do this if there is no compression yet, for
 * compatibility
 */
extern "C" void ProcessNewBasis() {
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    /* Also define/allocate the associate matrices */
    g_pMainWindow->m_oConstants.poPulsars[a].mdGReduceT = g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce[LO_TRANSPOSE];
    g_pMainWindow->m_oConstants.poPulsars[a].mdC.Initialize(
	g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce.m_pnDimSize[0],
	g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce.m_pnDimSize[0]);
    g_pMainWindow->m_oConstants.poPulsars[a].mdGCG.Initialize(
	g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1],
	g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1]);
  } // for a

  /* Also produce a full mdGReduce matrix */
  int nReduce=0, nIndex1=0, nIndex2=0;
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    nReduce += g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce.m_pnDimSize[1];
  } // for a

  g_pMainWindow->m_oData.mdGReduce.Initialize(g_pMainWindow->m_oConstants.n, nReduce);

  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    for(int i=0; i<g_pMainWindow->m_oConstants.poPulsars[a].nObservations; ++i) {
      nIndex2=0;
      for(int b=0; b<g_pMainWindow->m_oConstants.k; ++b) {
	for(int j=0; j<g_pMainWindow->m_oConstants.poPulsars[b].mdGReduce.m_pnDimSize[1]; ++j) {
	  if(a == b) {
	    g_pMainWindow->m_oData.mdGReduce[nIndex1][nIndex2] = double(g_pMainWindow->m_oConstants.poPulsars[a].mdGReduce[i][j]);
	  } else {
	    g_pMainWindow->m_oData.mdGReduce[nIndex1][nIndex2] = 0;
	  } /* if a == b */
	  nIndex2++;
	} /* for j */
      } /* for b */
      nIndex1++;
    } /* for i */
  } /* for a */

  g_pMainWindow->m_oData.mdGReduceT = g_pMainWindow->m_oData.mdGReduce[LO_TRANSPOSE];
  g_pMainWindow->m_oData.vdDataReduced = g_pMainWindow->m_oData.mdGReduceT * g_pMainWindow->m_oData.vdData;
} // ProcessNewBasis


/* This function adds a source that describes the tempo2 timing model parameters
 * or any other linear deterministic source
 *
 * This is sort of incompatible with the new HDF5 file format. This is described
 * internally on a per-pulsar level, though in the HDF5 file they are already
 * all lumped together. However, in the Python code we'll have to extract the
 * information from the "pulsar"-part of the data tree
 *
 * NOTE: it is not completely obvious to what extend this information is
 * relevant. The tempo2 model is required for the code to run, but since the
 * full design matrix and the nMarParameters are set correctly, these are always
 * dealt with.
 *
 * WARNING: Do not do fancy stuff. Just always marginalise over all parameters
 * and it will be fine. (Don't include parameters numerically)
 */
extern "C" void AddLinearSource(int nPsr, int nTMP) {
  int s = g_pMainWindow->m_oConstants.nSources;
  int si = -1;
    
  g_pMainWindow->AddSource(); // Now: s = nSources - 1

  /* oSourceTypes is a global variable already declared in corefunctions.h
   * Determine which index belongs to a deterministic source (is 0 now)
   */
  for(int k=0; strcmp(oSourceTypes[k].strID, "end") != 0; k++) {
    if(strcmp(oSourceTypes[k].strID, "det") == 0) {
      si = k;
    } // if strcmp
  } // for k

  // Standard initialisation for tempo2 timind model parameters
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.eID = oSourceTypes[si].eID;
  strcpy(g_pMainWindow->m_oConstants.poSources[s].oSourceType.strID, "det");
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nParameters = nTMP;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nMarTag = FITPARAMETER_ALL;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nWrapTag = oSourceTypes[si].nWrapTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nTag = DETSOURCE_TEMPO2;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nReduceTag = FITPARAMETER_NONE;

  // There is no correlation between pulsar. Work only on single pulsars
  g_pMainWindow->m_oConstants.poSources[s].eCorrelation = CID_SinglePulsar;

  // Set the scope
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    if(a == nPsr) {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = true;
    } else {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = false;
    } // if a
  } // for a

  strcpy(g_pMainWindow->m_oConstants.poSources[s].strScopeTag, "");
  // g_pMainWindow->m_oConstants.poSources[s].nFirstParIndex = 0;
  SetNumberOfParametersFromSources(g_pMainWindow->m_oConstants, g_pMainWindow->m_oData);
} // AddLinearSource

/* This function calculates the covariance matrix, given the stochastic model
 * parameters
 */
extern "C" void CovarianceMatrix(double *pdCovariance, double *pdPar) {
  if(g_pMainWindow) {
    // Set all the varying parameters
    for(int p=0; p<g_nDimSize; ++p) {
      g_pMainWindow->m_oParameters.pdPar[g_pnParameterIndex[p]] = pdPar[p];
    } // for p

    SetCoherenceMatrix(g_pMainWindow->m_oData, g_pMainWindow->m_oParameters, g_pMainWindow->m_oConstants);

    for(int i=0; i<g_pMainWindow->m_oData.mdC.m_pnDimSize[0]*g_pMainWindow->m_oData.mdC.m_pnDimSize[1]; ++i) {
        pdCovariance[i] = g_pMainWindow->m_oData.mdC.m_pdData[i];
    } // for i
  } else {
      printf("WARNING: g_MainWindow not set\n");
  }// if g_pMainWindow
} // CovarianceMatrix

/* This function calculates the compression matrix
 */
extern "C" int CompressionMatrix(double *pdCompressedBasis, double dFidelity, int nReducePoints) {
  CMatrix mdResult;
  ReduceData(g_pMainWindow->m_oData, g_pMainWindow->m_oConstants, mdResult, dFidelity, nReducePoints);

  /* The returned matrix is returned as a row-major matrix, because of memory
   * continuity. It will be transposed in the Python code */
  for(int i=0; i<mdResult.m_pnDimSize[0]; ++i) {
    for(int j=0; j<mdResult.m_pnDimSize[1]; ++j) {
      pdCompressedBasis[i+mdResult.m_pnDimSize[0]*j] = double(mdResult[i][j]);
    } // for j
  } // for i

  return mdResult.m_pnDimSize[1];
} // CompressionMatrix

/* This function calculates the unique elements of the covariance matrix of a
 * single (power-law) source at a specified number or bins. This is used for
 * initialising the interpolation cubic splines
 */
extern "C" void CalcInterpolData(int nPulsar, int nSource, int nInterpolationBins, double *pdData, double *pdGamma) {
  CalcSourceInterpolationBins(g_pMainWindow->m_oData, g_pMainWindow->m_oConstants, nPulsar, nSource, nInterpolationBins, pdData, pdGamma);
} // CalcInterpolData

/* This function intialises the GSL cubic spline interpolators for a specific
 * source
 */
extern "C" void SetSourceInterpolator(int nPulsar, int nSource, int nInterpolationBins, double *pdData, double *pdGamma) {
  SetSourceInterpolation(g_pMainWindow->m_oData, g_pMainWindow->m_oConstants, nPulsar, nSource, nInterpolationBins, pdData, pdGamma);
} // SetSourceInterpolator

extern "C" void SetSourceAmpAccelerator(int nPulsar, int nSource) {
  SetSourceAmpAccelerator(g_pMainWindow->m_oData, g_pMainWindow->m_oConstants, nPulsar, nSource);
} // SetSourceAmpAccelerator

/* Add an efac source, connected to the error bars.
 * Applies to the flag given by strFlag
 */
extern "C" void AddEfacSource(char strFlag[], int nPsr, double dMin, double dMax, double dStart, double dStepSize, int nVary, bool bRed) {
  int s = g_pMainWindow->m_oConstants.nSources;
  int si = -1;
    
  g_pMainWindow->AddSource(); // Now: s = nSources - 1

  /* oSourceTypes is a global variable already declared in corefunctions.h
   * Determine which index belongs to a deterministic source (is 0 now)
   */
  for(int k=0; strcmp(oSourceTypes[k].strID, "end") != 0; k++) {
    if(strcmp(oSourceTypes[k].strID, "err") == 0) {
      si = k;
    } // if strcmp
  } // for k

  // Standard initialisation for efac parameters
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.eID = oSourceTypes[si].eID;
  strcpy(g_pMainWindow->m_oConstants.poSources[s].oSourceType.strID, "err");
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nParameters = oSourceTypes[si].nParameters;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = (nVary ? FITPARAMETER(0) : FITPARAMETER_NONE);
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nMarTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nWrapTag = oSourceTypes[si].nWrapTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nTag = oSourceTypes[si].nTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nReduceTag = (bRed ? FITPARAMETER(0) : FITPARAMETER_NONE);

  // There is no correlation between pulsar. Work only on single pulsars
  g_pMainWindow->m_oConstants.poSources[s].eCorrelation = CID_SinglePulsar;

  // Set the tag/flag for this source
  strcpy(g_pMainWindow->m_oConstants.poSources[s].strScopeTag, strFlag);

  // Set the scope or do it with a tag
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    // Only do this if strFlag == ""
    if(a == nPsr && strcmp(strFlag, "") == 0) {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = true;
    } else {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = false;
    } // if a
  } // for a

  // Set all the boundary conditions
  g_pMainWindow->m_oConstants.poSources[s].pdPar[0] = dStart;
  g_pMainWindow->m_oConstants.poSources[s].pdParMinBound[0] = dMin;
  g_pMainWindow->m_oConstants.poSources[s].pdParMaxBound[0] = dMax;
  g_pMainWindow->m_oConstants.poSources[s].pdParStart[0] = dStart;
  g_pMainWindow->m_oConstants.poSources[s].pdParWidthMCMC[0] = dStepSize;
  g_pMainWindow->m_oConstants.poSources[s].pdParWidthFit[0] = dStepSize;
  g_pMainWindow->m_oConstants.poSources[s].pdParWidthPrior[0] = 0;
  g_pMainWindow->m_oConstants.poSources[s].pdParMeanPrior[0] = 0;

  SetNumberOfParametersFromSources(g_pMainWindow->m_oConstants, g_pMainWindow->m_oData);
} // AddEfacSource

/* Add an equad source
 * Applies to the flag given by strFlag
 */
extern "C" void AddEquadSource(char strFlag[], int nPsr, double dMin, double dMax, double dStart, double dStepSize, int nVary, bool bRed) {
  int s = g_pMainWindow->m_oConstants.nSources;
  int si = -1;
    
  g_pMainWindow->AddSource(); // Now: s = nSources - 1

  /* oSourceTypes is a global variable already declared in corefunctions.h
   * Determine which index belongs to a deterministic source (is 0 now)
   */
  for(int k=0; strcmp(oSourceTypes[k].strID, "end") != 0; k++) {
    if(strcmp(oSourceTypes[k].strID, "wit") == 0) {
      si = k;
    } // if strcmp
  } // for k

  // Standard initialisation for equad parameters
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.eID = oSourceTypes[si].eID;
  strcpy(g_pMainWindow->m_oConstants.poSources[s].oSourceType.strID, "wit");
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nParameters = oSourceTypes[si].nParameters;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = (nVary ? FITPARAMETER(0) : FITPARAMETER_NONE);
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nMarTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nWrapTag = oSourceTypes[si].nWrapTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nTag = oSourceTypes[si].nTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nReduceTag = (bRed ? FITPARAMETER(0) : FITPARAMETER_NONE);

  // There is no correlation between pulsar. Work only on single pulsars
  g_pMainWindow->m_oConstants.poSources[s].eCorrelation = CID_SinglePulsar;

  // Set the tag/flag for this source
  strcpy(g_pMainWindow->m_oConstants.poSources[s].strScopeTag, strFlag);

  // Set the scope or do it with a tag
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    // Only do this if strFlag == ""
    if(a == nPsr && strcmp(strFlag, "") == 0) {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = true;
    } else {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = false;
    } // if a
  } // for a

  // Set all the boundary conditions
  g_pMainWindow->m_oConstants.poSources[s].pdPar[0] = dStart;
  g_pMainWindow->m_oConstants.poSources[s].pdParMinBound[0] = dMin;
  g_pMainWindow->m_oConstants.poSources[s].pdParMaxBound[0] = dMax;
  g_pMainWindow->m_oConstants.poSources[s].pdParStart[0] = dStart;
  g_pMainWindow->m_oConstants.poSources[s].pdParWidthMCMC[0] = dStepSize;
  g_pMainWindow->m_oConstants.poSources[s].pdParWidthFit[0] = dStepSize;
  g_pMainWindow->m_oConstants.poSources[s].pdParWidthPrior[0] = 0;
  g_pMainWindow->m_oConstants.poSources[s].pdParMeanPrior[0] = 0;

  SetNumberOfParametersFromSources(g_pMainWindow->m_oConstants, g_pMainWindow->m_oData);
} // AddEquadSource


/* Add a red noise source (power-law)
 * Applies to a single pulsar, given by nPsr
 *
 * strFlag: the tag for the TOAs it applies to. If == "", use nPsr
 * nPsr: The pulsar it applies to
 * pdMin: minimum bound for all three parameters
 * ... etc.
 *
 * TODO: Perhaps also implement it so that it can be applied to a tag?
 */
extern "C" void AddRedNoiseSource(char strFlag[], int nPsr, double *pdMin, double *pdMax, double *pdStart, double *pdStepSize, long *pnVary, bool bRed) {
  int s = g_pMainWindow->m_oConstants.nSources;
  int si = -1;
    
  g_pMainWindow->AddSource(); // Now: s = nSources - 1

  /* oSourceTypes is a global variable already declared in corefunctions.h
   * Determine which index belongs to a deterministic source (is 0 now)
   */
  for(int k=0; strcmp(oSourceTypes[k].strID, "end") != 0; k++) {
    if(strcmp(oSourceTypes[k].strID, "pow") == 0) {
      si = k;
    } // if strcmp
  } // for k

  // Standard initialisation for equad parameters
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.eID = oSourceTypes[si].eID;
  strcpy(g_pMainWindow->m_oConstants.poSources[s].oSourceType.strID, "pow");
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nParameters = oSourceTypes[si].nParameters;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nMarTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nWrapTag = oSourceTypes[si].nWrapTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nTag = oSourceTypes[si].nTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nReduceTag = (bRed ? FITPARAMETER(0) : FITPARAMETER_NONE);

  // There is no correlation between pulsar. Work only on single pulsars
  g_pMainWindow->m_oConstants.poSources[s].eCorrelation = CID_SinglePulsar;

  // Set the tag/flag for this source
  strcpy(g_pMainWindow->m_oConstants.poSources[s].strScopeTag, strFlag);

  // Set the scope or do it with a tag
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    // Only do this if strFlag == ""
    if(a == nPsr && strcmp(strFlag, "") == 0) {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = true;
    } else {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = false;
    } // if a
  } // for a

  // Set all the boundary conditions
  for(int p=0; p<oSourceTypes[si].nParameters; ++p) {
    g_pMainWindow->m_oConstants.poSources[s].pdPar[p] = pdStart[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParMinBound[p] = pdMin[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParMaxBound[p] = pdMax[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParStart[p] = pdStart[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthMCMC[p] = pdStepSize[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthFit[p] = pdStepSize[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthPrior[p] = 0;
    g_pMainWindow->m_oConstants.poSources[s].pdParMeanPrior[p] = 0;

    if(pnVary[p]) {
      g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag | FITPARAMETER(p);
    } // if pnVary
  } // for p

  SetNumberOfParametersFromSources(g_pMainWindow->m_oConstants, g_pMainWindow->m_oData);
} // AddRedNoiseSource


/* Add a dispersion measure variation source (power-law in dmv)
 * Applies to a single pulsar, given by nPsr
 *
 * strFlag: the tag for the TOAs it applies to. If == "", use nPsr
 * nPsr: The pulsar it applies to
 * pdMin: minimum bound for all three parameters
 * ... etc.
 */
extern "C" void AddDMVSource(char strFlag[], int nPsr, double *pdMin, double *pdMax, double *pdStart, double *pdStepSize, long *pnVary, bool bRed) {
  int s = g_pMainWindow->m_oConstants.nSources;
  int si = -1;
    
  g_pMainWindow->AddSource(); // Now: s = nSources - 1

  /* oSourceTypes is a global variable already declared in corefunctions.h
   * Determine which index belongs to a deterministic source (is 0 now)
   */
  for(int k=0; strcmp(oSourceTypes[k].strID, "end") != 0; k++) {
    if(strcmp(oSourceTypes[k].strID, "dmv") == 0) {
      si = k;
    } // if strcmp
  } // for k

  // Standard initialisation for equad parameters
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.eID = oSourceTypes[si].eID;
  strcpy(g_pMainWindow->m_oConstants.poSources[s].oSourceType.strID, "dmv");
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nParameters = oSourceTypes[si].nParameters;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nMarTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nWrapTag = oSourceTypes[si].nWrapTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nTag = oSourceTypes[si].nTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nReduceTag = (bRed ? FITPARAMETER(0) : FITPARAMETER_NONE);

  // There is no correlation between pulsar. Work only on single pulsars
  g_pMainWindow->m_oConstants.poSources[s].eCorrelation = CID_SinglePulsar;

  // Set the tag/flag for this source
  strcpy(g_pMainWindow->m_oConstants.poSources[s].strScopeTag, strFlag);

  // Set the scope or do it with a tag
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    // Only do this if strFlag == ""
    if(a == nPsr && strcmp(strFlag, "") == 0) {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = true;
    } else {
      g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = false;
    } // if a
  } // for a

  // Set all the boundary conditions
  for(int p=0; p<oSourceTypes[si].nParameters; ++p) {
    g_pMainWindow->m_oConstants.poSources[s].pdPar[p] = pdStart[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParMinBound[p] = pdMin[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParMaxBound[p] = pdMax[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParStart[p] = pdStart[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthMCMC[p] = pdStepSize[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthFit[p] = pdStepSize[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthPrior[p] = 0;
    g_pMainWindow->m_oConstants.poSources[s].pdParMeanPrior[p] = 0;

    if(pnVary[p]) {
      g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag | FITPARAMETER(p);
    } // if pnVary
  } // for p

  SetNumberOfParametersFromSources(g_pMainWindow->m_oConstants, g_pMainWindow->m_oData);
} // AddDMVSource


/* Add a GWB source (power-law)
 */
extern "C" void AddGWBSource(double *pdMin, double *pdMax, double *pdStart, double *pdStepSize, long *pnVary, bool bRed) {
  int s = g_pMainWindow->m_oConstants.nSources;
  int si = -1;
    
  g_pMainWindow->AddSource(); // Now: s = nSources - 1

  /* oSourceTypes is a global variable already declared in corefunctions.h
   * Determine which index belongs to a deterministic source (is 0 now)
   */
  for(int k=0; strcmp(oSourceTypes[k].strID, "end") != 0; k++) {
    if(strcmp(oSourceTypes[k].strID, "pow") == 0) {
      si = k;
    } // if strcmp
  } // for k

  // Standard initialisation for equad parameters
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.eID = oSourceTypes[si].eID;
  strcpy(g_pMainWindow->m_oConstants.poSources[s].oSourceType.strID, "pow");
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nParameters = oSourceTypes[si].nParameters;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nMarTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nWrapTag = oSourceTypes[si].nWrapTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nTag = oSourceTypes[si].nTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nReduceTag = (bRed ? FITPARAMETER(0) : FITPARAMETER_NONE);

  // There is no correlation between pulsar. Work only on single pulsars
  g_pMainWindow->m_oConstants.poSources[s].eCorrelation = CID_GR;

  // Set the scope to none (do it with a tag)
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = true;
  } // for a

  // Set the tag/flag for this source
  strcpy(g_pMainWindow->m_oConstants.poSources[s].strScopeTag, "");

  // Set all the boundary conditions
  for(int p=0; p<oSourceTypes[si].nParameters; ++p) {
    g_pMainWindow->m_oConstants.poSources[s].pdPar[p] = pdStart[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParMinBound[p] = pdMin[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParMaxBound[p] = pdMax[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParStart[p] = pdStart[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthMCMC[p] = pdStepSize[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthFit[p] = pdStepSize[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthPrior[p] = 0;
    g_pMainWindow->m_oConstants.poSources[s].pdParMeanPrior[p] = 0;

    if(pnVary[p]) {
      g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag | FITPARAMETER(p);
    } // if pnVary
  } // for p

  SetNumberOfParametersFromSources(g_pMainWindow->m_oConstants, g_pMainWindow->m_oData);
} // AddGWBSource

/* Add an anisotropic GWB source (power-law)
 *
 * nMode: 0 = 0, 0
 *        1 = 1, -1
 *        2 = 1, 0
 *        3 = 1, 1
 */
extern "C" void AddDipoleGWBSource(int nMode, double *pdMin, double *pdMax, double *pdStart, double *pdStepSize, long *pnVary, bool bRed) {
  int s = g_pMainWindow->m_oConstants.nSources;
  int si = -1;
    
  g_pMainWindow->AddSource(); // Now: s = nSources - 1

  /* oSourceTypes is a global variable already declared in corefunctions.h
   * Determine which index belongs to a deterministic source (is 0 now)
   */
  for(int k=0; strcmp(oSourceTypes[k].strID, "end") != 0; k++) {
    if(strcmp(oSourceTypes[k].strID, "agw") == 0) {
      si = k;
    } // if strcmp
  } // for k

  // Standard initialisation for equad parameters
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.eID = oSourceTypes[si].eID;
  strcpy(g_pMainWindow->m_oConstants.poSources[s].oSourceType.strID, "agw");
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nParameters = oSourceTypes[si].nParameters;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nMarTag = FITPARAMETER_NONE;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nWrapTag = oSourceTypes[si].nWrapTag;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nTag = nMode;
  g_pMainWindow->m_oConstants.poSources[s].oSourceType.nReduceTag = (bRed ? FITPARAMETER(0) : FITPARAMETER_NONE);

  // There is no correlation between pulsar. Work only on single pulsars
  g_pMainWindow->m_oConstants.poSources[s].eCorrelation = CID_Uniform;

  // Set the scope to all (don't do it with a tag)
  for(int a=0; a<g_pMainWindow->m_oConstants.k; ++a) {
    g_pMainWindow->m_oConstants.poSources[s].pbScope[a] = true;
  } // for a

  // Set the tag/flag for this source
  strcpy(g_pMainWindow->m_oConstants.poSources[s].strScopeTag, "");

  // Set all the boundary conditions
  for(int p=0; p<oSourceTypes[si].nParameters; ++p) {
    g_pMainWindow->m_oConstants.poSources[s].pdPar[p] = pdStart[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParMinBound[p] = pdMin[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParMaxBound[p] = pdMax[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParStart[p] = pdStart[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthMCMC[p] = pdStepSize[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthFit[p] = pdStepSize[p];
    g_pMainWindow->m_oConstants.poSources[s].pdParWidthPrior[p] = 0;
    g_pMainWindow->m_oConstants.poSources[s].pdParMeanPrior[p] = 0;

    if(pnVary[p]) {
      g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag = g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag | FITPARAMETER(p);
    } // if pnVary
  } // for p

  SetNumberOfParametersFromSources(g_pMainWindow->m_oConstants, g_pMainWindow->m_oData);
} // AddDipoleGWBSource

/* This function does some extra necessary stuff to prepare for likelihood calls
 *
 * It also processes the reduced basis, so call after everything is read-in
 */
extern "C" void CompleteInitWrapper() {
  // Set the number of parameters
  SetNumberOfParametersFromSources(g_pMainWindow->m_oConstants, g_pMainWindow->m_oData);

  // Is this necessary? Beats me...
  g_pMainWindow->m_oData.mdC.Initialize(g_pMainWindow->m_oConstants.n,
      g_pMainWindow->m_oConstants.n);

  // Get the number of dimensions
  g_nDimSize = NumberOfVaryingParameters(g_pMainWindow->m_oConstants);

  // Set the index translation array
  if(g_pnParameterIndex) {
    delete[] g_pnParameterIndex;
    g_pnParameterIndex = NULL;
  } // if g_pnParameterIndex
  g_pnParameterIndex = new int[g_nDimSize];
  SetParameterIndexArray(g_pMainWindow->m_oConstants, g_pnParameterIndex);

  // Set all the parameter boundary values (start, width, bounds, etc)
  StartParametersFromSources(g_pMainWindow->m_oParameters,
      g_pMainWindow->m_oConstants.poSources,
      g_pMainWindow->m_oConstants.nSources);

  // Set _all_ the parameter values to their start values
  ParametersFromStartParameters(g_pMainWindow->m_oParameters,
      g_pMainWindow->m_oConstants.poSources,
      g_pMainWindow->m_oConstants.nSources);

  // Set the design matrix M (internally)
  SetBlockMatrix(g_pMainWindow->m_oData, g_pMainWindow->m_oConstants);

  // Process this non-reduced basis
  ProcessNewBasis();
} // CompleteInitWrapper





/* This function returns the number of dimensions in the MCMC
 */
extern "C" int GetParameterDimensions() {
  if(g_pMainWindow == NULL) {
    return 0;
  } else {
    return g_nDimSize;
  } // if g_pMainWindow
} // GetParameterDimensions

/* Obtain the boundary conditions, like min, max, start and width
 */
extern "C" void GetBoundaryConditions(double *pdMin, double *pdMax, double *pdStart, double *pdWidth) {
  if(g_pMainWindow != NULL) {
    for(int p=0; p<g_nDimSize; ++p) {
      pdMin[p] = g_pMainWindow->m_oParameters.pdParMinBound[g_pnParameterIndex[p]];
      pdMax[p] = g_pMainWindow->m_oParameters.pdParMaxBound[g_pnParameterIndex[p]];
      pdStart[p] = g_pMainWindow->m_oParameters.pdParStart[g_pnParameterIndex[p]];
      pdWidth[p] = g_pMainWindow->m_oParameters.pdParWidthMCMC[g_pnParameterIndex[p]];
          // * g_pMainWindow->m_oConstants.dGlobalMCMCWidthFactor;
    } // for p
  } // if g_pMainWindow
} // GetBoundaryConditions

/* Obtain the amplitude parameters (parameters that indicate the intensity of a
 * signal). This includes efac, equad/white, power-law amp... but not the
 * spectral indices.
 *
 * the integer array is 0 for a parameter that is not an amplitude, and 1 for an
 * amplitude parameter. Length of given array should be g_nDimSize
 */
extern "C" void GetAmplitudeIndices(double *pdIndices) {
  if(g_pMainWindow != NULL) {
    for(int p=0; p<g_nDimSize; ++p) {
      if(IsAmplitudeParameter(g_pMainWindow->m_oConstants, p)) {
	pdIndices[p] = 1.0;
      } else {
	pdIndices[p] = 0.0;
      } // if p
    } // for p
  } // if g_pMainWindow
} // GetAmplitudeIndices

/* This function wraps the loglikelihood times prior
 */
extern "C" double WrapLL(double *pdPar) {
  double dReturnValue=0;

  if(g_pMainWindow) {
    // Set all the varying parameters
    for(int p=0; p<g_nDimSize; ++p) {
      g_pMainWindow->m_oParameters.pdPar[g_pnParameterIndex[p]] = pdPar[p];
    } // for p

    dReturnValue = -LogLikelihoodTimesPrior(g_pMainWindow->m_oData,
	g_pMainWindow->m_oParameters,
	g_pMainWindow->m_oConstants);
  } // if g_pMainWindow

  return dReturnValue;
} // WrapLL

/* This function wraps the loglikelihood times prior, with respect to hypercube
 * parameters (on the interval 0, 1
 */
extern "C" double WrapHCLL(double *pdPar) {
  double dReturnValue=0;

  if(g_pMainWindow) {
    // Set all the varying parameters
    for(int p=0; p<g_nDimSize; ++p) {
      g_pMainWindow->m_oParameters.pdPar[g_pnParameterIndex[p]] =
	g_pMainWindow->m_oParameters.pdParMinBound[g_pnParameterIndex[p]] +
	(
	  g_pMainWindow->m_oParameters.pdParMaxBound[g_pnParameterIndex[p]] -
	  g_pMainWindow->m_oParameters.pdParMinBound[g_pnParameterIndex[p]]
	) * pdPar[p];
    } // for p

    dReturnValue = -LogLikelihoodTimesPrior(g_pMainWindow->m_oData,
	g_pMainWindow->m_oParameters,
	g_pMainWindow->m_oConstants);
  } // if g_pMainWindow

  return dReturnValue;
} // WrapLL

/* This function outputs all the relevant information about the model
 */
extern "C" void PrintModel() {
  printf("nSources: %i\n", g_pMainWindow->m_oConstants.nSources);
  if(g_pMainWindow) {
    for(int s=0; s<g_pMainWindow->m_oConstants.nSources; ++s) {
      printf("source[%2i]\n", s);
      printf("  eID:          %i\n", g_pMainWindow->m_oConstants.poSources[s].oSourceType.eID);
      printf("  eCorrelation: %i\n", g_pMainWindow->m_oConstants.poSources[s].eCorrelation);
      printf("  strID:        %s\n", g_pMainWindow->m_oConstants.poSources[s].oSourceType.strID);
      printf("  nParameters:  %i\n", g_pMainWindow->m_oConstants.poSources[s].oSourceType.nParameters);
      printf("  nFitTag:      %lu\n", g_pMainWindow->m_oConstants.poSources[s].oSourceType.nFitTag);
      printf("  nMarTag:      %lu\n", g_pMainWindow->m_oConstants.poSources[s].oSourceType.nMarTag);
      printf("  nWrapTag:     %lu\n", g_pMainWindow->m_oConstants.poSources[s].oSourceType.nWrapTag);
      printf("  nTag:         %li\n", g_pMainWindow->m_oConstants.poSources[s].oSourceType.nTag);
      printf("  nReduceTag:   %lu\n", g_pMainWindow->m_oConstants.poSources[s].oSourceType.nReduceTag);
      printf("  eCorr:        %i\n", g_pMainWindow->m_oConstants.poSources[s].eCorrelation);
      printf("  strScopeTag:  %s\n", g_pMainWindow->m_oConstants.poSources[s].strScopeTag);
      printf("  nFirstParI:   %i\n", g_pMainWindow->m_oConstants.poSources[s].nFirstParIndex);
    } // for s
  } // if g_pMainWindow
} // PrintModel

/* Set whether or not we will use the reduced basis formalism when calculating
 * the log-likelihood. This way, even for a compressed dataset, also the full
 * likelihood can be evaluated
 */
extern "C" void SetUseReducedBasis(bool bUseReducedBasis) {
  g_pMainWindow->m_oConstants.bUseReducedBasis = bUseReducedBasis;
} // SetUseReducedBasis

/* Clean up all memory */
extern "C" void FinishWrapper() {
  if(g_pMainWindow) {
    delete g_pMainWindow;
    g_pMainWindow = NULL;
  } // if g_pMainWindow

  g_nDimSize = 0;

  if(g_pnParameterIndex) {
    delete g_pnParameterIndex;
    g_pnParameterIndex = NULL;
  } // if g_pnParameterIndex
} // FinishWrapper

