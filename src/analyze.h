/* analyze.h -- Analysis header for PTA calculations program
 *
 * Rutger van Haasteren 15 August 2007 haasteren@strw.leidenuniv.nl
 *
 * Copyright (C) 2006-2008 Rutger van Haasteren.
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

#ifndef __ANALYZE_H__
#define __ANALYZE_H__

#include "corefunctions.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

struct MCMCDataPoint {
  SParametersType oParameters;
  CMatrix mdCXiInv;
  CVector vdChi;
  double dLogLik;
};

class CMCMCDataPoint {
public:
  CMCMCDataPoint() { bDefined = false; }
  ~CMCMCDataPoint() {
    if(bDefined) {
      delete[] pdPar;
      if(bTMPars) {
	delete[] pdCXiInv;
	delete[] pdChi;
      } /* if bTMPars */
    } /* if bDefined */
  } /* ~CMCMCDataPoint */

  void Initialize(int nParameters, int nTMParameters, bool bReadTMPars=true) {
    bDefined = true;
    bTMPars = bReadTMPars;
    pdPar = new double[nParameters];
    if(bTMPars) {
      pdChi = new double[nTMParameters];
      pdCXiInv = new double[nTMParameters*nTMParameters];
    } /* if bTMPars */
  } /* Initialize */

  bool bDefined;
  bool bTMPars;
  double *pdPar;
  double *pdCXiInv;
  double *pdChi;
  double dLogLik;
};

struct MCMCDataEnsemble {
  SParametersType oParameters;
  CMatrix mdCXiInv;
  CMatrix mdChi;
  CVector vdLogLik;
  double dLogLik;
};

class CMCMCDataEnsemble {
public:
  CMCMCDataEnsemble() { bDefined = false; }
  ~CMCMCDataEnsemble() { if(bDefined) { delete[] pdPar; delete[] pdCXiInv; delete[] pdChi; delete[] pdLogLik; }}

  void Initialize(int nDataSets, int nParameters, int nTMParameters) { bDefined = true; pdPar = new double[nParameters]; pdChi = new double[nDataSets*nTMParameters]; pdCXiInv = new double[nTMParameters*nTMParameters]; pdLogLik = new double[nDataSets];}

  bool bDefined;
  double *pdPar;
  double *pdCXiInv;
  double *pdChi;
  double *pdLogLik;
  double dLogLik;
};


void PowerSpectrum(SDataType &oData, SConstantsType &oConstants, int nPulsar=0);

double CalcResidualsRms(SDataType &oData, SConstantsType &oConstants, int nCalcRmsPulsar=0);
double CalcResidualsChisq(SDataType &oData, SConstantsType &oConstants, int nCalcChisqPulsar=0);

// The Markov chain takes samples from a Gaussian distribution around the current sample
void MCMCGauss(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName);		// Markov Chain Monte Carlo
void EnsembleSampler(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName);		// Ensemble-sampler

// MCMC Read/Write functions (Deprecated)
//bool WriteMCMCData(int nMCMCSteps, int nParameters, MCMCDataPoint *pDat, const char *strFileName);			// Write the data to a datafile 
//bool WriteMCMCData(int nProcesses, int nMCMCSteps, int nParameters, MCMCDataPoint **ppDat, const char *strFileName);			// Write the data to a datafile
//bool ReadMCMCData(int &nMCMCSteps, int nParameters, MCMCDataPoint *&pDat, const char *strFileName);			// Read the data from a datafile

// MCMC Read/Write functions
bool WriteMCMCDataFileStart(int nParameters, int nTMParameters, const char *strFileName);
bool WriteMCMCDataFileAppend(int nMCMCSteps, int nParameters, int nTMParameters, MCMCDataPoint *pDat, const char *strFileName, bool bCalcTMPars=true);
bool ReadMCMCDataFile(int &nMCMCSteps, int nParameters, int nTMParameters, CMCMCDataPoint *&pDat, const char *strFileName, bool bCalcTMPars=true);
bool WriteMCMCEnsembleDataFileStart(int nDataSets, int nObservations, int nParameters, int nTMParameters, double *pdData, const char *strFileName);
bool WriteMCMCEnsembleDataFileAppend(int nMCMCSteps, int nDataSets, int nParameters, int nTMParameters, MCMCDataEnsemble *pDat, const char *strFileName);
bool ReadMCMCEnsembleDataFile(int &nMCMCSteps, int &nDataSets, int nParameters, int nTMParameters, CMCMCDataEnsemble *&pDat, const char *strFileName, CMatrix &mdDataSets);
bool ReadMCMCEnsembleDataFileSets(int &nDataSets, int &nMCMCSteps, const char *strFileName, CMatrix &mdDataSets);
void MCMCEnsembleImportanceSample(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nSkipSamples=0);

void Calculate1DMCMCIntegrationOnTheFly(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nPlotParameterNumber, CVector &vdX, CVector &vdY);
void Calculate1DMCMCIntegration(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nPlotParameterNumber, CVector &vdX, CVector &vdY, CVector &vdYErr, bool bCalcErr);
void IntegrateMCMCData(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nPlotParameterNumber=0, bool bCalcErr=true);
//void Integrate3DMCMCData(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nPulsarNumber=-1, int nUseParameters=3);
void Calculate3DMCMCIntegrationOnTheFly(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ);
void Calculate3DMCMCIntegration(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ);
void Calculate3DMCMCEnsembleIntegration(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, int nDataSet, CVector &vdX, CVector &vdY, CMatrix &mdZ);
void Calculate3DMCMCEnsembleIntegrationOnTheFly(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, int nDataSet, CVector &vdX, CVector &vdY, CMatrix &mdZ);
void Calculate3DMCMCEnsembleMLDRIntegrationOnTheFly(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ);
void FindParameterBounds(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameter1, int nParameter2, bool &bMarPar1, bool &bMarPar2, int &nP1, int &nP2, double &dMinBound1, double &dMinBound2, double &dMaxBound1, double &dMaxBound2, CMatrix &mdCXi);
void Calculate3DMCMCIntegrationIncMar(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, int &nPR1, int &nPR2, CVector &vdX, CVector &vdY, CMatrix &mdZ, CMatrix &mdZML, CMatrix &mdZT2);
void Calculate2DEnsembleMLPoints(SDataType &oData, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY);
void Calculate1DEnsembleIntegrationMarParameters(SDataType &oData, SConstantsType &oConstants, const char *strFileName);
void Calculate1DEnsembleIntegrationMCMCParameters(SDataType &oData, SConstantsType &oConstants, const char *strFileName);
void Calculate3DMCMCIntegrationSkymap(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ);
void Integrate3DMCMCData(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName, int nParameter1=0, int nParameter2=1);
void CalculateMCMCEvidence(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, const char *strFileName);

void CreateClockCorrectedResiduals(SDataType &oData, SConstantsType &oConstants);

double MCMCEvidence(SDataType &oData, SConstantsType &oConstants, CMCMCDataPoint *pDat, int nMCMCSteps, int nParameters);
double Integral(int nDimensions, int nMCMCSteps, double **ppdX, double *pdZ, double *pdMu, double dEnlargement, CMatrix &mdInvC, double dLogDetC, bool bPrintResults);
void CalcCovariance(int nDimensions, int nMCMCSteps, double **ppdX, double *pdZ, double *pdMu, CMatrix &mdInvC, double &dEnlargement, double &dLogDetC, bool bPrintResults);
int PointsInEllipsoid(int nDimensions, int nMCMCSteps, double **ppdX, double *pdZ, double *pdMu, CMatrix &mdInvC, double dEnlargement);
double CalcEnlargement(int nDimensions, int nMCMCSteps, double **ppdX, double *pdZ, double *pdMu, CMatrix &mdInvC, bool bPrintResults);


void RecordTime(bool bShow, const char *strComment, const int nPosition);			// Print timing data on screen

double LogLikelihood(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants);	// Calculate the LogLikelihood value
void LogLikelihoodEnsemble(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants);
double LogLikelihoodSets(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants);	// Calculate the LogLikelihood value


void Calculate1DPlot(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameterNumber, CVector &vdX, CVector &vdY);
void MakePlot(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameterNumber=0);		// Make a plot of the likelihood function

void Calculate3DPlot(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameter1, int nParameter2, CVector &vdX, CVector &vdY, CMatrix &mdZ);		// Make a plot of the likelihood function
void Make3DPlot(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameter1=0, int nParameter2=1);		// Make a plot of the likelihood function

void CalcSourceInterpolationBins(SDataType &oData, SConstantsType &oConstants, int nPulsar, int nSource, int nInterpolationBins, double *pdAllData, double *pdGamma);
void ReduceData(SDataType &oData, SConstantsType &oConstants, CMatrix &mdResult, double dFidelity=0.99, int nReducePoints=-1);
void SetSourceInterpolation(SDataType &oData, SConstantsType &oConstants, int nPulsar, int nSource, int nInterpolationBins, double *pdAllData, double *pdGamma);
void SetSourceAmpAccelerator(SDataType &oData, SConstantsType &oConstants, int nPulsar, int nSource);
void CalcPlIntMatrix(gsl_interp_accel **ppAc, gsl_spline **ppSp, double dGamma, CMatrix &mdGCplG);

#if 0
void WriteReduceMatrix(const char strFileName[], CMatrix &mdMat);
void ReadReduceMatrix(const char strFileName[], CMatrix &mdMat);
void CreatePlInterpolation(SDataType &oData, SConstantsType &oConstants);
void CreatePlInterpolTest(SDataType &oData, SConstantsType &oConstants);
void CreateGWBInterpolation(SDataType &oData, SConstantsType &oConstants);
void CreateClockInterpolation(SDataType &oData, SConstantsType &oConstants);
void ReadPlInterpolation(SDataType &oData, SConstantsType &oConstants, int nPulsar, const char strFileName[]);
void ReadGWBInterpolation(SDataType &oData, SConstantsType &oConstants, const char strFileName[]);
void ReadClockInterpolation(SDataType &oData, SConstantsType &oConstants, const char strFileName[]);
#endif

#ifdef HAVE_MPI
// Messages that can be passed in evaluation msg
#define MPI_BAN_EVAL_SUCCESS	1
#define MPI_BAN_EVAL_ELERROR	2
#define MPI_BAN_EVAL_NPOSDEF	4
#define MPI_BAN_EVAL_OTHER		8

// General messages that can be passed
#define MPI_BAN_MSG_START		1
#define MPI_BAN_MSG_QUIT		2
#define MPI_BAN_MSG_QUIT_DELAYED	4

// Tags that represent the content of the messages
#define MPI_BAN_TAG_EVAL		1
#define MPI_BAN_TAG_MSG		2

// For the MPI client/server model, we need a datatype to pass parameters and
// succes/fail flags
struct SMPIParametersType  {
  double pdPar[MAX_PARAMETERS];
  double dLogLik;
  int nStatus;
};

void MakePlotServer(MPI_Datatype *pmpiParameters, SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameterNumber);
void MPIAnalyzeClient(MPI_Datatype *pmpiParameters, SDataType &oData, SConstantsType &oConstants);
double LogLikelihood(double *pdFitPar, SDataType *poData=NULL, SParametersType *poParameters=NULL, SConstantsType *poConstants=NULL);	// Calculate the LogLikelihood value
void MCMCGaussServer(MPI_Datatype *pmpiParameters, SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nParameterNumber);
void MPIStopClients();
double LogPrior(double *pdFitPar, SParametersType *poParameters=NULL, SConstantsType *poConstants=NULL);
double LogLikelihoodTimesPrior(double *pdFitPar, SDataType *poData=NULL, SParametersType *poParameters=NULL, SConstantsType *poConstants=NULL);
#endif // HAVE_MPI

double LorentzianPrior(double dParameter, double dMean, double dWidth);
double LogPrior(SParametersType &oParameters, SConstantsType &oConstants);
double LogLikelihoodTimesPrior(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants);

// endif __ANALYZE_H__
#endif
