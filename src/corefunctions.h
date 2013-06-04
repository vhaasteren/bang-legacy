/* corefunctions.h -- Core header for PTA calculations program
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


#ifndef __COREFUNCTIONS_H__
#define __COREFUNCTIONS_H__

#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_cdf.h>
//#include "tempo2.h"

// After how many steps we should save the MCMC data (also because of memory
// considerations
#define MAX_MCMC_BUFFER			2000
#define CLOUD_SIZE_PER_DIM		5

// General ban constants
// TODO: partially got rid of this non-dynamically allocated nonsense
//       only MAX_PARAMETERS and MAX_PARAMETERS_PER_SOURCE are left in code
//       High values for these does not seem to be an issue
#define MAX_PULSARS			40
#define MAX_SOURCES			200
#define MAX_PARAMETERS_PER_SOURCE	1000
#define MAX_PARAMETERS			10000
#define MAX_OBSERVATIONS_PER_PULSAR	1000

// Which parameter we want to fit for
// -2: ALL
// -1: 0
// x : x (for all other values)
// Combine with bitwise OR: |
#define FITPARAMETER(x)			(x == -1 ? 0 : ( x == -2 ? ~0 : 1 << x ))
#define FITPARAMETER_ALL		~0
#define FITPARAMETER_NONE		0

// Which deterministic source we are dealing with
// -2: ALL known
// -1: None
// x : x (for all other values)
// Combine with bitwise OR: |
#define DETSOURCENUMBER		9

#define DETSOURCE(x)			(x == -1 ? 0 : ( x == -2 ? ~0 : 1 << x ))
#define DETSOURCE_ALL		~0
#define DETSOURCE_NONE		0

#define DETSOURCE_CONSTANT	DETSOURCE(0)
#define DETSOURCE_LINEAR	DETSOURCE(1)
#define DETSOURCE_QUADRATIC	DETSOURCE(2)
#define DETSOURCE_CUBIC		DETSOURCE(3)
#define DETSOURCE_QUARTIC	DETSOURCE(4)
#define DETSOURCE_SINE		DETSOURCE(5)
#define DETSOURCE_GWMEM		DETSOURCE(6)
#define DETSOURCE_TEMPO2	DETSOURCE(7)
#define DETSOURCE_BHBINARY	DETSOURCE(8)

// QSD = Constant + Linear + Quadratic
#define DETSOURCE_QSD		7


// The default name of the configfile (parameters.conf)
const char g_strConfigFile[] = "parameters.conf";
const char g_strBaseDir[] = "./";


// Types of sources (Use NoSource to quantify fitting for parameters like
// pulsar position)
//
// Types of sources
// ----------------
//
// SID_Deterministic:	A deterministic source, this causes extra residuals to
// 			be subtracted
// SID_Nothing:		Not a source
// SID_White:		White noise
// SID_PowerLaw:	A power-law (for GWB). The cut-off frequency is not a
//			variable parameter by default
// SID_Exponential:	Exponential spectrum
// SID_Lorentzian:	A Lorentzian spectrum
// SID_PowerLaw_OneParameter: See power-law, with only an amplitude
// SID_Errorbars:	Error-bars are a form of white-noise. Use one parameter
//			as an overall scaling factor
// SID_NoSource:	Parameters that are not really a source, but these can
//			be adjusted and fit for anyway. Example is the
//			correlation (Zeta) between 2 pulsars to calculate a
//			Bayesian version of the Jenet et al. statistic.
// SID_NonStationary	A non-stationary Gaussian source (modelled on
// 			exponential source) - for testing purposes
// SID_DMVar:       DM variations, modelled as a power-law
// SID_DipoleGWB:   Anisotropic, Dipole, GWB, power-law (nTag = mode)
enum ESourceID {SID_Deterministic=0, SID_Nothing, SID_White, SID_PowerLaw, SID_Exponential, SID_Lorentzian, SID_PowerLaw_OneParameter, SID_Errorbars, SID_NoSource, SID_NonStationary, SID_DMVar, SID_DipoleGWB};

// Types of correlations
enum ECorrelationID {CID_SinglePulsar=0, CID_All, CID_Uniform, CID_GR, CID_Metric_Breathing, CID_Metric_Shear, CID_Metric_Longitudinal};

// The names in the configuration file for the above enumeration types (more
// user-friendly)
const char strCorrelationID[][4] = {
  {"sin"},
  {"all"},
  {"uni"},
  {"gr"},
  {"br"},
  {"sh"},
  {"lon"},
  {"end"}
};


// Predefined deterministic sources. This list is intended to let the user
// specify these with the help of strings in a configuration file. This also
// prevents having to specify each source separately
struct SDetSourceType {
  int nID;			// Type of source - integer code
  char strID[4];		// Type of source - 3 character code
  int nParameters;		// Amount of parameters
  bool bLinear;			// Whether this source is linear
}; // struct SDetSourceType


// Types of deterministic sources
// ------------------------------
//
// DETSOURCE_CONSTANT:	A constant offset added to all resdiduals (sec)
// DETSOURCE_LINEAR:	A linear term, with zero-point at zero (sec / yr)
// DETSOURCE_QUADRATIC: A quadratic term, with zero-point at zero (sec / yr^2)
// DETSOURCE_QUBIC:	A qubic term, with zero-point at zero (sec / yr^3)
// DETSOURCE_QUARTIC:	A quartic term, with zero-point at zero (sec / yr^4)
// DETSOURCE_QSD:	A quadratic spindown term: con + lin + qua
// DETSOURCE_SINE:	A sinusoidal term, with an amplitude, frequency and
// 			phase. Thus a non-linear source
// DETSOURCE_GWMEM:	A signal coming from GW memory
// 			p[0]: time the gwmem shockwave hits the earth
// 			p[1]: the size of the metric perturbation h
// 			p[2]: raj (rad) of source of gwmem shockwave
// 			p[3]: decl (rad) of source of gwmem shockwave
// 			p[4]: polarisation of the gwmem shockwave (Eulerian
// 			angle, rad)
// DETSOURCE_TEMPO2:	All tempo2 deterministic sources (only usable in the
//                      plugin). Note that the number of parameters of this
//                      source is seet in tempo2_compat.cpp
const SDetSourceType oDetSourceTypes[] = {
  {DETSOURCE_CONSTANT,		"con",		1, 	true},
  {DETSOURCE_LINEAR,		"lin",		1, 	true},
  {DETSOURCE_QUADRATIC,		"qua",		1, 	true},
  {DETSOURCE_CUBIC,		"qub",		1, 	true},
  {DETSOURCE_QUARTIC,		"qrt",		1, 	true},
  {DETSOURCE_QSD,		"qsd",		3, 	true},
  {DETSOURCE_SINE,		"sin",		3, 	false},
  {DETSOURCE_GWMEM,		"mem",		5, 	false},
  {DETSOURCE_TEMPO2,		"tem",		MAX_PARAMETERS_PER_SOURCE, 	true},
  {DETSOURCE_BHBINARY,		"bhb",		MAX_PARAMETERS_PER_SOURCE, 	false},
  {0, "end", 0, false}
};


// Predefined source-spectra to minimize confusion (i.e. number of parameters)
// eID:   determines the type of source of residuals (deterministic, white
//        spectrum, etc.
// strID: Same as eID, but now with chars (for in configfile)
// nParameters: The number of parameters that fully describe this source
// nFitTag:	The parameters that must be fit for. Uses binary notation, so
//              use the FITPARAMETER macro to set this. Example:
//              Only parameter 1: nFitTag = FITPARAMETER(0)    (1 = 0001)
//              Only parameter 3: nFitTag = FITPARAMETER(2)    (4 = 0100)
//              Parameter 1,2,3 : nFitTag = FITPARAMETER(0) |
//                                          FITPARAMETER(1) |
//                                          FITPARAMETER(2)    (7 = 0111)
struct SSourceType {
  ESourceID eID;		// Type of source - integer code
  char strID[4];		// Type of source - 3 character code
  int nParameters;		// Amount of parameters for this source
  unsigned long nFitTag;		// Which parameters to fit for (binary checks)
  unsigned long nMarTag;		// Which deterministic parameters we want to analytically
 				// marginalize (binary checks)
  unsigned long nWrapTag;	// Which parameters are wrapped (dMin <=> dMax)
  long nTag;			// Optional tag (i.e. what deterministic source this is, or which anisotropic mode (Sigma_{lm}))
  unsigned long nReduceTag;	// Reduce the data to be sensitive to this source
};


// Constant with source types
const SSourceType oSourceTypes[] = {
  {SID_Deterministic,		"det",	1, FITPARAMETER_NONE, FITPARAMETER_ALL,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_Nothing,			"not",	0, FITPARAMETER_ALL,  FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_White,			"wit",	1, FITPARAMETER_ALL,  FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_PowerLaw,		"pow",	3, FITPARAMETER(0) |
    					   FITPARAMETER(1),   FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_Exponential,		"exp",	2, FITPARAMETER_ALL,  FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_Lorentzian,		"lor",	2, FITPARAMETER_ALL,  FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_PowerLaw_OneParameter,	"po1",	1, FITPARAMETER_ALL,  FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_Errorbars,		"err",	1, FITPARAMETER_NONE, FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_NoSource,		"nsr",	0, FITPARAMETER_NONE, FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_NonStationary,		"nss",	3, FITPARAMETER(0) |
    					   FITPARAMETER(1),   FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_DMVar,	    	"dmv",	3, FITPARAMETER(0) |
    					   FITPARAMETER(1),   FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {SID_DipoleGWB,		"dgw",	7, FITPARAMETER(0) | FITPARAMETER(1) | FITPARAMETER(3) | FITPARAMETER(4) | FITPARAMETER(5) | FITPARAMETER(6),
    					   FITPARAMETER_NONE,  DETSOURCE_NONE, DETSOURCE_NONE, FITPARAMETER_NONE},
  {(ESourceID)0, "end", 0, 0, 0, 0, 0, 0}
};

// Acceleration of amplitude-only signals.
struct SAmplitudeAcc {
  bool bSet;				// Initialised yes/no
  int nPulsar;				// Which pulsar it applies to
  CMatrix mdGCaG;			// The matrix for unit-signal
}; // SAmplitudeAcc

// Interpolation information, can be used for interpolation of 1 parameter
struct SInterpolAcc {
  bool bSet;				// Initialised yes/no
  int nPulsar;				// Which pulsar it applies to (-1 = all)
  int nUniques;				// Number of unique elements in GCG
  CMatrix mdGCiG;			// The matrix for unit-signal
  gsl_interp_accel **ppIntAccel;	// Accelerator
  gsl_spline **ppIntSpline;		// Cubic spline interpolator
}; // SInterpol

// Sources should be defined with this structure
struct SSource {
  SSourceType oSourceType;		// The type of source
  ECorrelationID eCorrelation;		// Correlation type in case of Stochastic source
//  bool pbScope[MAX_PULSARS];		// Which pulsars are effected by this source
  bool *pbScope;		// Which pulsars are effected by this source
  char strScopeTag[128];		// Which TOAs are effected by this source
  int nFirstParIndex;			// Index of the first parameter belonging to this
  					// source

  SAmplitudeAcc oAmpAcc;		// In case it is efac/equad
  SInterpolAcc oIntAcc;			// In case it is a power-law

  double pdPar[MAX_PARAMETERS_PER_SOURCE];
  double pdParMinBound[MAX_PARAMETERS_PER_SOURCE];
  double pdParMaxBound[MAX_PARAMETERS_PER_SOURCE];
  double pdParStart[MAX_PARAMETERS_PER_SOURCE];
  double pdParWidthMCMC[MAX_PARAMETERS_PER_SOURCE];
  double pdParWidthFit[MAX_PARAMETERS_PER_SOURCE];
  double pdParWidthPrior[MAX_PARAMETERS_PER_SOURCE];
  double pdParMeanPrior[MAX_PARAMETERS_PER_SOURCE];
};


// All parameters needed when calculating another loglik value (the 2 + 2N values)
struct SParametersType {
  double pdPar[MAX_PARAMETERS];
  double pdParMinBound[MAX_PARAMETERS];
  double pdParMaxBound[MAX_PARAMETERS];
  double pdParStart[MAX_PARAMETERS];
  double pdParWidthMCMC[MAX_PARAMETERS];
  double pdParWidthFit[MAX_PARAMETERS];
  double pdParWidthPrior[MAX_PARAMETERS];
  double pdParMeanPrior[MAX_PARAMETERS];
};


// Pulsarsars should be defined with this structure
struct SPulsar {
  // Name and position
  char strPulsarName[80];		// The name of the pulsar
  char strRightAscension[80];		// Right ascension
  char strDeclination[80];		// Declination

  // Converted to workable units (mathematical/physical definition)
  double dPhi;
  double dTheta;

  // The actual data
  int nObservations;			// Number of observations
//  double pdTOA[MAX_OBSERVATIONS_PER_PULSAR];
//  double pdResiduals[MAX_OBSERVATIONS_PER_PULSAR];
//  double pdDeltaResiduals[MAX_OBSERVATIONS_PER_PULSAR];
//  bool pbFlagSet[MAX_OBSERVATIONS_PER_PULSAR];
//  char pstrFlags[MAX_OBSERVATIONS_PER_PULSAR][128];
  double *pdTOA;
  double *pdFreq;
  double *pdResiduals;
  double *pdDeltaResiduals;

  // The tags
  bool *pbFlagSet;
  char **pstrFlags;

  // The tempo2 parameter derivatives
  int nTempo2Parameters;
  CMatrix mdTempo2ParameterDerivative;
  CMatrix mdGReduce;			/* the compression matrix */
  CMatrix mdGReduceT;			/* the compression matrix */
  CMatrix mdC;				/* The single-pulsar cov. matrix */
  CMatrix mdGCG;			/* The single-pulsar cov. red. mat */

//  CMatrix mdGCeqG;			/* The C-diagonal eq cov. red. mat */
//  CMatrix mdGCefG;			/* The C-diagonal ef cov. red. mat */
//  CMatrix mdGCplG;			/* The C-power-law cov. red. mat */
//  int nParEfacIndex;
//  int nParEquadIndex;
//  int nParPlIndex;

  /* The GSL interpolation thingies */
//  bool bInterpolatorsSet;
//  int nUniques;
//  gsl_interp_accel **ppIntAccel;
//  gsl_spline **ppIntSpline;

  char pstrTempo2Descriptions[MAX_PARAMETERS_PER_SOURCE][80];
  double pdTempo2Multiplication[MAX_PARAMETERS_PER_SOURCE];
  double pdTempo2Value[MAX_PARAMETERS_PER_SOURCE];

  // Extra (not-needed) stuff
  double dPeriod;		// The period ( / frecuency)
};


// All constants needed in the pta calculation
struct SConstantsType {
  char strVersion[16];	// The version of this Package
  char strConfigVersion[16];	// The version of the configfile (paramerters.conf)

  int k;				// Number of pulsars
  int l;				// Number of datapoints per pulsar (when generating)
  int m;				// Amount of datasets
  int n;				// Total amount of datapoints.
  int nParameters;			// Total amount of parameters.
  int nFitParameters;			// Total amount of fit parameters.
  int nMarParameters;			// Total amount of fit parameters.
  int nSources;
  bool bUseReducedBasis;		// Use the reduced basis for likelihood evals
  bool bCalcTMPars;			// Calculate the ML TM parameters while mcmc
  double dReducedFidelity;		// How much information is kept (0 - 1)
  bool bCorAmpAccel;			// Whether we accalerate the correlated evals
  bool bUseInstrTags;			// Whether we use the instrtags (bug workaround)
  int nInterpolationBins;		// Number of bins used for interpolating the s.i.

//  SPulsar poPulsars[MAX_PULSARS];	// All the pulsars and their TOAs
  SPulsar *poPulsars;	// All the pulsars and their TOAs
//  SSource poSources[MAX_SOURCES];	// All the sources and their info
  SSource *poSources;	// All the sources and their info

  // Info for when generating data
  bool bUnevenlySampled;
  double dMeasureTimeInterval;	// Interval between measurements (sec - average)

  // Default pulsar parameters
  SPulsar oDefaultPulsar;
  SSource oDefaultPulsarSource;

  // General MCMC parameters
  int nMCMCSteps;		// Number of mcmc steps used in calculation
  int nPlotPoints;		// Number of points used in plot when plotting
  int nSmallPlotPoints;		// Number of points used in plot when plotting
  int nBurnInSteps;		// Number of mcmc steps before actually recording (burn-in time)
  int nAcceptanceRate;          // Acceptance rate in (%) of new points while sampling
  int nBootstrapAttempts;	// Number of trials in the bootstrap error-estimation
  int nEnsembleDataSets;	// Number of data sets in the ensemble
  double dEnsembleScalar;	// Scale size of ensemble sampler distribution
  double dGlobalMCMCWidthFactor;// Factor to multiply the width in MCMC runs with
				// for easy adjustment to the optimal acceptance rate

  int nSpindownStrength;	// The speed with which the pulsars follow qsd curves
  int nCorrectOrder;		// Number of correction terms
  int nCholeskyBits;		// Number of bits used in cholesky decomposition (precision variable)

  char strBaseDir[160];		// Base directory (./ or $TEMPO2/bayesian/)
  char strDataDir[160];		// Where all datafiles are stored
  char strTempoDataDir[160];	// Where all the tempo datafiles are stored
  char strResidualsFile[160];	// Name of the datafile
  char strParametersConf[160];	// Name of the config file
  char strAnglesFile[160];	// Name of the angles file
};



// All matices/data needed in calculating anything
struct SDataType {
  CMatrix mdGeometric;		// The geometrical part of the coherence matrix (geo subspace)
  CMatrix mdC;			// The total covariance matrix
  CMatrix mdCCor;			// The total covariance matrix
  CMatrix mdX;			// The data-matrix

  CMatrix mdEVC;		// The eigenvectors of mdC
  CMatrix mdCholeskyC;		// The Cholesky factor of mdC
  CMatrix mdInvC;		// The inverse of C
  CMatrix mdF;			// The QSD corrected version of InvC

  CMatrix mdGReduce;		/* the compression matrix */
  CMatrix mdGReduceT;		/* the compression matrix */
  CMatrix mdGCG;		/* The single-pulsar cov. red. mat */
  CMatrix mdGCGCorr;		/* The correlated GCG signal for fast eval */

  /* The GSL interpolation thingies */
//  bool bGWBInterpolatorsSet;
//  int nGWBUniques;
//  gsl_interp_accel **ppGWBIntAccel;
//  gsl_spline **ppGWBIntSpline;
//  bool bClockInterpolatorsSet;
//  int nClockUniques;
//  gsl_interp_accel **ppClockIntAccel;
//  gsl_spline **ppClockIntSpline;


  CMatrix mdBlock;		// The block matrix needed when removing functional forms
  CMatrix mdBlockTrans;		// The block matrix needed when removing functional forms
  CMatrix mdGenBlock;		// The block matrix needen when removing/adding qsd
  CMatrix mdAveBlock;		// The block matrix needen when removing/adding qsd

  CMatrix mdMarBlock;
  CMatrix mdFitBlock;

  CVector vdData;		// The residuals
  CVector vdDataReduced;
  CVector vdDataErr;		// The error bars on the data
  CVector vdLambdaC;		// The eigenvalues of mdC
  CVector vdLinearParameters;	// Linear parameters

  CVector vdChi;		// The best-fit values
  CMatrix mdCXiInv;

  // To do many dataset mcmc's
  int nDataSets;
  CMatrix mdDataSets;
  CMatrix mdChi;
  CVector vdLL;

  /* To avoid using static variables in LogLikelihoodReducedBasis
   * for their meaning, see LogLikelihoodReducedBasis
   */
  bool bRedLikFirstRun;
  int nRedLikGWBPar;
  int nRedLikAllPar;
  int nRedLikClockPar;
};

bool SourceWorksOnResidual(SConstantsType &oConstants, int s, int a, int i);
bool SourceWorksOnPulsar(SConstantsType &oConstants, int s, int a);
bool SourceWorks(SConstantsType &oConstants, int s);

void FitParametersFromVector(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec);
void VectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec);
void MinVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec);
void MaxVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec);
void WidthMCMCVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec);
void WidthFitVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec);
void WidthPriorVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec);
void MeanPriorVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec);

void ParametersToVector(SParametersType &oParameters, CVector &vdTheta);
void VectorToParameters(SParametersType &oParameters, CVector &vdTheta);
void ParametersFromSources(SParametersType &oParameters, SSource *poSources, int nSources=0);
void StartParametersFromSources(SParametersType &oParameters, SSource *poSources, int nSources=0);
void ParametersFromStartParameters(SParametersType &oParameters, SSource *poSources, int nSources=0);
int NumberOfVaryingParameters(SConstantsType &oConstants);
bool IsAmplitudeParameter(SConstantsType &oConstants, int nP);
void SetParameterIndexArray(SConstantsType &oConstants, int *pnParameterIndex);
double GenerateValidProposal(double dValue, double dWidth, double dMin, double dMax, gsl_rng *rng, bool bWrap=false);

void ReadGlobalConstants(SConstantsType &oConstants, SDataType &oData, const char *strParametersConf=NULL, const char *strResidualsFile=NULL, const char *strAnglesFile=NULL, const char *strDataDir=NULL);				// Read some constants from file
void ReadSourceConstants(SConstantsType &oConstants, SDataType &oData, const char *strParametersConf=NULL, const char *strResidualsFile=NULL, const char *strAnglesFile=NULL, const char *strDataDir=NULL);				// Read some constants from file
void SetNumberOfParametersFromSources(SConstantsType &oConstants, SDataType &oData);		// Name says it all
void SetGeometricPart(SDataType &oData, SConstantsType &oConstants);		// Set the geometric part of the coherence matrix
void CopySpecificPulsarsToObjects(SConstantsType &oConstants, SDataType &oData, SConstantsType &oNewConstants, SDataType &oNewData, bool *pbSelection);
void SetDataVectorFromPulsars(SConstantsType &oConstants, SDataType &oData);

void SetBlockMatrix(SDataType &oData, SConstantsType &oConstants);			// Set the block matrix (for quadratic spindown)
void ResidualsFromDetSources(SConstantsType &oConstants, SParametersType &oParameters, SDataType &oData, CVector *pvdReturnResiduals=NULL);
double ResidualFromDetSource(double *pdPar, SSource &oSource, int nPulsar, int nObsIndex, SDataType &oData, SConstantsType &oConstants);
double ParameterDerivative(int nPulsar, int nObsIndex, int nParameterTag, int nParameter, SDataType &oData, SConstantsType &oConstants);

double SineResidual(double *pdPar, int nPulsar, int nObsIndex, SConstantsType &oConstants);
double GWMemResidual(double *pdPar, int nPulsar, int nObsIndex, SConstantsType &oConstants);
double BHBinaryResidual(double *pdPar, int nPulsar, int nObsIndex, SConstantsType &oConstants);

void GenerateResiduals(SDataType &oData, SConstantsType &oConstants, bool bEraseDeltaResiduals=true);
void GenerateEnsembleResiduals(SDataType &oData, SConstantsType &oConstants, int nEnsembles);
void GenerateRadiometerEnsembleResiduals(SDataType &oData, SConstantsType &oConstants, int nEnsembles);
void CheckGenerateResiduals(SDataType &oData, SConstantsType &oConstants);
void GenerateDataSets(SDataType &oData, SConstantsType &oConstants);
bool WriteResiduals(SDataType &oData, SConstantsType &oConstants);			// Write the constants & data
bool ReadResiduals(SDataType &oData, SConstantsType &oConstants, SParametersType &oParameters);	// Read the data
bool ReadAppendAngles(SDataType &oData, SConstantsType &oConstants);
bool ReadAngles(SDataType &oData, SConstantsType &oConstants);
bool WriteAngles(SDataType &oData, SConstantsType &oConstants);
void GeneratePulsarAngles(SDataType &oData, SConstantsType &oConstants);
void PrintParameters(SConstantsType &oConstants);
void PrintParameters(SConstantsType &oConstants, SParametersType &oParameters);
void PrintResiduals(SConstantsType & oConstants);


void SetCoherenceMatrix(SDataType &oData, SConstantsType &oConstants, bool bTest=true); // Set the coherence matrix (data generation)
void SetCoherenceMatrix(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants); // Set the coherence matrix (estimation)
void SetCoherenceMatrixNonReduceSources(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants);
void SetCoherenceMatrixNonClockSources(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants);
void SetCoherenceMatrixClockSources(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants);
void SetCoherenceMatrixReduceSources(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants);
void SetCoherenceMatrixSinglePulsarNoCorrSig(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nPulsar);
void SetCoherenceMatrixCorrSig(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants);
void SetCoherenceMatrixOneSource(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, CMatrix &mdC, int nSource, int nPulsar=-1);


void SetCoherenceMatrixPulsarEfac(SConstantsType &oConstants, int nPulsar, CMatrix &mdC, double dEfac);
void SetCoherenceMatrixPulsarEquad(SConstantsType &oConstants, int nPulsar, CMatrix &mdC, double dEquad);
int GetSourceEquad(SConstantsType &oConstants, int nPulsar);
int GetSourceEfac(SConstantsType &oConstants, int nPulsar);
int GetSourcePowerLaw(SConstantsType &oConstants, int nPulsar);
int GetSourceGRCorr(SConstantsType &oConstants, int nPulsar=-1);
int GetSourceUniCorr(SConstantsType &oConstants, int nPulsar=-1);


void SaveCoherenceMatrix(SDataType &oData, SConstantsType &oConstants); // Save the coherence matrix (data generation)
void LoadCoherenceMatrix(SDataType &oData, SConstantsType &oConstants); // Load the coherence matrix (data generation)


void PrintTrueParameters(SDataType &oData, SConstantsType &oConstants);


// The function that sets the spectrum:
double SpectrumIntegral(ESourceID eID, double dTi, double dTj, CVector &vdParameters, bool bDMUnits=false);

// Example function for angular correlations (Chiara Mingarelli's work)
// double ChiaraExampleCorrelation(int np1, int np2, SConstantsType &oConstants, CVector &vdParameters);
double AnisotropicGWBCorrelation(int np1, int np2, int nMode, SConstantsType &oConstants, CVector &vdParameters);


// endif __COREFUNCTIONS_H__
#endif
