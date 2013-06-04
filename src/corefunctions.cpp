/* corefunctions.cpp -- Core cpp-file for PTA bayesian analysis program
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
// #include "moremath.h"
#include "linal.h"
#include "linalfunc.h"
#include "banfunc.h"
#include "corefunctions.h"
#include "filefunctions.h"


//#include "tempo2_compat.h"

using namespace std;

#define NO_DEFAULT_SOURCES


/* This function reads all parameters needed in the pta problem from a file
 * (usually this is parameters.conf).
 *
 * All data is stored in the oConstants & oData structures.
 *
 *  Last update: 2008-04-11
 *
 * */
void ReadGlobalConstants(SConstantsType &oConstants, SDataType &oData, const char *strParametersConf, const char *strResidualsFile, const char *strAnglesFile, const char *strDataDir) {				// Read some constants from file
#if 0
  string strTemp, strTemp2;
  char strBuf[160], strBuf2[160];
  int nMaxPulsarParameters=0, nMaxMeasurementsPerPulsar=0;
/*************************************
** Read the constants from a datafile
*************************************/
  if(strParametersConf && strlen(strParametersConf) != 0)
    strcpy(strBuf, strParametersConf);
  else
    sprintf(strBuf, "%s/%s", oConstants.strBaseDir, g_strConfigFile);

  printf("Opening %s...\n", strBuf);
  PrintStatus("Reading global constants...");
  ConfigFile parameters(strBuf);
  // First check whether the configfile is compatible
  parameters.readInto(strTemp, "VERSION", string("XXX"));
  strcpy(oConstants.strConfigVersion, strTemp.data());
  if(oConstants.strConfigVersion[0] == 'X') {
    PrintFailed();
    return;
  } // if 
  strcpy(oConstants.strVersion, VERSION);

  // Now read the file/directory names
  if(! parameters.readInto(strTemp, "strDataDir")) { //, "../data/workdata/");
    strTemp = "../data/workdata/";
    printf("WARNING: parameters.conf did not have: \"%s\"\n", "strDataDir");
    printf("         using: %s\n", strTemp.c_str());
  }
  sprintf(oConstants.strDataDir, "%s/%s", oConstants.strBaseDir, strTemp.data());

  if(! parameters.readInto(strTemp, "strResidualsFile")) { //, "../data/workdata/");, "residuals.dat");
    strTemp = "residuals.dat";
    printf("WARNING: parameters.conf did not have: \"%s\"\n", "strResidualsFile");
    printf("         using: %s\n", strTemp.c_str());
  }

  sprintf(oConstants.strResidualsFile, "%s/%s", oConstants.strDataDir, strTemp.data());

  if(! parameters.readInto(strTemp, "strTempoDataDir")) { //, "../data/workdata/tempo");
    strTemp = "../data/workdata/tempo";
    printf("WARNING: parameters.conf did not have: \"%s\"\n", "strTempoDataDir");
    printf("         using: %s\n", strTemp.c_str());
  }
  sprintf(oConstants.strTempoDataDir, "%s/%s", oConstants.strBaseDir, strTemp.data());
  if(! parameters.readInto(strTemp, "strAnglesFile")) { //, "../data/workdata/");, "angles.dat");
    strTemp = "angles.dat";
    printf("WARNING: parameters.conf did not have: \"%s\"\n", "strAnglesFile");
    printf("         using: %s\n", strTemp.c_str());
  }
  sprintf(oConstants.strAnglesFile, "%s/%s", oConstants.strDataDir, strTemp.data());
//  if(! parameters.readInto(strTemp, "strParametersConf")) { // "../data/parameters.conf
//    strTemp = "../data/parameters.conf";
//    printf("WARNING: parameters.conf did not have: \"%s\"\n", "strParametersConf");
//    printf("         using: %s\n", strTemp.c_str());
//  }
//  sprintf(oConstants.strParametersConf, "%s/%s", oConstants.strDataDir, strTemp.data());

  // Check whether we have a custom data-dir:
  if(strDataDir && strlen(strDataDir) != 0)
    strcpy(oConstants.strDataDir, strDataDir);
  if(strAnglesFile && strlen(strAnglesFile) != 0)
    strcpy(oConstants.strAnglesFile, strAnglesFile);
  if(strResidualsFile && strlen(strResidualsFile) != 0)
    strcpy(oConstants.strResidualsFile, strResidualsFile);
  if(strParametersConf && strlen(strParametersConf) != 0)
    strcpy(oConstants.strParametersConf, strParametersConf);
  else
    sprintf(oConstants.strParametersConf, "%s/%s", oConstants.strBaseDir, g_strConfigFile);


  // Now read the global parameters
//  parameters.readInto(oConstants.k, "k");
  parameters.readInto(strTemp, "k");
  oConstants.k = atoi(strTemp.c_str());
//  if(k > MAX_PULSARS) printf("WARNING: k > MAX_PULSARS\n");
//  parameters.readInto(oConstants.l, "l");
  parameters.readInto(strTemp, "l");
  oConstants.l = atoi(strTemp.c_str());
//  parameters.readInto(oConstants.m, "m");
  parameters.readInto(strTemp, "m");
  oConstants.m = atoi(strTemp.c_str());
  parameters.readInto(strTemp, "bUnevenlySampled", string("true"));
  if(strcmp(strTemp.data(), "false") == 0)
    oConstants.bUnevenlySampled = false;
  else
    oConstants.bUnevenlySampled = true;
  oConstants.nParameters = 0;
  oConstants.nFitParameters = 0;
//  parameters.readInto(oConstants.dMeasureTimeInterval, "dMeasureTimeInterval");
  parameters.readInto(strTemp, "dMeasureTimeInterval");
  oConstants.dMeasureTimeInterval = atof(strTemp.c_str());

  if(! parameters.readInto(strTemp, "bUseReducedBasis")) {
    strTemp = "false";
    printf("WARNING: parameters.conf did not have: \"%s\"\n", "bUseReducedBasis");
    printf("         using: %s\n", strTemp.c_str());
  }
  if(strcmp(strTemp.data(), "false") == 0)
    oConstants.bUseReducedBasis = false;
  else
    oConstants.bUseReducedBasis = true;
//  if(! parameters.readInto(oConstants.dReducedFidelity, "dReducedFidelity")) {
  if(! parameters.readInto(strTemp, "dReducedFidelity")) {
    oConstants.dReducedFidelity = 0.999;
  } else {
    oConstants.dReducedFidelity = atof(strTemp.c_str());
  }
  if(! parameters.readInto(strTemp, "bCorAmpAccel")) {
    strTemp = "false";
  }
  if(strcmp(strTemp.data(), "false") == 0)
    oConstants.bCorAmpAccel = false;
  else
    oConstants.bCorAmpAccel = true;
  if(! parameters.readInto(strTemp, "bCalcTMPars")) {
    strTemp = "true";
  }
  if(strcmp(strTemp.data(), "false") == 0)
    oConstants.bCalcTMPars = false;
  else
    oConstants.bCalcTMPars = true;
  if(! parameters.readInto(strTemp, "bUseInstrTags")) {
    strTemp = "true";
  }
  if(strcmp(strTemp.data(), "true") == 0)
    oConstants.bUseInstrTags = true;
  else
    oConstants.bUseInstrTags = false;
//  if(! parameters.readInto(oConstants.nInterpolationBins, "nInterpolationBins")) {
  if(! parameters.readInto(strTemp, "nInterpolationBins")) {
    oConstants.nInterpolationBins = 1000;
  } else {
    oConstants.nInterpolationBins = atoi(strTemp.c_str());
  }

  for(int a=0; a<oConstants.k; a++) {
    // Read the name of the pulsar
    sprintf(strBuf, "strPulsarName[%i]", a);
    sprintf(strBuf2, "Pulsar-%i", a);
    strTemp2 = strBuf2;
    parameters.readInto(strTemp, strBuf, strTemp2);
    strcpy(oConstants.poPulsars[a].strPulsarName, strTemp.c_str());
  } // for a


  PrintSuccess();

  ReadSourceConstants(oConstants, oData, strParametersConf, strResidualsFile, strAnglesFile, strDataDir);

  PrintStatus("Reading remaining constants...");

  oConstants.oDefaultPulsar.nObservations = oConstants.l;
  oConstants.oDefaultPulsar.dPhi = 0;
  oConstants.oDefaultPulsar.dTheta = 0;
  for(int i=0; i<MAX_OBSERVATIONS_PER_PULSAR; i++) {
    oConstants.oDefaultPulsar.pdTOA[i] = 0;
    oConstants.oDefaultPulsar.pbFlagSet[i] = false;
  }

  // Now check how many datapoints all pulsars have
  oConstants.n = oConstants.k * oConstants.l;

//  parameters.readInto(oConstants.nPlotPoints, "nPlotPoints");
  parameters.readInto(strTemp, "nPlotPoints");
  oConstants.nPlotPoints = atoi(strTemp.c_str());
//  parameters.readInto(oConstants.nSmallPlotPoints, "nSmallPlotPoints");
  parameters.readInto(strTemp, "nSmallPlotPoints");
  oConstants.nSmallPlotPoints = atoi(strTemp.c_str());
//  parameters.readInto(oConstants.nMCMCSteps, "nMCMCSteps");
  parameters.readInto(strTemp, "nMCMCSteps");
  oConstants.nMCMCSteps = atoi(strTemp.c_str());
//  parameters.readInto(oConstants.nBurnInSteps, "nBurnInSteps");
  parameters.readInto(strTemp, "nBurnInSteps");
  oConstants.nBurnInSteps = atoi(strTemp.c_str());
//  parameters.readInto(oConstants.nAcceptanceRate, "nAcceptanceRate");
  parameters.readInto(strTemp, "nAcceptanceRate");
  oConstants.nAcceptanceRate = atoi(strTemp.c_str());
  parameters.readInto(oConstants.nBootstrapAttempts, "nBootstrapAttempts");
  parameters.readInto(strTemp, "nBootstrapAttempts");
  oConstants.nBootstrapAttempts = atoi(strTemp.c_str());
//  parameters.readInto(oConstants.nEnsembleDataSets, "nEnsembleDataSets");
  parameters.readInto(strTemp, "nEnsembleDataSets");
  oConstants.nEnsembleDataSets = atoi(strTemp.c_str());
//  if(! parameters.readInto(oConstants.dEnsembleScalar, "dEnsembleScalar")) {
  if(! parameters.readInto(strTemp, "dEnsembleScalar")) {
    oConstants.dEnsembleScalar = 2.0;
    printf("WARNING: parameters.conf did not have: \"%s\"\n", "dEnsembleScalar");
    printf("         using: %f\n", oConstants.dEnsembleScalar);
  } else {
    oConstants.dEnsembleScalar = atof(strTemp.c_str());
  }
//  if(! parameters.readInto(oConstants.dGlobalMCMCWidthFactor, "dGlobalMCMCWidthFactor")) {
  if(! parameters.readInto(strTemp, "dGlobalMCMCWidthFactor")) {
    oConstants.dGlobalMCMCWidthFactor = 1.0;
  } else {
    oConstants.dGlobalMCMCWidthFactor = atof(strTemp.c_str());
  }
//  parameters.readInto(oConstants.nCholeskyBits, "nCholeskyBits");
  parameters.readInto(strTemp, "nCholeskyBits");
  oConstants.nCholeskyBits = atoi(strTemp.c_str());

  PrintSuccess();
/*************************************
** Read the constants from a datafile
*************************************/


  // First set the measure times. Original MATLAB statement says it all:
  // vdMeasureTimes = cumsum(2*rand(l,1)); % Cumulate random numbers as measure times
  CVector vdTemp;
  vdTemp.Initialize(oConstants.l);
//  parameters.readInto(oConstants.m, "m");
  parameters.readInto(strTemp, "m");
  oConstants.m = atoi(strTemp.c_str());

  vdTemp.Randomize();
//  vdTemp *= 2;
  for(int i=0; i<oConstants.l; i++) {
    if(i != 0)
      oConstants.oDefaultPulsar.pdTOA[i] = oConstants.oDefaultPulsar.pdTOA[i-1] + oConstants.dMeasureTimeInterval*double(vdTemp[i]);
    else
      oConstants.oDefaultPulsar.pdTOA[i] = oConstants.dMeasureTimeInterval + 76;
  }  // for

  for(int a=0; a<oConstants.k; a++) {
    strcpy(strBuf, oConstants.poPulsars[a].strPulsarName);
    oConstants.poPulsars[a] = oConstants.oDefaultPulsar;
    strcpy(oConstants.poPulsars[a].strPulsarName, strBuf);

    if(oConstants.bUnevenlySampled) {
      vdTemp.Randomize();
//  vdTemp *= 2;
      oConstants.poPulsars[a].nObservations = oConstants.l;
      for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	oConstants.poPulsars[a].pbFlagSet[i] = false;
	if(i != 0)
	  oConstants.poPulsars[a].pdTOA[i] = oConstants.poPulsars[a].pdTOA[i-1] + oConstants.dMeasureTimeInterval*double(vdTemp[i]);
	else
	  oConstants.poPulsars[a].pdTOA[i] = oConstants.dMeasureTimeInterval + 76;
      }  // for
    }
  } // for a

#if 0
    oData.mdGeometric.Initialize(oConstants.k,oConstants.k);
    oData.mdC.Initialize(oConstants.n, oConstants.n);
    oData.mdInvC.Initialize(oConstants.n, oConstants.n);
#endif
#endif
} // ReadGlobalConstants



/* This function reads all parameters belonging to specific sources. Usually
 * everything is read from a 'parameters.conf' file.
 *
 * Examples of sources are:
 * - Quadratic spindown
 * - GWB
 * - Pulsar timing noise
 *
 *  Last update: 2008-10-23 */
void ReadSourceConstants(SConstantsType &oConstants, SDataType &oData, const char *strParametersConf, const char *strResidualsFile, const char *strAnglesFile, const char *strDataDir) {
#if 0
  // TODO: Clear these up
  string strTemp;
  char strBuf[160];
  int nMaxPulsarParameters=0, nMaxMeasurementsPerPulsar=0;
  int nPulsar, nMarParameter, nDetSourceType;
  bool bFound=false;
  CNumber ndTemp;

  if(strParametersConf && strlen(strParametersConf) != 0)
    strcpy(strBuf, strParametersConf);
  else
    strcpy(strBuf, oConstants.strParametersConf);

  PrintStatus("Reading source constants...");

  ConfigFile parameters(strBuf);


  // First read default source constants
  if(! parameters.readInto(strTemp, "eDefaultSourceType")) {
    strTemp = "wit";
    printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
    printf("         using: %s\n", strTemp.c_str());
  }
  for(int i=0; strcmp(oSourceTypes[i].strID, "end") != 0; i++) {
    if(strcmp(oSourceTypes[i].strID, strTemp.c_str()) == 0) {
      oConstants.oDefaultPulsarSource.oSourceType.eID = oSourceTypes[i].eID;
      strcpy(oConstants.oDefaultPulsarSource.oSourceType.strID, oSourceTypes[i].strID);
      oConstants.oDefaultPulsarSource.oSourceType.nParameters = oSourceTypes[i].nParameters;
      oConstants.oDefaultPulsarSource.oSourceType.nFitTag = oSourceTypes[i].nFitTag;
      oConstants.oDefaultPulsarSource.oSourceType.nMarTag = oSourceTypes[i].nMarTag;
      oConstants.oDefaultPulsarSource.oSourceType.nWrapTag = oSourceTypes[i].nWrapTag;
      oConstants.oDefaultPulsarSource.oSourceType.nTag = oSourceTypes[i].nTag;
      oConstants.oDefaultPulsarSource.oSourceType.nReduceTag = oSourceTypes[i].nReduceTag;
    } // if strcmp
  } // for i

  strcpy(strBuf, "eDefaultCorrelation");
  if(! parameters.readInto(strTemp, strBuf)) {
    strTemp = "sin";
    printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
    printf("         using: %s\n", strTemp.c_str());
  }
  for(int i=0; strcmp(strCorrelationID[i], "end") != 0; i++) {
    if(strcmp(strCorrelationID[i], strTemp.c_str()) == 0) {
      oConstants.oDefaultPulsarSource.eCorrelation = (ECorrelationID) i;
    } // if strcmp

    switch(oConstants.oDefaultPulsarSource.eCorrelation) {
      case CID_SinglePulsar:
	for(int a=0; a<oConstants.k; a++)
	  oConstants.oDefaultPulsarSource.pbScope[a] = false;
	oConstants.oDefaultPulsarSource.pbScope[0] = true;
	strcpy(oConstants.oDefaultPulsarSource.strScopeTag, "");
	break;
      case CID_All:
      case CID_Uniform:
      case CID_GR:
      case CID_Metric_Breathing:
      case CID_Metric_Shear:
      case CID_Metric_Longitudinal:
	for(int a=0; a<oConstants.k; a++)
	  oConstants.oDefaultPulsarSource.pbScope[a] == true;
	strcpy(oConstants.oDefaultPulsarSource.strScopeTag, "");
	break;
      default:
	for(int a=0; a<oConstants.k; a++)
	  oConstants.oDefaultPulsarSource.pbScope[a] = false;
	oConstants.oDefaultPulsarSource.pbScope[0] = true;
	strcpy(oConstants.oDefaultPulsarSource.strScopeTag, "");
	break;
    } // switch eCorrelation
  } // for i

  for(int i=0; i<oConstants.oDefaultPulsarSource.oSourceType.nParameters; i++) {
    sprintf(strBuf, "dDefaultPar[%i]", i);
//    if(! parameters.readInto(oConstants.oDefaultPulsarSource.pdPar[i], strBuf)) {
    if(! parameters.readInto(strTemp, strBuf)) {
      oConstants.oDefaultPulsarSource.pdPar[i] = 0;
      printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
      printf("         using: %f\n", oConstants.oDefaultPulsarSource.pdPar[i]);
    } else {
      oConstants.oDefaultPulsarSource.pdPar[i] = atof(strTemp.c_str());
    }// if parameters.readInto

    sprintf(strBuf, "dDefaultParMinBound[%i]", i);
//    if(! parameters.readInto(oConstants.oDefaultPulsarSource.pdParMinBound[i], strBuf)) {
    if(! parameters.readInto(strTemp, strBuf)) {
      oConstants.oDefaultPulsarSource.pdParMinBound[i] = 0;
      printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
      printf("         using: %f\n", oConstants.oDefaultPulsarSource.pdParMinBound[i]);
    } else {
      oConstants.oDefaultPulsarSource.pdParMinBound[i] = atof(strTemp.c_str());
    }// if parameters.readInto

    sprintf(strBuf, "dDefaultParMaxBound[%i]", i);
//    if(! parameters.readInto(oConstants.oDefaultPulsarSource.pdParMaxBound[i], strBuf)) {
    if(! parameters.readInto(strTemp, strBuf)) {
      oConstants.oDefaultPulsarSource.pdParMaxBound[i] = 0;
      printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
      printf("         using: %f\n", oConstants.oDefaultPulsarSource.pdParMaxBound[i]);
    } else {
      oConstants.oDefaultPulsarSource.pdParMaxBound[i] = atof(strTemp.c_str());
    }// if parameters.readInto

    sprintf(strBuf, "dDefaultParStart[%i]", i);
//    if(! parameters.readInto(oConstants.oDefaultPulsarSource.pdParStart[i], strBuf)) {
    if(! parameters.readInto(strTemp, strBuf)) {
      oConstants.oDefaultPulsarSource.pdParStart[i] = 0;
      printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
      printf("         using: %f\n", oConstants.oDefaultPulsarSource.pdParStart[i]);
    } else {
      oConstants.oDefaultPulsarSource.pdParStart[i] = atof(strTemp.c_str());
    }// if parameters.readInto

    sprintf(strBuf, "dDefaultParWidthMCMC[%i]", i);
//    if(! parameters.readInto(oConstants.oDefaultPulsarSource.pdParWidthMCMC[i], strBuf)) {
    if(! parameters.readInto(strTemp, strBuf)) {
      oConstants.oDefaultPulsarSource.pdParWidthMCMC[i] = 0;
      printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
      printf("         using: %f\n", oConstants.oDefaultPulsarSource.pdParWidthMCMC[i]);
    } else {
      oConstants.oDefaultPulsarSource.pdParWidthMCMC[i] = atof(strTemp.c_str());
    }// if parameters.readInto

    sprintf(strBuf, "dDefaultParWidthFit[%i]", i);
//    if(! parameters.readInto(oConstants.oDefaultPulsarSource.pdParWidthFit[i], strBuf)) {
    if(! parameters.readInto(strTemp, strBuf)) {
      oConstants.oDefaultPulsarSource.pdParWidthFit[i] = 0;
      printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
      printf("         using: %f\n", oConstants.oDefaultPulsarSource.pdParWidthFit[i]);
    } else {
      oConstants.oDefaultPulsarSource.pdParWidthFit[i] = atof(strTemp.c_str());
    }// if parameters.readInto

    sprintf(strBuf, "dDefaultParWidthPrior[%i]", i);
//    if(! parameters.readInto(oConstants.oDefaultPulsarSource.pdParWidthPrior[i], strBuf)) {
    if(! parameters.readInto(strTemp, strBuf)) {
      oConstants.oDefaultPulsarSource.pdParWidthPrior[i] = 0;
      printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
      printf("         using: %f\n", oConstants.oDefaultPulsarSource.pdParWidthPrior[i]);
    } else {
      oConstants.oDefaultPulsarSource.pdParWidthPrior[i] = atof(strTemp.c_str());
    }// if parameters.readInto

    sprintf(strBuf, "dDefaultParMeanPrior[%i]", i);
//    if(! parameters.readInto(oConstants.oDefaultPulsarSource.pdParMeanPrior[i], strBuf)) {
    if(! parameters.readInto(strTemp, strBuf)) {
      oConstants.oDefaultPulsarSource.pdParMeanPrior[i] = 0;
      printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
      printf("         using: %f\n", oConstants.oDefaultPulsarSource.pdParMeanPrior[i]);
    } else {
      oConstants.oDefaultPulsarSource.pdParMeanPrior[i] = atof(strTemp.c_str());
    } // if parameters.readInto

    sprintf(strBuf, "bDefaultParFit[%i]", i);
    if(parameters.readInto(strTemp, strBuf)) {
      if(strcmp(strTemp.c_str(), "true") == 0) {
	oConstants.oDefaultPulsarSource.oSourceType.nFitTag =
	  oConstants.oDefaultPulsarSource.oSourceType.nFitTag | FITPARAMETER(i);
      } else {
	oConstants.oDefaultPulsarSource.oSourceType.nFitTag =
	  oConstants.oDefaultPulsarSource.oSourceType.nFitTag & (~FITPARAMETER(i));
      } // if strcmp
    } // if parameters.readInto

    sprintf(strBuf, "bDefaultParMar[%i]", i);
    if(parameters.readInto(strTemp, strBuf)) {
      if(strcmp(strTemp.c_str(), "true") == 0) {
	oConstants.oDefaultPulsarSource.oSourceType.nMarTag =
	  oConstants.oDefaultPulsarSource.oSourceType.nMarTag | FITPARAMETER(i);
      } else {
	oConstants.oDefaultPulsarSource.oSourceType.nMarTag =
	  oConstants.oDefaultPulsarSource.oSourceType.nMarTag & (~FITPARAMETER(i));
      } // if strcmp
    } // if parameters.readInto

    sprintf(strBuf, "bDefaultParWrap[%i]", i);
    if(parameters.readInto(strTemp, strBuf)) {
      if(strcmp(strTemp.c_str(), "true") == 0) {
	oConstants.oDefaultPulsarSource.oSourceType.nWrapTag =
	  oConstants.oDefaultPulsarSource.oSourceType.nWrapTag | FITPARAMETER(i);
      } else {
	oConstants.oDefaultPulsarSource.oSourceType.nWrapTag =
	  oConstants.oDefaultPulsarSource.oSourceType.nWrapTag & (~FITPARAMETER(i));
      } // if strcmp
    } // if parameters.readInto

    sprintf(strBuf, "bDefaultParReduce[%i]", i);
    if(parameters.readInto(strTemp, strBuf)) {
      if(strcmp(strTemp.c_str(), "true") == 0) {
	oConstants.oDefaultPulsarSource.oSourceType.nReduceTag =
	  oConstants.oDefaultPulsarSource.oSourceType.nReduceTag | FITPARAMETER(i);
      } else {
	oConstants.oDefaultPulsarSource.oSourceType.nReduceTag =
	  oConstants.oDefaultPulsarSource.oSourceType.nReduceTag & (~FITPARAMETER(i));
      } // if strcmp
    } // if parameters.readInto
  } // for i

  // Now read all the sources
  oConstants.nSources = 0;
  for(int i=0; ; i++) {
    sprintf(strBuf, "eSourceType[%i]", i);
    if(! parameters.readInto(strTemp, strBuf)) {
      break;
    } // if parameters.readInto

    // We have another source. Initialize
    oConstants.nSources++;
    for(int j=0; j<MAX_PULSARS; j++)
      oConstants.poSources[i].pbScope[j] = false;
    strcpy(oConstants.poSources[i].strScopeTag, "");

    // Determine the source type
    for(int k=0; strcmp(oSourceTypes[k].strID, "end") != 0; k++) {
      if(strcmp(oSourceTypes[k].strID, strTemp.c_str()) == 0) {
	oConstants.poSources[i].oSourceType.eID = oSourceTypes[k].eID;
	strcpy(oConstants.poSources[i].oSourceType.strID, oSourceTypes[k].strID);
	oConstants.poSources[i].oSourceType.nParameters = oSourceTypes[k].nParameters;
	oConstants.poSources[i].oSourceType.nFitTag = oSourceTypes[k].nFitTag;
	oConstants.poSources[i].oSourceType.nMarTag = oSourceTypes[k].nMarTag;
	oConstants.poSources[i].oSourceType.nWrapTag = oSourceTypes[k].nWrapTag;
	oConstants.poSources[i].oSourceType.nTag = oSourceTypes[k].nTag;
	oConstants.poSources[i].oSourceType.nReduceTag = oSourceTypes[k].nReduceTag;
      } // if strcmp
    } // for k

    if(strcmp(oConstants.poSources[i].oSourceType.strID, oSourceTypes[(int)SID_Deterministic].strID) == 0) {
      // This is a deterministic source
      sprintf(strBuf, "eDetSourceType[%i]", i);
      if(! parameters.readInto(strTemp, strBuf)) {
	oConstants.poSources[i].oSourceType.eID = SID_Nothing;
	oConstants.poSources[i].oSourceType.nParameters = 0;
	strcpy(oConstants.poSources[i].oSourceType.strID, oSourceTypes[SID_Nothing].strID);
      } else {
	for(int k=0; strcmp(oDetSourceTypes[k].strID, "end") != 0; k++) {
	  if(strcmp(oDetSourceTypes[k].strID, strTemp.c_str()) == 0) {
	    oConstants.poSources[i].oSourceType.nTag = oDetSourceTypes[k].nID;
	    oConstants.poSources[i].oSourceType.nParameters = oDetSourceTypes[k].nParameters;
	    if(oDetSourceTypes[k].bLinear) {
	      oConstants.poSources[i].oSourceType.nFitTag = FITPARAMETER_NONE;
	      oConstants.poSources[i].oSourceType.nMarTag = FITPARAMETER_ALL;
	    } else {
	      oConstants.poSources[i].oSourceType.nFitTag = FITPARAMETER_ALL;
	      oConstants.poSources[i].oSourceType.nMarTag = FITPARAMETER_NONE;
	    } //  if bLinear
	  } // if strcmp
	} // for k
      } // if ! readInto
    } // if strcmp strID

    // Determine the correlation type
    sprintf(strBuf, "eCorrelation[%i]", i);
    if(! parameters.readInto(strTemp, strBuf)) {
      strTemp = "sin";
      printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
      printf("         using: %s\n", strTemp.c_str());
    } // if parameters.readInto
    for(int k=0; strcmp(strCorrelationID[k], "end") != 0; k++) {
      if(strcmp(strCorrelationID[k], strTemp.c_str()) == 0) {
	oConstants.poSources[i].eCorrelation = (ECorrelationID) k;
      } // if strcmp
    } // for k


    switch(oConstants.poSources[i].eCorrelation) {
      case CID_All:
      case CID_Uniform:
      case CID_GR:
      case CID_Metric_Breathing:
      case CID_Metric_Shear:
      case CID_Metric_Longitudinal:
	for(int a=0; a<oConstants.k; a++)
	  oConstants.poSources[i].pbScope[a] = true;
	break;
      case CID_SinglePulsar:
      default:
	// Set the pulsars this source influences
	for(int a=0; a<oConstants.k; a++)
	  oConstants.poSources[i].pbScope[a] = false;
	sprintf(strBuf, "strScopeTag[%i]", i);
	if(! parameters.readInto(strTemp, strBuf)) {
	  sprintf(strBuf, "nPulsar[%i]", i);

//	  if(! parameters.readInto(nPulsar, strBuf)) {
	  if(! parameters.readInto(strTemp, strBuf)) {
	    sprintf(strBuf, "strPulsar[%i]", i);
	    if(! parameters.readInto(strTemp, strBuf)) {
	      nPulsar = -1;
	    } else {
	      nPulsar = -1;
	      for(int a=0; a<oConstants.k; a++) {
		if(strcmp(strTemp.c_str(), oConstants.poPulsars[a].strPulsarName) == 0) {
		  nPulsar = a;
		  break;
		} // if strcmp
	      } // for a
	    } // if parameters
	  } else {
	    nPulsar = atoi(strTemp.c_str());
	  }// if nPulsar
	  if(nPulsar >= 0)
	    oConstants.poSources[i].pbScope[nPulsar] = true;
	} else {
	  strcpy(oConstants.poSources[i].strScopeTag, strTemp.c_str());
	} // if strScopeTag
	break;
    } // switch eCorrelation


    if(oConstants.poSources[i].oSourceType.eID != SID_Deterministic ||
	oConstants.poSources[i].oSourceType.nTag != DETSOURCE_TEMPO2) {
      // Now read the parameters of this source
      for(int j=0; j<oConstants.poSources[i].oSourceType.nParameters; j++) {
	sprintf(strBuf, "dPar[%i][%i]", i, j);
//	if(! parameters.readInto(oConstants.poSources[i].pdPar[j], strBuf)) {
	if(! parameters.readInto(strTemp, strBuf)) {
	  oConstants.poSources[i].pdPar[j] = 0;
	  printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
	  printf("         using: %f\n", oConstants.poSources[i].pdPar[j]);
	} else {
	  oConstants.poSources[i].pdPar[j] = atof(strTemp.c_str());
	} // if parameters.readInto

	sprintf(strBuf, "dParMinBound[%i][%i]", i, j);
//	if(! parameters.readInto(oConstants.poSources[i].pdParMinBound[j], strBuf)) {
	if(! parameters.readInto(strTemp, strBuf)) {
	  oConstants.poSources[i].pdParMinBound[j] = 0;
	  printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
	  printf("         using: %f\n", oConstants.poSources[i].pdParMinBound[j]);
	} else {
	  oConstants.poSources[i].pdParMinBound[j] = atof(strTemp.c_str());
	}// if parameters.readInto

	sprintf(strBuf, "dParMaxBound[%i][%i]", i, j);
//	if(! parameters.readInto(oConstants.poSources[i].pdParMaxBound[j], strBuf)) {
	if(! parameters.readInto(strTemp, strBuf)) {
	  oConstants.poSources[i].pdParMaxBound[j] = 0;
	  printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
	  printf("         using: %f\n", oConstants.poSources[i].pdParMaxBound[j]);
	} else {
	  oConstants.poSources[i].pdParMaxBound[j] = atof(strTemp.c_str());
	}// if parameters.readInto

	sprintf(strBuf, "dParStart[%i][%i]", i, j);
//	if(! parameters.readInto(oConstants.poSources[i].pdParStart[j], strBuf)) {
	if(! parameters.readInto(strTemp, strBuf)) {
	  oConstants.poSources[i].pdParStart[j] = 0;
	  printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
	  printf("         using: %f\n", oConstants.poSources[i].pdParStart[j]);
	} else {
	  oConstants.poSources[i].pdParStart[j] = atof(strTemp.c_str());
	}// if parameters.readInto

	sprintf(strBuf, "dParWidthMCMC[%i][%i]", i, j);
//	if(! parameters.readInto(oConstants.poSources[i].pdParWidthMCMC[j], strBuf)) {
	if(! parameters.readInto(strTemp, strBuf)) {
	  oConstants.poSources[i].pdParWidthMCMC[j] = 0;
	  printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
	  printf("         using: %f\n", oConstants.poSources[i].pdParWidthMCMC[j]);
	} else {
	  oConstants.poSources[i].pdParWidthMCMC[j] = atof(strTemp.c_str());
	}// if parameters.readInto

	sprintf(strBuf, "dParWidthFit[%i][%i]", i, j);
//	if(! parameters.readInto(oConstants.poSources[i].pdParWidthFit[j], strBuf)) {
	if(! parameters.readInto(strTemp, strBuf)) {
	  oConstants.poSources[i].pdParWidthFit[j] = 0;
	  printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
	  printf("         using: %f\n", oConstants.poSources[i].pdParWidthFit[j]);
	} else {
	  oConstants.poSources[i].pdParWidthFit[j] = atof(strTemp.c_str());
	}// if parameters.readInto

	sprintf(strBuf, "dParWidthPrior[%i][%i]", i, j);
//	if(! parameters.readInto(oConstants.poSources[i].pdParWidthPrior[j], strBuf)) {
	if(! parameters.readInto(strTemp, strBuf)) {
	  oConstants.poSources[i].pdParWidthPrior[j] = 0;
	  printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
	  printf("         using: %f\n", oConstants.poSources[i].pdParWidthPrior[j]);
	} else {
	  oConstants.poSources[i].pdParWidthPrior[j] = atof(strTemp.c_str());
	}// if parameters.readInto

	sprintf(strBuf, "dParMeanPrior[%i][%i]", i, j);
//	if(! parameters.readInto(oConstants.poSources[i].pdParMeanPrior[j], strBuf)) {
	if(! parameters.readInto(strTemp, strBuf)) {
	  oConstants.poSources[i].pdParMeanPrior[j] = 0;
	  printf("WARNING: parameters.conf did not have: \"%s\"\n", strBuf);
	  printf("         using: %f\n", oConstants.poSources[i].pdParMeanPrior[j]);
	} else {
	  oConstants.poSources[i].pdParMeanPrior[j] = atof(strTemp.c_str());
	}// if parameters.readInto

	sprintf(strBuf, "bParFit[%i][%i]", i, j);
	if(parameters.readInto(strTemp, strBuf)) {
	  if(strcmp(strTemp.c_str(), "true") == 0) {
	    oConstants.poSources[i].oSourceType.nFitTag =
	      oConstants.poSources[i].oSourceType.nFitTag | FITPARAMETER(j);
	  } else {
	    oConstants.poSources[i].oSourceType.nFitTag =
	      oConstants.poSources[i].oSourceType.nFitTag & (~FITPARAMETER(j));
	  } // if strcmp
	} // if parameters.readInto

	sprintf(strBuf, "bParMar[%i][%i]", i, j);
	if(parameters.readInto(strTemp, strBuf)) {
	  if(strcmp(strTemp.c_str(), "true") == 0) {
	    oConstants.poSources[i].oSourceType.nMarTag =
	      oConstants.poSources[i].oSourceType.nMarTag | FITPARAMETER(j);
	  } else {
	    oConstants.poSources[i].oSourceType.nMarTag =
	      oConstants.poSources[i].oSourceType.nMarTag & (~FITPARAMETER(j));
	  } // if strcmp
	} // if parameters.readInto

	sprintf(strBuf, "bParWrap[%i][%i]", i, j);
	if(parameters.readInto(strTemp, strBuf)) {
	  if(strcmp(strTemp.c_str(), "true") == 0) {
	    oConstants.poSources[i].oSourceType.nWrapTag =
	      oConstants.poSources[i].oSourceType.nWrapTag | FITPARAMETER(j);
	  } else {
	    oConstants.poSources[i].oSourceType.nWrapTag =
	      oConstants.poSources[i].oSourceType.nWrapTag & (~FITPARAMETER(j));
	  } // if strcmp
	} // if parameters.readInto

	sprintf(strBuf, "bParReduce[%i][%i]", i, j);
	if(parameters.readInto(strTemp, strBuf)) {
	  if(strcmp(strTemp.c_str(), "true") == 0) {
	    oConstants.poSources[i].oSourceType.nReduceTag =
	      oConstants.poSources[i].oSourceType.nReduceTag | FITPARAMETER(j);
	  } else {
	    oConstants.poSources[i].oSourceType.nReduceTag =
	      oConstants.poSources[i].oSourceType.nReduceTag & (~FITPARAMETER(j));
	  } // if strcmp
	} // if parameters.readInto
      } // for j
    } else { // if nTag != TEMPO2
      // We set the amount of parameters in the plugin. So ignore if it doesn't
      // happen
      oConstants.poSources[i].oSourceType.nParameters = 0;

      // Still read the Fit & Mar parameters for the first one
      sprintf(strBuf, "bParFit[%i][0]", i);
      if(parameters.readInto(strTemp, strBuf)) {
	if(strcmp(strTemp.c_str(), "true") == 0) {
	  oConstants.poSources[i].oSourceType.nFitTag =
	    oConstants.poSources[i].oSourceType.nFitTag | FITPARAMETER(0);
	} else {
	  oConstants.poSources[i].oSourceType.nFitTag =
	    oConstants.poSources[i].oSourceType.nFitTag & (~FITPARAMETER(0));
	} // if strcmp
      } // if parameters.readInto

      sprintf(strBuf, "bParMar[%i][0]", i);
      if(parameters.readInto(strTemp, strBuf)) {
	if(strcmp(strTemp.c_str(), "true") == 0) {
	  oConstants.poSources[i].oSourceType.nMarTag =
	    oConstants.poSources[i].oSourceType.nMarTag | FITPARAMETER(0);
	} else {
	  oConstants.poSources[i].oSourceType.nMarTag =
	    oConstants.poSources[i].oSourceType.nMarTag & (~FITPARAMETER(0));
	} // if strcmp
      } // if parameters.readInto

      if(oConstants.poSources[i].oSourceType.nMarTag & FITPARAMETER(0))
	oConstants.poSources[i].oSourceType.nFitTag =
	  oConstants.poSources[i].oSourceType.nFitTag & (~FITPARAMETER(0));
    } // if nTag != TEMPO2
  } // for i

  // Put the default sources in per pulsar
  for(int k=0; k<oConstants.k; k++) {
    bFound = false;
    for(int i=0; i<oConstants.nSources; i++) {
      if(oConstants.poSources[i].pbScope[k] == true && oConstants.poSources[i].eCorrelation == CID_SinglePulsar && oConstants.poSources[i].oSourceType.eID != SID_Deterministic)
	bFound = true;
    } // for i

    // If there is no non-deterministic source for _only_ this pulsar, add the
    // default one
    if(! bFound) {
      // We need to add a source for pulsar k
#ifndef NO_DEFAULT_SOURCES
      oConstants.nSources++;
      for(int j=0; j<oConstants.k; j++) {
	oConstants.poSources[oConstants.nSources-1] = oConstants.oDefaultPulsarSource;
	// TODO: Check whether the strcpy is needed
	strcpy(oConstants.poSources[oConstants.nSources-1].oSourceType.strID,
	  oConstants.oDefaultPulsarSource.oSourceType.strID);
	for(int k2=0; k2<oConstants.k; k2++)
	  oConstants.poSources[oConstants.nSources-1].pbScope[k2] = k2==k ? true : false;
      } // for j
#endif // NO_DEFAULT_SOURCES
    } // if bFound
  } // for k

  // Now add the default deterministic sources per pulsar
  for(int k=0; k<oConstants.k; k++) {
    // We will not add the default deterministic source to this pulsar if there
    // already is a quadratic spindown term or a tempo2 term
    bFound = false;
    for(int i=0; i<oConstants.nSources; i++) {
      if(oConstants.poSources[i].pbScope[k] == true &&
	  oConstants.poSources[i].oSourceType.eID == SID_Deterministic &&
	  (oConstants.poSources[i].oSourceType.nTag == DETSOURCE_QSD ||
	   oConstants.poSources[i].oSourceType.nTag == DETSOURCE_TEMPO2))
	bFound = true;
    } // for i

    // If there is no deterministic source for this pulsar, add a qsd source
    if(! bFound) {
#ifndef NO_DEFAULT_SOURCES
      // We need to add a source for pulsar k
      oConstants.nSources++;
      for(int j=0; j<oConstants.k; j++) {
	oConstants.poSources[oConstants.nSources-1].oSourceType = oSourceTypes[SID_Deterministic];
	// TODO: Check whether the strcpy is needed
	strcpy(oConstants.poSources[oConstants.nSources-1].oSourceType.strID,
	  oSourceTypes[SID_Deterministic].strID);
	oConstants.poSources[oConstants.nSources-1].eCorrelation = CID_SinglePulsar;

	// Set the detsource type
	for(nDetSourceType = 0; strcmp(oDetSourceTypes[nDetSourceType].strID, "end") != 0; nDetSourceType++) {
          if(strcmp(oDetSourceTypes[nDetSourceType].strID, "qsd") == 0) {
	    oConstants.poSources[oConstants.nSources-1].oSourceType.nTag =
	      oDetSourceTypes[nDetSourceType].nID;
	    oConstants.poSources[oConstants.nSources-1].oSourceType.nParameters =
	      oDetSourceTypes[nDetSourceType].nParameters;
	    oConstants.poSources[oConstants.nSources-1].oSourceType.nWrapTag =
	      FITPARAMETER_NONE;
	    oConstants.poSources[oConstants.nSources-1].oSourceType.nReduceTag =
	      FITPARAMETER_NONE;
	    if(oDetSourceTypes[nDetSourceType].bLinear) {
	      oConstants.poSources[oConstants.nSources-1].oSourceType.nFitTag =
		FITPARAMETER_NONE;
	      oConstants.poSources[oConstants.nSources-1].oSourceType.nMarTag =
		FITPARAMETER_ALL;
	    } else {
	      oConstants.poSources[oConstants.nSources-1].oSourceType.nFitTag =
		FITPARAMETER_ALL;
	      oConstants.poSources[oConstants.nSources-1].oSourceType.nMarTag =
		FITPARAMETER_NONE;
	    } // if bLinear
	    break;
	  } // if strcmp
	} // for nDetSourceType

	// Set the scope
	for(int k2=0; k2<oConstants.k; k2++)
	  oConstants.poSources[oConstants.nSources-1].pbScope[k2] = k2==k ? true : false;

	// Now set the Boundaries and values
	// TODO: what are reasonable values?
	for(int p=0; p<oConstants.poSources[oConstants.nSources-1].oSourceType.nParameters; p++) {
	  oConstants.poSources[oConstants.nSources-1].pdPar[p] = 1E-13;
	  oConstants.poSources[oConstants.nSources-1].pdParMinBound[p] = -10E-13;
	  oConstants.poSources[oConstants.nSources-1].pdParMaxBound[p] = 10E-13;
	  oConstants.poSources[oConstants.nSources-1].pdParStart[p] = 1E-13;
	  oConstants.poSources[oConstants.nSources-1].pdParWidthMCMC[p] = 0.1E-13;
	  oConstants.poSources[oConstants.nSources-1].pdParWidthFit[p] = 0.1E-13;
	  oConstants.poSources[oConstants.nSources-1].pdParWidthPrior[p] = 0E-13;
	} // for p
      } // for j
#endif // NO_DEFAULT_SOURCES
    } // if bFound
  } // for k
  
  SetNumberOfParametersFromSources(oConstants, oData);

  PrintSuccess();
#endif
} // ReadSourceConstants


/* This function first checks how many parameters there are and sets the
 * appropriate variables in the oConstants struct. Then it checks whether we are
 * marginalising or fitting these parameters (mutally exclusive)
 * */
void SetNumberOfParametersFromSources(SConstantsType &oConstants, SDataType &oData) {
  int nMarParameter, nDetSourceType;
  oConstants.nParameters=0;
  oConstants.nFitParameters=0;
  oConstants.nMarParameters=0;
  for(int s=0; s<oConstants.nSources; s++) {
    oConstants.nParameters += oConstants.poSources[s].oSourceType.nParameters;
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p) && oConstants.poSources[s].oSourceType.eID == SID_Deterministic) {
	// Check what the DetSource ID is
	nDetSourceType = 0;
	while(strcmp(oDetSourceTypes[nDetSourceType].strID, "end") != 0) {
	  if(oConstants.poSources[s].oSourceType.nTag == oDetSourceTypes[nDetSourceType].nID)
	    break;
	  nDetSourceType++;
	} // while nDetSourceType

	// Check whether this parameter can be marginalised
	if(oDetSourceTypes[nDetSourceType].bLinear) {
	  oConstants.nMarParameters++;
	  // If this parameter is marginalized, then we will not fit
	  oConstants.poSources[s].oSourceType.nFitTag =
	    oConstants.poSources[s].oSourceType.nFitTag & (~FITPARAMETER(p));
	} else {
	  // This source is _not_ linear!!!
	  printf("WARNING: Source %i is not linear. Set to 'fit'\n", s);
	  oConstants.poSources[s].oSourceType.nMarTag =
	    oConstants.poSources[s].oSourceType.nMarTag & (~FITPARAMETER(p));
	  oConstants.poSources[s].oSourceType.nFitTag =
	    oConstants.poSources[s].oSourceType.nFitTag | FITPARAMETER(p);
	} // if bLinear
      } // if nMarTag
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p))
	oConstants.nFitParameters++;
    } // for p
    if(s==0) {
      oConstants.poSources[s].nFirstParIndex=0;
    } else {
      oConstants.poSources[s].nFirstParIndex =
	oConstants.poSources[s-1].nFirstParIndex +
	oConstants.poSources[s-1].oSourceType.nParameters;
       	oConstants.poSources[s].nFirstParIndex;
    } // if
  } // for s

  // In order to incorporate guadratic spindown, al pulsars must have a different
  // set of qsd parameters. 3 parameters*pulsars to realize the quadratic
  // Now this has been generalised: there are nMarParameter parameters
  oData.vdLinearParameters.Initialize(oConstants.nMarParameters);
  nMarParameter = 0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p) && oConstants.poSources[s].oSourceType.eID == SID_Deterministic) {
	oData.vdLinearParameters[nMarParameter] = oConstants.poSources[s].pdPar[p];
	nMarParameter++;
      } // if nMarTag
    } // for p
  } // for s
} // SetNumberOfParametersFromSources


/* This function calculates the residuals for _all_ deterministic sources. Thus
 * also the sources that are marginalised. This should not matter, because those
 * residuals are removed anyway. There is slight overhead though, because we
 * shouldn't need to calculate those resiudals
 * */
void ResidualsFromDetSources(SConstantsType &oConstants, SParametersType &oParameters, SDataType &oData, CVector *pvdReturnResiduals) {
  int nIndex, nDetParameter, nIndexParameter;

  if(! pvdReturnResiduals) return;
  if(! pvdReturnResiduals->Defined()) return;

  for(nIndex=0; nIndex<pvdReturnResiduals->m_pnDimSize[0]; nIndex++)
    (*pvdReturnResiduals)[nIndex] = 0;

  nIndex=0; nDetParameter=0; nIndexParameter=0;
  for(int a=0; a<oConstants.k; a++) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
      nDetParameter=0;
      nIndexParameter=0;
      for(int s=0; s<oConstants.nSources; s++) {
	switch(oConstants.poSources[s].oSourceType.eID) {
	  case SID_Deterministic:
//	    if(oConstants.poSources[s].pbScope[a]) {
	    if(SourceWorksOnResidual(oConstants, s, a, i)) {
	      (*pvdReturnResiduals)[nIndex] += ResidualFromDetSource(
		  &oParameters.pdPar[oConstants.poSources[s].nFirstParIndex],
		  oConstants.poSources[s], a, i, oData, oConstants);
//	    } // if pbScope
	    } // if SourceWorksOnResidual
	    nIndexParameter += oConstants.poSources[s].oSourceType.nParameters;
	    nDetParameter += oConstants.poSources[s].oSourceType.nParameters;
	    break;
	  default:
	    nIndexParameter += oConstants.poSources[s].oSourceType.nParameters;
	    break;
	} // switch eID
      } // for s
      nIndex++;
    } // for i
  } // for a
  return;
} // ResidualsFromDetSources


/* This function returns the residuals belonging to some deterministic source,
 * linear or non-linear.
 * */
double ResidualFromDetSource(double *pdPar, SSource &oSource, int nPulsar, int nObsIndex, SDataType &oData, SConstantsType &oConstants) {
  double dReturnValue=0;
  switch(oSource.oSourceType.nTag) {
    case DETSOURCE_QSD:
      dReturnValue += pdPar[0]*ParameterDerivative(nPulsar, nObsIndex, oSource.oSourceType.nTag, 0, oData, oConstants);
      dReturnValue += pdPar[1]*ParameterDerivative(nPulsar, nObsIndex, oSource.oSourceType.nTag, 1, oData, oConstants);
      dReturnValue += pdPar[2]*ParameterDerivative(nPulsar, nObsIndex, oSource.oSourceType.nTag, 2, oData, oConstants);
      break;
    case DETSOURCE_TEMPO2:
      for(int nPar=0; nPar<oSource.oSourceType.nParameters; nPar++) {
	dReturnValue += pdPar[nPar]*ParameterDerivative(nPulsar, nObsIndex, oSource.oSourceType.nTag, nPar, oData, oConstants);
      } // for nPar
      break;
    case DETSOURCE_SINE:
      dReturnValue += SineResidual(pdPar, nPulsar, nObsIndex, oConstants);
      break;
    case DETSOURCE_GWMEM:
      dReturnValue += GWMemResidual(pdPar, nPulsar, nObsIndex, oConstants);
      break;
    case DETSOURCE_BHBINARY:
      dReturnValue += BHBinaryResidual(pdPar, nPulsar, nObsIndex, oConstants);
      break;
    case DETSOURCE_CONSTANT:
    case DETSOURCE_LINEAR:
    case DETSOURCE_QUADRATIC:
    case DETSOURCE_CUBIC:
    case DETSOURCE_QUARTIC:
    default:
      dReturnValue += pdPar[0]*ParameterDerivative(nPulsar, nObsIndex, oSource.oSourceType.nTag, 0, oData, oConstants);
      break;
      break;
      // Add more functions here, like the residuals function that
      // Aleksandar is developing
  } // switch nParameterTag
  return dReturnValue;
} // ResidualFromDetSource

/* This function returns the parameter derivative, used in the Block matrix that
 * converts the deterministic (linear) parameters to timing residuals
 *
 * TODO: Make the default option scan for individual detsources and add them
 *       recursively
 * */
double ParameterDerivative(int nPulsar, int nObsIndex, int nParameterTag, int nParameter, SDataType &oData, SConstantsType &oConstants) {
  double dReturnValue=0;
  switch(nParameterTag) { 
    case DETSOURCE_CONSTANT:
      dReturnValue += 1;
      break;
    case DETSOURCE_LINEAR:
      dReturnValue += (oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR);
      break;
    case DETSOURCE_QUADRATIC:
      dReturnValue += (oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR) *
	(oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR);
      break;
    case DETSOURCE_CUBIC:
      dReturnValue += (oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR) *
	(oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR) *
	(oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR);
      break;
    case DETSOURCE_QUARTIC:
      dReturnValue += (oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR) *
	(oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR) *
	(oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR) *
	(oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR);
      break;
    case DETSOURCE_QSD:
      switch(nParameter) {
	case 0:
	  dReturnValue += 1;
	  break;
	case 1:
	  dReturnValue += (oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR);
	  break;
	case 2:
	  dReturnValue += (oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR) *
	    (oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR);
	default:
	  break;
      } // switch nParameter
      break;
    case DETSOURCE_TEMPO2:
//      dReturnValue += oConstants.poPulsars[nPulsar].ppdTempo2ParamterDerivative[nObsIndex][nParameter];
	dReturnValue += double(oConstants.poPulsars[nPulsar].mdTempo2ParameterDerivative[nObsIndex][nParameter]);
      break;
    default:
      break;
  } // switch nParameterTag
  return dReturnValue;
} // ParameterDerivative

/* This function returns the residuals for the sinusoidal deterministic source.
 * This is just an example of how non-linear deterministic sources can be
 * implemented. The sinusoidal deterministic source has 3 parameters: the
 * amplitude, frequency and phase.
 *
 * TODO: The phase is a periodic parameter. There is no support for these kinds
 *       of parameters yet. Might be interesting to make.
 * */
double SineResidual(double *pdPar, int nPulsar, int nObsIndex, SConstantsType &oConstants) {
  return pdPar[0] * sin(pdPar[1]*oConstants.poPulsars[nPulsar].pdTOA[nObsIndex]/SPERYEAR + pdPar[2]);
} // SineResidual

/* This function returns the residuals for the gw-memory effect of BH mergers.
 * It is just a linear function with a variable starting point
 *
 * First rotate the pulsar angles so that the gw-memory source is moving along
 * the z-direction. Then rotate the polarisation angle so that it's along the
 * x-axis.
 * */
double GWMemResidual(double *pdPar, int nPulsar, int nObsIndex, SConstantsType &oConstants) {
  CMatrix mdR1, mdR2, mdR3, mdR;
  CVector vdUnit;
  double dCosTwoPhi, dCosTheta;
  mdR1.Initialize(3,3);
  mdR2.Initialize(3,3);
  mdR3.Initialize(3,3);
  vdUnit.Initialize(3);

  // Rotate along the azimuthal angle
  mdR1[0][0] = cos(pdPar[2]);		mdR1[0][1] = sin(pdPar[2]);	mdR1[0][2] = 0;
  mdR1[1][0] = -sin(pdPar[2]);		mdR1[1][1] = cos(pdPar[2]);	mdR1[1][2] = 0;
  mdR1[2][0] = 0;			mdR1[2][1] = 0;			mdR1[2][2] = 1;

  // Rotate the polar angle to zero (set the direction of gw along the z-axis)
  mdR2[0][0] = cos(M_PI_2-pdPar[3]);	mdR2[0][1] = 0;		mdR2[0][2] = -sin(M_PI_2 - pdPar[3]);
  mdR2[1][0] = 0;			mdR2[1][1] = 1;		mdR2[1][2] = 0;
  mdR2[2][0] = sin(M_PI_2-pdPar[3]);	mdR2[2][1] = 0;		mdR2[2][2] = cos(M_PI_2 - pdPar[3]);

  // Rotate so that the gwmem polarisation is in the x-direction
  mdR3[0][0] = cos(pdPar[4]);	mdR3[0][1] = sin(pdPar[4]);	mdR3[0][2] = 0;
  mdR3[1][0] = -sin(pdPar[4]);	mdR3[1][1] = cos(pdPar[4]);	mdR3[1][2] = 0;
  mdR3[2][0] = 0;		mdR3[2][1] = 0;			mdR3[2][2] = 1;

  // The unit vector that represents the pulsar position
  vdUnit[0] = cos(oConstants.poPulsars[nPulsar].dPhi)*sin(oConstants.poPulsars[nPulsar].dTheta);
  vdUnit[1] = sin(oConstants.poPulsars[nPulsar].dPhi)*sin(oConstants.poPulsars[nPulsar].dTheta);
  vdUnit[2] = cos(oConstants.poPulsars[nPulsar].dTheta);

  // Calculate the total rotation matrix
  mdR = mdR3 * mdR2 * mdR1;

  // Calculate the new position of the pulsar
  vdUnit = mdR * vdUnit;

  // Calculate the necessary geometric factors
  dCosTheta = double(vdUnit[2]);
  dCosTwoPhi = (fabs(dCosTheta) == 1.0) ?
    0:
    2*(double(vdUnit[0])*double(vdUnit[0])/(1-double(vdUnit[2])*double(vdUnit[2]))) - 1;

  return oConstants.poPulsars[nPulsar].pdTOA[nObsIndex] > pdPar[0] ?
    0.5*dCosTwoPhi*(1+dCosTheta)*
    pdPar[1]*(oConstants.poPulsars[nPulsar].pdTOA[nObsIndex] - pdPar[0]) :
    0;
} // GWMemResidual

/* This function returns the residuals for a massive black-hole binary
 *
 * Function directly copied from paper by KJ Lee et al. (arXiv:1103.0115)
 *
 * Parameters used in the function:
 *   pdPar[0]: h_0, the GW amplitude in metric perturbation. typically ~10^{-15}
 *   pdPar[1]: w_g, the GW angular frequency at the observer [Hz]. typ.~10^{-7} Hz
 *   pdPar[2]: lambda, the ecliptic longitude of the GW source [0, 2pi]
 *   pdPar[3]: beta, the ecliptic lattitude of the GW source [0, pi]
 *   pdPar[4]: i, the inclination of the GW binary orbit [0, pi]
 *   pdPar[5]: phi, direction of the binary ascending node (polarisation) [0,pi]
 *   pdPar[6]: lambda_p, the ecliptic longitude of the pulsar [0, 2pi]
 *   pdPar[7]: beta_p, the ecliptic lattitude of the pulsar [0, pi]
 *   pdPar[8]: Dpsr, the distance to pulsar 1 [kpc]. typ. 1 kpc
 *      ..
 *   pdPar[MAX_PARAMETERS_PER_SOURCE]: the distance to pulsar MPPS-7 [kpc].
 * */
double BHBinaryResidual(double *pdPar, int nPulsar, int nObsIndex, SConstantsType &oConstants) { // BHBinaryResidual
  const double dC = 299792458;		// Speed of light in m/s
  const double dPc = 3.08568025e16;	// Parsec in meters

  double dH, dW_g, dDpsr, dLambda_p, dBeta_p, dLambda, dBeta, dI, dPhi,
	 dDeltaPhi, dB1, dB2, dCosTheta, dT,
	 dCoefficient, dTemp, dResidual;

  // Set the parameters given to this funciton
  dH = pdPar[0];			// GW amplitude
  dW_g = pdPar[1];			// GW angular frequency in Hz
  dLambda = pdPar[2];			// Ecliptic longitude in radians (source)
  dBeta = 0.5*M_PI-pdPar[3];		// Ecliptic lattitude in radians (source)
  dI = pdPar[4];			// Inclination of binary in radians
  dPhi = pdPar[5];			// The polarisation of the GW in radians

  // Although the location of the pulsar in principle is a variable, it is
  // probably better determined otherwise. Use the already known information?
  dLambda_p = pdPar[6];			// Ecliptic longitude in radians (pulsar)
  dBeta_p = 0.5*M_PI-pdPar[7];		// Ecliptic lattitude in radians (pulsar)
  dLambda_p = oConstants.poPulsars[nPulsar].dPhi;
  dBeta_p = 0.5*M_PI-oConstants.poPulsars[nPulsar].dTheta;
  //dBeta_p = oConstants.poPulsars[nPulsar].dTheta;

//  printf("%e   %e   %e   %e\n", dLambda_p, dBeta_p, dLambda, dBeta);

  // We are only working with one particular pulsar at the moment. Use that one
  // Check whether the parameter is in the parameter range
  if(8 + nPulsar < MAX_PARAMETERS_PER_SOURCE)
    dDpsr = pdPar[8+nPulsar] * dPc * 1000;	// Distance to pulsar in meters
  else {
    printf("WARNING: distance parameter out of range (MAX_PARAMETERS_PER_SOURCE);\n");
    return 0.0;
  } //

  // Some auxiliary values based on the parameters
  dCosTheta = cos(dBeta) * cos(dBeta_p) * cos(dLambda - dLambda_p) +
    sin(dBeta) * sin(dBeta_p);
  dDeltaPhi = dW_g * dDpsr * (1 - dCosTheta) / dC;
  dB1 = (1+sin(dBeta)*sin(dBeta))*cos(dBeta_p)*cos(dBeta_p)*cos(2*(dLambda - dLambda_p))
      - sin(2*dBeta)*sin(2*dBeta_p)*cos(dLambda-dLambda_p)
      + (2-3*cos(dBeta_p)*cos(dBeta_p))*cos(dBeta)*cos(dBeta);
  dB2 = 2*cos(dBeta)*sin(2*dBeta_p)*sin(dLambda-dLambda_p)
      - 2*sin(dBeta)*cos(dBeta_p)*cos(dBeta_p)*sin(2*(dLambda-dLambda_p));
  dT = oConstants.poPulsars[nPulsar].pdTOA[nObsIndex];

  // Now to calculate the residual
  dCoefficient = dH * sin(0.5*dDeltaPhi) / (2*dW_g*(1-dCosTheta));
  dTemp = ( dB1*cos(2*dPhi) + dB2*sin(2*dPhi) )
          *cos(dW_g*dT - 0.5*dDeltaPhi) * (1 + cos(dI)*cos(dI))
	+ 2*(dB2*cos(2*dPhi)-dB1*sin(2*dPhi))
	  *sin(dW_g*dT - 0.5*dDeltaPhi) * cos(dI);
  dResidual = dCoefficient * dTemp;
//  printf("%e  %e  %e  %e  %e\n", dCosTheta, dDeltaPhi, dB1, dB2, dResidual);
  return dResidual;
} // BHBinaryResidual

/* Sets the matrix that converts the linear deterministic parameters to timing
 * residuals. This matrix is also used for analytically marginalising those
 * parameters.
 *
 * There are 3 matrices: Ave, Gen, and Block.  Ave can be used to average the
 * residuals to zero, Gen is used to generate the data (nowadays equal to
 * 'block') Block will contain all deterministic linear parameters
 *
 * TODO: Nowadays it isn't really a block matrix anymore. So perhaps change the
 * name? Make it: designmatrix
 * */
void SetBlockMatrix(SDataType &oData, SConstantsType &oConstants) {
  int nIndex=0;
  int nMarParameter=0;

  oData.mdAveBlock.Initialize(oConstants.n, oConstants.k);
  oData.mdGenBlock.Initialize(oConstants.n, oConstants.nMarParameters);
  oData.mdBlock.Initialize(oConstants.n, oConstants.nMarParameters);

  // For use in LogLikelihood
  oData.mdCXiInv.Initialize(oData.mdBlock.m_pnDimSize[1], oData.mdBlock.m_pnDimSize[1]);
  oData.vdChi.Initialize(oData.mdBlock.m_pnDimSize[1]);

  // Now construct the matrix M which consists of functions of measure times
  nIndex=0;
  for(int a=0; a<oConstants.k; a++) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
      for(int b=0; b<oConstants.k; b++) {
        for(int j=0; j<1; j++) {
          oData.mdAveBlock[nIndex][j+b] = 0;
          if(a==b) // We need to fill this spot
            oData.mdAveBlock[nIndex][j+b] = pow(oConstants.poPulsars[a].pdTOA[i], j);
        } // for j
      } // for b
      nIndex++;
    } // for i
  } // for a

  nIndex=0;
  for(int a=0; a<oConstants.k; a++) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
      for(int d=0; d<oConstants.nMarParameters; d++) {
	oData.mdBlock[nIndex][d] = 0;
      } // for d
      nIndex++;
    } // for i
  } // for a

  // Now construct the matrix M which consists of all the deterministic sources
  nIndex=0; nMarParameter=0;
  for(int a=0; a<oConstants.k; a++) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
      nMarParameter=0;
      for(int s=0; s<oConstants.nSources; s++) {
	switch(oConstants.poSources[s].oSourceType.eID) {
	  case SID_Deterministic:
	    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
	      if(oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p)) {
//		if(oConstants.poSources[s].pbScope[a]) {
		if(SourceWorksOnResidual(oConstants, s, a, i)) {
		  oData.mdBlock[nIndex][nMarParameter] +=
		      ParameterDerivative(a, i, oConstants.poSources[s].oSourceType.nTag, p, oData, oConstants);
//		} // if pbScope
		} // if SourceWorksOnResidual
		nMarParameter++;
	      } // if nMarTag
	    } // for p
	    break;
	  default:
	    break;
	} // switch eID
      } // for s
      nIndex++;
    } // for i
  } // for a

  oData.mdBlockTrans = oData.mdBlock[LO_TRANSPOSE];
  oData.mdGenBlock = oData.mdBlock;
} // SetBlockMatrix


/* This function prepares new objects of the struct SConstantsType and SDataType
 *
 * In these new objects, only the pulsars and sources of specifically specified
 * pulsars of the original objects are present (provided by pbSelection). This
 * allows one to fit for just a specific set of pulsars, or just for one pulsar,
 * which can be handy for preprocessing the data before doing a MCMC on the
 * whole dataset.
 * */
void CopySpecificPulsarsToObjects(
    SConstantsType &oConstants, SDataType &oData,
    SConstantsType &oNewConstants, SDataType &oNewData,
    bool *pbSelection) {
  int *pnNewIndices;
  bool bCopiedSource, bTemp;
  pnNewIndices = new int[oConstants.k];

  // First copy all kinds of general stuff (plotpoints, datadirs, etc)
  oNewConstants = oConstants;

  // Find out how many and which new pulsars we are processing
  oNewConstants.k = 0;
  oNewConstants.n = 0;
  for(int p=0; p<oConstants.k; p++) {
    pnNewIndices[p] = -1;
    if(pbSelection[p]) {
      // We add this pulsar
      pnNewIndices[p] = oNewConstants.k;
      oNewConstants.poPulsars[oNewConstants.k] = oConstants.poPulsars[p];
      oNewConstants.n += oNewConstants.poPulsars[oNewConstants.k].nObservations;
      oNewConstants.k++;
    } // if pbSelection
  } // for p

  // Now only select the sources we need
  oNewConstants.nSources = 0;
  oNewConstants.nParameters = 0;
  oNewConstants.nFitParameters = 0;
  oNewConstants.nMarParameters = 0;
  for(int s=0; s<oConstants.nSources; s++) {
    bCopiedSource = false;
    for(int p=0; p<oConstants.k; p++) {
      bTemp = false;
      for(int i=0; i<oConstants.poPulsars[p].nObservations; i++) {
	if(SourceWorksOnResidual(oConstants, s, p, i))
	  bTemp = true;
      } // for i
//      if(oConstants.poSources[s].pbScope[p] && pnNewIndices[p] != -1) {
      if(bTemp && pnNewIndices[p] != -1) {
	// This source includes this pulsar, and it is a selected pulsar
	if(! bCopiedSource) {
	  oNewConstants.poSources[oNewConstants.nSources] = oConstants.poSources[s];
	  oNewConstants.poSources[oNewConstants.nSources].nFirstParIndex =
	    oNewConstants.poSources[oNewConstants.nSources].oSourceType.nParameters;
	  oNewConstants.nParameters +=
	    oNewConstants.poSources[oNewConstants.nSources].oSourceType.nParameters;
	  oNewConstants.nSources++;
	  for(int a=0; a<oNewConstants.k; a++)
  	    oNewConstants.poSources[oNewConstants.nSources-1].pbScope[a] = false;
	  bCopiedSource = true;
	} // if bCopiedSource

	// Set the right pulsar it applies to
	if(oConstants.poSources[s].pbScope[p])
	  oNewConstants.poSources[oNewConstants.nSources-1].pbScope[pnNewIndices[p]] = true;
	strcpy(oNewConstants.poSources[oNewConstants.nSources-1].strScopeTag, oConstants.poSources[s].strScopeTag);
      } // if bTemp
    } // for p
  } // for s


  // Set the correct number of parameters
  SetNumberOfParametersFromSources(oNewConstants, oNewData);

  // Set the right data structures
  SetDataVectorFromPulsars(oNewConstants, oNewData);
  SetBlockMatrix(oNewData, oNewConstants);
  SetGeometricPart(oNewData, oNewConstants);
  oNewData.mdC.Initialize(oNewConstants.n, oNewConstants.n);
//  SetCoherenceMatrix(oNewData, oNewConstants);

  delete[] pnNewIndices;
} // CopySpecificPulsarsToObjects

void SetDataVectorFromPulsars(SConstantsType &oConstants, SDataType &oData) {
  int nIndex=0;
  oData.vdData.Initialize(oConstants.n);
  oData.vdDataErr.Initialize(oConstants.n);

  for(int p=0; p<oConstants.k; p++) {
    for(int i=0; i<oConstants.poPulsars[p].nObservations; i++) {
      oData.vdData[nIndex] = oConstants.poPulsars[p].pdResiduals[i];
      oData.vdDataErr[nIndex] = oConstants.poPulsars[p].pdDeltaResiduals[i];
      nIndex++;
    } // for i
  } // for p
} // SetDataVectorFromPulsars


/* Generates a set of pulsar angles, evenly distributed over the sky. Whether
 * those pulsars are on the northern or southern hemisphere is not taken into
 * account.
 *
 * TODO: It is very unrealistic to generate them evenly distributed. Usually
 * only one part of the hemisphere is visible
 */
void GeneratePulsarAngles(SDataType &oData, SConstantsType &oConstants) {
// Generate random angles for each pulsar
  CVector vdTemp1, vdTemp, vdPhi, vdTheta;

  vdPhi.Initialize(oConstants.k);
  vdTheta.Initialize(oConstants.k);
  vdTemp.Initialize(oConstants.k);

  vdPhi.Randomize();
  vdPhi *= 2*M_PI;
  vdTemp.Randomize();
  vdTemp1 = vdTemp*2-1;
  for(int a=0; a<oConstants.k; a++) vdTemp[a] = 0.5*M_PI;
  vdTheta = vdTemp - ArcSin(vdTemp1);

  for(int a=0; a<oConstants.k; a++) {
    oConstants.poPulsars[a].dPhi = double(vdPhi[a]);
    oConstants.poPulsars[a].dTheta = double(vdTheta[a]);
  } // for a
  return;
} // GeneratePulsarAngles



/* Sets the geometrical part of the coherence matrix: the pulsar correlation
 * coefficients (this is only correlation due to GWB).
 *
 * Basically this is: 1.5xlog(x)x -x/4 +1/2 + delta(x)/2
 * */
void SetGeometricPart(SDataType &oData, SConstantsType &oConstants) {
// Generate random angles for each pulsar and calculate the angles between
// the pulsars in matrix geometric. (Geometric dependant correlation matrix)

  double dPhi, dTheta, dKsi, dPsi, dRz;
  oData.mdGeometric.Initialize(oConstants.k, oConstants.k);

  for(int a=0; a<oConstants.k; a++) {
    for(int b=0; b<oConstants.k; b++) {
      dPhi = oConstants.poPulsars[a].dPhi;
      dTheta = oConstants.poPulsars[a].dTheta;
      dKsi = oConstants.poPulsars[b].dPhi;
      dPsi = oConstants.poPulsars[b].dTheta;

#if 0
      double dkp[3], ddec=M_PI_2 - dTheta, dra = dPhi;
      double dkg[3];
      dkp[0] = sin(dTheta)*cos(dPhi);
      dkp[1] = sin(dTheta)*sin(dPhi);
      dkp[2] = cos(dTheta);
      dkg[0] = sin(dPsi)*cos(dKsi);
      dkg[1] = sin(dPsi)*sin(dKsi);
      dkg[2] = cos(dPsi);
      dRz = sin(dTheta)*sin(dPsi)*cos(dPhi-dKsi)+cos(dTheta)*cos(dPsi);

      fprintf(stderr, "kp: %lf %lf %lf\n", dkp[0], dkp[1], dkp[2]);
      fprintf(stderr, "dec, ra: %lf %lf\n", ddec, dra);
      fprintf(stderr, "cos(dec), cos(ra), sin(ra): %lf %lf %lf\n",
	cos(ddec), cos(dra), sin(dra));
      fprintf(stderr, "compare: %e  %e\n", dRz,
	  dkp[0]*dkg[0]+dkp[1]*dkg[1]+dkp[2]*dkg[2]);
#endif


      // Rotate vector r1 (phi, theta) to the northpole with two rotations
      // and then the rz of vector r2 (ksi, psi) is the cos of the angle
      // between the two pulsars
      // Or, just take the inner-product of the two (duh). That's dRz
      dRz = sin(dTheta)*sin(dPsi)*cos(dPhi-dKsi)+cos(dTheta)*cos(dPsi);
      // Compute the geometric dependent elements of the correlation matrix
      if(dRz>=1) {
        oData.mdGeometric[a][b] = 1.0/2.0;
//        oData.mdGeometric[a][b] = 1.0/3.0;
//        oData.mdGeometric[a][b] = 0;
      } else {
        oData.mdGeometric[a][b] = 0.75*(1-dRz)*log(0.5*(1-dRz))-(1.0/8.0)*(1-dRz)+1.0/2.0;
//        oData.mdGeometric[a][b] = 0.125 * (dRz + 3);
//        oData.mdGeometric[a][b] = 0.5*(1-dRz)*log(0.5*(1-dRz))-(1.0/12.0)*(1-dRz)+1.0/3.0;
//        oData.mdGeometric[a][b] = 0;
      } // if
      if(a==b) { // This part is for the self-correlation of a pulsar (the noise-term)
        oData.mdGeometric[a][b] += 1.0/2.0;
//        oData.mdGeometric[a][b] += 1.0/3.0;
      } // if a==b
    } // for b
  } // for a
} // SetGeometricPart


/* This function calculates the full Coherence matrix of the pulsar timing
 * array data. Both GWB and pulsar noise is taken into account.
 *
 * This version should be used to generate the data: it does not get its values
 * for the parameters from a SParameters struct, but from the SConstants struct.
 * */
void SetCoherenceMatrix(SDataType &oData, SConstantsType &oConstants, bool bTest) {
  double dTau, dTemp, dCor;
  bool bCheck;
  int nIndex1, nIndex2;
  CLinearObject *cloTemp;
  CVector vdParameters, vdParametersNull;
  vdParameters.Initialize(MAX_PARAMETERS_PER_SOURCE);

  PrintStatus("Setting coherence values...");

#if 0
  time_t tTimeNow, tTimePrevious;
  double dTimeSeconds;
  tTimePrevious = clock();
#endif

  if(! oData.mdC.Defined()) oData.mdC.Initialize(oConstants.n, oConstants.n);
  if(! oData.mdInvC.Defined()) oData.mdInvC.Initialize(oConstants.n, oConstants.n);

  // TODO: This is a lot of overhead for each calculation
  for(int i=0; i<oConstants.n; i++)
    for(int j=0; j<oConstants.n; j++)
      oData.mdC[i][j] = 0;

  // Fill the coherence matrix for all sources
  // TODO: Add error-bar handling
  // Continue here
  for(int s=0; s<oConstants.nSources; s++) {
    switch(oConstants.poSources[s].oSourceType.eID) {
      case SID_White:
      case SID_PowerLaw:
      case SID_Exponential:
      case SID_Lorentzian:
      case SID_PowerLaw_OneParameter:
      case SID_Errorbars:
      case SID_DMVar:
      case SID_DipoleGWB:
      case SID_NonStationary:
        for(int i=0; i<oConstants.poSources[s].oSourceType.nParameters; i++)
          vdParameters[i] = oConstants.poSources[s].pdPar[i];
        
        switch(oConstants.poSources[s].oSourceType.eID) {
          case SID_Errorbars:
            break;
          case SID_DipoleGWB:
            /* We need to do this separately, because the example 'Chiara'
             * source is actually a power-law, though with extra angular
             * parameters. "SpectrumIntegral" does not know about this*/
            SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
            break;
          case SID_DMVar:
            SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
            break;
          case SID_NonStationary:
            vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
            SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
            break;
          default:
            SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
            break;
        } // switch

	nIndex1 = nIndex2 = 0;
	for(int a=0; a<oConstants.k; a++) {
#if 0 // See what sources apply
	if(SourceWorksOnPulsar(oConstants, s, a))
	    printf("Source: %i\n", s);
#endif
	  for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	    for(int b=0; b<oConstants.k; b++) {
	      for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
		if(SourceWorksOnResidual(oConstants, s, a, i) &&
		    SourceWorksOnResidual(oConstants, s, b, j)) {
		  switch(oConstants.poSources[s].eCorrelation) {
		    case CID_SinglePulsar:
		      dCor = 1;
		      break;
		    case CID_GR:
		      dCor = double(oData.mdGeometric[a][b]);
		      break;
		    case CID_All:
		      dCor = a == b ? 1 : 0;
		      break;
		    case CID_Uniform:
		      dCor = 1;
		      break;
		    case CID_Metric_Breathing:
		    case CID_Metric_Shear:
		    case CID_Metric_Longitudinal:
		    default:
		      dCor = 1;
		      break;
		  } // switch eCorrelation
		  switch(oConstants.poSources[s].oSourceType.eID) {
		    case SID_Errorbars:
		      if(a==b && i==j)
			dTemp = exp(oConstants.poSources[s].pdPar[0])*
			  exp(oConstants.poSources[s].pdPar[0])*
			  oConstants.poPulsars[a].pdDeltaResiduals[i]*
			  oConstants.poPulsars[a].pdDeltaResiduals[i];
		      else
			dTemp = 0;
		      break;
		    case SID_White:
		      if(a==b && i==j) {
			dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				    oConstants.poPulsars[b].pdTOA[j]);
			dTemp = SpectrumIntegral(
			    oConstants.poSources[s].oSourceType.eID,
			    oConstants.poPulsars[a].pdTOA[i],
			    oConstants.poPulsars[b].pdTOA[j],
			    vdParametersNull); 
		      } else {
			dTemp = 0;
		      } // if a==b && i==j
		      break;
            case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[a].pdFreq[i]*oConstants.poPulsars[a].pdFreq[i]*
                  oConstants.poPulsars[b].pdFreq[j]*oConstants.poPulsars[b].pdFreq[j]);
              break;
            case SID_DipoleGWB:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull); 
              dTemp *= AnisotropicGWBCorrelation(a, b,
                      oConstants.poSources[s].oSourceType.nTag,
                      oConstants, vdParameters);
              break;
		    default:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
			  oConstants.poSources[s].oSourceType.eID,
			  oConstants.poPulsars[a].pdTOA[i],
			  oConstants.poPulsars[b].pdTOA[j],
			  vdParametersNull); 
		      break;
		  } // switch
		  oData.mdC[nIndex1][nIndex2] += dCor*dTemp;
		} // if SourceWorksOnResidual
		nIndex2++;
	      } // for j
	    } // for b
	    nIndex1++;
	    nIndex2 = 0;
	  } // for i
	} // for a

	break;
      case SID_Deterministic: // Do nothing
      case SID_Nothing:
      default:
	break;
    } // switch eID
  } // for s

  // Print diagonal elements for posdef check
#if 0
  printf("\n");
  for(int i=0; i<oData.mdC.m_pnDimSize[0]; i++) {
//    printf("%e  ", double(oData.mdC[i][i]));
    printf("%i: %e\n", i, double(oData.mdC[i][i]));
  } // for i
  printf("\n");
#endif

#if 0
  printf("\n");
  for(int i=0; i<oData.mdC.m_pnDimSize[0]; i++) {
      printf("%i: ", i);
    for(int j=0; j<oData.mdC.m_pnDimSize[0]; j++) {
      printf("%e  ", double(oData.mdC[i][j]));
    } // for j
    printf("\n");
  } // for i
  printf("\n");
#endif

  PrintSuccess();

#if 0
  for(int i=0; i<oData.mdC.m_pnDimSize[0]; i++) {
    if(double(oData.mdC[i][i]) == 0.0) {
      printf("i = %i\n", i);
    }
  } // for i
#endif

#if 0
  tTimeNow = clock();
  dTimeSeconds = double(tTimeNow - tTimePrevious)/CLOCKS_PER_SEC;
  printf("Time needed to fill matrix: %f\n", dTimeSeconds);
#endif

  if(bTest) {
    bCheck = true;
    double dLogDet;
//    CMatrix mdTemp;
    PrintStatus("Checking cholesky on mdC...");
//    try { LogDetChol(oData.mdC); } catch (ELinearError err) { bCheck = false;} 
//    try { mdTemp = oData.mdC.InverseChol(&dLogDet); } catch (ELinearError err) { bCheck = false;} 
    try { oData.mdCholeskyC = oData.mdC.Cholesky(); } catch (ELinearError err) { bCheck = false;} 
//    mdTemp = oData.mdC;
    if(bCheck) PrintSuccess(); else PrintFailed();

#if 0
    PrintStatus("Calculating eigenvectors/values...");
    cloTemp = oData.mdC.Eigen("V");
    oData.vdLambdaC = *((CVector *)&cloTemp[0]);
    oData.mdEVC = *((CMatrix *)&cloTemp[1]);
    delete[] cloTemp;
    PrintSuccess();

    bCheck = true;
    PrintStatus("Checking whether matrix is posdef by eigenvalues...");
    for(int i=0; i<oData.vdLambdaC.m_pnDimSize[0]; i++) {
      if(double(oData.vdLambdaC[i]) < 0)
	bCheck = false;
    } // for i
    if(bCheck)
      PrintSuccess();
    else
      PrintFailed();
  //  printf("Min ev: %f   Max ev: %f\n", Min(oData.vdLambdaC), Max(oData.vdLambdaC));
#endif
  } // bTest

  // Now do the new covariance stuff
#if 0
  {
    CMatrix mdC, mdNewC, mdM, mdMTrans, mdTemp1, mdTemp2, mdB, mdI;
    mdC = oData.mdC;
    mdM = oData.mdBlock;
    mdMTrans = oData.mdBlock[LO_TRANSPOSE];
    FILE *pFile;

    mdI = mdC;
    for(int i=0; i<mdC.m_pnDimSize[0]; i++) {
      for(int j=0; j<mdC.m_pnDimSize[1]; j++) {
	mdI[i][j] = ((i == j) ? 1.0 : 0.0);
      } // for j
    } // for i

    mdTemp1 = mdMTrans * mdM;
    mdTemp1.InvertChol();
    mdTemp2 = mdM * mdTemp1;
    mdTemp1 = mdTemp2 * mdMTrans;
    mdB = mdI - mdTemp1;

    mdTemp1 = mdC * mdB;
    mdNewC = mdB * mdTemp1;

    if(! (pFile=fopen("covariance.txt", "w+")) ) {
      printf("Unable to open file\n");
      return;
    }

    for(int i=0; i<mdC.m_pnDimSize[0]; i++) {
      for(int j=0; j<mdC.m_pnDimSize[1]; j++) {
	fprintf(pFile, "%.18e  %.18e  %.18e \n",
	    oConstants.poPulsars[0].pdTOA[i],
	    oConstants.poPulsars[0].pdTOA[j],
	    double(mdNewC[i][j]));
//	    double(mdC[i][j]));
      } // for j
      fprintf(pFile, "\n");
    } // for i
    fclose(pFile);
  }
#endif
} // SetCoherenceMatrix


/* This function calculates the full Coherence matrix of the pulsar timing
 * array data. Both GWB and pulsar noise is taken into account.
 *
 * This version should be used in the analysis routines: it does gets its values
 * for the parameters from a SParameters struct.
 * */
void SetCoherenceMatrix(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dTau, dTemp, dCor;
  int nIndex1, nIndex2;
  CVector vdParameters, vdParametersNull;
  vdParameters.Initialize(3);

  if(! oData.mdC.Defined()) oData.mdC.Initialize(oConstants.n, oConstants.n);
  if(! oData.mdInvC.Defined()) oData.mdInvC.Initialize(oConstants.n, oConstants.n);

  // TODO: This is a lot of overhead for each calculation
  for(int i=0; i<oConstants.n; i++)
    for(int j=0; j<oConstants.n; j++)
      oData.mdC[i][j] = 0;

  // Fill the coherence matrix for all sources
  // TODO: Add error-bar handling
  // Continue here
  for(int s=0; s<oConstants.nSources; s++) {
    switch(oConstants.poSources[s].oSourceType.eID) {
      case SID_White:
      case SID_PowerLaw:
      case SID_Exponential:
      case SID_Lorentzian:
      case SID_PowerLaw_OneParameter:
      case SID_Errorbars:
      case SID_DMVar:
      case SID_DipoleGWB:
      case SID_NonStationary:
	for(int i=0; i<oConstants.poSources[s].oSourceType.nParameters; i++)
	  vdParameters[i] = oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+i];
	switch(oConstants.poSources[s].oSourceType.eID) {
	  case SID_Errorbars:
	    break;
      case SID_DipoleGWB:
        /* We need to do this separately, because the example 'Chiara'
         * source is actually a power-law, though with extra angular
         * parameters. "SpectrumIntegral" does not know about this*/
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
        break;
      case SID_DMVar:
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
        break;
	  case SID_NonStationary:
	    vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
        break;
	  default:
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	} // switch

	nIndex1 = nIndex2 = 0;
	for(int a=0; a<oConstants.k; a++) {
	  for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	    for(int b=0; b<oConstants.k; b++) {
	      for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
		if(SourceWorksOnResidual(oConstants, s, a, i) &&
		    SourceWorksOnResidual(oConstants, s, b, j)) {
		  switch(oConstants.poSources[s].eCorrelation) {
		    case CID_SinglePulsar:
		      dCor = 1;
		      break;
		    case CID_GR:
		      dCor = double(oData.mdGeometric[a][b]);
		      break;
		    case CID_All:
		      dCor = a == b ? 1 : 0;
		      break;
		    case CID_Uniform:
		      dCor = 1;
		      break;
		    case CID_Metric_Breathing:
		    case CID_Metric_Shear:
		    case CID_Metric_Longitudinal:
		    default:
		      dCor = 1;
		      break;
		  } // switch eCorrelation
		  switch(oConstants.poSources[s].oSourceType.eID) {
		    case SID_Errorbars:
		      if(a==b && i==j)
			dTemp = exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  oConstants.poPulsars[a].pdDeltaResiduals[i]*
			  oConstants.poPulsars[a].pdDeltaResiduals[i];
		      else
			dTemp = 0;
		      break;
		    case SID_White:
		      if(a==b && i==j) {
			dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				    oConstants.poPulsars[b].pdTOA[j]);
			dTemp = SpectrumIntegral(
			    oConstants.poSources[s].oSourceType.eID,
			    oConstants.poPulsars[a].pdTOA[i],
			    oConstants.poPulsars[b].pdTOA[j],
			    vdParametersNull); 
		      } else {
			dTemp = 0;
		      } // if a==b && i==j
		      break;
            case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[a].pdFreq[i]*oConstants.poPulsars[a].pdFreq[i]*
                  oConstants.poPulsars[b].pdFreq[j]*oConstants.poPulsars[b].pdFreq[j]);
              break;
            case SID_DipoleGWB:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull); 
              dTemp *= AnisotropicGWBCorrelation(a, b,
                      oConstants.poSources[s].oSourceType.nTag,
                      oConstants, vdParameters);
              break;
		    default:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
			  oConstants.poSources[s].oSourceType.eID,
			  oConstants.poPulsars[a].pdTOA[i],
			  oConstants.poPulsars[b].pdTOA[j],
			  vdParametersNull); 
		      break;
		  } // switch
		  oData.mdC[nIndex1][nIndex2] += dCor*dTemp;
		} // if SourceWorksOnResidual
		nIndex2++;
	      } // for j
	    } // for b
	    nIndex1++;
	    nIndex2 = 0;
	  } // for i
	} // for a

	break;
      case SID_Deterministic: // Do nothing
      case SID_Nothing:
      default:
	break;
    } // switch eID
  } // for s

#if 0
  printf("\n");
  for(int i=0; i<oData.mdC.m_pnDimSize[0]; i++) {
    printf("%e  ", double(oData.mdC[i][i]));
  } // for i
  printf("\n");
#endif
} // SetCoherenceMatrix

/* This function calculates the Coherence matrix of the non-clock sources of
 * the pulsar timing array data. Both GWB and pulsar noise is taken into
 * account.
 *
 * This version should be used in the analysis routines: it does gets its values
 * for the parameters from a SParameters struct.
 * */
void SetCoherenceMatrixNonClockSources(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dTau, dTemp, dCor;
  int nIndex1, nIndex2;
  CVector vdParameters, vdParametersNull;
  vdParameters.Initialize(3);

  if(! oData.mdC.Defined()) oData.mdC.Initialize(oConstants.n, oConstants.n);
  if(! oData.mdInvC.Defined()) oData.mdInvC.Initialize(oConstants.n, oConstants.n);

  // TODO: This is a lot of overhead for each calculation
  for(int i=0; i<oConstants.n; i++)
    for(int j=0; j<oConstants.n; j++)
      oData.mdC[i][j] = 0;

  // Fill the coherence matrix for all sources
  for(int s=0; s<oConstants.nSources; s++) {
    if(oConstants.poSources[s].eCorrelation == CID_Uniform)
      continue;
    switch(oConstants.poSources[s].oSourceType.eID) {
      case SID_White:
      case SID_PowerLaw:
      case SID_Exponential:
      case SID_Lorentzian:
      case SID_PowerLaw_OneParameter:
      case SID_Errorbars:
      case SID_DMVar:
      case SID_DipoleGWB:
      case SID_NonStationary:
	for(int i=0; i<oConstants.poSources[s].oSourceType.nParameters; i++)
	  vdParameters[i] = oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+i];
	switch(oConstants.poSources[s].oSourceType.eID) {
	  case SID_Errorbars:
	    break;
      case SID_DMVar:
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
        break;
      case SID_DipoleGWB:
        /* We need to do this separately, because the example 'Chiara'
         * source is actually a power-law, though with extra angular
         * parameters. "SpectrumIntegral" does not know about this*/
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
        break;
	  case SID_NonStationary:
	    vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
        break;
	  default:
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	} // switch

	nIndex1 = nIndex2 = 0;
	for(int a=0; a<oConstants.k; a++) {
	  for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	    for(int b=0; b<oConstants.k; b++) {
	      for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
		if(SourceWorksOnResidual(oConstants, s, a, i) &&
		    SourceWorksOnResidual(oConstants, s, b, j)) {
		  switch(oConstants.poSources[s].eCorrelation) {
		    case CID_SinglePulsar:
		      dCor = 1;
		      break;
		    case CID_GR:
		      dCor = double(oData.mdGeometric[a][b]);
		      break;
		    case CID_All:
		      dCor = a == b ? 1 : 0;
		      break;
		    case CID_Uniform:
		      dCor = 1;
		      break;
		    case CID_Metric_Breathing:
		    case CID_Metric_Shear:
		    case CID_Metric_Longitudinal:
		    default:
		      dCor = 1;
		      break;
		  } // switch eCorrelation
		  switch(oConstants.poSources[s].oSourceType.eID) {
		    case SID_Errorbars:
		      if(a==b && i==j)
			dTemp = exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  oConstants.poPulsars[a].pdDeltaResiduals[i]*
			  oConstants.poPulsars[a].pdDeltaResiduals[i];
		      else
			dTemp = 0;
		      break;
		    case SID_White:
		      if(a==b && i==j) {
			dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				    oConstants.poPulsars[b].pdTOA[j]);
			dTemp = SpectrumIntegral(
			    oConstants.poSources[s].oSourceType.eID,
			    oConstants.poPulsars[a].pdTOA[i],
			    oConstants.poPulsars[b].pdTOA[j],
			    vdParametersNull); 
		      } else {
			dTemp = 0;
		      } // if a==b && i==j
		      break;
            case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[a].pdFreq[i]*oConstants.poPulsars[a].pdFreq[i]*
                  oConstants.poPulsars[b].pdFreq[j]*oConstants.poPulsars[b].pdFreq[j]);
              break;
            case SID_DipoleGWB:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull); 
              dTemp *= AnisotropicGWBCorrelation(a, b,
                      oConstants.poSources[s].oSourceType.nTag,
                      oConstants, vdParameters);
              break;
		    default:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
			  oConstants.poSources[s].oSourceType.eID,
			  oConstants.poPulsars[a].pdTOA[i],
			  oConstants.poPulsars[b].pdTOA[j],
			  vdParametersNull); 
		      break;
		  } // switch
		  oData.mdC[nIndex1][nIndex2] += dCor*dTemp;
		} // if SourceWorksOnResidual
		nIndex2++;
	      } // for j
	    } // for b
	    nIndex1++;
	    nIndex2 = 0;
	  } // for i
	} // for a

	break;
      case SID_Deterministic: // Do nothing
      case SID_Nothing:
      default:
	break;
    } // switch eID
  } // for s

#if 0
  printf("\n");
  for(int i=0; i<oData.mdC.m_pnDimSize[0]; i++) {
    printf("%e  ", double(oData.mdC[i][i]));
  } // for i
  printf("\n");
#endif
} // SetCoherenceMatrixNonClockSources

/* This function calculates the Coherence matrix of the clock sources of
 * the pulsar timing array data. These are all sources uniformly correlated.
 *
 * This version should be used in the analysis routines: it does gets its values
 * for the parameters from a SParameters struct.
 * */
void SetCoherenceMatrixClockSources(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dTau, dTemp, dCor;
  int nIndex1, nIndex2;
  CVector vdParameters, vdParametersNull;
  vdParameters.Initialize(3);

  if(! oData.mdC.Defined()) oData.mdC.Initialize(oConstants.n, oConstants.n);
  if(! oData.mdInvC.Defined()) oData.mdInvC.Initialize(oConstants.n, oConstants.n);

  // TODO: This is a lot of overhead for each calculation
  for(int i=0; i<oConstants.n; i++)
    for(int j=0; j<oConstants.n; j++)
      oData.mdC[i][j] = 0;

  // Fill the coherence matrix for all sources
  for(int s=0; s<oConstants.nSources; s++) {
    if(oConstants.poSources[s].eCorrelation != CID_Uniform)
      continue;
    switch(oConstants.poSources[s].oSourceType.eID) {
      case SID_White:
      case SID_PowerLaw:
      case SID_Exponential:
      case SID_Lorentzian:
      case SID_PowerLaw_OneParameter:
      case SID_Errorbars:
      case SID_DMVar:
      case SID_DipoleGWB:
      case SID_NonStationary:
	for(int i=0; i<oConstants.poSources[s].oSourceType.nParameters; i++)
	  vdParameters[i] = oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+i];
	switch(oConstants.poSources[s].oSourceType.eID) {
	  case SID_Errorbars:
	    break;
      case SID_DMVar:
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
        break;
      case SID_DipoleGWB:
        /* We need to do this separately, because the example 'Chiara'
         * source is actually a power-law, though with extra angular
         * parameters. "SpectrumIntegral" does not know about this*/
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
        break;
	  case SID_NonStationary:
	    vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	  default:
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	} // switch

	nIndex1 = nIndex2 = 0;
	for(int a=0; a<oConstants.k; a++) {
	  for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	    for(int b=0; b<oConstants.k; b++) {
	      for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
		if(SourceWorksOnResidual(oConstants, s, a, i) &&
		    SourceWorksOnResidual(oConstants, s, b, j)) {
		  switch(oConstants.poSources[s].eCorrelation) {
		    case CID_SinglePulsar:
		      dCor = 1;
		      break;
		    case CID_GR:
		      dCor = double(oData.mdGeometric[a][b]);
		      break;
		    case CID_All:
		      dCor = a == b ? 1 : 0;
		      break;
		    case CID_Uniform:
		      dCor = 1;
		      break;
		    case CID_Metric_Breathing:
		    case CID_Metric_Shear:
		    case CID_Metric_Longitudinal:
		    default:
		      dCor = 1;
		      break;
		  } // switch eCorrelation
		  switch(oConstants.poSources[s].oSourceType.eID) {
		    case SID_Errorbars:
		      if(a==b && i==j)
			dTemp = exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  oConstants.poPulsars[a].pdDeltaResiduals[i]*
			  oConstants.poPulsars[a].pdDeltaResiduals[i];
		      else
			dTemp = 0;
		      break;
		    case SID_White:
		      if(a==b && i==j) {
			dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				    oConstants.poPulsars[b].pdTOA[j]);
			dTemp = SpectrumIntegral(
			    oConstants.poSources[s].oSourceType.eID,
			    oConstants.poPulsars[a].pdTOA[i],
			    oConstants.poPulsars[b].pdTOA[j],
			    vdParametersNull); 
		      } else {
			dTemp = 0;
		      } // if a==b && i==j
		      break;
            case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[a].pdFreq[i]*oConstants.poPulsars[a].pdFreq[i]*
                  oConstants.poPulsars[b].pdFreq[j]*oConstants.poPulsars[b].pdFreq[j]);
              break;
            case SID_DipoleGWB:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull); 
              dTemp *= AnisotropicGWBCorrelation(a, b,
                      oConstants.poSources[s].oSourceType.nTag,
                      oConstants, vdParameters);
              break;
		    default:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
			  oConstants.poSources[s].oSourceType.eID,
			  oConstants.poPulsars[a].pdTOA[i],
			  oConstants.poPulsars[b].pdTOA[j],
			  vdParametersNull); 
		      break;
		  } // switch
		  oData.mdC[nIndex1][nIndex2] += dCor*dTemp;
		} // if SourceWorksOnResidual
		nIndex2++;
	      } // for j
	    } // for b
	    nIndex1++;
	    nIndex2 = 0;
	  } // for i
	} // for a

	break;
      case SID_Deterministic: // Do nothing
      case SID_Nothing:
      default:
	break;
    } // switch eID
  } // for s

#if 0
  printf("\n");
  for(int i=0; i<oData.mdC.m_pnDimSize[0]; i++) {
    printf("%e  ", double(oData.mdC[i][i]));
  } // for i
  printf("\n");
#endif
} // SetCoherenceMatrixClockSources




/* This function calculates the Coherence matrix of the non-reduced sources of
 * the pulsar timing array data. Both GWB and pulsar noise is taken into
 * account.
 *
 * This version should be used in the analysis routines: it does gets its values
 * for the parameters from a SParameters struct.
 * */
void SetCoherenceMatrixNonReduceSources(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dTau, dTemp, dCor;
  int nIndex1, nIndex2;
  CVector vdParameters, vdParametersNull;
  vdParameters.Initialize(3);

  if(! oData.mdC.Defined()) oData.mdC.Initialize(oConstants.n, oConstants.n);
  if(! oData.mdInvC.Defined()) oData.mdInvC.Initialize(oConstants.n, oConstants.n);

  // TODO: This is a lot of overhead for each calculation
  for(int i=0; i<oConstants.n; i++)
    for(int j=0; j<oConstants.n; j++)
      oData.mdC[i][j] = 0;

  // Fill the coherence matrix for all sources
  for(int s=0; s<oConstants.nSources; s++) {
    if(oConstants.poSources[s].oSourceType.nReduceTag & FITPARAMETER(0))
      continue;
    switch(oConstants.poSources[s].oSourceType.eID) {
      case SID_White:
      case SID_PowerLaw:
      case SID_Exponential:
      case SID_Lorentzian:
      case SID_PowerLaw_OneParameter:
      case SID_Errorbars:
      case SID_DMVar:
      case SID_DipoleGWB:
      case SID_NonStationary:
	for(int i=0; i<oConstants.poSources[s].oSourceType.nParameters; i++)
	  vdParameters[i] = oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+i];
	switch(oConstants.poSources[s].oSourceType.eID) {
	  case SID_Errorbars:
	    break;
      case SID_DMVar:
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
        break;
      case SID_DipoleGWB:
        /* We need to do this separately, because the example 'Chiara'
         * source is actually a power-law, though with extra angular
         * parameters. "SpectrumIntegral" does not know about this*/
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
        break;
	  case SID_NonStationary:
	    vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	  default:
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	} // switch

	nIndex1 = nIndex2 = 0;
	for(int a=0; a<oConstants.k; a++) {
	  for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	    for(int b=0; b<oConstants.k; b++) {
	      for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
		if(SourceWorksOnResidual(oConstants, s, a, i) &&
		    SourceWorksOnResidual(oConstants, s, b, j)) {
		  switch(oConstants.poSources[s].eCorrelation) {
		    case CID_SinglePulsar:
		      dCor = 1;
		      break;
		    case CID_GR:
		      dCor = double(oData.mdGeometric[a][b]);
		      break;
		    case CID_All:
		      dCor = a == b ? 1 : 0;
		      break;
		    case CID_Uniform:
		      dCor = 1;
		      break;
		    case CID_Metric_Breathing:
		    case CID_Metric_Shear:
		    case CID_Metric_Longitudinal:
		    default:
		      dCor = 1;
		      break;
		  } // switch eCorrelation
		  switch(oConstants.poSources[s].oSourceType.eID) {
		    case SID_Errorbars:
		      if(a==b && i==j)
			dTemp = exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  oConstants.poPulsars[a].pdDeltaResiduals[i]*
			  oConstants.poPulsars[a].pdDeltaResiduals[i];
		      else
			dTemp = 0;
		      break;
		    case SID_White:
		      if(a==b && i==j) {
			dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				    oConstants.poPulsars[b].pdTOA[j]);
			dTemp = SpectrumIntegral(
			    oConstants.poSources[s].oSourceType.eID,
			    oConstants.poPulsars[a].pdTOA[i],
			    oConstants.poPulsars[b].pdTOA[j],
			    vdParametersNull); 
		      } else {
			dTemp = 0;
		      } // if a==b && i==j
		      break;
            case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[a].pdFreq[i]*oConstants.poPulsars[a].pdFreq[i]*
                  oConstants.poPulsars[b].pdFreq[j]*oConstants.poPulsars[b].pdFreq[j]);
              break;
            case SID_DipoleGWB:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull); 
              dTemp *= AnisotropicGWBCorrelation(a, b,
                      oConstants.poSources[s].oSourceType.nTag,
                      oConstants, vdParameters);
              break;
		    default:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
			  oConstants.poSources[s].oSourceType.eID,
			  oConstants.poPulsars[a].pdTOA[i],
			  oConstants.poPulsars[b].pdTOA[j],
			  vdParametersNull); 
		      break;
		  } // switch
		  oData.mdC[nIndex1][nIndex2] += dCor*dTemp;
		} // if SourceWorksOnResidual
		nIndex2++;
	      } // for j
	    } // for b
	    nIndex1++;
	    nIndex2 = 0;
	  } // for i
	} // for a

	break;
      case SID_Deterministic: // Do nothing
      case SID_Nothing:
      default:
	break;
    } // switch eID
  } // for s

#if 0
  printf("\n");
  for(int i=0; i<oData.mdC.m_pnDimSize[0]; i++) {
    printf("%e  ", double(oData.mdC[i][i]));
  } // for i
  printf("\n");
#endif
} // SetCoherenceMatrixNonReduceSources

/* This function calculates the Coherence matrix of the reduceable sources of
 * the pulsar timing array data. Both GWB and pulsar noise is taken into
 * account.
 *
 * This version should be used in the analysis routines: it does gets its values
 * for the parameters from a SParameters struct.
 * */
void SetCoherenceMatrixReduceSources(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dTau, dTemp, dCor;
  int nIndex1, nIndex2;
  CVector vdParameters, vdParametersNull;
  vdParameters.Initialize(3);

  if(! oData.mdC.Defined()) oData.mdC.Initialize(oConstants.n, oConstants.n);
  if(! oData.mdInvC.Defined()) oData.mdInvC.Initialize(oConstants.n, oConstants.n);

  // TODO: This is a lot of overhead for each calculation
  for(int i=0; i<oConstants.n; i++)
    for(int j=0; j<oConstants.n; j++)
      oData.mdC[i][j] = 0;

  // Fill the coherence matrix for all sources
  for(int s=0; s<oConstants.nSources; s++) {
    if(! (oConstants.poSources[s].oSourceType.nReduceTag & FITPARAMETER(0)))
      continue;
    switch(oConstants.poSources[s].oSourceType.eID) {
      case SID_White:
      case SID_PowerLaw:
      case SID_Exponential:
      case SID_Lorentzian:
      case SID_PowerLaw_OneParameter:
      case SID_Errorbars:
      case SID_DMVar:
      case SID_DipoleGWB:
      case SID_NonStationary:
	for(int i=0; i<oConstants.poSources[s].oSourceType.nParameters; i++)
	  vdParameters[i] = oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+i];
	switch(oConstants.poSources[s].oSourceType.eID) {
	  case SID_Errorbars:
	    break;
      case SID_DMVar:
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
        break;
      case SID_DipoleGWB:
        /* We need to do this separately, because the example 'Chiara'
         * source is actually a power-law, though with extra angular
         * parameters. "SpectrumIntegral" does not know about this*/
        SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
        break;
	  case SID_NonStationary:
	    vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	  default:
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	} // switch

	nIndex1 = nIndex2 = 0;
	for(int a=0; a<oConstants.k; a++) {
	  for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	    for(int b=0; b<oConstants.k; b++) {
	      for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
		if(SourceWorksOnResidual(oConstants, s, a, i) &&
		    SourceWorksOnResidual(oConstants, s, b, j)) {
		  switch(oConstants.poSources[s].eCorrelation) {
		    case CID_SinglePulsar:
		      dCor = 1;
		      break;
		    case CID_GR:
		      dCor = double(oData.mdGeometric[a][b]);
		      break;
		    case CID_All:
		      dCor = a == b ? 1 : 0;
		      break;
		    case CID_Uniform:
		      dCor = 1;
		      break;
		    case CID_Metric_Breathing:
		    case CID_Metric_Shear:
		    case CID_Metric_Longitudinal:
		    default:
		      dCor = 1;
		      break;
		  } // switch eCorrelation
		  switch(oConstants.poSources[s].oSourceType.eID) {
		    case SID_Errorbars:
		      if(a==b && i==j)
			dTemp = exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  oConstants.poPulsars[a].pdDeltaResiduals[i]*
			  oConstants.poPulsars[a].pdDeltaResiduals[i];
		      else
			dTemp = 0;
		      break;
		    case SID_White:
		      if(a==b && i==j) {
			dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				    oConstants.poPulsars[b].pdTOA[j]);
			dTemp = SpectrumIntegral(
			    oConstants.poSources[s].oSourceType.eID,
			    oConstants.poPulsars[a].pdTOA[i],
			    oConstants.poPulsars[b].pdTOA[j],
			    vdParametersNull); 
		      } else {
			dTemp = 0;
		      } // if a==b && i==j
		      break;
            case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[a].pdFreq[i]*oConstants.poPulsars[a].pdFreq[i]*
                  oConstants.poPulsars[b].pdFreq[j]*oConstants.poPulsars[b].pdFreq[j]);
              break;
            case SID_DipoleGWB:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull); 
              dTemp *= AnisotropicGWBCorrelation(a, b,
                      oConstants.poSources[s].oSourceType.nTag,
                      oConstants, vdParameters);
              break;
		    default:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
			  oConstants.poSources[s].oSourceType.eID,
			  oConstants.poPulsars[a].pdTOA[i],
			  oConstants.poPulsars[b].pdTOA[j],
			  vdParametersNull); 
		      break;
		  } // switch
		  oData.mdC[nIndex1][nIndex2] += dCor*dTemp;
		} // if SourceWorksOnResidual
		nIndex2++;
	      } // for j
	    } // for b
	    nIndex1++;
	    nIndex2 = 0;
	  } // for i
	} // for a

	break;
      case SID_Deterministic: // Do nothing
      case SID_Nothing:
      default:
	break;
    } // switch eID
  } // for s

#if 0
  printf("\n");
  for(int i=0; i<oData.mdC.m_pnDimSize[0]; i++) {
    printf("%e  ", double(oData.mdC[i][i]));
  } // for i
  printf("\n");
#endif
} // SetCoherenceMatrixReduceSources


/* This function calculates the Coherence matrix of a single pulsar, without any
 * correlated signals, in order to speed up the calculations
 *
 * This version should be used in the analysis routines: it does gets its values
 * for the parameters from a SParameters struct.
 * */
void SetCoherenceMatrixSinglePulsarNoCorrSig(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, int nPulsar) {
  double dTau, dTemp, dCor;
  int nIndex1, nIndex2;
  CVector vdParameters, vdParametersNull;
  vdParameters.Initialize(3);

  // TODO: This is a lot of overhead for each calculation
  for(int i=0; i<oConstants.poPulsars[nPulsar].nObservations; ++i)
    for(int j=0; j<oConstants.poPulsars[nPulsar].nObservations; ++j)
      oConstants.poPulsars[nPulsar].mdC[i][j] = 0;

  // Fill the covariance matrix for all sources working on this pulsar
  for(int s=0; s<oConstants.nSources; s++) {
    if(SourceWorksOnPulsar(oConstants, s, nPulsar) &&
	oConstants.poSources[s].eCorrelation == CID_SinglePulsar) {
      switch(oConstants.poSources[s].oSourceType.eID) {
	case SID_White:
	case SID_PowerLaw:
	case SID_Exponential:
	case SID_Lorentzian:
	case SID_PowerLaw_OneParameter:
	case SID_Errorbars:
    case SID_DMVar:
    case SID_DipoleGWB:
	case SID_NonStationary:
	  for(int i=0; i<oConstants.poSources[s].oSourceType.nParameters; i++)
	    vdParameters[i] = oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+i];
	  switch(oConstants.poSources[s].oSourceType.eID) {
	    case SID_Errorbars:
	      break;
        case SID_DMVar:
          SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
          break;
        case SID_DipoleGWB:
          /* We need to do this separately, because the example 'Chiara'
           * source is actually a power-law, though with extra angular
           * parameters. "SpectrumIntegral" does not know about this*/
          SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
          break;
	    case SID_NonStationary:
	      vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
	      SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	      break;
	    default:
	      SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	      break;
	  } // switch

	  for(int i=0; i<oConstants.poPulsars[nPulsar].nObservations; ++i) {
	    for(int j=0; j<oConstants.poPulsars[nPulsar].nObservations; ++j) {
	      switch(oConstants.poSources[s].oSourceType.eID) {
		case SID_Errorbars:
		  if(i==j)
		    dTemp = exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
		      exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
		      oConstants.poPulsars[nPulsar].pdDeltaResiduals[i]*
		      oConstants.poPulsars[nPulsar].pdDeltaResiduals[i];
		  else
		    dTemp = 0;
		  break;
		case SID_White:
		  if(i==j) {
		    dTau = fabs(oConstants.poPulsars[nPulsar].pdTOA[i] -
			      oConstants.poPulsars[nPulsar].pdTOA[j]);
		    dTemp = SpectrumIntegral(
			oConstants.poSources[s].oSourceType.eID,
			oConstants.poPulsars[nPulsar].pdTOA[i],
			oConstants.poPulsars[nPulsar].pdTOA[j],
			vdParametersNull); 
		  } else {
		    dTemp = 0;
		  } // if i==j
		  break;
        case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[nPulsar].pdTOA[i] -
				  oConstants.poPulsars[nPulsar].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[nPulsar].pdTOA[i],
                      oConstants.poPulsars[nPulsar].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[nPulsar].pdFreq[i]*oConstants.poPulsars[nPulsar].pdFreq[i]*
                  oConstants.poPulsars[nPulsar].pdFreq[j]*oConstants.poPulsars[nPulsar].pdFreq[j]);
          break;
        case SID_DipoleGWB:
		      dTau = fabs(oConstants.poPulsars[nPulsar].pdTOA[i] -
				  oConstants.poPulsars[nPulsar].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[nPulsar].pdTOA[i],
                      oConstants.poPulsars[nPulsar].pdTOA[j],
                      vdParametersNull); 
              dTemp *= AnisotropicGWBCorrelation(nPulsar, nPulsar,
                      oConstants.poSources[s].oSourceType.nTag,
                      oConstants, vdParameters);
          break;
		default:
		  dTau = fabs(oConstants.poPulsars[nPulsar].pdTOA[i] -
			      oConstants.poPulsars[nPulsar].pdTOA[j]);
		  dTemp = SpectrumIntegral(
		      oConstants.poSources[s].oSourceType.eID,
		      oConstants.poPulsars[nPulsar].pdTOA[i],
		      oConstants.poPulsars[nPulsar].pdTOA[j],
		      vdParametersNull); 
		  break;
	      } // switch
	      oConstants.poPulsars[nPulsar].mdC[i][j] += dTemp;
	    } /* for j */
	  } /* for i */
	  break;
	case SID_Deterministic: // Do nothing
	case SID_Nothing:
	default:
	  break;
      } // switch eID
    } /* if SourceWorksOnPulsar */
  } /* for s */
} /* SetCoherenceMatrixSinglePulsarNoCorrSig */

/* This function calculates the Coherence matrix of a single pulsar, without any
 * correlated signals, in order to speed up the calculations
 *
 * This version should be used in the analysis routines: it does gets its values
 * for the parameters from a SParameters struct.
 * */
void SetCoherenceMatrixCorrSig(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants) {
  double dTau, dTemp, dCor;
  int nIndex1, nIndex2;
  CVector vdParameters, vdParametersNull;
  vdParameters.Initialize(3);

  if(! oData.mdCCor.Defined()) oData.mdCCor.Initialize(oConstants.n, oConstants.n);

  // TODO: This is a lot of overhead for each calculation
  for(int i=0; i<oConstants.n; i++)
    for(int j=0; j<oConstants.n; j++)
      oData.mdCCor[i][j] = 0;

  // Fill the coherence matrix for all sources
  for(int s=0; s<oConstants.nSources; s++) {
    if(oConstants.poSources[s].eCorrelation == CID_SinglePulsar)
      continue;

    switch(oConstants.poSources[s].oSourceType.eID) {
      case SID_White:
      case SID_PowerLaw:
      case SID_Exponential:
      case SID_Lorentzian:
      case SID_PowerLaw_OneParameter:
      case SID_Errorbars:
      case SID_DMVar:
      case SID_DipoleGWB:
      case SID_NonStationary:
	for(int i=0; i<oConstants.poSources[s].oSourceType.nParameters; i++)
	  vdParameters[i] = oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+i];
	switch(oConstants.poSources[s].oSourceType.eID) {
	  case SID_Errorbars:
	    break;
          case SID_DMVar:
            SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
            break;
          case SID_DipoleGWB:
            /* We need to do this separately, because the example 'Chiara'
             * source is actually a power-law, though with extra angular
             * parameters. "SpectrumIntegral" does not know about this*/
            SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
            break;
	  case SID_NonStationary:
	    vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	  default:
  	    SpectrumIntegral(oConstants.poSources[s].oSourceType.eID, 0, 0, vdParameters);
	    break;
	} // switch

	nIndex1 = nIndex2 = 0;
	for(int a=0; a<oConstants.k; a++) {
	  for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	    for(int b=0; b<oConstants.k; b++) {
	      switch(oConstants.poSources[s].eCorrelation) {
		case CID_GR:
		  dCor = double(oData.mdGeometric[a][b]);
		  break;
		case CID_All:
		  dCor = a == b ? 1 : 0;
		  break;
		case CID_Uniform:
		case CID_SinglePulsar:
		  dCor = 1;
		  break;
		case CID_Metric_Breathing:
		case CID_Metric_Shear:
		case CID_Metric_Longitudinal:
		default:
		  dCor = 1;
		  break;
	      } // switch eCorrelation

	      for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
		if(SourceWorksOnResidual(oConstants, s, a, i) &&
		    SourceWorksOnResidual(oConstants, s, b, j)) {
		  switch(oConstants.poSources[s].oSourceType.eID) {
		    case SID_Errorbars:
		      if(a==b && i==j)
			dTemp = exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  exp(oParameters.pdPar[oConstants.poSources[s].nFirstParIndex])*
			  oConstants.poPulsars[a].pdDeltaResiduals[i]*
			  oConstants.poPulsars[a].pdDeltaResiduals[i];
		      else
			dTemp = 0;
		      break;
		    case SID_White:
		      if(a==b && i==j) {
			dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				    oConstants.poPulsars[b].pdTOA[j]);
			dTemp = SpectrumIntegral(
			    oConstants.poSources[s].oSourceType.eID,
			    oConstants.poPulsars[a].pdTOA[i],
			    oConstants.poPulsars[b].pdTOA[j],
			    vdParametersNull); 
		      } else {
			dTemp = 0;
		      } // if a==b && i==j
		      break;
            case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[a].pdFreq[i]*oConstants.poPulsars[a].pdFreq[i]*
                  oConstants.poPulsars[b].pdFreq[j]*oConstants.poPulsars[b].pdFreq[j]);
              break;
            case SID_DipoleGWB:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull); 
              dTemp *= AnisotropicGWBCorrelation(a, b,
                      oConstants.poSources[s].oSourceType.nTag,
                      oConstants, vdParameters);
              break;
		    default:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
			  oConstants.poSources[s].oSourceType.eID,
			  oConstants.poPulsars[a].pdTOA[i],
			  oConstants.poPulsars[b].pdTOA[j],
			  vdParametersNull); 
		      break;
		  } // switch
		  oData.mdCCor[nIndex1][nIndex2] += dCor*dTemp;
		} // if SourceWorksOnResidual
		nIndex2++;
	      } // for j
	    } // for b
	    nIndex1++;
	    nIndex2 = 0;
	  } // for i
	} // for a
	break;
      case SID_Deterministic: // Do nothing
      case SID_Nothing:
      default:
	break;
    } // switch eID
  } /* for s */
} /* SetCoherenceMatrixCorrSig */


/* This function calculates the covariance matrix of one specific source. For
 * only one pulsar, or for the whole array, based on the nPulsar parameter
 * */
void SetCoherenceMatrixOneSource(SDataType &oData, SParametersType &oParameters, SConstantsType &oConstants, CMatrix &mdC, int nSource, int nPulsar) {
  double dTau, dTemp, dCor;
  int nIndex1, nIndex2;
  CVector vdParameters, vdParametersNull;
  vdParameters.Initialize(3);

  if(nPulsar < 0) {
    // Calculate covariance for whole array
    if(! mdC.Defined())
      mdC.Initialize(oConstants.n, oConstants.n);

    for(int i=0; i<oConstants.n; i++)
      for(int j=0; j<oConstants.n; j++)
	mdC[i][j] = 0;

    switch(oConstants.poSources[nSource].oSourceType.eID) {
      case SID_White:
      case SID_PowerLaw:
      case SID_Exponential:
      case SID_Lorentzian:
      case SID_PowerLaw_OneParameter:
      case SID_Errorbars:
      case SID_DMVar:
      case SID_DipoleGWB:
      case SID_NonStationary:
	for(int i=0; i<oConstants.poSources[nSource].oSourceType.nParameters; i++)
	  vdParameters[i] = oParameters.pdPar[oConstants.poSources[nSource].nFirstParIndex+i];
	switch(oConstants.poSources[nSource].oSourceType.eID) {
	  case SID_Errorbars:
	    break;
          case SID_DMVar:
            SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
            break;
          case SID_DipoleGWB:
            /* We need to do this separately, because the example 'Chiara'
             * source is actually a power-law, though with extra angular
             * parameters. "SpectrumIntegral" does not know about this*/
            SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
            break;
	  case SID_NonStationary:
	    vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
  	    SpectrumIntegral(oConstants.poSources[nSource].oSourceType.eID, 0, 0, vdParameters);
	    break;
	  default:
  	    SpectrumIntegral(oConstants.poSources[nSource].oSourceType.eID, 0, 0, vdParameters);
	    break;
	} // switch

	nIndex1 = nIndex2 = 0;
	for(int a=0; a<oConstants.k; a++) {
	  for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	    for(int b=0; b<oConstants.k; b++) {
	      switch(oConstants.poSources[nSource].eCorrelation) {
		case CID_GR:
		  dCor = double(oData.mdGeometric[a][b]);
		  break;
		case CID_All:
		  dCor = a == b ? 1 : 0;
		  break;
		case CID_Uniform:
		case CID_SinglePulsar:
		  dCor = 1;
		  break;
		case CID_Metric_Breathing:
		case CID_Metric_Shear:
		case CID_Metric_Longitudinal:
		default:
		  dCor = 1;
		  break;
	      } // switch eCorrelation

	      for(int j=0; j<oConstants.poPulsars[b].nObservations; j++) {
		if(SourceWorksOnResidual(oConstants, nSource, a, i) &&
		    SourceWorksOnResidual(oConstants, nSource, b, j)) {
		  switch(oConstants.poSources[nSource].oSourceType.eID) {
		    case SID_Errorbars:
		      if(a==b && i==j)
			dTemp = exp(oParameters.pdPar[oConstants.poSources[nSource].nFirstParIndex])*
			  exp(oParameters.pdPar[oConstants.poSources[nSource].nFirstParIndex])*
			  oConstants.poPulsars[a].pdDeltaResiduals[i]*
			  oConstants.poPulsars[a].pdDeltaResiduals[i];
		      else
			dTemp = 0;
		      break;
		    case SID_White:
		      if(a==b && i==j) {
			dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				    oConstants.poPulsars[b].pdTOA[j]);
			dTemp = SpectrumIntegral(
			    oConstants.poSources[nSource].oSourceType.eID,
			    oConstants.poPulsars[a].pdTOA[i],
			    oConstants.poPulsars[b].pdTOA[j],
			    vdParametersNull); 
		      } else {
			dTemp = 0;
		      } // if a==b && i==j
		      break;
            case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[a].pdFreq[i]*oConstants.poPulsars[a].pdFreq[i]*
                  oConstants.poPulsars[b].pdFreq[j]*oConstants.poPulsars[b].pdFreq[j]);
              break;
            case SID_DipoleGWB:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[a].pdTOA[i],
                      oConstants.poPulsars[b].pdTOA[j],
                      vdParametersNull); 
              dTemp *= AnisotropicGWBCorrelation(a, b,
                      oConstants.poSources[nSource].oSourceType.nTag,
                      oConstants, vdParameters);
              break;
		    default:
		      dTau = fabs(oConstants.poPulsars[a].pdTOA[i] -
				  oConstants.poPulsars[b].pdTOA[j]);
		      dTemp = SpectrumIntegral(
			  oConstants.poSources[nSource].oSourceType.eID,
			  oConstants.poPulsars[a].pdTOA[i],
			  oConstants.poPulsars[b].pdTOA[j],
			  vdParametersNull); 
		      break;
		  } // switch
		  mdC[nIndex1][nIndex2] += dCor*dTemp;
		} // if SourceWorksOnResidual
		nIndex2++;
	      } // for j
	    } // for b
	    nIndex1++;
	    nIndex2 = 0;
	  } // for i
	} // for a
	break;
      case SID_Deterministic: // Do nothing
      case SID_Nothing:
      default:
	break;
    } // switch eID
  } else {
    // Calculate covariance for one pulsar
    mdC.Initialize(oConstants.poPulsars[nPulsar].nObservations,
	oConstants.poPulsars[nPulsar].nObservations);

    for(int i=0; i<mdC.m_pnDimSize[0]; i++)
      for(int j=0; j<mdC.m_pnDimSize[0]; j++)
	mdC[i][j] = 0;

    switch(oConstants.poSources[nSource].oSourceType.eID) {
      case SID_White:
      case SID_PowerLaw:
      case SID_Exponential:
      case SID_Lorentzian:
      case SID_PowerLaw_OneParameter:
      case SID_Errorbars:
      case SID_DMVar:
      case SID_DipoleGWB:
      case SID_NonStationary:
	for(int i=0; i<oConstants.poSources[nSource].oSourceType.nParameters; i++)
	  vdParameters[i] = oParameters.pdPar[oConstants.poSources[nSource].nFirstParIndex+i];
	switch(oConstants.poSources[nSource].oSourceType.eID) {
	  case SID_Errorbars:
	    break;
          case SID_DMVar:
            SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters, true);
            break;
          case SID_DipoleGWB:
            /* We need to do this separately, because the example 'Chiara'
             * source is actually a power-law, though with extra angular
             * parameters. "SpectrumIntegral" does not know about this*/
            SpectrumIntegral(SID_PowerLaw, 0, 0, vdParameters);
            break;
	  case SID_NonStationary:
	    vdParameters[2] = oConstants.poPulsars[0].pdTOA[0];
  	    SpectrumIntegral(oConstants.poSources[nSource].oSourceType.eID, 0, 0, vdParameters);
	    break;
	  default:
  	    SpectrumIntegral(oConstants.poSources[nSource].oSourceType.eID, 0, 0, vdParameters);
	    break;
	} // switch

	for(int i=0; i<oConstants.poPulsars[nPulsar].nObservations; i++) {
	  switch(oConstants.poSources[nSource].eCorrelation) {
	    case CID_GR:
	      dCor = double(oData.mdGeometric[nPulsar][nPulsar]);
	      break;
	    case CID_All:
	      dCor = 1;
	      break;
	    case CID_Uniform:
	    case CID_SinglePulsar:
	      dCor = 1;
	      break;
	    case CID_Metric_Breathing:
	    case CID_Metric_Shear:
	    case CID_Metric_Longitudinal:
	    default:
	      dCor = 1;
	      break;
	  } // switch eCorrelation

	  for(int j=0; j<oConstants.poPulsars[nPulsar].nObservations; j++) {
	    if(SourceWorksOnResidual(oConstants, nSource, nPulsar, i) &&
		SourceWorksOnResidual(oConstants, nSource, nPulsar, j)) {
	      switch(oConstants.poSources[nSource].oSourceType.eID) {
		case SID_Errorbars:
		  if(i==j)
		    dTemp = exp(oParameters.pdPar[oConstants.poSources[nSource].nFirstParIndex])*
		      exp(oParameters.pdPar[oConstants.poSources[nSource].nFirstParIndex])*
		      oConstants.poPulsars[nPulsar].pdDeltaResiduals[i]*
		      oConstants.poPulsars[nPulsar].pdDeltaResiduals[i];
		  else
		    dTemp = 0;
		  break;
		case SID_White:
		  if(i==j) {
		    dTau = fabs(oConstants.poPulsars[nPulsar].pdTOA[i] -
				oConstants.poPulsars[nPulsar].pdTOA[j]);
		    dTemp = SpectrumIntegral(
			oConstants.poSources[nSource].oSourceType.eID,
			oConstants.poPulsars[nPulsar].pdTOA[i],
			oConstants.poPulsars[nPulsar].pdTOA[j],
			vdParametersNull); 
		  } else {
		    dTemp = 0;
		  } // if a==b && i==j
		  break;
        case SID_DMVar:
		      dTau = fabs(oConstants.poPulsars[nPulsar].pdTOA[i] -
				  oConstants.poPulsars[nPulsar].pdTOA[j]);
		      dTemp = SpectrumIntegral(
                      SID_PowerLaw,
                      oConstants.poPulsars[nPulsar].pdTOA[i],
                      oConstants.poPulsars[nPulsar].pdTOA[j],
                      vdParametersNull, true) / (DM_K * DM_K *
                  oConstants.poPulsars[nPulsar].pdFreq[i]*oConstants.poPulsars[nPulsar].pdFreq[i]*
                  oConstants.poPulsars[nPulsar].pdFreq[j]*oConstants.poPulsars[nPulsar].pdFreq[j]);
          break;
        case SID_DipoleGWB:
          dTau = fabs(oConstants.poPulsars[nPulsar].pdTOA[i] -
              oConstants.poPulsars[nPulsar].pdTOA[j]);
          dTemp = SpectrumIntegral(
                  SID_PowerLaw,
                  oConstants.poPulsars[nPulsar].pdTOA[i],
                  oConstants.poPulsars[nPulsar].pdTOA[j],
                  vdParametersNull); 
          dTemp *= AnisotropicGWBCorrelation(nPulsar, nPulsar,
                  oConstants.poSources[nSource].oSourceType.nTag,
                  oConstants, vdParameters);
          break;
		default:
		  dTau = fabs(oConstants.poPulsars[nPulsar].pdTOA[i] -
			      oConstants.poPulsars[nPulsar].pdTOA[j]);
		  dTemp = SpectrumIntegral(
		      oConstants.poSources[nSource].oSourceType.eID,
		      oConstants.poPulsars[nPulsar].pdTOA[i],
		      oConstants.poPulsars[nPulsar].pdTOA[j],
		      vdParametersNull); 
		  break;
	      } // switch
	      mdC[i][j] += dCor*dTemp;
	    } // if SourceWorksOnResidual
	  } // for j
	} // for i
	break;
      case SID_Deterministic: // Do nothing
      case SID_Nothing:
      default:
	break;
    } // switch eID
  } // if nPulsar
} /* SetCoherenceMatrixOneSource */


/* This matrix sets the EFAC covariance matrix for one pulsar */
void SetCoherenceMatrixPulsarEfac(SConstantsType &oConstants, int nPulsar, CMatrix &mdC, double dEfac) {
  if(! mdC.Defined()) throw ELENotDefined;
  if(! mdC.m_pnDimSize[0] == oConstants.poPulsars[nPulsar].nObservations) throw ELEDimensionMisMatch;
  if(! mdC.m_pnDimSize[1] == oConstants.poPulsars[nPulsar].nObservations) throw ELEDimensionMisMatch;

  for(int i=0; i<oConstants.poPulsars[nPulsar].nObservations; ++i) {
    for(int j=0; j<oConstants.poPulsars[nPulsar].nObservations; ++j) {
      if(i == j && dEfac >= 0) {
	mdC[i][j] = dEfac * dEfac * oConstants.poPulsars[nPulsar].pdDeltaResiduals[i] * oConstants.poPulsars[nPulsar].pdDeltaResiduals[i];
      } else {
	mdC[i][j] = 0;
      } /* if i == j */
    } /* for j */
  } /* for i */
} // SetCoherenceMatrixPulsarEfac

/* This matrix sets the EQUAD covariance matrix for one pulsar */
void SetCoherenceMatrixPulsarEquad(SConstantsType &oConstants, int nPulsar, CMatrix &mdC, double dEquad) {
  if(! mdC.Defined()) throw ELENotDefined;
  if(! mdC.m_pnDimSize[0] == oConstants.poPulsars[nPulsar].nObservations) throw ELEDimensionMisMatch;
  if(! mdC.m_pnDimSize[1] == oConstants.poPulsars[nPulsar].nObservations) throw ELEDimensionMisMatch;

  for(int i=0; i<oConstants.poPulsars[nPulsar].nObservations; ++i) {
    for(int j=0; j<oConstants.poPulsars[nPulsar].nObservations; ++j) {
      if(i == j && dEquad >= 0) {
	mdC[i][j] = dEquad * dEquad;
      } else {
	mdC[i][j] = 0;
      } /* if i == j */
    } /* for j */
  } /* for i */
} // SetCoherenceMatrixPulsarEfac

/* This function returns the parameter number of the efac for a pulsar */
int GetSourceEfac(SConstantsType &oConstants, int nPulsar) {
  int nEfacNumber = -1;

  for(int s=0; s<oConstants.nSources; ++s) {
    if(oConstants.poSources[s].oSourceType.eID == SID_Errorbars &&
      SourceWorksOnPulsar(oConstants, s, nPulsar)) {

      nEfacNumber = oConstants.poSources[s].nFirstParIndex;
      break;
    } /* if eID */
  } /* for s */

  return nEfacNumber;
} // GetSourceEfac

/* This function returns the parameter number of the efac for a pulsar */
int GetSourceEquad(SConstantsType &oConstants, int nPulsar) {
  int nEquadNumber = -1;

  for(int s=0; s<oConstants.nSources; ++s) {
    if(oConstants.poSources[s].oSourceType.eID == SID_White &&
      SourceWorksOnPulsar(oConstants, s, nPulsar)) {

      nEquadNumber = oConstants.poSources[s].nFirstParIndex;
      break;
    } /* if eID */
  } /* for s */

  return nEquadNumber;
} // GetSourceEquad

/* This function returns the parameter number of the efac for a pulsar */
int GetSourcePowerLaw(SConstantsType &oConstants, int nPulsar) {
  int nPlNumber = -1;

  for(int s=0; s<oConstants.nSources; ++s) {
    if(oConstants.poSources[s].oSourceType.eID == SID_PowerLaw &&
      SourceWorksOnPulsar(oConstants, s, nPulsar) &&
      oConstants.poSources[s].eCorrelation == CID_SinglePulsar) {

      nPlNumber = oConstants.poSources[s].nFirstParIndex;
      break;
    } /* if eID */
  } /* for s */

  return nPlNumber;
} // GetSourcePowerLaw

/* This function returns the parameter number of GR correlated signal */
int GetSourceGRCorr(SConstantsType &oConstants, int nPulsar) {
  int nPlNumber = -1;

  if(nPulsar == -1) {
    for(int s=0; s<oConstants.nSources; ++s) {
      if(oConstants.poSources[s].oSourceType.eID == SID_PowerLaw &&
	oConstants.poSources[s].eCorrelation == CID_GR) {

	nPlNumber = oConstants.poSources[s].nFirstParIndex;
	break;
      } /* if eID */
    } /* for s */
  } else {
    for(int s=0; s<oConstants.nSources; ++s) {
      if(oConstants.poSources[s].oSourceType.eID == SID_PowerLaw &&
	SourceWorksOnPulsar(oConstants, s, nPulsar) &&
	oConstants.poSources[s].eCorrelation == CID_GR) {

	nPlNumber = oConstants.poSources[s].nFirstParIndex;
	break;
      } /* if eID */
    } /* for s */
  } /* if nPulsar */

  return nPlNumber;
} // GetSourceGRCorr

/* This function returns the parameter number of uniformly correlated signal,
 * like a clock-signal
 * */
int GetSourceUniCorr(SConstantsType &oConstants, int nPulsar) {
  int nPlNumber = -1;

  if(nPulsar == -1) {
    for(int s=0; s<oConstants.nSources; ++s) {
      if(oConstants.poSources[s].oSourceType.eID == SID_PowerLaw &&
	oConstants.poSources[s].eCorrelation == CID_Uniform) {

	nPlNumber = oConstants.poSources[s].nFirstParIndex;
	break;
      } /* if eID */
    } /* for s */
  } else {
    for(int s=0; s<oConstants.nSources; ++s) {
      if(oConstants.poSources[s].oSourceType.eID == SID_PowerLaw &&
	SourceWorksOnPulsar(oConstants, s, nPulsar) &&
	oConstants.poSources[s].eCorrelation == CID_Uniform) {

	nPlNumber = oConstants.poSources[s].nFirstParIndex;
	break;
      } /* if eID */
    } /* for s */
  } /* if nPulsar */

  return nPlNumber;
} // GetSourceUniCorr


/* This function calculates the full Coherence matrix of the pulsar timing
 * array data, and saves is to disk. Then an external program can calculate the
 * eigenvectors/values
 *
 * Note: Deprecated*/
void SaveCoherenceMatrix(SDataType &oData, SConstantsType &oConstants) {
#if 0
  double dTau;
  bool bCheck;
  CMatrix mdTemp;
  FILE *pFile;
  char strFileName[160], strMsg[160];

  PrintStatus("Allocating memory...");
  oData.mdC.Initialize(oConstants.k*oConstants.l, oConstants.k*oConstants.l);
  oData.mdInvC.Initialize(oConstants.k*oConstants.l, oConstants.k*oConstants.l);
  oData.mdSpectral.Initialize(oConstants.l, oConstants.l);
  mdTemp = oData.mdSpectral;
  PrintSuccess();

  PrintStatus("Setting coherence values...");
  // First set the GWB part of the coherence matrix
  for(int i=0; i<oConstants.l; i++) {
    for(int j=0; j<oConstants.l; j++) {
      dTau = fabs(double(oData.vdMeasureTimes[i])-double(oData.vdMeasureTimes[j]));
      oData.mdSpectral[i][j] = GWBSpectrumIntegral(dTau, oConstants.dAmp, oConstants.dExp, oConstants.dCutoffFrequency);
      for(int a=0; a<oConstants.k; a++) {
        for(int b=0; b<oConstants.k; b++) {
          oData.mdC[a+i*oConstants.k][b+j*oConstants.k] = double(oData.mdGeometric[a][b])*double(oData.mdSpectral[i][j]);
        } // for b
      } // for a
    } // for j
  } // for i

  // Secondly set the PN part of the coherence matrix
  for(int a=0; a<oConstants.k; a++) {
    for(int i=0; i<oConstants.l; i++) {
      for(int j=0; j<oConstants.l; j++) {
        dTau = fabs(double(oData.vdMeasureTimes[i])-double(oData.vdMeasureTimes[j]));
        mdTemp[i][j] = PNSpectrumIntegral(dTau, double(oConstants.vdAmpi[a]), double(oConstants.vdExpi[a]), oConstants.dCutoffFrequency);
//        mdTemp[i][j] = GWBSpectrumIntegral(dTau, double(oConstants.vdAmpi[a]), double(oConstants.vdExpi[a])+2, oConstants.dCutoffFrequency);
        oData.mdC[a+i*oConstants.k][a+j*oConstants.k] += double(mdTemp[i][j]);
      } // for j
    } // for i
  } // for a
  PrintSuccess();

  // Now save the coherence matrix to disk
  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcat(strFileName, "coherencematrix.txt");

  strcpy(strMsg, "Writing coherence matrix to '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  try {
    if(! (pFile = fopen(strFileName, "w+")) ) throw 1;

    for(int i=0; i<oConstants.k*oConstants.l; i++) {
      for(int j=0; j<oConstants.k*oConstants.l; j++) {
	fprintf(pFile, "%.16e   ", double(oData.mdC[i][j]));
      } // for j
      fprintf(pFile, "\n");
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return ;
  } // try
  PrintSuccess();
#else
  PrintStatus("Skipping save-coherence matrix");
  PrintSuccess();
#endif
} // SaveCoherenceMatrix



/* This function reads a file with eigenvalues and a file with eigenvectors of
 * the coherence matrix. Also a file with the coherencematrix itself is read.
 *
 * Note: deprecated*/
void LoadCoherenceMatrix(SDataType &oData, SConstantsType &oConstants) {
#if 0
  double dTau;
  bool bCheck;
  CMatrix mdTemp;
  FILE *pFile;
  char strFileName[160], strMsg[160];
  double dBuf;

  PrintStatus("Allocating memory...");
  oData.mdC.Initialize(oConstants.k*oConstants.l, oConstants.k*oConstants.l);
  oData.mdEVC.Initialize(oConstants.k*oConstants.l, oConstants.k*oConstants.l);
  oData.mdInvC.Initialize(oConstants.k*oConstants.l, oConstants.k*oConstants.l);
  oData.mdSpectral.Initialize(oConstants.l, oConstants.l);
  oData.vdLambdaC.Initialize(oConstants.k*oConstants.l);
  mdTemp = oData.mdSpectral;
  PrintSuccess();

  // Now load the coherence matrix from disk
  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcat(strFileName, "coherencematrix.txt");

  strcpy(strMsg, "Loading coherence matrix from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");
  PrintStatus(strMsg);
  try {
    if(! (pFile = fopen(strFileName, "r+")) ) throw 1;

    for(int i=0; i<oConstants.k*oConstants.l; i++) {
      for(int j=0; j<oConstants.k*oConstants.l; j++) {
	fscanf(pFile, "%e", &dBuf);
        oData.mdC[i][j] = dBuf;
      } // for j
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return ;
  } // try
  PrintSuccess();


  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcat(strFileName, "eigenvectors.txt");

  strcpy(strMsg, "Loading eigenvector-matrix from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");
  PrintStatus(strMsg);
  try {
    if(! (pFile = fopen(strFileName, "r+")) ) throw 1;

    for(int i=0; i<oConstants.k*oConstants.l; i++) {
      for(int j=0; j<oConstants.k*oConstants.l; j++) {
	fscanf(pFile, "%e", &dBuf);
        oData.mdEVC[i][j] = dBuf;
      } // for j
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return ;
  } // try
  PrintSuccess();


  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcat(strFileName, "eigenvalues.txt");

  strcpy(strMsg, "Loading eigenvalue-vector from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");
  PrintStatus(strMsg);
  try {
    if(! (pFile = fopen(strFileName, "r+")) ) throw 1;

    for(int i=0; i<oConstants.k*oConstants.l; i++) {
      fscanf(pFile, "%e", &dBuf);
      oData.vdLambdaC[i] = dBuf;
    } // for i

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return ;
  } // try
  PrintSuccess();
#else
  PrintStatus("Skipping save-coherence matrix");
  PrintSuccess();
#endif

} // LoadCoherenceMatrix




/* This function Calculates a specified amount of datapoints for the pta
 * problem many times, and checks whether the coherence matrix is then
 * reproduced.
 *
 * Last update: 2009-04-13
 *
 * Note: Only for debugging purposes.
*/
void CheckGenerateResiduals(SDataType &oData, SConstantsType &oConstants) {
  CVector vdRand(oConstants.k*oConstants.l), vdTemp, vdTemp2, vdTemp3;
  CVector vdLambdaCPos;
  CMatrix mdXX;
  int nIterations=oConstants.m;
  int nWarning=0;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, (unsigned int)time(NULL));

  PrintStatus("Generating/checking many realisations of datapoints...");

  // Correct for negative eigenvalues
  vdLambdaCPos = oData.vdLambdaC;
  for(int i=0; i<oConstants.n; i++) {
    if(double(vdLambdaCPos[i]) < 0) {
      vdLambdaCPos[i] = 0;
      nWarning++;
    } // if vdTemp
  } // for i

  mdXX = oData.mdC;
  for(int i=0; i<oConstants.n; i++) {
    for(int j=0; j<oConstants.n; j++) {
      mdXX[i][j] = 0.0;
    } // for j
  } // for i

  for(int a=0; a<nIterations; a++) {
    // vdRand.Randomize();
    for(int b=0; b<oConstants.n; b++)
      vdRand[b] = gsl_ran_ugaussian(r);
    vdTemp = oData.mdEVC*(Sqrt(vdLambdaCPos)&&vdRand);

    for(int i=0; i<oConstants.n; i++) {
      for(int j=0; j<oConstants.n; j++) {
	mdXX[i][j] += double(vdTemp[i])*double(vdTemp[j]);
      } // for j
    } // for i
  } // for a


  for(int i=0; i<oConstants.n; i++) {
    for(int j=0; j<oConstants.n; j++) {
      mdXX[i][j] /= nIterations;
    } // for j
  } // for i
  
  vdTemp.Initialize(oConstants.n*oConstants.n);
  vdTemp2 = vdTemp;
  vdTemp3 = vdTemp;
  for(int i=0; i<oConstants.n; i++) {
    for(int j=0; j<oConstants.n; j++) {
      vdTemp[i+oConstants.n*j] = double(oData.mdC[i][j]);
      vdTemp2[i+oConstants.n*j] = double(mdXX[i][j]);
      vdTemp3[i+oConstants.n*j] = double(oData.mdC[i][j]) - double(mdXX[i][j]);
    } // for j
  } // for i

  WritePlot("testdata.txt", vdTemp, vdTemp2, vdTemp3);

  PrintSuccess();

  if(nWarning > 0)
    printf("WARNING: Corrected for %i negative eigenvalues!\n", nWarning);


#if 0
  bool bCheck;
  CMatrix mdD;
  CVector vdLambdaNew, vdEig;
//  CLinearObject *cloTemp;
//  cloTemp = oData.mdC.Eigen("V");
//  oData.vdLambdaC = *((CVector *)&cloTemp[0]);
//  oData.mdEVC = *((CMatrix *)&cloTemp[1]);
//  delete[] cloTemp;

  bCheck = true;
  PrintStatus("Checking cholesky on mdC...");
  try { LogDetChol(oData.mdC); } catch (ELinearError err) { bCheck = false;} 
  if(bCheck) PrintSuccess(); else PrintFailed();

  mdD = oData.mdC;
  vdEig = oData.vdLambdaC;
  vdLambdaNew = oData.vdLambdaC;
  for(int i=0; i<oConstants.n; i++) {
    for(int j=0; j<oConstants.n; j++) {
      mdD[i][j] = 0;
      if(i==j) {
	for(int k=0; k<oConstants.k; k++) {
	  vdEig[k] = double(oData.mdEVC[k][i]);
	} // for k
//	mdD[i][j] = double(oData.vdLambdaC[i]);
        mdD[i][j] = double(vdEig * oData.mdC * vdEig) / double(vdEig * vdEig);
      } // if i
    } // for j
  } // for i
  oData.mdX = oData.mdEVC * mdD * oData.mdEVC[LO_TRANSPOSE];
#endif

  gsl_rng_free(r);
  return;
} // CheckGenerateResiduals



/* This function Calculates a specified number of datapoints for the pta
 * problem.
 *
 * After the generation, it averages the residuals to zero (per pulsar)
 *
 * Last update: 2008-11-20
*/
void GenerateResiduals(SDataType &oData, SConstantsType &oConstants, bool bEraseDeltaResiduals) {
  CVector vdRand(oConstants.n), vdTemp, vdTemp2, vdTemp3, vdQSDFit, vdLambdaPos;
  CMatrix mdTemp;
  SParametersType oTempParameters;
  int nWarning=0;
  int nIndex;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, (unsigned int)time(NULL));
//  gsl_rng_set(r, (unsigned int)5);
  SetBlockMatrix(oData, oConstants);

  vdLambdaPos = oData.vdLambdaC;
#if 0
  { // This is based on the correct derivation
    double dVariance=0;

    CMatrix mdBlock, mdBlockTrans, mdC, mdInvC, mdB, mdBTrans,
	    mdTemp1, mdTemp2, mdTemp3, mdI;
    int nIndex1, nIndex2;

    // Use a different mdC here
    nIndex1 = 0;
    nIndex2 = 0;
    mdC.Initialize(oConstants.n, oConstants.n);
    mdInvC.Initialize(oConstants.n, oConstants.n);
    mdI.Initialize(oConstants.n, oConstants.n);
    for(int a=0; a<oConstants.k; a++ ) {
      for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	for(int b=0; b<oConstants.k; b++ ) {
	  for(int j=0; j<oConstants.poPulsars[a].nObservations; j++) {
	    if(nIndex1 == nIndex2) {
	      mdC[nIndex1][nIndex2] = gsl_pow_2(oConstants.poPulsars[a].pdDeltaResiduals[i]);
	      mdInvC[nIndex1][nIndex2] = 1.0 / gsl_pow_2(oConstants.poPulsars[a].pdDeltaResiduals[i]);
	      mdI[nIndex1][nIndex2] = 1.0;
	    } else {
	      mdC[nIndex1][nIndex2] = 0;
	      mdInvC[nIndex1][nIndex2] = 0;
	      mdI[nIndex1][nIndex2] = 0.0;
	    } // if nIndex1
	    nIndex2++;
	  } // for j
	} // for b
	nIndex1++;
	nIndex2 = 0;
      } // for i
    } // for a

    mdBlock = oData.mdBlock;
    mdBlockTrans = oData.mdBlock[LO_TRANSPOSE];

    mdTemp1 = mdInvC * mdBlock;
    mdTemp2 = mdBlockTrans * mdTemp1;
    mdTemp2.InvertChol();
    mdTemp1 = mdBlockTrans * mdInvC;
    mdTemp3 = mdTemp2 * mdTemp1;
    mdTemp2 = mdBlock * mdTemp3;
    mdB = mdI - mdTemp2;
    mdBTrans = mdB[LO_TRANSPOSE];

    mdTemp1 = oData.mdC * mdBTrans;
    mdTemp2 = mdB * mdTemp1;

    dVariance = Trace(mdTemp2);

    printf("\n");
    printf("Rms-theory = %e\n", sqrt(dVariance / mdC.m_pnDimSize[0]));
  }
#endif

#if 0
  for(int i=0; i<oConstants.n; i++) {
    if(double(vdLambdaPos[i]) < 0) {
      vdLambdaPos[i] = 0;
      nWarning++;
    } // if vdLambdaPos
  } // for i

  if(nWarning > 0)
    printf("WARNING: Correcting for %i negative eigenvalues!\n", nWarning);

  PrintStatus("Generating datapoints...");

  for(int b=0; b<oConstants.n; b++)
    vdRand[b] = gsl_ran_ugaussian(r);
  vdTemp = oData.mdEVC*(Sqrt(vdLambdaPos)&&vdRand);
#endif

  PrintStatus("Generating datapoints using Cholesky...");

  for(int b=0; b<oConstants.n; b++)
    vdRand[b] = gsl_ran_ugaussian(r);
  vdTemp = oData.mdCholeskyC * vdRand;

//  WritePlot("../data/workdata/eigenvalues.txt", oData.vdLambdaC, vdLambdaPos);

  // Adding some quadratic spindown to the data
  vdTemp2 = vdTemp;
  ParametersFromSources(oTempParameters, oConstants.poSources, oConstants.nSources);
  ResidualsFromDetSources(oConstants, oTempParameters, oData, &vdTemp2);

  vdTemp3 = vdTemp + vdTemp2;

  if(0 && oConstants.nMarParameters > 0) {
//  if(oConstants.nMarParameters > 0) {
    // Fitting for quadratici spindown next. As that is what observers do.
    // Here we only remove the constant component (average to zero)
    mdTemp = oData.mdGenBlock[LO_TRANSPOSE] * oData.mdGenBlock;
    mdTemp.InvertChol();
    vdQSDFit = mdTemp * (oData.mdGenBlock[LO_TRANSPOSE] * vdTemp3);
    oData.vdData = vdTemp3 - oData.mdGenBlock * vdQSDFit;
  } else {
    oData.vdData = vdTemp3;
  } // if nMarParameters

  PrintSuccess();

  gsl_rng_free(r);

  // Now convert this set to the residuals parameters of the pulsar structure
  nIndex=0;
  for(int a=0; a<oConstants.k; a++) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
      oConstants.poPulsars[a].pdResiduals[i] = double(oData.vdData[nIndex]);
      if(bEraseDeltaResiduals)
	oConstants.poPulsars[a].pdDeltaResiduals[i] = 0;
      nIndex++;
    } // for i
  } // for a
  return;
} // GenerateResiduals

/* This function generates an ensemble of datasets.
 *
 * Last update: 2011-08-04
*/
void GenerateEnsembleResiduals(SDataType &oData, SConstantsType &oConstants, int nEnsembles) {
  CMatrix mdRand(oConstants.n, nEnsembles), mdTemp;
  CVector vdDetSources;
  SParametersType oTempParameters;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, (unsigned int)time(NULL));

  SetBlockMatrix(oData, oConstants);
  vdDetSources.Initialize(oConstants.n);
  PrintStatus("Generating ensemble of datapoints using Cholesky...");

  // Generate a matrix with random elements
  for(int i=0; i<oConstants.n; i++) {
    for(int j=0; j<nEnsembles; j++) {
      mdRand[i][j] = gsl_ran_ugaussian(r);
    } // for j
  } // for i

  // Construct the stochastic datasets as columns of mdTemp
  mdTemp = oData.mdCholeskyC * mdRand;

  // Generate the deterministic source data
  ParametersFromSources(oTempParameters, oConstants.poSources, oConstants.nSources);
  ResidualsFromDetSources(oConstants, oTempParameters, oData, &vdDetSources);

  // Now fill the mdDataSets matrix
  oData.mdDataSets.Initialize(oConstants.n, nEnsembles);
  for(int i=0; i<oConstants.n; i++) {
    for(int j=0; j<nEnsembles; j++) {
      oData.mdDataSets[i][j] = double(mdTemp[i][j]) + double(vdDetSources[i]);
    } // for j
  } // for i
  oData.nDataSets = nEnsembles;

  PrintSuccess();
  gsl_rng_free(r);
} // GenerateEnsembleResiduals

/* This function generates an ensemble of datasets, all as independent
 * observations of the same realisation. The first (0-th) set are the true
 * arrival times. All the following datasets are the first set, with radiometer
 * noise superimposed on it.
 * The first set should be used as the true realisations. The second set should
 * be used as the importance-sampling kernel.
 */
void GenerateRadiometerEnsembleResiduals(SDataType &oData, SConstantsType &oConstants, int nEnsembles) {
  CMatrix mdRand(oConstants.n, 1), mdTemp;
  CVector vdDetSources;
  SParametersType oTempParameters;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, (unsigned int)time(NULL));
  int nIndex;

  SetBlockMatrix(oData, oConstants);
  vdDetSources.Initialize(oConstants.n);
  PrintStatus("Generating ensemble of datapoints from one set...");

  // Generate a matrix with random elements from a Gaussian distribution
  for(int i=0; i<oConstants.n; ++i) {
    mdRand[i][0] = gsl_ran_ugaussian(r);
  } // for i

  // Construct the stochastic datasets as columns of mdTemp
  mdTemp = oData.mdCholeskyC * mdRand;

  // Generate the deterministic source data
  ParametersFromSources(oTempParameters, oConstants.poSources, oConstants.nSources);
  ResidualsFromDetSources(oConstants, oTempParameters, oData, &vdDetSources);

  // Now fill the mdDataSets matrix
  oData.mdDataSets.Initialize(oConstants.n, nEnsembles);
  for(int i=0; i<oConstants.n; ++i) {
    for(int j=0; j<nEnsembles; ++j) {
      oData.mdDataSets[i][j] = double(mdTemp[i][0]) + double(vdDetSources[i]);
    } // for j
  } // for i

  // Add radiometer noise
  nIndex=0;
  for(int a=0; a<oConstants.k; ++a) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; ++i) {
      // Do not add radiometer noise to the original set
      for(int j=1; j<nEnsembles; j++) {
	oData.mdDataSets[nIndex][j] += oConstants.poPulsars[a].pdDeltaResiduals[i] *
	  gsl_ran_ugaussian(r);
      } // for j
      nIndex++;
    } // for i
  } // for a

  oData.nDataSets = nEnsembles;

  PrintSuccess();
  gsl_rng_free(r);
} // GenerateRadiometerEnsembleResiduals

/* This function generates many datasets, and computes the corresponding data
 * matrix for use in the loglikelihood function 
 * Last update: 2007-10-21
 *
 * Note: deprecated
*/
void GenerateDataSets(SDataType &oData, SConstantsType &oConstants) {
  CVector vdRand(oConstants.k*oConstants.l), vdTemp, vdTemp2, vdTemp3;
  CMatrix mdXX;
  int nIterations=oConstants.m;
  gsl_rng * r = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(r, (unsigned int)time(NULL));

  PrintStatus("Generating many datasets...");

  mdXX = oData.mdC;
  for(int i=0; i<oConstants.n; i++) {
    for(int j=0; j<oConstants.n; j++) {
      mdXX[i][j] = 0.0;
    } // for j
  } // for i

  for(int a=0; a<nIterations; a++) {
    // vdRand.Randomize();
    for(int b=0; b<oConstants.n; b++)
      vdRand[b] = gsl_ran_ugaussian(r);
    vdTemp = oData.mdEVC*(Sqrt(oData.vdLambdaC)&&vdRand);

    for(int i=0; i<oConstants.n; i++) {
      for(int j=0; j<oConstants.n; j++) {
	mdXX[i][j] += double(vdTemp[i])*double(vdTemp[j]);
      } // for j
    } // for i
  } // for a

  oData.mdX = mdXX;

  PrintSuccess();
  gsl_rng_free(r);
  return;
} // GenerateDataSets




/* Writes the timing residual data from a datafile: oConstants.strDataDir /
 * residuals.dat
 *
 * TODO: Other parameters, like the cutoff frequency should also be stored in
 * this file. One cannot rely on the file 'parameters.conf' for this, it
 * should be self-contained.
 * THIS STATEMENT IS NOT TRUE ANYMORE. RESIDUALS ARE RESIDUALS. NOTHING MORE.
 * HANDLE PARAMETERS SEPARATE FROM THE DATA!
 *
 * Note: MCMC parameters and plotting parameters that do not interfere with
 * the validity of the residuals can be stored in 'parameters.conf' and
 * should be omitted here.
 * */
bool WriteResiduals(SDataType &oData, SConstantsType &oConstants) {
  FILE *pFile;
  char strFileName[160], strMsg[160];

  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcpy(strFileName, oConstants.strResidualsFile);

  strcpy(strMsg, "Writing residuals to '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);
  if(! oData.mdX.Defined() ) oData.mdX.Initialize(oConstants.k*oConstants.l, oConstants.k*oConstants.l);
  if(! oData.vdData.Defined() ) oData.vdData.Initialize(oConstants.k*oConstants.l);

  try {
    if(! (pFile = fopen(strFileName, "wb+")) ) throw 1;

    if(! fwrite(oConstants.strVersion , sizeof(char) , 16 , pFile) ) throw 2;
    if(! fwrite(&oConstants.k , sizeof(int) , 1 , pFile) ) throw 3;

    // For each pulsar, write the stuff
    for(int a=0; a<oConstants.k; a++) {
      if(! fwrite(&oConstants.poPulsars[a].nObservations , sizeof(int) , 1 , pFile) ) throw 3;
      if(! fwrite(&oConstants.poPulsars[a].dPhi , sizeof(double) , 1 , pFile) ) throw 4;
      if(! fwrite(&oConstants.poPulsars[a].dTheta , sizeof(double) , 1 , pFile) ) throw 5;
      if(! fwrite(oConstants.poPulsars[a].pdTOA , sizeof(double) , oConstants.poPulsars[a].nObservations , pFile) ) throw 6;
    } // for k

    if(! fwrite(oData.vdData.m_pdData , sizeof(double) , oConstants.n , pFile) ) throw 15;

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    printf("Error number: %i\n", nError);
    return false;
  } // try

  PrintSuccess();
  return true;
} // WriteResiduals

/* Reads the timing residual data from a datafile: oConstants.strDataDir /
 * residuals.dat
 *
 * Todo: Other parameters, like the cutoff frequency should also be stored in
 * this file. On cannot rely on the file 'parameters.conf' for this, it
 * should be self-contained.
 * THIS STATEMENT IS NOT TRUE ANYMORE. RESIDUALS ARE RESIDUALS. NOTHING MORE.
 * HANDLE PARAMETERS SEPARATE FROM THE DATA!
 *
 * Note: MCMC parameters and plotting parameters that do not interfere with
 * the validity of the residuals can be stored in 'parameters.conf' and
 * should be omitted here.
 *
 * TODO: The vectors of SConstantsType are not initialized here. So what if we
 * have a larger parameters space for instance?
 * SOLUTION: There is redundant data. This file should be ONLY for residuals. If
 * there is extra data here (larger number of pulsars), read ONLY the relevant
 * residuals.
 * EDIT: This comment is very very outdated. TODO: Switch completely to text
 *       files like tempo2 does
 * */
bool ReadResiduals(SDataType &oData, SConstantsType &oConstants, SParametersType &oParameters) {
  FILE *pFile;
  char strFileName[160], strMsg[160], strFileVersion[16];
  int nIndex;

  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcpy(strFileName, oConstants.strResidualsFile);

  strcpy(strMsg, "Reading residuals from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  try {
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    if(! fread(strFileVersion , sizeof(char) , 16 , pFile) ) throw 2;

    if(strcmp(strFileVersion,oConstants.strVersion ) != 0)
      ; // wrong version

    if(! fread(&oConstants.k , sizeof(int) , 1 , pFile) ) throw 2;
    oConstants.l = 0;
    oConstants.n = 0;

    for(int a=0; a<oConstants.k; a++) {
      if(! fread(&oConstants.poPulsars[a].nObservations , sizeof(int) , 1 , pFile) ) throw 3;
      if(! fread(&oConstants.poPulsars[a].dPhi , sizeof(double) , 1 , pFile) ) throw 4;
      if(! fread(&oConstants.poPulsars[a].dTheta , sizeof(double) , 1 , pFile) ) throw 5;
      if(! fread(oConstants.poPulsars[a].pdTOA , sizeof(double) , oConstants.poPulsars[a].nObservations , pFile) ) throw 6;

      if(oConstants.poPulsars[a].nObservations > oConstants.l)
	oConstants.l = oConstants.poPulsars[a].nObservations;
      oConstants.n += oConstants.poPulsars[a].nObservations;
    } // for a


    oData.vdData.Initialize(oConstants.n);
    oData.mdC.Initialize(oConstants.n,oConstants.n);
    oData.mdX.Initialize(oConstants.n,oConstants.n);

    if(! fread(oData.vdData.m_pdData , sizeof(double) , oConstants.n , pFile) ) throw 13;

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return false;
  } // try

  PrintSuccess();

  // Now convert this set to the residuals parameters of the pulsar structure
  nIndex=0;
  for(int a=0; a<oConstants.k; a++) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
      oConstants.poPulsars[a].pdResiduals[i] = double(oData.vdData[nIndex]);
      oConstants.poPulsars[a].pdDeltaResiduals[i] = 0;
      nIndex++;
    } // for i
  } // for a

  // TODO: Why was SetBlockMatrix here again? Remove it! Ugly programming!
  SetBlockMatrix(oData, oConstants);

  // TODO: Why are we also initializing oParameters? See above!
  for(int s=0; s<oConstants.nSources; s++) {
    for(int i=0; i<oConstants.poSources[s].oSourceType.nParameters; i++) {
      oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+i] =
	oConstants.poSources[s].pdPar[i];
    } // for i
  } // for s
  return true;
} // ReadResiduals




/* Writes the pulsar angles (positions) from a datafile: oConstants.strDataDir /
 * oConstants.strAnglesFile
 *
 * TODO: This shouldn't be handled this way. Specify all pulsar details at one
 * place somewhere else (like in a .par file??)
 * */
bool WriteAngles(SDataType &oData, SConstantsType &oConstants) {
  FILE *pFile;
  char strFileName[160], strMsg[160];

  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcpy(strFileName, oConstants.strAnglesFile);

  strcpy(strMsg, "Writing pulsar angles to '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  try {
    if(! (pFile = fopen(strFileName, "wb+")) ) throw 1;

    if(! fwrite(oConstants.strVersion , sizeof(char) , 16 , pFile) ) throw 2;
    if(! fwrite(&oConstants.k , sizeof(int) , 1 , pFile) ) throw 3;

    // Now write the pulsar angles to this file
    for(int a=0; a<oConstants.k; a++) {
      if(! fwrite(&oConstants.poPulsars[a].dPhi , sizeof(double) , 1 , pFile) ) throw 4;
      if(! fwrite(&oConstants.poPulsars[a].dTheta , sizeof(double) , 1 , pFile) ) throw 5;
    } // for a

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    printf("Error number: %i\n", nError);
    return false;
  } // try

  PrintSuccess();
  return true;
} // WriteAngles


/* Read the pulsar angles (positions) from a datafile: oConstants.strDataDir /
 * oConstants.strAnglesFile
 *
 * */
bool ReadAngles(SDataType &oData, SConstantsType &oConstants) {
  FILE *pFile;
  char strFileName[160], strMsg[160], strBuf[160];
  int k;

  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcpy(strFileName, oConstants.strAnglesFile);

  strcpy(strMsg, "Reading pulsar angles from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  try {
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    if(! fread(strBuf, sizeof(char) , 16 , pFile) ) throw 2;
    if(strcmp(strBuf, oConstants.strVersion) != 0) throw 3;
    if(! fread(&k , sizeof(int) , 1 , pFile) ) throw 4;

    if(k < oConstants.k) throw 6;;

    for(int a=0; a<k; a++) {
      if(! fread(&oConstants.poPulsars[a].dPhi , sizeof(double) , 1 , pFile) ) throw 4;
      if(! fread(&oConstants.poPulsars[a].dTheta , sizeof(double) , 1 , pFile) ) throw 5;
    } // for a

    if(fclose(pFile) ) throw 0;
  } catch(int nError) {
    PrintFailed();
    return false;
  } // try

  PrintSuccess();
  return true;
} // ReadAngles


/* Read the pulsar angles (positions) from a datafile: oConstants.strDataDir /
 * oConstants.strAnglesFile
 *
 * */
bool ReadAppendAngles(SDataType &oData, SConstantsType &oConstants) {
  FILE *pFile;
  char strFileName[160], strMsg[160], strBuf[160];
  int k, nMin=0;
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, (unsigned int)time(NULL));

  strcpy(strFileName, oConstants.strDataDir);
  if(strFileName[strlen(strFileName)-1] != '/') strcat(strFileName, "/");
  strcpy(strFileName, oConstants.strAnglesFile);

  strcpy(strMsg, "Reading pulsar angles from '");
  strcat(strMsg, strFileName);
  strcat(strMsg, "'...");

  PrintStatus(strMsg);

  try {
    if(! (pFile = fopen(strFileName, "rb+")) ) throw 1;

    if(! fread(strBuf, sizeof(char) , 16 , pFile) ) throw 2;
    if(strcmp(strBuf, oConstants.strVersion) != 0) throw 3;
    if(! fread(&k , sizeof(int) , 1 , pFile) ) throw 4;

    nMin = (k < oConstants.k) ? k : oConstants.k;

    for(int a=0; a<oConstants.k; a++) {
      if(a < k) {
	if(! fread(&oConstants.poPulsars[a].dPhi , sizeof(double) , 1 , pFile) ) throw 4;
	if(! fread(&oConstants.poPulsars[a].dTheta , sizeof(double) , 1 , pFile) ) throw 5;
      } else {
	oConstants.poPulsars[a].dPhi = gsl_rng_uniform(rng)*2*M_PI;
	oConstants.poPulsars[a].dTheta = asin(gsl_rng_uniform(rng)*2-1);
      } // if a
    } // for a

    if(fclose(pFile) ) throw 0;

  } catch(int nError) {
    PrintFailed();

    gsl_rng_free(rng);
    return false;
  } // try

  PrintSuccess();
  return true;
} // ReadAppendAngles


/* This function shows the true parameters of the GWB and all pulsar noises
 * on-screen. First the file 'residuals.dat' is read. */
void PrintTrueParameters(SDataType &oData, SConstantsType &oConstants) {
#if 0
  SParametersType oParameters;
  if(ReadResiduals(oData, oConstants, oParameters)) {
    printf("Values of the parameters:\n");
       printf("GW:      Amp: %7.3f E-15 yr^(1/2)    Exp: %5.3f\n", oConstants.dAmp, oConstants.dExp);
    for(int a=0; a<oConstants.k; a++) {
      if(double(oConstants.vdParametersPerPulsar[a]) == 1.0)
	printf("PN[%2i]:  Amp: %7.3f ns\n", a, double(oConstants.mdPar[a][0]));
      else
  	printf("PN[%2i]:  Amp: %7.3f ns               Exp: %5.3f\n", a, double(oConstants.mdPar[a][0]), double(oConstants.mdPar[a][1]));
    } // for a
  } // if ReadResiduals
#endif
} // PrintTrueParameters

void PrintParameters(SConstantsType &oConstants) {
  char *strBuf;
  int nDetSourceType;
  strBuf = new char[oConstants.k<150?150:oConstants.k];
  printf("Sources for configuration: k=%i, l=%i, n=%i\n",
      oConstants.k, oConstants.l, oConstants.n);
  printf("strPulsarName[0]: %s\n", oConstants.poPulsars[0].strPulsarName);
  printf("Default pulsar source: \"%s\", correlation \"%s\", %i parameter(s)\n",
	oConstants.oDefaultPulsarSource.oSourceType.strID,
	strCorrelationID[(int) oConstants.oDefaultPulsarSource.eCorrelation],
	oConstants.oDefaultPulsarSource.oSourceType.nParameters);

  for(int p=0; p<oConstants.oDefaultPulsarSource.oSourceType.nParameters; p++)
    strBuf[p] = (oConstants.oDefaultPulsarSource.oSourceType.nFitTag & FITPARAMETER(p)) ? '1' : '0';
  strBuf[oConstants.oDefaultPulsarSource.oSourceType.nParameters] = '\0';
  printf("   nFitTag:    %s\n", strBuf);
  for(int p=0; p<oConstants.oDefaultPulsarSource.oSourceType.nParameters; p++)
    strBuf[p] = (oConstants.oDefaultPulsarSource.oSourceType.nMarTag & FITPARAMETER(p)) ? '1' : '0';
  strBuf[oConstants.oDefaultPulsarSource.oSourceType.nParameters] = '\0';
  printf("   nMarTag:    %s\n", strBuf);
  for(int p=0; p<oConstants.oDefaultPulsarSource.oSourceType.nParameters; p++)
    strBuf[p] = (oConstants.oDefaultPulsarSource.oSourceType.nWrapTag & FITPARAMETER(p)) ? '1' : '0';
  strBuf[oConstants.oDefaultPulsarSource.oSourceType.nParameters] = '\0';
  printf("   nWrapTag:   %s\n", strBuf);
  for(int p=0; p<oConstants.oDefaultPulsarSource.oSourceType.nParameters; p++)
    strBuf[p] = (oConstants.oDefaultPulsarSource.oSourceType.nReduceTag & FITPARAMETER(p)) ? '1' : '0';
  strBuf[oConstants.oDefaultPulsarSource.oSourceType.nParameters] = '\0';
  printf("   nReduceTag: %s\n", strBuf);
  for(int p=0; p<oConstants.oDefaultPulsarSource.oSourceType.nParameters; p++) {
    printf("   Param[%i]:  %e\n", p, oConstants.oDefaultPulsarSource.pdPar[p]);
    printf("      MinBound:   %e\n", oConstants.oDefaultPulsarSource.pdParMinBound[p]);
    printf("      MaxBound:   %e\n", oConstants.oDefaultPulsarSource.pdParMaxBound[p]);
    printf("      Start:      %e\n", oConstants.oDefaultPulsarSource.pdParStart[p]);
    printf("      WidthMCMC:  %e\n", oConstants.oDefaultPulsarSource.pdParWidthMCMC[p]);
    printf("      WidthFit:   %e\n", oConstants.oDefaultPulsarSource.pdParWidthFit[p]);
    printf("      WidthPrior: %e\n", oConstants.oDefaultPulsarSource.pdParWidthPrior[p]);
    printf("      MeanPrior:  %e\n", oConstants.oDefaultPulsarSource.pdParMeanPrior[p]);
  } // for p
  printf("\n");
  for(int s=0; s<oConstants.nSources; s++) {
    printf("Source[%i]: \"%s\", correlation \"%s\", %i parameter(s) (first: %3i)\n",
	s, oConstants.poSources[s].oSourceType.strID,
	strCorrelationID[(int) oConstants.poSources[s].eCorrelation],
	oConstants.poSources[s].oSourceType.nParameters,
	oConstants.poSources[s].nFirstParIndex);
    if(oConstants.poSources[s].oSourceType.eID == SID_Deterministic) {
      nDetSourceType = 0;
      while(oConstants.poSources[s].oSourceType.nTag != oDetSourceTypes[nDetSourceType].nID)
	nDetSourceType++;
      printf("   strID:   \"%s\"\n", oDetSourceTypes[nDetSourceType].strID);
    } // if eID
    for(int a=0; a<oConstants.k; a++)
      strBuf[a] = oConstants.poSources[s].pbScope[a] ? '1' : '0';
    strBuf[oConstants.k] = '\0';
    printf("   pbScope: %s\n", strBuf);
    printf("   strScopeTag: %s\n", oConstants.poSources[s].strScopeTag);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++)
      strBuf[p] = (oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) ? '1' : '0';
    strBuf[oConstants.poSources[s].oSourceType.nParameters] = '\0';
    printf("   nFitTag:    %s\n", strBuf);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++)
      strBuf[p] = (oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p)) ? '1' : '0';
    strBuf[oConstants.poSources[s].oSourceType.nParameters] = '\0';
    printf("   nMarTag:    %s\n", strBuf);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++)
      strBuf[p] = (oConstants.poSources[s].oSourceType.nWrapTag & FITPARAMETER(p)) ? '1' : '0';
    strBuf[oConstants.poSources[s].oSourceType.nParameters] = '\0';
    printf("   nWrapTag:   %s\n", strBuf);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++)
      strBuf[p] = (oConstants.poSources[s].oSourceType.nReduceTag & FITPARAMETER(p)) ? '1' : '0';
    strBuf[oConstants.poSources[s].oSourceType.nParameters] = '\0';
    printf("   nReduceTag: %s\n", strBuf);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      printf("   Param[%i]:   %e\n", p, oConstants.poSources[s].pdPar[p]);
      printf("      MinBound:   %e\n", oConstants.poSources[s].pdParMinBound[p]);
      printf("      MaxBound:   %e\n", oConstants.poSources[s].pdParMaxBound[p]);
      printf("      Start:      %e\n", oConstants.poSources[s].pdParStart[p]);
      printf("      WidthMCMC:  %e\n", oConstants.poSources[s].pdParWidthMCMC[p]);
      printf("      WidthFit:   %e\n", oConstants.poSources[s].pdParWidthFit[p]);
      printf("      WidthPrior: %e\n", oConstants.poSources[s].pdParWidthPrior[p]);
      printf("      MeanPrior:  %e\n", oConstants.poSources[s].pdParMeanPrior[p]);
    } // for p
    printf("\n");
  } // for s
  printf("\n");
  delete[] strBuf;
  return ;
} // PrintParameters

void PrintParameters(SConstantsType &oConstants, SParametersType &oParameters) {
  char *strBuf;
  strBuf = new char[oConstants.k+1];
  printf("Sources for configuration: k=%i, l=%i, n=%i\n",
      oConstants.k, oConstants.l, oConstants.n);
  printf("Default pulsar source: \"%s\", correlation \"%s\", %i parameter(s)\n",
	oConstants.oDefaultPulsarSource.oSourceType.strID,
	strCorrelationID[(int) oConstants.oDefaultPulsarSource.eCorrelation],
	oConstants.oDefaultPulsarSource.oSourceType.nParameters);

  for(int s=0; s<oConstants.nSources; s++) {
    printf("Source[%i]: \"%s\", correlation \"%s\", %i parameter(s)\n",
	s, oConstants.poSources[s].oSourceType.strID,
	strCorrelationID[(int) oConstants.poSources[s].eCorrelation],
	oConstants.poSources[s].oSourceType.nParameters);
    for(int a=0; a<oConstants.k; a++)
      strBuf[a] = oConstants.poSources[s].pbScope[a] ? '1' : '0';
    strBuf[oConstants.k] = '\0';
    printf("   pbScope: %s\n", strBuf);
    printf("   strScopeTag: %s\n", oConstants.poSources[s].strScopeTag);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++)
      strBuf[p] = (oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) ? '1' : '0';
    strBuf[oConstants.poSources[s].oSourceType.nParameters] = '\0';
    printf("   nFitTag:    %s\n", strBuf);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++)
      strBuf[p] = (oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p)) ? '1' : '0';
    strBuf[oConstants.poSources[s].oSourceType.nParameters] = '\0';
    printf("   nMarTag:    %s\n", strBuf);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++)
      strBuf[p] = (oConstants.poSources[s].oSourceType.nWrapTag & FITPARAMETER(p)) ? '1' : '0';
    strBuf[oConstants.poSources[s].oSourceType.nParameters] = '\0';
    printf("   nWrapTag:   %s\n", strBuf);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++)
      strBuf[p] = (oConstants.poSources[s].oSourceType.nReduceTag & FITPARAMETER(p)) ? '1' : '0';
    strBuf[oConstants.poSources[s].oSourceType.nParameters] = '\0';
    printf("   nReduceTag: %s\n", strBuf);
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      printf("   Param[%i]:   %e\n", p, oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+p]);
      printf("      MinBound:   %e\n", oConstants.poSources[s].pdParMinBound[p]);
      printf("      MaxBound:   %e\n", oConstants.poSources[s].pdParMaxBound[p]);
      printf("      Start:      %e\n", oConstants.poSources[s].pdParStart[p]);
      printf("      WidthMCMC:  %e\n", oConstants.poSources[s].pdParWidthMCMC[p]);
      printf("      WidthFit:   %e\n", oConstants.poSources[s].pdParWidthFit[p]);
      printf("      WidthPrior: %e\n", oConstants.poSources[s].pdParWidthPrior[p]);
      printf("      MeanPrior: %e\n", oConstants.poSources[s].pdParMeanPrior[p]);
    } // for p
    printf("\n");
  } // for s
  printf("\n");
  delete[] strBuf;
  return ;
} // PrintParameters

void PrintResiduals(SConstantsType & oConstants) {
  for(int a=0; a<oConstants.k; ++a) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; ++i) {
      printf("%i: %e   %e   %e\n", i,
	  oConstants.poPulsars[a].pdTOA[i],
	  oConstants.poPulsars[a].pdResiduals[i],
	  oConstants.poPulsars[a].pdDeltaResiduals[i]);
    } // for i
  } // for a
  return ;
} // PrintResiduals


/* This function calculates the value of an integral of the form
 * \int_{0}^{\infty}df  cos(2*pi*f*t) * S(f)
 *
 * Parameters: dTau = t [s]
 *             nSpectrumType = [0,1,2,3,4,...]
 *             vdParameters
 *
 * Assumptions: dFc * t < 1, dExp < -1
 *
 * The form of S(f) is given by the eID parameter. Examples:
 * White noise - vdParameters[0] = log(Amplitude), log(A) [log(s)]
 *
 * Power Law   - vdParameters[0] = log(Amplitude), log(A) [log(s)]
 *               vdParameters[1] = Exponent, Gamma [] (units of vHL2012)
 *               vdParameters[2] = Cut-off frequency, fL [0.5 yr^-1]
 *
 * Exponential - vdParameters[0] = log(Amplitude), log(A) [log(s)]
 *               vdParameters[1] = Exponent, Gamma [0.5 yr^-1]
 *
 * Lorentzian  - vdParameters[0] = Amplitude, log(A) [s]
 *               vdParameters[1] = Typical frequency, f0 [0.5 yr^-1]
 *
 *
 * To allow for fast computation of the GWB covariance matrix it is also possible
 * to allow for an empty (undefinded) vdParameters vector. Then, static
 * variables are used to evaluate the integral as fast as possible (then only
 * dTau is treated as a variable, the rest is assumed the same as before).
 * Usage:
 * Initialize with normal function call (and defined vdParameters)
 * Fast computation is used when vdParameters is undefined
 * */
double SpectrumIntegral(ESourceID eID, double dTi, double dTj, CVector &vdParameters, bool bDMUnits) {
  double dTau = fabs(dTj - dTi), dTij = 2*M_PI*dTau, dReturnValue=0;
  double dSum, dPrevSum;
  const int nMaxGeometricOrder=80;
  static double dGeometricCoefficients[nMaxGeometricOrder];
  static double dN, dAp, dAwn, dAe, dAl, dG, dF0, dFl, dExp;
  static double dAns, dGns, dTns;
  static double dN5, dN6, dGammaSine5, dGammaSine6;
  static double dGammaSine=0;
  static bool bInitialised=false, bInitialisedN=false;
  static bool bInitialised5=false, bInitialised6=false;
  static bool bInitialised0=false, bInitialised1=false;
  static bool bFirstRun=true;
  unsigned int n=0;

  if(bFirstRun) {
    // This is the first run of the function; initialise the vector
    dGeometricCoefficients[0] = -1.0/(1+dG);
    for(int i=1; i<nMaxGeometricOrder; i++) {
      dGeometricCoefficients[i] = double(gsl_pow_int(-1, i))/gsl_sf_fact(2*i);
    } // for i
    bFirstRun = false;
  } // bFirstRun

  switch(eID) {
    case SID_Nothing: // Nothing
      dReturnValue = 0;
      break;
    case SID_White: // White noise
      if(vdParameters.Defined()) {
  	dAwn = exp(double(vdParameters[0]));
      } // if vdParameters
//      if(fabs(dTij) <= 0.01 * 3600 * 24)
      if(fabs(dTij) == 0.0)
	dReturnValue = dAwn * dAwn;
      else
	dReturnValue = 0.0;
      break;
    case SID_PowerLaw: // Power-law
      if(vdParameters.Defined()) {
        if(bDMUnits) {
          dAp = exp(double(vdParameters[0]));
          dG = double(vdParameters[1]) - 2.0; // The -2.0 is to change definition of units
          dFl = double(vdParameters[2]) / SPERYEAR;
          dN = gsl_pow_int(dAp,2)/pow(dFl*SPERYEAR, 1+dG);
        } else {
          dAp = exp(double(vdParameters[0])) * 1E-15 * sqrt(SPERYEAR / 3.0);
          dG = double(vdParameters[1]) - 2.0; // The -2.0 is to change definition of units
          dFl = double(vdParameters[2]) / SPERYEAR;

          dN = gsl_pow_int(dAp,2)*SPERYEAR/
              (pow(dFl*SPERYEAR, 1+dG)*gsl_pow_int(2*M_PI,2));
        }

	// First the Fl dependent part: the Generalized Hypergeometrical Function
	dSum = -1.0/(1+dG);
	dPrevSum = 0;
	while(dPrevSum != dSum) {
	  n++;
	  dPrevSum = dSum;
	  dSum += gsl_pow_int(-1, n)*gsl_pow_int(dFl*dTij, 2*n)/(gsl_sf_fact(2*n)*(2*n-1-dG));

	  if(n > 100) {
	    printf("Too many iterations\n");
	    break;
	  }
	} // while

	// Make sure the gamma function is possible
	if(dG ==  3.0) dG -= 0.0001;
	if(dG ==  1.0) dG += 0.0001;
	if(dG ==  0.0) dG += 0.0001;
	if(dG == -1.0) dG += 0.0001;

	dGammaSine = gsl_sf_gamma(-1-dG)*sin(-M_PI_2*dG);
	dReturnValue = dN*(dGammaSine*pow(dFl*dTij, dG+1) - dSum);
	bInitialised = true;
      } else {
	// Assume that bInitialised == true
	if(! bInitialised) {
	  printf("SpectrumIntegral not initialised!\n");
	  return 0;
	} // if bInitialised

	// First the Fl dependent part: the Generalized Hypergeometrical Function
	dSum = -1.0/(1+dG);
	for(int i=1; i<nMaxGeometricOrder; i++) {
	  dPrevSum = dSum;
	  dSum += gsl_pow_int(dFl*dTij, 2*i)*dGeometricCoefficients[i]/(2*i-1-dG);
	  if(dPrevSum == dSum) break;
	} // for i

	dReturnValue = dN*(dGammaSine*pow(dFl*dTij, dG+1) - dSum);
      } // if vdParameters
      break;
    case SID_Exponential: // Exponential
      if(vdParameters.Defined()) {
	dAe = exp(double(vdParameters[0]));
	dG = double(vdParameters[1]) * SPERYEAR;
      } // if vdParameters

      // Check whether we have white noise
      if(dG == 0.0) {
	if(dTau == 0.0) dReturnValue = gsl_pow_int(dAe,2);
	else {
	  printf("dG=0.0");
	  dReturnValue = 0.0;
	} // if dTau
      } else {
	if(dExp < 0) {
	  printf("N.E.");
	} // if dExp
	dReturnValue = gsl_pow_int(dAe*dG,2)/(gsl_pow_int(dG,2) + gsl_pow_int(dTij,2));
      } // if dExp
      break;
    case SID_Lorentzian: // Lorentzian
      if(vdParameters.Defined()) {
	dAl = exp(double(vdParameters[0]));
	dF0 = double(vdParameters[1]);
      } // if vdParameters

      if(dF0 <= 0.0) {
	printf("ZN.E.");
      } // if dExp
      // Value of the returnvalue
      dReturnValue = dAl*dAl*exp(-fabs(dF0*dTij));
      break;
    case SID_NonStationary: // Non-Stationary red-noise (testing purposes)
      if(vdParameters.Defined()) {
	dAns = exp(double(vdParameters[0]));
	dGns = double(vdParameters[1]); // * SPERYEAR;
	dTns = double(vdParameters[2]);
      } // if vdParameters

      dTi -= dTns;
      dTj -= dTns;

      if(dGns <= 0.0) {
	if(dTau == 0.0) {
	  dReturnValue = gsl_pow_int(dAns, 2);
	} else {
	  dReturnValue = 0.0;
	} // if dTau
      } else {
	if(dTij != 0) {
	  dReturnValue = dAns*dAns * (1 + (2*M_PI / dGns) * (
		atan(dTj * dGns / fabs(dTij)) - atan(dTi * dGns / fabs(dTij))
	      ));
	} else {
	  dReturnValue = gsl_pow_int(dAns, 2);
	} // if dTij
//	printf("%e  ", dReturnValue);
      } // if dG
      break;
    default:
      printf("SpectrumIntegral(): default value called\n");
      break;
  } // switch nSpectrumType
  return dReturnValue;
} // SpectrumIntegral()


void FitParametersFromVector(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec) {
  int nIndex=0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+p] = double(vdVec[nIndex]);
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
  return ;
} // FitParametersFromVector

void VectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec) {
  int nIndex=0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	vdVec[nIndex] = oParameters.pdPar[oConstants.poSources[s].nFirstParIndex+p];
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
  return ;
} // VectorFromFitParameters

void MinVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec) {
  int nIndex=0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	vdVec[nIndex] = oParameters.pdParMinBound[oConstants.poSources[s].nFirstParIndex+p];
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
  return ;
} // MinVectorFromFitParameters

void MaxVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec) {
  int nIndex=0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	vdVec[nIndex] = oParameters.pdParMaxBound[oConstants.poSources[s].nFirstParIndex+p];
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
  return ;
} // MaxVectorFromFitParameters

void WidthMCMCVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec) {
  int nIndex=0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	vdVec[nIndex] = oParameters.pdParWidthMCMC[oConstants.poSources[s].nFirstParIndex+p];
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
  return ;
} // WidthMCMCVectorFromFitParameters

void WidthFitVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec) {
  int nIndex=0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	vdVec[nIndex] = oParameters.pdParWidthFit[oConstants.poSources[s].nFirstParIndex+p];
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
  return ;
} // WidthFitVectorFromFitParameters

void WidthPriorVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec) {
  int nIndex=0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	vdVec[nIndex] = oParameters.pdParWidthPrior[oConstants.poSources[s].nFirstParIndex+p];
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
  return ;
} // WidthPriorVectorFromFitParameters

void MeanPriorVectorFromFitParameters(SParametersType &oParameters, SConstantsType &oConstants, CVector &vdVec) {
  int nIndex=0;
  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	vdVec[nIndex] = oParameters.pdParMeanPrior[oConstants.poSources[s].nFirstParIndex+p];
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
  return ;
} // MeanPriorVectorFromFitParameters

void ParametersFromSources(SParametersType &oParameters, SSource *poSources, int nSources) {
  for(int s=0; s<nSources; s++)
    for(int p=0; p<poSources[s].oSourceType.nParameters; p++)
      oParameters.pdPar[poSources[s].nFirstParIndex+p] = poSources[s].pdPar[p];
  return ;
} // ParametersFromSources

void StartParametersFromSources(SParametersType &oParameters, SSource *poSources, int nSources) {
  for(int s=0; s<nSources; s++) {
    for(int i=0; i<poSources[s].oSourceType.nParameters; i++) {
      oParameters.pdParStart[poSources[s].nFirstParIndex+i] =
	poSources[s].pdParStart[i];
      oParameters.pdParMinBound[poSources[s].nFirstParIndex+i] =
	poSources[s].pdParMinBound[i];
      oParameters.pdParMaxBound[poSources[s].nFirstParIndex+i] =
	poSources[s].pdParMaxBound[i];
      oParameters.pdParWidthMCMC[poSources[s].nFirstParIndex+i] =
	poSources[s].pdParWidthMCMC[i];
      oParameters.pdParWidthFit[poSources[s].nFirstParIndex+i] =
	poSources[s].pdParWidthFit[i];
      oParameters.pdParWidthPrior[poSources[s].nFirstParIndex+i] =
	poSources[s].pdParWidthPrior[i];
      oParameters.pdParMeanPrior[poSources[s].nFirstParIndex+i] =
	poSources[s].pdParMeanPrior[i];
    } // for i
  } // for s
} // StartParametersFromSources

void ParametersFromStartParameters(SParametersType &oParameters, SSource *poSources, int nSources) {
  for(int s=0; s<nSources; s++) {
    for(int i=0; i<poSources[s].oSourceType.nParameters; i++) {
      oParameters.pdPar[poSources[s].nFirstParIndex+i] =
	oParameters.pdParStart[poSources[s].nFirstParIndex+i];
    } // for i
  } // for s
} // AllParametersFromSources

/* This function converts the oParameters struct to a vector */
void ParametersToVector(SParametersType &oParameters, CVector &vdTheta) {
  for(int i=0; i<vdTheta.m_pnDimSize[0]; i++)
    vdTheta[i] = oParameters.pdPar[i];
} // ParametersToVector

/* This function converts a vector to an oParameters struct */
void VectorToParameters(SParametersType &oParameters, CVector &vdTheta) {
  for(int i=0; i<vdTheta.m_pnDimSize[0]; i++)
    oParameters.pdPar[i] = double(vdTheta[i]);
} // VectorToParameters

bool SourceWorksOnResidual(SConstantsType &oConstants, int s, int a, int i) {
  bool bReturnValue = false;

  if(oConstants.poSources[s].pbScope[a] == true)
    bReturnValue = true;

  if(oConstants.bUseInstrTags) {
    if(oConstants.poPulsars[a].pbFlagSet[i]) {
      if(strcmp(oConstants.poSources[s].strScopeTag, oConstants.poPulsars[a].pstrFlags[i])==0) {
	bReturnValue = true;
      } // if pstrFlags
    } // if pbFlagSet
  } // if bUseInstrTags

  return bReturnValue;
} // SourceWorksOnResidual

bool SourceWorksOnPulsar(SConstantsType &oConstants, int s, int a) {
  bool bReturnValue = false;

  if(oConstants.poSources[s].pbScope[a] == true) {
    bReturnValue = true;
  }

  if(oConstants.bUseInstrTags) {
    for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
      if(oConstants.poPulsars[a].pbFlagSet[i] == true) {
	if(strcmp(oConstants.poSources[s].strScopeTag, oConstants.poPulsars[a].pstrFlags[i])==0) {
	  bReturnValue = true;
	  break;
	} // if pstrFlags
      }	// if pbFlagSet
    } // for i
  } // if bUseInstrTags

  return bReturnValue;
} // SourceWorksOnPulsar

bool SourceWorks(SConstantsType &oConstants, int s) {
  bool bReturnValue = false;

  for(int a=0; a<oConstants.k; a++) {
    if(oConstants.poSources[s].pbScope[a] == true) {
      bReturnValue = true;
      return true;
    } // if pbScope

    if(oConstants.bUseInstrTags) {
      for(int i=0; i<oConstants.poPulsars[a].nObservations; i++) {
	if(oConstants.poPulsars[a].pbFlagSet[i]) {
	  if(strcmp(oConstants.poSources[s].strScopeTag, oConstants.poPulsars[a].pstrFlags[i])==0) {
	    bReturnValue = true;
	    return true;
	  } // if pstrFlags
	} // if pbFlagSet
      } // for i
    } // if bUseInstrTags
  } // for a

  return bReturnValue;
} // SourceWorks

/* This function returns the number of parameters that is varying in the
 * sampling (e.g. the MCMC)
 * */
int NumberOfVaryingParameters(SConstantsType &oConstants) {
  int nVaryingParameters=0;

  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(SourceWorks(oConstants, s) && (oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p))) {
	nVaryingParameters++;
      } // if nFitTag
    } // for p
  } // for s

  return nVaryingParameters;
} // NumberOfVaryingParameters

bool IsAmplitudeParameter(SConstantsType &oConstants, int nP) {
  int nIndex=0;
  bool bReturnValue=false;

  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(SourceWorks(oConstants, s) && (oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p))) {
	if(nIndex == nP && p==0) {
	  bReturnValue = true;
	} // if nIndex
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
  return bReturnValue;
} // IsAmplitudeParameter

/* The parameters as in SParametersType oParameters objects are not numbered
 * according to varying parameters in an MCMC. In applications, like the
 * python wrappers, there therefore needs to be a mapping from varying
 * parameters to oParameters objects. This function creates that mapping
 */
void SetParameterIndexArray(SConstantsType &oConstants, int *pnParameterIndex) {
  int nIndex=0;

  for(int s=0; s<oConstants.nSources; s++) {
    for(int p=0; p<oConstants.poSources[s].oSourceType.nParameters; p++) {
      if(SourceWorks(oConstants, s) && (oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p))) {
	pnParameterIndex[nIndex] = oConstants.poSources[s].nFirstParIndex+p;
	nIndex++;
      } // if nFitTag
    } // for p
  } // for s
} // SetIndexArray

double GenerateValidProposal(double dValue, double dWidth, double dMin, double dMax, gsl_rng *rng, bool bWrap) {
  double dProposal;

  for(;;) { // Keep checking for a good value
    dProposal = dValue + dWidth * gsl_ran_ugaussian(rng);

    if(bWrap) {
      if(dProposal < dMin)
	dProposal = dMax - (dMin - dProposal);
      if(dProposal > dMax)
	dProposal = dMin + (dProposal - dMax);
    } // if bWrap

    if(dProposal > dMin && dProposal < dMax) {
      return dProposal;
    } // if dProposal
  } // for

  return 0;
} // GenerateValidProposal

/* This function is a template for how an arbitrary correlation function can be
 * implemented in the code. Requirements:
 * The SID code is set to SID_DipoleGWB (or something compatible)
 * The CID code is set to CID_Uniform (so that the correlation is always used)
 *
 * Input:
 * vdParameters: CVector object that contains the parameters of the signal
 *
 * Through oConstants, available:
 * dPhi1, dPhi2: The RA angle of both pulsars in radians (azimuthal angle)
 * dTheta1, dTheta2: pi/2 - declination (radians) (polar angle)
 *
 * nMode: 0 = 0, 0
 *        1 = 1, -1
 *        2 = 1, 0
 *        3 = 1, 1
 *
 * return: the correlation coefficient
 */
double AnisotropicGWBCorrelation(int np1, int np2, int nMode, SConstantsType &oConstants, CVector &vdParameters) {
    double dRz, dCorr, dPhi1, dPhi2, dTheta1, dTheta2;

    dPhi1 = oConstants.poPulsars[np1].dPhi;
    dPhi2 = oConstants.poPulsars[np2].dPhi;
    dTheta1 = oConstants.poPulsars[np1].dTheta;
    dTheta2 = oConstants.poPulsars[np2].dTheta;

    // Example: calculate the H&D correlation curve
    dRz = sin(dTheta1)*sin(dTheta2)*cos(dPhi1-dPhi2)+cos(dTheta1)*cos(dTheta2);

    // Here's an example of calculating the H&D coefficient. This one of course
    // does not depend on extra parameters. If the signal is a power-law (3
    // parameters), then the extra parameters can be found on spot 3 and
    // further...
    // double extrapar1 = vdParameters[3];
    // double extrapar2 = vdParameters[4];
    //
    // now extrapar1 and extrapar2 can be used to evaluate the correlation coeff
    //
    // The length of the vdParameters vector can be checked with:
    // len = vdParameters.m_pnDimSize[0]
    //
    // You can check that the vector exists with: vdParameters.Defined()

    // Do this check to prevent non-possible log evaluations
    if(dRz >= 1) {
        dCorr = 0.5;
    } else {
        dCorr = 0.75*(1-dRz)*log(0.5*(1-dRz))-(1.0/8.0)*(1-dRz)+0.5;
    } // if dRz

    // Do we need to add an auto-correlation term?
    if(np1 == np2) {
        dCorr += 0.5;
    }

    return dCorr;
} // AnisotropicGWBCorrelation
