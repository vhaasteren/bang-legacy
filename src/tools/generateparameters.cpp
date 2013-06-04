/* generateparameters.cpp -- Generate parameters.conf automatically

   Rutger van Haasteren 31 Januari 2008 haasteren@strw.leidenuniv.nl

   Copyright (C) 2008 Rutger van Haasteren.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main(int argc, char **argv) {
  FILE *pFile;
  bool bSetK=false, bSetL=false, bSetT=false, bSetN=false,
       bSetI=false, bSetO=false, bSetA=false, bSetG=false,
       bSetD=false, ba=false, br=false, bSetn=false, bSetB=false,
       bSetM=false;
  int nGWBSpectrum=1;
  double dK=10, dL=100, dI=0.025, dT=550, dN=100,
	 dA=1, dG=2.330001;
  char cPreviousOption, strFileName[160], strDataDir[160],
       strResidualsFile[160], strAnglesFile[160];

  strcpy(strFileName, "parameters.conf");
  strcpy(strDataDir, "../data/workdata/");
  strcpy(strResidualsFile, "resdiduals.dat");
  strcpy(strAnglesFile, "angles.dat");

  if(argc < 2) {
    printf("Usage is: generateparameters [-switches]\n");
    printf("\n");
    printf("               Switches: -i X  Set parameter i to X (delta t) DEPRECATED!!!\n");
    printf("                         -k X  Set parameter k to X (amount of pulsars)\n");
    printf("                         -l X  Set parameter l to X (datapoints per pulsar)\n");
    printf("                         -n X  Set parameter n to X (GWB spectrum)\n");
    printf("                         -N X  Set parameter N to X (pulsar timing noise)\n");
    printf("                         -M X  Set parameter N to 50*(1.2)^X (pulsar timing noise)\n");
    printf("                         -T X  Set parameter T to X (duration experiment)\n");
    printf("                         -A X  Set parameter A to X (GWB amplitude)\n");
    printf("                         -B X  Set parameter A to (1.3)^X/100 (GWB amplitude)\n");
    printf("                         -G X  Set parameter G to X (GWB exponent)\n");
    printf("                         -o X  Write to file X\n");
    printf("                         -d X  Set datadir to X\n");
    printf("                         -r X  Set residuals file to X\n");
    printf("                         -a X  Set angles file to X\n");
    printf("\n");
    return 1;
  } // if argc


  for(int i=1; i<argc; i++) {
    if(argv[i][0] == '-') {
      // process the switches
      for(int j=1; j<strlen(argv[i]); j++) {
        switch(argv[i][j]) {
        case 'A':
          bSetA=true;
          break;
	case 'B':
	  bSetB=true;
	  break;
        case 'G':
          bSetG=true;
          break;
        case 'i':
          bSetI=true;
          break;
        case 'k':
          bSetK=true;
          break;
        case 'l':
	  bSetL=true;
          break;
        case 'n':
	  bSetn=true;
          break;
        case 'N':
	  bSetN=true;
          break;
        case 'M':
	  bSetM=true;
          break;
        case 'T':
	  bSetT=true;
          break;
	case 'o':
	  bSetO=true;
	  break;
	case 'd':
	  bSetD=true;
	  break;
	case 'a':
	  ba=true;
	  break;
	case 'r':
	  br=true;
	  break;
        default:
	  printf("Option %s not recognized!\n", &argv[i][j]);
          break;
        } // switch
	cPreviousOption = argv[i][j];
      } // for j
    } else {
      // Set computer name
      switch(cPreviousOption) {
	case 'A':
	  dA = atof(argv[i]);
	  break;
	case 'B':
	  dA = 0.01*pow(1.2, atof(argv[i]));
	  break;
	case 'G':
	  dG = atof(argv[i]);
	  break;
	case 'i':
	  dI = atof(argv[i]);
	  break;
	case 'k':
	  dK = atof(argv[i]);
	  break;
	case 'l':
	  dL = atof(argv[i]);
	  break;
        case 'n':
	  nGWBSpectrum = atoi(argv[i]);
          break;
	case 'N':
	  dN = atof(argv[i]);
	  break;
	case 'M':
	  dN = 50*pow(1.3, atof(argv[i]));
	  break;
	case 'T':
	  dT = atof(argv[i]);
	  break;
	case 'o':
	  // Use this argument as a new data-dir
	  strcpy(strFileName, argv[i]); 
	  break;
	case 'd':
	  // Use this argument as a new data-dir
	  strcpy(strDataDir, argv[i]); 
	  break;
	case 'a':
	  // Use this argument as a new data-dir
	  strcpy(strAnglesFile, argv[i]); 
	  break;
	case 'r':
	  // Use this argument as a new data-dir
	  strcpy(strResidualsFile, argv[i]); 
	  break;
	case ' ':
	default:
	  break;
      } // switch
    } // if argv
  } // for i



  if(! (pFile = fopen(strFileName, "w")) ) throw 1;

  fprintf(pFile, "# parameters.conf -- Parameter values for the pta project\n");
  fprintf(pFile, "# 31 Januari 2008 haasteren@strw.leidenuniv.nl\n");
  fprintf(pFile, "# Copyright (C) 2006-2008 Rutger van Haasteren.\n");
  fprintf(pFile, "# \n");
  fprintf(pFile, "# Automatically generated by generateparameters\n");
  fprintf(pFile, "# \n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Versions of pta that work with this data-file\n");
  fprintf(pFile, "VERSION = 0.86\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Basic dimensions of the pta analysis\n");
  fprintf(pFile, "# TODO: Add possibility of various data-points per pulsar, and unevenly sampled\n");
  fprintf(pFile, "# data\n");
  fprintf(pFile, "k = %i 				# Number of pulsars\n", int(dK));
  if(bSetI && bSetT)
    dL = (dT * 7) / (dI);
  fprintf(pFile, "l = %i				# Number of datapoints per pulsar\n", int(dL));
  fprintf(pFile, "m = 1				# Number of datasets (testing purposes)\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Measuring times are randomized. What should be the average interval?\n");
  fprintf(pFile, "# TODO: Remove the need for this\n");

  if(bSetI) {
    fprintf(pFile, "dMeasureTimeInterval = %f	# Interval between measurements ([] = 2yr)\n", (dI / (730.0)) );
  } else {
    fprintf(pFile, "dMeasureTimeInterval = %f	# Interval between measurements ([] = 2yr)\n", (dT/ (104 * dL)) );
  } // if bSetI

  fprintf(pFile, "\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Global system parameters\n");
  fprintf(pFile, "# GWB spectrum parameters\n");
  fprintf(pFile, "nGWBSpectrum = %i		# Spectrum of GWB\n", nGWBSpectrum);
  fprintf(pFile, "dAmp = %f			# Amplitude of GW  ([] = 1E-15 yr^(1/2) )\n", dA);
  fprintf(pFile, "dExp = %f			# Exponent of GW Spectrum\n", dG);
  fprintf(pFile, "dCutoffFrequency = 0.05  	# The cutoff frequency ([] = 0.5 Yr^-1) (100 weeks)\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Extra parameters of the GWB\n");
  fprintf(pFile, "dAmpMinBound = %f\n", dA * 0.01);
  fprintf(pFile, "dAmpMaxBound = %f\n", dA * 20);
  fprintf(pFile, "dAmpStart = 0.8\n");
  fprintf(pFile, "dAmpWidthMCMC = 0.02\n");
  fprintf(pFile, "dAmpWidthFit = 0.01\n");
  fprintf(pFile, "dExpMinBound = 1.01\n");
  fprintf(pFile, "dExpMaxBound = 5.01\n");
  fprintf(pFile, "dExpStart = 2.23001\n");
  fprintf(pFile, "dExpWidthMCMC = 0.02\n");
  fprintf(pFile, "dExpWidthFit = 0.01\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Deterministic sources\n");
  fprintf(pFile, "# TODO: Generalize this\n");
  fprintf(pFile, "# Pulsars have quadratic spindown. This effect is randomized with this gauge\n");
  fprintf(pFile, "nSpindownStrength = 0 		# Quadratic spindown strength\n");
  fprintf(pFile, "nCorrectOrder = 3		# Number of correction terms\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Pulsar parameters can be specified per pulsar, or by default. If for a\n");
  fprintf(pFile, "# particular pulsar no parameters are specified, the default values will apply.\n");
  fprintf(pFile, "# TODO: Read pulsar-parameter-data from a data-file\n");
  fprintf(pFile, "#\n");
  fprintf(pFile, "# The following spectra are available:\n");
  fprintf(pFile, "# 0: White noise\n");
  fprintf(pFile, "# 1: Power law\n");
  fprintf(pFile, "# 2: Exponential\n");
  fprintf(pFile, "# 3: Lorentzian (0 = amp, 1 = f0)\n");
  fprintf(pFile, "# 4: White noise + Lorentzian (0 = wnamp, 1 = lamp, 2 = f0)\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Default pulsar parameters\n");
  fprintf(pFile, "nDefaultParametersPerPulsar = 1	# Amount of parameters per pulsar\n");
  fprintf(pFile, "nDefaultPulsarSpectrum = 0	# White Noise\n");
  fprintf(pFile, "# nDefaultMeasurementsPerPulsar = 100 # Amount of datapoints/measurements per pulsar (= 'l') Not used *YET*\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "dDefaultPar[0] = %f		# First parameter (pulsar noise amplitude)\n", dN);
  fprintf(pFile, "dDefaultParMinBound[0] = 0.01	# '' '' minimum bound for searching for max\n");
  fprintf(pFile, "dDefaultParMaxBound[0] = 990	# '' '' maximum bound for searching for max\n");
  fprintf(pFile, "dDefaultParStart[0] = 80	# '' '' start point for ML-search and MCMC\n");
  fprintf(pFile, "dDefaultParWidthMCMC[0] = 2.5	# '' '' typical width (step-size of MCMC)\n");
  fprintf(pFile, "dDefaultParWidthFit[0] = 2.5	# '' '' typical width (ML Gaussian fit width)\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "dDefaultPar[1] = 1		# Second parameter (pulsar noise exponent/amplitude)\n");
  fprintf(pFile, "dDefaultParMinBound[1] = 0.1	# '' '' minimum bound for searching for max\n");
  fprintf(pFile, "dDefaultParMaxBound[1] = 10	# '' '' maximum bound for searching for max\n");
  fprintf(pFile, "dDefaultParStart[1] = 0.7	# '' '' start point for ML-search and MCMC\n");
  fprintf(pFile, "dDefaultParWidthMCMC[1] = 0.03	# '' '' typical width (step-size of MCMC)\n");
  fprintf(pFile, "dDefaultParWidthFit[1] = 0.03	# '' '' typical width (ML Gaussian fit width)\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "dDefaultPar[2] = 1		# Third parameter (pulsar noise exponent)\n");
  fprintf(pFile, "dDefaultParMinBound[2] = 0.1	# '' '' minimum bound for searching for max\n");
  fprintf(pFile, "dDefaultParMaxBound[2] = 10	# '' '' maximum bound for searching for max\n");
  fprintf(pFile, "dDefaultParStart[2] = 0.7	# '' '' start point for ML-search and MCMC\n");
  fprintf(pFile, "dDefaultParWidthMCMC[2] = 0.03	# '' '' typical width (step-size of MCMC)\n");
  fprintf(pFile, "dDefaultParWidthFit[2] = 0.03	# '' '' typical width (ML Gaussian fit width)\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Pulsar parameters can be provided explicitly by leaving out the 'Default' and\n");
  fprintf(pFile, "# providing an extra [x] pulsar-number before the paramter-number.\n");
  fprintf(pFile, "#\n");
  fprintf(pFile, "# Example (for use, replace 'x' with the pulsar-number)\n");
  fprintf(pFile, "# nParametersPerPulsar[x] = 2	# Amount of parameters per pulsar\n");
  fprintf(pFile, "# nPulsarSpectrum[x] = 0	# White noise\n");
  fprintf(pFile, "# dPar[x][0] = 100		# First parameter (pulsar noise amplitude)\n");
  fprintf(pFile, "# dParMinBound[x][0] = 0.01	# '' '' minimum bound for searching for max\n");
  fprintf(pFile, "# dParMaxBound[x][0] = 990	# '' '' maximum bound for searching for max\n");
  fprintf(pFile, "# dParStart[x][0] = 80	# '' '' start point for ML-search and MCMC\n");
  fprintf(pFile, "# dParWidthMCMC[x][0] = 2.5	# '' '' typical width (step-size of MCMC)\n");
  fprintf(pFile, "# dParWidthFit[x][0] = 2.5	# '' '' typical width (ML Gaussian fit width)\n");
  fprintf(pFile, "# dPar[x][1] = 1		# First parameter (pulsar noise amplitude)\n");
  fprintf(pFile, "# dParMinBound[x][1] = 0.1	# '' '' minimum bound for searching for max\n");
  fprintf(pFile, "# dParMaxBound[x][1] = 10	# '' '' maximum bound for searching for max\n");
  fprintf(pFile, "# dParStart[x][1] = 0.7	# '' '' start point for ML-search and MCMC\n");
  fprintf(pFile, "# dParWidthMCMC[x][1] = 0.03	# '' '' typical width (step-size of MCMC)\n");
  fprintf(pFile, "# dParWidthFit[x][1] = 0.03	# '' '' typical width (ML Gaussian fit width)\n");
  fprintf(pFile, "#\n");
  fprintf(pFile, "# ---------------------------------\n");
  fprintf(pFile, "# The parameters for pulsar 1\n");
  fprintf(pFile, "# ---------------------------------\n");
  fprintf(pFile, "# Write some explicit parameters here... \n");
  fprintf(pFile, "\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# General MC parameters\n");
  fprintf(pFile, "# TODO: Implement acceptance-rate + ETA\n");
  fprintf(pFile, "nSmallPlotPoints = 200		# Number of points used in plot when plotting\n");
  fprintf(pFile, "nPlotPoints = 2000		# Number of points used in plot when plotting\n");
  fprintf(pFile, "nMCMCSteps = 65000		# Number of mcmc steps used in the calculation\n");
  fprintf(pFile, "nBurnInSteps = 100		# Number of mcmc steps (burn-in time)\n");
  fprintf(pFile, "nAcceptanceRate = 24		# Acceptance rate in (%) of new points while sampling\n");
  fprintf(pFile, "nBootstrapAttempts=1000		# Amount of trials in the bootstrap error-estimation\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# The amount of bits used in the Cholesky decomposition (GMP lib) (deprecated)\n");
  fprintf(pFile, "# TODO: Remove this\n");
  fprintf(pFile, "nCholeskyBits = 64		# Amount of bits used in calculation\n");
  fprintf(pFile, "\n");
  fprintf(pFile, "# Location of all the data-files generated by the program-chain where all the\n");
  fprintf(pFile, "# data-files are stored\n");
  fprintf(pFile, "# This file (parameters.conf) is located at the super-directory of strDataDir\n");
  fprintf(pFile, "strDataDir = %s\n", strDataDir);
  fprintf(pFile, "strResidualsFile = %s\n", strResidualsFile);
  fprintf(pFile, "strAnglesFile = %s\n", strAnglesFile);
  fprintf(pFile, "\n");

  fclose(pFile);
  return 0;
} // main
