/* banfunc.h -- extension functions for the pta program ban

   Rutger van Haasteren 15 August 2007 haasteren@strw.leidenuniv.nl

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


// Minimum version of the configfile that is compatible with this release
#define CONFIG_VERSION_MIN	0.20


// We want to be able to use the linal classes
#include "linal.h"
#include "config.h"

#ifndef __BANFUNC_H__
#define __BANFUNC_H__


// Check were the definition of pi really is
#ifndef PI
  #define PI	3.14159265358979323846264338327950288419716939937510
#endif

#ifndef NSPERYEAR
  #define NSPERYEAR	31557600000000000.0
  #define SPERYEAR	31557600.0
#endif

// The dispersion constant in seconds
#ifndef DM_K
//  #define DM_K 4.15e-3  // Units of Kj's paper (check these)
  #define DM_K   241.0    // GHz^-2 cm^-3 pc s^-1
#endif


int KbHit();

bool BANFuncInitialize(int nProcessRank=0);                   // Initialise some constantes, functions and random numbers

void InitProgressBar(const char *strMsg=NULL);
void DrawProgressBar(int nPercent=0);
void FinishProgressBar();

// For status report
void PrintStatus(const char *str);
void PrintFailed();
void PrintSuccess();

// For (MCMC) updates
void PrintUpdate(const char *str);
void FinishUpdate();

int n_min(int *pnArray, int nLength);
double d_min(double *pdArray, int nLength);
void d_swap(double &x, double &y);
void QuickSort(double *pdBufferX, double *pdBufferY, double *pdBufferZ, double *pdBufferComp, int left, int right);
void QuickSort(double *pdBuffer, int left, int right);

void QuickSort(double **ppdBuffer, double *pdBufferZ, int nDimensions, int left, int right);

// __BANFUNC_H__
#endif
