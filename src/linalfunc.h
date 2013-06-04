/* linalfunc.h -- extension functions for the linear algebra classes

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

// We want to be able to use the linal classes
#include "linal.h"

#ifndef __LINALFUNC_H__
#define __LINALFUNC_H__

// Write data to file so it can be plotted
void WritePlot(const char *strFileName, CVector &vecx);
void WritePlot(const char *strFileName, CVector &vecx, CVector &vecy);
void WritePlot(const char *strFileName, CVector &vecx, CVector &vecy, CVector &vecz);
void WritePlot(const char *strFileName, CVector &vecx, CVector &vec1, CVector &vec2, CVector &vec3);
void WritePlot(const char *strFileName, CVector &vecx, CVector &vec1, CVector &vec2, CVector &vec3, CVector &vec4);
void WritePlot(const char *strFileName, CVector &vecx, CVector &vec1, CVector &vec2, CVector &vec3, CVector &vec4, CVector &vec5);
bool WriteMatrix(const char *strFileName, CMatrix &mat);

double InvErf(double d);                      // The inverse error function
CNumber &InvErf(CNumber &num);
CVector &InvErf(CVector &vec);
CMatrix &InvErf(CMatrix &mat);

double Sqrt(double d);				// The square root function 
CNumber &Sqrt(CNumber &num);
CVector &Sqrt(CVector &vec);
CMatrix &Sqrt(CMatrix &mat);

// This is quick 'n dirty. Extent to general case
CVector &ArcSin(CVector &vec);
CVector &Sin(CVector &vec);
CVector &Fabs(CVector &vec);
CVector &Log(CVector &vec);
CVector &Exp(CVector &vec);
CVector &Exp_1(CVector &vec);		// Means Exp(-vec)

double Min(CVector &vec);
double Max(CVector &vec);
double Min(CMatrix &mat);
double Max(CMatrix &mat);
double Sum(CVector &vec);
double Sum(CMatrix &mat);
double Prod(CVector &vec);
double Prod(CMatrix &mat);
double Trace(CMatrix &mat);
double Det(CMatrix &mat);
double LogDetChol(CMatrix &mat);	// Log(Det(mat))   (= Trace(Log(mat)) )   <- Use Cholesky
double LogDetCholFloat(CMatrix &mat);	// Log(Det(mat))   (= Trace(Log(mat)) )   <- Use Cholesky
double LogDet(CMatrix &mat);		// Log(Det(mat))   (= Trace(Log(mat)) )   <- Don't use Cholesky

// GMP versions of Cholesky...
CVector &CholeskyGmp(CMatrix &mat, int nBits);
double LogDetCholGmp(CMatrix &mat, int nBits);

// Remove a polynomial from the data
CVector &PolynomialRemoved(CVector &vdX, CVector &vdY, int nOrder);


// For debugging purposes
void PrintMatrix(CMatrix &mat);
void PrintVector(CVector &vec);

// __LINALFUNC_H__
#endif
