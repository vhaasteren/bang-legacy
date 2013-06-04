/* linalfunc.cpp -- extension functions for the linear algebra classes

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


#include "config.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "moremath.h"
#include "linalfunc.h"


// Write data to datafile
void WritePlot(const char *strFileName, CVector &vecx) {
  FILE *file;
  if(! ( vecx.Defined() )) {throw ELENotDefined; }

  if(! (file=fopen(strFileName, "w+")) ) {
    // unable to open file
    printf("Unable to open file: %s\n", strFileName);
    return;
  } // if
  fprintf(file, "# File %s automatically created by linalfunc\n", strFileName);
  fprintf(file, "# Written by Rutger van Haasteren\n");
//              "12345678901234567890123456789012345678901234567890
  for(int i=0; i<vecx.m_pnDimSize[0]; i++) {
    fprintf(file, "%.18e \n", double(vecx[i]));
  } // for i

  fclose(file);
  return;
} // WritePlot

// Write data to datafile
void WritePlot(const char *strFileName, CVector &vecx, CVector &vecy) {
  FILE *file;
  if(! ( vecx.Defined() && vecy.Defined() )) {throw ELENotDefined; }
  if(vecx.m_pnDimSize[0] != vecy.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  if(! (file=fopen(strFileName, "w+")) ) {
    // unable to open file
    printf("Unable to open file: %s\n", strFileName);
    return;
  } // if
  fprintf(file, "# File %s automatically created by linalfunc\n", strFileName);
  fprintf(file, "# Written by Rutger van Haasteren\n");
//              "12345678901234567890123456789012345678901234567890
  for(int i=0; i<vecx.m_pnDimSize[0]; i++) {
    fprintf(file, "%.18e %.18e \n", double(vecx[i]), double(vecy[i]));
  } // for i

  fclose(file);
  return;
} // WritePlot

// Write data to datafile
void WritePlot(const char *strFileName, CVector &vecx, CVector &vecy, CVector &vecz) {
  FILE *file;
  if(! ( vecx.Defined() && vecy.Defined() && vecz.Defined())) {throw ELENotDefined; }
  if(vecx.m_pnDimSize[0] != vecy.m_pnDimSize[0] || vecx.m_pnDimSize[0] != vecz.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  if(! (file=fopen(strFileName, "w+")) ) {
    // unable to open file
    printf("Unable to open file: %s\n", strFileName);
    return;
  } // if
  fprintf(file, "# File %s automatically created by linalfunc\n", strFileName);
  fprintf(file, "# Written by Rutger van Haasteren\n");
//              "12345678901234567890123456789012345678901234567890
  for(int i=0; i<vecx.m_pnDimSize[0]; i++) {
    fprintf(file, "%.18e %.18e %.18e\n", double(vecx[i]), double(vecy[i]), double(vecz[i]));
  } // for i

  fclose(file);
  return;
} // WritePlot

void WritePlot(const char *strFileName, CVector &vecx, CVector &vec1, CVector &vec2, CVector &vec3) {
  FILE *file;
  if(! ( vecx.Defined() && vec1.Defined() && vec2.Defined() && vec3.Defined())) {throw ELENotDefined; }
  if(vecx.m_pnDimSize[0] != vec1.m_pnDimSize[0] || vecx.m_pnDimSize[0] != vec2.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }
  if(vecx.m_pnDimSize[0] != vec3.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }
  if(! (file=fopen(strFileName, "w+")) ) {
    // unable to open file
    printf("Unable to open file: %s\n", strFileName);
    return;
  } // if
  fprintf(file, "# File %s automatically created by linalfunc\n", strFileName);
  fprintf(file, "# Written by Rutger van Haasteren\n");
//              "12345678901234567890123456789012345678901234567890123456789012
  for(int i=0; i<vecx.m_pnDimSize[0]; i++) {
    fprintf(file, "%.18e %.18e %.18e %.18e\n", double(vecx[i]), double(vec1[i]), double(vec2[i]), double(vec3[i]));
  } // for i

  fclose(file);
  return;
} // WritePlot

void WritePlot(const char *strFileName, CVector &vecx, CVector &vec1, CVector &vec2, CVector &vec3, CVector &vec4) {
  FILE *file;
  if(! ( vecx.Defined() && vec1.Defined() && vec2.Defined() && vec3.Defined() && vec4.Defined())) {throw ELENotDefined; }
  if(vecx.m_pnDimSize[0] != vec1.m_pnDimSize[0] || vecx.m_pnDimSize[0] != vec2.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }
  if(vecx.m_pnDimSize[0] != vec3.m_pnDimSize[0] || vecx.m_pnDimSize[0] != vec4.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }
  if(! (file=fopen(strFileName, "w+")) ) {
    // unable to open file
    printf("Unable to open file: %s\n", strFileName);
    return;
  } // if
  fprintf(file, "# File %s automatically created by linalfunc\n", strFileName);
  fprintf(file, "# Written by Rutger van Haasteren\n");
//              "12345678901234567890123456789012345678901234567890123456789012
  for(int i=0; i<vecx.m_pnDimSize[0]; i++) {
    fprintf(file, "%.18e %.18e %.18e %.18e %.18e\n", double(vecx[i]), double(vec1[i]), double(vec2[i]), double(vec3[i]), double(vec4[i]));
  } // for i

  fclose(file);
  return;
} // WritePlot

void WritePlot(const char *strFileName, CVector &vecx, CVector &vec1, CVector &vec2, CVector &vec3, CVector &vec4, CVector &vec5) {
  FILE *file;
  if(! ( vecx.Defined() && vec1.Defined() && vec2.Defined() && vec3.Defined() && vec4.Defined() && vec5.Defined())) {throw ELENotDefined; }
  if(vecx.m_pnDimSize[0] != vec1.m_pnDimSize[0] || vecx.m_pnDimSize[0] != vec2.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }
  if(vecx.m_pnDimSize[0] != vec3.m_pnDimSize[0] || vecx.m_pnDimSize[0] != vec4.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }
  if(vecx.m_pnDimSize[0] != vec5.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  if(! (file=fopen(strFileName, "w+")) ) {
    // unable to open file
    printf("Unable to open file: %s\n", strFileName);
    return;
  } // if
  fprintf(file, "# File %s automatically created by linalfunc\n", strFileName);
  fprintf(file, "# Written by Rutger van Haasteren\n");
//              "12345678901234567890123456789012345678901234567890123456789012
  for(int i=0; i<vecx.m_pnDimSize[0]; i++) {
    fprintf(file, "%.18e %.18e %.18e %.18e %.18e %.18e\n", double(vecx[i]), double(vec1[i]), double(vec2[i]), double(vec3[i]), double(vec4[i]), double(vec5[i]));
  } // for i

  fclose(file);
  return;
} // WritePlot


// Write a whole matrix to a datafile (dat, tab-separated values)
bool WriteMatrix(const char *strFileName, CMatrix &mat) {
  FILE *file;
  if(! mat.Defined() ) {throw ELENotDefined; }

  if(! (file=fopen(strFileName, "w+")) ) {
    // unable to open file
//    printf("Unable to open file: %s\n", strFileName);
    return false;
  } // if

//  fprintf(file, "# File %s automatically created by linalfunc\n", strFileName);
//  fprintf(file, "# Written by Rutger van Haasteren\n");
  for(int i=0; i<mat.m_pnDimSize[0]; i++) {
    for(int j=0; j<mat.m_pnDimSize[1]; j++) {
      if(! fprintf(file, "%.30e ", double(mat[i][j]) ) ) return false;
      if(j!=mat.m_pnDimSize[1]-1) {
        if(! fprintf(file, "\t") ) return false;
      }
    } // for j
    fprintf(file, "\n");
  } // for i

  if( fclose(file) ) return false;
  return true;
} // WriteMatrix


double InvErf(double d) { return inverf(d); } // InvErf

// De inverse error function for numbers
CNumber &InvErf(CNumber &num) {
  CNumber *pnumInvErf;
  if(! num.Defined()) {throw ELENotDefined; }
  
  pnumInvErf = new CNumber(InvErf(double(num)));	// Make a deep copy of the original
  pnumInvErf->SetAllocated(true);		// We allocated the memory space

  if(num.Allocated()) {delete &num; } // This should call for the destructor
  return *pnumInvErf;
} // InvErf


// De inverse error function for vectors
CVector &InvErf(CVector &vec) {
  CVector *pvecInvErf;
  if(! vec.Defined()) {throw ELENotDefined; }
  
  pvecInvErf = new CVector(vec.m_pnDimSize[0]); // Make a deep copy of the original
  pvecInvErf->SetAllocated(true);		// We allocated the memory space
  
  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    (*pvecInvErf)[i] = InvErf(double(vec[i]));
  }

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return *pvecInvErf;
} // InvErf

// De inverse error function for matrices
CMatrix &InvErf(CMatrix &mat) {
  CMatrix *pmatInvErf;
  if(! mat.Defined()) {throw ELENotDefined; }
  
  pmatInvErf = new CMatrix(mat.m_pnDimSize[0], mat.m_pnDimSize[1]);		// Make a deep copy of the original
  pmatInvErf->SetAllocated(true);		// We allocated the memory space
  
  for(int i=0; i<mat.m_pnDimSize[0]; i++) {
    for(int j=0; j<mat.m_pnDimSize[1]; j++) {
      (*pmatInvErf)[i][j] = InvErf(double(mat[i][j]));
    }
  }

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return *pmatInvErf;
} // InvErf


// The square root function 
double Sqrt(double d) { return sqrt(d);}

// The square root function for numbers
CNumber &Sqrt(CNumber &num) {
  CNumber *pnumSqrt;
  if(! num.Defined()) {throw ELENotDefined; }
  
  pnumSqrt = new CNumber(sqrt(double(num)));	// Make a deep copy of the original
  pnumSqrt->SetAllocated(true);		// We allocated the memory space

  if(num.Allocated()) {delete &num; } // This should call for the destructor
  return *pnumSqrt;
} // Sqrt

// The square root function for vectors
CVector &Sqrt(CVector &vec) {
  CVector *pvecSqrt;
  if(! vec.Defined()) {throw ELENotDefined; }
  
  pvecSqrt = new CVector(vec.m_pnDimSize[0]);		// Make a deep copy of the original
  pvecSqrt->SetAllocated(true);		// We allocated the memory space
  
  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    (*pvecSqrt)[i] = sqrt(double(vec[i]));
  }

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return *pvecSqrt;
} // Sqrt

// The square root function for matrices
CMatrix &Sqrt(CMatrix &mat) {
  CMatrix *pmatSqrt;
  if(! mat.Defined()) {throw ELENotDefined; }
  
  pmatSqrt = new CMatrix(mat.m_pnDimSize[0],mat.m_pnDimSize[1]);		// Make a deep copy of the original
  pmatSqrt->SetAllocated(true);		// We allocated the memory space
  
  for(int i=0; i<mat.m_pnDimSize[0]; i++) {
    for(int j=0; j<mat.m_pnDimSize[1]; j++) {
      (*pmatSqrt)[i][j] = sqrt(double(mat[i][j]));
    }
  }

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return *pmatSqrt;
} // Sqrt

// The arcsin function for vectors
CVector &ArcSin(CVector &vec) {
  CVector *pvecArcSin;
  if(! vec.Defined()) {throw ELENotDefined; }
  pvecArcSin = new CVector(vec.m_pnDimSize[0]);			// Make a deep copy of the original
  pvecArcSin->SetAllocated(true);		// We allocated the memory space

  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    (*pvecArcSin)[i] = asin(double(vec[i]));
  }

  if(vec.Allocated()) { delete &vec; } // This should call for the destructor
  return *pvecArcSin;
} // ArcSin


// The sin function for vectors
CVector &Sin(CVector &vec) {
  CVector *pvecSin;
  if(! vec.Defined()) {throw ELENotDefined; }
  pvecSin = new CVector(vec.m_pnDimSize[0]);			// Make a deep copy of the original
  pvecSin->SetAllocated(true);		// We allocated the memory space

  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    (*pvecSin)[i] = sin(double(vec[i]));
  }

  if(vec.Allocated()) { delete &vec; } // This should call for the destructor
  return *pvecSin;
} // Sin


// The arcsin function for vectors
CVector &Fabs(CVector &vec) {
  CVector *pvecArcSin;
  if(! vec.Defined()) {throw ELENotDefined; }
  
  pvecArcSin = new CVector(vec.m_pnDimSize[0]);		// Make a deep copy of the original
  pvecArcSin->SetAllocated(true);		// We allocated the memory space

  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    (*pvecArcSin)[i] = fabs(double(vec[i]));
  }

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return *pvecArcSin;
} // Fabs

// The arcsin function for vectors
CVector &Log(CVector &vec) {
  CVector *pvecArcSin;
  if(! vec.Defined()) {throw ELENotDefined; }
  
  pvecArcSin = new CVector(vec.m_pnDimSize[0]);		// Make a deep copy of the original
  pvecArcSin->SetAllocated(true);		// We allocated the memory space
  
  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    (*pvecArcSin)[i] = log(double(vec[i]));
  }

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return *pvecArcSin;
} // Log

// The exp function for vectors
CVector &Exp(CVector &vec) {
  CVector *pvecExp;
  if(! vec.Defined()) {throw ELENotDefined; }
  
  pvecExp = new CVector(vec.m_pnDimSize[0]);		// Make a deep copy of the original
  pvecExp->SetAllocated(true);		// We allocated the memory space
  
  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    (*pvecExp)[i] = exp(double(vec[i]));
  }

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return *pvecExp;
} // Exp

// The 1/exp function for vectors
CVector &Exp_1(CVector &vec) {
  CVector *pvecExp;
  if(! vec.Defined()) {throw ELENotDefined; }
  
  pvecExp = new CVector(vec.m_pnDimSize[0]);		// Make a deep copy of the original
  pvecExp->SetAllocated(true);		// We allocated the memory space
  
  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    (*pvecExp)[i] = exp(0.0-double(vec[i]));
  }

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return *pvecExp;
} // Exp_1


double Sum(CVector &vec) {
  double dSum=0;
  if(! vec.Defined()) {throw ELENotDefined; }

  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    dSum += double(vec[i]);
  } // for i

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return dSum;
} // Sum

double Sum(CMatrix &mat) {
  double dSum=0;
  if(! mat.Defined()) {throw ELENotDefined; }

  for(int i=0; i<mat.m_pnDimSize[0]; i++) {
    for(int j=0; j<mat.m_pnDimSize[1]; j++) {
      dSum += double(mat[i][j]);
    } // for j
  } // for i

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dSum;
} // Sum

double Prod(CVector &vec) {
  double dProd=1;
  if(! vec.Defined()) {throw ELENotDefined; }

  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    dProd *= double(vec[i]);
  } // for i

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return dProd;
} // Prod

double Prod(CMatrix &mat) {
  double dProd=1;
  if(! mat.Defined()) {throw ELENotDefined; }

  for(int i=0; i<mat.m_pnDimSize[0]; i++) {
    for(int j=0; j<mat.m_pnDimSize[1]; j++) {
      dProd *= double(mat[i][j]);
    } // for j
  } // for i

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dProd;
} // Prod


double Trace(CMatrix &mat) {
  double dTrace=0;
  if(! mat.Defined()) {throw ELENotDefined; }
  if(mat.m_pnDimSize[0] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

  for(int i=0; i<mat.m_pnDimSize[0]; i++) {
    dTrace += double(mat[i][i]);
  } // for i

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dTrace;
} // Trace

double Det(CMatrix &mat) {
  double dDet=1;
  CVector vdEigenvalues;
  if(! mat.Defined()) {throw ELENotDefined; }
  if(mat.m_pnDimSize[0] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

  vdEigenvalues = mat.Eigenvalues();

  for(int i=0; i<vdEigenvalues.m_pnDimSize[0]; i++) {
    dDet *= double(vdEigenvalues[i]);
  } // for i

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dDet;
} // Det


// Given a positive-definite symmetric matrix
// a[1..n][1..n], this routine constructs its
// Cholesky decomposition, A = L · LT. On input,
// only the upper triangle of a need be given;
// it is not modified. The Cholesky factor L is
// returned in the lower triangle of a, except
// for its diagonal elements which are returned
// in p[1..n].
CVector &CholeskyGmp(CMatrix &mat, int nBits) {
  // Set the precision of the calculation:
/*  mpf_set_default_prec(nBits);		// Use 64 bit calculation
  mpf_class *mpfIn;			// The input matrix
  mpf_class *mpfOut;			// The output vector
  mpf_class mpfSum;			// Temporary variable
  CVector *vdOut;			// Output vector

  int i, j, k;				// Counting indices

  if(! mat.Defined()) {throw ELENotDefined; }
  if(mat.m_pnDimSize[0] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

  // Initialize the variables/matrices
  mpfIn = new mpf_class[int(mat.m_pnDimSize[0])*int(mat.m_pnDimSize[1])];
  mpfOut = new mpf_class[int(mat.m_pnDimSize[0])];
  for(i=0; i<mat.m_pnDimSize[0]; i++) {
    for(j=0; j<mat.m_pnDimSize[0]; j++) {
      mpfIn[i + mat.m_pnDimSize[0]*j] = double(mat[i][j]);
    } // for j
    mpfOut[i] = 0.0;
  } // for i
  mpfSum = 0.0;

//PrintMatrix(mat);

  for(i=0; i<mat.m_pnDimSize[0]; i++) {
    for(j=i; j<mat.m_pnDimSize[1]; j++) {
      mpfSum = mpfIn[i+mat.m_pnDimSize[0]*j];
//      printf("mpfSum = %e    ", mpfSum.get_d());
      for(k=i-1; k>=0; k--) {
        mpfSum -= mpfIn[i+mat.m_pnDimSize[0]*k]*mpfIn[j+mat.m_pnDimSize[0]*k];
//        printf("  mpfSum -= (mpfIn[%i+%i*%i]*mpfIn[%i+%i*%i] = mpfIn[%i][%i] = %e)\n", i, mat.m_pnDimSize[0], k, j, mat.m_pnDimSize[0], k, i, k, mpfIn[i+mat.m_pnDimSize[0]*k].get_d()*mpfIn[j+mat.m_pnDimSize[0]*k].get_d());
      } // for k
//      printf("mpfSum[%i][%i] = %e\n", i, j, mpfSum.get_d());
      if(i==j) {
        if(mpfSum <= 0.0) { 		// mpf_cmp_d(mpfSum, 0.0) < 0.0   ===    mpfSum < 0.0
          // The matrix wasn't positive definite according to Cholesky!
//          printf("Values:  i=%i  j=%i  \n", i, j);
          throw(ELELapack);
        } // if mpfSum < 0
        mpfOut[i] = sqrt(mpfSum);
      } else { // if i==j
        mpfIn[j + mat.m_pnDimSize[0]*i] = mpfSum / mpfOut[i];
      } // if i==j
    } // for j
  } // for i

  // Set the output variable
  vdOut = new CVector(mat.m_pnDimSize[0]);
  vdOut->SetAllocated(true);
  for(i=0; i<mat.m_pnDimSize[0]; i++) {
    (*vdOut)[i] = mpfOut[i].get_d();
  } // for i

  delete[] mpfIn;
  delete[] mpfOut;*/
  return (*(new CVector));
}  // choldc




// Log(Det(mat))  (=Trace(Log(mat)) ) <- Use Cholesky and gmp library!
double LogDetCholGmp(CMatrix &mat, int nBits){
  double dLogDet=0;
  CVector vdFact;
  if(! mat.Defined()) {throw ELENotDefined; }
  if(mat.m_pnDimSize[0] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

// Use Cholesky
  vdFact = CholeskyGmp(mat, nBits);

  for(int i =0; i<vdFact.m_pnDimSize[0]; i++) {
    dLogDet += log(double(vdFact[i]));
  } // for i
  dLogDet *= 2;

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dLogDet;
} // LogDetCholGmp


// Log(Det(mat))   (= Trace(Log(mat)) )  <- Use Cholesky
double LogDetChol(CMatrix &mat) {
  double dLogDet=0;
//  CVector vdEigenvalues;
  CMatrix mdFact;
  if(! mat.Defined()) {throw ELENotDefined; }
  if(mat.m_pnDimSize[0] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

// Use Cholesky
  mdFact = mat.Cholesky();

  for(int i =0; i<mdFact.m_pnDimSize[0]; i++) {
    if(double(mdFact[i][i]) > 0) dLogDet += log(double(mdFact[i][i]));
  } // for i
  dLogDet *= 2;

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dLogDet;
} // LogDetChol


// Log(Det(mat))   (= Trace(Log(mat)) )  <- Use Cholesky
double LogDetCholFloat(CMatrix &mat) {
  double dLogDet=0;
//  CVector vdEigenvalues;
  CMatrix mdFact;
  if(! mat.Defined()) {throw ELENotDefined; }
  if(mat.m_pnDimSize[0] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

// Use Cholesky
  mdFact = mat.CholeskySingle();

  for(int i =0; i<mdFact.m_pnDimSize[0]; i++) {
    if(double(mdFact[i][i]) > 0) dLogDet += log(double(mdFact[i][i]));
  } // for i
  dLogDet *= 2;

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dLogDet;
} // LogDet


// Log(Det(mat))   (= Trace(Log(mat)) )   <- Don't use Cholesky
/******************************** Careful! If matrix isn't pos.def.  ********************************/
/*********************************  This function is doesn't work!  *********************************/
double LogDet(CMatrix &mat) {
  double dLogDet=0;
//  CVector vdEigenvalues;
  CMatrix mdFact;
  if(! mat.Defined()) {throw ELENotDefined; }
  if(mat.m_pnDimSize[0] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

// Use Cholesky
  mdFact = mat.LUFact();

  for(int i =0; i<mdFact.m_pnDimSize[0]; i++) {
    if(double(mdFact[i][i]) > 0) dLogDet += log(fabs(double(mdFact[i][i])));
  } // for i
  dLogDet *= 2;

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dLogDet;
} // LogDet


double Min(CVector &vec) {
  double dMin;
  if(! vec.Defined()) {throw ELENotDefined; }
  
  dMin = vec[0];
  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    if(vec[i] < dMin) dMin = double(vec[i]);
  } // for i

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return dMin;
} // Min

double Max(CVector &vec) {
  double dMax;
  if(! vec.Defined()) {throw ELENotDefined; }
  
  dMax = vec[0];
  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    if(vec[i] > dMax) dMax = double(vec[i]);
  } // for i

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  return dMax;
} // Max

double Min(CMatrix &mat) {
  double dMin;
  if(! mat.Defined()) {throw ELENotDefined; }
  
  dMin = mat[0];
  for(int i=0; i<mat.m_pnDimSize[0]; i++) {
    if(mat[i] < dMin) dMin = double(mat[i]);
  } // for i

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dMin;
} // Min

double Max(CMatrix &mat) {
  double dMax;
  if(! mat.Defined()) {throw ELENotDefined; }
  
  dMax = mat[0];
  for(int i=0; i<mat.m_pnDimSize[0]; i++) {
    if(mat[i] > dMax) dMax = double(mat[i]);
  } // for i

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  return dMax;
} // Max


// Remove a polynomial from the data
CVector &PolynomialRemoved(CVector &vdX, CVector &vdY, int nOrder) {
  CVector *pvdPolyReduced, vdSubtract;
  CMatrix mdM, mdTemp;

  if(! vdX.Defined() || ! vdY.Defined()) {throw ELENotDefined; }
  if(vdX.m_pnDimSize[0] != vdY.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  pvdPolyReduced = new CVector(vdY.m_pnDimSize[0]); // Make a deep copy of the original
  pvdPolyReduced->SetAllocated(true);		// We allocated the memory space

  // Set the polynomial reduction matrix
  mdM.Initialize(vdX.m_pnDimSize[0], nOrder+1);
  for(int i=0; i<vdX.m_pnDimSize[0]; i++) {
    for(int j=0; j<=nOrder; j++) {
      mdM[i][j] = gsl_pow_int(double(vdX[i]), j);
    } // for j
  } // for i

  // Set the projection matrix
  mdTemp = mdM[LO_TRANSPOSE] * mdM;
  mdTemp.Invert();

  // Calculate the reduced vector
  vdSubtract = mdM * (mdTemp * (mdM[LO_TRANSPOSE] * vdY));
  *pvdPolyReduced = vdY - vdSubtract;

  if(vdX.Allocated()) { delete &vdX; } // This should call for the destructor
  if(vdY.Allocated()) { delete &vdY; } // This should call for the destructor
  return *pvdPolyReduced;
} // PolynomialRemoved



// For debugging purposes
void PrintMatrix(CMatrix &mat) {
  if(! mat.Defined()) {throw ELENotDefined; }

  putchar('\n');
  for(int i=0; i<mat.m_pnDimSize[0]; i++) {
    for(int j=0; j<mat.m_pnDimSize[1]; j++) {
      printf("%e  ", double(mat[i][j]));
    } // for j
    printf("\n");
  } // for i
} // PrintMatrix

void PrintVector(CVector &vec) {
  if(! vec.Defined()) {throw ELENotDefined; }

  putchar('\n');
  for(int i=0; i<vec.m_pnDimSize[0]; i++) {
    printf("%e  ", double(vec[i]));
  } // for i
  printf("\n");
} // PrintVector



