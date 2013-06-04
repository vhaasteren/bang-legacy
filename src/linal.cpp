/* linal.cpp -- linear algebra wrappers for libraries as gsl and lapack

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


// ******************************************************************************
// ***************** This todo list is /not/ up-to-date *************************
// ******************************************************************************


// Todo: add exception handling
//       add column handling for vectors
//           -tensor product
//
// Todo: Check why this gives segmentation fault:
//       oData.vdLambda[nCol] = vdLambdaAlpha[i]*vdLambdaTau[j];
//        +--> in CNumber::operator= the num.Defined() segfaults. Workaround: use double() operators
//
//       and
//       vdTempData = Sqrt(oData.vdLambda*2)&&InvErf(vdRand*2-1)
//
//
// Doc: Never allow for statements like: A*B;


#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_rng.h>

#if HAVE_LIBCBLAS_ATLAS == 1

extern "C" {
#if HAVE_LIBLAPACK_NETLIB == 1
  #include <f2c.h>
  #include <blaswrap.h>
#endif
  #include <cblas.h>
}

#elif HAVE_LIBCBLAS_MKL == 1

#include "mkl.h"

#elif HAVE_LIBBLAS_GSL == 1

#include <gsl/gsl_blas.h>

#elif HAVE_LIBCBLAS_GSL == 1

#include <gsl/gsl_cblas.h>

#elif HAVE_LIBCBLAS_VECLIB == 1

#include <vecLib/vecLib.h>

#elif HAVE_LIBCBLAS_CBLAS == 1

#include <cblas.h>

#endif

#if HAVE_LIBLAPACK_NETLIB == 1
// Include the lapack/atlas libraries for fast computation
extern "C" {
  #include <f2c.h>
  #include <clapack.h>
}

#elif HAVE_LIBLAPACK_CLAPACK == 1

#include <f2c.h>
#include <clapack.h>

#elif HAVE_LIBLAPACK_MKL == 1

#include "mkl.h"
#define integer int
#define real float
#define doublereal double


#elif HAVE_LIBLAPACK_GSL // Use the GSL routines

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#elif HAVE_LIBLAPACK_VECLIB == 1

#define integer int
#define real float
#define doublereal double

#elif HAVE_LIBLAPACK_LAPACK

#define integer int
#define real float
#define doublereal double

/* Declarations used in this file (don't use the headers) */
extern "C" integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1,  integer *n2, integer *n3, integer *n4);
extern "C" int dpotrf_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info);
extern "C" int dpotri_(char *uplo, integer *n, doublereal *a, integer *lda, integer *info);
extern "C" int dsyevd_(char *jobz, char *uplo, integer *n, doublereal *a, integer *lda, doublereal *w, doublereal *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
extern "C" int spotrf_(char *uplo, integer *n, real *a, integer *lda, integer *info);
extern "C" int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *info);
extern "C" int dgetrf_(integer *m, integer *n, doublereal *a, integer *lda, integer *ipiv, integer *info);
extern "C" int dgetri_(integer *n, doublereal *a, integer *lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);
extern "C" int dgeqp3_(integer *m, integer *n, doublereal *a, integer *lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork, integer *info);
extern "C" int dorgqr_(integer *m, integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

#endif



#include "linal.h"
#include "linalfunc.h"


CLinearObject::CLinearObject() {
  // Constructor of the linear object
  // No default data for this object
  InitializeNull();
}  // CLinearObject

CLinearObject::CLinearObject(int n) {
  InitializeNull();
} // CLinearObject

CLinearObject::CLinearObject(int n, int m) {
  InitializeNull();
} // CLinearObject

CLinearObject::CLinearObject(int nDim, int *pnDimSize) {
  InitializeNull();
} // CLinearObject


CLinearObject::~CLinearObject() {
  // Destructor of the linear object
  DeleteData();
} // ~CLinearObject


//#define __USE_GSL_RNG
void CLinearObject::Randomize() {
  int nElements=1;
#ifdef __USE_GSL_RNG
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
#endif

  if(! m_bDefined) {throw ELENotDefined; }

  // initialize random generator
  // Should be done only once (now done in LinalFuncInitialize() in linalfunc.cpp)
  // srand ( time(NULL) );

  // The original type of object determines what has to be done
  switch(m_eType) {
  case ENumber:
#ifdef __USE_GSL_RNG
    m_pdData[0] = gsl_rng_uniform(rng); //double(rand())/double(RAND_MAX);
#else
    m_pdData[0] = double(rand())/double(RAND_MAX);
#endif
    break;
  case EVector:
    for(int i=0; i<m_pnDimSize[0]; i++) {
#ifdef __USE_GSL_RNG
      m_pdData[i] = gsl_rng_uniform(rng); //double(rand())/double(RAND_MAX);
#else
      m_pdData[i] = double(rand())/double(RAND_MAX);
#endif
    } // for
    break;
  case EMatrix:
  case ETensor:
    for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
      nElements = nElements*m_pnDimSize[i];
    } // for
    m_pdData = new double[nElements];
    for(int i=0; i<nElements; i++) {
#ifdef __USE_GSL_RNG
      m_pdData[i] = gsl_rng_uniform(rng); //double(rand())/double(RAND_MAX);
#else
      m_pdData[i] = double(rand())/double(RAND_MAX);
#endif
    } // for
    break;
  case ELinearObject:
  default:
    // This should never happen
    printf("Serious error: CLinearObject::Randomize\n");
    throw ELEWrongClassType;
    break;
  }  // switch
#ifdef __USE_GSL_RNG
  gsl_rng_free(rng);
#endif
} // Randomize
#undef __USE_GSL_RNG

// Functions to initialize a linear object, needed when
// an object is made, but not initialized (like CMatrix mat;)
void CLinearObject::Initialize(int n) {
  int *pnDimSize;
  double *pdData;
  // If this object exists, delete all data
  if(m_bDefined) { DeleteData(); }

  pnDimSize = new int[1];
  pnDimSize[0] = n;
  pdData = new double[n];

  this->SetRawData(pdData, 1, pnDimSize, false, EVector);
} // Initialize

void CLinearObject::Initialize(int n, int m){
  int *pnDimSize;
  double *pdData;
  // If this object exists, delete all data
  if(m_bDefined) { DeleteData(); }

  pnDimSize = new int[2];
  pnDimSize[0] = n;
  pnDimSize[1] = m;
  
  pdData = new double[int(n*m)];

  this->SetRawData(pdData, 2, pnDimSize, false, EMatrix);
} // Initialize


// Unset all data in this object (total reset)
void CLinearObject::InitializeNull() {
  m_bDefined = false;
  DeleteData();
  m_eType = ELinearObject;
} // InitializeNull

double &CLinearObject::operator[](int i) {
  if(! m_bDefined) {throw ELENotDefined; }
  return m_pdData[i];
} // operator[]

&CLinearObject::operator double() {
  if(! m_bDefined) {throw ELENotDefined; }
  return m_pdData[0];
} // operator double

// Add a number (element-wise)
CLinearObject &CLinearObject::operator+=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] += d;
  } // for
  return *this;
} // operator+

// Subtract a number (element-wise)
CLinearObject &CLinearObject::operator-=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] -= d;
  } // for
  return *this;
} // operator-

// Multiply with a number (element-wise)
CLinearObject &CLinearObject::operator*=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] = m_pdData[i]*d;
  } // for
  return *this;
} // operator*

// Devide by a number (element-wise)
CLinearObject &CLinearObject::operator/=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] = m_pdData[i]/d;
  } // for
  return *this;
} // operator/

// Add a number (element-wise)
CLinearObject &CLinearObject::operator+(double d) {
  CLinearObject *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] += d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CLinearObject;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i] + d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator+

// Add a number (element-wise)
CLinearObject &CLinearObject::operator-(double d) {
  CLinearObject *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] -= d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CLinearObject;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i] - d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator-

// Add a number (element-wise)
CLinearObject &CLinearObject::operator*(double d) {
  CLinearObject *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] = m_pdData[i]*d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CLinearObject;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i]*d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator*

// Add a number (element-wise)
CLinearObject &CLinearObject::operator/(double d) {
  CLinearObject *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] = m_pdData[i]/d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CLinearObject;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i]/d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator/


void CLinearObject::SetRawData(double *pdData, int nDimensions, int *pnDimSize, bool bAllocated, ELinearClassType eType) {
  // Set the raw data of a Linear Object. No deletion of previous data for speed issues (can be changed)
  // DeleteData()
  m_nDimensions = nDimensions;
  m_pnDimSize = pnDimSize;
  m_pdData = pdData;
  m_bDefined = true;
  m_bAllocated = bAllocated;
  m_eType = eType;
} // SetRawData

void CLinearObject::DeleteData() {
//  int nIndexLength;

  // Delete all the data
  if(m_bDefined) {
    delete[] m_pnDimSize;
    delete[] m_pdData;
  } // if
  m_pdData = NULL;
  m_pnDimSize = NULL;
  m_nDimensions = 0;
  m_bDefined = false;
  m_bAllocated = false;
  m_dFlags = 0.0;
} // DeleteData

// This function makes a deep copy of another linear object. If this object was
// obtained by matrix algebra, then the original has to be destroyed. Then no
// copying is done, just movement of data. Further this function allows for flag
// notation and is therefore the heart of the entire Class Structure
//
// Todo: make flags possible for target, not only the source
//       after this function, the flags must always be reset
//       Chaining (a[flag] = b[flag] = c[flag] therefore is not possible
void CLinearObject::CopyFromLinearObject(CLinearObject &linob) {
  int nElements = 1;

  DeleteData();
  CopySingleMembers(linob);
  if(linob.m_bDefined) {
    if(linob.m_bAllocated) { // linob must be deleted
      m_pdData = linob.m_pdData;
      m_pnDimSize = linob.m_pnDimSize;
      m_bAllocated = false;
      linob.m_pdData = NULL;
      linob.m_pnDimSize = NULL;
      linob.m_bDefined = false;
      linob.m_bAllocated = false;
      linob.m_dFlags = 0.0;
      delete &linob;
    } else { // make a deep copy
      // The original type of object determines what has to be done
      switch(linob.m_eType) {
      case ENumber:
          m_pdData = new double[1];
          m_pnDimSize = new int[1];
          m_pdData[0] = linob.m_pdData[0];
          m_pnDimSize[0] = linob.m_pnDimSize[0];
        break;
      case EVector:
          m_pdData = new double[linob.m_pnDimSize[0]];
          m_pnDimSize = new int[1];
          m_pnDimSize[0] = linob.m_pnDimSize[0];
          for(int i=0; i<m_pnDimSize[0]; i++) {
            m_pdData[i] = linob.m_pdData[i];
          } // for
        break;
      case EMatrix:
      case ETensor:
          // First find out how many elements there are:
          m_pnDimSize = new int[linob.m_nDimensions];
          for(int i=0; i<linob.m_nDimensions; i++) {        // Go through each dimension
            nElements = nElements*linob.m_pnDimSize[i];
            m_pnDimSize[i] = linob.m_pnDimSize[i];
          } // for
          m_pdData = new double[nElements];
          for(int i=0; i<nElements; i++) {          // Copy all elements
            m_pdData[i] = linob.m_pdData[i];
          } // for
        break;
      case ELinearObject:
      default:
          // This should never happen
          printf("Serious error: CNumber::operator=\n");
          throw ELEWrongClassType;
        break;
      }  // switch
    }  // if

    // Reset the original linear object
    linob.ResetFlags();
  } else {
    m_pdData = NULL;
    m_pnDimSize = NULL;
  } // if
} // CopyFromLinearObject

int CLinearObject::GetElements() {
  int nElements=1;
  if(! m_bDefined) {throw ELENotDefined; }

  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  return nElements;
} // GetElements

// Todo: What to do with the flag parameter? Now only the source can
//       be flagged. Allow also for target.
void CLinearObject::CopySingleMembers(CLinearObject &linob) {
  m_bDefined = linob.m_bDefined;
  m_nDimensions = linob.m_nDimensions;
  m_eType = linob.m_eType;
  m_dFlags = 0.0;
} // CopySingleMembers

void CLinearObject::ResetFlags() {
  // Reset the flags for the object
  switch(int(m_dFlags)) {
  default:
    m_dFlags = 0;
    break;
  } // switch
} // ResetFlags

void CLinearObject::SetFlags(double d) {
  // Set the flags for the object
  ResetFlags();
  switch(int(d)) {
  default:
    m_dFlags = d;
    break;
  } // switch
} // SetFlags

double *CLinearObject::GetPD() {
  double *pd=NULL;

  if(! m_bDefined) {throw ELENotDefined; }

  // The original type of object determines what has to be done
  switch(m_eType) {
  case EMatrix:
    pd = &(m_pdData[m_indexj*m_pnDimSize[0]+m_indexi]);
    break;
  case EVector:
    pd = &(m_pdData[m_indexi]);
    break;
  case ENumber:
    pd = m_pdData;
    break;
  case ETensor:
  case ELinearObject:
  default:
    throw ELEWrongClassType;
    break;
  }  // switch
  return (pd);
} // CLinearObject::GetPD


CNumber::CNumber() {
  // Constructor
  InitializeNull();
  m_eType = ENumber;

  m_nDimensions = 0;
  m_pnDimSize = new int[1];
  m_pnDimSize[0] = 1;
  m_pdData = new double[1];
  m_bDefined = true;
  m_bAllocated = false;
} // CNumber


CNumber::CNumber(double d) {
  // Constructor
  m_eType = ENumber;

  m_nDimensions = 0;
  m_pnDimSize = new int[1];
  m_pnDimSize[0] = 1;
  m_pdData = new double[1];
  m_pdData[0] = d;
  m_bDefined = true;
  m_bAllocated = false;
} // CNumber

// Copy constructor for CNumber
CNumber::CNumber(CNumber &num) {
  CopyFromLinearObject(num);
} // CNumber


CNumber::~CNumber() {
  // Destructor
  DeleteData();
} // ~CNumber

// Here, the CNumber could also be a casted version of CVector or CMatrix
double &CNumber::operator=(double d) {
  // Assignment operator
  double *pd=NULL;

  if(! m_bDefined) {throw ELENotDefined; }

  // The original type of object determines what has to be done
  switch(m_eType) {
  case EMatrix:
      m_pdData[m_indexj*m_pnDimSize[0]+m_indexi] = d;
      pd = &(m_pdData[m_indexj*m_pnDimSize[0]+m_indexi]);
    break;
  case EVector:
      m_pdData[m_indexi] = d;
      pd = &(m_pdData[m_indexi]);
    break;
  case ENumber:
    m_pdData[0] = d;
    pd = m_pdData;
    break;
  case ETensor:
  case ELinearObject:
  default:
      // This should never happen
      printf("Serious error: CNumber::operator=\n");
      throw ELEWrongClassType;
    break;
  }  // switch
  return (*pd);
} // operator=

double &CNumber::operator=(int n) {
  // Assignment operator
  return this->operator=(double(n));
} // operator=

CNumber &CNumber::operator=(CNumber &num) {
  if(! num.Defined()) {throw ELENotDefined; }

  this->operator=(double(num));
  if(num.Allocated()) { delete &num; }
  return *this;
} // operator=

&CNumber::operator double() {
  return (*this->GetPD());
} // operator double

CNumber &CNumber::operator*(CNumber &num) {
  CNumber *newnum;
  double *pdNewData;
  if(! (m_bDefined && num.Defined())) {throw ELENotDefined; }


  switch(m_bAllocated) {
  case true:
    newnum = (CNumber *)this;
    pdNewData = new double[1];
    pdNewData[0] = double(*this)*double(num);
    delete[] m_pdData;
    newnum->SetRawData(pdNewData, 0, m_pnDimSize, true, ENumber);
    break;
  case false:
  default:
    CNumber *newnum=new CNumber(double(*this)*double(num));
    newnum->m_bAllocated = true;
    break;
  } // switch

  if(num.Allocated()) {delete &num; } // This should call for the destructor

  return *newnum;
} // operator *

CVector &CNumber::operator*(CVector &vec) {
  CVector *newvec;
  double *pdNewData;
  int *pnDimSize;

  if(! (m_bDefined && vec.Defined())) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // Delete the target data to allow M*N*O
    pdNewData = new double[int(vec.m_pnDimSize[0])];

    for(int i=0; i<vec.m_pnDimSize[0]; i++) {
      pdNewData[i] = m_pdData[0]*double(vec.m_pdData[i]);
    } // for

    newvec = (CVector *)this;

    pnDimSize = m_pnDimSize;
    pnDimSize[0] = vec.m_pnDimSize[0];
    
    delete[] m_pdData;
    newvec->SetRawData(pdNewData, 1, pnDimSize, true, EVector);

    break;
  case false:
    newvec = new CVector;
    pdNewData = new double[int(vec.m_pnDimSize[0])];

    for(int i=0; i<vec.m_pnDimSize[0]; i++) {
      pdNewData[i] = m_pdData[0]*double(vec.m_pdData[i]);
    } // for

    pnDimSize = new int[int(vec.m_pnDimSize[0])];
    pnDimSize[0] = vec.m_pnDimSize[0];
    newvec->SetRawData(pdNewData, 1, pnDimSize, true, EVector);

    break;
  default:
    break;
  } // switch

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  
  return *newvec;
} // operator *

CMatrix &CNumber::operator*(CMatrix &mat) {
  CMatrix *newmat;
  double *pdNewData;
  int *pnDimSize;

  if(! (m_bDefined && mat.Defined())) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // Delete the target data to allow M*N*O
    pdNewData = new double[int(mat.m_pnDimSize[0]*mat.m_pnDimSize[1])];

    for(int i=0; i<int(mat.m_pnDimSize[0]*mat.m_pnDimSize[1]); i++) {
      pdNewData[i] = m_pdData[0]*mat.m_pdData[i];
    } // for

    newmat = (CMatrix *)this;

    pnDimSize = new int[2];
    pnDimSize[0] = mat.m_pnDimSize[0];
    pnDimSize[1] = mat.m_pnDimSize[1];
    
    delete[] m_pnDimSize;
    delete[] m_pdData;
    newmat->SetRawData(pdNewData, 2, pnDimSize, true, EMatrix);

    break;
  case false:
    newmat = new CMatrix;
    pdNewData = new double[int(mat.m_pnDimSize[0]*mat.m_pnDimSize[1])];

    for(int i=0; i<int(mat.m_pnDimSize[0]*mat.m_pnDimSize[1]); i++) {
      pdNewData[i] = m_pdData[0]*mat.m_pdData[i];
    } // for

    pnDimSize = new int[2];
    pnDimSize[0] = mat.m_pnDimSize[0];
    pnDimSize[1] = mat.m_pnDimSize[1];

    newmat->SetRawData(pdNewData, 2, pnDimSize, true, EMatrix);
   break;
  default:
    break;
  } // switch

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  
  return *newmat;
} // operator *

// Add a number (element-wise)
CNumber &CNumber::operator+=(double d) {
//  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
/*  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] += d;
  } // for*/
  (*this->GetPD()) += d;
  return *this;
} // operator+

// Subtract a number (element-wise)
CNumber &CNumber::operator-=(double d) {
//  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
/*  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] -= d;
  } // for*/
  (*this->GetPD()) -= d;

  return *this;
} // operator-

// Multiply with a number (element-wise)
CNumber &CNumber::operator*=(double d) {
//  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
/*  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] = m_pdData[i]*d;
  } // for*/
  (*this->GetPD()) *= d;
  return *this;
} // operator*

// Devide by a number (element-wise)
CNumber &CNumber::operator/=(double d) {
//  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
/*  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] = m_pdData[i]/d;
  } // for*/
  (*this->GetPD()) /= d;
  return *this;
} // operator/

// Add a number (element-wise)
CNumber &CNumber::operator+(double d) {
  CNumber *newob;
  double *pdNewData;
  int *pnDimSize;

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
   (*this->GetPD()) += d;
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CNumber;
    pdNewData = new double[1];

    pdNewData[0] = (*this->GetPD()) + d;

    pnDimSize = new int[1];
    pnDimSize[0] = 1;

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, ENumber);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator+

// Add a number (element-wise)
CNumber &CNumber::operator-(double d) {
  CNumber *newob;
  double *pdNewData;
  int *pnDimSize;

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
   (*this->GetPD()) -= d;
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CNumber;
    pdNewData = new double[1];

    pdNewData[0] = (*this->GetPD()) - d;

    pnDimSize = new int[1];
    pnDimSize[0] = 1;

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, ENumber);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator-

// Add a number (element-wise)
CNumber &CNumber::operator*(double d) {
  CNumber *newob;
  double *pdNewData;
  int *pnDimSize;

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
   (*this->GetPD()) *= d;
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CNumber;
    pdNewData = new double[1];

    pdNewData[0] = (*this->GetPD()) * d;

    pnDimSize = new int[1];
    pnDimSize[0] = 1;

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, ENumber);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator*

// Add a number (element-wise)
CNumber &CNumber::operator/(double d) {
  CNumber *newob;
  double *pdNewData;
  int *pnDimSize;

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
   (*this->GetPD()) /= d;
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CNumber;
    pdNewData = new double[1];

    pdNewData[0] = (*this->GetPD()) / d;

    pnDimSize = new int[1];
    pnDimSize[0] = 1;

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, ENumber);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator/


CVector::CVector() {
  // Constructor
  InitializeNull();
  m_eType = EVector;
} // CVector

CVector::CVector(int n) {
  // Constructor
  m_eType = EVector;

  m_nDimensions = 1;
  m_pnDimSize = new int[1];
  m_pnDimSize[0] = n;
  m_pdData = new double[n];
  m_bDefined = true;
  m_bAllocated = false;
} // CVector

CVector::CVector(CVector &vec) {
  // Copy constructor
  CopyFromLinearObject(vec);
} // CVector

CVector::~CVector() {
  // Destructor
  DeleteData();
} // ~CVector

CNumber &CVector::operator[](int i) {
  // Index operator

  if(! m_bDefined) {throw ELENotDefined; }

  // The original type of object determines what has to be done
  switch(m_eType) {
  case EMatrix:
      if(i >= m_pnDimSize[1]) {throw ELEBadIndex;}
      m_indexj = i;
    break;
  case EVector:
      if(i >= m_pnDimSize[0]) {throw ELEBadIndex;}
    m_indexi = i;
    break;
  case ETensor:
  case ELinearObject:
  case ENumber:
  default:
      throw ELEWrongClassType;
    break;
  }  // switch
  return *((CNumber *)this);
}

CVector &CVector::operator=(CVector &vec) {
  CopyFromLinearObject(vec);
  return *this;
} // operator=


/*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/


// Todo: Use function SetRawData();
CNumber &CVector::operator*(CVector &vec) { // * operator
  // Use the lapack/atlas libraries for fast matrix multiplication
  CNumber *newnum;
  double *pdNewData;
  int *pnDimSize;

#ifdef HAVE_LIBCBLAS
/*  integer inCols;
  doublereal doAlpha, doBeta;

  inCols = (integer)1;
  doAlpha = (doublereal)1.0;
  doBeta = (doublereal)0.0;*/
#else // HAVE_LIBCBLAS
  gsl_matrix_view mat1, mat2, mat3;
#endif  // HAVE_LIBCBLAS

  if(! (m_bDefined && vec.m_bDefined)) {throw ELENotDefined; }
  if(m_pnDimSize[0] != vec.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // Delete the target data to allow M*N*O
    pdNewData = new double[1];


#ifdef HAVE_LIBCBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,1,1,m_pnDimSize[0],1.0,m_pdData,1,vec.m_pdData,vec.m_pnDimSize[0],0.0,pdNewData,vec.m_pnDimSize[0]);
//    dgemm_("n","n",&inCols,&inCols,(integer *)&m_pnDimSize[0],&doAlpha,(doublereal *)m_pdData,(integer *)&inCols,(doublereal *)vec.m_pdData,&inCols,&doBeta,(doublereal *)pdNewData,&inCols);


#else // HAVE_LIBCBLAS
    // USE THE GNU SCIENTIFIC LIBRARY
    mat1 = gsl_matrix_view_array(m_pdData, 1, m_pnDimSize[0]);
    mat2 = gsl_matrix_view_array(vec.m_pdData, m_pnDimSize[0], 1);
    mat3 = gsl_matrix_view_array(pdNewData, 1, 1);

    // C = \alpha*A*B + \beta*C
    gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&mat1.matrix,&mat2.matrix,0,&mat3.matrix);
    // use the gnu scientific library
#endif // HAVE_LIBCBLAS


    delete[] m_pdData;
    newnum = (CNumber *)this;

    pnDimSize = m_pnDimSize;
    pnDimSize[0] = 1;
    newnum->SetRawData(pdNewData, 0, pnDimSize, true, ENumber);

    break;
  case false:
    newnum = new CNumber;
    pdNewData = new double[1];

#ifdef HAVE_LIBCBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,1,1,m_pnDimSize[0],1.0,m_pdData,1,vec.m_pdData,vec.m_pnDimSize[0],0.0,pdNewData,vec.m_pnDimSize[0]);

//    dgemm_("n","n",&inCols,&inCols,(integer *)&m_pnDimSize[0],&doAlpha,(doublereal *)m_pdData,(integer *)&inCols,(doublereal *)vec.m_pdData,&inCols,&doBeta,(doublereal *)pdNewData,&inCols);

#else // HAVE_LIBCBLAS
    // USE THE GNU SCIENTIFIC LIBRARY
    mat1 = gsl_matrix_view_array(m_pdData, 1, m_pnDimSize[0]);
    mat2 = gsl_matrix_view_array(vec.m_pdData, m_pnDimSize[0], 1);
    mat3 = gsl_matrix_view_array(pdNewData, 1, 1);

    // C = \alpha*A*B + \beta*C
    gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&mat1.matrix,&mat2.matrix,0,&mat3.matrix);
    // use the gnu scientific library

#endif // HAVE_LIBCBLAS

    pnDimSize = new int[1];
    pnDimSize[0] = 1;
    newnum->SetRawData(pdNewData, 0, pnDimSize, true, ENumber);

    break;
  default:
    break;
  } // switch

  if(vec.m_bAllocated) {delete &vec; } // This should call for the destructor
  
  return *newnum;
} // operator*

CVector &CVector::operator*(CMatrix &mat) { // * operator
  // Use the lapack/atlas libraries for fast matrix multiplication
  CVector *newvec;
  double *pdNewData;
  int *pnDimSize;

#ifdef HAVE_LIBCBLAS
/*  integer inCols;
  doublereal doAlpha, doBeta;

  inCols = (integer)1;
  doAlpha = (doublereal)1.0;
  doBeta = (doublereal)0.0;*/
#else // HAVE_LIBCBLAS
  gsl_matrix_view mat1, mat2, mat3;
#endif // HAVE_LIBCBLAS

  if(! (m_bDefined && mat.Defined())) {throw ELENotDefined; }
  if(m_pnDimSize[0] != mat.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // Delete the target data to allow M*N*O
    pdNewData = new double[int(mat.m_pnDimSize[1])];

#ifdef HAVE_LIBCBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,1,mat.m_pnDimSize[1],m_pnDimSize[0],1.0,m_pdData,1,mat.m_pdData,mat.m_pnDimSize[0],0.0,pdNewData,1);
//    dgemm_("n","n",&inCols,(integer *)&mat.m_pnDimSize[1],(integer *)&m_pnDimSize[0],&doAlpha,(doublereal *)m_pdData,(integer *)&m_pnDimSize[0],(doublereal *)mat.m_pdData,(integer *)&mat.m_pnDimSize[1],&doBeta,(doublereal *)pdNewData,(integer *)&mat.m_pnDimSize[1]);

#else // HAVE_LIBCBLAS
    // USE THE GNU SCIENTIFIC LIBRARY
    mat1 = gsl_matrix_view_array(m_pdData, 1, m_pnDimSize[0]);
    mat2 = gsl_matrix_view_array(mat.m_pdData, mat.m_pnDimSize[1], mat.m_pnDimSize[0]);
    mat3 = gsl_matrix_view_array(pdNewData, 1, mat.m_pnDimSize[1]);

    // C = \alpha*A*B + \beta*C
    gsl_blas_dgemm (CblasNoTrans,CblasTrans,1,&mat1.matrix,&mat2.matrix,0,&mat3.matrix);
    // use the gnu scientific library

#endif // HAVE_LIBCBLAS

    newvec = this;

    delete[] m_pdData;

    newvec->SetRawData(pdNewData, 1, m_pnDimSize, true, EVector);
    break;
  case false:
    newvec = new CVector;
    pdNewData = new double[int(mat.m_pnDimSize[1])];

#ifdef HAVE_LIBCBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,1,mat.m_pnDimSize[1],m_pnDimSize[0],1.0,m_pdData,1,mat.m_pdData,mat.m_pnDimSize[0],0.0,pdNewData,1);
//    dgemm_("n","n",&inCols,(integer *)&mat.m_pnDimSize[1],(integer *)&m_pnDimSize[0],&doAlpha,(doublereal *)m_pdData,(integer *)&m_pnDimSize[0],(doublereal *)mat.m_pdData,(integer *)&mat.m_pnDimSize[1],&doBeta,(doublereal *)pdNewData,(integer *)&mat.m_pnDimSize[1]);

#else // HAVE_LIBCBLAS

    // USE THE GNU SCIENTIFIC LIBRARY
    mat1 = gsl_matrix_view_array(m_pdData, 1, m_pnDimSize[0]);
    mat2 = gsl_matrix_view_array(mat.m_pdData, mat.m_pnDimSize[1], mat.m_pnDimSize[0]);
    mat3 = gsl_matrix_view_array(pdNewData, 1, mat.m_pnDimSize[1]);

    // C = \alpha*A*B + \beta*C
    gsl_blas_dgemm (CblasNoTrans,CblasTrans,1,&mat1.matrix,&mat2.matrix,0,&mat3.matrix);
    // use the gnu scientific library
#endif // HAVE_LIBCBLAS


    pnDimSize = new int[1];
    pnDimSize[0] = m_pnDimSize[0];
    newvec->SetRawData(pdNewData, 1, pnDimSize, true, EVector);

    break;
  default:
    break;
  } // switch

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  
  return *newvec;
} // operator*

// Element-wise multiplication
CVector &CVector::operator&&(CVector &vec) {
  CVector *newvec;
  double *pdNewData;
  int *pnDimSize;

  if(! (m_bDefined && vec.Defined())) {throw ELENotDefined; }
  if(! (m_pnDimSize[0] == vec.m_pnDimSize[0])) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // Delete the target data to allow M*N*O
    pdNewData = new double[int(m_pnDimSize[0])];

// Maybe this can go through some kind of lapack routine?
    for(int i=0; i<m_pnDimSize[0]; i++) {
        pdNewData[i] = m_pdData[i]*vec.m_pdData[i];
    }

    newvec = this;
    delete[] m_pdData;

    newvec->SetRawData(pdNewData, 1, m_pnDimSize, true, EVector);

    break;
  case false:
    newvec = new CVector;
    pdNewData = new double[int(m_pnDimSize[0])];

    for(int i=0; i<m_pnDimSize[0]; i++) {
      pdNewData[i] = m_pdData[i]*vec.m_pdData[i];
    }

    pnDimSize = new int[1];
    pnDimSize[0] = m_pnDimSize[0];
    newvec->SetRawData(pdNewData, 1, pnDimSize, true, EMatrix);

    break;
  default:
    break;
  } // switch

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  
  return *newvec;
} // operator.*


// Add a number (element-wise)
CVector &CVector::operator+=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] += d;
  } // for
  return *this;
} // operator+

// Subtract a number (element-wise)
CVector &CVector::operator-=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] -= d;
  } // for
  return *this;
} // operator-

// Multiply with a number (element-wise)
CVector &CVector::operator*=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] = m_pdData[i]*d;
  } // for
  return *this;
} // operator*

// Devide by a number (element-wise)
CVector &CVector::operator/=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] = m_pdData[i]/d;
  } // for
  return *this;
} // operator/


// Add a number (element-wise)
CVector &CVector::operator+(double d) {
  CVector *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] += d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CVector;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i] + d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator+

// Add a number (element-wise)
CVector &CVector::operator-(double d) {
  CVector *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] -= d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CVector;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i] - d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator-

// Add a number (element-wise)
CVector &CVector::operator*(double d) {
  CVector *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] = m_pdData[i]*d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CVector;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i]*d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator*

// Add a number (element-wise)
CVector &CVector::operator/(double d) {
  CVector *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] = m_pdData[i]/d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CVector;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i]/d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator/


// Add a number (element-wise)
CVector &CVector::operator+(CVector &vec) {
  CVector *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! (m_bDefined && vec.Defined())) {throw ELENotDefined; }
  if(m_pnDimSize[0] != vec.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] += double(vec.m_pdData[i]);
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CVector;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = double(m_pdData[i]) + double(vec.m_pdData[i]);
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  if(vec.m_bAllocated) {delete &vec; } // This should call for the destructor
  return *newob;
} // operator+

// Subtract a number (element-wise)
CVector &CVector::operator-(CVector &vec) {
  CVector *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! (m_bDefined && vec.Defined())) {throw ELENotDefined; }
  if(m_pnDimSize[0] != vec.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] -= vec.m_pdData[i];
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CVector;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i] - vec.m_pdData[i];
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  if(vec.m_bAllocated) {delete &vec; } // This should call for the destructor
  return *newob;
} // operator-


CMatrix::CMatrix() {
  // Constructor
  InitializeNull();
  m_eType = EMatrix;
} // CMatrix

CMatrix::CMatrix(int n, int m) {
  // Constructor
  m_eType = EMatrix;

  m_nDimensions = 2;
  m_pnDimSize = new int[2];
  m_pnDimSize[0] = n;
  m_pnDimSize[1] = m;
  m_pdData = new double[n*m];
  m_bDefined = true;
  m_bAllocated = false;
} // CMatrix

CMatrix::CMatrix(CMatrix &mat) {
  CopyFromLinearObject(mat);
} // CMatrix

CMatrix::~CMatrix() {
  // Destructor
  DeleteData();
} // ~CMatrix

CVector &CMatrix::operator[](int i) {
  // Index operator

  if(! m_bDefined) {throw ELENotDefined; }

  // The original type of object determines what has to be done
  switch(m_eType) {
  case EMatrix:
      if(i >= m_pnDimSize[0]) {throw ELEBadIndex;}
      m_indexi = i;
    break;
  case ETensor:
  case ELinearObject:
  case ENumber:
  case EVector:
  default:
      throw ELEWrongClassType;
    break;
  }  // switch
  return *((CVector *)this);
}  // operator[]

// Todo: this must be adjusted for generality. Now just quick 'n dirty
CMatrix &CMatrix::operator[](double d) {
  // Special functions index operator
  CMatrix *pmReturn = this, *pmTemp;
  double *pdNewData;
  int *pnDimSize;

  if(! m_bDefined) {throw ELENotDefined; }

  // The passed double determines the function that must be used
  SetFlags(d);  // Is nu de enige vereiste
  switch(int(d)) {
  case int(LO_ALL):
    // All elements must be assigned
    break;
  case int(LO_TRANSPOSE):
    // Return the transpose of this matrix
    pmTemp = new CMatrix;

    pdNewData = new double[int(m_pnDimSize[0]*m_pnDimSize[1])];

    // Transpose the data
    for(int i=0; i<m_pnDimSize[0]; i++) {
      for(int j=0; j<m_pnDimSize[1]; j++) {
        pdNewData[j+i*m_pnDimSize[1]] = m_pdData[i+j*m_pnDimSize[0]];
      } // for j
    } // for i
    
    pnDimSize = new int[2];
    pnDimSize[0] = m_pnDimSize[1];
    pnDimSize[1] = m_pnDimSize[0];
    
    pmTemp->SetRawData(pdNewData, 1, pnDimSize, true, EMatrix);
    pmReturn = pmTemp;
    break;
  default:
    break;
  } // switch
  return *pmReturn;
} // operator[]

CMatrix &CMatrix::operator=(CMatrix &mat) {
  CopyFromLinearObject(mat);
  return *this;
} // operator=

// Multiply a matrix with a matrix. Memory allocation problem solved
CMatrix &CMatrix::operator*(CMatrix &mat) {
  // Use the lapack/atlas libraries for fast matrix multiplication
  CMatrix *newmat;
  double *pdNewData;
  int *pnDimSize;

#ifdef HAVE_LIBCBLAS
/*  integer inCols;
  doublereal doAlpha, doBeta;

  inCols = (integer)1;
  doAlpha = (doublereal)1.0;
  doBeta = (doublereal)0.0;*/
#else // HAVE_LIBCBLAS
  gsl_matrix_view mat1, mat2, mat3;
#endif // HAVE_LIBCBLAS

  if(! (m_bDefined && mat.m_bDefined)) {throw ELENotDefined; }
  if(m_pnDimSize[1] != mat.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // Delete the target data to allow M*N*O
    pdNewData = new double[int(m_pnDimSize[0]*mat.m_pnDimSize[1])];

#ifdef HAVE_LIBCBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m_pnDimSize[0],mat.m_pnDimSize[1],m_pnDimSize[1],1.0,m_pdData,m_pnDimSize[0],mat.m_pdData,mat.m_pnDimSize[0],0.0,pdNewData,m_pnDimSize[0]);
//    dgemm_("n","n",(integer *)&m_pnDimSize[0],(integer *)&mat.m_pnDimSize[1],(integer *)&m_pnDimSize[1],&doAlpha,(doublereal *)m_pdData,(integer *)&m_pnDimSize[1],(doublereal *)mat.m_pdData,(integer *)&mat.m_pnDimSize[1],&doBeta,(doublereal *)pdNewData,(integer *)&mat.m_pnDimSize[1]);

#else // HAVE_LIBCBLAS
    // USE THE GNU SCIENTIFIC LIBRARY
    mat2 = gsl_matrix_view_array(m_pdData, m_pnDimSize[1], m_pnDimSize[0]);
    mat1 = gsl_matrix_view_array(mat.m_pdData, mat.m_pnDimSize[1], mat.m_pnDimSize[0]);
    mat3 = gsl_matrix_view_array(pdNewData, mat.m_pnDimSize[1], m_pnDimSize[0]);

    // C = \alpha*A*B + \beta*C
    gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&mat1.matrix,&mat2.matrix,0,&mat3.matrix);
    // use the gnu scientific library
#endif // HAVE_LIBCBLAS


    newmat = this;
    delete[] m_pdData;

    pnDimSize = m_pnDimSize;
    pnDimSize[1] = mat.m_pnDimSize[1];
    newmat->SetRawData(pdNewData, 2, m_pnDimSize, true, EMatrix);

    break;
  case false:
    newmat = new CMatrix;
    pdNewData = new double[int(m_pnDimSize[0]*mat.m_pnDimSize[1])];

#ifdef HAVE_LIBCBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m_pnDimSize[0],mat.m_pnDimSize[1],m_pnDimSize[1],1.0,m_pdData,m_pnDimSize[0],mat.m_pdData,mat.m_pnDimSize[0],0.0,pdNewData,m_pnDimSize[0]);
//    dgemm_("n","n",(integer *)&m_pnDimSize[0],(integer *)&mat.m_pnDimSize[1],(integer *)&m_pnDimSize[1],&doAlpha,(doublereal *)m_pdData,(integer *)&m_pnDimSize[1],(doublereal *)mat.m_pdData,(integer *)&mat.m_pnDimSize[1],&doBeta,(doublereal *)pdNewData,(integer *)&mat.m_pnDimSize[1]);

#else // HAVE_LIBCBLAS
    // USE THE GNU SCIENTIFIC LIBRARY
    mat2 = gsl_matrix_view_array(m_pdData, m_pnDimSize[1], m_pnDimSize[0]);
    mat1 = gsl_matrix_view_array(mat.m_pdData, mat.m_pnDimSize[1], mat.m_pnDimSize[0]);
    mat3 = gsl_matrix_view_array(pdNewData, mat.m_pnDimSize[1], m_pnDimSize[0]);

    // C = \alpha*A*B + \beta*C
    gsl_blas_dgemm (CblasNoTrans,CblasNoTrans,1,&mat1.matrix,&mat2.matrix,0,&mat3.matrix);
    // use the gnu scientific library
#endif // HAVE_LIBCBLAS

    pnDimSize = new int[2];
    pnDimSize[0] = m_pnDimSize[0];
    pnDimSize[1] = mat.m_pnDimSize[1];
    newmat->SetRawData(pdNewData, 2, pnDimSize, true, EMatrix);

    break;
  default:
    break;
  } // switch

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  
  return *newmat;
} // operator*

// Multiply with the inverse of a matrix (division)
CMatrix &CMatrix::operator/(CMatrix &mat) {
  CMatrix matInverse;
  
  // These restrictions must be more strict
  if(! (m_bDefined && mat.m_bDefined)) {throw ELENotDefined; }
  if(m_pnDimSize[1] != mat.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  matInverse = mat.Inverse();
  return ((*this) * matInverse);
} // operator /

// Multiply a matrix with a vector
CVector &CMatrix::operator*(CVector &vec) {
  // Use the lapack/atlas libraries for fast matrix multiplication
  CVector *newvec;
  double *pdNewData;
  int *pnDimSize;

#ifdef HAVE_LIBCBLAS
/*  integer inCols;
  doublereal doAlpha, doBeta;

  inCols = (integer)1;
  doAlpha = (doublereal)1.0;
  doBeta = (doublereal)0.0;*/
#else // HAVE_LIBCBLAS
  gsl_matrix_view mat1, mat2, mat3;
#endif // HAVE_LIBCBLAS

  if(! (m_bDefined && vec.Defined())) {throw ELENotDefined; }
  if(m_pnDimSize[1] != vec.m_pnDimSize[0]) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // Delete the target data to allow M*N*O
    pdNewData = new double[int(m_pnDimSize[0])];

#ifdef HAVE_LIBCBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m_pnDimSize[0],1,m_pnDimSize[1],1.0,m_pdData,m_pnDimSize[0],vec.m_pdData,vec.m_pnDimSize[0],0.0,pdNewData,m_pnDimSize[0]);
//    dgemm_("n","n",(integer *)&m_pnDimSize[0],&inCols,(integer *)&m_pnDimSize[1],&doAlpha,(doublereal *)m_pdData,(integer *)&m_pnDimSize[1],(doublereal *)vec.m_pdData,&inCols,&doBeta,(doublereal *)pdNewData,&inCols);

#else // HAVE_LIBCBLAS
    // USE THE GNU SCIENTIFIC LIBRARY
    mat1 = gsl_matrix_view_array(m_pdData, m_pnDimSize[1], m_pnDimSize[0]);
    mat2 = gsl_matrix_view_array(vec.m_pdData, vec.m_pnDimSize[0], 1);
    mat3 = gsl_matrix_view_array(pdNewData, m_pnDimSize[0], 1);

    // C = \alpha*A*B + \beta*C
    gsl_blas_dgemm (CblasTrans,CblasNoTrans,1,&mat1.matrix,&mat2.matrix,0,&mat3.matrix);
    // use the gnu scientific library
#endif // HAVE_LIBCBLAS


    newvec = (CVector *)this;
    
    pnDimSize = new int[1];
    pnDimSize[0] = m_pnDimSize[0];

    delete[] m_pnDimSize;
    delete[] m_pdData;
    
    newvec->SetRawData(pdNewData, 1, pnDimSize, true, EVector);

    break;
  case false:
    newvec = new CVector;
    pdNewData = new double[int(m_pnDimSize[0])];

#ifdef HAVE_LIBCBLAS
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m_pnDimSize[0],1,m_pnDimSize[1],1.0,m_pdData,m_pnDimSize[0],vec.m_pdData,vec.m_pnDimSize[0],0.0,pdNewData,m_pnDimSize[0]);
//    dgemm_("n","n",(integer *)&m_pnDimSize[0],&inCols,(integer *)&m_pnDimSize[1],&doAlpha,(doublereal *)m_pdData,(integer *)&m_pnDimSize[1],(doublereal *)vec.m_pdData,&inCols,&doBeta,(doublereal *)pdNewData,&inCols);

#else // HAVE_LIBCBLAS
    // USE THE GNU SCIENTIFIC LIBRARY
    mat1 = gsl_matrix_view_array(m_pdData, m_pnDimSize[1], m_pnDimSize[0]);
    mat2 = gsl_matrix_view_array(vec.m_pdData, vec.m_pnDimSize[0], 1);
    mat3 = gsl_matrix_view_array(pdNewData, m_pnDimSize[0], 1);

    // C = \alpha*A*B + \beta*C
    gsl_blas_dgemm (CblasTrans,CblasNoTrans,1,&mat1.matrix,&mat2.matrix,0,&mat3.matrix);
    // use the gnu scientific library
#endif // HAVE_LIBCBLAS


    pnDimSize = new int[1];
    pnDimSize[0] = m_pnDimSize[0];

    newvec->SetRawData(pdNewData, 1, pnDimSize, true, EVector);

    break;
  default:
    break;
  } // switch

  if(vec.Allocated()) {delete &vec; } // This should call for the destructor
  
  return *newvec;
} // operator*

// Multiply a matrix with another matrix (element-wise)
CMatrix &CMatrix::operator&&(CMatrix &mat) {
  CMatrix *newmat;
  double *pdNewData;
  int *pnDimSize;

  if(! (m_bDefined && mat.Defined())) {throw ELENotDefined; }
  if(! ((m_pnDimSize[0] == mat.m_pnDimSize[0])&&(m_pnDimSize[1] == mat.m_pnDimSize[1]))) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // Delete the target data to allow M*N*O
    pdNewData = new double[int(m_pnDimSize[0]*m_pnDimSize[1])];

// Maybe this can go through some kind of lapack routine?
    for(int i=0; i<m_pnDimSize[0]; i++) {
      for(int j=0; j<m_pnDimSize[1]; j++) {
        pdNewData[j+m_pnDimSize[0]*i] = m_pdData[j+m_pnDimSize[0]*i]*mat.m_pdData[j+m_pnDimSize[0]*i];
      }
    }

    newmat = this;
    delete[] m_pdData;

    newmat->SetRawData(pdNewData, 2, m_pnDimSize, true, EMatrix);

    break;
  case false:
    newmat = new CMatrix;
    pdNewData = new double[int(m_pnDimSize[0]*m_pnDimSize[1])];

    for(int i=0; i<m_pnDimSize[0]; i++) {
      for(int j=0; j<m_pnDimSize[1]; j++) {
        pdNewData[j+m_pnDimSize[0]*i] = m_pdData[j+m_pnDimSize[0]*i]*mat.m_pdData[j+m_pnDimSize[0]*i];
      }
    }

    pnDimSize = new int[2];
    pnDimSize[0] = m_pnDimSize[0];
    pnDimSize[1] = m_pnDimSize[1];
    newmat->SetRawData(pdNewData, 2, pnDimSize, true, EMatrix);

    break;
  default:
    break;
  } // switch

  if(mat.Allocated()) {delete &mat; } // This should call for the destructor
  
  return *newmat;
} // operator .*

// & operator (tensor product)
CMatrix &CMatrix::operator&(CMatrix &mat) {
  // Use the lapack/atlas libraries for fast matrix multiplication
  CMatrix *newmat;
  double *pdNewData, *pdMatData;
  int *pnDimSize, k, l, i, ii, j, jj, nRow, nCol;

  if(! (m_bDefined && mat.m_bDefined)) {throw ELENotDefined; }
  if(m_pnDimSize[0] != m_pnDimSize[1] || mat.m_pnDimSize[0] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }
  k = m_pnDimSize[0];
  l = mat.m_pnDimSize[0];
  pdMatData = mat.GetRawData();

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // Delete the target data to allow M&N&O
    pdNewData = new double[int(m_pnDimSize[0]*m_pnDimSize[1]*mat.m_pnDimSize[0]*mat.m_pnDimSize[1])];

    // There has to be a lapack/atlas routine that does the same. Which one?
    // cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m_pnDimSize[0],mat.m_pnDimSize[1],m_pnDimSize[1],1.0,m_pdData,m_pnDimSize[1],mat.m_pdData,mat.m_pnDimSize[1],0.0,pdNewData,mat.m_pnDimSize[1]);
    // Hier komt ie..

    // Base vectors are ordered first by i, then by j (steps of k)
    for(i=0; i<k; i++) {
        for(j=0; j<l; j++) {
            // Now construct the ij'th part of the tensor product
            nCol = i + j*k;   // number of column we are computing
            for(ii=0; ii<k; ii++) {
                for(jj=0; jj<l; jj++) {
                    nRow = ii + jj*k;
                    pdNewData[nCol+nRow*k*l] = m_pdData[i+k*ii]*pdMatData[j+l*jj];
                }
            } // go to the next element of the column
        }
    }

    newmat = this;
    delete[] m_pdData;

    pnDimSize = m_pnDimSize;
    pnDimSize[0] = pnDimSize[1] = k*l;
    newmat->SetRawData(pdNewData, 2, m_pnDimSize, true, EMatrix);

    break;
  case false:
    newmat = new CMatrix;
    pdNewData = new double[int(m_pnDimSize[0]*m_pnDimSize[1]*mat.m_pnDimSize[0]*mat.m_pnDimSize[1])];

//    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,m_pnDimSize[0],mat.m_pnDimSize[1],m_pnDimSize[1],1.0,m_pdData,m_pnDimSize[1],mat.m_pdData,mat.m_pnDimSize[1],0.0,pdNewData,mat.m_pnDimSize[1]);

    // Base vectors are ordered first by i, then by j (steps of k)
    for(i=0; i<k; i++) {
        for(j=0; j<l; j++) {
            // Now construct the ij'th part of the tensor product
            nCol = i + j*k;   // number of column we are computing
            for(ii=0; ii<k; ii++) {
                for(jj=0; jj<l; jj++) {
                    nRow = ii + jj*k;
                    pdNewData[nCol+nRow*k*l] = m_pdData[i+k*ii]*pdMatData[j+l*jj];
                }
            } // go to the next element of the column
        }
    }

    pnDimSize = new int[2];
    pnDimSize[0] = pnDimSize[1] = k*l;
    newmat->SetRawData(pdNewData, 2, pnDimSize, true, EMatrix);

    break;
  default:
    break;
  } // switch

  if(mat.m_bAllocated) {delete &mat; } // This should call for the destructor
  
  return *newmat;
} // operator&

/* Calculates the inverse of a matrix using LU-factorization
   The original matrix is overwritten by this function */
CMatrix &CMatrix::Invert(double *pdLogDet) {
  if(! m_bDefined) {throw ELENotDefined; }
  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

#ifdef HAVE_LIBLAPACK

  // Invert this matrix using the clapack/atlas routine
  integer info, *ipiv, ispec=1, N1, N2, N3, N4, k;
  doublereal *pdWorkSpace;
#if HAVE_LIBLAPACK_MKL
  const char *name="dgetri", *opts="";
#else
  char *name="dgetri", *opts="";
#endif
  long namelen = 6, optslen=0;
  
  N1 = N2 = integer(m_pnDimSize[0]);
  N3 = N4 = integer(-1);

  ipiv = new integer[m_pnDimSize[0]];
#if HAVE_LIBLAPACK_MKL == 1
  k = m_pnDimSize[0]*int(ilaenv_(&ispec, name, opts, &N1, &N2, &N3, &N4));
#elif HAVE_LIBLAPACK_NETLIB == 1
  k = m_pnDimSize[0]*int(ilaenv_(&ispec, name, opts, &N1, &N2, &N3, &N4, namelen, optslen));
#elif HAVE_LIBLAPACK_LAPACK == 1
//  k = m_pnDimSize[0]*int(ilaenv_(&ispec, name, opts, &N1, &N2, &N3, &N4, namelen, optslen));
  k = m_pnDimSize[0]*int(ilaenv_(&ispec, name, opts, &N1, &N2, &N3, &N4));
#endif

  pdWorkSpace = new doublereal[unsigned(k)];

//  info = clapack_dgetrf(CblasColMajor, m_pnDimSize[0], m_pnDimSize[0], m_pdData, m_pnDimSize[0], ipiv);
  N1 = (integer)m_pnDimSize[0];
  dgetrf_(&N1, &N1, (doublereal *)m_pdData, &N1, ipiv, &info);
  if(int(info) != 0) {throw ELELapack; }

  if(pdLogDet != NULL) { // We need to calculate the determinant as well
    *pdLogDet = 0;
    for(int i=0; i<m_pnDimSize[0]; i++) {
      *pdLogDet += log(double(m_pdData[i+m_pnDimSize[0]*i]));
    } // for i
//    (*pdLogDet)*=2;
  } // if

//  info = clapack_dgetri(CblasColMajor, m_pnDimSize[0], m_pdData, m_pnDimSize[0], ipiv);
  dgetri_(&N1, (doublereal *)m_pdData, &N1, ipiv, pdWorkSpace, &k, &info);
  if(int(info) != 0) {throw ELELapack; }

  delete[] ipiv;
  delete[] pdWorkSpace;

#elif HAVE_LIBLAPACK_GSL == 1

  // USE THE GNU SCIENTIFIC LIBRARY
  CMatrix matLU = *this;
  int nSign;
  gsl_matrix_view mLU = gsl_matrix_view_array(matLU.m_pdData, matLU.m_pnDimSize[1], matLU.m_pnDimSize[0]);
  gsl_matrix_view mINV = gsl_matrix_view_array(m_pdData, m_pnDimSize[1], m_pnDimSize[0]);
  gsl_permutation *pPerm = gsl_permutation_alloc(m_pnDimSize[0]);

  // Decompose in Upper & Lower diagonal factors
  gsl_linalg_LU_decomp(&mLU.matrix, pPerm, &nSign);
  // Compute the determinant if needed
  if(pdLogDet != NULL) *pdLogDet = gsl_linalg_LU_lndet(&mLU.matrix);
  // Invert the matrix, using the LU decomposition
  gsl_linalg_LU_invert(&mLU.matrix, pPerm, &mINV.matrix);

  gsl_permutation_free(pPerm);
  // use the gnu scientific library

#endif // HAVE_LIBLAPACK

  return *this;
} // Invert

/* Calculates the inverse of a matrix using LU-factorization
   The original matrix is overwritten by this function
   
   This matrix assumes a positive definite, symmetric square matrix. LOOK OUT!
   DO NOT USE OTHERWISE
   */
CMatrix &CMatrix::InvertSVD() {
  if(! m_bDefined) {throw ELENotDefined; }
  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

#ifdef HAVE_LIBLAPACK
  CMatrix mdTemp, mdU, mdV, mdVTrans, mdSigma;
  CVector vdS;

  mdTemp = *this;

  mdTemp.SVD(mdU, mdV, vdS);

//  PrintVector(vdS);

  mdSigma.Initialize(mdU.m_pnDimSize[0], mdU.m_pnDimSize[1]);
  for(int i=0; i<mdU.m_pnDimSize[0]; ++i) {
    for(int j=0; j<mdU.m_pnDimSize[1]; ++j) {
      if(i == j && double(vdS[i]) > 0) {
	mdSigma[i][j] = 1.0 / double(vdS[i]);
      } else {
	mdSigma[i][j] = 0.0;
      } /* if i == j */
    } /* for j */
  } /* for i */
  mdVTrans = mdV[LO_TRANSPOSE];
  mdV = mdSigma * mdVTrans;
  *this = mdU * mdV;

#elif HAVE_LIBLAPACK_GSL == 1

  // USE THE GNU SCIENTIFIC LIBRARY
  // use the gnu scientific library

#endif // HAVE_LIBLAPACK

  return *this;
} // InvertSVD



/* Calculates the inverse of a matrix using Cholesky-decomposition
   The original matrix is overwritten by this function */
CMatrix &CMatrix::InvertChol(double *pdLogDet) {
  if(! m_bDefined) {throw ELENotDefined; }
  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

#ifdef HAVE_LIBLAPACK

  // Invert this matrix using the clapack/atlas routine for Cholesky
  doublereal *pdWorkSpace = (doublereal *) m_pdData;
  integer info, dim;
  char strBuf[2];
  strcpy(strBuf, "U");

  dim = (integer) (m_pnDimSize[0]);
  dpotrf_(strBuf, &dim, pdWorkSpace, &dim, &info);
  if(int(info) != 0) { printf("Cholesky: %d\n", (int) info); throw ELELapack; }

  if(pdLogDet != NULL) { // We need to calculate the determinant as well
    *pdLogDet = 0;
    for(int i=0; i<m_pnDimSize[0]; i++) {
      *pdLogDet += log(double(pdWorkSpace[i+int(dim)*i]));
    } // for i
    (*pdLogDet)*=2;
  } // if

  dpotri_(strBuf, &dim, pdWorkSpace, &dim, &info);
  if(int(info) != 0) { printf("Cholesky-Solve: %d\n", (int) info); throw ELELapack; }

  for(int i=1; i<m_pnDimSize[0]; i++) {
    for(int j=0; j<i; j++) {
      pdWorkSpace[i+int(dim)*j] = pdWorkSpace[j+int(dim)*i];
    } // for j
  } // for i

#elif HAVE_LIBLAPACK_GSL == 1

  // USE THE GNU SCIENTIFIC LIBRARY

  // Copy current matrix in matChol, *this will be linked to mInv wich will contain the inverse
  int info;
  CMatrix matChol = *this;
  gsl_matrix_view mChol = gsl_matrix_view_array(matChol.m_pdData, matChol.m_pnDimSize[0], matChol.m_pnDimSize[1]);
  gsl_matrix_view mInv = gsl_matrix_view_array(m_pdData, m_pnDimSize[0], m_pnDimSize[1]);
  gsl_vector_view vInvColumn;

  // Calculate the Cholesky decomposition
  info = gsl_linalg_cholesky_decomp(&mChol.matrix);
  if(info != 0) { printf("Cholesky: %d\n", (int) info); throw ELELapack; }

  if(pdLogDet != NULL) { // We need to calculate the determinant as well
    *pdLogDet = 0;
    for(int i=0; i<m_pnDimSize[0]; i++) {
      *pdLogDet += log(gsl_matrix_get(&mChol.matrix,i,i) );
    } // for i
    (*pdLogDet)*=2;
  } // if


  // Use the Cholesky decomposition to calculate the inverse column-by-column
  gsl_matrix_set_identity(&mInv.matrix);
  for(int i=0; i<m_pnDimSize[1]; i++) {
    vInvColumn = gsl_matrix_column(&mInv.matrix, i);
    info = gsl_linalg_cholesky_svx(&mChol.matrix,&vInvColumn.vector);
    if(info != 0) { printf("Cholesky-Solve: %d\n", (int) info); throw ELELapack; }
  } // for i


#endif // HAVE_LIBLAPACK

  return *this;
} // InvertChol

/* Calculates the inverse of a matrix using LU-factorization
   The original matrix is not overwritten by this function */
CMatrix &CMatrix::Inverse(double *pdLogDet) {
  if(! m_bDefined) {throw ELENotDefined; }
  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

#ifdef HAVE_LIBLAPACK

  // Invert this matrix using the clapack/atlas routine
  integer info, *ipiv, ispec=1, N1, N2, N3, N4, k;
  doublereal *pdWorkSpace;
#if HAVE_LIBLAPACK_MKL
  const char *name="dgetri", *opts="";
#else
  char *name="dgetri", *opts="";
#endif
  long namelen = 6, optslen=0;
  CMatrix *pmatInverse;
  
  pmatInverse = new CMatrix(*this);
  pmatInverse->SetAllocated(true);

  N1 = N2 = integer(pmatInverse->m_pnDimSize[0]);
  N3 = N4 = integer(-1);

  ipiv = new integer[pmatInverse->m_pnDimSize[0]];
#if HAVE_LIBLAPACK_MKL == 1
  k = pmatInverse->m_pnDimSize[0]*int(ilaenv_(&ispec, name, opts, &N1, &N2, &N3, &N4));
#elif HAVE_LIBLAPACK_NETLIB == 1
  k = pmatInverse->m_pnDimSize[0]*int(ilaenv_(&ispec, name, opts, &N1, &N2, &N3, &N4, namelen, optslen));
#elif HAVE_LIBLAPACK_LAPACK
  k = pmatInverse->m_pnDimSize[0]*int(ilaenv_(&ispec, name, opts, &N1, &N2, &N3, &N4));
//  k = pmatInverse->m_pnDimSize[0]*int(ilaenv_(&ispec, name, opts, &N1, &N2, &N3, &N4, namelen, optslen));
#endif
  pdWorkSpace = new doublereal[unsigned(k)];

//  info = clapack_dgetrf(CblasColMajor, m_pnDimSize[0], m_pnDimSize[0], m_pdData, m_pnDimSize[0], ipiv);
  N1 = (integer)pmatInverse->m_pnDimSize[0];
  dgetrf_(&N1, &N1, (doublereal *)pmatInverse->m_pdData, &N1, ipiv, &info);
  if(int(info) != 0) {throw ELELapack; }

  if(pdLogDet != NULL) { // We need to calculate the determinant as well
    *pdLogDet = 0;
    for(int i=0; i<m_pnDimSize[0]; i++) {
      *pdLogDet += log(double(pmatInverse->m_pdData[i+m_pnDimSize[0]*i]));
    } // for i
//    (*pdLogDet)*=2;
  } // if

//  info = clapack_dgetri(CblasColMajor, m_pnDimSize[0], m_pdData, m_pnDimSize[0], ipiv);
  dgetri_(&N1, (doublereal *)pmatInverse->m_pdData, &N1, ipiv, pdWorkSpace, &k, &info);
  if(int(info) != 0) {throw ELELapack; }

  delete[] ipiv;
  delete[] pdWorkSpace;

#elif HAVE_LIBLAPACK_GSL == 1

  // USE THE GNU SCIENTIFIC LIBRARY
  CMatrix matLU = *this, *pmatInverse;
  pmatInverse = new CMatrix(*this);
  pmatInverse->SetAllocated(true);

  int nSign;
  gsl_matrix_view mLU = gsl_matrix_view_array(matLU.m_pdData, matLU.m_pnDimSize[1], matLU.m_pnDimSize[0]);
  gsl_matrix_view mINV = gsl_matrix_view_array(pmatInverse->m_pdData, m_pnDimSize[1], m_pnDimSize[0]);
  gsl_permutation *pPerm = gsl_permutation_alloc(m_pnDimSize[0]);

  // Decompose in Upper & Lower diagonal factors
  gsl_linalg_LU_decomp(&mLU.matrix, pPerm, &nSign);
  // Compute the determinant if needed
  if(pdLogDet != NULL) *pdLogDet = gsl_linalg_LU_lndet(&mLU.matrix);
  // Invert the matrix, using the LU decomposition
  gsl_linalg_LU_invert(&mLU.matrix, pPerm, &mINV.matrix);

  gsl_permutation_free(pPerm);
  // use the gnu scientific library
#endif // HAVE_LIBLAPACK


  return *pmatInverse;
} // Inverse


/* Calculates the inverse of a matrix using Cholesky-decomposition
   The original matrix is not overwritten by this function */
CMatrix &CMatrix::InverseChol(double *pdLogDet) {
  if(! m_bDefined) {throw ELENotDefined; }
  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

  // The inverse of the matrix will be allocated here
  CMatrix *pmatInverse;
  
  pmatInverse = new CMatrix(*this);
  pmatInverse->SetAllocated(true);

#ifdef HAVE_LIBLAPACK

  // Invert this matrix using the clapack/atlas routine for Cholesky
  doublereal *pdWorkSpace = (doublereal *) pmatInverse->m_pdData;
  integer info, dim;
  char strBuf[2];
  strcpy(strBuf, "U");

  dim = (integer) (pmatInverse->m_pnDimSize[0]);
  dpotrf_(strBuf, &dim, pdWorkSpace, &dim, &info);
  if(int(info) != 0) { printf("Cholesky: %d\n", (int) info); throw ELELapack; }

  if(pdLogDet != NULL) { // We need to calculate the determinant as well
    *pdLogDet = 0;
    for(int i=0; i<pmatInverse->m_pnDimSize[0]; i++) {
      *pdLogDet += log(double(pdWorkSpace[i+int(dim)*i]));
    } // for i
    (*pdLogDet)*=2;
  } // if

  dpotri_(strBuf, &dim, pdWorkSpace, &dim, &info);
  if(int(info) != 0) {
    printf("Cholesky-Solve: %d\n", (int) info);
    printf("U[%d][%d]: %f\n", (int) info, (int) info, pdWorkSpace[(int) info + int(dim)*int(info)]);
    throw ELELapack;
  }

  for(int i=1; i<pmatInverse->m_pnDimSize[0]; i++) {
    for(int j=0; j<i; j++) {
      pdWorkSpace[i+int(dim)*j] = pdWorkSpace[j+int(dim)*i];
    } // for j
  } // for i

#elif HAVE_LIBLAPACK_GSL == 1

  // USE THE GNU SCIENTIFIC LIBRARY

  // Copy current matrix in matChol, *this will be linked to mInv wich will contain the inverse
  int info;
  CMatrix matChol = *this;
  gsl_matrix_view mChol = gsl_matrix_view_array(matChol.m_pdData, matChol.m_pnDimSize[0], matChol.m_pnDimSize[1]);
  gsl_matrix_view mInv = gsl_matrix_view_array(pmatInverse->m_pdData, pmatInverse->m_pnDimSize[0], pmatInverse->m_pnDimSize[1]);
  gsl_vector_view vInvColumn;

  // Calculate the Cholesky decomposition
  info = gsl_linalg_cholesky_decomp(&mChol.matrix);
  if(info != 0) { printf("Cholesky: %d\n", (int) info); throw ELELapack; }

  if(pdLogDet != NULL) { // We need to calculate the determinant as well
    *pdLogDet = 0;
    for(int i=0; i<pmatInverse->m_pnDimSize[0]; i++) {
      *pdLogDet += log(gsl_matrix_get(&mChol.matrix,i,i) );
    } // for i
    (*pdLogDet)*=2;
  } // if


  // Use the Cholesky decomposition to calculate the inverse column-by-column
  gsl_matrix_set_identity(&mInv.matrix);
  for(int i=0; i<m_pnDimSize[1]; i++) {
    vInvColumn = gsl_matrix_column(&mInv.matrix, i);
    info = gsl_linalg_cholesky_svx(&mChol.matrix,&vInvColumn.vector);
    if(info != 0) { printf("Cholesky-Solve: %d\n", (int) info); throw ELELapack; }
  } // for i


#endif // HAVE_LIBLAPACK

  return *pmatInverse;
} // InverseChol


CVector &CMatrix::Eigenvalues() {
// pcJob == "V"/"N". Calculate Eigenvalues/Eigenvectors for "V". Only eigenvalues for "N"
// Output is CLinearObject array of CVector * or CVector * and CMatrix for "N" or "V" option.
// Make sure to dynamically delete the output of this function
// CLinearObject *CMatrix::Eigen(const char *pcJob) {
  CVector *cvEigenValues = new CVector;
  CLinearObject *plo;

  // Compute the eigenvalues
  plo = Eigen("N");
  *cvEigenValues = *((CVector *)&plo[0]);
  delete[] plo;
  return *cvEigenValues;
} // Eigenvalues

CMatrix &CMatrix::Eigenvectors() {
  CMatrix *cmEigenVectors = new CMatrix;
  CLinearObject *plo;

  // Compute the eigenvalues
  plo = Eigen("V");
  *cmEigenVectors = *((CMatrix *)&plo[1]);
  delete[] plo;
  return *cmEigenVectors;
} // Eigenvectors


// In Eigen wordt deze lapack driver gebruikt
/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   

    Purpose   
    =======   

    DSYEVD computes all eigenvalues and, optionally, eigenvectors of a   
    real symmetric matrix A. If eigenvectors are desired, it uses a   
    divide and conquer algorithm.   

    The divide and conquer algorithm makes very mild assumptions about   
    floating point arithmetic. It will work on machines with a guard   
    digit in add/subtract, or on those binary machines without guard   
    digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or   
    Cray-2. It could conceivably fail on hexadecimal or decimal machines   
    without guard digits, but we know of none.   

    Because of large use of BLAS of level 3, DSYEVD needs N**2 more   
    workspace than DSYEVX.   

    Arguments   
    =========   

    JOBZ    (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only;   
            = 'V':  Compute eigenvalues and eigenvectors.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the   
            leading N-by-N upper triangular part of A contains the   
            upper triangular part of the matrix A.  If UPLO = 'L',   
            the leading N-by-N lower triangular part of A contains   
            the lower triangular part of the matrix A.   
            On exit, if JOBZ = 'V', then if INFO = 0, A contains the   
            orthonormal eigenvectors of the matrix A.   
            If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')   
            or the upper triangle (if UPLO='U') of A, including the   
            diagonal, is destroyed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    W       (output) DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, the eigenvalues in ascending order.   

    WORK    (workspace/output) DOUBLE PRECISION array,   
                                           dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If N <= 1,               LWORK must be at least 1.   
            If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.   
            If JOBZ = 'V' and N > 1, LWORK must be at least   
                                                  1 + 6*N + 2*N**2.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    IWORK   (workspace/output) INTEGER array, dimension (LIWORK)   
            On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.   

    LIWORK  (input) INTEGER   
            The dimension of the array IWORK.   
            If N <= 1,                LIWORK must be at least 1.   
            If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.   
            If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.   

            If LIWORK = -1, then a workspace query is assumed; the   
            routine only calculates the optimal size of the IWORK array,   
            returns this value as the first entry of the IWORK array, and   
            no error message related to LIWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   





            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the algorithm failed to converge; i   
                  off-diagonal elements of an intermediate tridiagonal   
                  form did not converge to zero
*/

// pcJob == "V"/"N". Calculate Eigenvalues/Eigenvectors for "V". Only eigenvalues for "N"
// Output is CLinearObject array of CVector * or CVector * and CMatrix for "N" or "V" option.
// Make sure to dynamically delete the output of this function
CLinearObject *CMatrix::Eigen(const char *pcJob) {
  CMatrix *pEVMat;
  CVector *pEVVec;
  CLinearObject *ploOutput;

  if(! m_bDefined) {throw ELENotDefined; }
  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }


#ifdef HAVE_LIBLAPACK

  // Compute the eigenvalues of the matrix using the clapack/atlas routine
  double *pdEigenVectors, *pdEigenValues, *pdWorkSpace;
  int n, n2, nWork, nWork2, *pnDimSize;
//  char *strUpper="U";
  integer N1, N2, N3, info, *pnWorkSpace;
  char strJob[2], strUpper[2];
  strcpy(strJob, pcJob);
  strcpy(strUpper, "U");

  n2 = int(m_pnDimSize[0]*m_pnDimSize[1]);
  n = int(m_pnDimSize[0]);
  nWork = 1 + 6*n + 2*n2;
  nWork2 = 3 + 5*n;
  pdEigenVectors = new double[n2];
  pdWorkSpace = new double[nWork];
  pnWorkSpace = new integer[nWork2];
  pdEigenValues = new double[m_pnDimSize[0]];
  for(int i=0; i<n2; i++) { pdEigenVectors[i] = m_pdData[i]; }

  N1 = (integer)m_pnDimSize[0]; N2 = nWork; N3 = nWork2;
  dsyevd_(strJob, strUpper, &N1, (doublereal *) pdEigenVectors,&N1, (doublereal *)pdEigenValues, (doublereal *)pdWorkSpace, &N2, pnWorkSpace, &N3, &info);

  if(int(info) != 0) { throw ELELapack; }
  switch(pcJob[0]) {
  case 'V': // Compute eigenvalues and eigenvectors
    ploOutput = new CLinearObject[2];
    pEVVec = (CVector *)&ploOutput[0];
    pEVMat = (CMatrix *)&ploOutput[1];
    
    pnDimSize = new int[1];
    pnDimSize[0] = n;
    pEVVec->SetRawData(pdEigenValues, 1, pnDimSize, false, EVector);
    
    pnDimSize = new int[2];
    pnDimSize[0] = pnDimSize[1] = n;
    pEVMat->SetRawData(pdEigenVectors, 2, pnDimSize, false, EMatrix);

    break;
  case 'N': // Compute eigenvalues only
    ploOutput = new CLinearObject[1];
    pEVVec = (CVector *)&ploOutput[0];
//    pEVMat = (CMatrix *)ploOutput[1];
    
    pnDimSize = new int[1];
    pnDimSize[0] = n;
    pEVVec->SetRawData(pdEigenValues, 1, pnDimSize, false, EVector);
    
//    pnDimSize = new int[2];
//    pnDimSize[0] = pnDimSize[1] = n;
//   pEVMat->SetRawData(pdEigenVectors, 2, pnDimSize, true, EMatrix);

    delete[] pdEigenVectors;
    break;
  default:
    break;
  } // switch
  
  delete[] pnWorkSpace;
  delete[] pdWorkSpace;

#elif HAVE_LIBLAPACK_GSL == 1
  // USE THE GNU SCIENTIFIC LIBRARY
  // Transpose this matrix, since the GSL stores diffently from LAPACK
  CMatrix mA = (*this)[LO_TRANSPOSE];
  gsl_eigen_symmv_workspace *wv;
  gsl_eigen_symm_workspace *w;
  gsl_matrix_view m;
  gsl_matrix_view mEV;
  gsl_vector_view vEV;

  switch(pcJob[0]) {
  case 'V': // Compute eigenvalues and eigenvectors
    // Allocate memory
    ploOutput = new CLinearObject[2];
    pEVVec = (CVector *)&ploOutput[0];
    pEVMat = (CMatrix *)&ploOutput[1];

    pEVVec->Initialize(mA.m_pnDimSize[0]);
    pEVVec->SetAllocated(false);
    pEVMat->Initialize(mA.m_pnDimSize[0], mA.m_pnDimSize[1]);
    pEVMat->SetAllocated(false);

    // Set the gsl matrix/vector views
    wv = gsl_eigen_symmv_alloc(mA.m_pnDimSize[0]);
    m = gsl_matrix_view_array(mA.m_pdData, mA.m_pnDimSize[0], mA.m_pnDimSize[1]);
    mEV = gsl_matrix_view_array(pEVMat->m_pdData, mA.m_pnDimSize[0], mA.m_pnDimSize[1]);
    vEV = gsl_vector_view_array(pEVVec->m_pdData, mA.m_pnDimSize[0]);

    // Compute the eigenvalues/vectors
    gsl_eigen_symmv(&m.matrix, &vEV.vector, &mEV.matrix, wv);

    // We need to transpose the eigenvector matrix for compatibility with lapack
    // routines.
    gsl_matrix_transpose(&mEV.matrix);

    // Free the workspace
    gsl_eigen_symmv_free(wv);

    // Sort the eigenvalues/vectors
    // gsl_eigen_symmv_sort(&vEV.vector, &mEV.matrix, GSL_EIGEN_SORT_ABS_ASC);
    break;
  case 'N': // Compute eigenvalues only
    // Allocate memory
    ploOutput = new CLinearObject[1];
    pEVVec = (CVector *)&ploOutput[0];

    pEVVec->Initialize(mA.m_pnDimSize[0]);
    pEVVec->SetAllocated(false);

    // Set the gsl vector views
    w = gsl_eigen_symm_alloc(mA.m_pnDimSize[0]);
    m = gsl_matrix_view_array(mA.m_pdData, mA.m_pnDimSize[0], mA.m_pnDimSize[1]);
    vEV = gsl_vector_view_array(pEVVec->m_pdData, mA.m_pnDimSize[0]);

    // Compute the eigenvalues
    gsl_eigen_symm(&m.matrix, &vEV.vector, w);

    // Free the workspace
    gsl_eigen_symm_free(w);

    // Sort the eigenvalues/vectors
    // gsl_eigen_symm_sort(&vEV.vector, GSL_EIGEN_SORT_ABS_ASC);
    break;
  default:
    break;
  } // switch

  // use the gnu scientific library
#endif // HAVE_LIBLAPACK
  
  return ploOutput;
} // Eigen




/*
    =======   

    DPOTRF computes the Cholesky factorization of a real symmetric   
    positive definite matrix A.   

    The factorization has the form   
       A = U**T * U,  if UPLO = 'U', or   
       A = L  * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the block version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading   
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower   
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U**T*U or A = L*L**T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   
*/
// Compute the Cholesky factorization of a real-symmetric matrix
CMatrix &CMatrix::Cholesky() {
  CMatrix *cmCholesky = new CMatrix(*this);
  cmCholesky->SetAllocated(true);

  if(! m_bDefined) {throw ELENotDefined; }
  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

#ifdef HAVE_LIBLAPACK

  doublereal *pdWorkSpace = (doublereal *) cmCholesky->m_pdData;
  integer info, dim;
  char strBuf[2];
  strcpy(strBuf, "L");

  dim = (integer) (m_pnDimSize[0]);
  dpotrf_(strBuf, &dim, pdWorkSpace, &dim, &info);
  if(int(info) != 0) { printf("Cholesky: %d\n", (int) info); throw ELELapack; }

  // Lower-diagonal cholesky. Set other elements equal to zero
  for(int i=0; i<m_pnDimSize[0]; i++) {
    for(int j=i+1; j<m_pnDimSize[0]; j++) {
      pdWorkSpace[i + m_pnDimSize[0]*j] = 0;
    } // for j
  } // for i

#elif HAVE_LIBLAPACK_GSL == 1 // HAVE_LIBLAPACK
  // USE THE GNU SCIENTIFIC LIBRARY
  int info;
  gsl_matrix_view m = gsl_matrix_view_array(cmCholesky->m_pdData, cmCholesky->m_pnDimSize[0], cmCholesky->m_pnDimSize[1]);
  info = gsl_linalg_cholesky_decomp(&m.matrix);
  if(info != 0) { printf("Cholesky: %d\n", (int) info); throw ELELapack; }

  // Lower-diagonal cholesky. Set other elements equal to zero
  for(int i=0; i<m_pnDimSize[0]; i++) {
    for(int j=i+1; j<m_pnDimSize[0]; j++) {
      gsl_matrix_set(&m.matrix, j, i, 0);
    } // for j
  } // for i

  // use the gnu scientific library
#endif // HAVE_LIBLAPACK

  // Return the factorization
  return *cmCholesky;
}  // Cholesky


// Compute the Cholesky factorization of a real-symmetric matrix
CMatrix &CMatrix::CholeskySingle() {
  CMatrix *cmCholesky = new CMatrix(*this);
  cmCholesky->SetAllocated(true);

  if(! m_bDefined) {throw ELENotDefined; }
  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }


#ifdef HAVE_LIBLAPACK

//  doublereal *pdWorkSpace = (doublereal *) cmCholesky->m_pdData;
  real *pfWorkSpace;
  integer info, dim;
  char strBuf[2];
  strcpy(strBuf, "L");

  pfWorkSpace = new real[int(m_pnDimSize[0]*m_pnDimSize[1])];

  for(int i=0; i<m_pnDimSize[0]*m_pnDimSize[1]; i++) pfWorkSpace[i] = real(m_pdData[i]);

  dim = (integer) (m_pnDimSize[0]);
  spotrf_(strBuf, &dim, pfWorkSpace, &dim, &info);
  if(int(info) != 0) { printf("\n\n\n\n\n\nCholesky: %d\n", (int) info); throw ELELapack; }

  delete[] pfWorkSpace;

#elif HAVE_LIBLAPACK_GSL == 1

  // USE THE GNU SCIENTIFIC LIBRARY
  int info;
  gsl_matrix_view m = gsl_matrix_view_array(cmCholesky->m_pdData, cmCholesky->m_pnDimSize[0], cmCholesky->m_pnDimSize[1]);

  info = gsl_linalg_cholesky_decomp(&m.matrix);
  if(info != 0) { printf("Cholesky: %d\n", (int) info); throw ELELapack; }

  // Lower-diagonal cholesky. Set other elements equal to zero
  for(int i=0; i<m_pnDimSize[0]; i++) {
    for(int j=i+1; j<m_pnDimSize[0]; j++) {
      gsl_matrix_set(&m.matrix, j, i, 0);
    } // for j
  } // for i
  // use the gnu scientific library
#endif // HAVE_LIBLAPACK

  // Return the factorization
  return *cmCholesky;
}  // Cholesky


// Return the Q-factor of QR decomposition
CMatrix &CMatrix::QRDecompQ() {
  CMatrix *cmCholesky = new CMatrix(*this);
  cmCholesky->SetAllocated(true);

  if(! m_bDefined) {throw ELENotDefined; }
//  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

#if 0

#ifdef HAVE_LIBLAPACK

//  doublereal *pdWorkSpace = (doublereal *) cmCholesky->m_pdData;
  real *pfWorkSpace;
  integer info, dim;
  char strBuf[2];
  strcpy(strBuf, "U");

  pfWorkSpace = new real[int(m_pnDimSize[0]*m_pnDimSize[1])];

  for(int i=0; i<m_pnDimSize[0]*m_pnDimSize[1]; i++) pfWorkSpace[i] = real(m_pdData[i]);

  dim = (integer) (m_pnDimSize[0]);
  spotrf_(strBuf, &dim, pfWorkSpace, &dim, &info);
  if(int(info) != 0) { printf("\n\n\n\n\n\nCholesky: %d\n", (int) info); throw ELELapack; }

  delete[] pfWorkSpace;

#elif HAVE_LIBLAPACK_GSL == 1

  // USE THE GNU SCIENTIFIC LIBRARY
  int info;
  gsl_matrix_view m = gsl_matrix_view_array(cmCholesky->m_pdData, cmCholesky->m_pnDimSize[0], cmCholesky->m_pnDimSize[1]);
  info = gsl_linalg_cholesky_decomp(&m.matrix);

  if(info != 0) { printf("Cholesky: %d\n", (int) info); throw ELELapack; }
  // use the gnu scientific library

#endif // HAVE_LIBLAPACK

#endif
  // Return the factorization
  return *cmCholesky;
}  // QRDecompQ


/* This member function calculates the singular value decomposition of the
 * matrix *this.
 * M = U Sigma Vt
 * M is an mxn matrix (this matrix)
 * U is an mxm orthogonal matrix
 * V is an nxn orthogonal matrix (Vt is the transpose of V)
 * Sigma is a diagonal mxn matrix.
 *
 * Vector S vdS is the diagonal elements of Sigma
 * */
void CMatrix::SVD(CMatrix &mdU, CMatrix &mdV, CVector &vdS) {
  int nLWork, nInfo, nLdVt, nLdA, nLdU, m, n;
  double *pdWork=NULL, dWorkOpt;
  CMatrix mdTemp;
  char strBuf[4];
  strcpy(strBuf, "All");

  m = this->m_pnDimSize[0];
  n = this->m_pnDimSize[1];

  nLdA = m;
  nLdU = m;
  nLdVt = n;
  mdU.Initialize(m,m);
  mdTemp.Initialize(n,n);
  vdS.Initialize(n < m ? n : m);

  nLWork = -1;

#ifdef HAVE_LIBLAPACK

  // Find the optimal workspace
  dgesvd_(strBuf, strBuf, &m, &n, this->m_pdData, &nLdA, vdS.m_pdData, mdU.m_pdData, &nLdU, mdTemp.m_pdData, &nLdVt, &dWorkOpt, &nLWork, &nInfo);

  nLWork = (int) dWorkOpt;
  pdWork = new double[nLWork];
  dgesvd_(strBuf, strBuf, &m, &n, this->m_pdData, &nLdA, vdS.m_pdData, mdU.m_pdData, &nLdU, mdTemp.m_pdData, &nLdVt, pdWork, &nLWork, &nInfo);

  mdV = mdTemp[LO_TRANSPOSE];

  delete[] pdWork;
#endif // HAVE_LIBLAPACK
} // SVD


// Compute the LU-factorization of a real general matrix
CMatrix &CMatrix::LUFact() {
  CMatrix *cmLUFact = new CMatrix(*this);
  cmLUFact->SetAllocated(true);

  if(! m_bDefined) {throw ELENotDefined; }
  if(! m_pnDimSize[0] == m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

#ifdef HAVE_LIBLAPACK

  integer info, *ipiv, N;
  ipiv = new integer[cmLUFact->m_pnDimSize[0]];

  N = (integer)cmLUFact->m_pnDimSize[0];

  dgetrf_(&N, &N, (doublereal *)cmLUFact->m_pdData, &N, ipiv, &info);
  if(int(info) != 0) {throw ELELapack; }

#elif HAVE_LIBLAPACK_GSL == 1

  // USE THE GNU SCIENTIFIC LIBRARY
  int info, nSign;
  gsl_matrix_view m = gsl_matrix_view_array(cmLUFact->m_pdData, cmLUFact->m_pnDimSize[0], cmLUFact->m_pnDimSize[1]);
  gsl_permutation *pPerm = gsl_permutation_alloc(cmLUFact->m_pnDimSize[0]);

  // Decompose in Upper & Lower diagonal factors
  info = gsl_linalg_LU_decomp(&m.matrix, pPerm, &nSign);

  gsl_permutation_free(pPerm);

  if(info != 0) {throw ELELapack; }

  // use the gnu scientific library
#endif // HAVE_LIBLAPACK

  // Return the factorization
  return *cmLUFact;
} // LUFact()




// Add a number (element-wise)
CMatrix &CMatrix::operator+=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] += d;
  } // for
  return *this;
} // operator+

// Subtract a number (element-wise)
CMatrix &CMatrix::operator-=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] -= d;
  } // for
  return *this;
} // operator-

// Multiply with a number (element-wise)
CMatrix &CMatrix::operator*=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] = m_pdData[i]*d;
  } // for
  return *this;
} // operator*

// Devide by a number (element-wise)
CMatrix &CMatrix::operator/=(double d) {
  int nElements=1;  // The amount of data-elements in this object
  if(! m_bDefined) {throw ELENotDefined; }
  
  for(int i=0; i<m_nDimensions; i++) {        // Go through each dimension
    nElements = nElements*m_pnDimSize[i];
  } // for
  for(int i=0; i<nElements; i++) {          // Calculate all new elements
    m_pdData[i] = m_pdData[i]/d;
  } // for
  return *this;
} // operator/

// Add a number (element-wise)
CMatrix &CMatrix::operator+(double d) {
  CMatrix *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] += d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CMatrix;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i] + d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator+

// Subtract a number (element-wise)
CMatrix &CMatrix::operator-(double d) {
  CMatrix *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] -= d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CMatrix;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i] - d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator-

// Multiply with a number (element-wise)
CMatrix &CMatrix::operator*(double d) {
  CMatrix *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] = m_pdData[i]*d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CMatrix;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i]*d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator*

// Devide by a number (element-wise)
CMatrix &CMatrix::operator/(double d) {
  CMatrix *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! m_bDefined) {throw ELENotDefined; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] = m_pdData[i]/d;
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CMatrix;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i]/d;
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  return *newob;
} // operator/

// Add a number (element-wise)
CMatrix &CMatrix::operator+(CMatrix &mat) {
  CMatrix *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! (m_bDefined && mat.Defined())) {throw ELENotDefined; }
  if(m_pnDimSize[0] != mat.m_pnDimSize[0] || m_pnDimSize[1] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] += mat.m_pdData[i];
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CMatrix;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i] + mat.m_pdData[i];
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  if(mat.m_bAllocated) {delete &mat; } // This should call for the destructor
  return *newob;
} // operator+

// Subtract a number (element-wise)
CMatrix &CMatrix::operator-(CMatrix &mat) {
  CMatrix *newob;
  double *pdNewData;
  int *pnDimSize;
  int nElements=this->GetElements();

  if(! (m_bDefined && mat.Defined())) {throw ELENotDefined; }
  if(m_pnDimSize[0] != mat.m_pnDimSize[0] || m_pnDimSize[1] != mat.m_pnDimSize[1]) {throw ELEDimensionMisMatch; }

  // Check whether the this object is also allocated. Ass is left to right!!!!!
  switch(m_bAllocated) {
  case true:
    // If this object was allocated, then we just substitute the data
    for(int i=0; i<nElements; i++) {
      m_pdData[i] -= mat.m_pdData[i];
    } // for
    newob = this;
    break;
  case false:
    // If this object was not allocated, we must make a copy that is allocated
    newob = new CMatrix;
    pdNewData = new double[nElements];

    for(int i=0; i<nElements; i++) {
      pdNewData[i] = m_pdData[i] - mat.m_pdData[i];
    } // for

    pnDimSize = new int[m_nDimensions];
    for(int i=0; i<m_nDimensions; i++) {
      pnDimSize[i] = m_pnDimSize[i];
    }

    newob->SetRawData(pdNewData, m_nDimensions, pnDimSize, true, m_eType);
   break;
  default:
    break;
  } // switch

  if(mat.m_bAllocated) {delete &mat; } // This should call for the destructor
  return *newob;
} // operator-

/* Returns a random square orthogonal matrix Q by using the QR
   decomposition on a random square matrix A=QR. */
CMatrix &RandomOrthogonalMatrix(int n) {
  CMatrix matA(n, n), matR(n,n), *pmatQ = new CMatrix(n,n);
  CVector vecTau(n);
  pmatQ->SetAllocated(true);

  // Generate a random matrix (elements between 0 and 1)
  matA.Randomize();

#ifdef HAVE_LIBLAPACK
  integer M, N, LDA, *JPVT, LWORK, NB, NB1, NB2, info, ispec=(integer)1, N1, N2, N3, N4;
  doublereal *A, *TAU, *WORK;
#if HAVE_LIBLAPACK_MKL
  const char *name1="dgeqp3", *name2="dorgqr" , *opts="";
#else
  char *name1="dgeqp3", *name2="dorgqr" , *opts="";
#endif
  long namelen = 6, optslen=0;

  M = (integer) n;
  N = (integer) n;
  A = (doublereal *) matA.m_pdData;
  LDA = (integer) n;
  N1 = M; N2=N; N3=(integer)-1; N4=(integer)-1;

  // Compute the optimal blocksize
#if HAVE_LIBLAPACK_MKL == 1
  NB1 = ilaenv_(&ispec, name1, opts, &N1, &N2, &N3, &N4);
  NB2 = ilaenv_(&ispec, name2, opts, &N1, &N2, &N3, &N4);
#elif HAVE_LIBLAPACK_NETLIB == 1
  NB1 = ilaenv_(&ispec, name1, opts, &N1, &N2, &N3, &N4, namelen, optslen);
  NB2 = ilaenv_(&ispec, name2, opts, &N1, &N2, &N3, &N4, namelen, optslen);
#elif HAVE_LIBLAPACK_LAPACK
  NB1 = ilaenv_(&ispec, name1, opts, &N1, &N2, &N3, &N4);
  NB2 = ilaenv_(&ispec, name2, opts, &N1, &N2, &N3, &N4);
//  NB1 = ilaenv_(&ispec, name1, opts, &N1, &N2, &N3, &N4, namelen, optslen);
//  NB2 = ilaenv_(&ispec, name2, opts, &N1, &N2, &N3, &N4, namelen, optslen);
#endif

  if(NB1 > NB2) NB = NB1; else NB = NB2;

  JPVT = new integer[int(N)];
  TAU = new doublereal[int(N)];
  LWORK = integer(2)*N+(N+integer(1))*NB;
  WORK = new doublereal[int(LWORK)];

  dgeqp3_(&M, &N, A, &LDA, JPVT, TAU, WORK, &LWORK, &info);
  if(int(info) != 0) {throw ELELapack; }

  dorgqr_(&M, &N, &N, A, &LDA, TAU, WORK, &LWORK, &info);
  if(int(info) != 0) {throw ELELapack; }

  // This is inefficient!
  *pmatQ = matA;

  delete[] JPVT;
  delete[] WORK;
  delete[] TAU;
#elif HAVE_LIBLAPACK_GSL == 1
  gsl_matrix_view mA = gsl_matrix_view_array(matA.m_pdData, matA.m_pnDimSize[0], matA.m_pnDimSize[1]);
  gsl_matrix_view mR = gsl_matrix_view_array(matR.m_pdData, matR.m_pnDimSize[0], matR.m_pnDimSize[1]);
  gsl_matrix_view mQ = gsl_matrix_view_array(pmatQ->m_pdData, pmatQ->m_pnDimSize[0], pmatQ->m_pnDimSize[1]);
  gsl_vector_view vTau = gsl_vector_view_array(vecTau.m_pdData, vecTau.m_pnDimSize[0]);

  // Use the QR decomposition
  gsl_linalg_QR_decomp(&mA.matrix, &vTau.vector);
  gsl_linalg_QR_unpack(&mA.matrix, &vTau.vector, &mQ.matrix, &mR.matrix);
#endif

  return (*pmatQ);
} // RandomOrthogonalMatrix



/*******************************************************************************
*  Copyright (C) 2009-2011 Intel Corporation. All Rights Reserved.
*  The information and material ("Material") provided below is owned by Intel
*  Corporation or its suppliers or licensors, and title to such Material remains
*  with Intel Corporation or its suppliers or licensors. The Material contains
*  proprietary information of Intel or its suppliers and licensors. The Material
*  is protected by worldwide copyright laws and treaty provisions. No part of
*  the Material may be copied, reproduced, published, uploaded, posted,
*  transmitted, or distributed in any way without Intel's prior express written
*  permission. No license under any patent, copyright or other intellectual
*  property rights in the Material is granted to or conferred upon you, either
*  expressly, by implication, inducement, estoppel or otherwise. Any license
*  under such intellectual property rights must be express and approved by Intel
*  in writing.
*
********************************************************************************
*/
/*
   LAPACKE_dgesvd Example.
   =======================

   Program computes the singular value decomposition of a general
   rectangular matrix A:

     8.79   9.93   9.83   5.45   3.16
     6.11   6.91   5.04  -0.27   7.98
    -9.15  -7.93   4.86   4.85   3.01
     9.57   1.64   8.83   0.74   5.80
    -3.49   4.02   9.80  10.00   4.27
     9.84   0.15  -8.99  -6.02  -5.31

   Description.
   ============

   The routine computes the singular value decomposition (SVD) of a real
   m-by-n matrix A, optionally computing the left and/or right singular
   vectors. The SVD is written as

   A = U*SIGMA*VT

   where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
   diagonal elements, U is an m-by-m orthogonal matrix and VT (V transposed)
   is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
   are the singular values of A; they are real and non-negative, and are
   returned in descending order. The first min(m, n) columns of U and V are
   the left and right singular vectors of A.

   Note that the routine returns VT, not V.

   Example Program Results.
   ========================

 LAPACKE_dgesvd (row-major, high-level) Example Program Results

 Singular values
  27.47  22.64   8.56   5.99   2.01

 Left singular vectors (stored columnwise)
  -0.59   0.26   0.36   0.31   0.23
  -0.40   0.24  -0.22  -0.75  -0.36
  -0.03  -0.60  -0.45   0.23  -0.31
  -0.43   0.24  -0.69   0.33   0.16
  -0.47  -0.35   0.39   0.16  -0.52
   0.29   0.58  -0.02   0.38  -0.65

 Right singular vectors (stored rowwise)
  -0.25  -0.40  -0.69  -0.37  -0.41
   0.81   0.36  -0.25  -0.37  -0.10
  -0.26   0.70  -0.22   0.39  -0.49
   0.40  -0.45   0.25   0.43  -0.62
  -0.22   0.14   0.59  -0.63  -0.44
*/

#ifdef HAVE_LIBLAPACK

#define min(a,b) ((a)>(b)?(b):(a))

/* Parameters */
#define M 6
#define N 5
#define LDA M
#define LDU M
#define LDVT N

#ifndef MKL_INT
#define MKL_INT int
#endif


/* Auxiliary routine: printing a matrix */
void print_matrix(const char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
        MKL_INT i, j;
        printf( "\n %s\n", desc );
//        for( i = 0; i < m; i++ ) {
        for( i = 0; i < n; i++ ) {
//                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
                for( j = 0; j < m; j++ ) printf( " %6.2f", a[i*lda+j] );
                printf( "\n" );
        }
}


/* Main program */
int testingsvd() {
/* Locals */
  // This is column-major. Need row-major. So transpose it! M6 N5 LDAM LDUM
  // LDVTN
        int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info, lwork;
        double wkopt;
        double* work;
        /* Local arrays */
        double s[N], u[LDU*M], vt[LDVT*N];
	char strBuf[4];
        double a[LDA*N] = {
            8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
            9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
            9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
            5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
            3.16,  7.98,  3.01,  5.80,  4.27, -5.31
        };
	strcpy(strBuf, "All");

        print_matrix( "Matrix m", m, n, a, lda );
        /* Executable statements */
        printf( " DGESVD Example Program Results\n" );
        /* Query and allocate the optimal workspace */
        lwork = -1;
        dgesvd_(strBuf, strBuf, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info);
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        /* Compute SVD */
        dgesvd_(strBuf, strBuf, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
        /* Check for convergence */
        if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge.\n" );
		return 1;
        }
	/* Print matrix a */
        print_matrix( "Matrix m", m, n, a, lda );
        /* Print singular values */
        print_matrix( "Singular values", 1, n, s, 1 );
        /* Print left singular vectors */
        print_matrix( "Left singular vectors (stored columnwise)", m, n, u, ldu );
        /* Print right singular vectors */
        print_matrix( "Right singular vectors (stored rowwise)", n, n, vt, ldvt );
        /* Free workspace */
        free( (void*)work );
	return 0;
} /* End of LAPACKE_dgesvd Example */

#endif
