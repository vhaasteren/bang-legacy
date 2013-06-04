/* linal.h   -- linear algebra wrappers for libraries as gsl and lapack

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



#ifndef NULL
  #define NULL      0
#endif

#ifndef __LINAL_H__
#define __LINAL_H__

// Define the constants that can be passed as indices to make
// certain operations possible. Stored in m_dFlags
// Note: not all these are/might ever be usable
#define LO_ALL          1.0
#define LO_TRANSPOSE    2.0
#define LO_INVERSE      4.0
#define LO_SQRT         8.0
#define LO_LOG          16.0



// Enumeration type to identify the class type
enum ELinearClassType {
  ELinearObject, ENumber, EVector, EMatrix, ETensor
};

// Enumeration type to classify different errors
enum ELinearError {
  ELENotDefined, ELEWrongClassType, ELEBadIndex, ELEDimensionMisMatch, ELELapack
};

// Available classes in in this header file:
class CLinearObject;
class CNumber;
class CVector;
class CMatrix;


// The superclass of all linear objects
class CLinearObject {
public:
  // public members
  double *m_pdData;     // The actual data array
  int m_nDimensions;        // The amount of dimensions
  int *m_pnDimSize;     // The size of a dimension

  // public functions
  CLinearObject();      // Constructor
  CLinearObject(int n);     // Constructor
  CLinearObject(int n, int m);  // Constructor
  CLinearObject(int nDim, int *pnDimSize);  // Constructor
  ~CLinearObject();     // Destructor
  void Randomize();	// Give all elements a random value between 0 and 1

// Functions to initialize a linear object
// Todo: maybe add variable arguments to this function?
  void Initialize(int n);
  void Initialize(int n, int m);

  
// Should be used only by functions of derived classes
  void SetRawData(double *pdData, int nDimensions, int *pnDimSize, bool bAllocated, ELinearClassType eType);
  inline void SetAllocated(bool bAllocated) { if(! m_bDefined){ throw ELENotDefined; } m_bAllocated = bAllocated; }
  inline double *GetRawData() { return m_pdData; }
  int GetElements();	// Return the amount of elements in the object

  void DeleteData();        // Delete all data in the object
  
  inline bool Defined() {return m_bDefined; }
  inline bool Allocated() {return m_bAllocated; }
  inline ELinearClassType GetType() {return m_eType; }

  // public operators
  double &operator[](int i);    // Make index notation possible
  &operator double();       	// Conversion to double
  CLinearObject &operator+=(double d);	// Add a number (element-wise)
  CLinearObject &operator+=(int n) { return this->operator+=(double(n)); }
  CLinearObject &operator-=(double d);	// Subtract a number (element-wise)
  CLinearObject &operator-=(int n) { return this->operator-=(double(n)); }
  CLinearObject &operator*=(double d);	// Multiply with a number (element-wise)
  CLinearObject &operator*=(int n) { return this->operator*=(double(n)); }
  CLinearObject &operator/=(double d);	// Devide by a number (element-wise)
  CLinearObject &operator/=(int n) { return this->operator/=(double(n)); }

  CLinearObject &operator+(double d);	// Add a number (element-wise)
  CLinearObject &operator+(int n) { return this->operator+(double(n)); }
  CLinearObject &operator-(double d);	// Add a number (element-wise)
  CLinearObject &operator-(int n) { return this->operator-(double(n)); }
  CLinearObject &operator*(double d);	// Add a number (element-wise)
  CLinearObject &operator*(int n) { return this->operator*(double(n)); }
  CLinearObject &operator/(double d);	// Add a number (element-wise)
  CLinearObject &operator/(int n) { return this->operator/(double(n)); }

private:
  // private members
protected:
  // protected members
  bool m_bDefined;      // Whether or not the data is defined
  bool m_bAllocated;        // Whether or not the object is allocated for multiplication etc.

  // For assignmetruent purposes:
  int m_indexi;         // Keep track of what has been called for
  int m_indexj;         // by the index operators
  double m_dFlags;      // The flags that might be set

  ELinearClassType m_eType; // The type of the class

  // protected functions
  void InitializeNull();    // Reset all data (only for constructor use)
  inline void Initialize() { this->InitializeNull(); }  // For consistency
  void CopyFromLinearObject(CLinearObject &linob);  // Make deep copies of objects
  void CopySingleMembers(CLinearObject &linob);     // Make deep copies of objects
  void ResetFlags();        // Reset/undo all flag operations
  void SetFlags(double d=0.0);  // Set a flag operation
  double *GetPD();		// Get the pd of the element that is indexed (like a[i][j])
};


// Derived class for numbers
class CNumber: public CLinearObject {
public:
  // public members

  // public functions
  CNumber();            // Constructor
  CNumber(double d);        // Constuctor
  CNumber(CNumber &num);    // Copy-constructor
  ~CNumber();           // Destructor

//  double &operator[](int i);  // Make index notation possible
  double &operator=(double d);  // Assign a double to a CNumber
  double &operator=(int n); // Assign a int to a CNumber
  CNumber &operator=(CNumber &num); // Deep-copy assignment operator
  &operator double();       // Conversion to double
  CNumber &operator*(CNumber &num); // * operator
  CVector &operator*(CVector &vec); // * operator
  CMatrix &operator*(CMatrix &mat); // * operator

// Wrapper functions/operators
  CNumber &operator+=(double d);		// Add a number (element-wise)
  CNumber &operator+=(int n) { return this->operator+=(double(n)); }
  CNumber &operator-=(double d);		// Subtract a number (element-wise)
  CNumber &operator-=(int n) { return this->operator-=(double(n)); }
  CNumber &operator*=(double d);		// Multiply with a number (element-wise)
  CNumber &operator*=(int n) { return this->operator*=(double(n)); }
  CNumber &operator/=(double d);		// Devide by a number (element-wise)
  CNumber &operator/=(int n) { return this->operator/=(double(n)); }

  CNumber &operator+(double d);	// Add a number (element-wise)
  CNumber &operator+(int n) { return this->operator+(double(n)); }
  CNumber &operator-(double d);	// Add a number (element-wise)
  CNumber &operator-(int n) { return this->operator-(double(n)); }
  CNumber &operator*(double d);	// Add a number (element-wise)
  CNumber &operator*(int n) { return this->operator*(double(n)); }
  CNumber &operator/(double d);	// Add a number (element-wise)
  CNumber &operator/(int n) { return this->operator/(double(n)); }
private:
  // private members
protected:
  // protected members
};

// Derived class for vectors
class CVector: public CLinearObject {
public:
  // public members

  // public functions
  CVector();            // Constructor
  CVector(int n);       // Constructor
  CVector(CVector &vec);    // Copy-constructor
  ~CVector();           // Destructor

  CNumber &operator[](int i);   // Make index notation possible
  CVector &operator=(CVector &vec); // Deep-copy assignment operator
  CNumber &operator*(CVector &vec); // * operator
  CVector &operator*(CMatrix &mat); // * operator
  CVector &operator&&(CVector &vec);// && operator (element-wise multiplication)

// Wrapper functions/operators
  CVector &operator+=(double d);		// Add a number (element-wise)
  CVector &operator+=(int n) { return this->operator+=(double(n)); }
  CVector &operator-=(double d);		// Subtract a number (element-wise)
  CVector &operator-=(int n) { return this->operator-=(double(n)); }
  CVector &operator*=(double d);		// Multiply with a number (element-wise)
  CVector &operator*=(int n) { return this->operator*=(double(n)); }
  CVector &operator/=(double d);		// Devide by a number (element-wise)
  CVector &operator/=(int n) { return this->operator/=(double(n)); }

  CVector &operator+(double d);	// Add a number (element-wise)
  CVector &operator+(int n) { return this->operator+(double(n)); }
  CVector &operator-(double d);	// Add a number (element-wise)
  CVector &operator-(int n) { return this->operator-(double(n)); }
  CVector &operator*(double d);	// Add a number (element-wise)
  CVector &operator*(int n) { return this->operator*(double(n)); }
  CVector &operator/(double d);	// Add a number (element-wise)
  CVector &operator/(int n) { return this->operator/(double(n)); }

// These operators must be extended more generally. This is just quick 'n dirty
  CVector &operator+(CVector &vec);	// Add a number (element-wise)
  CVector &operator-(CVector &vec);	// Add a number (element-wise)
private:
  // private members
protected:
  // protected members
};


// Derived class for matrices
class CMatrix: public CLinearObject {
public:
  // public members

  // public functions
  CMatrix();            // Constructor
  CMatrix(int n, int m);    // Constructor
  CMatrix(CMatrix &mat);    // Copy-constructor
  ~CMatrix();           // Destructor

  CVector &operator[](int i);   // Make index notation possible
  CMatrix &operator[](double d);// Make special functions possible
  CMatrix &operator=(CMatrix &mat); // Deep-copy assignment operator
  CMatrix &operator*(CMatrix &mat); // * operator
  CMatrix &operator/(CMatrix &mat); // / operator
  CVector &operator*(CVector &vec); // * operator
  CMatrix &operator&&(CMatrix &mat);// && operator (element-wise multiplication)
  
  CMatrix &operator&(CMatrix &mat); // & operator (tensor product)

  CMatrix &Invert(double *pdLogDet=NULL);        // Invert this matrix using LU-factorization
  CMatrix &InvertSVD();                          // Invert with SVD (POS DEF SyM !!!)
  CMatrix &InvertChol(double *pdLogDet=NULL);    // ''  but using Cholesky
  CMatrix &Inverse(double *pdLogDet=NULL);       // Return the inverse of this matrix (no change to original) using LU
  CMatrix &InverseChol(double *pdLogDet=NULL);   // '' but using Cholesky

  CVector &Eigenvalues();   // Compute the eigenvalues
  CMatrix &Eigenvectors();  // Compute the eigenvectors
  CLinearObject *Eigen(const char *pcJob="N");   // Compute the eigenvectors/values

  CMatrix &Cholesky();		// Compute the Cholesky factorization of real-symm. matrix
  void SVD(CMatrix &mdU, CMatrix &mdV, CVector &vdS);	// The Singular Value Decomposition
  CMatrix &LUFact();		// Compute the LU-factorization of real general matrix

  // A function for debugging the precision of the package:
  CMatrix &CholeskySingle();	// Compute the single precision version of Cholesky
  CMatrix &QRDecompQ();		// Return the Q-factor of QR decomposition

// Wrapper functions/operators
  CMatrix &operator+=(double d);		// Add a number (element-wise)
  CMatrix &operator+=(int n) { return this->operator+=(double(n)); }
  CMatrix &operator-=(double d);		// Subtract a number (element-wise)
  CMatrix &operator-=(int n) { return this->operator-=(double(n)); }
  CMatrix &operator*=(double d);		// Multiply with a number (element-wise)
  CMatrix &operator*=(int n) { return this->operator*=(double(n)); }
  CMatrix &operator/=(double d);		// Devide by a number (element-wise)
  CMatrix &operator/=(int n) { return this->operator/=(double(n)); }

  CMatrix &operator+(double d);	// Add a number (element-wise)
  CMatrix &operator+(int n) { return this->operator+(double(n)); }
  CMatrix &operator-(double d);	// Add a number (element-wise)
  CMatrix &operator-(int n) { return this->operator-(double(n)); }
  CMatrix &operator*(double d);	// Add a number (element-wise)
  CMatrix &operator*(int n) { return this->operator*(double(n)); }
  CMatrix &operator/(double d);	// Add a number (element-wise)
  CMatrix &operator/(int n) { return this->operator/(double(n)); }

// These operators must be extended more generally. This is just quick 'n dirty
  CMatrix &operator+(CMatrix &mat);	// Add a number (element-wise)
  CMatrix &operator-(CMatrix &mat);	// Add a number (element-wise)
private:
  // private members
protected:
  // protected members
};

CMatrix &RandomOrthogonalMatrix(int n);		// Generate a random orthogonal square matrix

// endif __LINAL_H__
#endif

