/* correlation.cpp -- Whitening cpp file for PTA calculations program
 *
 * Rutger van Haasteren 16 April 2008 haasteren@strw.leidenuniv.nl
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

#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include "correlation.h"
#include "linal.h"
#include "linalfunc.h"
#include "banfunc.h"
#include "corefunctions.h"


/* This function calculates the Maximum likelihood value for the correlation,
 * given the spectral index */
double MLCorrelation(CVector &vdSatA, CVector &vdResA, CVector &vdSatB, CVector &vdResB, double dGamma) {
  SCorrelation oCor;
  double dML;

  oCor.vdAX = vdSatA;
  oCor.vdAY = vdResA;
  oCor.vdBX = vdSatB;
  oCor.vdBY = vdResB;

#if 0
  double dMean;
    dML = 0; dMean = 0;
    for(int i=0; i<vdResA.m_pnDimSize[0]; i++) {
      dML += double(vdResA[i]) * double(vdResB[i]);
      dMean += gsl_sf_pow_int(double(vdResA[i]), 2) * 
	gsl_sf_pow_int(double(vdResB[i]), 2);
    } // for i
    dML = dML / sqrt(vdResA.m_pnDimSize[0]*dMean);
#endif

  Brent1DOptimum(oCor, -0.9, 0, 0.9, 0.000001, dML);
  return dML;
} // MLCorrelation


/* We will see what this function does next... */
double PowerLaw(double dT, double dA, double dExp, double dFc) {
  double dTij, dReturnValue=0;
  double dG, dFl, dN, dSum, dPrevSum, dAp;
  double dGammaSine;
  unsigned int n=0;

//	dAp = dA * sqrt(1 / 3.0) * 1E-15;
	dAp = dA * sqrt(SPERYEAR / 3.0) * 1E-15;
	dG = dExp;
//	dFl = dFc;
	dFl = dFc / SPERYEAR;
	dTij = 2*M_PI * dT; /// (2 * 365);

//	dN = gsl_pow_int(dAp,2)*pow(2*NSPERYEAR/dFl, 1+dG)
//	  /gsl_pow_int(2*M_PI,2);
//	dN = gsl_pow_int(NSPERYEAR,2)*
//	  gsl_pow_int(dAp,2)*pow(2/dFl, 1+dG)/gsl_pow_int(2*M_PI,2);
//	dN = gsl_pow_int(dAp,2)*pow(SPERYEAR/dFl, 1+dG)/
//	  (SPERYEAR*gsl_pow_int(2*M_PI,2));
	dN = gsl_pow_int(dAp,2)*SPERYEAR/
	  (pow(dFl*SPERYEAR, 1+dG)*gsl_pow_int(2*M_PI,2));

	// First the Fl dependent part: the Generalized Hypergeometrical Function
	dSum = -1.0/(1+dG);
	dPrevSum = 0;
	while(dPrevSum != dSum) {
	  n++;
	  dPrevSum = dSum;
	  dSum += gsl_pow_int(-1, n)*gsl_pow_int(dFl*dTij, 2*n)/(gsl_sf_fact(2*n)*(2*n-1-dG));

//	  printf("Iteration %i\n", n);
	  if(n > 100) {
	    printf("Too many iterations\n");
	    break;
	  }
	} // while

	dGammaSine = gsl_sf_gamma(-1-dG)*sin(-M_PI_2*dG);
	dReturnValue = dN*(dGammaSine*pow(dFl*dTij, dG+1) - dSum);
  return dReturnValue;
} // PowerLaw

/* This function calculates f(oCor, zeta, fl=0.005), given vector x, y and stuff
 *
 * Now only does white noise */
double CorrelationLikelihood(SCorrelation &oCor, double dZeta, double dFl) {
  CMatrix mdC, mdG, mdM, mdTemp;
  CVector vdX, vdY, vdQSD;
  int nIndex, nSize, nSize1;
  double dTimeCor, dSpaceCor, dLogDetC, dLogDetMCM;
  double dReturnValue, dRMS;
  static bool bDrawn=false;

  /* Perform sanity checks */
  if(! oCor.vdAX.Defined() || ! oCor.vdAY.Defined()) { throw ELENotDefined; }
  if(! oCor.vdBX.Defined() || ! oCor.vdBY.Defined()) { throw ELENotDefined; }
  if(oCor.vdAX.m_pnDimSize[0] != oCor.vdAY.m_pnDimSize[0]) {
    throw ELEDimensionMisMatch;}
  if(oCor.vdBX.m_pnDimSize[0] != oCor.vdBY.m_pnDimSize[0]) {
    throw ELEDimensionMisMatch;}

  /* Allocate memory */
  nSize = oCor.vdAX.m_pnDimSize[0] + oCor.vdBX.m_pnDimSize[0];
  nSize1 = oCor.vdAX.m_pnDimSize[0];
  mdC.Initialize(nSize, nSize);
  mdG.Initialize(2, 2);
  mdM.Initialize(nSize, 6);
  vdX.Initialize(nSize);
  vdY.Initialize(nSize);

  /* Fill the correlation vectors */
  nIndex = 0;
//  dRMS = sqrt(oCor.vdAY * oCor.vdAY);
  dRMS = 1;
  for(int i=0; i<oCor.vdAY.m_pnDimSize[0]; i++) {
    vdX[nIndex] = double(oCor.vdAX[i]);
    vdY[nIndex] = double(oCor.vdAY[i]) / dRMS;
    nIndex++;
  } // for i
//  dRMS = sqrt(oCor.vdBY * oCor.vdBY);
  dRMS = 1;
  for(int i=0; i<oCor.vdBY.m_pnDimSize[0]; i++) {
    vdX[nIndex] = double(oCor.vdBX[i]);
    vdY[nIndex] = double(oCor.vdBY[i]) / dRMS;
    nIndex++;
  } // for i

  /* Fill the correlation matrices */
  mdG[0][0] = 1; mdG[1][1] = 1;
  mdG[0][1] = dZeta; mdG[1][0] = dZeta;
  for(int i=0; i<nSize; i++) {
    for(int j=0; j<nSize; j++) {
      /* White noise assumption in dTimeCor */
//      dTimeCor = double(vdX[i]) == double(vdX[j]) ? 1 : 0;
//      dTimeCor = (fabs(double(vdX[i]) - double(vdX[j])) < 1.7 ? 1 : 0);
//      dTimeCor = PowerLaw(fabs(double(vdX[i])-double(vdX[j]))/(2*365), 1, 2.330001, 0.05); //
      dTimeCor = PowerLaw(fabs(double(vdX[i])-double(vdX[j])), 1, 2.330001, 0.2); //
      dSpaceCor = double(mdG[i < nSize1 ? 0 : 1][j < nSize1 ? 0 : 1]);
      mdC[i][j] = dTimeCor * dSpaceCor;

      /* add uncorrelated white noise */
//      if( i==j )
//	mdC[i][j] += 25;
    } // for j
  } // for i

  /* Fill the quadratic spindown matrix */
  for(int i=0; i<nSize; i++) {
    for(int a=0; a<2; a++) {
      for(int j=0; j<3; j++) {
	mdM[i][j + a*3] = 0;
	if(a == (i<nSize1 ? 0 : 1))
	  mdM[i][j + a*3] += gsl_sf_pow_int(double(vdX[i][j+a*3])/2*365, j);
      } // for j
    } // for a
  } // for i

//  printf("mdC.InvertChol(&dLogDetC);\n");
  /* Calculate the intermediate matrices */
//  mdC.InvertChol(&dLogDetC);
  mdC.InvertChol(&dLogDetC);
//  mdTemp = mdM[LO_TRANSPOSE] * mdC * mdM;
//  printf("mdTemp.InvertChol(&dLogDetMCM);\n");
//  mdTemp.InvertChol(&dLogDetMCM);
//  vdY *= 1E9;

  double dTemp=0;
  dReturnValue = 0;
  dReturnValue = 0.5 * dLogDetC;
  dTemp += (vdY * mdC * vdY);
//  dTemp -= vdY * (mdC * (mdM * (mdTemp * (mdM[LO_TRANSPOSE] * (mdC * vdY)))));
  dReturnValue += ((vdY.m_pnDimSize[0] - 1)/2.0) * gsl_sf_log(dTemp);

  /* Now without integrating, with white noise... */
//  dReturnValue = 0.5 * dLogDetC;
//  printf("dLogDetC: %.15f\n", dLogDetC);
//  dReturnValue = 0;
//  dReturnValue += 0.5 * (vdY * mdC * vdY);
//  printf("vdY*mdC*vdY: %.25f\n", 2*dReturnValue);
//  dTemp -= vdY * (mdC * (mdM * (mdTemp * (mdM[LO_TRANSPOSE] * (mdC * vdY)))));
//  dReturnValue += ((vdY.m_pnDimSize[0] - 1)/2.0) * gsl_sf_log(dTemp);

  return dReturnValue;
} // CorrelationLikelihood


/* This function speaks for itself: shift 3 values
 * */
inline void shft3(double &a, double &b, double &c, const double d) {
  a=b;
  b=c;
  c=d;
}

inline double SIGN(const double &a, const float &b)
        {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}


/* 1 Dimensional optimization algorithm - Brent
 * Takes as input a direction (vector) and a starting position in a parameter
 * space, and then minimizes a given function in the given direction.
 *
 * Quote NR:
 * Given a function f, and given a bracketing triplet of abscissas ax, bx, cx
 * (such that bx is between ax and cx, and f(bx) is less than both f(ax) and
 * f(cx)), this routine isolates the minimum to a fractional precision of about
 * tol using Brents's method. The abscissa of the minimum is returned as xmin,
 * and the minimum function value is returned as brent, the returned function
 * value */
double Brent1DOptimum(SCorrelation &oCor, const double ax, const double bx,
    const double cx, const double tol, double &xmin) {
  // Here, ITMAX is the maximum allowed number of iterations, CGOLD is the
  // golden ratio; ZEPS is a small number that portects against trying to
  // achieve fractional accuracy for a minimum that happens to be exactly zero.
  const int ITMAX=100;
  const double CGOLD=0.3819660;
  const double ZEPS = 1.0e-10;
  int iter;
  double a,b,d=0.0,etemp,fu,fv,fw,fx;
  double p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=CorrelationLikelihood(oCor, x);
  for (iter=0;iter<ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=CorrelationLikelihood(oCor, u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      shft3(v,w,x,u);
      shft3(fv,fw,fx,fu);
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    }
  }
//  nrerror("Too many iterations in brent");
  xmin=x;
  return fx;
} // Brent1DOptimum
