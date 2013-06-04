/* moremath.cpp -- some mathematical functions that do not exist in the
                 math liberary (the are also available in the GSL)

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

#include "moremath.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

// De inverse error function
double inverf(double d) {
  double tol, z_old, z_new, s;
  int iterations=0, maxit;

  tol = 3.5 * 2.2204E-16;	// Machine precision - tolerance
  maxit = 100;

  if(d == -1) return GSL_NEGINF;
  if(d == 1) return GSL_POSINF;
  if(d < -1 || d > 1) return GSL_NAN;

  s = M_SQRTPI/2;
  z_old = 1;
  z_new = sqrt(-log(1-fabs(d))) * GSL_SIGN(d);
  while(fabs(gsl_sf_erf(z_new) - d) > tol * fabs(d)) {
    z_old = z_new;
    z_new = z_old - (gsl_sf_erf(z_old) - d) * exp(gsl_pow_2(z_old)) * s;
    if(iterations++ > maxit) {
      printf("iterations: %i  maxit: %i \n", iterations, maxit);
      GSL_ERROR("erfinv: iteration limit exceeded", GSL_ETOL);
      break;
    } // if iterations
  } // while
  return z_new;
} // inverf


// The log of the gamma function (numerical recipes -> corrected)
double gammln(double d) {
  double tmp, ser;
  tmp  = d + 4.5;
  tmp -= (d - 0.5) * log(tmp);

  ser = 1.000000000190015
      + (76.18009172947146     / d)
      - (86.50532032941677     / (d + 1.0))
      + (24.01409824083091     / (d + 2.0))
      - (1.231739572450155     / (d + 3.0))
      + (0.1208650973866179e-2 / (d + 4.0))
      - (0.5395239384953e-5    / (d + 5.0));
  return (log(2.5066282746310005 * ser) - tmp);
} // gammln

// The practical gamma function
double gamm(double d) {
  double dReturnValue;

  if(d > 0)
    dReturnValue = exp(gammln(d));
  else
    dReturnValue = -M_PI/(d*exp(gammln(-d))*sin(M_PI*d));

  return dReturnValue;
} // gamma


// Return the sign of the parameter
int fsign(double d) {
  return int(d/fabs(d));
} // fsign

