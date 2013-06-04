/* moremath.h -- some mathematical functions that do not exist in the
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


#ifndef __MOREMATH_H__
#define __MOREMATH_H__

double inverf(double d);                      // The inverse error function
int fsign(double d);                          // Return the sign of the parameter
double gammln(double d);			// The log of the gamma function
double gamm(double d);

#endif // __PTAFUNC_H__
