/* banfunc.cpp -- extension functions for the pta program

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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <termios.h>
#include <unistd.h>

#include "banfunc.h"


/* Some variable that handle the state of the line-image.
 *
 * These handle progress-bars and status reports
 *
 * Phases:
 * 0: Nothing is going on
 * 1: ProgressBar initialised
 * 2: Statusreport initialised
 * */
const int nScreenWidth=80;
const char strRotateImage[] = "-\\|/";

static int nRotateChar=0;
static int nShowPhase=0;
static int nBarStart=0;
static int nProcessRank=0;
static char strImage[160];

/* This function is a linux-version of the DOS kbhit (keyboard-hit) function. It
 * returns an integer representing the input from the keyboard, non-blocking
 *
 * Returnvalue:
 *   -1: No input present
 *  >=0: ASCII value of the input
 *
 *  Last update: 31-10-2008
 * */
int KbHit() {
  struct termios oOldTerm, oNewTerm;
  int nCh;

  tcgetattr(STDIN_FILENO, &oOldTerm);
  oNewTerm = oOldTerm;

  oNewTerm.c_cc[VMIN] = 0;
  oNewTerm.c_cc[VTIME] = 1;
  oNewTerm.c_lflag &= ~(ICANON | ECHO);
  tcsetattr(STDIN_FILENO, TCSANOW, &oNewTerm);

  nCh = getchar();
  tcsetattr(STDIN_FILENO, TCSANOW, &oOldTerm);

  return nCh;
} // kbhit


// Initialize some constants, functions and random numbers. Also set the
// processrank. If it != 0, we will not print output.
bool BANFuncInitialize(int nProcessRank) {
  // Initialize random number generator
  srand ( time(NULL) );

  // Set nProcessRank
  ::nProcessRank = nProcessRank;

  return true;
} // initialize


/* The progressbar look like this:
 * "[==== 67x max ===>   ] xxx %  / "
 * Which is in total 79 bytes long. Universal:
 * "[==== nScreenWidth-13x max ===>   ] xxx %  / "
 *
 * */
void InitProgressBar(const char *strMsg) {
  if(nShowPhase != 0) fprintf(stderr, "\n");

  nShowPhase = 1;
  nRotateChar = 0;

  if(strMsg == NULL) { // We don't have to print any msg
    nBarStart = 0;
  } else { // First print a Msg, then the progressbar
    strcpy(strImage, strMsg);
    nBarStart = strlen(strImage) + 1;
    strImage[nBarStart - 1] = ' ';
  } // if strMsgt

  strImage[nBarStart] = '[';
  for(int i=nBarStart+1; i<nScreenWidth-12; i++)
    strImage[i] = ' ';
  if(! ::nProcessRank) fprintf(stderr, "[?25l"); 
  sprintf(strImage + nScreenWidth-12, "]   0 %%  %c \r", *(strRotateImage+nRotateChar));
  return;
} // InitProgressBar()

void DrawProgressBar(int nPercent) {
  int nMax;

  if(nShowPhase != 1) {
    if(! ::nProcessRank) fprintf(stderr, "\n");
    InitProgressBar();
    nShowPhase = 1;
  } // if nShowPhase

  nRotateChar = (nRotateChar + 1) % 4;

//  nMax = int(double((nPercent*(nScreenWidth-13)))/100.0);
  nMax = nBarStart + 1 + int(double((nPercent*(nScreenWidth-nBarStart-12)))/100.0);
  for(int i=nBarStart + 1; i<nMax; i++) {
    if(i == nMax - 1 && nPercent != 100)
      strImage[i] = '>';
    else
      strImage[i] = '=';
  } // for i
  sprintf(strImage + nScreenWidth-12, "] %3i %%  %c \r", nPercent, *(strRotateImage+nRotateChar));

  if(! ::nProcessRank) fprintf(stderr, strImage);
  return;
} // DrawProgressBar

void FinishProgressBar() {
  int nMax;

  if(nShowPhase != 1) fprintf(stderr, "\n");

  nShowPhase = 0;
  nRotateChar = (nRotateChar + 1) % 4;

//  nMax = nScreenWidth-13;
  nMax = nBarStart + 1 + nScreenWidth-nBarStart-12;
  for(int i=nBarStart + 1; i<nMax; i++) {
    strImage[i] = '=';
  } // for i
  sprintf(strImage + nScreenWidth-12, "] 100 %%    \n");
  if(! ::nProcessRank) fprintf(stderr, strImage);
  if(! ::nProcessRank) fprintf(stderr, "[?25h"); 
  return;
} // FinishProgressBar


void PrintStatus(const char *str) {
  if(nShowPhase != 0) printf("\n");
  nShowPhase = 2;

  strcpy(strImage, str);
  for(int i=strlen(strImage); i<nScreenWidth-10; i++)
    strImage[i] = ' ';

  sprintf(strImage+nScreenWidth-10, "[      ]\r");

  if(! ::nProcessRank) fprintf(stderr, strImage);
  return;
} // PrintStatus

void PrintFailed() {
  if(nShowPhase != 2) {
    PrintStatus(" ");
  } // if nShowPhase
  sprintf(strImage+nScreenWidth-10, "[FAILED]\n");

  if(! ::nProcessRank) fprintf(stderr, strImage);
  nShowPhase = 0;
} // PrintFailed


void PrintSuccess() {
  if(nShowPhase != 2) {
    PrintStatus(" ");
  } // if nShowPhase
  sprintf(strImage+nScreenWidth-10, "[  OK  ]\n");

  if(! ::nProcessRank) fprintf(stderr, strImage);
  nShowPhase = 0;
} // PrintSuccess

void PrintUpdate(const char *str) {
//  if(nShowPhase != 0) printf("\n");
  nShowPhase = 2;

  strcpy(strImage, str);
  for(int i=strlen(strImage); i<nScreenWidth; i++)
    strImage[i] = ' ';

  sprintf(strImage+nScreenWidth-1, "\r");

  if(! ::nProcessRank) fprintf(stderr, strImage);
  return;
} // PrintUpdate


void FinishUpdate() {
  if(nShowPhase != 2) {
    PrintStatus(" ");
  } // if nShowPhase
  sprintf(strImage+nScreenWidth-10, "[  OK  ]\n");

  if(! ::nProcessRank) fprintf(stderr, strImage);
  nShowPhase = 0;
} // FinishUpdate

int n_min(int *pnArray, int nLength) {
  int nMinIndex=0;
  for(int i=0; i<nLength; i++) {
    if(pnArray[i] < pnArray[nMinIndex])
      nMinIndex = i;
  } // for i
  return pnArray[nMinIndex];
} // d_min

double d_min(double *pdArray, int nLength) {
  int nMinIndex=0;
  for(int i=0; i<nLength; i++) {
    if(pdArray[i] < pdArray[nMinIndex])
      nMinIndex = i;
  } // for i
  return pdArray[nMinIndex];
} // d_min

void d_swap(double &x, double &y) {
  double dTemp;
  dTemp = x;
  x = y;
  y = dTemp;
} // d_swap()


/* This function is a 'parallel' quicksort implementation. It sorts using a
 * recursive devide & conquer algorithm, making the total required running time
 * of order n*log(n).
 *
 * Input/Output: pfBufferX, pfBufferY, pfBufferZ
 *   These are the arrays that need to be sorted. Of course they all need to
 *   exist between the left & right indices. On return, one of them
 *   (pfBufferComp) is sorted, they others have changed with it (thus not
 *   necessarily sorted)
 *
 * Input: pfBufferComp
 *   This should be one of the XYZ buffers above. This is the one that is used
 *   as the comparison buffer.
 *
 * Input: left, right
 *   The first & last index of the buffers that need to be sorted.
 *
 * Example:
 *   double x[N], y[N], z[N];
 *   for(i=0; i<N; i++) { x[i] = N-i; y[i] = N*i; z[i] = x[i]*y[i]; }
 *   QuickSort(x, y, z, x, 0, N-1);
 *
 * Author: Rutger van Haasteren <haasteren@strw.leidenuniv.nl>
 * Date: 15-08-2007
 * */
void QuickSort(double *pdBufferX, double *pdBufferY, double *pdBufferZ, double *pdBufferComp, int left, int right) {
//  left is the lower index, right is the upper index
//  of the region of array a that is to be sorted
  int i=left, j=right;
  double pivot = pdBufferComp[(left+right)/2];

  //  partition
  while(i <= j) {    
    while (pdBufferComp[i]<pivot) i++; 
    while (pdBufferComp[j]>pivot) j--;
    if (i<=j) {
      d_swap(pdBufferX[i], pdBufferX[j]);
      d_swap(pdBufferY[i], pdBufferY[j]);
      d_swap(pdBufferZ[i], pdBufferZ[j]);
      i++; j--;
    } // if i<=j
  } // while i<=j

  //  recursion
  if (left<j) QuickSort(pdBufferX, pdBufferY, pdBufferZ, pdBufferComp, left, j);
  if (i<right) QuickSort(pdBufferX, pdBufferY, pdBufferZ, pdBufferComp, i, right);
} // QuickSort()

/* This function is a 'regular' quicksort implementation. It sorts using a
 * recursive devide & conquer algorithm, making the total required running time
 * of order n*log(n).
 *
 * Input/Output: pfBuffer
 *   These are the arrays that need to be sorted.
 *
 * Input: left, right
 *   The first & last index of the buffers that need to be sorted.
 *
 * Example:
 *   double x[N];
 *   for(i=0; i<N; i++) { x[i] = N-i; }
 *   QuickSort(x, 0, N-1);
 *
 * Author: Rutger van Haasteren <haasteren@strw.leidenuniv.nl>
 * Date: 26-01-2009
 * */
void QuickSort(double *pdBuffer, int left, int right) {
//  left is the lower index, right is the upper index
//  of the region of array a that is to be sorted
  int i=left, j=right;
  double pivot = pdBuffer[(left+right)/2];

  //  partition
  while(i <= j) {    
    while (pdBuffer[i]<pivot) i++; 
    while (pdBuffer[j]>pivot) j--;
    if (i<=j) {
      d_swap(pdBuffer[i], pdBuffer[j]);
      i++; j--;
    } // if i<=j
  } // while i<=j

  //  recursion
  if (left<j) QuickSort(pdBuffer, left, j);
  if (i<right) QuickSort(pdBuffer, i, right);
} // QuickSort()

void QuickSort(double **ppdBuffer, double *pdBufferZ, int nDimensions, int left, int right) {
//void Evidence::QuickSort(double *pdBufferX, double *pdBufferY, double *pdBufferZ, int *pnBufferI, int *pnBufferComp, int left, int right) {
//  left is the lower index, right is the upper index
//  of the region of array a that is to be sorted
  int i=left, j=right, d;
  double pivot = pdBufferZ[(left+right)/2];

  //  partition
  while(i <= j) {    
    while (double(pdBufferZ[i])<pivot) i++; 
    while (double(pdBufferZ[j])>pivot) j--;
    if (i<=j) {
      d_swap(pdBufferZ[i], pdBufferZ[j]);
      for(d=0; d<nDimensions; d++)
	d_swap(ppdBuffer[i][d], ppdBuffer[j][d]);
      i++; j--;
    } // if i<=j
  } // while i<=j

  //  recursion
  if (left<j) QuickSort(ppdBuffer, pdBufferZ, nDimensions, left, j);
  if (i<right) QuickSort(ppdBuffer, pdBufferZ, nDimensions, i, right);
} // QuickSort()
