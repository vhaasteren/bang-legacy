/* iugui.cpp -- Interactive universal graphical user interface
 *
 * Rutger van Haasteren 22 Novermber 2008 haasteren@strw.leidenuniv.nl
 *
 * Copyright (C) 2008-2009 Rutger van Haasteren.
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

 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.  */

// TODO: Implement the NeedUpdate stuff
// TODO: Use malloc instead of new
// TODO: It's possible to add a single widget several times to a list / window.
//       This could cause bugs.

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "iugui.h"

#ifdef HAVE_PGPLOT
#include <cpgplot.h>
#endif // HAVE_PGPLOT


CIWidgetList::CIWidgetList() {
  m_nSize = 0;
  m_peFirst = NULL;
  m_peLast = NULL;
} // CIWidgetList

CIWidgetList::~CIWidgetList() {
  // Remove all the list elements
  CIWidgetListElement *peCurrent, *peNext;
  peCurrent = m_peFirst;

  while(peCurrent != NULL) {
    peNext = peCurrent->GetNext();
    delete peCurrent;
    m_nSize--;
    peCurrent = peNext;
    peNext = NULL;
  } // while peCurrent

  if(m_nSize != 0) printf("WARNING: Serious error in ~CIWidgetList()\n");
  m_nSize = 0;
} // ~CIWidgetList

void CIWidgetList::AddWidget(CIWidget *pwWidget) {
  CIWidgetListElement *peNewElement;
  peNewElement = new CIWidgetListElement(pwWidget);

  // Set the new couplings
  if(m_peLast != NULL) {
    // The last element exists
    m_peLast->SetNext(peNewElement);
    peNewElement->SetPrevious(m_peLast);
  } else {
    // This is the first item in the list
    m_peFirst = peNewElement;
    peNewElement->SetPrevious(NULL);
  } // if m_peLast

  // We have a new last element
  peNewElement->SetNext(NULL);
  m_peLast = peNewElement;
  m_nSize++;
} // AddWidget

void CIWidgetList::RemoveWidget(CIWidget *pwWidget) {
  CIWidgetListElement *peCurrent;
  peCurrent = GetListElement(pwWidget);
  if(peCurrent == NULL) {
    printf("WARNING: List is empty\n");
    return ;
  } // if m_peFirst

  // Found the widget, now remove it from the list
  if(peCurrent->GetNext() != NULL && peCurrent->GetPrevious() != NULL) {
    peCurrent->GetNext()->SetPrevious(peCurrent->GetPrevious());
    peCurrent->GetPrevious()->SetNext(peCurrent->GetNext());
  } else if(peCurrent->GetNext()) {
    peCurrent->GetNext()->SetPrevious(NULL);
    m_peFirst = peCurrent->GetNext();
  } else if(peCurrent->GetPrevious()) {
    peCurrent->GetPrevious()->SetNext(NULL);
    m_peLast = peCurrent->GetPrevious();
  } else {
    m_peFirst = NULL;
    m_peLast = NULL;
  } // if GetNext

  delete peCurrent;
  m_nSize--;
  return ;
} // RemoveWidget

int CIWidgetList::GetSize() {
  return m_nSize;
} // GetSize

bool CIWidgetList::HasWidget(CIWidget *pwWidget) {
  CIWidget *pwCurrent = NULL;

  if(! m_peFirst)
    return false;

  pwCurrent = m_peFirst->GetWidget();
  while(pwCurrent) {
    if(pwWidget == pwCurrent)
      return true;
    pwCurrent = GetNext(pwCurrent);
  } // while pwCurrent
  return false;
} // HasWidget

CIWidget *CIWidgetList::GetFirst() {
  if(m_peFirst)
    return m_peFirst->GetWidget();
  else
    return NULL;
} // GetFirst

CIWidget *CIWidgetList::GetLast() {
  if(m_peLast)
    return m_peLast->GetWidget();
  else
    return NULL;
} // GetLast


CIWidget *CIWidgetList::GetNext(CIWidget *pwWidget) {
  CIWidgetListElement *peCurrent;
  peCurrent = GetListElement(pwWidget);
  if(peCurrent == NULL) {
    printf("WARNING: List is empty\n");
    return NULL;
  } // if m_peFirst

  if(! peCurrent->GetNext())
    return NULL;
  else
    return  peCurrent->GetNext()->GetWidget();
} // GetNext

CIWidget *CIWidgetList::GetPrevious(CIWidget *pwWidget) {
  CIWidgetListElement *peCurrent;
  peCurrent = GetListElement(pwWidget);
  if(peCurrent == NULL) {
    printf("WARNING: List is empty\n");
    return NULL;
  } // if m_peFirst

  if(! peCurrent->GetPrevious() )
    return NULL;
  else
    return  peCurrent->GetPrevious()->GetWidget();
} // GetPrevious


CIWidgetListElement *CIWidgetList::GetListElement(CIWidget *pwWidget) {
  CIWidgetListElement *peCurrent = m_peFirst;
  if(m_peFirst == NULL) {
    printf("WARNING: List is empty\n");
    return NULL;
  } // if m_peFirst

  while(peCurrent->GetWidget() != pwWidget) {
    peCurrent = peCurrent->GetNext();
    if(peCurrent == NULL) {
      printf("WARNING: List element not found\n");
      return NULL;
    } // if m_peCurrent
  } // while GetWidget

  return peCurrent;
} // GetListElement


CIWidgetListElement::CIWidgetListElement(CIWidget *pwWidget) {
  m_pwWidget = pwWidget;
  m_peNext = NULL;
  m_pePrevious = NULL;
} // CIWidgetListElement

CIWidgetListElement::~CIWidgetListElement() {
  m_pwWidget = NULL;
  m_peNext = NULL;
  m_pePrevious = NULL;
} // ~CIWidgetListElement

CIWidget *CIWidgetListElement::GetWidget() {
  return m_pwWidget;
} // GetWidget

CIWidgetListElement *CIWidgetListElement::GetNext() {
  return m_peNext;
} // GetNext

CIWidgetListElement *CIWidgetListElement::GetPrevious() {
  return m_pePrevious;
} // GetPrevious

CIWidgetListElement *CIWidgetListElement::SetNext(CIWidgetListElement *peNext) {
  m_peNext = peNext;
  return peNext;
} // SetNext

CIWidgetListElement *CIWidgetListElement::SetPrevious(CIWidgetListElement *pePrevious) {
  m_pePrevious = pePrevious;
  return pePrevious;
} // SetPrevious


CIWindow::CIWindow() {
  // In the case of MPI we cannot call this anymore. TODO: Demand from the user
  // to initialise after creating an object
  m_bInitialised = false;
  m_bAlreadyPapped = false;
//  Initialise("/xs"); // Device /xs = standard output device
} // CIWindow

CIWindow::~CIWindow() {
  if(m_bInitialised)
    DeInitialise();
} // ~CIWindow

void CIWindow::DeInitialise() {
#ifdef HAVE_PGPLOT
  cpgend();
#endif // HAVE_PGPLOT
  m_bInitialised = false;
} // DeInitialise

void CIWindow::Initialise(const char strDevice[]) {
//  const char strDevice[] = "/xs";	// Standard output device
  const char strBkgrdColour[] = "black";
  const char strLineColour[] = "white";
  const char strButtonColour[] = "green";
  const char strButtonTextColour[] = "black";
  const char strTextColour[] = "white";
  const char strActiveColour[] = "red";
  const char strSigma1Colour[] = "green";
  const char strSigma2Colour[] = "orange";
  const char strSigma3Colour[] = "red";

  int nStatus;				// A status return variable
  m_fAspectRatio = 0.75;		// Aspect ratio. 1 = square
  m_fWidth = 0.0;			// Width of view in inches. 0 = big
  m_fFontSize = 1.0;			// Font size (vertical). 1 = normal
  m_nFontType = 1;			// Normal font (1 - 4)
  m_nLineWidth = 1;			// Normal line width (1 - 201)

  // Save the device we are working with
  strcpy(m_strDeviceFile, strDevice);

  // Begin PGPLOT window (only 1 subdevision)
#ifdef HAVE_PGPLOT
  cpgbeg(0, strDevice, 1, 1);
#endif // HAVE_PGPLOT

  // Set the size of the view surface
  if(! m_bAlreadyPapped) {
#ifdef HAVE_PGPLOT
    cpgpap(m_fWidth, m_fAspectRatio);
#endif // HAVE_PGPLOT
    m_bAlreadyPapped = true;
  } // if m_bAlreadyPapped

  // Set the viewport (normalised device coordinates)
#ifdef HAVE_PGPLOT
  cpgsvp(0, 1, 0, 1);

  // Set attributes
  cpgscf(m_nFontType);
  cpgslw(m_nLineWidth);

  // Set colour attributes
  cpgscrn(IUGUI_COLOR_BKGRD, strBkgrdColour, &nStatus);
  cpgscrn(IUGUI_COLOR_LINE, strLineColour, &nStatus);
  cpgscrn(IUGUI_COLOR_BUTTON, strButtonColour, &nStatus);
  cpgscrn(IUGUI_COLOR_BUTTONTEXT, strButtonTextColour, &nStatus);
  cpgscrn(IUGUI_COLOR_ACTIVE, strActiveColour, &nStatus);
  cpgscrn(IUGUI_COLOR_TEXT, strTextColour, &nStatus);
  cpgscrn(IUGUI_COLOR_SIGMA1, strSigma1Colour, &nStatus);
  cpgscrn(IUGUI_COLOR_SIGMA2, strSigma2Colour, &nStatus);
  cpgscrn(IUGUI_COLOR_SIGMA3, strSigma3Colour, &nStatus);

  // Set pgplot to non-interactie (we do the interaction ourselves)
  cpgask(0);
#endif // HAVE_PGPLOT

  m_dSizeX = 1;
  m_dSizeY = 1;
  m_bNeedsUpdate = true;
  m_bDone = false;
  m_bInitialised = true;
} // Initialise

/* This function processes signals from the user. Basically all signals present
 * in this design are keyboard strokes and mouse clicks. This function handles
 * those and, if necessary, passes the signals to the children widgets
 * */
void CIWindow::ProcessSignals() {
  float fMouseX, fMouseY;
  float fXMinNor, fXMaxNor, fYMinNor, fYMaxNor;
  float fXMinWorld, fXMaxWorld, fYMinWorld, fYMaxWorld;
  float fXMinFullNor, fXMaxFullNor, fYMinFullNor, fYMaxFullNor;
  char cKey;
  CIWidget *pwTarget = NULL;

  // Request the coordinates (normalised device) of the viewport
#ifdef HAVE_PGPLOT
//  cpgqvp(0, &fXMinNor, &fXMaxNor, &fYMinNor, &fYMaxNor);
  // Set the viewport (normalised device coordinates)
  cpgsvp(0, 1, 0, 1);

  // Request the window coordinates
//  cpgqwin(&fXMinWorld, &fXMaxWorld, &fYMinWorld, &fYMaxWorld);

  // Request total window coordinates (normalised wrt full device)
  cpgqvsz(0, &fXMinFullNor, &fXMaxFullNor, &fYMinFullNor, &fYMaxFullNor);

  // Find out whether there has been clicked/typed
  // Left mouse: 65, right mouse: 88
//  cpgcurs(&fMouseX,&fMouseY,&nKey);
  cpgband(0, 0, 0, 0, &fMouseX, &fMouseY, &cKey);
#endif // HAVE_PGPLOT

//  printf("cKey = %i\n", (int)cKey);
  switch(cKey) {
    case (char)65: // Left mouse button
    case (char)88: // Right mouse button
      // First check for an active widget
      pwTarget = GetWidget(fMouseX, fMouseY);
      if(pwTarget) {
  	pwTarget->OnClicked(fMouseX, fMouseY);
      }
      else
	OnClicked(fMouseX, fMouseY);
      break;
    case 'q':
    case 'Q':
      m_bDone = true;
      break;
    case 's':
    case 'S': // Make a screenshot
      ScreenShot("iugui.ps/ps");
      break;
    default:
      OnKey(cKey);
      break;
  } // switch nKey
} // CIWindow

void CIWindow::ScreenShot(const char strDevice[]) {
  if(! m_bInitialised) return;
#ifdef HAVE_PGPLOT
  cpgend();
#endif // HAVE_PGPLOT
  Initialise(strDevice);
  Update();
#ifdef HAVE_PGPLOT
  cpgend();
#endif // HAVE_PGPLOT

  Initialise("/xs");
  Update();
  return;
} // ScreenShot

/* This function copies the device string we are working with to the argument
 * */
void CIWindow::GetDevice(char *strDevice) {
  strcpy(strDevice, m_strDeviceFile);
} // GetDevice

float CIWindow::GetFontSize() {
  return m_fFontSize;
} // GetFontSize

bool CIWindow::NeedsUpdate() {
  return m_bNeedsUpdate;
} // NeedsUpdate

/* This function draws all the widgets on the screen. It simply loops over all
 * the widgets.
 * */
void CIWindow::Update() {
  CIWidget *pwCurrent = NULL;

#ifdef HAVE_PGPLOT
  cpgeras();
#endif // HAVE_PGPLOT
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    if(pwCurrent->Visible()) {
      pwCurrent->Update();
      pwCurrent->UpdateChildren();
    } // if Visible
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  pwCurrent = GetActiveWidget();
  if(pwCurrent) {
    if(pwCurrent->Visible()) {
      pwCurrent->Update();
      pwCurrent->UpdateChildren();
    } // if Visible
  } // if pwCurrent

//  m_bNeedsUpdate = false;
  return ;
} // Update

void CIWindow::AddWidget(CIWidget *pwWidget) {
  CIWidget *pwCompetition=NULL;
  pwCompetition = GetWidgetInRange(pwWidget->GetPosX(), pwWidget->GetPosY(), pwWidget->GetSizeX(), pwWidget->GetSizeY());
  if(pwCompetition) {
    printf("WARNING: Widget cannot be placed there: competing widget\n");
    return;
  } // if pwCompetition

  m_wlChildren.AddWidget(pwWidget);
  pwWidget->SetParent(this);
} // AddWidget

CIWidget *CIWindow::GetActiveWidget() {
  CIWidget *pwCurrent = m_wlChildren.GetFirst();
  CIWidget *pwActive = NULL,
	   *pwActiveChild = NULL;

  while(pwCurrent) {
    if(pwCurrent->Active()) {
      pwActive = pwCurrent;
      break;
    } // if Active
    pwActiveChild = pwCurrent->GetActiveWidget();
    if(pwActiveChild) {
      pwActive = pwActiveChild;
      break;
    } // if pwActiveChild
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  return pwActive;
} // GetActiveWidget

CIWidget *CIWindow::GetWidget(double dPosX, double dPosY) {
  CIWidget *pwCurrent = NULL;
  CIWidget *pwTarget = NULL;

  // First check for the active widget, that one is always on top
  pwCurrent = GetActiveWidget();
  if(pwCurrent) {
    if(pwCurrent->GetActivePosX() < dPosX &&
	pwCurrent->GetActivePosX() + pwCurrent->GetActiveSizeX() > dPosX &&
	pwCurrent->GetActivePosY() < dPosY &&
	pwCurrent->GetActivePosY() + pwCurrent->GetActiveSizeY() > dPosY &&
	pwCurrent->Visible()) {
      return pwCurrent;
    } // if Active
  } // if pwCurrent

  // It wasn't the active widget. Check all the others
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    if(pwCurrent->GetPosX() < dPosX &&
	pwCurrent->GetPosX() + pwCurrent->GetSizeX() > dPosX &&
	pwCurrent->GetPosY() < dPosY &&
	pwCurrent->GetPosY() + pwCurrent->GetSizeY() > dPosY &&
	pwCurrent->Visible()) {
      pwTarget = pwCurrent;
      break;
    } // if Active
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  return pwTarget;
} // GetWidget

CIWidget *CIWindow::GetNonActiveWidget(double dPosX, double dPosY) {
  CIWidget *pwCurrent = NULL;
  CIWidget *pwTarget = NULL;

  // We don't check for the active widget
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    if(pwCurrent->GetPosX() < dPosX &&
	pwCurrent->GetPosX() + pwCurrent->GetSizeX() > dPosX &&
	pwCurrent->GetPosY() < dPosY &&
	pwCurrent->GetPosY() + pwCurrent->GetSizeY() > dPosY &&
	pwCurrent->Visible()) {
      pwTarget = pwCurrent;
      break;
    } // if Active
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  return pwTarget;
} // GetNonActiveWidget

CIWidget *CIWindow::GetWidgetInRange(double dPosX, double dPosY, double dSizeX, double dSizeY) {
  CIWidget *pwCurrent = NULL;
  CIWidget *pwTarget = NULL;

  // Check all the widgets
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    if(pwCurrent->InRange(dPosX, dPosY) ||
	pwCurrent->InRange(dPosX + dSizeX, dPosY) ||
	pwCurrent->InRange(dPosX, dPosY + dSizeY) ||
	pwCurrent->InRange(dPosX + dSizeX, dPosY + dSizeY) &&
	pwCurrent->Visible()) {
      pwTarget = pwCurrent;
      break;
    } // if InRange
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  return pwTarget;
} // GetWidgetInRange

bool CIWindow::HasWidget(CIWidget *pwWidget) {
  return m_wlChildren.HasWidget(pwWidget);
} // HasWidget

void CIWindow::RemoveWidget(CIWidget *pwWidget) {
  m_wlChildren.RemoveWidget(pwWidget);
} // RemoveWidget

void CIWindow::SetChildrenNotActive() {
  CIWidget *pwCurrent = NULL;

  // Check all the widgets
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    pwCurrent->SetActive(false);
    pwCurrent->SetChildrenNotActive();
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent
  return ;
} // SetChildrenNotActive

/* This is a virtual function, supposed to be overloaded by the derived class.
 * It responds to mouse clicks
 * */
void CIWindow::OnClicked(double dPosX, double dPosY) {
  SetChildrenNotActive();
} // OnClicked

/* This is a virtual function, supposed to be overloaded by the derived class.
 * It responds to keyboard keys
 * */
void CIWindow::OnKey(char cKey) {
} // OnKey

bool CIWindow::Done() {
  return m_bDone;
} // Done

void CIWindow::SetDone(bool bDone) {
  m_bDone = bDone;
} // SetDone

CIWidget::CIWidget(double dSizeX, double dSizeY) {
  Initialise();

  if(dSizeX == -1)
    m_dSizeX = 1;
  else
    m_dSizeX = dSizeX;
  if(dSizeY == -1)
    m_dSizeY = 1;
  else
    m_dSizeY = dSizeY;
} // CIWidget

CIWidget::~CIWidget() {
  return ;
} // ~CIWidget

void CIWidget::Initialise() {
  m_dSizeX = 1;
  m_dSizeY = 1;
  m_dPosX = 0;
  m_dPosY = 0;
  m_dActivePosX = 0;
  m_dActivePosY = 0;
  m_dActiveSizeX = 0;
  m_dActiveSizeY = 0;

  m_bVisible = true;
  m_bActive = false;
  m_bRespectAbsolutePosition = true;

  m_pwParentWindow = NULL;
  m_pwParentWidget = NULL;
} // Initialise

/* This function should only be called when writing directly to a device, not
 * when updating a widget contained in a window. Example: writing a plot to a
 * file. This function initialises the writing session.
 * */
void CIWidget::DeviceInit(const char strDevice[]) {
  const char strBkgrdColour[] = "white";
  const char strLineColour[] = "black";
  const char strButtonColour[] = "green";
  const char strButtonTextColour[] = "black";
  const char strTextColour[] = "black";
  const char strActiveColour[] = "red";
  const char strSigma1Colour[] = "green";
  const char strSigma2Colour[] = "orange";
  const char strSigma3Colour[] = "red";
  int nStatus;

  // Begin PGPLOT window (only 1 subdevision)
#ifdef HAVE_PGPLOT
  cpgbeg(0, strDevice, 1, 1);

  // Set the viewport (normalised device coordinates)
  cpgsvp(0, 1, 0, 1);

  // Set colour attributes
  cpgscrn(IUGUI_COLOR_BKGRD, strBkgrdColour, &nStatus);
  cpgscrn(IUGUI_COLOR_LINE, strLineColour, &nStatus);
  cpgscrn(IUGUI_COLOR_BUTTON, strButtonColour, &nStatus);
  cpgscrn(IUGUI_COLOR_BUTTONTEXT, strButtonTextColour, &nStatus);
  cpgscrn(IUGUI_COLOR_ACTIVE, strActiveColour, &nStatus);
  cpgscrn(IUGUI_COLOR_TEXT, strTextColour, &nStatus);
  cpgscrn(IUGUI_COLOR_SIGMA1, strSigma1Colour, &nStatus);
  cpgscrn(IUGUI_COLOR_SIGMA2, strSigma2Colour, &nStatus);
  cpgscrn(IUGUI_COLOR_SIGMA3, strSigma3Colour, &nStatus);

  // Set pgplot to non-interactie (we do the interaction ourselves)
  cpgask(0);
#endif // HAVE_PGPLOT

  m_bRespectAbsolutePosition = false;
} // DeviceInit

/* This function should only be called when writing directly to a device, not
 * when updating a widget contained in a window. Example: writing a plot to a
 * file. This function ends the writing session.
 * */
void CIWidget::DeviceDeInit() {
#ifdef HAVE_PGPLOT
  cpgend();
#endif // HAVE_PGPLOT
  m_bRespectAbsolutePosition = true;
} // DeviceDeInit

/* This is a virtual function, supposed to be overloaded by the derived class.
 * It responds to mouse clicks
 * */
void CIWidget::OnClicked(double dPosX, double dPosY) {
  if(m_pwParentWindow);
    m_pwParentWindow->SetChildrenNotActive();
} // OnClicked

/* This is a virtual function, supposed to be overloaded by the derived class.
 * It responds to keyboard keys
 * */
void CIWidget::OnKey(char cKey) {
} // OnKey

double CIWidget::GetSizeX() {
  return m_dSizeX;
} // GetSizeX

double CIWidget::GetSizeY() {
  return m_dSizeY;
} // GetSizeX

double CIWidget::GetPosX() {
  return m_dPosX;
} // GetPosX

double CIWidget::GetPosY() {
  return m_dPosY;
} // GetPosX

double CIWidget::GetActivePosX() {
  return m_dActivePosX;
} // GetSizeX

double CIWidget::GetActivePosY() {
  return m_dActivePosY;
} // GetSizeX

double CIWidget::GetActiveSizeX() {
  return m_dActiveSizeX;
} // GetSizeX

double CIWidget::GetActiveSizeY() {
  return m_dActiveSizeY;
} // GetSizeX

void CIWidget::SetActive(bool bActive) {
  m_bActive = bActive;
} // SetActive

bool CIWidget::Active() {
  return m_bActive;
} // Active

void CIWidget::SetChildrenNotActive() {
  CIWidget *pwCurrent = NULL;

  // Check all the widgets
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    pwCurrent->SetActive(false);
    pwCurrent->SetChildrenNotActive();
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent
  return ;
} // SetChildrenNotActive

void CIWidget::SetVisible(bool bVisible) {
  CIWidget *pwCurrent = NULL;
  m_bVisible = bVisible;

  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    pwCurrent->SetVisible(bVisible);
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent
} // SetVisible

bool CIWidget::Visible() {
  return m_bVisible;
} // Visible

bool CIWidget::InRange(double dPosX, double dPosY) {
  bool bReturnValue = false;

  if(GetPosX() < dPosX &&
      GetPosX() + GetSizeX() > dPosX &&
      GetPosY() < dPosY &&
      GetPosY() + GetSizeY() > dPosY){
    bReturnValue = true;
  } // if in range

  return bReturnValue;
} // InRange

bool CIWidget::InActiveRange(double dPosX, double dPosY) {
  bool bReturnValue = false;

  if(GetPosX() < dPosX &&
      GetPosX() + GetActiveSizeX() > dPosX &&
      GetPosY() < dPosY &&
      GetPosY() + GetActiveSizeY() > dPosY){
    bReturnValue = true;
  } // if Active

  return bReturnValue;
} // InActiveRange

/* This function adds a child widget
 *
 * TODO: This only checks for competition on the window-level, not on the child
 * level. Child widgets naturally compete on the window-level, by design
 * */
void CIWidget::AddWidget(CIWidget *pwWidget) {
  CIWidget *pwCompetition;
  pwCompetition = GetWidgetInRange(pwWidget->GetPosX(), pwWidget->GetPosY(), pwWidget->GetSizeX(), pwWidget->GetSizeY());
  if(pwCompetition) {
    printf("WARNING: Widget cannot be placed there: competing widget\n");
    return;
  } // if pwCompetition

  m_wlChildren.AddWidget(pwWidget);
  pwWidget->SetParent(this);
  pwWidget->SetParent(m_pwParentWindow);
} // AddWidget

CIWidget *CIWidget::GetActiveWidget() {
  CIWidget *pwCurrent = m_wlChildren.GetFirst();
  CIWidget *pwActive = NULL,
	   *pwActiveChild = NULL;

  while(pwCurrent) {
    if(pwCurrent->Active()) {
      pwActive = pwCurrent;
      break;
    } // if Active
    pwActiveChild = pwCurrent->GetActiveWidget();
    if(pwActiveChild) {
      pwActive = pwActiveChild;
      break;
    } // if pwActiveChild
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  return pwActive;
} // GetActiveWidget

CIWidget *CIWidget::GetWidget(double dPosX, double dPosY) {
  CIWidget *pwCurrent = NULL;
  CIWidget *pwTarget = NULL;

  // First check for the active widget, that one is always on top
  pwCurrent = GetActiveWidget();
  if(pwCurrent) {
    if(pwCurrent->GetActivePosX() < dPosX &&
	pwCurrent->GetActivePosX() + pwCurrent->GetActiveSizeX() > dPosX &&
	pwCurrent->GetActivePosY() < dPosY &&
	pwCurrent->GetActivePosY() + pwCurrent->GetActiveSizeY() > dPosY &&
	pwCurrent->Visible()) {
      return pwCurrent;
    } // if Active
  } // if pwCurrent

  // It wasn't the active widget. Check all the others
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    if(pwCurrent->GetPosX() < dPosX &&
	pwCurrent->GetPosX() + pwCurrent->GetSizeX() > dPosX &&
	pwCurrent->GetPosY() < dPosY &&
	pwCurrent->GetPosY() + pwCurrent->GetSizeY() > dPosY &&
	pwCurrent->Visible()) {
      pwTarget = pwCurrent;
      break;
    } // if Active
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  return pwTarget;
} // GetWidget

CIWidget *CIWidget::GetNonActiveWidget(double dPosX, double dPosY) {
  CIWidget *pwCurrent = NULL;
  CIWidget *pwTarget = NULL;

  // We don't check for the active widget
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    if(pwCurrent->GetPosX() < dPosX &&
	pwCurrent->GetPosX() + pwCurrent->GetSizeX() > dPosX &&
	pwCurrent->GetPosY() < dPosY &&
	pwCurrent->GetPosY() + pwCurrent->GetSizeY() > dPosY &&
	pwCurrent->Visible()) {
      pwTarget = pwCurrent;
      break;
    } // if Active
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  return pwTarget;
} // GetNonActiveWidget

CIWidget *CIWidget::GetWidgetInRange(double dPosX, double dPosY, double dSizeX, double dSizeY) {
  CIWidget *pwCurrent = NULL;
  CIWidget *pwTarget = NULL;

  // Check all the widgets
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    if(pwCurrent->InRange(dPosX, dPosY) ||
	pwCurrent->InRange(dPosX + dSizeX, dPosY) ||
	pwCurrent->InRange(dPosX, dPosY + dSizeY) ||
	pwCurrent->InRange(dPosX + dSizeX, dPosY + dSizeY) &&
	pwCurrent->Visible()) {
      pwTarget = pwCurrent;
      break;
    } // if InRange
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  return pwTarget;
} // GetWidgetInRange

bool CIWidget::HasWidget(CIWidget *pwWidget) {
  return m_wlChildren.HasWidget(pwWidget);
} // HasWidget

void CIWidget::RemoveWidget(CIWidget *pwWidget) {
  m_wlChildren.RemoveWidget(pwWidget);
} // RemoveWidget

void CIWidget::SetParent(CIWindow *pwWindow) {
  m_pwParentWindow = pwWindow;
} // SetParent

void CIWidget::SetParent(CIWidget *pwWidget) {
  m_pwParentWidget = pwWidget;
} // SetParent

/* This is a virtual function, supposed to be overloaded by the derived class.
 * It is supposed to draw this widget to the screen
 * */
void CIWidget::Update() {
} // Update

/* This function updates all the children of this widget
 * */
void CIWidget::UpdateChildren() {
  CIWidget *pwCurrent = NULL;
  pwCurrent = m_wlChildren.GetFirst();
  while(pwCurrent) {
    if(pwCurrent->Visible()) {
      pwCurrent->Update();
      pwCurrent->UpdateChildren();
    } // if Visible
    pwCurrent = m_wlChildren.GetNext(pwCurrent);
  } // while pwCurrent

  pwCurrent = GetActiveWidget();
  if(pwCurrent) {
    if(pwCurrent->Visible()) {
      pwCurrent->Update();
      pwCurrent->UpdateChildren();
    } // if Visible
  } // if pwCurrent
} // Update


CIButton::CIButton() {
  m_pOnClicked = NULL;
  strcpy(m_strText, "");
} // CIButton


CIButton::~CIButton() {
  m_pOnClicked = NULL;
  strcpy(m_strText, "");
} // ~CIButton

/* This function connects the given function to the signal handler for the click
 * of a mouse button. When this connection is made, an arbitrary function can be
 * called when the mouse button is pressed */
void CIButton::ConnectOnClicked(void (*pFunc)(CIWindow *pwThis)) {
  m_pOnClicked = pFunc;
} // Settext

void CIButton::SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY) {
  m_dPosX = dPosX;
  m_dPosY = dPosY;
  m_dSizeX = dSizeX;
  m_dSizeY = dSizeY;

  m_dActivePosX = dPosX;
  m_dActivePosY = dPosY;
  m_dActiveSizeX = dSizeX;
  m_dActiveSizeY = dSizeY;
} // SetPos

void CIButton::SetText(const char strText[]) {
  strcpy(m_strText, strText);
} // Settext

/* This function draws the button to the screen, and puts some text on it
 * */
void CIButton::Update() {
#ifdef HAVE_PGPLOT
  cpgsci(IUGUI_COLOR_BUTTON);
  cpgrect(m_dPosX, m_dPosX+m_dSizeX, m_dPosY, m_dPosY+m_dSizeY);
  cpgsci(IUGUI_COLOR_BUTTONTEXT);
  if(strlen(m_strText) > 0)
    cpgtext(m_dPosX, m_dPosY+0.01, m_strText);
  cpgsci(IUGUI_COLOR_BKGRD);
#endif // HAVE_PGPLOT
} // Update

/* This function calls the signal handler, if it exists, for button clicks
 * */
void CIButton::OnClicked(double dPosX, double dPosY) {
  if(m_pOnClicked)
    m_pOnClicked(m_pwParentWindow);
  if(m_pwParentWindow)
    m_pwParentWindow->SetChildrenNotActive();
} // OnClicked

/* This function doesn't do anything yet, but it should respond to key-strokes
 * */
void CIButton::OnKey(char cKey) {
} // OnKey

/* The constructor of the CIGroup widget
 * */
CIGroup::CIGroup() {
  m_dPosX = 0;
  m_dSizeX = 1;
  m_dActivePosX = 0;
  m_dActiveSizeX = 1;
  m_dPosY = 0;
  m_dSizeY = 1;
  m_dActivePosY = 0;
  m_dActiveSizeY = 1;

  m_bDrawBorder = true;

  strcpy(m_strText, "");
  m_dTextWidth = 0;
  m_dTextPadding = 0.01;
} // CIGroup

CIGroup::~CIGroup() {
  m_bDrawBorder = false;

  strcpy(m_strText, "");
  m_dTextWidth = 0;
  m_dTextPadding = 0.01;
} // CIGroup

void CIGroup::SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY) {
  m_dPosX = dPosX;
  m_dPosY = dPosY;
  m_dSizeX = dSizeX;
  m_dSizeY = dSizeY;

  m_dActivePosX = dPosX;
  m_dActivePosY = dPosY;
  m_dActiveSizeX = dSizeX;
  m_dActiveSizeY = dSizeY;
} // SetPos

/* This function sets the flag whether to draw the border around the group or
 * not
 * */
void CIGroup::SetBorder(bool bShowBorder) {
  m_bDrawBorder = bShowBorder;
} // SetBorder

/* This function returns the state of the border (draw or don't draw)
 * */
bool CIGroup::GetBorder() {
  return m_bDrawBorder;
} // GetBorder

/* This function sets the info text of the group
 * */
void CIGroup::SetText(const char strText[]) {
  if(strText)
    strcpy(m_strText, strText);
  else
    strcpy(m_strText, "");

  m_dTextWidth = strlen(m_strText) * 0.0105;
} // SetText

/* This function sets the text width so that the line is interrupted where the
 * text is written.
 *
 * TODO: Make this more user-friendly
 * */
void CIGroup::SetTextWidth(double dTextWidth) {
  m_dTextWidth = dTextWidth;
} // SetTextWidth

/* This function adds a widget to the tab.
 * */
void CIGroup::AddWidget(CIWidget *pwWidget) {
  CIWidget *pwCompetition;
  pwCompetition = GetWidgetInRange(pwWidget->GetPosX(), pwWidget->GetPosY(), pwWidget->GetSizeX(), pwWidget->GetSizeY());
  if(pwCompetition) {
    // Check whether it is on the same tab
    if(m_wlChildren.HasWidget(pwWidget)) {
      printf("WARNING: Widget cannot be placed there: competing widget\n");
      return;
    } // if HasWidget
  } // if pwCompetition

  // Check whether widget fits in this container
  if(pwWidget->GetPosX() > m_dPosX &&
      pwWidget->GetPosX() + pwWidget->GetSizeX() < m_dPosX + m_dSizeX &&
      pwWidget->GetPosY() > m_dPosY &&
      pwWidget->GetPosY() + pwWidget->GetSizeY() < m_dPosY + m_dSizeY) {
    m_wlChildren.AddWidget(pwWidget);
    pwWidget->SetParent(this);
    pwWidget->SetParent(m_pwParentWindow);
#if 0
    if(nTab == m_nActiveTab)
      pwWidget->SetVisible(true);
    else
      pwWidget->SetVisible(false);
#endif
  } else {
    printf("WARNING: Widget cannot be placed there: doesn't fit in container\n");
  } // if GetPos
} // AddWidget

/* This function draws the line around the group, if the m_bDrawBorder flag is
 * set to true
 * */
void CIGroup::Update() {
  float pfLineX[7], pfLineY[7];
  float fPosX, fPosY;
  float fParentFontSize = m_pwParentWindow->GetFontSize();

  if(m_bDrawBorder) {
    pfLineX[0] = float(m_dPosX);
    pfLineX[1] = float(m_dPosX);
    pfLineX[2] = float(m_dPosX + m_dTextPadding);
    pfLineX[3] = float(m_dPosX + m_dTextPadding + m_dTextWidth);
    pfLineX[4] = float(m_dPosX + m_dSizeX);
    pfLineX[5] = float(m_dPosX + m_dSizeX);
    pfLineX[6] = float(m_dPosX);
    pfLineY[0] = float(m_dPosY);
    pfLineY[1] = float(m_dPosY + m_dSizeY);
    pfLineY[2] = float(m_dPosY + m_dSizeY);
    pfLineY[3] = float(m_dPosY + m_dSizeY);
    pfLineY[4] = float(m_dPosY + m_dSizeY);
    pfLineY[5] = float(m_dPosY);
    pfLineY[6] = float(m_dPosY);

    // Draw the box
#ifdef HAVE_PGPLOT
    cpgsci(IUGUI_COLOR_LINE);
    cpgline(3, pfLineX, pfLineY);
    cpgline(4, pfLineX+3, pfLineY+3);
    cpgsci(IUGUI_COLOR_BKGRD);
#endif // HAVE_PGPLOT
  } // if m_bDrawBorder

#ifdef HAVE_PGPLOT
  // Draw the text (if any)
  cpgsci(IUGUI_COLOR_TEXT);
  cpgsch(fParentFontSize); 
  if(strlen(m_strText) > 0)
    cpgtext(m_dPosX+m_dTextPadding, m_dPosY+m_dSizeY-0.001, m_strText);
  cpgsch(fParentFontSize); 
  cpgsci(IUGUI_COLOR_BKGRD);
#endif // HAVE_PGPLOT
} // Update


/* The user has clicked in the area of the group. Check whether it was on a
 * child widget
 * */
void CIGroup::OnClicked(double dPosX, double dPosY) {
  CIWidget *pwTarget;
  double dButtonPosX;

  pwTarget = GetWidget(dPosX, dPosY);
  if(pwTarget) {
    pwTarget->OnClicked(dPosX, dPosY);
  } else {
    // What to do if not clicked on a widget in the group? Set here
    if(m_pwParentWindow)
      m_pwParentWindow->SetChildrenNotActive();
  } // if pwTarget
} // OnClicked

/* Key pressed
 * */
void CIGroup::OnKey(char cKey) {
  ;
} // OnKey


CITab::CITab() {
  m_nTabs = 0;
  m_nActiveTab = 0;
  m_dMenuHeight = 0.03;
  m_dMenuItemWidth = 0.12;
  m_dMenuItemWidthSpace = 0.01;
  m_dPosX = 0;
  m_dSizeX = 1;
  m_dActivePosX = 0;
  m_dActiveSizeX = 1;
  m_dPosY = 0;
  m_dSizeY = 1;
  m_dActivePosY = 0;
  m_dActiveSizeY = 1;
} // CITab

CITab::~CITab() {
  m_nTabs = 0;
  m_nActiveTab = 0;
  m_dMenuHeight = 0.03;
  m_dMenuItemWidth = 0.09;
  m_dMenuItemWidthSpace = 0.01;
} // CITab

void CITab::SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY) {
  m_dPosX = dPosX;
  m_dPosY = dPosY;
  m_dSizeX = dSizeX;
  m_dSizeY = dSizeY;

  m_dActivePosX = dPosX;
  m_dActivePosY = dPosY;
  m_dActiveSizeX = dSizeX;
  m_dActiveSizeY = dSizeY;
} // SetPos

void CITab::SetTabs(int nTabs) {
  m_nTabs = nTabs;
} // SetTabs

/* This function returns the number of tabs
 * */
int CITab::GetTabs() {
  return m_nTabs;
} // GetTabs

void CITab::SetText(int nTab, const char strText[]) {
  strcpy(m_strTabText[nTab], strText);
} // SetText

void CITab::SetActiveTab(int nTab) {
  CIWidget *pwCurrent = NULL;
  m_nActiveTab = nTab;
  for(int i=0; i<m_nTabs; i++) {
    pwCurrent = m_wlTabChildren[i].GetFirst();
    while(pwCurrent) {
      if(i == m_nActiveTab)
	pwCurrent->SetVisible(true);
      else
	pwCurrent->SetVisible(false);
      pwCurrent = m_wlTabChildren[i].GetNext(pwCurrent);
    } // while pwCurrent
  } // for i
} // SetActiveTab

/* This function returns the active tab id
 * */
int CITab::GetActiveTab() {
  return m_nActiveTab;
} // GetActiveTab

void CITab::SetVisible(bool bVisible) {
  CIWidget *pwCurrent = NULL;
  m_bVisible = bVisible;

  for(int i=0; i<m_nTabs; i++) {
    pwCurrent = m_wlTabChildren[i].GetFirst();
    while(pwCurrent) {
      if(i == m_nActiveTab)
	pwCurrent->SetVisible(true);
      else
	pwCurrent->SetVisible(false);
      pwCurrent = m_wlTabChildren[i].GetNext(pwCurrent);
    } // while pwCurrent
  } // for i
} // SetVisible

/* This function changes the menu item width
 * */
void CITab::SetMenuItemWidth(double dMenuItemWidth) {
  m_dMenuItemWidth = dMenuItemWidth;
} // SetMenuItemWidth

/* This function adds a widget to the tab.
 *
 * TODO: This function doesn't work well enough. It only checks for one
 * competing widget, but there can be many on different tabs. Make the check per
 * tab. SOLVED: Checks are only for _visible_ widgets.
 * */
void CITab::AddWidget(int nTab, CIWidget *pwWidget) {
  CIWidget *pwCompetition;
  pwCompetition = GetWidgetInRange(pwWidget->GetPosX(), pwWidget->GetPosY(), pwWidget->GetSizeX(), pwWidget->GetSizeY());
  if(pwCompetition) {
    // Check whether it is on the same tab
    if(m_wlTabChildren[nTab].HasWidget(pwWidget)) {
      printf("WARNING: Widget cannot be placed there: competing widget\n");
      return;
    } // if HasWidget
  } // if pwCompetition

  // Check whether widget fits in this container
  if(pwWidget->GetPosX() > m_dPosX &&
      pwWidget->GetPosX() + pwWidget->GetSizeX() < m_dPosX + m_dSizeX &&
      pwWidget->GetPosY() > m_dPosY &&
      pwWidget->GetPosY() + pwWidget->GetSizeY() < m_dPosY + m_dSizeY - m_dMenuHeight) {
    m_wlChildren.AddWidget(pwWidget);
    m_wlTabChildren[nTab].AddWidget(pwWidget);
    pwWidget->SetParent(this);
    pwWidget->SetParent(m_pwParentWindow);
    if(nTab == m_nActiveTab)
      pwWidget->SetVisible(true);
    else
      pwWidget->SetVisible(false);
  } else {
    printf("WARNING: Widget cannot be placed there: doesn't fit in container\n");
  } // if GetPos
  SetActiveTab(m_nActiveTab);
} // AddWidget

void CITab::RemoveWidget(CIWidget *pwWidget) {
  m_wlChildren.RemoveWidget(pwWidget);
  for(int i=0; i<m_nTabs; i++) {
    if(m_wlTabChildren[i].HasWidget(pwWidget))
      m_wlTabChildren[i].RemoveWidget(pwWidget);
  } // for i
} // RemoveWidget

void CITab::Update() {
  float pfLineX[5], pfLineY[5];
  float fPosX, fPosY;
  pfLineX[0] = float(m_dPosX);
  pfLineX[1] = float(m_dPosX);
  pfLineX[2] = float(m_dPosX + m_dSizeX);
  pfLineX[3] = float(m_dPosX + m_dSizeX);
  pfLineX[4] = float(m_dPosX);
  pfLineY[0] = float(m_dPosY);
  pfLineY[1] = float(m_dPosY + m_dSizeY - m_dMenuHeight);
  pfLineY[2] = float(m_dPosY + m_dSizeY - m_dMenuHeight);
  pfLineY[3] = float(m_dPosY);
  pfLineY[4] = float(m_dPosY);

  // Draw the box
#ifdef HAVE_PGPLOT
  cpgsci(IUGUI_COLOR_LINE);
  cpgline(5, pfLineX, pfLineY);
  cpgsci(IUGUI_COLOR_BKGRD);
#endif // HAVE_PGPLOT

  // Draw the menu on top
  for(int i=0; i<m_nTabs; i++) {
    // Start position for button
    fPosX = float(m_dPosX + i*(m_dMenuItemWidth+m_dMenuItemWidthSpace));
    fPosY = float(m_dPosY + m_dSizeY - m_dMenuHeight);

    // Draw the rectangle
    if(i == m_nActiveTab) {
#ifdef HAVE_PGPLOT
      cpgsci(IUGUI_COLOR_BUTTON);
#endif // HAVE_PGPLOT
    } else {
#ifdef HAVE_PGPLOT
      cpgsci(IUGUI_COLOR_LINE);
#endif // HAVE_PGPLOT
    } // if m_nActiveTab
#ifdef HAVE_PGPLOT
    cpgrect(fPosX, fPosX + float(m_dMenuItemWidth), fPosY, fPosY + float(m_dMenuHeight));
#endif // HAVE_PGPLOT

#ifdef HAVE_PGPLOT
    cpgsci(IUGUI_COLOR_BUTTONTEXT);
    if(strlen(m_strTabText[i]) > 0)
      cpgtext(fPosX, fPosY+0.01, m_strTabText[i]);

    cpgsci(IUGUI_COLOR_BKGRD);
#endif // HAVE_PGPLOT
  } // for i
} // Update


void CITab::OnClicked(double dPosX, double dPosY) {
  CIWidget *pwTarget;
  double dButtonPosX;

  pwTarget = GetWidget(dPosX, dPosY);
  if(pwTarget) {
    pwTarget->OnClicked(dPosX, dPosY);
  } else {
    if(dPosY > m_dPosY + m_dSizeY - m_dMenuHeight &&
      dPosY < m_dPosY + m_dSizeY) {
      // Check whether it was on a tab-button
      for(int i=0; i<m_nTabs; i++) {
	dButtonPosX = m_dPosX + i*(m_dMenuItemWidth+m_dMenuItemWidthSpace);

	if(dPosX > dButtonPosX && dPosX < dButtonPosX + m_dMenuItemWidth) {
	  SetActiveTab(i);
	  break;
	} // if dPosX
      } // for i
    } // if dPos
    if(m_pwParentWindow)
      m_pwParentWindow->SetChildrenNotActive();
  } // if pwTarget
} // OnClicked

void CITab::OnKey(char cKey) {
  ;
} // OnKey


CIComboBox::CIComboBox() {
  m_nItems = 0;
  m_nCurrentItem = 0;
  m_pOnSelected = NULL;
  m_dPosX = 0;
  m_dSizeX = 0.1;
  m_dActivePosX = 0;
  m_dActiveSizeX = 0.1;
  m_dPosY = 0;
  m_dSizeY = 0.03;
  m_dActivePosY = 0;
  m_dActiveSizeY = 0.03;
  m_dItemHeight = 0.03;
  m_dItemHeightSpace = 0.005;

  for(int i=0; i<MAX_COMBOITEMS; i++) {
    strcpy(m_strItemText[i], "");
  } // for i
} // CIComboBox


CIComboBox::~CIComboBox() {
  m_nItems = 0;
  m_nCurrentItem = 0;
  m_pOnSelected = NULL;
} // ~CIButton

/* This function connects the given function to the signal handler for the
 * selection of an element of the combobox. When this connection is made, an
 * arbitrary function can be called when an element is selected
 * */
void CIComboBox::ConnectOnSelected(void (*pFunc)(CIWindow *pwThis)) {
  m_pOnSelected = pFunc;
} // ConnectOnSelected

void CIComboBox::SetPos(double dPosX, double dPosY, double dWidth) {
  m_dPosX = dPosX;
  m_dActivePosX = dPosX;
  m_dPosY = dPosY;
  m_dActivePosY = dPosY;

  if(dWidth != -1) {
    m_dSizeX = dWidth;
    m_dActiveSizeX = dWidth;
  } // if dWidth
} // SetPos

void CIComboBox::SetItems(int nItems) {
  if(nItems > MAX_COMBOITEMS) {
    printf("WARNING: Max combobox items is set to: %i\n", MAX_COMBOITEMS);
    m_nItems = MAX_COMBOITEMS;
  } else {
    m_nItems = nItems;
  } // if nItems
} // SetItems

void CIComboBox::SetItem(int nItem, char strItem[]) {
  if(nItem < MAX_COMBOITEMS)
    strcpy(m_strItemText[nItem], strItem);
} // SetItem

/* This function sets the current item to a new value. It must _not_ call the
 * signal handler as if it has been clicked, since that would lead to bugs.
 * */
void CIComboBox::SetCurrentItem(int nItem) {
  if(nItem < m_nItems)
    m_nCurrentItem = nItem;
} // SetCurrentItem

/* This function returns the index of the current item
 * */
int CIComboBox::GetCurrentItem() {
  return m_nCurrentItem;
} // GetCurrentItem

/* This function copies the text of item nItem into the parameter strText
 * */
void CIComboBox::CopyItemText(int nItem, char *strText) {
  if(! strText) {
    printf("WARNING: cannot copy to NULL pointer\n");
    return;
  } // if strText
  if(nItem >= m_nItems) {
    printf("WARNING: requested item out of range\n");
    return;
  } // if nItem

  strcpy(strText, m_strItemText[nItem]);
  return;
} // GetCurrentItem

/* This function draws the combobox to the screen
 * */
void CIComboBox::Update() {
  float pfLineX[5], pfLineY[5];
  float fPosX, fPosY;

  // First, draw the primary select-box
  pfLineX[0] = float(m_dPosX);
  pfLineX[1] = float(m_dPosX);
  pfLineX[2] = float(m_dPosX + m_dSizeX);
  pfLineX[3] = float(m_dPosX + m_dSizeX);
  pfLineX[4] = float(m_dPosX);
  pfLineY[0] = float(m_dPosY);
  pfLineY[1] = float(m_dPosY + m_dSizeY);
  pfLineY[2] = float(m_dPosY + m_dSizeY);
  pfLineY[3] = float(m_dPosY);
  pfLineY[4] = float(m_dPosY);

#ifdef HAVE_PGPLOT
  cpgsci(IUGUI_COLOR_BUTTON);
  cpgline(5, pfLineX, pfLineY);

  cpgsci(IUGUI_COLOR_LINE);
  if(strlen(m_strItemText[m_nCurrentItem]) > 0)
    cpgtext(m_dPosX, m_dPosY+0.01, m_strItemText[m_nCurrentItem]);
#endif // HAVE_PGPLOT

  switch(m_bActive) {
    case true:
      m_dActiveSizeY = m_dSizeY + m_nItems*(m_dItemHeightSpace+m_dItemHeight);
      m_dActivePosY = m_dPosY + m_dSizeY - m_dActiveSizeY;

      // Draw the active part of the select box -- fill
#ifdef HAVE_PGPLOT
      cpgsci(IUGUI_COLOR_LINE);
      cpgrect(m_dPosX, m_dPosX+m_dSizeX, m_dPosY + m_dSizeY - m_dActiveSizeY, m_dPosY);
#endif // HAVE_PGPLOT

      // Draw the active part of the select-box -- border
      pfLineX[0] = float(m_dPosX);
      pfLineX[1] = float(m_dPosX);
      pfLineX[2] = float(m_dPosX + m_dSizeX);
      pfLineX[3] = float(m_dPosX + m_dSizeX);
      pfLineX[4] = float(m_dPosX);
      pfLineY[0] = float(m_dPosY + m_dSizeY - m_dActiveSizeY);
      pfLineY[1] = float(m_dPosY);
      pfLineY[2] = float(m_dPosY);
      pfLineY[3] = float(m_dPosY + m_dSizeY - m_dActiveSizeY);
      pfLineY[4] = float(m_dPosY + m_dSizeY - m_dActiveSizeY);
#ifdef HAVE_PGPLOT
      cpgsci(IUGUI_COLOR_BUTTON);
      cpgline(5, pfLineX, pfLineY);
#endif // HAVE_PGPLOT

      // Now write the items
      for(int i=0; i<m_nItems; i++) {
#ifdef HAVE_PGPLOT
	cpgsci(IUGUI_COLOR_BKGRD);
	if(strlen(m_strItemText[i]) > 0)
	  cpgtext(m_dPosX, m_dPosY + 0.01 - (i+1)*(m_dItemHeight+m_dItemHeightSpace), m_strItemText[i]);
#endif // HAVE_PGPLOT

	// Draw the border between elements
	pfLineX[0] = float(m_dPosX);
	pfLineX[1] = float(m_dPosX + m_dSizeX);
	pfLineY[0] = float(m_dPosY - (i+1)*(m_dItemHeight+m_dItemHeightSpace));
	pfLineY[1] = float(m_dPosY - (i+1)*(m_dItemHeight+m_dItemHeightSpace));
#ifdef HAVE_PGPLOT
	cpgsci(IUGUI_COLOR_BUTTON);
	cpgline(2, pfLineX, pfLineY);
#endif // HAVE_PGPLOT
      } // for i
      break;
    case false:
    default:
      break;
  } // switch m_bActive
} // SetUpdate


/* This function calls the signal handler, if it exists, for button clicks
 * */
void CIComboBox::OnClicked(double dPosX, double dPosY) {
  double dOffset, dSingleOffset;
  int nProposal;
  if(m_bActive) {
    if(m_pwParentWindow)
      m_pwParentWindow->SetChildrenNotActive();
    m_bActive = false;

    // Check where was clicked
    dOffset = (dPosY - m_dActivePosY);
    dSingleOffset = dOffset - int(dOffset / (m_dItemHeight + m_dItemHeightSpace)) *
      (m_dItemHeight + m_dItemHeightSpace);
    if(dSingleOffset < m_dItemHeight) {
      // We clicked right on the mark
      nProposal = int(dOffset / (m_dItemHeight + m_dItemHeightSpace));
      if(nProposal >= 0 && nProposal < m_nItems) {
	// New item selected
	m_nCurrentItem = m_nItems - 1 - nProposal;

	// If set, call the callback function
	if(m_pOnSelected)
	  m_pOnSelected(m_pwParentWindow);
      } // if nProposal
    } // if dSingleOffset
  } else {
    if(m_pwParentWindow)
      m_pwParentWindow->SetChildrenNotActive();
    m_bActive = true;
  } // if m_bActive
} // OnClicked

/* This function doesn't do anything yet, but it should respond to key-strokes
 * */
void CIComboBox::OnKey(char cKey) {
} // OnKey



CICheckbox::CICheckbox() {
  m_pOnClicked = NULL;
  m_bChecked = false;
  m_dPadding = 0.0025;
  m_dBorderSize = 0.0025;
  m_dPosX = 0;
  m_dSizeX = 0.03;
  m_dActiveSizeX = 0.03;
  m_dPosY = 0;
  m_dSizeY = 0.03;
  m_dActiveSizeY = 0.03;
} // CICheckbox


CICheckbox::~CICheckbox() {
  m_pOnClicked = NULL;
  m_bChecked = false;
} // ~CICheckbox

/* This function connects the given function to the signal handler for the click
 * of a mouse button. When this connection is made, an arbitrary function can be
 * called when the mouse button is pressed */
void CICheckbox::ConnectOnClicked(void (*pFunc)(CIWindow *pwThis)) {
  m_pOnClicked = pFunc;
} // Settext

void CICheckbox::SetPos(double dPosX, double dPosY) {
  m_dPosX = dPosX;
  m_dPosY = dPosY;

  m_dActivePosX = dPosX;
  m_dActivePosY = dPosY;
} // SetPos

void CICheckbox::SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY) {
  m_dPosX = dPosX;
  m_dPosY = dPosY;
  m_dSizeX = dSizeX;
  m_dSizeY = dSizeY;

  m_dActivePosX = dPosX;
  m_dActivePosY = dPosY;
  m_dActiveSizeX = dSizeX;
  m_dActiveSizeY = dSizeY;
} // SetPos

void CICheckbox::SetChecked(bool bChecked) {
  m_bChecked = bChecked;
} // SetChecked

bool CICheckbox::GetChecked() {
  return m_bChecked;
} // GetChecked

/* This function draws the button to the screen, and puts some text on it
 * */
void CICheckbox::Update() {
#ifdef HAVE_PGPLOT
  cpgsci(IUGUI_COLOR_BUTTON);
  cpgrect(m_dPosX, m_dPosX+m_dSizeX, m_dPosY, m_dPosY+m_dSizeY);
  cpgsci(IUGUI_COLOR_BKGRD);
  cpgrect(m_dPosX+m_dBorderSize, m_dPosX+m_dSizeX-m_dBorderSize,
      m_dPosY+m_dBorderSize, m_dPosY+m_dSizeY-m_dBorderSize);
  if(m_bChecked) {
    cpgsci(IUGUI_COLOR_ACTIVE);
    cpgrect(m_dPosX+m_dPadding+m_dBorderSize, m_dPosX+m_dSizeX-m_dPadding-m_dBorderSize,
	m_dPosY+m_dPadding+m_dBorderSize, m_dPosY+m_dSizeY-m_dPadding-m_dBorderSize);
  } // if m_bChecked
  cpgsci(IUGUI_COLOR_BKGRD);
#endif // HAVE_PGPLOT
} // Update

/* This function calls the signal handler, if it exists, for button clicks
 * */
void CICheckbox::OnClicked(double dPosX, double dPosY) {
  m_bChecked = m_bChecked ? false : true;

  if(m_pOnClicked)
    m_pOnClicked(m_pwParentWindow);
  if(m_pwParentWindow)
    m_pwParentWindow->SetChildrenNotActive();
} // OnClicked

/* This function doesn't do anything yet, but it should respond to key-strokes
 * */
void CICheckbox::OnKey(char cKey) {
} // OnKey



CIPlot::CIPlot() {
  m_pOnClicked = NULL;
  m_dPosX = 0;
  m_dSizeX = 1;
  m_dActivePosX = 0;
  m_dActiveSizeX = 1;
  m_dPosY = 0;
  m_dSizeY = 1;
  m_dActivePosY = 0;
  m_dActiveSizeY = 1;
  strcpy(m_strTitle, "");
  strcpy(m_strXLabel, "");
  strcpy(m_strYLabel, "");

  m_nSymbol = 16;
  m_bUserDefinedView = false;

  m_bLogX = false;
  m_bLogY = false;

  SetNull();
} // CIPlot


CIPlot::~CIPlot() {
  m_pOnClicked = NULL;
  strcpy(m_strTitle, "");
  strcpy(m_strXLabel, "");
  strcpy(m_strYLabel, "");

  DeletePlot();
} // ~CIPlot

/* This function connects the given function to the signal handler for the click
 * of a mouse button. When this connection is made, an arbitrary function can be
 * called when the mouse button is pressed */
void CIPlot::ConnectOnClicked(void (*pFunc)(CIWindow *pwThis)) {
  m_pOnClicked = pFunc;
} // Settext

void CIPlot::SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY) {
  m_dPosX = dPosX;
  m_dPosY = dPosY;
  m_dSizeX = dSizeX;
  m_dSizeY = dSizeY;

  m_dActivePosX = dPosX;
  m_dActivePosY = dPosY;
  m_dActiveSizeX = dSizeX;
  m_dActiveSizeY = dSizeY;
} // SetPos

void CIPlot::SetPlot(int nPlotPoints, double *pdX, double *pdY) {
  DeletePlot();
  m_nPlotPoints = nPlotPoints;
  m_pfPlotX = new float[m_nPlotPoints];
  m_pfPlotY = new float[m_nPlotPoints];

  for(int i=0; i<m_nPlotPoints; i++) {
    if(m_bLogX)
      m_pfPlotX[i] = float(log10(pdX[i]));
    else
      m_pfPlotX[i] = float(pdX[i]);
    if(m_bLogY)
      m_pfPlotY[i] = float(log10(pdY[i]));
    else
      m_pfPlotY[i] = float(pdY[i]);
  } // for i
} // SetPlot

void CIPlot::SetXErr(int nPlotPoints, double *pdXErr) {
  if(m_nPlotPoints != nPlotPoints) {
    printf("WARNING: m_nPlotPoints != nPlotPoints\n");
    return;
  } // if m_nPlotPoints
  if(! m_pfPlotX) {
    printf("WARNING: can only set errors after defining points\n");
    return;
  } // if m_pfPlotX

  if(m_pfPlotXErr)
    delete[] m_pfPlotXErr;
  if(m_pfPlotXErr1)
    delete[] m_pfPlotXErr1;
  if(m_pfPlotXErr2)
    delete[] m_pfPlotXErr2;

  m_pfPlotXErr = new float[m_nPlotPoints];
  m_pfPlotXErr1 = new float[m_nPlotPoints];
  m_pfPlotXErr2 = new float[m_nPlotPoints];

  for(int i=0; i<m_nPlotPoints; i++) {
    m_pfPlotXErr[i] = float(pdXErr[i]);
    m_pfPlotXErr1[i] = m_pfPlotX[i] - 0.5*m_pfPlotXErr[i];
    m_pfPlotXErr2[i] = m_pfPlotX[i] + 0.5*m_pfPlotXErr[i];
  } // for i
} // SetPlotXErr

void CIPlot::SetYErr(int nPlotPoints, double *pdYErr) {
  if(m_nPlotPoints != nPlotPoints) {
    printf("WARNING: m_nPlotPoints != nPlotPoints\n");
    return;
  } // if m_nPlotPoints
  if(! m_pfPlotY) {
    printf("WARNING: can only set errors after defining points\n");
    return;
  } // if m_pfPlotY

  if(m_pfPlotYErr)
    delete[] m_pfPlotYErr;
  if(m_pfPlotYErr1)
    delete[] m_pfPlotYErr1;
  if(m_pfPlotYErr2)
    delete[] m_pfPlotYErr2;

  m_pfPlotYErr = new float[m_nPlotPoints];
  m_pfPlotYErr1 = new float[m_nPlotPoints];
  m_pfPlotYErr2 = new float[m_nPlotPoints];

  for(int i=0; i<m_nPlotPoints; i++) {
    m_pfPlotYErr[i] = float(pdYErr[i]);
    m_pfPlotYErr1[i] = m_pfPlotY[i] - 0.5*m_pfPlotYErr[i];
    m_pfPlotYErr2[i] = m_pfPlotY[i] + 0.5*m_pfPlotYErr[i];
  } // for i
} // SetYErr

void CIPlot::SetYFit(int nPlotPoints, double *pdYFit) {
  if(m_nPlotPoints != nPlotPoints) {
    printf("WARNING: m_nPlotPoints != nPlotPoints\n");
    return;
  } // if m_nPlotPoints
  if(! m_pfPlotY) {
    printf("WARNING: can only set errors after defining points\n");
    return;
  } // if m_pfPlotY

  if(m_pfPlotYFit)
    delete[] m_pfPlotYFit;

  m_pfPlotYFit = new float[m_nPlotPoints];

  for(int i=0; i<m_nPlotPoints; i++) {
    m_pfPlotYFit[i] = float(pdYFit[i]);
  } // for i
} // SetYFit

void CIPlot::SetLogScale(int nAxis, bool bLog) {
  switch(nAxis) {
    case 0:
      m_bLogX = bLog;
      break;
    case 1:
      m_bLogY = bLog;
      break;
    default:
      break;
  } // switch nAxis
} // SetLogScale

void CIPlot::DrawLine(bool bLine) {
  m_bDrawLine = bLine;
} // DrawLine

void CIPlot::UnsetYFit() {
  if(m_pfPlotYFit)
    delete[] m_pfPlotYFit;
  m_pfPlotYFit = NULL;
} // UnsetYFit

void CIPlot::SetTitle(const char strText[]) {
  strcpy(m_strTitle, strText);
} // SetTitle

void CIPlot::SetXLabel(const char strText[]) {
  strcpy(m_strXLabel, strText);
} // SetXLabel

void CIPlot::SetYLabel(const char strText[]) {
  strcpy(m_strYLabel, strText);
} // SetYLabel

void CIPlot::SetSymbol(int nSymbol) {
  m_nSymbol = nSymbol;
} // SetSymbol

void CIPlot::SetUserView(double dMinX, double dMaxX, double dMinY, double dMaxY) {
  m_dMinX = dMinX;
  m_dMaxX = dMaxX;
  m_dMinY = dMinY;
  m_dMaxY = dMaxY;
  m_bUserDefinedView = true;
} // SetUserView

void CIPlot::UnSetUserView() {
  m_bUserDefinedView = false;
} // UnSetUserView

void CIPlot::DeletePlot() {
  if(m_pfPlotX)
    delete[] m_pfPlotX;
  if(m_pfPlotY)
    delete[] m_pfPlotY;
  if(m_pfPlotXErr)
    delete[] m_pfPlotXErr;
  if(m_pfPlotXErr1)
    delete[] m_pfPlotXErr1;
  if(m_pfPlotXErr2)
    delete[] m_pfPlotXErr2;
  if(m_pfPlotYErr)
    delete[] m_pfPlotYErr;
  if(m_pfPlotYErr1)
    delete[] m_pfPlotYErr1;
  if(m_pfPlotYErr2)
    delete[] m_pfPlotYErr2;
  if(m_pfPlotYFit)
    delete[] m_pfPlotYFit;
  SetNull();
} // DeletePlot

/* This function saves the plot to a user-defined file or device */
void CIPlot::SavePlot(const char strDeviceName[]) {
  char strDeviceFile[160];
  m_pwParentWindow->GetDevice(strDeviceFile);
  m_pwParentWindow->DeInitialise();

  DeviceInit(strDeviceName);
  Update();
  DeviceDeInit();

  m_pwParentWindow->Initialise(strDeviceFile);
} // SavePlot

void CIPlot::SetNull() {
  m_nPlotPoints = 0;
  m_pfPlotX = NULL;
  m_pfPlotY = NULL;
  m_pfPlotXErr = NULL;
  m_pfPlotXErr1 = NULL;
  m_pfPlotXErr2 = NULL;
  m_pfPlotYErr = NULL;
  m_pfPlotYErr1 = NULL;
  m_pfPlotYErr2 = NULL;
  m_pfPlotYFit = NULL;
} // SetNull

/* This function draws the button to the screen, and puts some text on it
 * */
void CIPlot::Update() {
  double dMinX, dMaxX, dMinY, dMaxY;
  char strBoxX[10], strBoxY[10];

  if(! m_pfPlotX || ! m_pfPlotY)
    return;

  // Set the viewport and window
#ifdef HAVE_PGPLOT
  if(m_bRespectAbsolutePosition)
    cpgsvp(m_dPosX, m_dPosX+m_dSizeX, m_dPosY, m_dPosY+m_dSizeY);
  else
    cpgsvp(0.1, 0.9, 0.1, 0.9);
#endif // HAVE_PGPLOT

  // Check whether we need to set the user view here
  if(m_bUserDefinedView) {
    dMinX = m_dMinX;
    dMaxX = m_dMaxX;
    dMinY = m_dMinY;
    dMaxY = m_dMaxY;

    dMinX = 0.5 * (m_dMaxX + m_dMinX) - 0.55 * (m_dMaxX - m_dMinX);
    dMaxX = 0.5 * (m_dMaxX + m_dMinX) + 0.55 * (m_dMaxX - m_dMinX);
    dMinY = 0.5 * (m_dMaxY + m_dMinY) - 0.55 * (m_dMaxY - m_dMinY);
    dMaxY = 0.5 * (m_dMaxY + m_dMinY) + 0.55 * (m_dMaxY - m_dMinY);
  } else {
    dMinX = MinX();
    dMaxX = MaxX();
    dMinY = MinY();
    dMaxY = MaxY();

    dMinX = 0.5 * (MaxX() + MinX()) - 0.55 * (MaxX() - MinX());
    dMaxX = 0.5 * (MaxX() + MinX()) + 0.55 * (MaxX() - MinX());
    dMinY = 0.5 * (MaxY() + MinY()) - 0.55 * (MaxY() - MinY());
    dMaxY = 0.5 * (MaxY() + MinY()) + 0.55 * (MaxY() - MinY());
  } // if m_bUserDefinedView

  // Set the boxing stuff
  if(m_bLogX)
    strcpy(strBoxX, "BCNLST2");
  else
    strcpy(strBoxX, "BCNST2");

  if(m_bLogY)
    strcpy(strBoxY, "BCNLST2");
  else
    strcpy(strBoxY, "BCNST2");

#ifdef HAVE_PGPLOT
  // Set the environment
  cpgsci(IUGUI_COLOR_LINE);
  cpglab(m_strXLabel, m_strYLabel, m_strTitle);
  cpgswin(dMinX, dMaxX, dMinY, dMaxY);
  cpgbox(strBoxX, 0.0, 0, strBoxY, 0.0, 0);

  cpgpt(m_nPlotPoints, m_pfPlotX, m_pfPlotY, m_nSymbol);
  if(m_bDrawLine)
    cpgline(m_nPlotPoints, m_pfPlotX, m_pfPlotY);

  if(m_pfPlotXErr)
    cpgerrx(m_nPlotPoints, m_pfPlotXErr1, m_pfPlotXErr2, m_pfPlotY, 1);

  if(m_pfPlotYErr)
    cpgerry(m_nPlotPoints, m_pfPlotX, m_pfPlotYErr1, m_pfPlotYErr2, 1);

  if(m_pfPlotYFit) {
    cpgsci(IUGUI_COLOR_BUTTON);
    cpgpt(m_nPlotPoints, m_pfPlotX, m_pfPlotYFit, 0-m_nSymbol);
    cpgsci(IUGUI_COLOR_LINE);
  } // if m_pfPlotYFit

  cpgsvp(0.0,1.0,0.0,1.0);
  cpgswin(0,1,0,1);
#endif // HAVE_PGPLOT
} // Update

/* This function calls the signal handler, if it exists, for button clicks
 * */
void CIPlot::OnClicked(double dPosX, double dPosY) {
  if(m_pOnClicked)
    m_pOnClicked(m_pwParentWindow);
  if(m_pwParentWindow)
    m_pwParentWindow->SetChildrenNotActive();
} // OnClicked

/* This function doesn't do anything yet, but it should respond to key-strokes
 * */
void CIPlot::OnKey(char cKey) {
} // OnKey

float CIPlot::MinX() {
  float fMin = m_pfPlotX[0];
  for(int i=1; i<m_nPlotPoints; i++)
    if(fMin > m_pfPlotX[i])
      fMin = m_pfPlotX[i];
  return fMin;
} // MinX

float CIPlot::MaxX() {
  float fMax = m_pfPlotX[0];
  for(int i=1; i<m_nPlotPoints; i++)
    if(fMax < m_pfPlotX[i])
      fMax = m_pfPlotX[i];
  return fMax;
} // MaxX

float CIPlot::MinY() {
  float fMin = m_pfPlotY[0];
  for(int i=1; i<m_nPlotPoints; i++)
    if(fMin > m_pfPlotY[i])
      fMin = m_pfPlotY[i];
  return fMin;
} // MinY

float CIPlot::MaxY() {
  float fMax = m_pfPlotY[0];
  for(int i=1; i<m_nPlotPoints; i++)
    if(fMax < m_pfPlotY[i])
      fMax = m_pfPlotY[i];
  return fMax;
} // MaxY

CIContourPlot::CIContourPlot() {
  m_pOnClicked = NULL;
  m_dPosX = 0;
  m_dSizeX = 1;
  m_dActivePosX = 0;
  m_dActiveSizeX = 1;
  m_dPosY = 0;
  m_dSizeY = 1;
  m_dActivePosY = 0;
  m_dActiveSizeY = 1;
  strcpy(m_strTitle, "");
  strcpy(m_strXLabel, "");
  strcpy(m_strYLabel, "");

  m_nSymbol = 16;

  SetNull();
} // CIContourPlot


CIContourPlot::~CIContourPlot() {
  m_pOnClicked = NULL;
  strcpy(m_strTitle, "");
  strcpy(m_strXLabel, "");
  strcpy(m_strYLabel, "");

  DeletePlot();
} // ~CIContourPlot

/* This function connects the given function to the signal handler for the click
 * of a mouse button. When this connection is made, an arbitrary function can be
 * called when the mouse button is pressed */
void CIContourPlot::ConnectOnClicked(void (*pFunc)(CIWindow *pwThis)) {
  m_pOnClicked = pFunc;
} // Settext

void CIContourPlot::SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY) {
  m_dPosX = dPosX;
  m_dPosY = dPosY;
  m_dSizeX = dSizeX;
  m_dSizeY = dSizeY;

  m_dActivePosX = dPosX;
  m_dActivePosY = dPosY;
  m_dActiveSizeX = dSizeX;
  m_dActiveSizeY = dSizeY;
} // SetPos

void CIContourPlot::SetPlot(int nPlotPointsX, int nPlotPointsY, double *pdX, double *pdY, double **pdZ) {
  DeletePlot();
  m_nPlotPointsX = nPlotPointsX;
  m_nPlotPointsY = nPlotPointsY;

  m_pfPlotX = new float[m_nPlotPointsX];
  m_pfPlotY = new float[m_nPlotPointsY];
  m_pfPlotZ = new float[m_nPlotPointsX*m_nPlotPointsY];

  for(int i=0; i<m_nPlotPointsX; i++) {
    m_pfPlotX[i] = float(pdX[i]);
    for(int j=0; j<m_nPlotPointsY; j++) {
      m_pfPlotZ[i + m_nPlotPointsX*j] = pdZ[i][j];
    } // for j
  } // for i

  for(int i=0; i<m_nPlotPointsY; i++) {
    m_pfPlotY[i] = float(pdY[i]);
  } // for i
} // SetPlot

void CIContourPlot::SetPlot(int nPlotPointsX, int nPlotPointsY, double *pdX, double *pdY, double *pdZ) {
  DeletePlot();
  m_nPlotPointsX = nPlotPointsX;
  m_nPlotPointsY = nPlotPointsY;

  m_pfPlotX = new float[m_nPlotPointsX];
  m_pfPlotY = new float[m_nPlotPointsY];
  m_pfPlotZ = new float[m_nPlotPointsX*m_nPlotPointsY];

  for(int i=0; i<m_nPlotPointsX; i++) {
    m_pfPlotX[i] = float(pdX[i]);
  } // for i

  for(int i=0; i<m_nPlotPointsX*m_nPlotPointsY; i++) {
    m_pfPlotZ[i] = pdZ[i];
  } // for j

  for(int i=0; i<m_nPlotPointsY; i++) {
    m_pfPlotY[i] = float(pdY[i]);
  } // for i
} // SetPlot

void CIContourPlot::SetLevels(int nLevels, double *pdLevels) {
  m_nLevels = nLevels;

  m_pfLevels = new float[m_nLevels];

  for(int i=0; i<m_nLevels; i++)
    m_pfLevels[i] = float(pdLevels[i]);
} // SetLevels

void CIContourPlot::SetTitle(const char strText[]) {
  strcpy(m_strTitle, strText);
} // SetTitle

void CIContourPlot::SetXLabel(const char strText[]) {
  strcpy(m_strXLabel, strText);
} // SetXLabel

void CIContourPlot::SetYLabel(const char strText[]) {
  strcpy(m_strYLabel, strText);
} // SetYLabel

void CIContourPlot::DeletePlot() {
  if(m_pfPlotX)
    delete[] m_pfPlotX;
  if(m_pfPlotY)
    delete[] m_pfPlotY;
  if(m_pfPlotZ)
    delete[] m_pfPlotZ;
  if(m_pfLevels)
    delete[] m_pfLevels;

  SetNull();
} // DeletePlot

/* This function saves the plot to a user-defined file or device */
void CIContourPlot::SavePlot(const char strDeviceName[]) {
  char strDeviceFile[160];
  m_pwParentWindow->GetDevice(strDeviceFile);
  m_pwParentWindow->DeInitialise();
#ifdef HAVE_PGPLOT
  cpgbeg(0, strDeviceName, 1, 1);
#endif // HAVE_PGPLOT

  m_bRespectAbsolutePosition = false;
  Update();
  m_bRespectAbsolutePosition = true;

#ifdef HAVE_PGPLOT
  cpgend();
#endif // HAVE_PGPLOT
  m_pwParentWindow->Initialise(strDeviceFile);
  return;
} // SavePlot


void CIContourPlot::SetNull() {
  m_nPlotPointsX = 0;
  m_nPlotPointsY = 0;
  m_nLevels = 0;
  m_pfPlotX = NULL;
  m_pfPlotY = NULL;
  m_pfPlotZ = NULL;
  m_pfLevels = NULL;
} // SetNull

/* This function draws the plot to the screen, and puts some text on it
 * */
void CIContourPlot::Update() {
  float pfTr[6];
  float fTemp;
  pfTr[0] = MinX(); pfTr[1] = (MaxX()-MinX())/m_nPlotPointsX; pfTr[2] = 0.0;
  pfTr[3] = MinY(); pfTr[4] = 0.0; pfTr[5] = (MaxY()-MinY())/m_nPlotPointsY;
  if(! m_pfPlotX || ! m_pfPlotY || ! m_pfPlotZ || ! m_pfLevels)
    return;

  // Set the viewport and window
#ifdef HAVE_PGPLOT
  if(m_bRespectAbsolutePosition)
    cpgsvp(m_dPosX, m_dPosX+m_dSizeX, m_dPosY, m_dPosY+m_dSizeY);
  else
    cpgsvp(0.1, 0.9, 0.1, 0.9);

  // Set the environment
  cpgsci(IUGUI_COLOR_LINE);
  cpgswin(MinX(), MaxX(), MinY(), MaxY());
  cpgbox("BCNST", 0.0, 0, "BCNST", 0.0, 0);
  cpglab(m_strXLabel, m_strYLabel, m_strTitle);

  // Make the contour plot
  cpgbbuf();
  cpgsci(IUGUI_COLOR_SIGMA1);
  cpgcont(m_pfPlotZ, m_nPlotPointsX, m_nPlotPointsY, 1, m_nPlotPointsX, 1, m_nPlotPointsY, m_pfLevels+0, -1, pfTr);
  cpgsci(IUGUI_COLOR_SIGMA2);
  cpgcont(m_pfPlotZ, m_nPlotPointsX, m_nPlotPointsY, 1, m_nPlotPointsX, 1, m_nPlotPointsY, m_pfLevels+1, -1, pfTr);
  cpgsci(IUGUI_COLOR_SIGMA3);
  cpgcont(m_pfPlotZ, m_nPlotPointsX, m_nPlotPointsY, 1, m_nPlotPointsX, 1, m_nPlotPointsY, m_pfLevels+2, -1, pfTr);
  cpgebuf();

  cpgsvp(0.0,1.0,0.0,1.0);
  cpgswin(0,1,0,1);

  cpgsci(IUGUI_COLOR_LINE);
#endif // HAVE_PGPLOT
} // Update

/* This function calls the signal handler, if it exists, for button clicks
 * */
void CIContourPlot::OnClicked(double dPosX, double dPosY) {
  if(m_pOnClicked)
    m_pOnClicked(m_pwParentWindow);
  if(m_pwParentWindow)
    m_pwParentWindow->SetChildrenNotActive();
} // OnClicked

/* This function doesn't do anything yet, but it should respond to key-strokes
 * */
void CIContourPlot::OnKey(char cKey) {
} // OnKey

float CIContourPlot::MinX() {
  float fMin = m_pfPlotX[0];
  for(int i=1; i<m_nPlotPointsX; i++)
    if(fMin > m_pfPlotX[i])
      fMin = m_pfPlotX[i];
  return fMin;
} // MinX

float CIContourPlot::MaxX() {
  float fMax = m_pfPlotX[0];
  for(int i=1; i<m_nPlotPointsX; i++)
    if(fMax < m_pfPlotX[i])
      fMax = m_pfPlotX[i];
  return fMax;
} // MaxX

float CIContourPlot::MinY() {
  float fMin = m_pfPlotY[0];
  for(int i=1; i<m_nPlotPointsY; i++)
    if(fMin > m_pfPlotY[i])
      fMin = m_pfPlotY[i];
  return fMin;
} // MinY

float CIContourPlot::MaxY() {
  float fMax = m_pfPlotY[0];
  for(int i=1; i<m_nPlotPointsY; i++)
    if(fMax < m_pfPlotY[i])
      fMax = m_pfPlotY[i];
  return fMax;
} // MaxY

CIText::CIText() {
  strcpy(m_strText, "");

  m_fFontSize = 1.0; //m_pwParentWindow->GetFontSize();
} // CIText


CIText::~CIText() {
  strcpy(m_strText, "");
} // ~CIText


void CIText::SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY) {
  m_dPosX = dPosX;
  m_dPosY = dPosY;
  m_dSizeX = dSizeX;
  m_dSizeY = dSizeY;

  m_dActivePosX = dPosX;
  m_dActivePosY = dPosY;
  m_dActiveSizeX = dSizeX;
  m_dActiveSizeY = dSizeY;
} // SetPos

void CIText::SetText(const char strText[]) {
  strcpy(m_strText, strText);
} // Settext

void CIText::SetFontSize(float fFontSize) {
  m_fFontSize = fFontSize;
} // SetFontSize

/* This function draws the button to the screen, and puts some text on it
 * */
void CIText::Update() {
  float fParentFontSize = m_pwParentWindow->GetFontSize();

#ifdef HAVE_PGPLOT
  cpgsci(IUGUI_COLOR_TEXT);
  cpgsch(m_fFontSize); 
  if(strlen(m_strText) > 0)
    cpgtext(m_dPosX, m_dPosY+0.01, m_strText);
  cpgsch(fParentFontSize); 
  cpgsci(IUGUI_COLOR_BKGRD);
#endif // HAVE_PGPLOT
} // Update

