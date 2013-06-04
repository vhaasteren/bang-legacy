/* iugui.h -- Interactive universal graphical user interface
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


#ifndef __IUGUI_H__
#define __IUGUI_H__

#define IUGUI_COLOR_BKGRD	0
#define IUGUI_COLOR_LINE	1
#define IUGUI_COLOR_BUTTON	2
#define IUGUI_COLOR_BUTTONTEXT	3
#define IUGUI_COLOR_ACTIVE	4
#define IUGUI_COLOR_TEXT	5
#define IUGUI_COLOR_SIGMA1	6
#define IUGUI_COLOR_SIGMA2	7
#define IUGUI_COLOR_SIGMA3	8

#define MAX_TABS		25
#define MAX_COMBOITEMS		150

#include "config.h"
#include <stdio.h>

// We will define these classes in this header
class CIWindow;
class CIWidget;
class CIWidgetList;
class CIWidgetListElement;

// A linked-list of widgets (for children)
class CIWidgetList {
public:
  // Constructor and destructor
  CIWidgetList();
  ~CIWidgetList();

  void AddWidget(CIWidget *pwWidget);
  void RemoveWidget(CIWidget *pwWidget);

  int GetSize();
  bool HasWidget(CIWidget *pwWidget);

  // Member functions to retrieve members of the list
  CIWidget *GetFirst();
  CIWidget *GetLast();

  CIWidget *GetNext(CIWidget *pwWidget);
  CIWidget *GetPrevious(CIWidget *pwWidget);

  // These are probably not really needed. Don't waste time yet
//  void InsertBefore(CIWidget *pwWidgetElement, CIWidget *pwWidgetAddition);
//  void InsertAfter(CIWidget *pwWidgetElement, CIWidget *pwWidgetAddition);
private:
  // In principle this is not necessary, but keep track of it anayway
  int m_nSize;

  CIWidgetListElement *m_peFirst;
  CIWidgetListElement *m_peLast;
protected:
  CIWidgetListElement *GetListElement(CIWidget *pwWidget);
} ;

// A linked-list element (basically the widget itself)
class CIWidgetListElement {
public:
  // Constructor and destructor
  CIWidgetListElement(CIWidget *pwWidget=NULL);
  ~CIWidgetListElement();

  CIWidget *GetWidget();

  CIWidgetListElement *GetNext();
  CIWidgetListElement *GetPrevious();

  CIWidgetListElement *SetNext(CIWidgetListElement *peNext);
  CIWidgetListElement *SetPrevious(CIWidgetListElement *pePrevious);
private:
  CIWidget *m_pwWidget;

  CIWidgetListElement *m_peNext;
  CIWidgetListElement *m_pePrevious;
protected:
} ;

// The main-window class
class CIWindow {
public:
  // Constructor / Destructor
  CIWindow();
  ~CIWindow();

  // Initialise the window
  void Initialise(const char strDevice[]);
  void DeInitialise();

  // Process all the input and send signals
  void ProcessSignals();

  // Whether we need to update the screen
  bool NeedsUpdate();

  // Update the screen
  void Update();

  // Add widgets
  void AddWidget(CIWidget *pwWidget);
  CIWidget *GetActiveWidget();
  CIWidget *GetWidget(double dPosX, double dPosY);
  CIWidget *GetNonActiveWidget(double dPosX, double dPosY);
  CIWidget *GetWidgetInRange(double dPosX, double dPosY, double dSizeX, double dSizeY);

  bool HasWidget(CIWidget *pwWidget);
  void RemoveWidget(CIWidget *pwWidget);

  // Set all widgets to non-active
  void SetChildrenNotActive();

  // Check whether the window is done
  bool Done();
  void SetDone(bool bDone=true);
  void ScreenShot(const char strDevice[]);
  void GetDevice(char *strDevice);

  float GetFontSize();
protected:
  //Signal handlers:
  virtual void OnClicked(double dPosX, double dPosY);
  virtual void OnKey(char cKey);

  // Children widgets:
  CIWidgetList m_wlChildren;

  // Member variables
  char m_strDeviceFile[160];
  double m_dSizeX;
  double m_dSizeY;
  bool m_bNeedsUpdate;
  bool m_bInitialised;
  bool m_bDone;

  // Font variables
  float m_fAspectRatio;		// Aspect ratio. 1 = square
  float m_fWidth;		// Width of view in inches. 0=big
  float m_fFontSize;		// Font size (vertical). 1 = normal
  int m_nFontType;		// Normal font (1 - 4)
  int m_nLineWidth;		// Normal line width (1 - 201)

  // Member variable to prevent tempo2/PGPLOT bug
  bool m_bAlreadyPapped;
};

// The widget class
class CIWidget {
public:
  // Constructor / Destructor
  CIWidget(double dSizeX=-1, double dSizeY=-1);
  virtual ~CIWidget();

  // Initialise the widget
  void Initialise();
  void DeviceInit(const char strDevice[]);
  void DeviceDeInit();

  //Signal handlers:
  virtual void OnClicked(double dPosX, double dPosY);
  virtual void OnKey(char cKey);

  // Position stuff
  double GetPosX();
  double GetPosY();
  double GetSizeX();
  double GetSizeY();
  double GetActivePosX();
  double GetActivePosY();
  double GetActiveSizeX();
  double GetActiveSizeY();

  void SetActive(bool bActive=true);
  bool Active();
  void SetChildrenNotActive();
  virtual void SetVisible(bool bVisible=true);
  bool Visible();

  // Checks for range of widgets
  bool InRange(double dPosX, double dPosY);
  bool InActiveRange(double dPosX, double dPosY);

  // Add widgets
  virtual void AddWidget(CIWidget *pwWidget);
  CIWidget *GetActiveWidget();
  CIWidget *GetWidget(double dPosX, double dPosY);
  CIWidget *GetNonActiveWidget(double dPosX, double dPosY);
  CIWidget *GetWidgetInRange(double dPosX, double dPosY, double dSizeX, double dSizeY);

  bool HasWidget(CIWidget *pwWidget);
  void RemoveWidget(CIWidget *pwWidget);

  void SetParent(CIWindow *pwWindow);
  void SetParent(CIWidget *pwWidget);

  // Update the screen
  virtual void Update();
  void UpdateChildren();
protected:
  // The position in the parent widget/window
  double m_dPosX;
  double m_dPosY;

  // The size of this widget
  double m_dSizeX;
  double m_dSizeY;

  // The active dimensions of this widget
  double m_dActivePosX;
  double m_dActivePosY;
  double m_dActiveSizeX;
  double m_dActiveSizeY;

  // Whether this widget is visible/active
  bool m_bVisible;
  bool m_bActive;

  // Whether we use whole device or not (not fully implemented)
  bool m_bRespectAbsolutePosition;

  // The parents of this widget
  CIWindow *m_pwParentWindow;
  CIWidget *m_pwParentWidget;

  // Children widgets:
  CIWidgetList m_wlChildren;
private:
};

// The button class
class CIButton : public CIWidget {
public:
//  CIButton(double dSizeX=-1, double dSizeY=-1);
  CIButton();
  ~CIButton();

  void SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY);
  void SetText(const char strText[]);
  void ConnectOnClicked(void (*pFunc)(CIWindow *pwThis));

  void Update();
protected:
  //Signal handlers:
  void OnClicked(double dPosX, double dPosY);
  void OnKey(char cKey);

  // Call-back function
  void (*m_pOnClicked)(CIWindow *pwThis);

  // Member variables
  char m_strText[80];
};

// The group widget
class CIGroup : public CIWidget {
public:
  CIGroup();
  ~CIGroup();

  void SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY);
  void SetBorder(bool bShowBorder=true);
  bool GetBorder();

  void SetText(const char strText[]);
  void SetTextWidth(double dTextWidth);

  void AddWidget(CIWidget *pwWidget);

  void Update();
protected:
  //Signal handlers:
  void OnClicked(double dPosX, double dPosY);
  void OnKey(char cKey);

  bool m_bDrawBorder;
  char m_strText[80];
  double m_dTextWidth;
  double m_dTextPadding;
};


// The tab class
class CITab : public CIWidget {
public:
  CITab();
  ~CITab();

  void SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY);

  void SetTabs(int nTabs=1);
  int GetTabs();
  void SetText(int nTab, const char strText[]);
  void SetActiveTab(int nTab);
  int GetActiveTab();
  void SetVisible(bool bVisible=true);

  void SetMenuItemWidth(double dMenuItemWidth);

  void AddWidget(int nTab, CIWidget *pwWidget);
  void RemoveWidget(CIWidget *pwWidget);

  void Update();
protected:
  //Signal handlers:
  void OnClicked(double dPosX, double dPosY);
  void OnKey(char cKey);

  // Member variables
  double m_dMenuHeight;
  double m_dMenuItemWidth;
  double m_dMenuItemWidthSpace;
  int m_nTabs;
  int m_nActiveTab;
  char m_strTabText[MAX_TABS][80];
  CIWidgetList m_wlTabChildren[MAX_TABS];
};

// The combobox class
class CIComboBox : public CIWidget {
public:
//  CIButton(double dSizeX=-1, double dSizeY=-1);
  CIComboBox();
  ~CIComboBox();

  void SetPos(double dPosX, double dPosY, double dWidth=-1);
  void SetItems(int nItems);
  void SetItem(int nItem, char strItem[]);

  void SetCurrentItem(int nItem);
  int GetCurrentItem();
  void CopyItemText(int nItem, char *strText);

  // The user just selected a new item
  void ConnectOnSelected(void (*pFunc)(CIWindow *pwThis));

  void Update();
protected:
  //Signal handlers:
  void OnClicked(double dPosX, double dPosY);
  void OnKey(char cKey);

  // Call-back function
  void (*m_pOnSelected)(CIWindow *pwThis);

  // Member variables
  int m_nItems;
  int m_nCurrentItem;
  double m_dItemHeight;
  double m_dItemHeightSpace;

  char m_strItemText[MAX_COMBOITEMS][80];
};


// The button class
class CICheckbox : public CIWidget {
public:
//  CICheckbox(double dSizeX=-1, double dSizeY=-1);
  CICheckbox();
  ~CICheckbox();

  void SetPos(double dPosX, double dPosY);
  void SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY);
  void ConnectOnClicked(void (*pFunc)(CIWindow *pwThis));

  void SetChecked(bool bChecked=true);
  bool GetChecked();

  void Update();
protected:
  //Signal handlers:
  void OnClicked(double dPosX, double dPosY);
  void OnKey(char cKey);

  // Call-back function
  void (*m_pOnClicked)(CIWindow *pwThis);

  // Member variables
  bool m_bChecked;
  double m_dPadding;
  double m_dBorderSize;
};

/*
 * cpgeras();			Erase window
 * cpgsvp(0.0,1.0,0.8,1.0);	Set viewport
 * cpgswin(0,1,0,1);		Set world coordinates
 * */
// The CIPlot class (simple 1D plots)
class CIPlot : public CIWidget {
public:
  CIPlot();
  ~CIPlot();

  void SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY);

  void SetPlot(int nPlotPoints, double *pdX, double *pdY);
  void SetXErr(int nPlotPoints, double *pdXErr);
  void SetYErr(int nPlotPoints, double *pdYErr);
  void SetYFit(int nPlotPoints, double *pdYFit);
  void SetLogScale(int nAxis, bool bLog);
  void DrawLine(bool bLine);
  void UnsetYFit();

  void SetTitle(const char strText[]);
  void SetXLabel(const char strText[]);
  void SetYLabel(const char strText[]);

  void SetSymbol(int nSymbol);
  void SetUserView(double dMinX, double dMaxX, double dMinY, double dMaxY);
  void UnSetUserView();

  void DeletePlot();
  void SavePlot(const char strDeviceName[]);

  void ConnectOnClicked(void (*pFunc)(CIWindow *pwThis));

  void Update();
protected:
  void SetNull();

  //Signal handlers:
  void OnClicked(double dPosX, double dPosY);
  void OnKey(char cKey);

  float MinX();
  float MaxX();
  float MinY();
  float MaxY();

  // Call-back function
  void (*m_pOnClicked)(CIWindow *pwThis);

  // Member variables
  char m_strTitle[80];
  char m_strXLabel[80];
  char m_strYLabel[80];

  int m_nSymbol;
  bool m_bUserDefinedView;
  double m_dMinX;
  double m_dMaxX;
  double m_dMinY;
  double m_dMaxY;

  // Flags
  bool m_bLogX;
  bool m_bLogY;
  bool m_bDrawLine;

  // Plot stuff
  int m_nPlotPoints;
  float *m_pfPlotX;
  float *m_pfPlotY;
  float *m_pfPlotXErr;
  float *m_pfPlotXErr1;
  float *m_pfPlotXErr2;
  float *m_pfPlotYErr;
  float *m_pfPlotYErr1;
  float *m_pfPlotYErr2;

  float *m_pfPlotYFit;
};

/*
 * cpgeras();			Erase window
 * cpgsvp(0.0,1.0,0.8,1.0);	Set viewport
 * cpgswin(0,1,0,1);		Set world coordinates
 * */
// The CIContourPlot class (contour plots)
class CIContourPlot : public CIWidget {
public:
  CIContourPlot();
  ~CIContourPlot();

  void SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY);

  void SetPlot(int nPlotPointsX, int nPlotPointsY, double *pdX, double *pdY, double **pdZ);
  void SetPlot(int nPlotPointsX, int nPlotPointsY, double *pdX, double *pdY, double *pdZ);
  void SetLevels(int nLevels, double *pdLevels);

  void SetTitle(const char strText[]);
  void SetXLabel(const char strText[]);
  void SetYLabel(const char strText[]);

  void DeletePlot();
  void SavePlot(const char strDeviceName[]);

  void ConnectOnClicked(void (*pFunc)(CIWindow *pwThis));

  void Update();
protected:
  void SetNull();

  //Signal handlers:
  void OnClicked(double dPosX, double dPosY);
  void OnKey(char cKey);

  float MinX();
  float MaxX();
  float MinY();
  float MaxY();

  // Call-back function
  void (*m_pOnClicked)(CIWindow *pwThis);

  // Member variables
  char m_strTitle[80];
  char m_strXLabel[80];
  char m_strYLabel[80];
  int m_nSymbol;

  // Plot stuff
  int m_nPlotPointsX;
  int m_nPlotPointsY;
  int m_nLevels;
  float *m_pfPlotX;
  float *m_pfPlotY;
  float *m_pfPlotZ;
  float *m_pfLevels;
};


// The text class
class CIText : public CIWidget {
public:
  CIText();
  ~CIText();

  void SetPos(double dPosX, double dSizeX, double dPosY, double dSizeY);
  void SetText(const char strText[]);
  void SetFontSize(float fFontSize);
//  void ConnectOnClicked(void (*pFunc)(CIWindow *pwThis));

  void Update();
protected:
  //Signal handlers:
//  void OnClicked(double dPosX, double dPosY);
//  void OnKey(char cKey);

  // Call-back function
//  void (*m_pOnClicked)(CIWindow *pwThis);

  // Member variables
  char m_strText[80];

  float m_fFontSize;
};


#endif // __IUGUI_H__
