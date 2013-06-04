/* bangui.cpp -- iugui interface
 *
 * Rutger van Haasteren 23 December 2008 haasteren@strw.leidenuniv.nl
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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_spline.h>

#include "linal.h"
#include "linalfunc.h"
#include "corefunctions.h"
#include "analyze.h"
#include "banfunc.h"
//#include "tempo2_compat.h"

#include "bangui.h"
#include "filefunctions.h"


//BanWindow::BanWindow(int argc, char *argv[], bool bPlugin) {
BanWindow::BanWindow(int argc, char *argv[], bool bPlugin, const char *strParametersConf) {
  // Remember whether we are in a plugin or not
  m_bPlugin = bPlugin;

  // Set Random Number generator, read tabulas, etc.
  BANFuncInitialize(m_nProcessRank);

  // Set several default values for the member variables
  SetDefaultValues(bPlugin);

  // If we have a parameters-file, add it
  if(strParametersConf && strlen(strParametersConf) > 0) {
    strcpy(m_strParametersConf, strParametersConf);
  } // if strParametersConf


  // Some constants that are not used (fill anyway)
  strcpy(m_oConstants.strConfigVersion, "0.60");
  strcpy(m_oConstants.strVersion, "0.60");
  strcpy(m_oConstants.strDataDir, "notused");
  strcpy(m_oConstants.strDataDir, "residuals.dat");
  strcpy(m_oConstants.strTempoDataDir, "notused");
  strcpy(m_oConstants.strAnglesFile, "notused");
  m_oConstants.bUnevenlySampled = true;
  m_oConstants.nFitParameters = 0;
  m_oConstants.nParameters = 0;
  m_oConstants.dMeasureTimeInterval = 0;

  // Set all things to zero
  m_oConstants.k = 0;
  m_oConstants.l = 0;
  m_oConstants.m = 0;
  m_oConstants.nSources = 0;
  m_oConstants.poPulsars = NULL;
  m_oConstants.poSources = NULL;

  // TODO: set this one based on the presence of reduced basis data!
  m_oConstants.bUseReducedBasis = false;

  // TODO: data compression stuff should be defined in the HDF5 file
  m_oConstants.dReducedFidelity = 0.999;
  m_oConstants.bCorAmpAccel = false;
  m_oConstants.bCalcTMPars = false;
  m_oConstants.bUseInstrTags = true;  // NANOGrav hack, always true!
  m_oConstants.nInterpolationBins = 118;

  // The interpolation stuff (global)
#if 0
  m_oData.bGWBInterpolatorsSet = false;
  m_oData.ppGWBIntAccel = NULL;
  m_oData.ppGWBIntSpline = NULL;
  m_oData.bClockInterpolatorsSet = false;
  m_oData.ppClockIntAccel = NULL;
  m_oData.ppClockIntSpline = NULL;

  m_oData.nRedLikGWBPar = -1;;
  m_oData.nRedLikAllPar = -1;;
  m_oData.nRedLikClockPar = -1;;
#endif

  // Setting the static-variable replacements
  m_oData.bRedLikFirstRun = true;
} // BanWindow


BanWindow::~BanWindow() {
  if(m_oConstants.poPulsars) {
    for(int a=0; a<m_oConstants.k; ++a) {
      delete[] m_oConstants.poPulsars[a].pdTOA;
      delete[] m_oConstants.poPulsars[a].pdResiduals;
      delete[] m_oConstants.poPulsars[a].pdDeltaResiduals;
      delete[] m_oConstants.poPulsars[a].pdFreq;
      delete[] m_oConstants.poPulsars[a].pbFlagSet;

      for(int i=0; i<m_oConstants.poPulsars[a].nObservations; ++i) {
	delete[] m_oConstants.poPulsars[a].pstrFlags[i];
      } // for i
      delete[] m_oConstants.poPulsars[a].pstrFlags;
    } // for a

    delete[] m_oConstants.poPulsars;
    m_oConstants.poPulsars = NULL;
    m_oConstants.k = 0;
  } // if poPulsars

  if(m_oConstants.poSources) {
    for(int s=0; s<m_oConstants.nSources; ++s) {
      if(m_oConstants.poSources[s].pbScope) {
	delete[] m_oConstants.poSources[s].pbScope;
	m_oConstants.poSources[s].pbScope = NULL;
      } // if pbScope

      if(m_oConstants.poSources[s].oIntAcc.bSet) {
	for(int i=0; i<m_oConstants.poSources[s].oIntAcc.nUniques; ++i) {
	  gsl_spline_free(m_oConstants.poSources[s].oIntAcc.ppIntSpline[i]);
	  gsl_interp_accel_free(m_oConstants.poSources[s].oIntAcc.ppIntAccel[i]);
	} // for i
	delete[] m_oConstants.poSources[s].oIntAcc.ppIntSpline;
	delete[] m_oConstants.poSources[s].oIntAcc.ppIntAccel;
	m_oConstants.poSources[s].oIntAcc.ppIntSpline = NULL;
	m_oConstants.poSources[s].oIntAcc.ppIntAccel = NULL;
	m_oConstants.poSources[s].oIntAcc.bSet = false;
	m_oConstants.poSources[s].oIntAcc.nUniques = 0;
      } // if bSet
    } // for s

    delete[] m_oConstants.poSources;
    m_oConstants.poSources = NULL;
    m_oConstants.nSources = 0;
  } // if poPulsars
} // ~BanWindow


void BanWindow::InitializeWidgets() {
  int nPulsarStochasticSources, nPulsarDeterministicSources, nPulsar, nIndex,
      nType;

  // The main tab (the content of this widget is fixed)
  m_MainTab.SetPos(0.00, 1.00, 0.05, 0.95);
  m_MainTab.SetMenuItemWidth(0.11);
  m_MainTab.SetTabs(6);
  AddWidget(&m_MainTab);

  // Add button quit to main window
  m_ButtonQuit.SetPos(0.92, 0.06, 0.00, 0.028);
  m_ButtonQuit.ConnectOnClicked(BanWindow::OnClickedQuit);
  AddWidget(&m_ButtonQuit);

  // Add button load to main window
  m_ButtonLoad.SetPos(0.8, 0.10, 0.00, 0.028);
  m_ButtonLoad.ConnectOnClicked(BanWindow::OnClickedLoad);
  AddWidget(&m_ButtonLoad);

  // Add button load to main window
  m_ButtonPrint.SetPos(0.7, 0.08, 0.00, 0.028);
  m_ButtonPrint.ConnectOnClicked(BanWindow::OnClickedPrint);
  AddWidget(&m_ButtonPrint);

  m_Text_Copyright.SetPos(0.01, 0.08, 0.00, 0.028);
  m_Text_Copyright.SetFontSize(0.6);
  AddWidget(&m_Text_Copyright);

  InitializeWidgetsTab0();
  InitializeWidgetsTab1();
  InitializeWidgetsTab2();
  InitializeWidgetsTab3();
  InitializeWidgetsTab4();
  InitializeWidgetsTab5();
} // InitializeWidgets

/* This function Initialises all Widgets in Tab 0
 * */
void BanWindow::InitializeWidgetsTab0() {
  // Add pulsar combobox to tab 0
  m_Tab0_Combo_Pulsars.SetPos(0.13, 0.90, 0.15);
  m_Tab0_Combo_Pulsars.SetItems(m_oConstants.k);
  m_Tab0_Combo_Pulsars.ConnectOnSelected(BanWindow::OnPulsarTab0Selected);
  m_MainTab.AddWidget(0, &m_Tab0_Combo_Pulsars);

  // Add texts to tab 0
  m_Tab0_Text_Pulsars.SetPos(0.03, 0.07, 0.90, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_Pulsars);

  m_Tab0_Text_Datapoints_Description.SetPos(0.03, 0.09, 0.84, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_Datapoints_Description);
  m_Tab0_Text_Datapoints_Value.SetPos(0.13, 0.10, 0.84, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_Datapoints_Value);

  m_Tab0_Text_Timespan_Description.SetPos(0.03, 0.09, 0.80, 0.028);
  m_Tab0_Text_Timespan_Description.SetText("Timespan:");
  m_MainTab.AddWidget(0, &m_Tab0_Text_Timespan_Description);
  m_Tab0_Text_Timespan_Value.SetPos(0.13, 0.10, 0.80, 0.028);
  m_Tab0_Text_Timespan_Value.SetText("");
  m_MainTab.AddWidget(0, &m_Tab0_Text_Timespan_Value);

  m_Tab0_Text_RMS_Description.SetPos(0.03, 0.09, 0.76, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_RMS_Description);
  m_Tab0_Text_RMS_Value.SetPos(0.13, 0.10, 0.76, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_RMS_Value);

  m_Tab0_Text_Chisq_Description.SetPos(0.03, 0.09, 0.72, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_Chisq_Description);
  m_Tab0_Text_Chisq_Value.SetPos(0.13, 0.10, 0.72, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_Chisq_Value);

  m_Tab0_Text_RedChisq_Description.SetPos(0.03, 0.09, 0.68, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_RedChisq_Description);
  m_Tab0_Text_RedChisq_Value.SetPos(0.13, 0.10, 0.68, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_RedChisq_Value);

  m_Tab0_Button_MemMap.SetPos(0.03, 0.08, 0.26, 0.028);
  m_Tab0_Button_MemMap.ConnectOnClicked(BanWindow::OnClickedTab0MemMap);
//  m_MainTab.AddWidget(0, &m_Tab0_Button_MemMap);
  m_Tab0_Text_MemMap_Description.SetPos(0.13, 0.10, 0.26, 0.028);
//  m_MainTab.AddWidget(0, &m_Tab0_Text_MemMap_Description);

  m_Tab0_Button_SavePlot.SetPos(0.03, 0.08, 0.21, 0.028);
  m_Tab0_Button_SavePlot.ConnectOnClicked(BanWindow::OnClickedTab0SavePlot);
  m_MainTab.AddWidget(0, &m_Tab0_Button_SavePlot);

  m_Tab0_Button_Fit.SetPos(0.03, 0.08, 0.16, 0.028);
  m_Tab0_Button_Fit.ConnectOnClicked(BanWindow::OnFitTab0Selected);
  m_MainTab.AddWidget(0, &m_Tab0_Button_Fit);
  m_Tab0_Text_Fit_Description.SetPos(0.13, 0.10, 0.16, 0.028);
  m_MainTab.AddWidget(0, &m_Tab0_Text_Fit_Description);

  m_Tab0_Button_Plk.SetPos(0.03, 0.08, 0.11, 0.028);
  m_Tab0_Button_Plk.ConnectOnClicked(BanWindow::OnSelectedTab0Plk);
  if(m_bPlugin)
    m_MainTab.AddWidget(0, &m_Tab0_Button_Plk);

  m_Tab0_Text_Plk_Description.SetPos(0.13, 0.10, 0.11, 0.028);
  if(m_bPlugin)
    m_MainTab.AddWidget(0, &m_Tab0_Text_Plk_Description);

  // Add residuals plot to tab 0
  m_Tab0_Plot_Residuals.SetPos(0.33, 0.62, 0.18, 0.68);
  m_MainTab.AddWidget(0, &m_Tab0_Plot_Residuals);
} // InitializeWidgetsTab0

/* This function Initialises all Widgets in Tab 1
 * */
void BanWindow::InitializeWidgetsTab1() {
  int nPulsarStochasticSources, nPulsarDeterministicSources, nPulsar, nIndex, nType;

  // Check how many sources belong to the active pulsar
  nPulsarStochasticSources = 0;
  nPulsarDeterministicSources = 0;
  nIndex = 0;
  nPulsar = m_Tab0_Combo_Pulsars.GetCurrentItem();
  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      if(m_oConstants.poSources[s].oSourceType.eID == SID_Deterministic)
	nPulsarDeterministicSources++;
      else
  	nPulsarStochasticSources++;
    } // if SourceWorksOnPulsar
  } // for s

  m_Tab1_Combo_Pulsars.SetPos(0.03, 0.90, 0.15);
  m_Tab1_Combo_Pulsars.SetItems(m_oConstants.k);
  m_Tab1_Combo_Pulsars.ConnectOnSelected(BanWindow::OnPulsarTab1Selected);
  if(! m_MainTab.HasWidget(&m_Tab1_Combo_Pulsars))
    m_MainTab.AddWidget(1, &m_Tab1_Combo_Pulsars);

  m_Tab1_Text_Stochastic.SetPos(0.03, 0.24, 0.82, 0.028);
  if(! m_MainTab.HasWidget(&m_Tab1_Text_Stochastic))
    m_MainTab.AddWidget(1, &m_Tab1_Text_Stochastic);

  m_Tab1_StochasticTab.SetPos(0.03, 0.44, 0.4, 0.40);
  m_Tab1_StochasticTab.SetTabs(nPulsarStochasticSources);
  m_Tab1_StochasticTab.SetMenuItemWidth(0.04);
  if(! m_MainTab.HasWidget(&m_Tab1_StochasticTab))
    m_MainTab.AddWidget(1, &m_Tab1_StochasticTab);

  m_Tab1_Text_Deterministic.SetPos(0.53, 0.24, 0.82, 0.028);
  if(! m_MainTab.HasWidget(&m_Tab1_Text_Deterministic))
    m_MainTab.AddWidget(1, &m_Tab1_Text_Deterministic);

  m_Tab1_DeterministicTab.SetPos(0.53, 0.44, 0.4, 0.40);
  m_Tab1_DeterministicTab.SetTabs(nPulsarDeterministicSources);
  m_Tab1_DeterministicTab.SetMenuItemWidth(0.04);
  if(! m_MainTab.HasWidget(&m_Tab1_DeterministicTab))
    m_MainTab.AddWidget(1, &m_Tab1_DeterministicTab);

  // Calculation widgets
  m_Tab1_Button_PowellOptimum.SetPos(0.03, 0.14, 0.32, 0.028);
  m_Tab1_Button_PowellOptimum.ConnectOnClicked(BanWindow::OnClickedTab1PowellOptimum);
  m_MainTab.AddWidget(1, &m_Tab1_Button_PowellOptimum);

  SetModelWidgets();
} // InitializeWidgetsTab1

/* This function Initialises all Widgets in Tab 2
 * */
void BanWindow::InitializeWidgetsTab2() {
  // Add button RunMCMC to Tab 2
  m_Tab2_Button_RunMCMC.SetPos(0.04, 0.10, 0.91, 0.028);
  m_Tab2_Button_RunMCMC.ConnectOnClicked(BanWindow::OnClickedRunMCMC);
  m_MainTab.AddWidget(2, &m_Tab2_Button_RunMCMC);

  // Add Group 1D Integrate to Tab 2
  m_Tab2_Group_1D.SetPos(0.03, 0.25, 0.63, 0.24);
  m_Tab2_Group_1D.SetTextWidth(12*0.0105-0.007);
  m_MainTab.AddWidget(2, &m_Tab2_Group_1D);

  // Add elements of Tab 2 group 1D
  m_Tab2_Text_Pulsars1D_Description.SetPos(0.04, 0.08, 0.82, 0.028);
  m_Tab2_Group_1D.AddWidget(&m_Tab2_Text_Pulsars1D_Description);
  m_Tab2_Combo_Pulsars1D.SetPos(0.12, 0.82, 0.15);
  m_Tab2_Combo_Pulsars1D.SetItems(0);
  m_Tab2_Combo_Pulsars1D.ConnectOnSelected(BanWindow::OnSelectedTab2Pulsar1D);
  m_Tab2_Group_1D.AddWidget(&m_Tab2_Combo_Pulsars1D);

  m_Tab2_Text_Parameters1D_Description.SetPos(0.04, 0.08, 0.77, 0.028);
  m_Tab2_Group_1D.AddWidget(&m_Tab2_Text_Parameters1D_Description);
  m_Tab2_Combo_Parameter1D.SetPos(0.12, 0.77, 0.15);
  m_Tab2_Combo_Parameter1D.SetItems(0);
  m_Tab2_Group_1D.AddWidget(&m_Tab2_Combo_Parameter1D);

  m_Tab2_Text_CalcErr_Description.SetPos(0.04, 0.08, 0.72, 0.028);
  m_Tab2_Group_1D.AddWidget(&m_Tab2_Text_CalcErr_Description);
  m_Tab2_Check_CalcErr.SetPos(0.16, 0.72);
  m_Tab2_Group_1D.AddWidget(&m_Tab2_Check_CalcErr);

  m_Tab2_Button_Integrate1D.SetPos(0.04, 0.12, 0.65, 0.028);
  m_Tab2_Button_Integrate1D.ConnectOnClicked(BanWindow::OnClickedIntegrate1D);
  m_Tab2_Group_1D.AddWidget(&m_Tab2_Button_Integrate1D);


  // Add Group 2D Integrate to Tab 2
  m_Tab2_Group_2D.SetPos(0.03, 0.25, 0.31, 0.25);
  m_Tab2_Group_2D.SetTextWidth(12*0.0105-0.007);
  m_MainTab.AddWidget(2, &m_Tab2_Group_2D);

  // Add elements of Tab 2 group 2D
  m_Tab2_Text_Pulsars2D_Description.SetPos(0.04, 0.08, 0.51, 0.028);
  m_Tab2_Group_2D.AddWidget(&m_Tab2_Text_Pulsars2D_Description);
  m_Tab2_Combo_Pulsars2D.SetPos(0.12, 0.51, 0.15);
  m_Tab2_Combo_Pulsars2D.SetItems(0);
  m_Tab2_Combo_Pulsars2D.ConnectOnSelected(BanWindow::OnSelectedTab2Pulsar2D);
  m_Tab2_Group_2D.AddWidget(&m_Tab2_Combo_Pulsars2D);

  m_Tab2_Text_Parameters2D_1_Description.SetPos(0.04, 0.08, 0.46, 0.028);
  m_Tab2_Group_2D.AddWidget(&m_Tab2_Text_Parameters2D_1_Description);
  m_Tab2_Combo_Parameter2D_1.SetPos(0.12, 0.46, 0.15);
  m_Tab2_Combo_Parameter2D_1.SetItems(0);
  m_Tab2_Group_2D.AddWidget(&m_Tab2_Combo_Parameter2D_1);

  m_Tab2_Text_Parameters2D_2_Description.SetPos(0.04, 0.08, 0.41, 0.028);
  m_Tab2_Group_2D.AddWidget(&m_Tab2_Text_Parameters2D_2_Description);
  m_Tab2_Combo_Parameter2D_2.SetPos(0.12, 0.41, 0.15);
  m_Tab2_Combo_Parameter2D_2.SetItems(0);
  m_Tab2_Group_2D.AddWidget(&m_Tab2_Combo_Parameter2D_2);

  m_Tab2_Button_Integrate2D.SetPos(0.04, 0.12, 0.34, 0.028);
  m_Tab2_Button_Integrate2D.ConnectOnClicked(BanWindow::OnClickedIntegrate2D);
  m_Tab2_Group_2D.AddWidget(&m_Tab2_Button_Integrate2D);

  m_Tab2_Button_SavePlot.SetPos(0.04, 0.12, 0.25, 0.028);
  m_Tab2_Button_SavePlot.ConnectOnClicked(BanWindow::OnClickedTab2SavePlot);
  m_MainTab.AddWidget(2, &m_Tab2_Button_SavePlot);

  m_Tab2_Button_Evidence.SetPos(0.04, 0.12, 0.20, 0.028);
  m_Tab2_Button_Evidence.ConnectOnClicked(BanWindow::OnClickedTab2Evidence);
  m_MainTab.AddWidget(2, &m_Tab2_Button_Evidence);

  // Add likelihood plot to Tab 2
  m_Tab2_Plot_Results_1D.SetPos(0.38, 0.57, 0.18, 0.68);
//  m_MainTab.AddWidget(2, &m_Tab2_Plot_Results_1D);

  m_Tab2_Plot_Results_2D.SetPos(0.38, 0.57, 0.18, 0.68);
//  m_MainTab.AddWidget(2, &m_Tab2_Plot_Results_2D);
} // InitializeWidgetsTab2

/* This function Initialises all Widgets in Tab 3
 * */
void BanWindow::InitializeWidgetsTab3() {
  // Add Group 2D Plot to Tab 3
  m_Tab3_Group_2D.SetPos(0.03, 0.25, 0.62, 0.25);
  m_Tab3_Group_2D.SetTextWidth(7*0.0105-0.003);
  m_MainTab.AddWidget(3, &m_Tab3_Group_2D);

  // Add elements of Tab 3 group 2D
  m_Tab3_Text_Pulsars2D_Description.SetPos(0.04, 0.08, 0.82, 0.028);
  m_Tab3_Group_2D.AddWidget(&m_Tab3_Text_Pulsars2D_Description);
  m_Tab3_Combo_Pulsars2D.SetPos(0.12, 0.82, 0.15);
  m_Tab3_Combo_Pulsars2D.ConnectOnSelected(BanWindow::OnSelectedTab3Pulsar2D);
  m_Tab3_Group_2D.AddWidget(&m_Tab3_Combo_Pulsars2D);

  m_Tab3_Text_Parameters2D_1_Description.SetPos(0.04, 0.08, 0.77, 0.028);
  m_Tab3_Group_2D.AddWidget(&m_Tab3_Text_Parameters2D_1_Description);
  m_Tab3_Combo_Parameter2D_1.SetPos(0.12, 0.77, 0.15);
  m_Tab3_Group_2D.AddWidget(&m_Tab3_Combo_Parameter2D_1);

  m_Tab3_Text_Parameters2D_2_Description.SetPos(0.04, 0.08, 0.72, 0.028);
  m_Tab3_Group_2D.AddWidget(&m_Tab3_Text_Parameters2D_2_Description);
  m_Tab3_Combo_Parameter2D_2.SetPos(0.12, 0.72, 0.15);
  m_Tab3_Group_2D.AddWidget(&m_Tab3_Combo_Parameter2D_2);

  m_Tab3_Button_Plot2D.SetPos(0.04, 0.12, 0.65, 0.028);
  m_Tab3_Button_Plot2D.ConnectOnClicked(BanWindow::OnClickedTab3Plot2D);
  m_Tab3_Group_2D.AddWidget(&m_Tab3_Button_Plot2D);

  m_Tab3_Button_PlotP1.SetPos(0.17, 0.09, 0.65, 0.028);
  m_Tab3_Button_PlotP1.ConnectOnClicked(BanWindow::OnClickedTab3PlotP1);
  m_Tab3_Group_2D.AddWidget(&m_Tab3_Button_PlotP1);

  m_Tab3_Button_SavePlot.SetPos(0.04, 0.12, 0.56, 0.028);
  m_Tab3_Button_SavePlot.ConnectOnClicked(BanWindow::OnClickedTab3SavePlot);
  m_MainTab.AddWidget(3, &m_Tab3_Button_SavePlot);


  m_Tab3_Plot_1D.SetPos(0.38, 0.57, 0.18, 0.68);
//  m_MainTab.AddWidget(2, &m_Tab3_Plot_1D);
  m_Tab3_Plot_2D.SetPos(0.38, 0.57, 0.18, 0.68);
//  m_MainTab.AddWidget(2, &m_Tab3_Plot_2D);


  // Vertical shift: 38 -> 62-24 = 38
  // Add Group 2D ensemble to Tab 3
  m_Tab3_Group_2D_Ensemble.SetPos(0.03, 0.25, 0.24, 0.25);
  m_Tab3_Group_2D_Ensemble.SetTextWidth(7*0.0105-0.003);
  m_MainTab.AddWidget(3, &m_Tab3_Group_2D_Ensemble);

  // Add elements of Tab 3 group 2D ensemble
  m_Tab3_Text_Pulsars2D_Ensemble_Description.SetPos(0.04, 0.08, 0.44, 0.028);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Text_Pulsars2D_Ensemble_Description);
  m_Tab3_Combo_Pulsars2D_Ensemble.SetPos(0.12, 0.44, 0.15);
  m_Tab3_Combo_Pulsars2D_Ensemble.ConnectOnSelected(BanWindow::OnSelectedTab3Pulsar2DEnsemble);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Combo_Pulsars2D_Ensemble);

  m_Tab3_Text_Parameters2D_1_Ensemble_Description.SetPos(0.04, 0.08, 0.39, 0.028);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Text_Parameters2D_1_Ensemble_Description);
  m_Tab3_Combo_Parameter2D_1_Ensemble.SetPos(0.12, 0.39, 0.15);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Combo_Parameter2D_1_Ensemble);

  m_Tab3_Text_Parameters2D_2_Ensemble_Description.SetPos(0.04, 0.08, 0.34, 0.028);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Text_Parameters2D_2_Ensemble_Description);
  m_Tab3_Combo_Parameter2D_2_Ensemble.SetPos(0.12, 0.34, 0.15);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Combo_Parameter2D_2_Ensemble);

  m_Tab3_Button_Plot2D_MLEnsemble.SetPos(0.04, 0.12, 0.30, 0.028);
  m_Tab3_Button_Plot2D_MLEnsemble.ConnectOnClicked(BanWindow::OnClickedTab3Plot2DMLEnsemble);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Button_Plot2D_MLEnsemble);

  m_Tab3_Button_Plot2D_Ensemble.SetPos(0.04, 0.12, 0.26, 0.028);
  m_Tab3_Button_Plot2D_Ensemble.ConnectOnClicked(BanWindow::OnClickedTab3Plot2DEnsemble);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Button_Plot2D_Ensemble);

  m_Tab3_Button_Plot1D_M.SetPos(0.17, 0.04, 0.26, 0.028);
  m_Tab3_Button_Plot1D_M.ConnectOnClicked(BanWindow::OnClickedTab3Plot1DM);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Button_Plot1D_M);
  m_Tab3_Button_Plot1D_M2.SetPos(0.22, 0.04, 0.26, 0.028);
  m_Tab3_Button_Plot1D_M2.ConnectOnClicked(BanWindow::OnClickedTab3Plot1Dm);
  m_Tab3_Group_2D_Ensemble.AddWidget(&m_Tab3_Button_Plot1D_M2);

  // Add button RunMCMCEnsemble to Tab 3
  m_Tab3_Combo_DataSet_Ensemble.SetPos(0.04, 0.18, 0.15);
  m_MainTab.AddWidget(3, &m_Tab3_Combo_DataSet_Ensemble);

  m_Tab3_Button_RunMCMCEnsemble.SetPos(0.04, 0.14, 0.13, 0.028);
  m_Tab3_Button_RunMCMCEnsemble.ConnectOnClicked(BanWindow::OnClickedTab3RunMCMCEnsemble);
  m_MainTab.AddWidget(3, &m_Tab3_Button_RunMCMCEnsemble);

  m_Tab3_Button_Down.SetPos(0.20, 0.03, 0.18, 0.028);
  m_Tab3_Button_Down.ConnectOnClicked(BanWindow::OnClickedTab3Down);
  m_MainTab.AddWidget(3, &m_Tab3_Button_Down);
  m_Tab3_Button_Up.SetPos(0.24, 0.03, 0.18, 0.028);
  m_Tab3_Button_Up.ConnectOnClicked(BanWindow::OnClickedTab3Up);
  m_MainTab.AddWidget(3, &m_Tab3_Button_Up);
} // InitializeWidgetsTab3

/* This function Initialises all Widgets in Tab 4
 * */
void BanWindow::InitializeWidgetsTab4() {
  // Add pulsar texts to tab 0
  m_Tab4_Text_Pulsars.SetPos(0.03, 0.07, 0.90, 0.028);
  m_MainTab.AddWidget(4, &m_Tab4_Text_Pulsars);

  // Add pulsar combobox to tab 4
  m_Tab4_Combo_Pulsars.SetPos(0.13, 0.90, 0.15);
  m_Tab4_Combo_Pulsars.SetItems(m_oConstants.k);
  m_Tab4_Combo_Pulsars.ConnectOnSelected(BanWindow::OnSelectedTab4Pulsar);
  m_MainTab.AddWidget(4, &m_Tab4_Combo_Pulsars);

  // Add save plot button
  m_Tab4_Button_SavePlot.SetPos(0.04, 0.12, 0.25, 0.028);
  m_Tab4_Button_SavePlot.ConnectOnClicked(BanWindow::OnClickedTab4SavePlot);
  m_MainTab.AddWidget(4, &m_Tab4_Button_SavePlot);

  // Add Sigmaz plot to Tab 4
  m_Tab4_Plot_Sigmaz.SetPos(0.38, 0.57, 0.18, 0.68);
} // InitializeWidgetsTab4

/* This function Initialises all Widgets in Tab 5
 * */
void BanWindow::InitializeWidgetsTab5() {
  // Add Simulation button to tab 5
  m_Tab5_Button_Sensitivity.SetPos(0.03, 0.08, 0.16, 0.028);
  m_Tab5_Button_Sensitivity.ConnectOnClicked(BanWindow::OnClickedTab5Sensitivity);
  m_MainTab.AddWidget(5, &m_Tab5_Button_Sensitivity);
  m_Tab5_Button_Simulate.SetPos(0.03, 0.08, 0.11, 0.028);
  m_Tab5_Button_Simulate.ConnectOnClicked(BanWindow::OnClickedTab5Simulate);
  m_MainTab.AddWidget(5, &m_Tab5_Button_Simulate);

  // Add Sensitivity plot to Tab 5
  m_Tab5_Plot_Sensitivity.SetPos(0.38, 0.57, 0.18, 0.68);
  m_bHave_Sensitivity_Curve = false;

  // Add plot combobox to tab 5
  m_Tab5_Combo_Curve.SetPos(0.03, 0.21, 0.15);
  m_Tab5_Combo_Curve.SetItems(5);
  m_Tab5_Combo_Curve.ConnectOnSelected(BanWindow::OnSelectedTab5Curve);
  m_MainTab.AddWidget(5, &m_Tab5_Combo_Curve);
} // InitializeWidgetsTab5

/* This function sets the content of all Widgets in Tab 0
 * */
void BanWindow::SetWidgetsContentTab0() {
  char strBuf[80];
  int nPulsar;
  double dChisq, dRedChisq;

  m_Tab0_Combo_Pulsars.SetItems(m_oConstants.k);
  for(int a=0; a<m_oConstants.k; a++) {
    m_Tab0_Combo_Pulsars.SetItem(a, m_oConstants.poPulsars[a].strPulsarName);
  } // for a

  m_Tab0_Text_Pulsars.SetText("Pulsar:");
  m_Tab0_Text_Datapoints_Description.SetText("# of Obs:");
  m_Tab0_Text_Datapoints_Value.SetText("");
  m_Tab0_Text_Timespan_Description.SetText("Timespan:");
  m_Tab0_Text_Timespan_Value.SetText("");
  m_Tab0_Text_RMS_Description.SetText("RMS res:");
  m_Tab0_Text_RMS_Value.SetText("");
  m_Tab0_Text_Chisq_Description.SetText("Chisq:");
  m_Tab0_Text_Chisq_Value.SetText("");
  m_Tab0_Text_RedChisq_Description.SetText("R. Chisq:");
  m_Tab0_Text_RedChisq_Value.SetText("");
  m_Tab0_Button_MemMap.SetText("Skymap");
  m_Tab0_Text_MemMap_Description.SetText("of memory");
  m_Tab0_Button_SavePlot.SetText("Save Plot");
  m_Tab0_Button_Fit.SetText("Fit");
  m_Tab0_Text_Fit_Description.SetText("Plot fit");
  m_Tab0_Button_Plk.SetText("Run plk");
  m_Tab0_Text_Plk_Description.SetText("For this pulsar");

  m_Tab0_Plot_Residuals.SetTitle("Residuals");
  m_Tab0_Plot_Residuals.SetXLabel("TOA");
  m_Tab0_Plot_Residuals.SetYLabel("Timing residuals");

  if(m_bHasData) {
    SetResidualsPlot();
    nPulsar = m_Tab0_Combo_Pulsars.GetCurrentItem();

    sprintf(strBuf, "%i", m_oConstants.poPulsars[nPulsar].nObservations);
    m_Tab0_Text_Datapoints_Value.SetText(strBuf);

    sprintf(strBuf, "%.2f yr", (m_oConstants.poPulsars[nPulsar].pdTOA[
  m_oConstants.poPulsars[nPulsar].nObservations-1] - 
  m_oConstants.poPulsars[nPulsar].pdTOA[0]) / SPERYEAR);
    m_Tab0_Text_Timespan_Value.SetText(strBuf);

    sprintf(strBuf, "%.3e", CalcResidualsRms(m_oData, m_oConstants, nPulsar));
    m_Tab0_Text_RMS_Value.SetText(strBuf);

    dChisq = CalcResidualsChisq(m_oData, m_oConstants, nPulsar);
    dRedChisq = dChisq / (m_oConstants.poPulsars[nPulsar].nObservations - m_oConstants.poPulsars[nPulsar].nTempo2Parameters);

    sprintf(strBuf, "%.2f", dChisq);
    m_Tab0_Text_Chisq_Value.SetText(strBuf);

    sprintf(strBuf, "%.2f", dRedChisq);
    m_Tab0_Text_RedChisq_Value.SetText(strBuf);
  }
} // SetWidgetsContentTab0

/* This function sets the content of all Widgets in Tab 1
 * */
void BanWindow::SetWidgetsContentTab1() {
  m_Tab1_Combo_Pulsars.SetItems(m_oConstants.k);
  for(int a=0; a<m_oConstants.k; a++) {
    m_Tab1_Combo_Pulsars.SetItem(a, m_oConstants.poPulsars[a].strPulsarName);
  } // for a

  m_Tab1_Text_Stochastic.SetText("Stochastic sources:");

  m_Tab1_Text_Deterministic.SetText("Deterministic sources:");

  SetModelWidgets();

  m_Tab1_Button_PowellOptimum.SetText("PowellOptimum");
} // SetWidgetsContentTab1

/* This function sets the content of all Widgets in Tab 2
 * */
void BanWindow::SetWidgetsContentTab2() {
  m_Tab2_Button_RunMCMC.SetText("Run MCMC");
  m_Tab2_Group_1D.SetText("Integrate 1D");
  m_Tab2_Text_Pulsars1D_Description.SetText("Pulsar:");
  m_Tab2_Combo_Pulsars1D.SetItems(0);
  m_Tab2_Text_Parameters1D_Description.SetText("Par:");
  m_Tab2_Combo_Parameter1D.SetItems(0);
  m_Tab2_Text_CalcErr_Description.SetText("Bootstrap:");
  m_Tab2_Button_Integrate1D.SetText("Integrate 1D");
  m_Tab2_Group_2D.SetText("Integrate 2D");
  m_Tab2_Text_Pulsars2D_Description.SetText("Pulsar:");
  m_Tab2_Combo_Pulsars2D.SetItems(0);
  m_Tab2_Text_Parameters2D_1_Description.SetText("Par 1:");
  m_Tab2_Combo_Parameter2D_1.SetItems(0);
  m_Tab2_Text_Parameters2D_2_Description.SetText("Par 2:");
  m_Tab2_Combo_Parameter2D_2.SetItems(0);
  m_Tab2_Button_Integrate2D.SetText("Integrate 2D");
  m_Tab2_Button_SavePlot.SetText("Save Plot");
  m_Tab2_Button_Evidence.SetText("Evidence");

  m_Tab2_Plot_Results_1D.SetTitle("MCMC Likelihood analysis");
  m_Tab2_Plot_Results_1D.SetXLabel("Parameter");
  m_Tab2_Plot_Results_1D.SetYLabel("Likelihood / Max(Likelihood)");

  m_Tab2_Plot_Results_2D.SetTitle("MCMC Likelihood analysis");
  m_Tab2_Plot_Results_2D.SetXLabel("Parameter 1");
  m_Tab2_Plot_Results_2D.SetYLabel("Parameter 2");

  // Set the pulsars for the MCMC integration combo boxes
  m_Tab2_Combo_Pulsars1D.SetItems(m_oConstants.k);
  m_Tab2_Combo_Pulsars2D.SetItems(m_oConstants.k);
  for(int a=0; a<m_oConstants.k; a++) {
    m_Tab2_Combo_Pulsars1D.SetItem(a, m_oConstants.poPulsars[a].strPulsarName);
    m_Tab2_Combo_Pulsars2D.SetItem(a, m_oConstants.poPulsars[a].strPulsarName);
  } // for a

  // For the currently selected pulsars, set the parameters in the combo boxes
  Set1DMCMCParameterComboBoxes();
  Set2DMCMCParameterComboBoxes();
} // SetWidgetsContentTab2

/* This function sets the content of all Widgets in Tab 3
 * */
void BanWindow::SetWidgetsContentTab3() {
  m_Tab3_Group_2D.SetText("Plot 2D");
  m_Tab3_Text_Pulsars2D_Description.SetText("Pulsar:");
  m_Tab3_Combo_Pulsars2D.SetItems(0);
  m_Tab3_Text_Parameters2D_1_Description.SetText("Par 1:");
  m_Tab3_Combo_Parameter2D_1.SetItems(0);
  m_Tab3_Text_Parameters2D_2_Description.SetText("Par 2:");
  m_Tab3_Combo_Parameter2D_2.SetItems(0);
  m_Tab3_Button_Plot2D.SetText("Plot 2D");
  m_Tab3_Button_PlotP1.SetText("Plot P1");
  m_Tab3_Button_SavePlot.SetText("Save Plot");

  m_Tab3_Plot_2D.SetTitle("MCMC Likelihood analysis");
  m_Tab3_Plot_2D.SetXLabel("Parameter 1");
  m_Tab3_Plot_2D.SetYLabel("Parameter 2");

  m_Tab3_Group_2D_Ensemble.SetText("Ensemble plot");
  m_Tab3_Text_Pulsars2D_Ensemble_Description.SetText("Pulsar:");
  m_Tab3_Combo_Pulsars2D_Ensemble.SetItems(0);
  m_Tab3_Text_Parameters2D_1_Ensemble_Description.SetText("Par 1:");
  m_Tab3_Combo_Parameter2D_1_Ensemble.SetItems(0);
  m_Tab3_Text_Parameters2D_2_Ensemble_Description.SetText("Par 2:");
  m_Tab3_Combo_Parameter2D_2_Ensemble.SetItems(0);
  m_Tab3_Button_Plot2D_MLEnsemble.SetText("Ensemble ML");
  m_Tab3_Button_Plot2D_Ensemble.SetText("Ensemble 2D");
  m_Tab3_Button_Plot1D_M.SetText("1D M");
  m_Tab3_Button_Plot1D_M2.SetText("1D m");
  m_Tab3_Button_RunMCMCEnsemble.SetText("Ensemble Run");
  m_Tab3_Combo_DataSet_Ensemble.SetItems(0);
  m_Tab3_Button_Down.SetText("<");
  m_Tab3_Button_Up.SetText(">");

  // Set the pulsars for the MCMC integration combo boxes
  m_Tab3_Combo_Pulsars2D.SetItems(m_oConstants.k);
  m_Tab3_Combo_Pulsars2D_Ensemble.SetItems(m_oConstants.k);
  for(int a=0; a<m_oConstants.k; a++) {
    m_Tab3_Combo_Pulsars2D.SetItem(a, m_oConstants.poPulsars[a].strPulsarName);
    m_Tab3_Combo_Pulsars2D_Ensemble.SetItem(a, m_oConstants.poPulsars[a].strPulsarName);
  } // for a

  Set2DPlotParameterComboBoxes();
} // SetWidgetsContentTab3

/* This function sets the content of all Widgets in Tab 4
 * */
void BanWindow::SetWidgetsContentTab4() {
  char strBuf[160];

  m_Tab4_Text_Pulsars.SetText("Pulsar:");

  m_Tab4_Combo_Pulsars.SetItems(m_oConstants.k);
  for(int a=0; a<m_oConstants.k; a++) {
    m_Tab4_Combo_Pulsars.SetItem(a, m_oConstants.poPulsars[a].strPulsarName);
  } // for a

  m_Tab4_Button_SavePlot.SetText("Save Plot");

  sprintf(strBuf, "\\gs\\dz\\u for %s", m_oConstants.poPulsars[
	m_Tab4_Combo_Pulsars.GetCurrentItem()].strPulsarName);

  m_Tab4_Plot_Sigmaz.SetTitle(strBuf);
  m_Tab4_Plot_Sigmaz.SetXLabel("Timescale (yr)");
  m_Tab4_Plot_Sigmaz.SetYLabel("\\gs\\dz\\u");
  m_Tab4_Plot_Sigmaz.SetLogScale(0, true);
  m_Tab4_Plot_Sigmaz.SetLogScale(1, true);
  m_Tab4_Plot_Sigmaz.DrawLine(true);
} // SetWidgetsContentTab4

/* This function sets the content of all Widgets in Tab 5
 * */
void BanWindow::SetWidgetsContentTab5() {
  m_Tab5_Button_Sensitivity.SetText("Sense");
  m_Tab5_Button_Simulate.SetText("Simulate");

  m_Tab5_Combo_Curve.SetItems(5);
  m_Tab5_Combo_Curve.SetItem(0, "DS");
  m_Tab5_Combo_Curve.SetItem(1, "DS Theory");
  m_Tab5_Combo_Curve.SetItem(2, "Bayesian");
  m_Tab5_Combo_Curve.SetItem(3, "DST & B");
  m_Tab5_Combo_Curve.SetItem(4, "DST & DS");
} // SetWidgetsContentTab5

/* This function sets the widgets on the model tab. The amount of widgets etc.
 * can change after adding/removing sources, so this needs to be done in other
 * places than the Initialise as well.
 * */
void BanWindow::SetModelWidgets() {
  int nPulsarStochasticSources, nPulsarDeterministicSources, nPulsar, nIndex, nType,
      nSourceParameter;
  char strBuf[80];

  // Check how many sources belong to the active pulsar
  nPulsarStochasticSources = 0;
  nPulsarDeterministicSources = 0;
  nIndex = 0;
  nPulsar = m_Tab0_Combo_Pulsars.GetCurrentItem();
  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      if(m_oConstants.poSources[s].oSourceType.eID == SID_Deterministic)
	nPulsarDeterministicSources++;
      else
  	nPulsarStochasticSources++;
    } // if SourceWorksOnPulsar
  } // for s

  m_Tab1_StochasticTab.SetTabs(nPulsarStochasticSources);
  nIndex = 0;
  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      if(m_oConstants.poSources[s].oSourceType.eID != SID_Deterministic) {
	// We have another stochastic source for this pulsar. Set the tab text
	m_Tab1_StochasticTab.SetText(nIndex, m_oConstants.poSources[s].oSourceType.strID);

	// Set the parameter combobox
	m_Tab1_Text_StochasticParameters[nIndex].SetPos(0.05, 0.10, 0.72, 0.028);
	m_Tab1_Text_StochasticParameters[nIndex].SetText("Parameter:");
	if(! m_Tab1_StochasticTab.HasWidget(&m_Tab1_Text_StochasticParameters[nIndex]))
	  m_Tab1_StochasticTab.AddWidget(nIndex, &m_Tab1_Text_StochasticParameters[nIndex]);

	m_Tab1_Combo_StochasticParameters[nIndex].SetPos(0.16, 0.72, 0.15);
	m_Tab1_Combo_StochasticParameters[nIndex].SetItems(
	    m_oConstants.poSources[s].oSourceType.nParameters);
	for(int p=0; p<m_oConstants.poSources[s].oSourceType.nParameters; p++) {
	  GetParameterCode(m_oConstants.poSources[s], p, strBuf);
	  m_Tab1_Combo_StochasticParameters[nIndex].SetItem(p, strBuf);
	} // for p
	if(! m_Tab1_StochasticTab.HasWidget(&m_Tab1_Combo_StochasticParameters[nIndex])) {
	  m_Tab1_StochasticTab.AddWidget(nIndex, &m_Tab1_Combo_StochasticParameters[nIndex]);
	} // if HasWidget
	m_Tab1_Combo_StochasticParameters[nIndex].ConnectOnSelected(BanWindow::OnTab1SelectedStochasticParameter);
  
	if(strlen(m_oConstants.poSources[s].strScopeTag) > 0)
	  sprintf(strBuf, "Tag: %s", m_oConstants.poSources[s].strScopeTag);
	else
	  strcpy(strBuf, "");
	m_Tab1_Text_StochasticTag[nIndex].SetPos(0.32, 0.10, 0.72, 0.028);
	m_Tab1_Text_StochasticTag[nIndex].SetText(strBuf);
	if(! m_Tab1_StochasticTab.HasWidget(&m_Tab1_Text_StochasticTag[nIndex])) {
	  m_Tab1_StochasticTab.AddWidget(nIndex, &m_Tab1_Text_StochasticTag[nIndex]);
	} // if HasWidget

	nSourceParameter = m_Tab1_Combo_StochasticParameters[nIndex].GetCurrentItem();

	// Set the checkboxes
	m_Tab1_Text_StochasticFit_Description[nIndex].SetPos(0.05, 0.04, 0.42, 0.028);
	m_Tab1_Text_StochasticFit_Description[nIndex].SetText("Vary:");
	if(! m_Tab1_StochasticTab.HasWidget(&m_Tab1_Text_StochasticFit_Description[nIndex]))
	  m_Tab1_StochasticTab.AddWidget(nIndex, &m_Tab1_Text_StochasticFit_Description[nIndex]);
	m_Tab1_Check_StochasticFit[nIndex].SetPos(0.10, 0.42);
	m_Tab1_Check_StochasticFit[nIndex].SetChecked(bool(m_oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(nSourceParameter)));
	m_Tab1_Check_StochasticFit[nIndex].ConnectOnClicked(BanWindow::OnTab1ClickedStochasticCheckbox);
	if(! m_Tab1_StochasticTab.HasWidget(&m_Tab1_Check_StochasticFit[nIndex]))
	  m_Tab1_StochasticTab.AddWidget(nIndex, &m_Tab1_Check_StochasticFit[nIndex]);

	m_Tab1_Text_StochasticMarginalise_Description[nIndex].SetPos(0.16, 0.10, 0.42, 0.028);
	m_Tab1_Text_StochasticMarginalise_Description[nIndex].SetText("An. Marg.:");
	if(! m_Tab1_StochasticTab.HasWidget(&m_Tab1_Text_StochasticMarginalise_Description[nIndex]))
	  m_Tab1_StochasticTab.AddWidget(nIndex, &m_Tab1_Text_StochasticMarginalise_Description[nIndex]);
	m_Tab1_Check_StochasticMarginalise[nIndex].SetPos(0.27, 0.42);
	m_Tab1_Check_StochasticMarginalise[nIndex].SetChecked(bool(m_oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(nSourceParameter)));
	m_Tab1_Check_StochasticMarginalise[nIndex].ConnectOnClicked(BanWindow::OnTab1ClickedStochasticCheckbox);
	if(false) {		// Stochastic sources are never linear
	  if(! m_Tab1_StochasticTab.HasWidget(&m_Tab1_Check_StochasticMarginalise[nIndex])) {
	    m_Tab1_StochasticTab.AddWidget(nIndex, &m_Tab1_Check_StochasticMarginalise[nIndex]);
	  } // if HasWidget
	} else {
	  m_Tab1_Check_StochasticMarginalise[nIndex].SetChecked(false);
	} // if false

	// Set the delete button
	m_Tab1_Button_StochasticDelete[nIndex].SetPos(0.36, 0.08, 0.42, 0.028);
	m_Tab1_Button_StochasticDelete[nIndex].SetText("Delete");
	m_Tab1_Button_StochasticDelete[nIndex].ConnectOnClicked(
	    BanWindow::OnTab1ClickedStochasticDelete);
	if(! m_Tab1_StochasticTab.HasWidget(&m_Tab1_Button_StochasticDelete[nIndex]))
	  m_Tab1_StochasticTab.AddWidget(nIndex, &m_Tab1_Button_StochasticDelete[nIndex]);

	nIndex++;
      } // if eID
    } // if SourceWorksOnPulsar
  } // for s

  // The deterministic control part of tab1
  m_Tab1_DeterministicTab.SetTabs(nPulsarDeterministicSources);

  nIndex = 0;
  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      if(m_oConstants.poSources[s].oSourceType.eID == SID_Deterministic) {
	// We have another stochastic source for this pulsar. Set the tab text
	nType = 0;
	while(oDetSourceTypes[nType].nID != m_oConstants.poSources[s].oSourceType.nTag) {
	  nType++;
	} // while nType
	m_Tab1_DeterministicTab.SetText(nIndex, oDetSourceTypes[nType].strID);

	// Set the parameter combobox
	m_Tab1_Text_DeterministicParameters[nIndex].SetPos(0.55, 0.10, 0.72, 0.028);
	m_Tab1_Text_DeterministicParameters[nIndex].SetText("Parameter:");
	if(! m_Tab1_DeterministicTab.HasWidget(&m_Tab1_Text_DeterministicParameters[nIndex]))
	  m_Tab1_DeterministicTab.AddWidget(nIndex, &m_Tab1_Text_DeterministicParameters[nIndex]);
	m_Tab1_Combo_DeterministicParameters[nIndex].ConnectOnSelected(BanWindow::OnTab1SelectedDeterministicParameter);

	m_Tab1_Combo_DeterministicParameters[nIndex].SetPos(0.66, 0.72, 0.15);
	m_Tab1_Combo_DeterministicParameters[nIndex].SetItems(
	    m_oConstants.poSources[s].oSourceType.nParameters);
	for(int p=0; p<m_oConstants.poSources[s].oSourceType.nParameters; p++) {
	  GetParameterCode(m_oConstants.poSources[s], p, strBuf);
	  m_Tab1_Combo_DeterministicParameters[nIndex].SetItem(p, strBuf);
	} // for p
	if(! m_Tab1_DeterministicTab.HasWidget(&m_Tab1_Combo_DeterministicParameters[nIndex])) {
	  m_Tab1_DeterministicTab.AddWidget(nIndex, &m_Tab1_Combo_DeterministicParameters[nIndex]);
	} // if HasWidget

	if(strlen(m_oConstants.poSources[s].strScopeTag) > 0)
	  sprintf(strBuf, "Tag: %s", m_oConstants.poSources[s].strScopeTag);
	else
	  strcpy(strBuf, "");
	m_Tab1_Text_DeterministicTag[nIndex].SetPos(0.82, 0.10, 0.72, 0.028);
	m_Tab1_Text_DeterministicTag[nIndex].SetText(strBuf);
	if(! m_Tab1_DeterministicTab.HasWidget(&m_Tab1_Text_DeterministicTag[nIndex])) {
	  m_Tab1_DeterministicTab.AddWidget(nIndex, &m_Tab1_Text_DeterministicTag[nIndex]);
	} // if HasWidget

	nSourceParameter = m_Tab1_Combo_DeterministicParameters[nIndex].GetCurrentItem();

	// Set the checkboxes
	m_Tab1_Text_DeterministicFit_Description[nIndex].SetPos(0.55, 0.04, 0.42, 0.028);
	m_Tab1_Text_DeterministicFit_Description[nIndex].SetText("Vary:");
	if(! m_Tab1_DeterministicTab.HasWidget(&m_Tab1_Text_DeterministicFit_Description[nIndex]))
	  m_Tab1_DeterministicTab.AddWidget(nIndex, &m_Tab1_Text_DeterministicFit_Description[nIndex]);
	m_Tab1_Check_DeterministicFit[nIndex].SetPos(0.60, 0.42);
	m_Tab1_Check_DeterministicFit[nIndex].SetChecked(bool(m_oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(nSourceParameter)));
	m_Tab1_Check_DeterministicFit[nIndex].ConnectOnClicked(BanWindow::OnTab1ClickedDeterministicCheckbox);
	if(! m_Tab1_DeterministicTab.HasWidget(&m_Tab1_Check_DeterministicFit[nIndex]))
	  m_Tab1_DeterministicTab.AddWidget(nIndex, &m_Tab1_Check_DeterministicFit[nIndex]);

	m_Tab1_Text_DeterministicMarginalise_Description[nIndex].SetPos(0.67, 0.10, 0.42, 0.028);
	m_Tab1_Text_DeterministicMarginalise_Description[nIndex].SetText("An. Marg.:");
	if(! m_Tab1_DeterministicTab.HasWidget(&m_Tab1_Text_DeterministicMarginalise_Description[nIndex]))
	  m_Tab1_DeterministicTab.AddWidget(nIndex, &m_Tab1_Text_DeterministicMarginalise_Description[nIndex]);
	m_Tab1_Check_DeterministicMarginalise[nIndex].SetPos(0.77, 0.42);
	m_Tab1_Check_DeterministicMarginalise[nIndex].SetChecked(bool(m_oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(nSourceParameter)));
	m_Tab1_Check_DeterministicMarginalise[nIndex].ConnectOnClicked(BanWindow::OnTab1ClickedDeterministicCheckbox);
	if(oDetSourceTypes[nType].bLinear) {
	  // This is a linear source
	  if(! m_Tab1_DeterministicTab.HasWidget(&m_Tab1_Check_DeterministicMarginalise[nIndex]))
	    m_Tab1_DeterministicTab.AddWidget(nIndex, &m_Tab1_Check_DeterministicMarginalise[nIndex]);
	} else {
	  // This is a nonlinear source
	  if(m_Tab1_DeterministicTab.HasWidget(&m_Tab1_Check_DeterministicMarginalise[nIndex])) {
	    m_Tab1_DeterministicTab.RemoveWidget(&m_Tab1_Check_DeterministicMarginalise[nIndex]);
	  } // if HasWidget
	  m_Tab1_Check_DeterministicMarginalise[nIndex].SetChecked(false);
	} // if bLinear

	// Set the delete button
	m_Tab1_Button_DeterministicDelete[nIndex].SetPos(0.86, 0.08, 0.42, 0.028);
	m_Tab1_Button_DeterministicDelete[nIndex].SetText("Delete");
	m_Tab1_Button_DeterministicDelete[nIndex].ConnectOnClicked(
	    BanWindow::OnTab1ClickedDeterministicDelete);
	if(! m_Tab1_DeterministicTab.HasWidget(&m_Tab1_Button_DeterministicDelete[nIndex]))
	  m_Tab1_DeterministicTab.AddWidget(nIndex, &m_Tab1_Button_DeterministicDelete[nIndex]);

	nIndex++;
      } // if eID
    } // if SourceWorksOnPulsar
  } // for s
} // SetModelWidgets

/* This function (re)sets the content of all the widget according to the member
 * variables m_oConstants and m_oData
 * */
void BanWindow::SetWidgetsContent() {
  // The main tab
  m_MainTab.SetTabs(6);
  m_MainTab.SetText(0, "Pulsars");
  m_MainTab.SetText(1, "Model");
  m_MainTab.SetText(2, "MCMC");
  m_MainTab.SetText(3, "MCMC Adv.");
  m_MainTab.SetText(4, "Statistics");
  m_MainTab.SetText(5, "Simulate");

  // Add button quit to main window
  m_ButtonQuit.SetText("Quit");

  // Add button load to main window
  m_ButtonLoad.SetText("Load .dat");

  // Add button load to main window
  m_ButtonPrint.SetText("Print P");

  m_Text_Copyright.SetText("BANgui v.0.5 (Rutger van Haasteren)");

  SetWidgetsContentTab0();
  SetWidgetsContentTab1();
  SetWidgetsContentTab2();
  SetWidgetsContentTab3();
  SetWidgetsContentTab4();
  SetWidgetsContentTab5();
} // SetWidgetsContent

/* This function sets the correct parameters in the 1D combo boxes for a given
 * pulsar */
void BanWindow::Set1DMCMCParameterComboBoxes() {
  int nPulsar, nSources, nParameters;
  char strMsg[80];

  nPulsar = m_Tab2_Combo_Pulsars1D.GetCurrentItem();
  nSources = 0; nParameters = 0;
  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      // This source is acting on the pulsar
      nSources++;
      for(int p=0; p<m_oConstants.poSources[s].oSourceType.nParameters; p++) {
	if(m_oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	  // This parameter is included in the MCMC
	  nParameters++;
	  GetParameterCode(m_oConstants.poSources[s], p, strMsg);
	  m_Tab2_Combo_Parameter1D.SetItems(nParameters);
	  m_Tab2_Combo_Parameter1D.SetItem(nParameters-1, strMsg);
	} // if FITPARAMETER
      } // for p
    } // if SourceWorksOnPulsar
  } // for  s

  m_Tab2_Combo_Parameter1D.SetCurrentItem(0);
} // Set1DMCMCParameterComboBoxes

/* This function sets the correct parameters in the 1D combo boxes for a given
 * pulsar */
void BanWindow::Set2DMCMCParameterComboBoxes() {
  int nPulsar, nSources, nParameters;
  char strMsg[80];

  nPulsar = m_Tab2_Combo_Pulsars2D.GetCurrentItem();
  nSources = 0; nParameters = 0;
  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      // This source is acting on the pulsar
      nSources++;
      for(int p=0; p<m_oConstants.poSources[s].oSourceType.nParameters; p++) {
	if(m_oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	  // This parameter is included in the MCMC
	  nParameters++;
	  GetParameterCode(m_oConstants.poSources[s], p, strMsg);
	  m_Tab2_Combo_Parameter2D_1.SetItems(nParameters);
	  m_Tab2_Combo_Parameter2D_2.SetItems(nParameters);
	  m_Tab2_Combo_Parameter2D_1.SetItem(nParameters-1, strMsg);
	  m_Tab2_Combo_Parameter2D_2.SetItem(nParameters-1, strMsg);
	} // if FITPARAMETER
      } // for p
    } // if SourceWorksOnPulsar
  } // for  s

  m_Tab2_Combo_Parameter2D_1.SetCurrentItem(0);
  m_Tab2_Combo_Parameter2D_2.SetCurrentItem(0);
  if(nParameters >= 2)
    m_Tab2_Combo_Parameter2D_2.SetCurrentItem(1);
} // Set2DMCMCParameterComboBoxes


/* This function sets the correct parameters in the 2D combo boxes for a given
 * pulsar */
void BanWindow::Set2DPlotParameterComboBoxes() {
  int nPulsar, nSources, nParameters;
  char strMsg[80];

  // First do it for the MCMC plot
  nPulsar = m_Tab3_Combo_Pulsars2D.GetCurrentItem();
  nSources = 0; nParameters = 0;
  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      // This source is acting on the pulsar
      nSources++;
      for(int p=0; p<m_oConstants.poSources[s].oSourceType.nParameters; p++) {
	if(m_oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	  // This parameter is included in the MCMC
	  nParameters++;
	  GetParameterCode(m_oConstants.poSources[s], p, strMsg);
	  m_Tab3_Combo_Parameter2D_1.SetItems(nParameters);
	  m_Tab3_Combo_Parameter2D_2.SetItems(nParameters);
	  m_Tab3_Combo_Parameter2D_1.SetItem(nParameters-1, strMsg);
	  m_Tab3_Combo_Parameter2D_2.SetItem(nParameters-1, strMsg);
	} else if(m_oConstants.poSources[s].oSourceType.eID == SID_Deterministic &&
	    (m_oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p))) {
	  // This is a marginalisation parameter
	  nParameters++;
	  GetParameterCode(m_oConstants.poSources[s], p, strMsg);
	  m_Tab3_Combo_Parameter2D_1.SetItems(nParameters);
	  m_Tab3_Combo_Parameter2D_2.SetItems(nParameters);
	  m_Tab3_Combo_Parameter2D_1.SetItem(nParameters-1, strMsg);
	  m_Tab3_Combo_Parameter2D_2.SetItem(nParameters-1, strMsg);
	} // if FITPARAMETER
      } // for p
    } // if SourceWorksOnPulsar
  } // for  s

  m_Tab3_Combo_Parameter2D_1.SetCurrentItem(0);
  m_Tab3_Combo_Parameter2D_2.SetCurrentItem(0);
  if(nParameters >= 2)
    m_Tab3_Combo_Parameter2D_2.SetCurrentItem(1);

  // Do the same thing, but now for the Ensemble plot
  nPulsar = m_Tab3_Combo_Pulsars2D_Ensemble.GetCurrentItem();
  nSources = 0; nParameters = 0;
  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      // This source is acting on the pulsar
      nSources++;
      for(int p=0; p<m_oConstants.poSources[s].oSourceType.nParameters; p++) {
	if(m_oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	  // This parameter is included in the MCMC
	  nParameters++;
	  GetParameterCode(m_oConstants.poSources[s], p, strMsg);
	  m_Tab3_Combo_Parameter2D_1_Ensemble.SetItems(nParameters);
	  m_Tab3_Combo_Parameter2D_2_Ensemble.SetItems(nParameters);
	  m_Tab3_Combo_Parameter2D_1_Ensemble.SetItem(nParameters-1, strMsg);
	  m_Tab3_Combo_Parameter2D_2_Ensemble.SetItem(nParameters-1, strMsg);
	} else if(m_oConstants.poSources[s].oSourceType.eID == SID_Deterministic &&
	    (m_oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p))) {
	  // This is a marginalisation parameter
	  /*
	  nParameters++;
	  GetParameterCode(m_oConstants.poSources[s], p, strMsg);
	  m_Tab3_Combo_Parameter2D_1_Ensemble.SetItems(nParameters);
	  m_Tab3_Combo_Parameter2D_2_Ensemble.SetItems(nParameters);
	  m_Tab3_Combo_Parameter2D_1_Ensemble.SetItem(nParameters-1, strMsg);
	  m_Tab3_Combo_Parameter2D_2_Ensemble.SetItem(nParameters-1, strMsg);
	  */
	} // if FITPARAMETER
      } // for p
    } // if SourceWorksOnPulsar
  } // for  s

  m_Tab3_Combo_Parameter2D_1_Ensemble.SetCurrentItem(0);
  m_Tab3_Combo_Parameter2D_2_Ensemble.SetCurrentItem(0);
  if(nParameters >= 2)
    m_Tab3_Combo_Parameter2D_2_Ensemble.SetCurrentItem(1);
} // Set2DPlotParameterComboBoxes

/* This function looks up the correct descriptions of parameters for use in the
 * plots.
 *
 * Input: strParameter - A string that describes the parameter in question
 *        strPulsar    - A string with the pulsar name
 *
 * Output: strTitle, strLabel - The strings used in the plot
 * */
void BanWindow::GetParameterDescriptions(const char strParameter[], const char strPulsar[], char *strTitle, char *strLabel) {
  int nDescriptionValue, i;
  int nTempo2Parameter;
  char strBuf[80];
  char strNumber[3];
  char *pstrSin;

  // First check whether it is deterministic or not
  if(strncmp(oParameterDescriptions[0].strParameter, strParameter, 3) == 0) {
    ;
  } else {
    // It's stochastic, see if we can find it
    nDescriptionValue = -1; i = 0;

    // First check whether or not it's a tempo2 parameter
    if(strncmp(strParameter, "tem", 3) == 0) {
      // It´s a tempo2 parameter
      // Figure out which number it is
      strNumber[0] = strParameter[5];
      strNumber[1] = strParameter[6];
      strNumber[2] = '\0';
      if(strNumber[1] == ')')
	strNumber[1] = '\0';
      nTempo2Parameter = atoi(strNumber);
      sprintf(strTitle, "%s: likelihood plot of %s", strPulsar, m_oConstants.poPulsars[0].pstrTempo2Descriptions[nTempo2Parameter]);
      sprintf(strLabel, "%s", m_oConstants.poPulsars[0].pstrTempo2Descriptions[nTempo2Parameter]);
      nDescriptionValue = 0;
    } // if tempo2 parameter

    while(nDescriptionValue == -1 && strcmp(oParameterDescriptions[i].strParameter, "end") != 0) {
      if(strcmp(oParameterDescriptions[i].strParameter, strParameter) == 0) {
	// Found the parameter
	nDescriptionValue = i;
	strcpy(strBuf, strParameter);
	pstrSin = strstr(strBuf, "sin");
	if(! pstrSin) {
	  sprintf(strTitle, "Likelihood plot of %s",
	      oParameterDescriptions[nDescriptionValue].strDescription);
	} else {
	  sprintf(strTitle, "%s: likelihood plot of %s",
	      strPulsar,
	      oParameterDescriptions[nDescriptionValue].strDescription);
	} // if strstr
	sprintf(strLabel, "%s [%s]",
	    oParameterDescriptions[nDescriptionValue].strDescription,
	    oParameterDescriptions[nDescriptionValue].strUnit);
      } // if
      i++;
    } // while nDescriptionValue

    if(nDescriptionValue == -1) {
      // No parameter found
      sprintf(strTitle, "%s: likelihood plot of %s", strPulsar, strParameter);
      sprintf(strLabel, "%s []", strParameter);
    } // if nDescriptionValue
  } // if strncmp
} // GetParameterDescriptions

/* This function returns the parameter-code which is used in the combo-boxes in
 * the gui, given a source and parameter number. This parameter code is also
 * used to look up the descriptions in the pre-defined table.
 *
 * Input:  oSource    - The source structure
 *         nParameter - The parameter number of the source
 *
 * Output: strCode     - The code that describes the parameter
 *
 * */
void BanWindow::GetParameterCode(const SSource &oSource, int nParameter, char *strCode) {
  int nIndex=0;
  switch(oSource.oSourceType.eID) {
    case SID_Deterministic:
      while(strncmp(oDetSourceTypes[nIndex].strID, "end", 3) != 0) {
	if(oSource.oSourceType.nTag == oDetSourceTypes[nIndex].nID) {
	  // Found the detsource
	  break;
	} // if nTag
	nIndex++;
      } // while strncmp
      // Describe the tempo2 codes (Needs to be other than pulsar 0), (also in
      // tempo2_plugin_compat.cpp)
      if(strncmp(oDetSourceTypes[nIndex].strID, "tem", 3) == 0) {
	// Name the code properly
	sprintf(strCode, "%s", m_oConstants.poPulsars[0].pstrTempo2Descriptions[nParameter]);
      } else {
	sprintf(strCode, "%s (%i)", oDetSourceTypes[nIndex].strID, nParameter);
      }
      break;
    default:
      sprintf(strCode, "%s-%s (%i)",
	  oSource.oSourceType.strID,
	  strCorrelationID[(int)oSource.eCorrelation],
	  nParameter);
      break;
  } // switch eID
} // GetParameterCode

/* This function returns the absolute parameter number from the pulsar number
 * and pulsar parameter number (as read from the combo box) so it can be passed
 * to an integration function
 * */
int BanWindow::GetParameterFromPulsarSource(int nPulsar, int nParameter) {
  int nSources=0, nParameters=0, nReturnParameter=0;

  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      // This source is acting on the pulsar
      nSources++;
      for(int p=0; p<m_oConstants.poSources[s].oSourceType.nParameters; p++) {
	if(m_oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) {
	  // This parameter is included in the MCMC
	  nParameters++;
	  if(nParameters == nParameter+1)
	    nReturnParameter = m_oConstants.poSources[s].nFirstParIndex + p;
	} // if FITPARAMETER
      } // for p
    } // if SourceWorksOnPulsar
  } // for  s

  return nReturnParameter;
} // GetParameterFromPulsarSource


/* This function returns the absolute parameter number from the pulsar number
 * and pulsar parameter number (as read from the combo box) so it can be passed
 * to an integration function
 *
 * Marginalisation parameters included
 * */
int BanWindow::GetParameterFromPulsarSourceIncMar(int nPulsar, int nParameter) {
  int nSources=0, nParameters=0, nReturnParameter=0;

  for(int s=0; s<m_oConstants.nSources; s++) {
//    if(m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(m_oConstants, s, nPulsar)) {
      // This source is acting on the pulsar
      nSources++;
      for(int p=0; p<m_oConstants.poSources[s].oSourceType.nParameters; p++) {
	if((m_oConstants.poSources[s].oSourceType.nFitTag & FITPARAMETER(p)) ||
	   (m_oConstants.poSources[s].oSourceType.nMarTag & FITPARAMETER(p))) {
	  // This parameter is included in the MCMC
	  nParameters++;
	  if(nParameters == nParameter+1)
	    nReturnParameter = m_oConstants.poSources[s].nFirstParIndex + p;
	} // if FITPARAMETER
      } // for p
    } // if SourceWorksOnPulsar
  } // for  s

  return nReturnParameter;
} // GetParameterFromPulsarSourceIncMar


/* This function is called when an item is selected in the pulsar combobox
 * */
void BanWindow::OnPulsarTab0Selected(CIWindow *pwThis) {
  ((BanWindow *)pwThis)->m_Tab1_Combo_Pulsars.SetCurrentItem(
      ((BanWindow *)pwThis)->m_Tab0_Combo_Pulsars.GetCurrentItem());

  if(! ((BanWindow *)pwThis)->m_bHasData) return;
  ((BanWindow *)pwThis)->SetWidgetsContent();
} // OnPulsarTab0Selected

/* This function is called when an item is selected in the pulsar combobox
 * */
void BanWindow::OnPulsarTab1Selected(CIWindow *pwThis) {
  ((BanWindow *)pwThis)->m_Tab0_Combo_Pulsars.SetCurrentItem(
      ((BanWindow *)pwThis)->m_Tab1_Combo_Pulsars.GetCurrentItem());

  if(! ((BanWindow *)pwThis)->m_bHasData) return;
  ((BanWindow *)pwThis)->SetWidgetsContent();
} // OnPulsarTab1Selected

/* This function will react on the button click to load the residuals from a
 * file
 * */
void BanWindow::OnClickedLoad(CIWindow *pwThis) {
  ((BanWindow *)pwThis)->LoadPulsarData();

  OnPulsarTab0Selected(pwThis);
} // OnClickedLoad

void BanWindow::OnClickedQuit(CIWindow *pwThis) {
  pwThis->SetDone(true);
} // OnClickedQuit

/* This function will react on the button click to load the residuals from a
 * file
 * */
void BanWindow::OnClickedPrint(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  PrintParameters(pwBanThis->m_oConstants);
} // OnClickedPrint

/* This function is called when the SavePlot button is clicked on the residuals tab
 * */
void BanWindow::OnClickedTab0SavePlot(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  char strDeviceName[160];
  char strTimeString[80];
  time_t tNow;
  struct tm *ptiTimeInfo;

  time(&tNow);
  ptiTimeInfo = localtime(&tNow);
  strftime(strTimeString, 80, "%Y%m%d%H%M%S", ptiTimeInfo);
  printf("Plotting with timestamp %s\n", strTimeString);

  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab0_Plot_Residuals)) {
    sprintf(strDeviceName, "residuals-%s-%s.ps/ps", pwBanThis->m_oConstants.poPulsars[
	pwBanThis->m_Tab0_Combo_Pulsars.GetCurrentItem()].strPulsarName, strTimeString);
    pwBanThis->m_Tab0_Plot_Residuals.SavePlot(strDeviceName);
  } // if HasWidget

  if (1) {
    // Also write a residuals.txt file:
    CVector vdX, vdY, vdYErr;
    vdX.Initialize(pwBanThis->m_oConstants.poPulsars[0].nObservations);
    vdY.Initialize(pwBanThis->m_oConstants.poPulsars[0].nObservations);
    vdYErr.Initialize(pwBanThis->m_oConstants.poPulsars[0].nObservations);
    for(int i=0; i<pwBanThis->m_oConstants.poPulsars[0].nObservations; i++) {
      vdX[i] = pwBanThis->m_oConstants.poPulsars[0].pdTOA[i];
      vdY[i] = pwBanThis->m_oConstants.poPulsars[0].pdResiduals[i];
      vdYErr[i] = pwBanThis->m_oConstants.poPulsars[0].pdDeltaResiduals[i];
    } // for i
    WritePlot("plotresiduals.txt", vdX, vdY, vdYErr);
  } // if 0
} // OnClickedTab0SavePlot


/* This function is called when the MemMap-checkbox is clicked on the pulsar
 * tab. This function will make a sky-map of the pta sensitivity, and write it
 * to a datafile
 * */
void BanWindow::OnClickedTab0MemMap(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  // MakeArraySensitivityScalingPlot(pwBanThis->m_oConstants, "sensitivity.txt");
} // OnClickedTab0MemMap

/* This function is called when the fit-checkbox is clicked on the pulsar tab.
 * Either add or remove the line that shows the fit to the tempo2 fitting parameters.
 * */
void BanWindow::OnFitTab0Selected(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  SConstantsType oNewConstants;
  SDataType oNewData;
  CVector vdParameters, vdFitData, vdData;
  CMatrix mdPseudoInverse, mdTemp1, mdTemp2, mdTemp3, mdCInv;
  int nPulsar;
  bool *pbPulsars;

  nPulsar = pwBanThis->m_Tab1_Combo_Pulsars.GetCurrentItem();

  if(! pwBanThis->m_bHasData) return;
  if(pwBanThis->m_oConstants.poPulsars[nPulsar].nTempo2Parameters == 0) {
    printf("WARNING: Pulsar %i has nTempo2Parameters == 0, skipping fit\n", nPulsar);
    return;
  } // if nTempo2Parameters

  pbPulsars = new bool[pwBanThis->m_oConstants.k];
  for(int p=0; p<pwBanThis->m_oConstants.k; p++)
    pbPulsars[p] = false;
  pbPulsars[nPulsar] = true;

#if 1
      CopySpecificPulsarsToObjects(pwBanThis->m_oConstants, pwBanThis->m_oData,
	  oNewConstants, oNewData, pbPulsars);
      vdData = oNewData.vdData;
      pwBanThis->m_oData.mdC.Initialize(pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations, pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations);
      mdCInv.Initialize(pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations, pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations);
//      SetCoherenceMatrix(pwBanThis->m_oData, pwBanThis->m_oConstants);
//      mdCInv = oNewData.mdC.InverseChol();
      mdTemp1 = oNewData.mdBlock;

      for(int i=0;
	  i<pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations; i++) {
	for(int j=0;
	    j<pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations; j++) {
	  mdCInv[i][j] = 0;
	  if(i == j)
	    mdCInv[i][j] = 1.0 / (pwBanThis->m_oConstants.poPulsars[nPulsar].pdDeltaResiduals[i] *
		pwBanThis->m_oConstants.poPulsars[nPulsar].pdDeltaResiduals[i]);
	} // for j
      } // for i
//      PrintParameters(oNewConstants);
#else
      vdData.Initialize(pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations);
      mdTemp1.Initialize(pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations,
	  pwBanThis->m_oConstants.poPulsars[nPulsar].nTempo2Parameters);
      mdCInv.Initialize(pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations,
	  pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations);
      for(int i=0;
	  i<pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations; i++) {
	vdData[i] = pwBanThis->m_oConstants.poPulsars[nPulsar].pdResiduals[i];
	for(int j=0;
	    j<pwBanThis->m_oConstants.poPulsars[nPulsar].nTempo2Parameters; j++) {
	  mdTemp1[i][j] = double(
	      pwBanThis->m_oConstants.poPulsars[nPulsar].mdTempo2ParameterDerivative[i][j]);
	} // for j
	for(int j=0;
	    j<pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations; j++) {
	  mdCInv[i][j] = 0;
	  if(i == j)
	    mdCInv[i][j] = 1.0 / (pwBanThis->m_oConstants.poPulsars[nPulsar].pdDeltaResiduals[i] *
		pwBanThis->m_oConstants.poPulsars[nPulsar].pdDeltaResiduals[i]);
	} // for j
      } // for i
#endif

      // Now that we have the tempo2 parameter derivatives, find the parameters
      mdTemp2 = mdTemp1[LO_TRANSPOSE];
      mdPseudoInverse = mdTemp2 * mdCInv * mdTemp1;
      mdPseudoInverse.InvertChol();
      mdTemp3 = mdPseudoInverse *  mdTemp2;
      vdParameters = mdTemp3 * mdCInv * vdData;

      // Now calculate the parameters
//      vdFitData = mdTemp1 * vdParameters;
      vdFitData = mdTemp1 * (mdPseudoInverse * (mdTemp2 * (mdCInv * vdData)));

      // Add the fit to the plot
      pwBanThis->m_Tab0_Plot_Residuals.SetYFit(
	  pwBanThis->m_oConstants.poPulsars[nPulsar].nObservations, vdFitData.m_pdData);
//      break;
//    case false:
//      pwBanThis->m_Tab0_Plot_Residuals.UnsetYFit();
//      break;
//  } // switch GetChecked
  delete[] pbPulsars;
} // OnFitTab0Selected

/* This function is called when the fit-checkbox is clicked on the pulsar tab.
 * Either add or remove the line that shows the fit to the tempo2 fitting parameters.
 * */
void BanWindow::OnSelectedTab0Plk(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  pwBanThis->m_nPlkPulsar = pwBanThis->m_Tab1_Combo_Pulsars.GetCurrentItem();
} // OnSelectedTab0Plk


/* This function is called when the stochastic parameter is selected in the
 * model tab
 * */
void BanWindow::OnTab1SelectedStochasticParameter(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  pwBanThis->SetModelWidgets();
} // OnTab1SelectedStochasticParameter

/* This function is called when the checkbox of a stochastic source (de)selected
 * */
void BanWindow::OnTab1ClickedStochasticCheckbox(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  int nIndex, nSourceParameter, nSource, nSourceIndex, nPulsar;

  nIndex = 0;
  nPulsar = pwBanThis->m_Tab0_Combo_Pulsars.GetCurrentItem();
  for(int s=0; s<pwBanThis->m_oConstants.nSources; s++) {
//    if(pwBanThis->m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(pwBanThis->m_oConstants, s, nPulsar)) {
      if(pwBanThis->m_oConstants.poSources[s].oSourceType.eID != SID_Deterministic) {
	if(nIndex == pwBanThis->m_Tab1_StochasticTab.GetActiveTab()) {
	  nSource = s;
	  nSourceIndex = nIndex;
	}
	nIndex++;
      } // if SID_Deterministic
    } // if SourceWorksOnPulsar
  } // for s

  nSourceParameter = pwBanThis->m_Tab1_Combo_StochasticParameters[nSourceIndex].GetCurrentItem();

  switch(pwBanThis->m_Tab1_Check_StochasticFit[nSourceIndex].GetChecked()) {
    case true:
      if(! (pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag & FITPARAMETER(nSourceParameter))) {
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag | FITPARAMETER(nSourceParameter);
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag & ~FITPARAMETER(nSourceParameter);
	  pwBanThis->m_Tab1_Check_StochasticMarginalise[nSourceIndex].SetChecked(false);
      } // if nFitTag
      break;
    case false:
      if((pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag & FITPARAMETER(nSourceParameter))) {
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag & ~FITPARAMETER(nSourceParameter);
//	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag = 
//	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag | FITPARAMETER(nSourceParameter);
//	  pwBanThis->m_Tab1_Check_StochasticMarginalise[nSourceIndex].SetChecked(true);
      } // if nFitTag
    default:
      break;
  } // switch

  switch(pwBanThis->m_Tab1_Check_StochasticMarginalise[nSourceIndex].GetChecked()) {
    case true:
      if(! pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag & FITPARAMETER(nSourceParameter)) {
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag | FITPARAMETER(nSourceParameter);
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag & ~FITPARAMETER(nSourceParameter);
	  pwBanThis->m_Tab1_Check_StochasticFit[nSourceIndex].SetChecked(false);
      } // if nMarTag
      break;
    case false:
      if(pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag & FITPARAMETER(nSourceParameter)) {
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag & ~FITPARAMETER(nSourceParameter);
//	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag = 
//	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag | FITPARAMETER(nSourceParameter);
//	  pwBanThis->m_Tab1_Check_StochasticFit[nSourceIndex].SetChecked(true);
      } // if nMarTag
    default:
      break;
  } // switch

  pwBanThis->SetWidgetsContentTab2();
  pwBanThis->SetWidgetsContentTab3();
  SetNumberOfParametersFromSources(pwBanThis->m_oConstants, pwBanThis->m_oData);
} // OnTab1ClickedStochasticCheckbox

/* This function is called when a stochastic source is about to be deleted
 * */
void BanWindow::OnTab1ClickedStochasticDelete(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
} // OnTab1ClickedStochasticDelete

/* This function is called when the deterministic parameter is selected in the
 * model tab
 * */
void BanWindow::OnTab1SelectedDeterministicParameter(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  pwBanThis->SetModelWidgets();
} // OnTab1SelectedDeterministicParameter

/* This function is called when the checkbox of a deterministic source (de)selected
 * */
void BanWindow::OnTab1ClickedDeterministicCheckbox(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  int nIndex, nSourceParameter, nSource, nSourceIndex, nPulsar;

  nIndex = 0;
  nPulsar = pwBanThis->m_Tab0_Combo_Pulsars.GetCurrentItem();
  for(int s=0; s<pwBanThis->m_oConstants.nSources; s++) {
//    if(pwBanThis->m_oConstants.poSources[s].pbScope[nPulsar]) {
    if(SourceWorksOnPulsar(pwBanThis->m_oConstants, s, nPulsar)) {
      if(pwBanThis->m_oConstants.poSources[s].oSourceType.eID == SID_Deterministic) {
	if(nIndex == pwBanThis->m_Tab1_DeterministicTab.GetActiveTab()) {
	  nSource = s;
	  nSourceIndex = nIndex;
	}
	nIndex++;
      } // if SID_Deterministic
    } // if SourceWorksOnPulsar
  } // for s

  nSourceParameter = pwBanThis->m_Tab1_Combo_DeterministicParameters[nSourceIndex].GetCurrentItem();

  switch(pwBanThis->m_Tab1_Check_DeterministicFit[nSourceIndex].GetChecked()) {
    case true:
      if(! (pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag & FITPARAMETER(nSourceParameter))) {
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag | FITPARAMETER(nSourceParameter);
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag & ~FITPARAMETER(nSourceParameter);
	pwBanThis->m_Tab1_Check_DeterministicMarginalise[nSourceIndex].SetChecked(false);
      } // if nFitTag
      break;
    case false:
      if(pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag & FITPARAMETER(nSourceParameter)) {
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag & ~FITPARAMETER(nSourceParameter);
//	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag = 
//	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag | FITPARAMETER(nSourceParameter);
//	  pwBanThis->m_Tab1_Check_DeterministicMarginalise[nSourceIndex].SetChecked(true);
      } // if nFitTag
    default:
      break;
  } // switch

  switch(pwBanThis->m_Tab1_Check_DeterministicMarginalise[nSourceIndex].GetChecked()) {
    case true:
      if(! (pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag & FITPARAMETER(nSourceParameter))) {
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag | FITPARAMETER(nSourceParameter);
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag & ~FITPARAMETER(nSourceParameter);
	pwBanThis->m_Tab1_Check_DeterministicFit[nSourceIndex].SetChecked(false);
      } // if nMarTag
      break;
    case false:
      if(pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag & FITPARAMETER(nSourceParameter)) {
	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag = 
	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nMarTag & ~FITPARAMETER(nSourceParameter);
//	pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag = 
//	  pwBanThis->m_oConstants.poSources[nSource].oSourceType.nFitTag | FITPARAMETER(nSourceParameter);
//	  pwBanThis->m_Tab1_Check_DeterministicFit[nSourceIndex].SetChecked(true);
      } // if nMarTag
    default:
      break;
  } // switch

  pwBanThis->SetWidgetsContentTab2();
  pwBanThis->SetWidgetsContentTab3();
  SetNumberOfParametersFromSources(pwBanThis->m_oConstants, pwBanThis->m_oData);
} // OnTab1ClickedDeterministicCheckbox

/* This function is called when a deterministic source is about to be deleted
 * */
void BanWindow::OnTab1ClickedDeterministicDelete(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
} // OnTab1ClickedDeterministicDelete

/* This function is called when the user has requested to calculate the ML with
 * the Powell iteration technique
 * */
void BanWindow::OnClickedTab1PowellOptimum(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  CVector vdMax;

//  PowellOptimum(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, vdMax);
} // OnClickedTab1PowellOptimum

/* This function can be called from the plug-in. It should simply run the MCMC
 * chain */
void BanWindow::RunPluginMCMC() {
  BanWindow::OnClickedRunMCMC(this);
} // RunPluginMCMC

/* This function can be called from the plug-in. It should simply run the MCMC
 * chain */
void BanWindow::RunPluginPlot1D(int nParameter1) {
  BanWindow::OnClickedRunPlot1D(this, nParameter1);
} // RunPluginPlot1D

/* This function can be called from the plug-in. It should simply run the MCMC
 * chain */
void BanWindow::RunPluginPlot2D(int nParameter1, int nParameter2) {
  BanWindow::OnClickedRunPlot2D(this, nParameter1, nParameter2);
} // RunPluginPlot2D

/* This function can be called from the plug-in.
 * chain */
void BanWindow::RunPluginEnsemble() {
  BanWindow::OnClickedTab3RunMCMCEnsemble(this);
} // RunPluginEnsemble

/* This function can be called from the plug-in.
 * chain */
void BanWindow::RunPlugin1DEnsembleMCMC() {
  BanWindow::OnClickedTab3Plot1Dm(this);
} // RunPlugin1DEnsembleMCMC

/* This function can be called from the plug-in.
 * chain */
void BanWindow::RunPlugin1DEnsembleMar() {
  BanWindow::OnClickedTab3Plot1DM(this);
} // RunPlugin1DEnsembleMar

/* This function reduces the data
 */
void BanWindow::RunReduce() {
//  ReduceData(m_oData, m_oConstants);
} // RunReduce

/* This function runs the affine invariant sampler
 */
void BanWindow::RunAffine() {
  SetBlockMatrix(m_oData, m_oConstants);
  SetGeometricPart(m_oData, m_oConstants);
  m_oData.mdC.Initialize(m_oConstants.n, m_oConstants.n);

  char strFileName[240];
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(rng, (unsigned int)time(NULL));
  int nUnique = gsl_rng_uniform_int(rng, 10000);

  sprintf(strFileName, "./mcmcdata.%s-%i-%i.dat", m_strProcessorName, m_nProcessRank, nUnique);
  while(FileExists(strFileName)) {
    nUnique = gsl_rng_uniform_int(rng, 10000);
    sprintf(strFileName, "./mcmcdata.%s-%i-%i.dat", m_strProcessorName, m_nProcessRank, nUnique);
  } // while FileExists

  printf("Starting affine invariant sampler...\n");
  EnsembleSampler(m_oData,
      m_oParameters,
      m_oConstants, strFileName);
} // RunReduce

/* This function produces a file that can be use to interpolate a
 * powerlaw-signal in a data reduction scheme. Will be used instead of a
 * data reduction process
 */
void BanWindow::RunPlInterpolation() {
//  CreatePlInterpolation(m_oData, m_oConstants);
} // RunPlInterpolation

/* This function sets the number of pulsars (allocates memory)
 */
void BanWindow::AllocatePulsars(int nPsr, int *pnPsrObs) {
  if(m_oConstants.poPulsars) {
    for(int a=0; a<m_oConstants.k; ++a) {
      delete[] m_oConstants.poPulsars[a].pdTOA;
      delete[] m_oConstants.poPulsars[a].pdResiduals;
      delete[] m_oConstants.poPulsars[a].pdDeltaResiduals;
      delete[] m_oConstants.poPulsars[a].pdFreq;
      delete[] m_oConstants.poPulsars[a].pbFlagSet;

      for(int i=0; i<m_oConstants.poPulsars[a].nObservations; ++i) {
	delete[] m_oConstants.poPulsars[a].pstrFlags[i];
      } // for i
      delete[] m_oConstants.poPulsars[a].pstrFlags;
    } // for a

    delete[] m_oConstants.poPulsars;
  } // if poPulsars

  m_oConstants.k = nPsr;
  m_oConstants.poPulsars = new SPulsar[nPsr];
  for(int a=0; a<nPsr; ++a) {
    m_oConstants.poPulsars[a].nObservations = pnPsrObs[a];
    m_oConstants.poPulsars[a].pdTOA = new double[pnPsrObs[a]];
    m_oConstants.poPulsars[a].pdResiduals = new double[pnPsrObs[a]];
    m_oConstants.poPulsars[a].pdDeltaResiduals = new double[pnPsrObs[a]];
    m_oConstants.poPulsars[a].pdFreq = new double[pnPsrObs[a]];
    m_oConstants.poPulsars[a].pbFlagSet = new bool[pnPsrObs[a]];
    m_oConstants.poPulsars[a].pstrFlags = new char*[pnPsrObs[a]];
    for(int i=0; i<pnPsrObs[a]; ++i) {
      m_oConstants.poPulsars[a].pstrFlags[i] = new char[128];
    } // for i
  } // for a
} // SetPulsars

/* This function allocates memory for an extra source, and copies all the data
 * to the newly allocated memory chunk
 */
void BanWindow::AddSource() {
  SSource *poNewSources = new SSource[++m_oConstants.nSources];

  for(int s=0; s<m_oConstants.nSources-1; ++s) {
    // Copy the old array here
    poNewSources[s] = m_oConstants.poSources[s]; // This copies all but pbScope
/*
    // The reserved memory for pbScope is also transferred. Undeep, but that's
    // ok. So comment this out
    poNewSources[s].pbScope = new bool[m_oConstants.k];
    for(int a=0; a<m_oConstants.k; ++a) {
      poNewSources[s].pbScope[a] = m_oConstants.poSources[s].pbScope[a];
    } // for a
    delete[] m_oConstants.poSources[s].pbScope;
    m_oConstants.poSources[s].pbScope = NULL;
    */
  } // for s
  delete[] m_oConstants.poSources;
  m_oConstants.poSources = poNewSources;

  m_oConstants.poSources[m_oConstants.nSources-1].pbScope = new bool[m_oConstants.k];
  m_oConstants.poSources[m_oConstants.nSources-1].oAmpAcc.bSet = false;
  m_oConstants.poSources[m_oConstants.nSources-1].oIntAcc.bSet = false;
  m_oConstants.poSources[m_oConstants.nSources-1].oIntAcc.ppIntSpline = NULL;
  m_oConstants.poSources[m_oConstants.nSources-1].oIntAcc.ppIntAccel = NULL;
} // AddSource

/* This function is called when the Plot1D button is clicked
 * */
void BanWindow::OnClickedRunPlot1D(CIWindow *pwThis, int n1) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  int nParameter1;

  double pdLevels[3];
  CVector vdX, vdY;

  nParameter1 = n1;

  Calculate1DPlot(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, nParameter1, vdX, vdY);


  WritePlot("plot-1d.txt", vdX, vdY);
} // OnClickedRunPlot1D


/* This function is called when the Plot2D button is clicked
 * */
void BanWindow::OnClickedRunPlot2D(CIWindow *pwThis, int n1, int n2) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  int nParameter1, nParameter2;
  double dMax, dMin;
  double pdLevels[3];
  double dSigma1, dSigma2, dSigma3, dVolume, dCum;
  char strPulsarName[80], strXLabel[80], strYLabel[80],
       strTitle[80], strParameter1[80], strParameter2[80];
  int nStep, nP1, nP2, nPlotPoints;
  CVector vdX, vdY;
  CMatrix mdZ, mdZML, mdZT2, mdTemp;

  /*
  nParameter1 = pwBanThis->GetParameterFromPulsarSourceIncMar(
      pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem(),
      pwBanThis->m_Tab3_Combo_Parameter2D_1.GetCurrentItem());
  nParameter2 = pwBanThis->GetParameterFromPulsarSourceIncMar(
      pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem(),
      pwBanThis->m_Tab3_Combo_Parameter2D_2.GetCurrentItem());
      */
  nParameter1 = 474;
  nParameter2 = 475;
  nParameter1 = n1;
  nParameter2 = n2;


  Calculate3DPlot(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, nParameter1, nParameter2, vdX, vdY, mdZ);


#if 0
  // Calculate the sigma's for mdZ
  mdTemp = mdZ;
  dMin = double(mdZ[0][0]);
  dMax = dMin;
  for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
      mdTemp[i][j] = double(mdZ[i][j]);
      if(dMin > double(mdTemp[i][j]))
	dMin = double(mdTemp[i][j]);
      if(dMax < double(mdTemp[i][j]))
	dMax = double(mdTemp[i][j]);
    } // for j
  } // for i

  dVolume = 0;
  for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
      mdZ[i][j] = double(mdTemp[i][j]) / dMax;
      dVolume += double(mdZ[i][j]);
    } // for j
  } // for i

  mdTemp = mdZ;
  QuickSort(mdTemp.m_pdData, 0, mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1);


  dCum = 0;
  nStep = 0;
  for(int i=mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1; i >= 0; i--) {
    dCum += mdTemp.m_pdData[i];
    switch(nStep) {
      case 0:
	if(dCum > 0.68268949 * dVolume) {
	  nStep++;
	  dSigma1 = mdTemp.m_pdData[i];
	} // if
	break;
      case 1:
	if(dCum > 0.95449974 * dVolume) {
	  nStep++;
	  dSigma2 = mdTemp.m_pdData[i];
	} // if
	break;
      case 2:
	if(dCum > 0.99730024 * dVolume) {
	  nStep++;
	  dSigma3 = mdTemp.m_pdData[i];
	} // if
	break;
      default:
	break;
    } // switch
  } // for i
  printf("dSigma1: %f,  dSigma2: %f,  dSigma3: %f\n", dSigma1, dSigma2, dSigma3);

  pdLevels[0] = dSigma1;
  pdLevels[1] = dSigma2;
  pdLevels[2] = dSigma3;

  // Check whether we have some parameter descriptions available
  pwBanThis->m_Tab3_Combo_Parameter2D_1.CopyItemText(
      pwBanThis->m_Tab3_Combo_Parameter2D_1.GetCurrentItem(),
      strParameter1);
  pwBanThis->m_Tab3_Combo_Parameter2D_2.CopyItemText(
      pwBanThis->m_Tab3_Combo_Parameter2D_2.GetCurrentItem(),
      strParameter2);
  pwBanThis->m_Tab3_Combo_Pulsars2D.CopyItemText(
      pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem(), strPulsarName);
  pwBanThis->GetParameterDescriptions(
      strParameter1, strPulsarName, strTitle, strXLabel);
  pwBanThis->GetParameterDescriptions(
      strParameter2, strPulsarName, strTitle, strYLabel);

  pwBanThis->m_Tab3_Plot_2D.SetTitle(strTitle);
  pwBanThis->m_Tab3_Plot_2D.SetXLabel(strXLabel);
  pwBanThis->m_Tab3_Plot_2D.SetYLabel(strYLabel);

  pwBanThis->m_Tab3_Plot_2D.SetPlot(mdZ.m_pnDimSize[0], mdZ.m_pnDimSize[1], vdX.m_pdData, vdY.m_pdData, mdZ.m_pdData);
  pwBanThis->m_Tab3_Plot_2D.SetLevels(3, pdLevels);


  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_1D)) {
    pwBanThis->m_MainTab.RemoveWidget(&pwBanThis->m_Tab3_Plot_1D);
  } // if HasWidget

  if(! pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_2D)) {
    pwBanThis->m_MainTab.AddWidget(3, &pwBanThis->m_Tab3_Plot_2D);
  } // if HasWidget
#endif


  CVector vdXX, vdYY, vdZZ;
  int nCount=0;
  vdXX.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);
  vdYY.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);
  vdZZ.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);

  for(int i=0; i<mdZ.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdZ.m_pnDimSize[1]; j++) {
      vdXX[nCount] = double(vdX[i]);
      vdYY[nCount] = double(vdY[j]);
      vdZZ[nCount] = double(mdZ[i][j]);
      nCount++;
    } // for j
  } // for i

  WritePlot("plot-2d.txt", vdXX, vdYY, vdZZ);

} // Plot2D

/* This function is called when the RunMCMC button is clicked
 * */
void BanWindow::OnClickedRunMCMC(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  SetBlockMatrix(pwBanThis->m_oData, pwBanThis->m_oConstants);
  SetGeometricPart(pwBanThis->m_oData, pwBanThis->m_oConstants);
  pwBanThis->m_oData.mdC.Initialize(pwBanThis->m_oConstants.n, pwBanThis->m_oConstants.n);
//  SetCoherenceMatrix(pwBanThis->m_oData, pwBanThis->m_oConstants);
//  printf("0,0: %e\n1,1: %e\n", double(pwBanThis->m_oData.mdC[0][0]), double(pwBanThis->m_oData.mdC[1][1]));
//  PrintMatrix(pwBanThis->m_oData.mdC);
#ifdef HAVE_MPI
	if(pwBanThis->m_nProcesses >= 2)
	{
	  printf("Entering MPI-mcmc mode...\n");
	  printf("nProcessRank: %i of %i\n", pwBanThis->m_nProcessRank, pwBanThis->m_nProcesses);
	  MCMCGaussServer(&pwBanThis->m_MPI_Parameters,
	      pwBanThis->m_oData,
	      pwBanThis->m_oParameters,
	      pwBanThis->m_oConstants,
	      NULL);
	}
	else
#endif // HAVE_MPI
	  {
	    char strFileName[240];
	    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
	    gsl_rng_set(rng, (unsigned int)time(NULL));
	    int nUnique = gsl_rng_uniform_int(rng, 10000);

	    sprintf(strFileName, "./mcmcdata.%s-%i-%i.dat", pwBanThis->m_strProcessorName, pwBanThis->m_nProcessRank, nUnique);
	    while(FileExists(strFileName)) {
	      nUnique = gsl_rng_uniform_int(rng, 10000);
	      sprintf(strFileName, "./mcmcdata.%s-%i-%i.dat", pwBanThis->m_strProcessorName, pwBanThis->m_nProcessRank, nUnique);
	    } // while FileExists

	    if(0) {
	      printf("Starting affine invariant sampler...\n");
	      EnsembleSampler(pwBanThis->m_oData,
		  pwBanThis->m_oParameters,
		  pwBanThis->m_oConstants, strFileName);
	    } else {
	      printf("Entering normal mcmc mode...\n");
	      MCMCGauss(pwBanThis->m_oData,
		  pwBanThis->m_oParameters,
		  pwBanThis->m_oConstants, strFileName);
            } // if sampler
	  }
} // OnClickedLoad

/* This function is called when the user selects a new pulsar from the MCMC
 * Integration 1D Group in Tab 2
 * */
void BanWindow::OnSelectedTab2Pulsar1D(CIWindow *pwThis) {
  ((BanWindow*)pwThis)->Set1DMCMCParameterComboBoxes();
} // OnSelectedTab2Pulsar1D

/* This function is called when the user selects a new pulsar from the MCMC
 * Integration 1D Group in Tab 2
 * */
void BanWindow::OnSelectedTab2Pulsar2D(CIWindow *pwThis) {
  ((BanWindow*)pwThis)->Set2DMCMCParameterComboBoxes();
} // OnSelectedTab2Pulsar2D


/* This function is called when the user selects a new pulsar from the MCMC
 * Integration 1D Group in Tab 3
 * */
void BanWindow::OnSelectedTab3Pulsar2D(CIWindow *pwThis) {
  ((BanWindow*)pwThis)->Set2DPlotParameterComboBoxes();
} // OnSelectedTab3Pulsar2D

/* This function is called when the Integrate button is clicked on the MCMC tab
 * */
void BanWindow::OnClickedIntegrate1D(CIWindow *pwThis) {
  CVector vdX, vdY, vdYErr;
  int nParameter1;
  char strPulsarName[80], strXLabel[80], strYLabel[80], strTitle[80], strParameter[80];
  bool bCalcErr;
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  // Check which parameter to integrate and whether to bootstrap
  nParameter1 = pwBanThis->GetParameterFromPulsarSource(
      pwBanThis->m_Tab2_Combo_Pulsars1D.GetCurrentItem(),
      pwBanThis->m_Tab2_Combo_Parameter1D.GetCurrentItem());
  bCalcErr = pwBanThis->m_Tab2_Check_CalcErr.GetChecked();

  // We don't use the bootstrap in the "OnTheFly" method
  bCalcErr = false;

  // Calculate the integrated likelihood function
//  Calculate1DMCMCIntegration(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, "./mcmcdata.bangui.dat", nParameter1, vdX, vdY, vdYErr, bCalcErr);
  Calculate1DMCMCIntegrationOnTheFly(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, "./mcmcdata.bangui.dat", nParameter1, vdX, vdY);

  // Set the plot
  pwBanThis->m_Tab2_Plot_Results_1D.SetPlot(vdX.m_pnDimSize[0], vdX.m_pdData, vdY.m_pdData);

  // Set the errors, and write data to disk
  if(bCalcErr) {
    WritePlot("plotmcmcdata-1d.txt", vdX, vdY, vdYErr);
    pwBanThis->m_Tab2_Plot_Results_1D.SetYErr(vdX.m_pnDimSize[0], vdYErr.m_pdData);
  } else {
    WritePlot("plotmcmcdata-1d.txt", vdX, vdY);
  } // if bCalcErr

  // Check whether we have some parameter descriptions available
  pwBanThis->m_Tab2_Combo_Parameter1D.CopyItemText(
      pwBanThis->m_Tab2_Combo_Parameter1D.GetCurrentItem(),
      strParameter);
  pwBanThis->m_Tab2_Combo_Pulsars1D.CopyItemText(
      pwBanThis->m_Tab2_Combo_Pulsars1D.GetCurrentItem(), strPulsarName);
  pwBanThis->GetParameterDescriptions(
      strParameter, strPulsarName, strTitle, strXLabel);
  strcpy(strYLabel, "Likelihood / Max(Likelihood)");


  // Set the title and labels on the plot
  pwBanThis->m_Tab2_Plot_Results_1D.SetTitle(strTitle);
  pwBanThis->m_Tab2_Plot_Results_1D.SetXLabel(strXLabel);
  pwBanThis->m_Tab2_Plot_Results_1D.SetYLabel(strYLabel);

  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab2_Plot_Results_2D)) {
    pwBanThis->m_MainTab.RemoveWidget(&pwBanThis->m_Tab2_Plot_Results_2D);
  } // if HasWidget

  if(! pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab2_Plot_Results_1D)) {
    pwBanThis->m_MainTab.AddWidget(2, &pwBanThis->m_Tab2_Plot_Results_1D);
  } // if HasWidget
} // OnClickedIntegrate1D

/* This function is called when the Integrate3D button is clicked on the MCMC tab
 * */
void BanWindow::OnClickedIntegrate2D(CIWindow *pwThis) {
  int nParameter1, nParameter2;
  double dMax, dMin;
  double pdLevels[3];
  double dSigma1, dSigma2, dSigma3, dVolume, dCum;
  char strPulsarName[80], strXLabel[80], strYLabel[80],
       strTitle[80], strParameter1[80], strParameter2[80];
  int nStep;
  CVector vdX, vdY;
  CMatrix mdZ, mdTemp;
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  nParameter1 = pwBanThis->GetParameterFromPulsarSource(
      pwBanThis->m_Tab2_Combo_Pulsars2D.GetCurrentItem(),
      pwBanThis->m_Tab2_Combo_Parameter2D_1.GetCurrentItem());
  nParameter2 = pwBanThis->GetParameterFromPulsarSource(
      pwBanThis->m_Tab2_Combo_Pulsars2D.GetCurrentItem(),
      pwBanThis->m_Tab2_Combo_Parameter2D_2.GetCurrentItem());

//  Integrate3DMCMCData(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, "./mcmcdata.bangui.dat", nParameter1, nParameter2);
#if 0 // Make a skymap thingy
  if(nParameter1 == 2 && nParameter2 == 3)
    Calculate3DMCMCIntegrationSkymap(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, "./mcmcdata.bangui.dat", nParameter1, nParameter2, vdX, vdY, mdZ);
  else
#endif
  {
//    Calculate3DMCMCIntegration(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, "./mcmcdata.bangui.dat", nParameter1, nParameter2, vdX, vdY, mdZ);
    Calculate3DMCMCIntegrationOnTheFly(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, "./mcmcdata.bangui.dat", nParameter1, nParameter2, vdX, vdY, mdZ);
  } // if nParameter

  mdTemp = mdZ;
  dMin = double(mdZ[0][0]);
  dMax = dMin;
  for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
      mdTemp[i][j] = double(mdZ[i][j]);
      if(dMin > double(mdTemp[i][j]))
	dMin = double(mdTemp[i][j]);
      if(dMax < double(mdTemp[i][j]))
	dMax = double(mdTemp[i][j]);
    } // for j
  } // for i

  dVolume = 0;
  for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
      mdZ[i][j] = double(mdTemp[i][j]) / dMax;
      dVolume += double(mdZ[i][j]);
    } // for j
  } // for i

  mdTemp = mdZ;
  QuickSort(mdTemp.m_pdData, 0, mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1);

  dCum = 0;
  nStep = 0;
  for(int i=mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1; i >= 0; i--) {
    dCum += mdTemp.m_pdData[i];
    switch(nStep) {
      case 0:
	if(dCum > 0.68268949 * dVolume) {
	  nStep++;
	  dSigma1 = mdTemp.m_pdData[i];
	} // if
	break;
      case 1:
	if(dCum > 0.95449974 * dVolume) {
	  nStep++;
	  dSigma2 = mdTemp.m_pdData[i];
	} // if
	break;
      case 2:
	if(dCum > 0.99730024 * dVolume) {
	  nStep++;
	  dSigma3 = mdTemp.m_pdData[i];
	} // if
	break;
      default:
	break;
    } // switch
  } // for i
  printf("dSigma1: %f,  dSigma2: %f,  dSigma3: %f\n", dSigma1, dSigma2, dSigma3);

  pdLevels[0] = dSigma1;
  pdLevels[1] = dSigma2;
  pdLevels[2] = dSigma3;

  // Check whether we have some parameter descriptions available
  pwBanThis->m_Tab2_Combo_Parameter2D_1.CopyItemText(
      pwBanThis->m_Tab2_Combo_Parameter2D_1.GetCurrentItem(),
      strParameter1);
  pwBanThis->m_Tab2_Combo_Parameter2D_2.CopyItemText(
      pwBanThis->m_Tab2_Combo_Parameter2D_2.GetCurrentItem(),
      strParameter2);
  pwBanThis->m_Tab2_Combo_Pulsars2D.CopyItemText(
      pwBanThis->m_Tab2_Combo_Pulsars2D.GetCurrentItem(), strPulsarName);
  pwBanThis->GetParameterDescriptions(
      strParameter1, strPulsarName, strTitle, strXLabel);
  pwBanThis->GetParameterDescriptions(
      strParameter2, strPulsarName, strTitle, strYLabel);

  pwBanThis->m_Tab2_Plot_Results_2D.SetTitle(strTitle);
  pwBanThis->m_Tab2_Plot_Results_2D.SetXLabel(strXLabel);
  pwBanThis->m_Tab2_Plot_Results_2D.SetYLabel(strYLabel);

  pwBanThis->m_Tab2_Plot_Results_2D.SetPlot(mdZ.m_pnDimSize[0], mdZ.m_pnDimSize[1], vdX.m_pdData, vdY.m_pdData, mdZ.m_pdData);
  pwBanThis->m_Tab2_Plot_Results_2D.SetLevels(3, pdLevels);

  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab2_Plot_Results_1D)) {
    pwBanThis->m_MainTab.RemoveWidget(&pwBanThis->m_Tab2_Plot_Results_1D);
  } // if HasWidget

  if(! pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab2_Plot_Results_2D)) {
    pwBanThis->m_MainTab.AddWidget(2, &pwBanThis->m_Tab2_Plot_Results_2D);
  } // if HasWidget

#if 1
  CVector vdXX, vdYY, vdZZ;
  int nCount=0;
  vdXX.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);
  vdYY.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);
  vdZZ.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);

  for(int i=0; i<mdZ.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdZ.m_pnDimSize[1]; j++) {
      vdXX[nCount] = double(vdX[i]);
      vdYY[nCount] = double(vdY[j]);
      vdZZ[nCount] = double(mdZ[i][j]);
      nCount++;
    } // for j
  } // for i

  WritePlot("plotmcmcdata-2d.txt", vdXX, vdYY, vdZZ);
#endif
} // OnClickedIntegrate2D

/* This function is called when the SavePlot button is clicked on the MCMC tab
 * */
void BanWindow::OnClickedTab2SavePlot(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  char strDeviceName[160];
  char strTimeString[80];
  time_t tNow;
  struct tm *ptiTimeInfo;

  time(&tNow);
  ptiTimeInfo = localtime(&tNow);
  strftime(strTimeString, 80, "%Y%m%d%H%M%S", ptiTimeInfo);
  printf("Plotting with timestamp %s\n", strTimeString);

  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab2_Plot_Results_1D)) {
    sprintf(strDeviceName, "mcmc1d-%s.ps/ps", strTimeString);
    pwBanThis->m_Tab2_Plot_Results_1D.SavePlot(strDeviceName);
  } // if HasWidget

  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab2_Plot_Results_2D)) {
    sprintf(strDeviceName, "mcmc2d-%s.ps/ps", strTimeString);
    pwBanThis->m_Tab2_Plot_Results_2D.SavePlot(strDeviceName);
  } // if HasWidget
} // OnClickedTab2SavePlot

/* This function is called when the Evidence button is clicked on the Plot tab
 * */
void BanWindow::OnClickedTab2Evidence(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  CalculateMCMCEvidence(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, "./mcmcdata.bangui.dat");
} // OnClickedTab2Evidence

/* This function is called when the MakePlot button is clicked on the 1D Plot tab
 * */
void BanWindow::OnClickedTab3Plot2D(CIWindow *pwThis) {
  int nParameter1, nParameter2;
  double dMax, dMin;
  double pdLevels[3];
  double dSigma1, dSigma2, dSigma3, dVolume, dCum;
  char strPulsarName[80], strXLabel[80], strYLabel[80],
       strTitle[80], strParameter1[80], strParameter2[80];
  int nStep, nP1, nP2;
  CVector vdX, vdY;
  CMatrix mdZ, mdZML, mdZT2, mdTemp;
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  nParameter1 = pwBanThis->GetParameterFromPulsarSourceIncMar(
      pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem(),
      pwBanThis->m_Tab3_Combo_Parameter2D_1.GetCurrentItem());
  nParameter2 = pwBanThis->GetParameterFromPulsarSourceIncMar(
      pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem(),
      pwBanThis->m_Tab3_Combo_Parameter2D_2.GetCurrentItem());

  nP1 = -1;
  nP2 = -1;

  {
    Calculate3DMCMCIntegrationIncMar(pwBanThis->m_oData, pwBanThis->m_oParameters, pwBanThis->m_oConstants, "./mcmcdata.bangui.dat", nParameter1, nParameter2, nP1, nP2, vdX, vdY, mdZ, mdZML, mdZT2);
  } // if nParameter

  if(nP1 != -1) {
    // Parameter one is a marginalisation parameter
    for(int i=0; i<vdX.m_pnDimSize[0]; i++) {
//      printf("%e   %e   ", double(vdX[i]), pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem()].pdTempo2Multiplication[nP1]);
      vdX[i] = double(vdX[i])*pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem()].pdTempo2Multiplication[nP1];
//      vdX[i] = double(vdX[i]) + pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem()].pdTempo2Value[nP1];
//      printf("%e   ", double(vdX[i]));
    } // for i
  } // if nP1
  if(nP2 != -1) {
    // Parameter two is a marginalisation parameter
    for(int i=0; i<vdY.m_pnDimSize[0]; i++) {
      vdY[i] = double(vdY[i])*pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem()].pdTempo2Multiplication[nP2];
//      vdY[i] = double(vdY[i])+pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem()].pdTempo2Value[nP2];
    } // for i
  } // if nP2

  // Calculate the sigma's for mdZ
  mdTemp = mdZ;
  dMin = double(mdZ[0][0]);
  dMax = dMin;
  for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
      mdTemp[i][j] = double(mdZ[i][j]);
      if(dMin > double(mdTemp[i][j]))
	dMin = double(mdTemp[i][j]);
      if(dMax < double(mdTemp[i][j]))
	dMax = double(mdTemp[i][j]);
    } // for j
  } // for i

  dVolume = 0;
  for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
      mdZ[i][j] = double(mdTemp[i][j]) / dMax;
      dVolume += double(mdZ[i][j]);
    } // for j
  } // for i

  mdTemp = mdZ;
  QuickSort(mdTemp.m_pdData, 0, mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1);

  dCum = 0;
  nStep = 0;
  for(int i=mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1; i >= 0; i--) {
    dCum += mdTemp.m_pdData[i];
    switch(nStep) {
      case 0:
	if(dCum > 0.68268949 * dVolume) {
	  nStep++;
	  dSigma1 = mdTemp.m_pdData[i];
	} // if
	break;
      case 1:
	if(dCum > 0.95449974 * dVolume) {
	  nStep++;
	  dSigma2 = mdTemp.m_pdData[i];
	} // if
	break;
      case 2:
	if(dCum > 0.99730024 * dVolume) {
	  nStep++;
	  dSigma3 = mdTemp.m_pdData[i];
	} // if
	break;
      default:
	break;
    } // switch
  } // for i
  printf("dSigma1: %f,  dSigma2: %f,  dSigma3: %f\n", dSigma1, dSigma2, dSigma3);

  pdLevels[0] = dSigma1;
  pdLevels[1] = dSigma2;
  pdLevels[2] = dSigma3;

  // Check whether we have some parameter descriptions available
  pwBanThis->m_Tab3_Combo_Parameter2D_1.CopyItemText(
      pwBanThis->m_Tab3_Combo_Parameter2D_1.GetCurrentItem(),
      strParameter1);
  pwBanThis->m_Tab3_Combo_Parameter2D_2.CopyItemText(
      pwBanThis->m_Tab3_Combo_Parameter2D_2.GetCurrentItem(),
      strParameter2);
  pwBanThis->m_Tab3_Combo_Pulsars2D.CopyItemText(
      pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem(), strPulsarName);
  pwBanThis->GetParameterDescriptions(
      strParameter1, strPulsarName, strTitle, strXLabel);
  pwBanThis->GetParameterDescriptions(
      strParameter2, strPulsarName, strTitle, strYLabel);

  pwBanThis->m_Tab3_Plot_2D.SetTitle(strTitle);
  pwBanThis->m_Tab3_Plot_2D.SetXLabel(strXLabel);
  pwBanThis->m_Tab3_Plot_2D.SetYLabel(strYLabel);

  pwBanThis->m_Tab3_Plot_2D.SetPlot(mdZ.m_pnDimSize[0], mdZ.m_pnDimSize[1], vdX.m_pdData, vdY.m_pdData, mdZ.m_pdData);
  pwBanThis->m_Tab3_Plot_2D.SetLevels(3, pdLevels);


  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_1D)) {
    pwBanThis->m_MainTab.RemoveWidget(&pwBanThis->m_Tab3_Plot_1D);
  } // if HasWidget

  if(! pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_2D)) {
    pwBanThis->m_MainTab.AddWidget(3, &pwBanThis->m_Tab3_Plot_2D);
  } // if HasWidget

#if 1
  CVector vdXX, vdYY, vdZZ;
  int nCount=0;
  vdXX.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);
  vdYY.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);
  vdZZ.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);

  for(int i=0; i<mdZ.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdZ.m_pnDimSize[1]; j++) {
      vdXX[nCount] = double(vdX[i]);
      vdYY[nCount] = double(vdY[j]);
      vdZZ[nCount] = double(mdZ[i][j]);
      nCount++;
    } // for j
  } // for i

  WritePlot("plotmcmcdata-2d.txt", vdXX, vdYY, vdZZ);

  // Calculate the sigma's for mdZML
  if(mdZML.Defined()) {
    mdTemp = mdZML;
    dMin = double(mdZML[0][0]);
    dMax = dMin;
    for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
      for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
	mdTemp[i][j] = double(mdZML[i][j]);
	if(dMin > double(mdTemp[i][j]))
	  dMin = double(mdTemp[i][j]);
	if(dMax < double(mdTemp[i][j]))
	  dMax = double(mdTemp[i][j]);
      } // for j
    } // for i

    dVolume = 0;
    for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
      for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
	mdZML[i][j] = double(mdTemp[i][j]) / dMax;
	dVolume += double(mdZML[i][j]);
      } // for j
    } // for i

    mdTemp = mdZML;
    QuickSort(mdTemp.m_pdData, 0, mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1);

    dCum = 0;
    nStep = 0;
    for(int i=mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1; i >= 0; i--) {
      dCum += mdTemp.m_pdData[i];
      switch(nStep) {
	case 0:
	  if(dCum > 0.68268949 * dVolume) {
	    nStep++;
	    dSigma1 = mdTemp.m_pdData[i];
	  } // if
	  break;
	case 1:
	  if(dCum > 0.95449974 * dVolume) {
	    nStep++;
	    dSigma2 = mdTemp.m_pdData[i];
	  } // if
	  break;
	case 2:
	  if(dCum > 0.99730024 * dVolume) {
	    nStep++;
	    dSigma3 = mdTemp.m_pdData[i];
	  } // if
	  break;
	default:
	  break;
      } // switch
    } // for i
    printf("mdZML dSigma1: %f,  dSigma2: %f,  dSigma3: %f\n", dSigma1, dSigma2, dSigma3);

    nCount = 0;
    for(int i=0; i<mdZML.m_pnDimSize[0]; i++) {
      for(int j=0; j<mdZML.m_pnDimSize[1]; j++) {
	vdXX[nCount] = double(vdX[i]);
	vdYY[nCount] = double(vdY[j]);
	vdZZ[nCount] = double(mdZML[i][j]);
	nCount++;
      } // for j
    } // for i

    WritePlot("plotmcmcdata-2d-ML.txt", vdXX, vdYY, vdZZ);
  } // if Defined

  if(mdZT2.Defined()) {
    // Calculate the sigma's for mdZT2
    mdTemp = mdZT2;
    dMin = double(mdZT2[0][0]);
    dMax = dMin;
    for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
      for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
	mdTemp[i][j] = double(mdZT2[i][j]);
	if(dMin > double(mdTemp[i][j]))
	  dMin = double(mdTemp[i][j]);
	if(dMax < double(mdTemp[i][j]))
	  dMax = double(mdTemp[i][j]);
      } // for j
    } // for i

    dVolume = 0;
    for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
      for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
	mdZT2[i][j] = double(mdTemp[i][j]) / dMax;
	dVolume += double(mdZT2[i][j]);
      } // for j
    } // for i

    mdTemp = mdZT2;
    QuickSort(mdTemp.m_pdData, 0, mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1);

    dCum = 0;
    nStep = 0;
    for(int i=mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1; i >= 0; i--) {
      dCum += mdTemp.m_pdData[i];
      switch(nStep) {
	case 0:
	  if(dCum > 0.68268949 * dVolume) {
	    nStep++;
	    dSigma1 = mdTemp.m_pdData[i];
	  } // if
	  break;
	case 1:
	  if(dCum > 0.95449974 * dVolume) {
	    nStep++;
	    dSigma2 = mdTemp.m_pdData[i];
	  } // if
	  break;
	case 2:
	  if(dCum > 0.99730024 * dVolume) {
	    nStep++;
	    dSigma3 = mdTemp.m_pdData[i];
	  } // if
	  break;
	default:
	  break;
      } // switch
    } // for i
    printf("mdZT2 dSigma1: %f,  dSigma2: %f,  dSigma3: %f\n", dSigma1, dSigma2, dSigma3);

    nCount = 0;
    for(int i=0; i<mdZT2.m_pnDimSize[0]; i++) {
      for(int j=0; j<mdZT2.m_pnDimSize[1]; j++) {
	vdXX[nCount] = double(vdX[i]);
	vdYY[nCount] = double(vdY[j]);
	vdZZ[nCount] = double(mdZT2[i][j]);
	nCount++;
      } // for j
    } // for i

    WritePlot("plotmcmcdata-2d-T2.txt", vdXX, vdYY, vdZZ);
  } // if Defined
#endif
} // OnClickedTab3Plot2D

/* This function is called when the PlotP1 button is clicked on the 1D Plot tab
 * */
void BanWindow::OnClickedTab3PlotP1(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  CVector vdX, vdY;
  char strPulsarName[80], strXLabel[80], strYLabel[80],
       strTitle[80], strParameter[80];
  int nParameter;

  nParameter = pwBanThis->GetParameterFromPulsarSourceIncMar(
      pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem(),
      pwBanThis->m_Tab3_Combo_Parameter2D_1.GetCurrentItem());

  SetBlockMatrix(pwBanThis->m_oData, pwBanThis->m_oConstants);
  SetGeometricPart(pwBanThis->m_oData, pwBanThis->m_oConstants);
  pwBanThis->m_oData.mdC.Initialize(pwBanThis->m_oConstants.n, pwBanThis->m_oConstants.n);
  SetCoherenceMatrix(pwBanThis->m_oData, pwBanThis->m_oConstants, false);

#ifdef HAVE_MPI
  if(pwBanThis->m_nProcesses >= 2)
  {
    printf("Entering MPI-1DPlot mode...\n");
    printf("nProcessRank: %i of %i\n", pwBanThis->m_nProcessRank, pwBanThis->m_nProcesses);
  } else
#endif // HAVE_MPI
  {
    Calculate1DPlot(pwBanThis->m_oData,
      pwBanThis->m_oParameters,
      pwBanThis->m_oConstants, nParameter, vdX, vdY);

    pwBanThis->m_Tab3_Plot_1D.SetPlot(
	vdX.m_pnDimSize[0],
	vdX.m_pdData,
	vdY.m_pdData);
#if 1
    WritePlot("plotdata-1d.txt", vdX, vdY);
#endif

    // Check whether we have some parameter descriptions available
    pwBanThis->m_Tab3_Combo_Parameter2D_1.CopyItemText(
	pwBanThis->m_Tab3_Combo_Parameter2D_1.GetCurrentItem(),
	strParameter);
    pwBanThis->m_Tab3_Combo_Pulsars2D.CopyItemText(
	pwBanThis->m_Tab3_Combo_Pulsars2D.GetCurrentItem(), strPulsarName);
    pwBanThis->GetParameterDescriptions(
	strParameter, strPulsarName, strTitle, strXLabel);
    strcpy(strYLabel, "Likelihood / Max(Likelihood)");

    pwBanThis->m_Tab3_Plot_1D.SetTitle(strTitle);
    pwBanThis->m_Tab3_Plot_1D.SetXLabel(strXLabel);
    pwBanThis->m_Tab3_Plot_1D.SetYLabel(strYLabel);
  }

  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_2D)) {
    pwBanThis->m_MainTab.RemoveWidget(&pwBanThis->m_Tab3_Plot_2D);
  } // if HasWidget

  if(! pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_1D)) {
    pwBanThis->m_MainTab.AddWidget(3, &pwBanThis->m_Tab3_Plot_1D);
  } // if HasWidget
} // OnClickedTab3PlotP1

/* This function is called when the SavePlot button is clicked on the Plot tab
 * */
void BanWindow::OnClickedTab3SavePlot(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  char strDeviceName[160];
  char strTimeString[80];
  time_t tNow;
  struct tm *ptiTimeInfo;

  time(&tNow);
  ptiTimeInfo = localtime(&tNow);
  strftime(strTimeString, 80, "%Y%m%d%H%M%S", ptiTimeInfo);
  printf("Plotting with timestamp %s\n", strTimeString);

  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_1D)) {
    sprintf(strDeviceName, "plot1d-%s.ps/ps", strTimeString);
    pwBanThis->m_Tab3_Plot_1D.SavePlot(strDeviceName);
  } // if HasWidget
  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_2D)) {
    sprintf(strDeviceName, "plot2d-%s.ps/ps", strTimeString);
    pwBanThis->m_Tab3_Plot_2D.SavePlot(strDeviceName);
  } // if HasWidget
} // OnClickedTab3SavePlot


/* This function is called when the RunMCMCEnsemble button is clicked on tab 3
 * */
void BanWindow::OnClickedTab3RunMCMCEnsemble(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  FILE *pFile;
  int nMCMCSteps=0;

  try {
    pwBanThis->m_oData.mdC.Initialize(pwBanThis->m_oConstants.n, pwBanThis->m_oConstants.n);
    if(!pwBanThis->m_bPlugin)
      GeneratePulsarAngles(pwBanThis->m_oData, pwBanThis->m_oConstants);

    SetGeometricPart(pwBanThis->m_oData, pwBanThis->m_oConstants);
    SetCoherenceMatrix(pwBanThis->m_oData, pwBanThis->m_oConstants, true);

    // Check if mcmcensembledata-resume.dat already exists
    if(pFile = fopen("mcmcensemble-resume.dat", "r")) {
      fclose(pFile);
      // Yes, we should resume the calculation
      fprintf(stderr, "Using \"mcmcensemble-resume.dat\"\n");
      ReadMCMCEnsembleDataFileSets(pwBanThis->m_oData.nDataSets, nMCMCSteps, "mcmcensemble-resume.dat", pwBanThis->m_oData.mdDataSets);
    } else if(pFile = fopen("mcmcensemble-usedata.dat", "r")) {
      fclose(pFile);
      // Yes, we should use the data
      fprintf(stderr, "Using data of \"mcmcensemble-usedata.dat\"\n");
      ReadMCMCEnsembleDataFileSets(pwBanThis->m_oData.nDataSets, nMCMCSteps, "mcmcensemble-usedata.dat", pwBanThis->m_oData.mdDataSets);
      nMCMCSteps = 0;
    } else {
      // Generate datasets
//      GenerateEnsembleResiduals(pwBanThis->m_oData, pwBanThis->m_oConstants, pwBanThis->m_oConstants.nEnsembleDataSets);
      GenerateRadiometerEnsembleResiduals(pwBanThis->m_oData, pwBanThis->m_oConstants, pwBanThis->m_oConstants.nEnsembleDataSets);
    } // if pFile

    MCMCEnsembleImportanceSample(pwBanThis->m_oData, pwBanThis->m_oConstants, "mcmcdata.bangui.dat", nMCMCSteps);

    // Now save the residuals, and the MCMC
//    pwBanThis->m_bHasData = true;
//    pwBanThis->SetWidgetsContent();
//    pwBanThis->m_bSimulatedData = true;
  } catch (ELinearError err) {  // Error handling
    switch(err) {
    case ELENotDefined:
      printf("Not Defined\n");
      break;
    case ELEWrongClassType:
      printf("Wrong Class Type\n");
      break;
    case ELEBadIndex:
      printf("Wrong Index number\n");
      break;
    case ELEDimensionMisMatch:
      printf("Matrix Dimensions Do Not Match\n");
      break;
    case ELELapack:
      printf("Lapack error!\n");
      break;
    default:
      printf("Default error!\n");
      break;
    } // switch
  } // try
} // OnClickedTab3RunMCMCEnsemble

/* This function is called when the button 'down' for the DataSet selection is
 * clicked
 * */
void BanWindow::OnClickedTab3Down(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  int nCurrentItem;

  if(pwBanThis->m_oData.nDataSets > 0) {
    nCurrentItem = pwBanThis->m_Tab3_Combo_DataSet_Ensemble.GetCurrentItem();
    nCurrentItem--;
    if(nCurrentItem < 0) nCurrentItem += pwBanThis->m_oData.nDataSets;
    pwBanThis->m_Tab3_Combo_DataSet_Ensemble.SetCurrentItem(nCurrentItem);
  } // if nDataSets
} // OnClickedTab3Down

/* This function is called when the button 'up' for the DataSet selection is
 * clicked
 * */
void BanWindow::OnClickedTab3Up(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  int nCurrentItem;

  if(pwBanThis->m_oData.nDataSets > 0) {
    nCurrentItem = pwBanThis->m_Tab3_Combo_DataSet_Ensemble.GetCurrentItem();
    nCurrentItem++;
    if(nCurrentItem >= pwBanThis->m_oData.nDataSets) nCurrentItem -= pwBanThis->m_oData.nDataSets;
    pwBanThis->m_Tab3_Combo_DataSet_Ensemble.SetCurrentItem(nCurrentItem);
  } // if nDataSets
} // OnClickedTab3Up

/* This function is called when a pulsar is selected in the ensemble group
 * */
void BanWindow::OnSelectedTab3Pulsar2DEnsemble(CIWindow *pwThis) {
  ((BanWindow*)pwThis)->Set2DPlotParameterComboBoxes();
} // OnSelectedTab3Pulsar2DEnsemble


/* This function is called when the button EnsembleML is pressed in the ensemble group
 * */
void BanWindow::OnClickedTab3Plot2DMLEnsemble(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  int nParameter1, nParameter2;
  char strPulsarName[80], strXLabel[80], strYLabel[80],
       strTitle[80], strParameter1[80], strParameter2[80], strBuf[80];
  CVector vdX, vdY;

  // If we include marginalisation parameters, change these to
  // "GetParameterFromPulsarSourceIncMar".
  nParameter1 = pwBanThis->GetParameterFromPulsarSource(
      pwBanThis->m_Tab3_Combo_Pulsars2D_Ensemble.GetCurrentItem(),
      pwBanThis->m_Tab3_Combo_Parameter2D_1_Ensemble.GetCurrentItem());
  nParameter2 = pwBanThis->GetParameterFromPulsarSource(
      pwBanThis->m_Tab3_Combo_Pulsars2D_Ensemble.GetCurrentItem(),
      pwBanThis->m_Tab3_Combo_Parameter2D_2_Ensemble.GetCurrentItem());

  Calculate2DEnsembleMLPoints(pwBanThis->m_oData, pwBanThis->m_oConstants, "./mcmcensembledata.dat", nParameter1, nParameter2, vdX, vdY);

  // Set the number of datasets
  pwBanThis->m_Tab3_Combo_DataSet_Ensemble.SetItems(pwBanThis->m_oData.nDataSets);
  for(int i=0; i<pwBanThis->m_oData.nDataSets; i++) {
    sprintf(strBuf, "Set %i", i);
    pwBanThis->m_Tab3_Combo_DataSet_Ensemble.SetItem(i, strBuf);
  } // for i

  pwBanThis->m_Tab3_Combo_DataSet_Ensemble.SetItems(pwBanThis->m_oData.nDataSets);
  for(int i=0; i<pwBanThis->m_oData.nDataSets; i++) {
    sprintf(strBuf, "Set %i", i);
    pwBanThis->m_Tab3_Combo_DataSet_Ensemble.SetItem(i, strBuf);
  } // for i

  // Check whether we have some parameter descriptions available
  pwBanThis->m_Tab3_Combo_Parameter2D_1_Ensemble.CopyItemText(
      pwBanThis->m_Tab3_Combo_Parameter2D_1_Ensemble.GetCurrentItem(),
      strParameter1);
  pwBanThis->m_Tab3_Combo_Parameter2D_2_Ensemble.CopyItemText(
      pwBanThis->m_Tab3_Combo_Parameter2D_2_Ensemble.GetCurrentItem(),
      strParameter2);
  pwBanThis->m_Tab3_Combo_Pulsars2D_Ensemble.CopyItemText(
      pwBanThis->m_Tab3_Combo_Pulsars2D_Ensemble.GetCurrentItem(), strPulsarName);
  pwBanThis->GetParameterDescriptions(
      strParameter1, strPulsarName, strTitle, strXLabel);
  pwBanThis->GetParameterDescriptions(
      strParameter2, strPulsarName, strTitle, strYLabel);

  pwBanThis->m_Tab3_Plot_1D.SetTitle(strTitle);
  pwBanThis->m_Tab3_Plot_1D.SetXLabel(strXLabel);
  pwBanThis->m_Tab3_Plot_1D.SetYLabel(strYLabel);

  pwBanThis->m_Tab3_Plot_1D.SetPlot(vdX.m_pnDimSize[0], vdX.m_pdData, vdY.m_pdData);

  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_2D)) {
    pwBanThis->m_MainTab.RemoveWidget(&pwBanThis->m_Tab3_Plot_2D);
  } // if HasWidget
  if(! pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_1D)) {
    pwBanThis->m_MainTab.AddWidget(3, &pwBanThis->m_Tab3_Plot_1D);
  } // if HasWidget

#if 1
  WritePlot("plotmlensembledata-2d.txt", vdX, vdY);
#endif
} // OnClickedTab3Plot2DMLEnsemble

/* This function is called when the button Ensemble2D is pressed in the ensemble group
 * */
void BanWindow::OnClickedTab3Plot2DEnsemble(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  int nParameter1, nParameter2;
  double dMax, dMin;
  double pdLevels[3];
  double dSigma1, dSigma2, dSigma3, dVolume, dCum;
  char strPulsarName[80], strXLabel[80], strYLabel[80],
       strTitle[80], strParameter1[80], strParameter2[80], strBuf[80];
  int nStep, nDataSet;
  CVector vdX, vdY;
  CMatrix mdZ, mdTemp;

  int nMCMCSteps;

  // If we include marginalisation parameters, change these to
  // "GetParameterFromPulsarSourceIncMar".
  nParameter1 = pwBanThis->GetParameterFromPulsarSource(
      pwBanThis->m_Tab3_Combo_Pulsars2D_Ensemble.GetCurrentItem(),
      pwBanThis->m_Tab3_Combo_Parameter2D_1_Ensemble.GetCurrentItem());
  nParameter2 = pwBanThis->GetParameterFromPulsarSource(
      pwBanThis->m_Tab3_Combo_Pulsars2D_Ensemble.GetCurrentItem(),
      pwBanThis->m_Tab3_Combo_Parameter2D_2_Ensemble.GetCurrentItem());

  nDataSet = pwBanThis->m_Tab3_Combo_DataSet_Ensemble.GetCurrentItem();

  // Are we doing the MLDR variant, or regular importance re-sampling?
#if 0
  Calculate3DMCMCEnsembleIntegrationOnTheFly(pwBanThis->m_oData, pwBanThis->m_oConstants, "./mcmcensembledata.dat", nParameter1, nParameter2, nDataSet, vdX, vdY, mdZ);
#else
  Calculate3DMCMCEnsembleMLDRIntegrationOnTheFly(pwBanThis->m_oData, pwBanThis->m_oConstants, "./mcmcensembledata.dat", nParameter1, nParameter2, vdX, vdY, mdZ);
#endif

  pwBanThis->m_Tab3_Combo_DataSet_Ensemble.SetItems(pwBanThis->m_oData.nDataSets);
  for(int i=0; i<pwBanThis->m_oData.nDataSets; i++) {
    sprintf(strBuf, "Set %i", i);
    pwBanThis->m_Tab3_Combo_DataSet_Ensemble.SetItem(i, strBuf);
  } // for i

  mdTemp = mdZ;
  dMin = double(mdZ[0][0]);
  dMax = dMin;
  for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
      mdTemp[i][j] = double(mdZ[i][j]);
      if(dMin > double(mdTemp[i][j]))
	dMin = double(mdTemp[i][j]);
      if(dMax < double(mdTemp[i][j]))
	dMax = double(mdTemp[i][j]);
    } // for j
  } // for i

  dVolume = 0;
  for(int i=0; i<mdTemp.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdTemp.m_pnDimSize[1]; j++) {
      mdZ[i][j] = double(mdTemp[i][j]) / dMax;
      dVolume += double(mdZ[i][j]);
    } // for j
  } // for i

  mdTemp = mdZ;
  QuickSort(mdTemp.m_pdData, 0, mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1);

  dCum = 0;
  nStep = 0;
  for(int i=mdTemp.m_pnDimSize[0]*mdTemp.m_pnDimSize[1]-1; i >= 0; i--) {
    dCum += mdTemp.m_pdData[i];
    switch(nStep) {
      case 0:
	if(dCum > 0.68268949 * dVolume) {
	  nStep++;
	  dSigma1 = mdTemp.m_pdData[i];
	} // if
	break;
      case 1:
	if(dCum > 0.95449974 * dVolume) {
	  nStep++;
	  dSigma2 = mdTemp.m_pdData[i];
	} // if
	break;
      case 2:
	if(dCum > 0.99730024 * dVolume) {
	  nStep++;
	  dSigma3 = mdTemp.m_pdData[i];
	} // if
	break;
      default:
	break;
    } // switch
  } // for i
  printf("dSigma1: %f,  dSigma2: %f,  dSigma3: %f\n", dSigma1, dSigma2, dSigma3);

  pdLevels[0] = dSigma1;
  pdLevels[1] = dSigma2;
  pdLevels[2] = dSigma3;

  // Check whether we have some parameter descriptions available
  pwBanThis->m_Tab3_Combo_Parameter2D_1_Ensemble.CopyItemText(
      pwBanThis->m_Tab3_Combo_Parameter2D_1_Ensemble.GetCurrentItem(),
      strParameter1);
  pwBanThis->m_Tab3_Combo_Parameter2D_2_Ensemble.CopyItemText(
      pwBanThis->m_Tab3_Combo_Parameter2D_2_Ensemble.GetCurrentItem(),
      strParameter2);
  pwBanThis->m_Tab3_Combo_Pulsars2D_Ensemble.CopyItemText(
      pwBanThis->m_Tab3_Combo_Pulsars2D_Ensemble.GetCurrentItem(), strPulsarName);
  pwBanThis->GetParameterDescriptions(
      strParameter1, strPulsarName, strTitle, strXLabel);
  pwBanThis->GetParameterDescriptions(
      strParameter2, strPulsarName, strTitle, strYLabel);

  pwBanThis->m_Tab3_Plot_2D.SetTitle(strTitle);
  pwBanThis->m_Tab3_Plot_2D.SetXLabel(strXLabel);
  pwBanThis->m_Tab3_Plot_2D.SetYLabel(strYLabel);

  pwBanThis->m_Tab3_Plot_2D.SetPlot(mdZ.m_pnDimSize[0], mdZ.m_pnDimSize[1], vdX.m_pdData, vdY.m_pdData, mdZ.m_pdData);
  pwBanThis->m_Tab3_Plot_2D.SetLevels(3, pdLevels);


  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_1D)) {
    pwBanThis->m_MainTab.RemoveWidget(&pwBanThis->m_Tab3_Plot_1D);
  } // if HasWidget
  if(! pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab3_Plot_2D)) {
    pwBanThis->m_MainTab.AddWidget(3, &pwBanThis->m_Tab3_Plot_2D);
  } // if HasWidget

#if 1
  CVector vdXX, vdYY, vdZZ;
  int nCount=0;
  vdXX.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);
  vdYY.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);
  vdZZ.Initialize(mdZ.m_pnDimSize[0] * mdZ.m_pnDimSize[1]);

  for(int i=0; i<mdZ.m_pnDimSize[0]; i++) {
    for(int j=0; j<mdZ.m_pnDimSize[1]; j++) {
      vdXX[nCount] = double(vdX[i]);
      vdYY[nCount] = double(vdY[j]);
      vdZZ[nCount] = double(mdZ[i][j]);
      nCount++;
    } // for j
  } // for i

  WritePlot("plotmcmcensembledata-2d.txt", vdXX, vdYY, vdZZ);
#endif
} // OnClickedTab3Plot2DEnsemble

/* This function is called when the button is pressed
 * */
void BanWindow::OnClickedTab3Plot1DM(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  Calculate1DEnsembleIntegrationMarParameters(pwBanThis->m_oData, pwBanThis->m_oConstants, "mcmcensembledata.dat");
} // OnClickedTab3Plot1DM

/* This function is called when the button is pressed
 * */
void BanWindow::OnClickedTab3Plot1Dm(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  Calculate1DEnsembleIntegrationMCMCParameters(pwBanThis->m_oData, pwBanThis->m_oConstants, "mcmcensembledata.dat");
} // OnClickedTab3Plot1Dm

/* This function is called when an item in the Pulsar Combobox is selected in
 * tab 4
 * */
void BanWindow::OnSelectedTab4Pulsar(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  int nMaxPoints;
  double dTau;
  CVector vdSZ, vdT, vdX, vdY, vdYErr;

  nMaxPoints = int(log(pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].nObservations) / log(2.0)) - 2;
  nMaxPoints *= 5; // ************************************** // *************

  if(pwBanThis->m_bPlugin && pwBanThis->m_bHasData) {
    vdSZ.Initialize(nMaxPoints);
    vdT.Initialize(nMaxPoints);
    vdX.Initialize(pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].nObservations);
    vdY.Initialize(pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].nObservations);
    vdYErr.Initialize(pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].nObservations);
    for(int i=0; i<pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].nObservations; i++) {
      vdX[i] = pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].pdTOA[i];
      vdY[i] = pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].pdResiduals[i];
      vdYErr[i] = pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].pdDeltaResiduals[i];
    } // for i

    for(int i=0; i<nMaxPoints; i++) {
//      dTau = (pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].pdTOA[pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].nObservations-1] - pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].pdTOA[0]) / pow(2, nMaxPoints - i - 1);
      dTau = (pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].pdTOA[pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].nObservations-1] - pwBanThis->m_oConstants.poPulsars[pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].pdTOA[0]) / pow(2, 0.2*(nMaxPoints - i - 1)); // ***************************************** // *******************
      vdT[i] = dTau / SPERYEAR;
      vdSZ[i] = 0.0; //sigma_z(vdX, vdY, vdYErr, dTau);
    } // for i

    pwBanThis->m_Tab4_Plot_Sigmaz.SetPlot(vdT.m_pnDimSize[0], vdT.m_pdData, vdSZ.m_pdData);

    if(! pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab4_Plot_Sigmaz)) {
      pwBanThis->m_MainTab.AddWidget(4, &pwBanThis->m_Tab4_Plot_Sigmaz);
    } // if HasWidget
  } // if m_bPlugin && m_bHasData

  pwBanThis->SetWidgetsContentTab4();

#if 0 // Write sigma_z
    WritePlot("sigmaz.txt", vdT, vdSZ);
#endif // Write sigma_z
} // OnSelectedTab4Pulsar

/* This function is called when the SavePlot button is clicked on the residuals tab
 * */
void BanWindow::OnClickedTab4SavePlot(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  char strDeviceName[160];
  char strTimeString[80];
  time_t tNow;
  struct tm *ptiTimeInfo;

  time(&tNow);
  ptiTimeInfo = localtime(&tNow);
  strftime(strTimeString, 80, "%Y%m%d%H%M%S", ptiTimeInfo);
  printf("Plotting with timestamp %s\n", strTimeString);

  if(pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab4_Plot_Sigmaz)) {
    sprintf(strDeviceName, "sigmaz-%s-%s.ps/ps", pwBanThis->m_oConstants.poPulsars[
	pwBanThis->m_Tab4_Combo_Pulsars.GetCurrentItem()].strPulsarName, strTimeString);
    pwBanThis->m_Tab4_Plot_Sigmaz.SavePlot(strDeviceName);
  } // if HasWidget
} // OnClickedTab4SavePlot


/* This function is called when an item in the Curve Combobox is selected in
 * tab 5
 * */
void BanWindow::OnSelectedTab5Curve(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;

  if(pwBanThis->m_bPlugin && pwBanThis->m_bHasData && pwBanThis->m_bHave_Sensitivity_Curve) {

    pwBanThis->m_Tab5_Plot_Sensitivity.SetLogScale(0, true);
    switch(pwBanThis->m_Tab5_Combo_Curve.GetCurrentItem()) {
      case 0:
	pwBanThis->m_Tab5_Plot_Sensitivity.SetPlot(
	    pwBanThis->m_vdTab5_A.m_pnDimSize[0],
	    pwBanThis->m_vdTab5_A.m_pdData,
	    pwBanThis->m_vdTab5_DS.m_pdData);
	pwBanThis->m_Tab5_Plot_Sensitivity.DrawLine(false);
	pwBanThis->m_Tab5_Plot_Sensitivity.SetTitle("Sensitvity of Array");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetXLabel("GWB amplitude");
	//pwBanThis->m_Tab5_Plot_Sensitivity.SetYLabel("\\gs\\dz\\u");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetYLabel("Detection Significance");
	break;
      case 1:
	pwBanThis->m_Tab5_Plot_Sensitivity.SetPlot(
	    pwBanThis->m_vdTab5_A.m_pnDimSize[0],
	    pwBanThis->m_vdTab5_A.m_pdData,
	    pwBanThis->m_vdTab5_S.m_pdData);
	pwBanThis->m_Tab5_Plot_Sensitivity.DrawLine(true);
	pwBanThis->m_Tab5_Plot_Sensitivity.SetTitle("Sensitvity of Array");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetXLabel("GWB amplitude");
	//pwBanThis->m_Tab5_Plot_Sensitivity.SetYLabel("\\gs\\dz\\u");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetYLabel("Theoretical DS");
	break;
      case 2:
	pwBanThis->m_Tab5_Plot_Sensitivity.SetPlot(
	    pwBanThis->m_vdTab5_A.m_pnDimSize[0],
	    pwBanThis->m_vdTab5_A.m_pdData,
	    pwBanThis->m_vdTab5_AdA.m_pdData);
	pwBanThis->m_Tab5_Plot_Sensitivity.DrawLine(false);
	pwBanThis->m_Tab5_Plot_Sensitivity.SetTitle("Sensitvity of Array");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetXLabel("GWB amplitude");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetYLabel("Bayesian SNR");
	break;
      case 3:
	pwBanThis->m_Tab5_Plot_Sensitivity.SetPlot(
	    pwBanThis->m_vdTab5_A.m_pnDimSize[0],
	    pwBanThis->m_vdTab5_A.m_pdData,
	    pwBanThis->m_vdTab5_S.m_pdData);
	pwBanThis->m_Tab5_Plot_Sensitivity.SetYFit(
	    pwBanThis->m_vdTab5_A.m_pnDimSize[0],
	    pwBanThis->m_vdTab5_AdA.m_pdData);
	pwBanThis->m_Tab5_Plot_Sensitivity.DrawLine(true);
	pwBanThis->m_Tab5_Plot_Sensitivity.SetTitle("Sensitvity of Array");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetXLabel("GWB amplitude");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetYLabel("Theoretical DS");
	break;
      case 4:
	pwBanThis->m_Tab5_Plot_Sensitivity.SetPlot(
	    pwBanThis->m_vdTab5_A.m_pnDimSize[0],
	    pwBanThis->m_vdTab5_A.m_pdData,
	    pwBanThis->m_vdTab5_S.m_pdData);
	pwBanThis->m_Tab5_Plot_Sensitivity.SetYFit(
	    pwBanThis->m_vdTab5_A.m_pnDimSize[0],
	    pwBanThis->m_vdTab5_DS.m_pdData);
	pwBanThis->m_Tab5_Plot_Sensitivity.DrawLine(true);
	pwBanThis->m_Tab5_Plot_Sensitivity.SetTitle("Sensitvity of Array");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetXLabel("GWB amplitude");
	pwBanThis->m_Tab5_Plot_Sensitivity.SetYLabel("Theoretical DS");
	break;
      default:
	break;
    } // switch

    if(! pwBanThis->m_MainTab.HasWidget(&pwBanThis->m_Tab5_Plot_Sensitivity)) {
      pwBanThis->m_MainTab.AddWidget(5, &pwBanThis->m_Tab5_Plot_Sensitivity);
    } // if HasWidget
  } // if m_bPlugin && m_bHasData

  pwBanThis->SetWidgetsContentTab5();
} // OnSelectedTab5Curve


/* This function is called when the 'sense' button is clicked on the simulate tab
 * */
void BanWindow::OnClickedTab5Sensitivity(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
} // OnClickedTab5Simulate

/* This function is called when the Simulate button is clicked on the simulate tab
 * */
void BanWindow::OnClickedTab5Simulate(CIWindow *pwThis) {
  BanWindow *pwBanThis = (BanWindow *)pwThis;
  try {
    pwBanThis->m_oData.mdC.Initialize(pwBanThis->m_oConstants.n, pwBanThis->m_oConstants.n);
    if(!pwBanThis->m_bPlugin)
      GeneratePulsarAngles(pwBanThis->m_oData, pwBanThis->m_oConstants);

    SetGeometricPart(pwBanThis->m_oData, pwBanThis->m_oConstants);
    SetCoherenceMatrix(pwBanThis->m_oData, pwBanThis->m_oConstants, true);
    GenerateResiduals(pwBanThis->m_oData, pwBanThis->m_oConstants, false);
    SetBlockMatrix(pwBanThis->m_oData, pwBanThis->m_oConstants);
    pwBanThis->m_bHasData = true;
    pwBanThis->SetWidgetsContent();
    pwBanThis->m_bSimulatedData = true;
  } catch (ELinearError err) {  // Error handling
    switch(err) {
    case ELENotDefined:
      printf("Not Defined\n");
      break;
    case ELEWrongClassType:
      printf("Wrong Class Type\n");
      break;
    case ELEBadIndex:
      printf("Wrong Index number\n");
      break;
    case ELEDimensionMisMatch:
      printf("Matrix Dimensions Do Not Match\n");
      break;
    case ELELapack:
      printf("Lapack error!\n");
      break;
    default:
      printf("Default error!\n");
      break;
    } // switch
  } // try
} // OnClickedTab5Simulate

/* This function initializes mpi and stores the variables in the member
 * variables of the main  window
 * */
void BanWindow::InitMpi(int argc, char *argv[]) {
  int nNameLength;
  m_nProcesses = 1;
  m_nProcessRank = 0;
  strcpy(m_strProcessorName, "default");
} // InitMpi


/* This function registers an mpi datatype. This allows for communication
 * between different mpi processes.
 * */
void BanWindow::PrepareMPISend() {
  // Initialize the MPI datatypes for communication with the clients
#ifdef HAVE_MPI
  struct SMPIParametersType oMPIParameters;	// A representative variable for communication

  // length, displacement, and type arrays used to describe an MPI derived type
  // their size reflects the number of components in SMPIParametersType
  int          lena[3]; 			// Length of the members of SMPIParametersType
  MPI_Aint     loca[3];			// Relative memory location of the members 
  MPI_Datatype typa[3];			// Type indicator

  MPI_Aint     baseaddress;			// The base-address of the representation

  // a variable to hold the MPI type indicator for SMPIParametersType
//  MPI_Datatype m_MPI_Parameters;

  // set up the MPI description of SMPIParametersType
  MPI_Address(&oMPIParameters, &baseaddress);

  // oMPIParameters.pdPar has length of MAX_PARAMETERS doubles
  lena[0] = MAX_PARAMETERS;
  MPI_Address(&oMPIParameters.pdPar, &loca[0]); 
  loca[0] -= baseaddress;			// byte address relative to start of structure
  typa[0] = MPI_DOUBLE;

  // oMPIParameters.nLength has length of 1 integer
  lena[1] = 1;
  MPI_Address(&oMPIParameters.dLogLik, &loca[1]); 
  loca[1] -= baseaddress; 
  typa[1] = MPI_DOUBLE;

  // oMPIParameters.nLength has length of 1 integer
  lena[2] = 1;
  MPI_Address(&oMPIParameters.nStatus, &loca[2]); 
  loca[2] -= baseaddress; 
  typa[2] = MPI_INT;

  // Now define and commit the datatype 'm_MPI_Parameters'
  MPI_Type_struct(3, lena, loca, typa, &m_MPI_Parameters);
  MPI_Type_commit(&m_MPI_Parameters);
#endif // HAVE_MPI
} // PrepareMPISend


void BanWindow::SetDefaultValues(bool bPlugin) {
  strcpy(m_strParam, "kT");
  if(bPlugin) {
    sprintf(m_oOriginalConstants.strBaseDir, "%s/bayesian", getenv("TEMPO2"));
    sprintf(m_oConstants.strBaseDir, "%s/bayesian", getenv("TEMPO2"));
  }
  else
    strcpy(m_strBaseDir, g_strBaseDir);

  // Use the current directory from now on!
  sprintf(m_oOriginalConstants.strBaseDir, "./");
  sprintf(m_oConstants.strBaseDir, "./");
  // Use the current directory from now on!

  strcpy(m_strDataDir, "./");
  strcpy(m_strParametersConf, "");
  strcpy(m_strResidualsFile, "residuals.dat");
  strcpy(m_strAnglesFile, "angles.dat");

  // We do not have data yet
  m_bHasData = false;

  // We need to display this window...
  m_bNeedsUpdate = true;

  // We will not run the Plk plug-in ( < 0 )
  m_nPlkPulsar = -1;
} // SetDefaultValues

void BanWindow::ReadParametersFile() {
  // Read the constants set in the parameters file
  try {
    ReadGlobalConstants(m_oOriginalConstants, m_oOriginalData, m_strParametersConf, m_strResidualsFile, m_strAnglesFile, m_strDataDir);
//    ReadGlobalConstants(m_oConstants, m_oData, m_strParametersConf, m_strResidualsFile, m_strAnglesFile, m_strDataDir);
    m_oConstants = m_oOriginalConstants;
    m_oData = m_oOriginalData;
  } catch (ELinearError err) {  // Error handling
    switch(err) {
    case ELENotDefined:
      printf("Not Defined\n");
      break;
    case ELEWrongClassType:
      printf("Wrong Class Type\n");
      break;
    case ELEBadIndex:
      printf("Wrong Index number\n");
      break;
    case ELEDimensionMisMatch:
      printf("Matrix Dimensions Do Not Match\n");
      break;
    case ELELapack:
      printf("Lapack error!\n");
      break;
    default:
      printf("Default error!\n");
      break;
    } // switch
  } // try
} // ReadParametersFile

/* This function loads the pulsar timing residuals from a file
 * */
void BanWindow::LoadPulsarData() {
  ReadAngles(m_oData, m_oConstants);
  ReadResiduals(m_oData, m_oConstants, m_oParameters);
  SetBlockMatrix(m_oData, m_oConstants);
  SetGeometricPart(m_oData, m_oConstants);
  m_bHasData = true;
} // LoadResiduals

/* This function should be run after LoadTempoData
 * */
void BanWindow::SetMatrices(bool bSetWidgets) {
  SetBlockMatrix(m_oData, m_oConstants);
  SetGeometricPart(m_oData, m_oConstants);
  m_bHasData = true;

  if(bSetWidgets)
    SetWidgetsContent();
} // SetMatrices

/* This function resets the m_nPlkPulsar variable to -1
 * */
void BanWindow::PlkHasBeenRun() {
  m_nPlkPulsar = -1;
} // PlkHasBeenRun

/* This function returns the number of the pulsar that we need to run Plk on
 * */
int BanWindow::GetPlkPulsarNumber() {
  return m_nPlkPulsar;
} // PlkHasBeenRun

/* This function returns whether or not data has been simulated
 * */
bool BanWindow::HasSimulated() {
  bool bReturnValue = m_bSimulatedData;
  m_bSimulatedData = false;
  return bReturnValue;
} // HasSimulated


/* This function sets the correct residuals and errors in the plot of the first
 * tab
 * */
void BanWindow::SetResidualsPlot() {
  char strBuf[80];
  double *pdBuf = new double[m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].nObservations];
  for(int i=0; i<m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].nObservations; i++)
    pdBuf[i] = m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].pdTOA[i]/SPERYEAR;

  m_Tab0_Plot_Residuals.SetPlot(
      m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].nObservations,
//      m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].pdTOA,
      pdBuf,
      m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].pdResiduals);

  if(m_bPlugin) {
    for(int i=0; i<m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].nObservations; i++)
      pdBuf[i] = m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].pdDeltaResiduals[i]*2;

    m_Tab0_Plot_Residuals.SetYErr(
	m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].nObservations,
	pdBuf);
  } // if m_bPlugin

  sprintf(strBuf, "Residuals of pulsar %s",
      m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].strPulsarName);
  m_Tab0_Plot_Residuals.SetTitle(strBuf);
  m_Tab0_Plot_Residuals.SetXLabel("TOA (yr.)");
  m_Tab0_Plot_Residuals.SetYLabel("Timing residuals (sec.)");

# if 0
  {
    CVector vdX, vdY, vdYErr;
    vdX.Initialize(m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].nObservations);
    vdY.Initialize(m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].nObservations);
    vdYErr.Initialize(m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].nObservations);
    for(int i=0; i<m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].nObservations; i++) {
      vdX[i] = m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].pdTOA[i];
      vdY[i] = m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].pdResiduals[i];
      vdYErr[i] = m_oConstants.poPulsars[m_Tab0_Combo_Pulsars.GetCurrentItem()].pdDeltaResiduals[i];
    } // for i
    WritePlot("plotresiduals.txt", vdX, vdY, vdYErr);
  }
#endif

  delete[] pdBuf;
} // SetResidualsPlot

