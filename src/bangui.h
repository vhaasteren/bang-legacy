/* bangui.h -- iugui interface
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


#ifndef __BANGUI_H__
#define __BANGUI_H__

#include "config.h"

#include "iugui.h"
#include "linal.h"
#include "corefunctions.h"


// Define the descriptions for the parameters. These are hard-coded here since
// it has more to do with visualisation than calculations.
struct SParameterDescriptions {
  char strParameter[80];
  char strDescription[80];
  char strUnit[80];
  char strFilePiece[80];
};

const SParameterDescriptions oParameterDescriptions[] = {
  {"det-sin (0)", "Deterministic source", "", "det"},
  {"pow-gr (0)", "GWB amplitude", "yr\\u1/2\\d", "gwbamp"},
  {"pow-gr (1)", "GWB spectral index", "", "gwbsi"},
  {"pow-gr (2)", "GWB cut-off frequency", "yr\\u-1\\d", "gwbcof"},
  {"err-sin (0)", "EFAC value", "", "efac"},
  {"wit-sin (0)", "White noise amplitude", "sec", "wna"},
  {"exp-sin (0)", "Exponential red-noise amplitude", "sec", "erna"},
  {"exp-sin (1)", "Exponential red-noise redness", "yr", "ernr"},
  {"lor-sin (0)", "Lorentzian red-noise amplitude", "sec", "lrna"},
  {"lor-sin (1)", "Lorentzian red-noise redness", "sec\\u-1\\d", "lrnr"},
  {"pow-sin (0)", "Power-law red-noise amplitude", "yr\\u1/2\\d", "plrna"},
  {"pow-sin (1)", "Power-law red-noise spectral index", "", "plrnsi"},
  {"pow-sin (2)", "Power-law red-noise cut-off frequency", "yr\\u-1\\d", "plrncof"},
  {"tem", "tempo2 parameter", "", "tempo2"},
  {"mem", "gw", "", "gwmem"},
  {"con (0)", "Offset", "sec", "con"},
  {"lin (0)", "Linear", "", "lin"},
  {"qua (0)", "Quadratic", "sec\\u-1\\d", "qua"},
  {"qub (0)", "Qubic", "sec\\u-2\\d", "qub"},
  {"qrt (0)", "Quartic", "sec\\u-3\\d", "qrt"},
  {"qsd (0)", "QSD offset", "sec", "qsd0"},
  {"qsd (1)", "QSD linear", "", "qsd1"},
  {"qsd (2)", "QSD quadratic", "sec\\u-1\\d", "qsd2"},
  {"sin (0)", "Sine amplitude", "sec", "sin0"},
  {"sin (1)", "Sine frequency", "sec\\u-1\\d", "sin1"},
  {"sin (2)", "Sine phase", "", "sin2"},
  {"mem (0)", "GW-Memory time hits earth", "sec", "mem0"},
  {"mem (1)", "GW-Memory metric perturbation", "", "mem1"},
  {"mem (2)", "GW-Memory RAJ of source", "rad", "mem2"},
  {"mem (3)", "GW-Memory DECL of source", "rad", "mem3"},
  {"mem (4)", "GW-Memory polarisation", "rad", "mem4"},
  {"tem (0)", "Tempo2 parameter", "", "tem0"},
  {"tem (1)", "Tempo2 parameter", "", "tem1"},
  {"tem (2)", "Tempo2 parameter", "", "tem2"},
  {"tem (3)", "Tempo2 parameter", "", "tem3"},
  {"tem (4)", "Tempo2 parameter", "", "tem4"},
  {"tem (5)", "Tempo2 parameter", "", "tem5"},
  {"tem (6)", "Tempo2 parameter", "", "tem6"},
  {"tem (7)", "Tempo2 parameter", "", "tem7"},
  {"tem (8)", "Tempo2 parameter", "", "tem8"},
  {"tem (9)", "Tempo2 parameter", "", "tem9"},
  {"tem (10)", "Tempo2 parameter", "", "tem10"},
  {"tem (11)", "Tempo2 parameter", "", "tem11"},
  {"tem (12)", "Tempo2 parameter", "", "tem12"},
  {"tem (13)", "Tempo2 parameter", "", "tem13"},
  {"tem (14)", "Tempo2 parameter", "", "tem14"},
  {"tem (15)", "Tempo2 parameter", "", "tem15"},
  {"tem (16)", "Tempo2 parameter", "", "tem16"},
  {"tem (17)", "Tempo2 parameter", "", "tem17"},
  {"tem (18)", "Tempo2 parameter", "", "tem18"},
  {"tem (19)", "Tempo2 parameter", "", "tem19"},
  {"tem (20)", "Tempo2 parameter", "", "tem20"},
  {"tem (21)", "Tempo2 parameter", "", "tem21"},
  {"tem (22)", "Tempo2 parameter", "", "tem22"},
  {"end", "end", "end", "end"}
};

struct SDeterministicParameterDescriptions {
  char strParameter[80];
  char strDescription[80];
  char strUnit[80];
  char strFilePiece[80];
};

const SDeterministicParameterDescriptions oDeterministicParameterDescriptions[] = {
};


class BanWindow : public CIWindow {
public:
  BanWindow(int argc, char *argv[], bool bPlugin=false, const char *strParametersConf = "");
  ~BanWindow();

  void SetMatrices(bool bSetWidgets=true);
  int GetParameterFromPulsarSource(int nPulsar, int nParameter);
  int GetParameterFromPulsarSourceIncMar(int nPulsar, int nParameter);
  void PlkHasBeenRun();
  int GetPlkPulsarNumber();
  bool HasSimulated();

  void RunPluginMCMC();
  void RunPluginPlot1D(int nParameter1);
  void RunPluginPlot2D(int nParameter1, int nParameter2);
  void RunPluginEnsemble();
  void RunPlugin1DEnsembleMCMC();
  void RunPlugin1DEnsembleMar();
  void RunReduce();
  void RunAffine();
  void RunPlInterpolation();

  // Memory management stuff
  void AllocatePulsars(int nPsr, int *pnPsrObs); // Allocate memory for all pulsars
  void AddSource(); // Allocate memory for extra source, and copy all

  // The ban variables
  SDataType m_oData;
  SConstantsType m_oConstants;
  SParametersType m_oParameters;
  SDataType m_oOriginalData;
  SConstantsType m_oOriginalConstants;
  SParametersType m_oOriginalParameters;

#ifdef HAVE_MPI
  MPI_Datatype m_MPI_Parameters;
#endif

  // MPI member variables
  int m_nProcesses, m_nProcessRank;
#ifdef HAVE_MPI
  char m_strProcessorName[MPI_MAX_PROCESSOR_NAME];
#else
  char m_strProcessorName[80];
#endif

protected:
  void InitMpi(int argc, char *argv[]);
  void PrepareMPISend();

  void SetDefaultValues(bool bPlugin);
  void ReadParametersFile();

  void InitializeWidgets();
  void SetWidgetsContent();

  void InitializeWidgetsTab0();
  void InitializeWidgetsTab1();
  void InitializeWidgetsTab2();
  void InitializeWidgetsTab3();
  void InitializeWidgetsTab4();
  void InitializeWidgetsTab5();
  void SetWidgetsContentTab0();
  void SetWidgetsContentTab1();
  void SetWidgetsContentTab2();
  void SetWidgetsContentTab3();
  void SetWidgetsContentTab4();
  void SetWidgetsContentTab5();

  void SetModelWidgets();

  void SetResidualsPlot();
  void Set1DMCMCParameterComboBoxes();
  void Set2DMCMCParameterComboBoxes();
  void Set2DPlotParameterComboBoxes();

  void LoadPulsarData();

  void GetParameterDescriptions(const char strParameter[], const char strPulsar[], char *strTitle, char *strLabel);
  void GetParameterDescriptionsIncMar(const char strParameter[], const char strPulsar[], char *strTitle, char *strLabel);
  void GetParameterCode(const SSource &oSource, int nParameter, char *strCode);
protected:
  // Control member variables
  bool m_bHasData;		// Whether we have data loaded
  bool m_bPlugin;		// Whether we are in a tempo2 plug-in
  int m_nPlkPulsar;		// What pulsar to run plk on (default = -1, which means don't)
  bool m_bSimulatedData;	// Whether we just simulated some dataset

  // File/Directory names
  char m_strParam[160];
  char m_strBaseDir[160];
  char m_strDataDir[160];
  char m_strParametersConf[160];
  char m_strResidualsFile[160];
  char m_strAnglesFile[160];

protected:
  // The signal handlers for the window
  static void OnClickedQuit(CIWindow *pwThis);
  static void OnClickedLoad(CIWindow *pwThis);
  static void OnClickedPrint(CIWindow *pwThis);
  static void OnPulsarTab0Selected(CIWindow *pwThis);
  static void OnPulsarTab1Selected(CIWindow *pwThis);

  static void OnFitTab0Selected(CIWindow *pwThis);
  static void OnSelectedTab0Plk(CIWindow *pwThis);
  static void OnClickedTab0MemMap(CIWindow *pwThis);
  static void OnClickedTab0SavePlot(CIWindow *pwThis);

  static void OnTab1SelectedStochasticParameter(CIWindow *pwThis);
  static void OnTab1ClickedStochasticCheckbox(CIWindow *pwThis);
  static void OnTab1ClickedStochasticDelete(CIWindow *pwThis);

  static void OnTab1SelectedDeterministicParameter(CIWindow *pwThis);
  static void OnTab1ClickedDeterministicCheckbox(CIWindow *pwThis);
  static void OnTab1ClickedDeterministicDelete(CIWindow *pwThis);
  static void OnClickedTab1PowellOptimum(CIWindow *pwThis);

  static void OnClickedRunMCMC(CIWindow *pwThis);
  static void OnClickedRunPlot1D(CIWindow *pwThis, int n1);
  static void OnClickedRunPlot2D(CIWindow *pwThis, int n1, int n2);
  static void OnClickedIntegrate1D(CIWindow *pwThis);
  static void OnClickedIntegrate2D(CIWindow *pwThis);
  static void OnClickedTab2SavePlot(CIWindow *pwThis);
  static void OnClickedTab2Evidence(CIWindow *pwThis);
  static void OnSelectedTab2Pulsar1D(CIWindow *pwThis);
  static void OnSelectedTab2Pulsar2D(CIWindow *pwThis);

  static void OnClickedTab3Plot2D(CIWindow *pwThis);
  static void OnClickedTab3PlotP1(CIWindow *pwThis);
  static void OnClickedTab3SavePlot(CIWindow *pwThis);
  static void OnSelectedTab3Pulsar2D(CIWindow *pwThis);

  static void OnSelectedTab3Pulsar2DEnsemble(CIWindow *pwThis);
  static void OnClickedTab3Plot2DMLEnsemble(CIWindow *pwThis);
  static void OnClickedTab3Plot2DEnsemble(CIWindow *pwThis);
  static void OnClickedTab3Plot1DM(CIWindow *pwThis);
  static void OnClickedTab3Plot1Dm(CIWindow *pwThis);
  static void OnClickedTab3RunMCMCEnsemble(CIWindow *pwThis);
  static void OnClickedTab3Down(CIWindow *pwThis);
  static void OnClickedTab3Up(CIWindow *pwThis);

  static void OnSelectedTab4Pulsar(CIWindow *pwThis);
  static void OnClickedTab4SavePlot(CIWindow *pwThis);

  static void OnSelectedTab5Curve(CIWindow *pwThis);
  static void OnClickedTab5Sensitivity(CIWindow *pwThis);
  static void OnClickedTab5Simulate(CIWindow *pwThis);

  // The main window widgets
  CITab m_MainTab;
  CIButton m_ButtonLoad;
  CIButton m_ButtonQuit;
  CIButton m_ButtonPrint;
  CIText m_Text_Copyright;

  // Tab 0 (Pulsar tab)
  CIComboBox m_Tab0_Combo_Pulsars;
  CIText m_Tab0_Text_Pulsars;
  CIText m_Tab0_Text_Datapoints_Description;
  CIText m_Tab0_Text_Datapoints_Value;
  CIText m_Tab0_Text_Timespan_Description;
  CIText m_Tab0_Text_Timespan_Value;
  CIText m_Tab0_Text_RMS_Description;
  CIText m_Tab0_Text_RMS_Value;
  CIText m_Tab0_Text_Chisq_Description;
  CIText m_Tab0_Text_Chisq_Value;
  CIText m_Tab0_Text_RedChisq_Description;
  CIText m_Tab0_Text_RedChisq_Value;

  CIButton m_Tab0_Button_MemMap;
  CIText m_Tab0_Text_MemMap_Description;
  CIButton m_Tab0_Button_SavePlot;
  CIButton m_Tab0_Button_Fit;
  CIText m_Tab0_Text_Fit_Description;
  CIButton m_Tab0_Button_Plk;
  CIText m_Tab0_Text_Plk_Description;
  CIPlot m_Tab0_Plot_Residuals;

  // Tab 1 (Model tab)
  CIComboBox m_Tab1_Combo_Pulsars;

  // Tab 1 Stochastic part
  CIText m_Tab1_Text_Stochastic;
  CITab m_Tab1_StochasticTab;
  CIComboBox m_Tab1_Combo_StochasticParameters[MAX_PARAMETERS_PER_SOURCE];
  CIText m_Tab1_Text_StochasticParameters[MAX_PARAMETERS_PER_SOURCE];
  CIText m_Tab1_Text_StochasticTag[MAX_PARAMETERS_PER_SOURCE];
  CIButton m_Tab1_Button_StochasticDelete[MAX_PARAMETERS_PER_SOURCE];

  CIText m_Tab1_Text_StochasticFit_Description[MAX_PARAMETERS_PER_SOURCE];
  CICheckbox m_Tab1_Check_StochasticFit[MAX_PARAMETERS_PER_SOURCE];
  CIText m_Tab1_Text_StochasticMarginalise_Description[MAX_PARAMETERS_PER_SOURCE];
  CICheckbox m_Tab1_Check_StochasticMarginalise[MAX_PARAMETERS_PER_SOURCE];

  // Tab 1 Deterministic part
  CIText m_Tab1_Text_Deterministic;
  CITab m_Tab1_DeterministicTab;
  CIComboBox m_Tab1_Combo_DeterministicParameters[MAX_PARAMETERS_PER_SOURCE];
  CIText m_Tab1_Text_DeterministicParameters[MAX_PARAMETERS_PER_SOURCE];
  CIText m_Tab1_Text_DeterministicTag[MAX_PARAMETERS_PER_SOURCE];
  CIButton m_Tab1_Button_DeterministicDelete[MAX_PARAMETERS_PER_SOURCE];

  CIText m_Tab1_Text_DeterministicFit_Description[MAX_PARAMETERS_PER_SOURCE];
  CICheckbox m_Tab1_Check_DeterministicFit[MAX_PARAMETERS_PER_SOURCE];
  CIText m_Tab1_Text_DeterministicMarginalise_Description[MAX_PARAMETERS_PER_SOURCE];
  CICheckbox m_Tab1_Check_DeterministicMarginalise[MAX_PARAMETERS_PER_SOURCE];

  // Tab 1, calculations
  CIButton m_Tab1_Button_PowellOptimum;

  // Tab 2 (MCMC tab)
  CIButton m_Tab2_Button_RunMCMC;

  CIGroup m_Tab2_Group_1D;
  CIText m_Tab2_Text_Pulsars1D_Description;
  CIComboBox m_Tab2_Combo_Pulsars1D;
  CIText m_Tab2_Text_Parameters1D_Description;
  CIComboBox m_Tab2_Combo_Parameter1D;
  CIText m_Tab2_Text_CalcErr_Description;
  CICheckbox m_Tab2_Check_CalcErr;
  CIButton m_Tab2_Button_Integrate1D;

  CIGroup m_Tab2_Group_2D;
  CIText m_Tab2_Text_Pulsars2D_Description;
  CIComboBox m_Tab2_Combo_Pulsars2D;
  CIText m_Tab2_Text_Parameters2D_1_Description;
  CIComboBox m_Tab2_Combo_Parameter2D_1;
  CIText m_Tab2_Text_Parameters2D_2_Description;
  CIComboBox m_Tab2_Combo_Parameter2D_2;
  CIButton m_Tab2_Button_Integrate2D;

  CIButton m_Tab2_Button_SavePlot;
  CIButton m_Tab2_Button_Evidence;

  CIPlot m_Tab2_Plot_Results_1D;
  CIContourPlot m_Tab2_Plot_Results_2D;

  // Tab 3 (MCMC Adv. tab)
  CIGroup m_Tab3_Group_2D;
  CIText m_Tab3_Text_Pulsars2D_Description;
  CIComboBox m_Tab3_Combo_Pulsars2D;
  CIText m_Tab3_Text_Parameters2D_1_Description;
  CIComboBox m_Tab3_Combo_Parameter2D_1;
  CIText m_Tab3_Text_Parameters2D_2_Description;
  CIComboBox m_Tab3_Combo_Parameter2D_2;
  CIButton m_Tab3_Button_Plot2D;
  CIButton m_Tab3_Button_PlotP1;

  CIButton m_Tab3_Button_SavePlot;
  CIPlot m_Tab3_Plot_1D;
  CIContourPlot m_Tab3_Plot_2D;

  CIGroup m_Tab3_Group_2D_Ensemble;
  CIText m_Tab3_Text_Pulsars2D_Ensemble_Description;
  CIComboBox m_Tab3_Combo_Pulsars2D_Ensemble;
  CIText m_Tab3_Text_Parameters2D_1_Ensemble_Description;
  CIComboBox m_Tab3_Combo_Parameter2D_1_Ensemble;
  CIText m_Tab3_Text_Parameters2D_2_Ensemble_Description;
  CIComboBox m_Tab3_Combo_Parameter2D_2_Ensemble;
  CIButton m_Tab3_Button_Plot2D_MLEnsemble;
  CIButton m_Tab3_Button_Plot2D_Ensemble;
  CIButton m_Tab3_Button_Plot1D_M;
  CIButton m_Tab3_Button_Plot1D_M2;

  CIComboBox m_Tab3_Combo_DataSet_Ensemble;
  CIButton m_Tab3_Button_Down;
  CIButton m_Tab3_Button_Up;
  CIButton m_Tab3_Button_RunMCMCEnsemble;

  // Tab 4 (ML tab)
  CIText m_Tab4_Text_Pulsars;
  CIComboBox m_Tab4_Combo_Pulsars;
  CIButton m_Tab4_Button_SavePlot;
  CIPlot m_Tab4_Plot_Sigmaz;

  // Tab 5 (Simulation tab)
  CIButton m_Tab5_Button_Sensitivity;
  CIButton m_Tab5_Button_Simulate;

  CIComboBox m_Tab5_Combo_Curve;
  CIPlot m_Tab5_Plot_Sensitivity;

  bool m_bHave_Sensitivity_Curve;
  CVector m_vdTab5_A;
  CVector m_vdTab5_S;
  CVector m_vdTab5_DS;
  CVector m_vdTab5_AdA;
};

#endif // __BANGUI_H__

