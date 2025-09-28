//Jinfeng, 25.09.25
//Fitter to calculate the uncertainty from the selection of lifetime variable [M(JJ)]

// H1. T Series
#include "TString.h"
#include "TH2.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFile.h"
#include "TROOT.h"
// H2. C Series
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
// H3. RooFit pack
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooKeysPdf.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooFitResult.h"
#include "RooStepFunction.h"
#include "RooGenericPdf.h"
#include "TLatex.h"
#include "RooCrystalBall.h"
#include "RooFFTConvPdf.h"
#include "RooHist.h"

using namespace RooFit;

void Fit_Final(float FourMuonMass_LowCut_Value, float FourMuonMass_HighCut_Value, int MassNumber, float Standard_Value, TString Output_Name) {

	//Shape Parameters
	//For prompt J/psi on the ctau dimension

	float ctau_JJ_Gaus_Mean_Value = 0.00003;
	float ctau_JJ_Gaus_Sigma1_Value = 0.00227;
	float ctau_JJ_Gaus_Sigma2_Value = 0.0074;
	float ctau_JJ_Gaus_Frac_Value = 0.882;

	//For non-prompt J/psi on the ctau dimension
	float ctau_JJ_Convolution_Tau_Value = -27.30;
	float ctau_JJ_Convolution_Mean_Value = 0.00036;
	float ctau_JJ_Convolution_Sigma_Value = 0.00292;

	//For JMuMu on the ctau(J/psi) dimension 
	float ctau_JMuMu_J_Convolution_Tau_Value = -26.7;
	float ctau_JMuMu_J_Convolution_Mean_Value = -0.0032;
	float ctau_JMuMu_J_Convolution_Sigma_Value = 0.0032;
	float ctau_JMuMu_J_Gaus_Mean_Value = 0.00062;
	float ctau_JMuMu_J_Gaus_Sigma_Value = 0.00098;
	float ctau_JMuMu_J_Frac_Value = 0.086;

	//For JMuMu on the ctau(MuMu) dimension 
	float ctau_JMuMu_MuMu_Convolution_Tau_Value = -22.3;
	float ctau_JMuMu_MuMu_Convolution_Mean_Value = -0.0065;
	float ctau_JMuMu_MuMu_Convolution_Sigma_Value = 0.0102;

	//For prompt J/psi on the SigLxy dimension
	float SigLxy_JJ_Convolution_Tau_Value = -1.328;
	float SigLxy_JJ_Convolution_Mean_Value = 0.237;
	float SigLxy_JJ_Convolution_Sigma_Value = 0.191;

	//For non-prompt J/psi on the SigLxy dimension
	float SigLxy_JJ_Exp_Tau_Value = -0.07056;

	//For JMuMu on the SigLxy(J/psi) dimension 
	float SigLxy_JMuMu_J_Convolution_Tau_Value = -1.08;
	float SigLxy_JMuMu_J_Convolution_Mean_Value = 0.2591;
	float SigLxy_JMuMu_J_Convolution_Sigma_Value = 0.002;
	float SigLxy_JMuMu_J_Exp_Tau_Value = -0.0608;
	float SigLxy_JMuMu_J_Frac_Value = 0.191;

	//For JMuMu on the SigLxy(MuMu) dimension 
	float SigLxy_JMuMu_MuMu_Exp_Tau_Value = -0.0625;

	//For Jpsi on the mass dimension
	float J_Mean_Value = 3.09312;
	float J_Sigma1_Value = 0.0240;
	float J_Sigma2_Value = 0.059;
	float J_Alpha_Value = 1.72;
	float J_N_Value = 3.00;
	float J_Frac_Value = 0.743;

	//Mass Windows
	float LowCut = 2.95;
	float HighCut = 3.25;

	//Bin cut
	float FourMuonMass_LowCut = FourMuonMass_LowCut_Value;
	float FourMuonMass_HighCut = FourMuonMass_HighCut_Value;

	//Particle Massage
	const float J_Mass = 3.0969;

	//Set the variables
	RooRealVar J1Mass_("J1Mass_", "m(J/#psi_{1}) [GeV]", LowCut, HighCut);
	RooRealVar J2Mass_("J2Mass_", "m(J/#psi_{2}) [GeV]", LowCut, HighCut);
	RooRealVar x1_("x1_", "c#tau^{J/#psi1} [cm]", -0.05, 0.28);
	RooRealVar x2_("x2_", "c#tau^{J/#psi2} [cm]", -0.05, 0.28);
	RooRealVar LxyPVSig1_("LxyPVSig1_", "Significance of L_{xy}^{J/#psi1}PV ", 0, 60);
	RooRealVar LxyPVSig2_("LxyPVSig2_", "Significance of L_{xy}^{J/#psi2}PV ", 0, 60);

	RooRealVar FourMuonMass_("FourMuonMass_", "m(J/#psi_{1}J/#psi_{2}) [Gev]", 0, 107.5);

	RooRealVar Weight_sum("Weight_sum", "Weight_sum", +0, 5000);

	//1D PDF
	//Mass - Parameters - J/Psi 
	RooRealVar J_Mean("J_Mean", "J_Mean", J_Mean_Value);  //1

	RooRealVar J_Sigma1("J_Sigma1", "J_Sigma1", J_Sigma1_Value);                //2
	RooRealVar J_Sigma2("J_Sigma2", "J_Sigma2", J_Sigma2_Value);                //3

	RooRealVar J_Alpha("J_Alpha", "J_Alpha", J_Alpha_Value);                            //4
	RooRealVar J_N("J_N", "J_N", J_N_Value);                                        //5

	RooRealVar J_Frac("J_Frac", "J_Frac", J_Frac_Value);                             //6

	//Mass1 - Function - CBJ1 - double CB
	RooAbsPdf* CBJ1_1 = new RooCBShape("CBJ1_1", "CBJ1_1", J1Mass_, J_Mean, J_Sigma1, J_Alpha, J_N);
	RooAbsPdf* CBJ1_2 = new RooCBShape("CBJ1_2", "CBJ1_2", J1Mass_, J_Mean, J_Sigma2, J_Alpha, J_N);

	RooAddPdf CBJ1("CBJ1", "CBJ1", RooArgList(*CBJ1_1, *CBJ1_2), J_Frac);

	//Mass2 - Function - CBJ2 - double CB
	RooAbsPdf* CBJ2_1 = new RooCBShape("CBJ2_1", "CBJ2_1", J2Mass_, J_Mean, J_Sigma1, J_Alpha, J_N);
	RooAbsPdf* CBJ2_2 = new RooCBShape("CBJ2_2", "CBJ2_2", J2Mass_, J_Mean, J_Sigma2, J_Alpha, J_N);

	RooAddPdf CBJ2("CBJ2", "CBJ2", RooArgList(*CBJ2_1, *CBJ2_2), J_Frac);

	//Mass - Parameters - MuMu
	//RooRealVar J_Cheb_Co1("J_Cheb_Co1", "J_Cheb_Co1", J_Cheb_Co1_Value);
	//RooRealVar J_Cheb_Co2("J_Cheb_Co2", "J_Cheb_Co2", J_Cheb_Co2_Value);

	RooRealVar J_Cheb_Co1("J_Cheb_Co1", "J_Cheb_Co1", -2.00, -3.00, 0.10);
	RooRealVar J_Cheb_Co2("J_Cheb_Co2", "J_Cheb_Co2", 0.40, -0.20, 1.00);

	//Mass1 - Function - J1_Cheb - second order cheb
	RooAbsPdf* J1_Cheb = new RooChebychev("J1_Cheb", "J1_Cheb", J1Mass_, RooArgList(J_Cheb_Co1, J_Cheb_Co2));

	//Mass2 - Function - J2_Cheb - second order cheb
	RooAbsPdf* J2_Cheb = new RooChebychev("J2_Cheb", "J2_Cheb", J2Mass_, RooArgList(J_Cheb_Co1, J_Cheb_Co2));

	//ctau - Parameters - Prompt JJ
	RooRealVar ctau_JJ_Gaus_Mean("ctau_JJ_Gaus_Mean", "ctau_JJ_Gaus_Mean", ctau_JJ_Gaus_Mean_Value);

	RooRealVar ctau_JJ_Gaus_Sigma1("ctau_JJ_Gaus_Sigma1", "ctau_JJ_Gaus_Sigma1", ctau_JJ_Gaus_Sigma1_Value);
	RooRealVar ctau_JJ_Gaus_Sigma2("ctau_JJ_Gaus_Sigma2", "ctau_JJ_Gaus_Sigma2", ctau_JJ_Gaus_Sigma2_Value);

	RooRealVar ctau_JJ_Gaus_Frac("ctau_JJ_Gaus_Frac", "ctau_JJ_Gaus_Frac", ctau_JJ_Gaus_Frac_Value);

	//ctau1 - Function - ctau1_JJ_Gaus - double gaus
	RooAbsPdf* ctau1_JJ_Gaus1 = new RooGaussian("ctau1_JJ_Gaus1", "ctau1_JJ_Gaus1", x1_, ctau_JJ_Gaus_Mean, ctau_JJ_Gaus_Sigma1);
	RooAbsPdf* ctau1_JJ_Gaus2 = new RooGaussian("ctau1_JJ_Gaus2", "ctau1_JJ_Gaus2", x1_, ctau_JJ_Gaus_Mean, ctau_JJ_Gaus_Sigma2);

	RooAddPdf ctau1_JJ_Gaus("ctau1_JJ_Gaus", "ctau1_JJ_Gaus", RooArgList(*ctau1_JJ_Gaus1, *ctau1_JJ_Gaus2), ctau_JJ_Gaus_Frac);

	//ctau2 - Function - ctau1_JJ_Gaus - double gaus
	RooAbsPdf* ctau2_JJ_Gaus1 = new RooGaussian("ctau2_JJ_Gaus1", "ctau2_JJ_Gaus1", x2_, ctau_JJ_Gaus_Mean, ctau_JJ_Gaus_Sigma1);
	RooAbsPdf* ctau2_JJ_Gaus2 = new RooGaussian("ctau2_JJ_Gaus2", "ctau2_JJ_Gaus2", x2_, ctau_JJ_Gaus_Mean, ctau_JJ_Gaus_Sigma2);

	RooAddPdf ctau2_JJ_Gaus("ctau2_JJ_Gaus", "ctau2_JJ_Gaus", RooArgList(*ctau2_JJ_Gaus1, *ctau2_JJ_Gaus2), ctau_JJ_Gaus_Frac);

	//ctau - Parameters - Non-prompt JJ
	RooRealVar ctau_JJ_Convolution_Tau("ctau_JJ_Convolution_Tau", "ctau_JJ_Convolution_Tau", ctau_JJ_Convolution_Tau_Value);

	RooRealVar ctau_JJ_Convolution_Mean("ctau_JJ_Convolution_Mean", "ctau_JJ_Convolution_Mean", ctau_JJ_Convolution_Mean_Value);
	RooRealVar ctau_JJ_Convolution_Sigma("ctau_JJ_Convolution_Sigma", "ctau_JJ_Convolution_Sigma", ctau_JJ_Convolution_Sigma_Value);

	RooRealVar ctau_StepValue("ctau_StepValue", "ctau_StepValue", 1);
	RooRealVar ctau_LowLimit("ctau_LowLimit", "ctau_LowLimit", 0);
	RooRealVar ctau_HighLimit("ctau_HighLimit", "ctau_HighLimit", 0.28);

	RooAbsReal* ctau1_Step = new RooStepFunction("ctau1_Step", "ctau1_Step", x1_, RooArgList(ctau_StepValue), RooArgList(ctau_LowLimit, ctau_HighLimit));
	RooAbsReal* ctau2_Step = new RooStepFunction("ctau2_Step", "ctau2_Step", x2_, RooArgList(ctau_StepValue), RooArgList(ctau_LowLimit, ctau_HighLimit));

	//ctau1 - Function - ctau1_JJ_Convolution - Convolution
	RooAbsPdf* ctau1_JJ_Exp = new RooGenericPdf("ctau1_JJ_Exp", "ctau1_JJ_Exp", "ctau1_Step * exp(x1_ * ctau_JJ_Convolution_Tau)", RooArgSet(*ctau1_Step, x1_, ctau_JJ_Convolution_Tau));
	RooAbsPdf* ctau1_JJ_Gaus_ForConvolution = new RooGaussian("ctau1_JJ_Gaus_ForConvolution", "ctau1_JJ_Gaus_ForConvolution", x1_, ctau_JJ_Convolution_Mean, ctau_JJ_Convolution_Sigma);

	RooFFTConvPdf ctau1_JJ_Convolution("ctau1_JJ_Convolution", "ctau1_JJ_Convolution", x1_, *ctau1_JJ_Exp, *ctau1_JJ_Gaus_ForConvolution);

	//ctau2 - Function - ctau2_JJ_Convolution - Convolution
	RooAbsPdf* ctau2_JJ_Exp = new RooGenericPdf("ctau2_JJ_Exp", "ctau2_JJ_Exp", "ctau2_Step * exp(x2_ * ctau_JJ_Convolution_Tau)", RooArgSet(*ctau2_Step, x2_, ctau_JJ_Convolution_Tau));
	RooAbsPdf* ctau2_JJ_Gaus_ForConvolution = new RooGaussian("ctau2_JJ_Gaus_ForConvolution", "ctau2_JJ_Gaus_ForConvolution", x2_, ctau_JJ_Convolution_Mean, ctau_JJ_Convolution_Sigma);

	RooFFTConvPdf ctau2_JJ_Convolution("ctau2_JJ_Convolution", "ctau2_JJ_Convolution", x2_, *ctau2_JJ_Exp, *ctau2_JJ_Gaus_ForConvolution);

	//ctau - Parameters - JMuMu - J/psi
	RooRealVar ctau_JMuMu_J_Convolution_Tau("ctau_JMuMu_J_Convolution_Tau", "ctau_JMuMu_J_Convolution_Tau", ctau_JMuMu_J_Convolution_Tau_Value);

	RooRealVar ctau_JMuMu_J_Convolution_Mean("ctau_JMuMu_J_Convolution_Mean", "ctau_JMuMu_J_Convolution_Mean", ctau_JMuMu_J_Convolution_Mean_Value);
	RooRealVar ctau_JMuMu_J_Convolution_Sigma("ctau_JMuMu_J_Convolution_Sigma", "ctau_JMuMu_J_Convolution_Sigma", ctau_JMuMu_J_Convolution_Sigma_Value);

	RooRealVar ctau_JMuMu_J_Gaus_Mean("ctau_JMuMu_J_Gaus_Mean", "ctau_JMuMu_J_Gaus_Mean", ctau_JMuMu_J_Gaus_Mean_Value);
	RooRealVar ctau_JMuMu_J_Gaus_Sigma("ctau_JMuMu_J_Gaus_Sigma", "ctau_JMuMu_J_Gaus_Sigma", ctau_JMuMu_J_Gaus_Sigma_Value);

	RooRealVar ctau_JMuMu_J_Frac("ctau_JMuMu_J_Frac", "ctau_JMuMu_J_Frac", ctau_JMuMu_J_Frac_Value);

	//ctau1 - Function - ctau1_JMuMu_J_PDF - Convolution+Gaus
	RooAbsPdf* ctau1_JMuMu_J_Exp = new RooGenericPdf("ctau1_JMuMu_J_Exp", "ctau1_JMuMu_J_Exp", "ctau1_Step * exp(x1_ * ctau_JMuMu_J_Convolution_Tau)", RooArgSet(*ctau1_Step, x1_, ctau_JMuMu_J_Convolution_Tau));
	RooAbsPdf* ctau1_JMuMu_J_Gaus_ForConvolution = new RooGaussian("ctau1_JMuMu_J_Gaus_ForConvolution", "ctau1_JMuMu_J_Gaus_ForConvolution", x1_, ctau_JMuMu_J_Convolution_Mean, ctau_JMuMu_J_Convolution_Sigma);

	RooFFTConvPdf ctau1_JMuMu_J_Convolution("ctau1_JMuMu_J_Convolution", "ctau1_JMuMu_J_Convolution", x1_, *ctau1_JMuMu_J_Exp, *ctau1_JMuMu_J_Gaus_ForConvolution);

	RooAbsPdf* ctau1_JMuMu_J_Gaus = new RooGaussian("ctau1_JMuMu_J_Gaus", "ctau1_JMuMu_J_Gaus2", x1_, ctau_JMuMu_J_Gaus_Mean, ctau_JMuMu_J_Gaus_Sigma);
	RooAddPdf ctau1_JMuMu_PDF("ctau1_JMuMu_PDF", "ctau1_JMuMu_PDF", RooArgList(*ctau1_JMuMu_J_Gaus, ctau1_JMuMu_J_Convolution), ctau_JMuMu_J_Frac);

	//ctau2 - Function - ctau2_JMuMu_J_PDF - Convolution+Gaus
	RooAbsPdf* ctau2_JMuMu_J_Exp = new RooGenericPdf("ctau2_JMuMu_J_Exp", "ctau2_JMuMu_J_Exp", "ctau2_Step * exp(x2_ * ctau_JMuMu_J_Convolution_Tau)", RooArgSet(*ctau2_Step, x2_, ctau_JMuMu_J_Convolution_Tau));
	RooAbsPdf* ctau2_JMuMu_J_Gaus_ForConvolution = new RooGaussian("ctau2_JMuMu_J_Gaus_ForConvolution", "ctau2_JMuMu_J_Gaus_ForConvolution", x2_, ctau_JMuMu_J_Convolution_Mean, ctau_JMuMu_J_Convolution_Sigma);

	RooFFTConvPdf ctau2_JMuMu_J_Convolution("ctau2_JMuMu_J_Convolution", "ctau2_JMuMu_J_Convolution", x2_, *ctau2_JMuMu_J_Exp, *ctau2_JMuMu_J_Gaus_ForConvolution);

	RooAbsPdf* ctau2_JMuMu_J_Gaus = new RooGaussian("ctau2_JMuMu_J_Gaus", "ctau2_JMuMu_J_Gaus2", x2_, ctau_JMuMu_J_Gaus_Mean, ctau_JMuMu_J_Gaus_Sigma);
	RooAddPdf ctau2_JMuMu_PDF("ctau2_JMuMu_PDF", "ctau2_JMuMu_PDF", RooArgList(*ctau2_JMuMu_J_Gaus, ctau2_JMuMu_J_Convolution), ctau_JMuMu_J_Frac);

	//ctau - Parameters - JMuMu - MuMu
	RooRealVar ctau_JMuMu_MuMu_Convolution_Tau("ctau_JMuMu_MuMu_Convolution_Tau", "ctau_JMuMu_MuMu_Convolution_Tau", ctau_JMuMu_MuMu_Convolution_Tau_Value);

	RooRealVar ctau_JMuMu_MuMu_Convolution_Mean("ctau_JMuMu_MuMu_Convolution_Mean", "ctau_JMuMu_MuMu_Convolution_Mean", ctau_JMuMu_MuMu_Convolution_Mean_Value);
	RooRealVar ctau_JMuMu_MuMu_Convolution_Sigma("ctau_JMuMu_MuMu_Convolution_Sigma", "ctau_JMuMu_MuMu_Convolution_Sigma", ctau_JMuMu_MuMu_Convolution_Sigma_Value);

	//ctau1 - Function - ctau1_JMuMu_Convolution - Convolution
	RooAbsPdf* ctau1_JMuMu_MuMu_Exp = new RooGenericPdf("ctau1_JMuMu_MuMu_Exp", "ctau1_JMuMu_MuMu_Exp", "ctau1_Step * exp(x1_ * ctau_JMuMu_MuMu_Convolution_Tau)", RooArgSet(*ctau1_Step, x1_, ctau_JMuMu_MuMu_Convolution_Tau));
	RooAbsPdf* ctau1_JMuMu_MuMu_Gaus_ForConvolution = new RooGaussian("ctau1_JMuMu_MuMu_Gaus_ForConvolution", "ctau1_JMuMu_MuMu_Gaus_ForConvolution", x1_, ctau_JMuMu_MuMu_Convolution_Mean, ctau_JMuMu_MuMu_Convolution_Sigma);

	RooFFTConvPdf ctau1_JMuMu_Convolution("ctau1_JMuMu_Convolution", "ctau1_JMuMu_Convolution", x1_, *ctau1_JMuMu_MuMu_Exp, *ctau1_JMuMu_MuMu_Gaus_ForConvolution);

	//ctau2 - Function - ctau2_JMuMu_Convolution - Convolution
	RooAbsPdf* ctau2_JMuMu_MuMu_Exp = new RooGenericPdf("ctau2_JMuMu_MuMu_Exp", "ctau2_JMuMu_MuMu_Exp", "ctau2_Step * exp(x2_ * ctau_JMuMu_MuMu_Convolution_Tau)", RooArgSet(*ctau2_Step, x2_, ctau_JMuMu_MuMu_Convolution_Tau));
	RooAbsPdf* ctau2_JMuMu_MuMu_Gaus_ForConvolution = new RooGaussian("ctau2_JMuMu_MuMu_Gaus_ForConvolution", "ctau2_JMuMu_MuMu_Gaus_ForConvolution", x2_, ctau_JMuMu_MuMu_Convolution_Mean, ctau_JMuMu_MuMu_Convolution_Sigma);

	RooFFTConvPdf ctau2_JMuMu_Convolution("ctau2_JMuMu_Convolution", "ctau2_JMuMu_Convolution", x2_, *ctau2_JMuMu_MuMu_Exp, *ctau2_JMuMu_MuMu_Gaus_ForConvolution);

	//SigLxy - Step 
	RooRealVar LxySig_StepValue("LxySig_StepValue", "LxySig_StepValue", 1);
	RooRealVar LxySig_LowLimit("LxySig_LowLimit", "LxySig_LowLimit", 0);
	RooRealVar LxySig_HighLimit_Narrow("LxySig_HighLimit_Narrow", "LxySig_HighLimit_Narrow", 6);
	RooRealVar LxySig_HighLimit_Wide("LxySig_HighLimit_Wide", "LxySig_HighLimit_Wide", 60);

	RooAbsReal* LxySig1_Step_Narrow = new RooStepFunction("LxySig1_Step_Narrow", "LxySig1_Step_Narrow", LxyPVSig1_, RooArgList(LxySig_StepValue), RooArgList(LxySig_LowLimit, LxySig_HighLimit_Narrow));
	RooAbsReal* LxySig2_Step_Narrow = new RooStepFunction("LxySig2_Step_Narrow", "LxySig2_Step_Narrow", LxyPVSig2_, RooArgList(LxySig_StepValue), RooArgList(LxySig_LowLimit, LxySig_HighLimit_Narrow));
	RooAbsReal* LxySig1_Step_Wide = new RooStepFunction("LxySig1_Step_Wide", "LxySig1_Step_Wide", LxyPVSig1_, RooArgList(LxySig_StepValue), RooArgList(LxySig_LowLimit, LxySig_HighLimit_Wide));
	RooAbsReal* LxySig2_Step_Wide = new RooStepFunction("LxySig2_Step_Wide", "LxySig2_Step_Wide", LxyPVSig2_, RooArgList(LxySig_StepValue), RooArgList(LxySig_LowLimit, LxySig_HighLimit_Wide));

	//SigLxy - Parameters - Prompt JJ
	RooRealVar SigLxy_JJ_Convolution_Tau("SigLxy_JJ_Convolution_Tau", "SigLxy_JJ_Convolution_Tau", SigLxy_JJ_Convolution_Tau_Value);

	RooRealVar SigLxy_JJ_Convolution_Mean("SigLxy_JJ_Convolution_Mean", "SigLxy_JJ_Convolution_Mean", SigLxy_JJ_Convolution_Mean_Value);
	RooRealVar SigLxy_JJ_Convolution_Sigma("SigLxy_JJ_Convolution_Sigma", "SigLxy_JJ_Convolution_Sigma", SigLxy_JJ_Convolution_Sigma_Value);

	//SigLxy1 - Function - SigLxy1_JJ_Convolution - Convolution
	RooAbsPdf* LxySig1_JJ_ExpForConvolution = new RooGenericPdf("LxySig1_JJ_ExpForConvolution", "LxySig1_JJ_ExpForConvolution", "LxySig1_Step_Narrow * exp(LxyPVSig1_ * SigLxy_JJ_Convolution_Tau)", RooArgSet(*LxySig1_Step_Narrow, LxyPVSig1_, SigLxy_JJ_Convolution_Tau));
	RooAbsPdf* LxySig1_JJ_Gaus = new RooGaussian("LxySig1_JJ_Gaus", "LxySig1_JJ_Gaus", LxyPVSig1_, SigLxy_JJ_Convolution_Mean, SigLxy_JJ_Convolution_Sigma);

	RooFFTConvPdf LxySig1_JJ_Convolution("LxySig1_JJ_Convolution", "LxySig1_JJ_Convolution", LxyPVSig1_, *LxySig1_JJ_ExpForConvolution, *LxySig1_JJ_Gaus);

	//SigLxy2 - Function - SigLxy2_JJ_Convolution - Convolution
	RooAbsPdf* LxySig2_JJ_ExpForConvolution = new RooGenericPdf("LxySig2_JJ_ExpForConvolution", "LxySig2_JJ_ExpForConvolution", "LxySig2_Step_Narrow * exp(LxyPVSig2_ * SigLxy_JJ_Convolution_Tau)", RooArgSet(*LxySig2_Step_Narrow, LxyPVSig2_, SigLxy_JJ_Convolution_Tau));
	RooAbsPdf* LxySig2_JJ_Gaus = new RooGaussian("LxySig2_JJ_Gaus", "LxySig2_JJ_Gaus", LxyPVSig2_, SigLxy_JJ_Convolution_Mean, SigLxy_JJ_Convolution_Sigma);

	RooFFTConvPdf LxySig2_JJ_Convolution("LxySig2_JJ_Convolution", "LxySig2_JJ_Convolution", LxyPVSig2_, *LxySig2_JJ_ExpForConvolution, *LxySig2_JJ_Gaus);

	//SigLxy - Parameters - Non-prompt JJ
	RooRealVar SigLxy_JJ_Exp_Tau("SigLxy_JJ_Exp_Tau", "SigLxy_JJ_Exp_Tau", SigLxy_JJ_Exp_Tau_Value);

	//SigLxy1 - Function - SigLxy1_JJ_Exp - Exponent
	RooAbsPdf* LxySig1_JJ_Exp = new RooGenericPdf("LxySig1_JJ_Exp", "LxySig1_JJ_Exp", "LxySig1_Step_Wide * exp(LxyPVSig1_ * SigLxy_JJ_Exp_Tau)", RooArgSet(*LxySig1_Step_Wide, LxyPVSig1_, SigLxy_JJ_Exp_Tau));

	//SigLxy2 - Function - SigLxy2_JJ_Exp - Exponent
	RooAbsPdf* LxySig2_JJ_Exp = new RooGenericPdf("LxySig2_JJ_Exp", "LxySig2_JJ_Exp", "LxySig2_Step_Wide * exp(LxyPVSig2_ * SigLxy_JJ_Exp_Tau)", RooArgSet(*LxySig2_Step_Wide, LxyPVSig2_, SigLxy_JJ_Exp_Tau));

	//SigLxy - Parameters - JMuMu - J/psi
	RooRealVar SigLxy_JMuMu_J_Convolution_Tau("SigLxy_JMuMu_J_Convolution_Tau", "SigLxy_JMuMu_J_Convolution_Tau", SigLxy_JMuMu_J_Convolution_Tau_Value);
	RooRealVar SigLxy_JMuMu_J_Convolution_Mean("SigLxy_JMuMu_J_Convolution_Mean", "SigLxy_JMuMu_J_Convolution_Mean", SigLxy_JMuMu_J_Convolution_Mean_Value);
	RooRealVar SigLxy_JMuMu_J_Convolution_Sigma("SigLxy_JMuMu_J_Convolution_Sigma", "SigLxy_JMuMu_J_Convolution_Sigma", SigLxy_JMuMu_J_Convolution_Sigma_Value);

	RooRealVar SigLxy_JMuMu_J_Exp_Tau("SigLxy_JMuMu_J_Exp_Tau", "SigLxy_JMuMu_J_Exp_Tau", SigLxy_JMuMu_J_Exp_Tau_Value);

	RooRealVar SigLxy_JMuMu_J_Frac("SigLxy_JMuMu_J_Frac", "SigLxy_JMuMu_J_Frac", SigLxy_JMuMu_J_Frac_Value);

	//SigLxy1 - Function - SigLxy1_JMuMu_PDF - Convolution+Exp
	RooAbsPdf* LxySig1_JMuMu_ExpForConvolution = new RooGenericPdf("LxySig1_JMuMu_ExpForConvolution", "LxySig1_JMuMu_ExpForConvolution", "LxySig1_Step_Narrow * exp(LxyPVSig1_ * SigLxy_JMuMu_J_Convolution_Tau)", RooArgSet(*LxySig1_Step_Narrow, LxyPVSig1_, SigLxy_JMuMu_J_Convolution_Tau));
	RooAbsPdf* LxySig1_JMuMu_Gaus = new RooGaussian("LxySig1_JMuMu_Gaus", "LxySig1_JMuMu_Gaus", LxyPVSig1_, SigLxy_JMuMu_J_Convolution_Mean, SigLxy_JMuMu_J_Convolution_Sigma);

	RooFFTConvPdf LxySig1_JMuMu_Convolution("LxySig1_JMuMu_Convolution", "LxySig1_JMuMu_Convolution", LxyPVSig1_, *LxySig1_JMuMu_ExpForConvolution, *LxySig1_JMuMu_Gaus);

	RooAbsPdf* LxySig1_JMuMu_Exp = new RooGenericPdf("LxySig1_JMuMu_Exp", "LxySig1_JMuMu_Exp", "LxySig1_Step_Wide * exp(LxyPVSig1_ * SigLxy_JMuMu_J_Exp_Tau)", RooArgSet(*LxySig1_Step_Wide, LxyPVSig1_, SigLxy_JMuMu_J_Exp_Tau));

	RooAddPdf LxySig1_JMuMu_PDF("LxySig1_JMuMu_PDF", "LxySig1_JMuMu_PDF", RooArgList(LxySig1_JMuMu_Convolution, *LxySig1_JMuMu_Exp), SigLxy_JMuMu_J_Frac);

	//SigLxy2 - Function - SigLxy2_JMuMu_PDF - Convolution+Exp
	RooAbsPdf* LxySig2_JMuMu_ExpForConvolution = new RooGenericPdf("LxySig2_JMuMu_ExpForConvolution", "LxySig2_JMuMu_ExpForConvolution", "LxySig2_Step_Narrow * exp(LxyPVSig2_ * SigLxy_JMuMu_J_Convolution_Tau)", RooArgSet(*LxySig2_Step_Narrow, LxyPVSig2_, SigLxy_JMuMu_J_Convolution_Tau));
	RooAbsPdf* LxySig2_JMuMu_Gaus = new RooGaussian("LxySig2_JMuMu_Gaus", "LxySig2_JMuMu_Gaus", LxyPVSig2_, SigLxy_JMuMu_J_Convolution_Mean, SigLxy_JMuMu_J_Convolution_Sigma);

	RooFFTConvPdf LxySig2_JMuMu_Convolution("LxySig2_JMuMu_Convolution", "LxySig2_JMuMu_Convolution", LxyPVSig2_, *LxySig2_JMuMu_ExpForConvolution, *LxySig2_JMuMu_Gaus);

	RooAbsPdf* LxySig2_JMuMu_Exp = new RooGenericPdf("LxySig2_JMuMu_Exp", "LxySig2_JMuMu_Exp", "LxySig2_Step_Wide * exp(LxyPVSig2_ * SigLxy_JMuMu_J_Exp_Tau)", RooArgSet(*LxySig2_Step_Wide, LxyPVSig2_, SigLxy_JMuMu_J_Exp_Tau));

	RooAddPdf LxySig2_JMuMu_PDF("LxySig2_JMuMu_PDF", "LxySig2_JMuMu_PDF", RooArgList(LxySig2_JMuMu_Convolution, *LxySig2_JMuMu_Exp), SigLxy_JMuMu_J_Frac);

	//SigLxy - Parameters - JMuMu - J/psi
	RooRealVar SigLxy_JMuMu_MuMu_Exp_Tau("SigLxy_JMuMu_MuMu_Exp_Tau", "SigLxy_JMuMu_MuMu_Exp_Tau", SigLxy_JMuMu_MuMu_Exp_Tau_Value);

	//SigLxy1 - Function - SigLxy1_JMuMu_Exp - Exponent
	RooAbsPdf* LxySig1_JMuMu_MuMu_Exp = new RooGenericPdf("LxySig1_JMuMu_MuMu_Exp", "LxySig1_JMuMu_MuMu_Exp", "LxySig1_Step_Wide * exp(LxyPVSig1_ * SigLxy_JMuMu_MuMu_Exp_Tau)", RooArgSet(*LxySig1_Step_Wide, LxyPVSig1_, SigLxy_JMuMu_MuMu_Exp_Tau));

	//SigLxy2 - Function - SigLxy2_JMuMu_Exp - Exponent
	RooAbsPdf* LxySig2_JMuMu_MuMu_Exp = new RooGenericPdf("LxySig2_JMuMu_MuMu_Exp", "LxySig2_JMuMu_MuMu_Exp", "LxySig2_Step_Wide * exp(LxyPVSig2_ * SigLxy_JMuMu_MuMu_Exp_Tau)", RooArgSet(*LxySig2_Step_Wide, LxyPVSig2_, SigLxy_JMuMu_MuMu_Exp_Tau));

	//4D PDF
	RooProdPdf PDF_JJ_PP("PDF_JJ_PP", "PDF_JJ_PP", RooArgList(CBJ1, CBJ2, LxySig1_JJ_Convolution, LxySig2_JJ_Convolution));

	RooProdPdf PDF_JJ_NPP("PDF_JJ_NPP", "PDF_JJ_NPP", RooArgList(CBJ1, CBJ2, *LxySig1_JJ_Exp, LxySig2_JJ_Convolution));
	RooProdPdf PDF_JJ_PNP("PDF_JJ_PNP", "PDF_JJ_PNP", RooArgList(CBJ1, CBJ2, LxySig1_JJ_Convolution, *LxySig2_JJ_Exp));
	RooProdPdf PDF_JJ_NPNP("PDF_JJ_NPNP", "PDF_JJ_NPNP", RooArgList(CBJ1, CBJ2, *LxySig1_JJ_Exp, *LxySig2_JJ_Exp));

	RooProdPdf PDF_J1MuMu("PDF_J1MuMu", "PDF_J1MuMu", RooArgList(CBJ1, *J2_Cheb, LxySig1_JMuMu_PDF, *LxySig2_JMuMu_MuMu_Exp));
	RooProdPdf PDF_MuMuJ2("PDF_MuMuJ2", "PDF_MuMuJ2", RooArgList(*J1_Cheb, CBJ2, *LxySig1_JMuMu_MuMu_Exp, LxySig2_JMuMu_PDF));
	RooProdPdf PDF_MuMuMuMu("PDF_MuMuMuMu", "PDF_MuMuMuMu", RooArgList(*J1_Cheb, *J2_Cheb, *LxySig1_JMuMu_MuMu_Exp, *LxySig2_JMuMu_MuMu_Exp));

	//Yield
	RooRealVar N_JJ_PP("N_JJ_PP", "N_JJ_PP", 2000, 0, 1E6);

	RooRealVar N_JJ_NPP("N_JJ_NPP", "N_JJ_NPP", 500, 0, 1E4);
	RooRealVar N_JJ_NPNP("N_JJ_NPNP", "N_JJ_NPNP", 1000, 0, 1E6);

	RooRealVar N_JMuMu("N_JMuMu", "N_JMuMu", 0, 0, 5000);
	RooRealVar N_MuMuMuMu("N_MuMuMuMu", "N_MuMuMuMu", 0, 0, 300);

	//Fitting function
	RooArgList PDF_List(PDF_JJ_PP, PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP, PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu);
	RooArgList Yield_List(N_JJ_PP, N_JJ_NPP, N_JJ_NPP, N_JJ_NPNP, N_JMuMu, N_JMuMu, N_MuMuMuMu);

	RooAbsPdf* FittingFunction = new RooAddPdf("FittingFunction", "FittingFunction", PDF_List, Yield_List);

	//Mass windows cut, for insurance
	TString MassCut = "";
	MassCut = MassCut + "(J1Mass_ >" + LowCut + ") && (J1Mass_ < " + HighCut + ") && (J2Mass_ >" + LowCut + ") && (J2Mass_ < " + HighCut + ")";
	MassCut = MassCut + "&& (FourMuonMass_ > " + FourMuonMass_LowCut + ") && (FourMuonMass_ < " + FourMuonMass_HighCut + ")";

	//Data acquire
	TFile F_input("../../Sample/Data_MixWeighted", "read");
	TTree* Tree_ = (TTree*)gROOT->FindObject("WeightedTree");

	RooArgList Parameter_List(J1Mass_, J2Mass_, x1_, x2_, LxyPVSig1_, LxyPVSig2_, Weight_sum, FourMuonMass_);
	RooDataSet* Set = new RooDataSet("Set", "Set", Tree_, Parameter_List, MassCut, Weight_sum.GetName());

	//Fit
	RooFitResult* FitResult = FittingFunction->fitTo(*Set, Extended(kTRUE), Hesse(kTRUE), Strategy(1), NumCPU(4), Save(kTRUE), Minos(kFALSE));

	//Plot
	//Plot J1
	RooPlot* J1Frame = J1Mass_.frame(Title("M_{J/#psi1}"));

	Set->plotOn(J1Frame, MarkerColor(kBlack), LineColor(kBlack), Name("Dat"));
	FittingFunction->plotOn(J1Frame, LineColor(kBlue), Name("Total"));

	FittingFunction->plotOn(J1Frame, Components(PDF_JJ_PP), LineColor(kRed), Name("PP"));

	FittingFunction->plotOn(J1Frame, Components(PDF_JJ_NPP), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NPP"));
	FittingFunction->plotOn(J1Frame, Components(PDF_JJ_PNP), LineColor(kMagenta - 1), LineStyle(kDashed), Name("PNP"));
	FittingFunction->plotOn(J1Frame, Components(PDF_JJ_NPNP), LineColor(kMagenta - 5), LineStyle(kDashed), Name("NPNP"));

	FittingFunction->plotOn(J1Frame, Components(PDF_J1MuMu), LineColor(kCyan + 3), LineStyle(kDashed), Name("J1MuMu"));
	FittingFunction->plotOn(J1Frame, Components(PDF_MuMuJ2), LineColor(kCyan - 2), LineStyle(kDashed), Name("MuMuJ2"));
	FittingFunction->plotOn(J1Frame, Components(PDF_MuMuMuMu), LineColor(kCyan - 6), LineStyle(kDashed), Name("MuMuMuMu"));
	/*
	FittingFunction->plotOn(J1Frame, Components(RooArgSet(PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP)), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NP"));
	FittingFunction->plotOn(J1Frame, Components(RooArgSet(PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu)), LineColor(kOrange), LineStyle(kDashed), Name("Combi"));
	*/
	//Pull
	RooHist* J1_Pull = J1Frame->pullHist("Dat", "Total");
	RooPlot* J1_Pull_Frame = J1Mass_.frame(Title(" "));
	J1_Pull_Frame->addPlotable(J1_Pull, "P");

	//Legend
	TLegend legend(0.7, 0.5, 0.9, 0.9);
	//TLegend legend(0.7, 0.65, 0.9, 0.9);

	legend.AddEntry(J1Frame->findObject("Dat"), "Data(2016)", "LP");
	legend.AddEntry(J1Frame->findObject("Total"), "Total PDF", "L");

	legend.AddEntry(J1Frame->findObject("PP"), "J/#psi_{1}(P)J/#psi_{2}(P)", "L");

	legend.AddEntry(J1Frame->findObject("NPP"), "J/#psi_{1}(NP)J/#psi_{2}(P)", "L");
	legend.AddEntry(J1Frame->findObject("PNP"), "J/#psi_{1}(P)J/#psi_{2}(NP)", "L");
	legend.AddEntry(J1Frame->findObject("NPNP"), "J/#psi_{1}(NP)J/#psi_{2}(NP)", "L");

	legend.AddEntry(J1Frame->findObject("J1MuMu"), "J/#psi_{1}#mu^{+}#mu^{-}", "L");
	legend.AddEntry(J1Frame->findObject("MuMuJ2"), "#mu^{+}#mu^{-}J/#psi_{2}", "L");
	legend.AddEntry(J1Frame->findObject("MuMuMuMu"), "#mu^{+}#mu^{-}#mu^{+}#mu^{-}", "L");
	/*
	legend.AddEntry(J1Frame->findObject("PP"), "J/#psi_{1}J/#psi_{2}(P)", "L");
	legend.AddEntry(J1Frame->findObject("NP"), "J/#psi_{1}J/#psi_{2}(NP)", "L");
	legend.AddEntry(J1Frame->findObject("Combi"), "Combi", "L");
	*/
	//Draw
	TString OutputNumber = "";
	OutputNumber = OutputNumber + MassNumber;

	TCanvas* J1_Canvas = new TCanvas("J1_Canvas", "J1_Canvas", 3000, 3000);
	J1_Canvas->Divide(1, 2);
	J1_Canvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	J1Frame->Draw();
	legend.DrawClone();

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.5 * J1_Canvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.18, 0.85, "CMS");

	J1_Canvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	J1_Pull_Frame->GetYaxis()->SetTitle("Pull");
	J1_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	J1_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	J1_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	J1_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	J1_Pull_Frame->Draw();

	J1_Canvas->SaveAs("./Plot/J1_Result" + OutputNumber + ".pdf");
	delete J1_Canvas;

	//Plot J2
	RooPlot* J2Frame = J2Mass_.frame(Title("M_{J/#psi2}"));

	Set->plotOn(J2Frame, MarkerColor(kBlack), LineColor(kBlack), Name("Dat"));
	FittingFunction->plotOn(J2Frame, LineColor(kBlue), Name("Total"));

	FittingFunction->plotOn(J2Frame, Components(PDF_JJ_PP), LineColor(kRed), Name("PP"));

	FittingFunction->plotOn(J2Frame, Components(PDF_JJ_NPP), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NPP"));
	FittingFunction->plotOn(J2Frame, Components(PDF_JJ_PNP), LineColor(kMagenta - 1), LineStyle(kDashed), Name("PNP"));
	FittingFunction->plotOn(J2Frame, Components(PDF_JJ_NPNP), LineColor(kMagenta - 5), LineStyle(kDashed), Name("NPNP"));

	FittingFunction->plotOn(J2Frame, Components(PDF_J1MuMu), LineColor(kCyan + 3), LineStyle(kDashed), Name("J1MuMu"));
	FittingFunction->plotOn(J2Frame, Components(PDF_MuMuJ2), LineColor(kCyan - 2), LineStyle(kDashed), Name("MuMuJ2"));
	FittingFunction->plotOn(J2Frame, Components(PDF_MuMuMuMu), LineColor(kCyan - 6), LineStyle(kDashed), Name("MuMuMuMu"));
	/*
	FittingFunction->plotOn(J2Frame, Components(RooArgSet(PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP)), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NP"));
	FittingFunction->plotOn(J2Frame, Components(RooArgSet(PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu)), LineColor(kOrange), LineStyle(kDashed), Name("Combi"));
	*/
	//Pull
	RooHist* J2_Pull = J2Frame->pullHist("Dat", "Total");
	RooPlot* J2_Pull_Frame = J2Mass_.frame(Title(" "));
	J2_Pull_Frame->addPlotable(J2_Pull, "P");

	//Draw
	TCanvas* J2_Canvas = new TCanvas("J2_Canvas", "J2_Canvas", 3000, 3000);
	J2_Canvas->Divide(1, 2);
	J2_Canvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	J2Frame->Draw();
	legend.DrawClone();

	latex.SetTextSize(0.5 * J2_Canvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.18, 0.85, "CMS");

	J2_Canvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	J2_Pull_Frame->GetYaxis()->SetTitle("Pull");
	J2_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	J2_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	J2_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	J2_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	J2_Pull_Frame->Draw();

	J2_Canvas->SaveAs("./Plot/J2_Result" + OutputNumber + ".pdf");
	delete J2_Canvas;

	//Plot LxySig1
	RooPlot* LxySig1Frame = LxyPVSig1_.frame(Title("Significance of L_{xy}^{J/#psi1}PV"), Range(+0, 40));

	Set->plotOn(LxySig1Frame, MarkerColor(kBlack), LineColor(kBlack), Name("Dat"));
	FittingFunction->plotOn(LxySig1Frame, LineColor(kBlue), Name("Total"));

	FittingFunction->plotOn(LxySig1Frame, Components(PDF_JJ_PP), LineColor(kRed), Name("PP"));

	FittingFunction->plotOn(LxySig1Frame, Components(PDF_JJ_NPP), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NPP"));
	FittingFunction->plotOn(LxySig1Frame, Components(PDF_JJ_PNP), LineColor(kMagenta - 1), LineStyle(kDashed), Name("PNP"));
	FittingFunction->plotOn(LxySig1Frame, Components(PDF_JJ_NPNP), LineColor(kMagenta - 5), LineStyle(kDashed), Name("NPNP"));

	FittingFunction->plotOn(LxySig1Frame, Components(PDF_J1MuMu), LineColor(kCyan + 3), LineStyle(kDashed), Name("J1MuMu"));
	FittingFunction->plotOn(LxySig1Frame, Components(PDF_MuMuJ2), LineColor(kCyan - 2), LineStyle(kDashed), Name("MuMuJ2"));
	FittingFunction->plotOn(LxySig1Frame, Components(PDF_MuMuMuMu), LineColor(kCyan - 6), LineStyle(kDashed), Name("MuMuMuMu"));
	/*
	FittingFunction->plotOn(LxySig1Frame, Components(RooArgSet(PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP)), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NP"));
	FittingFunction->plotOn(LxySig1Frame, Components(RooArgSet(PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu)), LineColor(kOrange), LineStyle(kDashed), Name("Combi"));
	*/
	//Pull
	RooHist* LxySig1_Pull = LxySig1Frame->pullHist("Dat", "Total");
	LxySig1_Pull->SetMarkerStyle(20);
	LxySig1_Pull->SetMarkerSize(3);
	RooPlot* LxySig1_Pull_Frame = LxyPVSig1_.frame(Title(" "), Range(+0, 40));
	LxySig1_Pull_Frame->addPlotable(LxySig1_Pull, "P");

	//Draw
	TCanvas* LxySig1_Canvas = new TCanvas("LxySig1_Canvas", "LxySig1_Canvas", 3000, 3000);
	LxySig1_Canvas->Divide(1, 2);
	LxySig1_Canvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy();
	LxySig1Frame->Draw();
	legend.DrawClone();

	latex.SetTextSize(0.5 * LxySig1_Canvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.60, 0.85, "CMS");

	LxySig1_Canvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	LxySig1_Pull_Frame->GetYaxis()->SetTitle("Pull");
	LxySig1_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	LxySig1_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	LxySig1_Pull_Frame->GetYaxis()->SetTitleOffset(0.5);
	LxySig1_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	LxySig1_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	LxySig1_Pull_Frame->Draw();

	LxySig1_Canvas->SaveAs("./Plot/SigLxy1_Result" + OutputNumber + ".pdf");
	delete LxySig1_Canvas;

	//Plot LxySig2
	RooPlot* LxySig2Frame = LxyPVSig2_.frame(Title("Significance of L_{xy}^{J/#psi2}PV"), Range(+0, 40));

	Set->plotOn(LxySig2Frame, MarkerColor(kBlack), LineColor(kBlack), Name("Dat"));
	FittingFunction->plotOn(LxySig2Frame, LineColor(kBlue), Name("Total"));

	FittingFunction->plotOn(LxySig2Frame, Components(PDF_JJ_PP), LineColor(kRed), Name("PP"));

	FittingFunction->plotOn(LxySig2Frame, Components(PDF_JJ_NPP), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NPP"));
	FittingFunction->plotOn(LxySig2Frame, Components(PDF_JJ_PNP), LineColor(kMagenta - 1), LineStyle(kDashed), Name("PNP"));
	FittingFunction->plotOn(LxySig2Frame, Components(PDF_JJ_NPNP), LineColor(kMagenta - 5), LineStyle(kDashed), Name("NPNP"));

	FittingFunction->plotOn(LxySig2Frame, Components(PDF_J1MuMu), LineColor(kCyan + 3), LineStyle(kDashed), Name("J1MuMu"));
	FittingFunction->plotOn(LxySig2Frame, Components(PDF_MuMuJ2), LineColor(kCyan - 2), LineStyle(kDashed), Name("MuMuJ2"));
	FittingFunction->plotOn(LxySig2Frame, Components(PDF_MuMuMuMu), LineColor(kCyan - 6), LineStyle(kDashed), Name("MuMuMuMu"));
	/*
	FittingFunction->plotOn(LxySig2Frame, Components(RooArgSet(PDF_JJ_NPP, PDF_JJ_PNP, PDF_JJ_NPNP)), LineColor(kMagenta + 4), LineStyle(kDashed), Name("NP"));
	FittingFunction->plotOn(LxySig2Frame, Components(RooArgSet(PDF_J1MuMu, PDF_MuMuJ2, PDF_MuMuMuMu)), LineColor(kOrange), LineStyle(kDashed), Name("Combi"));
	*/
	//Pull
	RooHist* LxySig2_Pull = LxySig2Frame->pullHist("Dat", "Total");
	LxySig1_Pull->SetMarkerStyle(20);
	LxySig1_Pull->SetMarkerSize(3);
	RooPlot* LxySig2_Pull_Frame = LxyPVSig2_.frame(Title(" "), Range(+0, 40));
	LxySig2_Pull_Frame->addPlotable(LxySig2_Pull, "P");

	//Draw
	TCanvas* LxySig2_Canvas = new TCanvas("LxySig2_Canvas", "LxySig2_Canvas", 3000, 3000);
	LxySig2_Canvas->Divide(1, 2);
	LxySig2_Canvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy();
	LxySig2Frame->Draw();
	legend.DrawClone();

	latex.SetTextSize(0.5 * LxySig2_Canvas->GetTopMargin());
	latex.SetTextFont(62);
	latex.SetTextAlign(11);
	latex.DrawLatex(0.60, 0.85, "CMS");

	LxySig2_Canvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	LxySig2_Pull_Frame->GetYaxis()->SetTitle("Pull");
	LxySig2_Pull_Frame->GetYaxis()->SetTitleSize(0.1);
	LxySig2_Pull_Frame->GetYaxis()->SetLabelSize(0.1);
	LxySig2_Pull_Frame->GetYaxis()->SetTitleOffset(0.5);
	LxySig2_Pull_Frame->GetXaxis()->SetTitleSize(0.1);
	LxySig2_Pull_Frame->GetXaxis()->SetLabelSize(0.1);
	LxySig2_Pull_Frame->Draw();

	LxySig2_Canvas->SaveAs("./Plot/SigLxy2_Result" + OutputNumber + ".pdf");
	delete LxySig2_Canvas;

	//Print the result to the target file
	FILE* OutputFile = NULL;
	OutputFile = fopen(Output_Name + ".txt", "a");

	fprintf(OutputFile, "%i: %i\n", MassNumber, FitResult->status());
	fprintf(OutputFile, "PP: %f +/- %f\n", N_JJ_PP.getValV(), N_JJ_PP.getError());
	fprintf(OutputFile, "PNP: %f +/- %f\n", N_JJ_NPP.getValV(), N_JJ_NPP.getError());
	fprintf(OutputFile, "NPNP: %f +/- %f\n", N_JJ_NPNP.getValV(), N_JJ_NPNP.getError());
	fprintf(OutputFile, "JMuMu: %f +/- %f\n", N_JMuMu.getValV(), N_JMuMu.getError());
	fprintf(OutputFile, "MuMuMuMu: %f +/- %f\n", N_MuMuMuMu.getValV(), N_MuMuMuMu.getError());
	fprintf(OutputFile, "Sys: %f\n", (N_JJ_PP.getValV() - Standard_Value) / Standard_Value);

	fclose(OutputFile);
}
