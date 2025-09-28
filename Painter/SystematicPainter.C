//Jinfeng, 2025.09.25
//To draw the systematic uncertainties

#include "TGraphErrors.h"

void SystematicPainter()
{

	//Bin edge (add the "overflow edge" here)
	float JJMassBinEdge[] = { 7.5,17.5,27.5,37.5,47.5,57.5,67.5,77.5,87.5 };
	float DeltaYBinEdge[] = { 0.0,0.5,1.0,1.5,2.0,2.5,3.0 };
	float JJPtBinEdge[] = { 0,5,10,15,20,25,30,35,40,45 };
	float JJYBinEdge[] = { 0,0.4,0.8,1.2,1.6,2.0 };
	float DeltaPhiBinEdge[] = { 0.0,TMath::Pi() / 8,TMath::Pi() / 4,TMath::Pi() * 3 / 8,TMath::Pi() / 2,TMath::Pi() * 5 / 8,TMath::Pi() * 3 / 4,TMath::Pi() * 7 / 8,TMath::Pi() };
	float JPtBinEdge[] = { 10,15,20,25,30,35,40 };

	//Content 
	float JJMass_Lumi[] = { 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12 }; //Luminosity
	float JJMass_BR[] = { 1.106, 1.106, 1.106, 1.106, 1.106, 1.106, 1.106, 1.106 }; //Branch ratio
	float JJMass_Correction[] = { 14.67, 17.22, 16.08, 13.15, 14.88, 13.72, 18.95, 25.23 }; //Correction
	float JJMass_Ctau[] = { 2.727, 4.331, 6.181, 4.092, 5.872, 4.032, 0.899, 6.778 }; //Ctau shape
	float JJMass_SigLxy[] = { 0.253, 0.002, 3.811, 4.085, 12.166, 3.995, 16.180, 8.800 }; //Lifetime variable (SigLxy)
	float JJMass_Fitter[] = { 0.764, 0.455, 5.997, 1.076, 2.131, 4.204, 6.754, 3.209 }; //Fitter stability
	float JJMass_Fix[] = { 0.81, 0.15, 0.89, 0.54, 6.06, 1.61, 7.88, 4.83 }; //Parameter fix
	float JJMass_Total[] = { 15.007, 17.800, 18.687, 14.461, 21.126, 15.553, 28.467, 28.190}; //Total uncertainty

	float DeltaY_Lumi[] = { 0.12, 0.12, 0.12, 0.12, 0.12, 0.12 };
	float DeltaY_BR[] = { 1.106, 1.106, 1.106, 1.106, 1.106, 1.106 };
	float DeltaY_Correction[] = { 14.234, 20.327, 24.430, 21.076, 25.562, 36.862 };
	float DeltaY_Ctau[] = { 3.218, 5.078, 4.805, 5.168, 5.018, 4.116 };
	float DeltaY_SigLxy[] = { 0.659, 1.331, 2.856, 1.920, 7.014, 7.930 };
	float DeltaY_Fitter[] = { 0.571, 1.255, 2.490, 1.305, 2.568, 3.662 };
	float DeltaY_Fix[] = { 2.86, 7.26, 4.84, 5.51, 5.62, 3.71 };
	float DeltaY_Total[] = { 14.937, 22.277, 25.670, 22.536, 27.699, 38.302 };

	float JJPt_Lumi[] = { 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12};
	float JJPt_BR[] = { 1.106, 1.106, 1.106, 1.106, 1.106, 1.106, 1.106, 1.106, 1.106 };
	float JJPt_Correction[] = { 16.749, 16.217, 15.897, 14.735, 15.339, 16.012, 10.812, 15.768, 32.783 };
	float JJPt_Ctau[] = { 3.896, 3.624, 5.453, 7.301, 3.809, 4.752, 3.843, 2.718, 3.412 };
	float JJPt_SigLxy[] = { 2.983, 0.168, 3.125, 4.491, 0.397, 2.163, 3.288, 14.683, 4.536 };
	float JJPt_Fitter[] = { 1.002, 1.897, 1.006, 3.952, 1.086, 1.487, 1.420, 1.573, 0.720 };
	float JJPt_Fix[] = { 4.44, 21.3, 7.44, 12.2, 3.40, 5.67, 4.20, 2.45, 3.79 };
	float JJPt_Total[] = { 18.062, 27.105, 18.705, 21.361, 16.246, 17.868, 12.782, 21.939, 33.512 };

	float JJY_Lumi[] = { 0.12, 0.12, 0.12, 0.12, 0.12 };
	float JJY_BR[] = { 1.106, 1.106, 1.106, 1.106, 1.106 };
	float JJY_Correction[] = { 15.671, 15.609, 16.080, 12.235, 12.184 };
	float JJY_Ctau[] = { 3.562, 4.031, 4.235, 4.172, 6.253 };
	float JJY_SigLxy[] = { 0.770, 1.594, 0.799, 3.579, 8.584 };
	float JJY_Fitter[] = { 0.707, 1.074, 0.647, 0.458, 6.057 };
	float JJY_Fix[] = { 3.77, 3.91, 4.48, 6.02, 5.79 };
	float JJY_Total[] = { 16.578, 16.738, 17.288, 14.752, 18.240 };

	float DeltaPhi_Lumi[] = { 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12 };
	float DeltaPhi_BR[] = { 1.106, 1.106, 1.106, 1.106, 1.106, 1.106, 1.106, 1.106 };
	float DeltaPhi_Correction[] = { 18.547, 13.055, 21.782, 31.598, 19.444, 22.171, 24.587, 23.997 };
	float DeltaPhi_Ctau[] = { 3.047, 2.839, 3.832, 4.959, 5.074, 3.635, 4.058, 3.928 };
	float DeltaPhi_SigLxy[] = { 0.138, 3.109, 0.027, 5.703, 4.946, 1.548, 0.525, 0.144};
	float DeltaPhi_Fitter[] = { 0.992, 1.607, 1.737, 2.536, 5.864, 2.632, 2.403, 1.703 };
	float DeltaPhi_Fix[] = { 0.18, 0.97, 1.20, 0.503, 2.25, 1.93, 1.53, 0.263 };
	float DeltaPhi_Total[] = { 18.856, 13.889, 22.245, 32.611, 21.656, 22.841, 25.113, 24.403 };

	float JPt_Lumi[] = { 0.12, 0.12, 0.12, 0.12, 0.12, 0.12 };
	float JPt_BR[] = { 1.106, 1.106, 1.106, 1.106, 1.106, 1.106 };
	float JPt_Correction[] = { 16.657, 13.112, 13.076, 13, 13, 13 };
	float JPt_Ctau[] = { 4.127, 3.338, 3.646, 5.356, 16.242, 11.387 };
	float JPt_SigLxy[] = { 2.884, 1.7713, 8.510, 17.886, 16.945, 21.264 };
	float JPt_Fitter[] = { 0.885, 1.113, 3.409, 1.734, 13.1, 25.677 };
	float JPt_Fix[] = { 4.56, 5.58, 4.05, 5.45, 9.07, 6.96 };
	float JPt_Total[] = { 18.045, 14.826, 16.910, 23.485, 31.321, 38.207 };

	//Graph
	float JJMass_Center[8];
	float JJMass_XError[8];
	float JJMass_YError[8];

	for (int i = 0; i < 8; i++)
	{
		JJMass_Center[i] = (JJMassBinEdge[i + 1] + JJMassBinEdge[i]) / 2;
		JJMass_XError[i] = (JJMassBinEdge[i + 1] - JJMassBinEdge[i]) / 2;
		JJMass_YError[i] = 0; //No vertical error bar
	}

	auto JJMass_LumiGraph = new TGraphErrors(8, JJMass_Center, JJMass_Lumi, JJMass_XError, JJMass_YError);
	auto JJMass_BRGraph = new TGraphErrors(8, JJMass_Center, JJMass_BR, JJMass_XError, JJMass_YError);
	auto JJMass_CorrectionGraph = new TGraphErrors(8, JJMass_Center, JJMass_Correction, JJMass_XError, JJMass_YError);
	auto JJMass_CtauGraph = new TGraphErrors(8, JJMass_Center, JJMass_Ctau, JJMass_XError, JJMass_YError);
	auto JJMass_SigLxyGraph = new TGraphErrors(8, JJMass_Center, JJMass_SigLxy, JJMass_XError, JJMass_YError);
	auto JJMass_FitterGraph = new TGraphErrors(8, JJMass_Center, JJMass_Fitter, JJMass_XError, JJMass_YError);
	auto JJMass_FixGraph = new TGraphErrors(8, JJMass_Center, JJMass_Fix, JJMass_XError, JJMass_YError);
	auto JJMass_TotalGraph = new TGraphErrors(8, JJMass_Center, JJMass_Total, JJMass_XError, JJMass_YError);

	//
	float DeltaY_Center[6];
	float DeltaY_XError[6];
	float DeltaY_YError[6];

	for (int i = 0; i < 6; i++)
	{
		DeltaY_Center[i] = (DeltaYBinEdge[i + 1] + DeltaYBinEdge[i]) / 2;
		DeltaY_XError[i] = (DeltaYBinEdge[i + 1] - DeltaYBinEdge[i]) / 2;
		DeltaY_YError[i] = 0;
	}

	auto DeltaY_LumiGraph = new TGraphErrors(6, DeltaY_Center, DeltaY_Lumi, DeltaY_XError, DeltaY_YError);
	auto DeltaY_BRGraph = new TGraphErrors(6, DeltaY_Center, DeltaY_BR, DeltaY_XError, DeltaY_YError);
	auto DeltaY_CorrectionGraph = new TGraphErrors(6, DeltaY_Center, DeltaY_Correction, DeltaY_XError, DeltaY_YError);
	auto DeltaY_CtauGraph = new TGraphErrors(6, DeltaY_Center, DeltaY_Ctau, DeltaY_XError, DeltaY_YError);
	auto DeltaY_SigLxyGraph = new TGraphErrors(6, DeltaY_Center, DeltaY_SigLxy, DeltaY_XError, DeltaY_YError);
	auto DeltaY_FitterGraph = new TGraphErrors(6, DeltaY_Center, DeltaY_Fitter, DeltaY_XError, DeltaY_YError);
	auto DeltaY_FixGraph = new TGraphErrors(6, DeltaY_Center, DeltaY_Fix, DeltaY_XError, DeltaY_YError);
	auto DeltaY_TotalGraph = new TGraphErrors(6, DeltaY_Center, DeltaY_Total, DeltaY_XError, DeltaY_YError);

	//
	float JJPt_Center[9];
	float JJPt_XError[9];
	float JJPt_YError[9];

	for (int i = 0; i < 9; i++)
	{
		JJPt_Center[i] = (JJPtBinEdge[i + 1] + JJPtBinEdge[i]) / 2;
		JJPt_XError[i] = (JJPtBinEdge[i + 1] - JJPtBinEdge[i]) / 2;
		JJPt_YError[i] = 0;
	}

	auto JJPt_LumiGraph = new TGraphErrors(9, JJPt_Center, JJPt_Lumi, JJPt_XError, JJPt_YError);
	auto JJPt_BRGraph = new TGraphErrors(9, JJPt_Center, JJPt_BR, JJPt_XError, JJPt_YError);
	auto JJPt_CorrectionGraph = new TGraphErrors(9, JJPt_Center, JJPt_Correction, JJPt_XError, JJPt_YError);
	auto JJPt_CtauGraph = new TGraphErrors(9, JJPt_Center, JJPt_Ctau, JJPt_XError, JJPt_YError);
	auto JJPt_SigLxyGraph = new TGraphErrors(9, JJPt_Center, JJPt_SigLxy, JJPt_XError, JJPt_YError);
	auto JJPt_FitterGraph = new TGraphErrors(9, JJPt_Center, JJPt_Fitter, JJPt_XError, JJPt_YError);
	auto JJPt_FixGraph = new TGraphErrors(9, JJPt_Center, JJPt_Fix, JJPt_XError, JJPt_YError);
	auto JJPt_TotalGraph = new TGraphErrors(9, JJPt_Center, JJPt_Total, JJPt_XError, JJPt_YError);

	//
	float JJY_Center[5];
	float JJY_XError[5];
	float JJY_YError[5];

	for (int i = 0; i < 5; i++)
	{
		JJY_Center[i] = (JJYBinEdge[i + 1] + JJYBinEdge[i]) / 2;
		JJY_XError[i] = (JJYBinEdge[i + 1] - JJYBinEdge[i]) / 2;
		JJY_YError[i] = 0;
	}

	auto JJY_LumiGraph = new TGraphErrors(5, JJY_Center, JJY_Lumi, JJY_XError, JJY_YError);
	auto JJY_BRGraph = new TGraphErrors(5, JJY_Center, JJY_BR, JJY_XError, JJY_YError);
	auto JJY_CorrectionGraph = new TGraphErrors(5, JJY_Center, JJY_Correction, JJY_XError, JJY_YError);
	auto JJY_CtauGraph = new TGraphErrors(5, JJY_Center, JJY_Ctau, JJY_XError, JJY_YError);
	auto JJY_SigLxyGraph = new TGraphErrors(5, JJY_Center, JJY_SigLxy, JJY_XError, JJY_YError);
	auto JJY_FitterGraph = new TGraphErrors(5, JJY_Center, JJY_Fitter, JJY_XError, JJY_YError);
	auto JJY_FixGraph = new TGraphErrors(5, JJY_Center, JJY_Fix, JJY_XError, JJY_YError);
	auto JJY_TotalGraph = new TGraphErrors(5, JJY_Center, JJY_Total, JJY_XError, JJY_YError);

	//
	float DeltaPhi_Center[8];
	float DeltaPhi_XError[8];
	float DeltaPhi_YError[8];

	for (int i = 0; i < 8; i++)
	{
		DeltaPhi_Center[i] = (DeltaPhiBinEdge[i + 1] + DeltaPhiBinEdge[i]) / 2;
		DeltaPhi_XError[i] = (DeltaPhiBinEdge[i + 1] - DeltaPhiBinEdge[i]) / 2;
		DeltaPhi_YError[i] = 0;
	}

	auto DeltaPhi_LumiGraph = new TGraphErrors(8, DeltaPhi_Center, DeltaPhi_Lumi, DeltaPhi_XError, DeltaPhi_YError);
	auto DeltaPhi_BRGraph = new TGraphErrors(8, DeltaPhi_Center, DeltaPhi_BR, DeltaPhi_XError, DeltaPhi_YError);
	auto DeltaPhi_CorrectionGraph = new TGraphErrors(8, DeltaPhi_Center, DeltaPhi_Correction, DeltaPhi_XError, DeltaPhi_YError);
	auto DeltaPhi_CtauGraph = new TGraphErrors(8, DeltaPhi_Center, DeltaPhi_Ctau, DeltaPhi_XError, DeltaPhi_YError);
	auto DeltaPhi_SigLxyGraph = new TGraphErrors(8, DeltaPhi_Center, DeltaPhi_SigLxy, DeltaPhi_XError, DeltaPhi_YError);
	auto DeltaPhi_FitterGraph = new TGraphErrors(8, DeltaPhi_Center, DeltaPhi_Fitter, DeltaPhi_XError, DeltaPhi_YError);
	auto DeltaPhi_FixGraph = new TGraphErrors(8, DeltaPhi_Center, DeltaPhi_Fix, DeltaPhi_XError, DeltaPhi_YError);
	auto DeltaPhi_TotalGraph = new TGraphErrors(8, DeltaPhi_Center, DeltaPhi_Total, DeltaPhi_XError, DeltaPhi_YError);

	//
	float JPt_Center[6];
	float JPt_XError[6];
	float JPt_YError[6];

	for (int i = 0; i < 6; i++)
	{
		JPt_Center[i] = (JPtBinEdge[i + 1] + JPtBinEdge[i]) / 2;
		JPt_XError[i] = (JPtBinEdge[i + 1] - JPtBinEdge[i]) / 2;
		JPt_YError[i] = 0;
	}

	auto JPt_LumiGraph = new TGraphErrors(6, JPt_Center, JPt_Lumi, JPt_XError, JPt_YError);
	auto JPt_BRGraph = new TGraphErrors(6, JPt_Center, JPt_BR, JPt_XError, JPt_YError);
	auto JPt_CorrectionGraph = new TGraphErrors(6, JPt_Center, JPt_Correction, JPt_XError, JPt_YError);
	auto JPt_CtauGraph = new TGraphErrors(6, JPt_Center, JPt_Ctau, JPt_XError, JPt_YError);
	auto JPt_SigLxyGraph = new TGraphErrors(6, JPt_Center, JPt_SigLxy, JPt_XError, JPt_YError);
	auto JPt_FitterGraph = new TGraphErrors(6, JPt_Center, JPt_Fitter, JPt_XError, JPt_YError);
	auto JPt_FixGraph = new TGraphErrors(6, JPt_Center, JPt_Fix, JPt_XError, JPt_YError);
	auto JPt_TotalGraph = new TGraphErrors(6, JPt_Center, JPt_Total, JPt_XError, JPt_YError);

	//Draw
	TCanvas* JJMassCanvas = new TCanvas("JJMassCanvas", "JJMassCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	JJMass_TotalGraph->SetTitle("");
	JJMass_TotalGraph->GetYaxis()->SetTitle("Relative systematics (%)");
	JJMass_TotalGraph->GetXaxis()->SetTitle("M_{J/#psi_{1}J/#psi_{2}} (GeV)");
	JJMass_TotalGraph->GetYaxis()->SetRangeUser(0, 45);
	JJMass_TotalGraph->GetXaxis()->SetRangeUser(7.5, 87.5);
	JJMass_TotalGraph->SetLineWidth(6);
	JJMass_TotalGraph->SetLineColor(kRed);
	JJMass_TotalGraph->SetMarkerStyle(20);
	JJMass_TotalGraph->SetMarkerSize(5);
	JJMass_TotalGraph->SetMarkerColor(kRed);
	JJMass_TotalGraph->Draw("AP same");
	JJMass_LumiGraph->SetLineWidth(2);
	JJMass_LumiGraph->SetLineColor(kBlack);
	JJMass_LumiGraph->SetMarkerStyle(20);
	JJMass_LumiGraph->SetMarkerSize(5);
	JJMass_LumiGraph->SetMarkerColor(kBlack);
	JJMass_LumiGraph->Draw("P same");
	JJMass_BRGraph->SetLineWidth(4);
	JJMass_BRGraph->SetLineColor(kBlue - 1);
	JJMass_BRGraph->SetMarkerStyle(20);
	JJMass_BRGraph->SetMarkerSize(5);
	JJMass_BRGraph->SetMarkerColor(kBlue - 1);
	JJMass_BRGraph->Draw("P same");
	JJMass_CorrectionGraph->SetLineWidth(4);
	JJMass_CorrectionGraph->SetLineColor(kBlue + 2);
	JJMass_CorrectionGraph->SetMarkerStyle(20);
	JJMass_CorrectionGraph->SetMarkerSize(5);
	JJMass_CorrectionGraph->SetMarkerColor(kBlue + 2);
	JJMass_CorrectionGraph->Draw("P same");
	JJMass_CtauGraph->SetLineWidth(4);
	JJMass_CtauGraph->SetLineColor(kBlue);
	JJMass_CtauGraph->SetMarkerStyle(20);
	JJMass_CtauGraph->SetMarkerSize(5);
	JJMass_CtauGraph->SetMarkerColor(kBlue);
	JJMass_CtauGraph->Draw("P same");
	JJMass_SigLxyGraph->SetLineWidth(4);
	JJMass_SigLxyGraph->SetLineColor(kBlue - 7);
	JJMass_SigLxyGraph->SetMarkerStyle(20);
	JJMass_SigLxyGraph->SetMarkerSize(5);
	JJMass_SigLxyGraph->SetMarkerColor(kBlue - 7);
	JJMass_SigLxyGraph->Draw("P same");
	JJMass_FitterGraph->SetLineWidth(4);
	JJMass_FitterGraph->SetLineColor(kBlue - 10);
	JJMass_FitterGraph->SetMarkerStyle(20);
	JJMass_FitterGraph->SetMarkerSize(5);
	JJMass_FitterGraph->SetMarkerColor(kBlue - 10);
	JJMass_FitterGraph->Draw("P same");
	JJMass_FixGraph->SetLineWidth(4);
	JJMass_FixGraph->SetLineColor(kGray);
	JJMass_FixGraph->SetMarkerStyle(20);
	JJMass_FixGraph->SetMarkerSize(5);
	JJMass_FixGraph->SetMarkerColor(kGray);
	JJMass_FixGraph->Draw("P same");

	TLegend Legend_JJMass(0.50, 0.70, 0.89, 0.88);
	Legend_JJMass.SetNColumns(2);
	Legend_JJMass.AddEntry(JJMass_LumiGraph, "Luminosity", "LP");
	Legend_JJMass.AddEntry(JJMass_BRGraph, "BR", "LP");
	Legend_JJMass.AddEntry(JJMass_CorrectionGraph, "Correction", "LP");
	Legend_JJMass.AddEntry(JJMass_CtauGraph, "c#tau shape", "LP");
	Legend_JJMass.AddEntry(JJMass_SigLxyGraph, "Lifetime variable", "LP");
	Legend_JJMass.AddEntry(JJMass_FitterGraph, "Fitter stability", "LP");
	Legend_JJMass.AddEntry(JJMass_FixGraph, "Paramters fix", "LP");
	Legend_JJMass.AddEntry(JJMass_TotalGraph, "Total", "LP");
	Legend_JJMass.SetLineWidth(0);
	Legend_JJMass.DrawClone();

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.05);//0.5
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    // default is helvetica-italics
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.22, 0.80, "(a)");

	JJMassCanvas->SaveAs("./SystematicPlot/JJMass.png");

	//
	TCanvas* DeltaYCanvas = new TCanvas("DeltaYCanvas", "DeltaYCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	DeltaY_TotalGraph->SetTitle("");
	DeltaY_TotalGraph->GetYaxis()->SetTitle("Relative systematics (%)");
	DeltaY_TotalGraph->GetXaxis()->SetTitle("#Deltay(J/#psi_{1},J/#psi_{2})");
	DeltaY_TotalGraph->GetYaxis()->SetRangeUser(0, 55);
	DeltaY_TotalGraph->GetXaxis()->SetRangeUser(0, 3.0);
	DeltaY_TotalGraph->SetLineWidth(6);
	DeltaY_TotalGraph->SetLineColor(kRed);
	DeltaY_TotalGraph->SetMarkerStyle(20);
	DeltaY_TotalGraph->SetMarkerSize(5);
	DeltaY_TotalGraph->SetMarkerColor(kRed);
	DeltaY_TotalGraph->Draw("AP same");
	DeltaY_LumiGraph->SetLineWidth(2);
	DeltaY_LumiGraph->SetLineColor(kBlack);
	DeltaY_LumiGraph->SetMarkerStyle(20);
	DeltaY_LumiGraph->SetMarkerSize(5);
	DeltaY_LumiGraph->SetMarkerColor(kBlack);
	DeltaY_LumiGraph->Draw("P same");
	DeltaY_BRGraph->SetLineWidth(4);
	DeltaY_BRGraph->SetLineColor(kBlue - 1);
	DeltaY_BRGraph->SetMarkerStyle(20);
	DeltaY_BRGraph->SetMarkerSize(5);
	DeltaY_BRGraph->SetMarkerColor(kBlue - 1);
	DeltaY_BRGraph->Draw("P same");
	DeltaY_CorrectionGraph->SetLineWidth(4);
	DeltaY_CorrectionGraph->SetLineColor(kBlue + 2);
	DeltaY_CorrectionGraph->SetMarkerStyle(20);
	DeltaY_CorrectionGraph->SetMarkerSize(5);
	DeltaY_CorrectionGraph->SetMarkerColor(kBlue + 2);
	DeltaY_CorrectionGraph->Draw("P same");
	DeltaY_CtauGraph->SetLineWidth(4);
	DeltaY_CtauGraph->SetLineColor(kBlue);
	DeltaY_CtauGraph->SetMarkerStyle(20);
	DeltaY_CtauGraph->SetMarkerSize(5);
	DeltaY_CtauGraph->SetMarkerColor(kBlue);
	DeltaY_CtauGraph->Draw("P same");
	DeltaY_SigLxyGraph->SetLineWidth(4);
	DeltaY_SigLxyGraph->SetLineColor(kBlue - 7);
	DeltaY_SigLxyGraph->SetMarkerStyle(20);
	DeltaY_SigLxyGraph->SetMarkerSize(5);
	DeltaY_SigLxyGraph->SetMarkerColor(kBlue - 7);
	DeltaY_SigLxyGraph->Draw("P same");
	DeltaY_FitterGraph->SetLineWidth(4);
	DeltaY_FitterGraph->SetLineColor(kBlue - 10);
	DeltaY_FitterGraph->SetMarkerStyle(20);
	DeltaY_FitterGraph->SetMarkerSize(5);
	DeltaY_FitterGraph->SetMarkerColor(kBlue - 10);
	DeltaY_FitterGraph->Draw("P same");
	DeltaY_FixGraph->SetLineWidth(4);
	DeltaY_FixGraph->SetLineColor(kGray);
	DeltaY_FixGraph->SetMarkerStyle(20);
	DeltaY_FixGraph->SetMarkerSize(5);
	DeltaY_FixGraph->SetMarkerColor(kGray);
	DeltaY_FixGraph->Draw("P same");

	TLegend Legend_DeltaY(0.50, 0.70, 0.89, 0.88);
	Legend_DeltaY.SetNColumns(2);
	Legend_DeltaY.AddEntry(DeltaY_LumiGraph, "Luminosity", "LP");
	Legend_DeltaY.AddEntry(DeltaY_BRGraph, "BR", "LP");
	Legend_DeltaY.AddEntry(DeltaY_CorrectionGraph, "Correction", "LP");
	Legend_DeltaY.AddEntry(DeltaY_CtauGraph, "c#tau shape", "LP");
	Legend_DeltaY.AddEntry(DeltaY_SigLxyGraph, "Lifetime variable", "LP");
	Legend_DeltaY.AddEntry(DeltaY_FitterGraph, "Fitter stability", "LP");
	Legend_DeltaY.AddEntry(DeltaY_FitterGraph, "Parameter fix", "LP");
	Legend_DeltaY.AddEntry(DeltaY_TotalGraph, "Total", "LP");
	Legend_DeltaY.SetLineWidth(0);
	Legend_DeltaY.DrawClone();

	latex.SetTextSize(0.05);//0.5
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    // default is helvetica-italics
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.22, 0.80, "(b)");

	DeltaYCanvas->SaveAs("./SystematicPlot/DeltaY.png");

	//
	TCanvas* JJPtCanvas = new TCanvas("JJPtCanvas", "JJPtCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	JJPt_TotalGraph->SetTitle("");
	JJPt_TotalGraph->GetYaxis()->SetTitle("Relative systematics (%)");
	JJPt_TotalGraph->GetXaxis()->SetTitle("p_{T}(J/#psi_{1}J/#psi_{2}) (GeV)");
	JJPt_TotalGraph->GetYaxis()->SetRangeUser(0, 50);
	JJPt_TotalGraph->GetXaxis()->SetRangeUser(0, 45);
	JJPt_TotalGraph->SetLineWidth(6);
	JJPt_TotalGraph->SetLineColor(kRed);
	JJPt_TotalGraph->SetMarkerStyle(20);
	JJPt_TotalGraph->SetMarkerSize(5);
	JJPt_TotalGraph->SetMarkerColor(kRed);
	JJPt_TotalGraph->Draw("AP same");
	JJPt_LumiGraph->SetLineWidth(2);
	JJPt_LumiGraph->SetLineColor(kBlack);
	JJPt_LumiGraph->SetMarkerStyle(20);
	JJPt_LumiGraph->SetMarkerSize(5);
	JJPt_LumiGraph->SetMarkerColor(kBlack);
	JJPt_LumiGraph->Draw("P same");
	JJPt_BRGraph->SetLineWidth(4);
	JJPt_BRGraph->SetLineColor(kBlue - 1);
	JJPt_BRGraph->SetMarkerStyle(20);
	JJPt_BRGraph->SetMarkerSize(5);
	JJPt_BRGraph->SetMarkerColor(kBlue - 1);
	JJPt_BRGraph->Draw("P same");
	JJPt_CorrectionGraph->SetLineWidth(4);
	JJPt_CorrectionGraph->SetLineColor(kBlue + 2);
	JJPt_CorrectionGraph->SetMarkerStyle(20);
	JJPt_CorrectionGraph->SetMarkerSize(5);
	JJPt_CorrectionGraph->SetMarkerColor(kBlue + 2);
	JJPt_CorrectionGraph->Draw("P same");
	JJPt_CtauGraph->SetLineWidth(4);
	JJPt_CtauGraph->SetLineColor(kBlue);
	JJPt_CtauGraph->SetMarkerStyle(20);
	JJPt_CtauGraph->SetMarkerSize(5);
	JJPt_CtauGraph->SetMarkerColor(kBlue);
	JJPt_CtauGraph->Draw("P same");
	JJPt_SigLxyGraph->SetLineWidth(4);
	JJPt_SigLxyGraph->SetLineColor(kBlue - 7);
	JJPt_SigLxyGraph->SetMarkerStyle(20);
	JJPt_SigLxyGraph->SetMarkerSize(5);
	JJPt_SigLxyGraph->SetMarkerColor(kBlue - 7);
	JJPt_SigLxyGraph->Draw("P same");
	JJPt_FitterGraph->SetLineWidth(4);
	JJPt_FitterGraph->SetLineColor(kBlue - 10);
	JJPt_FitterGraph->SetMarkerStyle(20);
	JJPt_FitterGraph->SetMarkerSize(5);
	JJPt_FitterGraph->SetMarkerColor(kBlue - 10);
	JJPt_FitterGraph->Draw("P same");
	JJPt_FixGraph->SetLineWidth(4);
	JJPt_FixGraph->SetLineColor(kGray);
	JJPt_FixGraph->SetMarkerStyle(20);
	JJPt_FixGraph->SetMarkerSize(5);
	JJPt_FixGraph->SetMarkerColor(kGray);
	JJPt_FixGraph->Draw("P same");

	TLegend Legend_JJPt(0.50, 0.70, 0.89, 0.88);
	Legend_JJPt.SetNColumns(2);
	Legend_JJPt.AddEntry(JJPt_LumiGraph, "Luminosity", "LP");
	Legend_JJPt.AddEntry(JJPt_BRGraph, "BR", "LP");
	Legend_JJPt.AddEntry(JJPt_CorrectionGraph, "Correction", "LP");
	Legend_JJPt.AddEntry(JJPt_CtauGraph, "c#tau shape", "LP");
	Legend_JJPt.AddEntry(JJPt_SigLxyGraph, "Lifetime variable", "LP");
	Legend_JJPt.AddEntry(JJPt_FitterGraph, "Fitter stability", "LP");
	Legend_JJPt.AddEntry(JJPt_FixGraph, "Parameter fix", "LP");
	Legend_JJPt.AddEntry(JJPt_TotalGraph, "Total", "LP");
	Legend_JJPt.SetLineWidth(0);
	Legend_JJPt.DrawClone();

	latex.SetTextSize(0.05);//0.5
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    // default is helvetica-italics
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.22, 0.80, "(c)");

	JJPtCanvas->SaveAs("./SystematicPlot/JJPt.png");

	//
	TCanvas* JJYCanvas = new TCanvas("JJYCanvas", "JJYCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	JJY_TotalGraph->SetTitle("");
	JJY_TotalGraph->GetYaxis()->SetTitle("Relative systematics (%)");
	JJY_TotalGraph->GetXaxis()->SetTitle("|y(J/#psi_{1}J/#psi_{2})|");
	JJY_TotalGraph->GetYaxis()->SetRangeUser(0, 35);
	JJY_TotalGraph->GetXaxis()->SetRangeUser(0, 2.0);
	JJY_TotalGraph->SetLineWidth(6);
	JJY_TotalGraph->SetLineColor(kRed);
	JJY_TotalGraph->SetMarkerStyle(20);
	JJY_TotalGraph->SetMarkerSize(5);
	JJY_TotalGraph->SetMarkerColor(kRed);
	JJY_TotalGraph->Draw("AP same");
	JJY_LumiGraph->SetLineWidth(2);
	JJY_LumiGraph->SetLineColor(kBlack);
	JJY_LumiGraph->SetMarkerStyle(20);
	JJY_LumiGraph->SetMarkerSize(5);
	JJY_LumiGraph->SetMarkerColor(kBlack);
	JJY_LumiGraph->Draw("P same");
	JJY_BRGraph->SetLineWidth(4);
	JJY_BRGraph->SetLineColor(kBlue - 1);
	JJY_BRGraph->SetMarkerStyle(20);
	JJY_BRGraph->SetMarkerSize(5);
	JJY_BRGraph->SetMarkerColor(kBlue - 1);
	JJY_BRGraph->Draw("P same");
	JJY_CorrectionGraph->SetLineWidth(4);
	JJY_CorrectionGraph->SetLineColor(kBlue + 2);
	JJY_CorrectionGraph->SetMarkerStyle(20);
	JJY_CorrectionGraph->SetMarkerSize(5);
	JJY_CorrectionGraph->SetMarkerColor(kBlue + 2);
	JJY_CorrectionGraph->Draw("P same");
	JJY_CtauGraph->SetLineWidth(4);
	JJY_CtauGraph->SetLineColor(kBlue);
	JJY_CtauGraph->SetMarkerStyle(20);
	JJY_CtauGraph->SetMarkerSize(5);
	JJY_CtauGraph->SetMarkerColor(kBlue);
	JJY_CtauGraph->Draw("P same");
	JJY_SigLxyGraph->SetLineWidth(4);
	JJY_SigLxyGraph->SetLineColor(kBlue - 7);
	JJY_SigLxyGraph->SetMarkerStyle(20);
	JJY_SigLxyGraph->SetMarkerSize(5);
	JJY_SigLxyGraph->SetMarkerColor(kBlue - 7);
	JJY_SigLxyGraph->Draw("P same");
	JJY_FitterGraph->SetLineWidth(4);
	JJY_FitterGraph->SetLineColor(kBlue - 10);
	JJY_FitterGraph->SetMarkerStyle(20);
	JJY_FitterGraph->SetMarkerSize(5);
	JJY_FitterGraph->SetMarkerColor(kBlue - 10);
	JJY_FitterGraph->Draw("P same");
	JJY_FixGraph->SetLineWidth(4);
	JJY_FixGraph->SetLineColor(kGray);
	JJY_FixGraph->SetMarkerStyle(20);
	JJY_FixGraph->SetMarkerSize(5);
	JJY_FixGraph->SetMarkerColor(kGray);
	JJY_FixGraph->Draw("P same");

	TLegend Legend_JJY(0.50, 0.70, 0.89, 0.88);
	Legend_JJY.SetNColumns(2);
	Legend_JJY.AddEntry(JJY_LumiGraph, "Luminosity", "LP");
	Legend_JJY.AddEntry(JJY_BRGraph, "BR", "LP");
	Legend_JJY.AddEntry(JJY_CorrectionGraph, "Correction", "LP");
	Legend_JJY.AddEntry(JJY_CtauGraph, "c#tau shape", "LP");
	Legend_JJY.AddEntry(JJY_SigLxyGraph, "Lifetime variable", "LP");
	Legend_JJY.AddEntry(JJY_FitterGraph, "Fitter stability", "LP");
	Legend_JJY.AddEntry(JJY_FixGraph, "Parameter fix", "LP");
	Legend_JJY.AddEntry(JJY_TotalGraph, "Total", "LP");
	Legend_JJY.SetLineWidth(0);
	Legend_JJY.DrawClone();

	latex.SetTextSize(0.05);//0.5
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    // default is helvetica-italics
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.22, 0.80, "(d)");

	JJYCanvas->SaveAs("./SystematicPlot/JJY.png");

	//
	TCanvas* DeltaPhiCanvas = new TCanvas("DeltaPhiCanvas", "DeltaPhiCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	DeltaPhi_TotalGraph->SetTitle("");
	DeltaPhi_TotalGraph->GetYaxis()->SetTitle("Relative systematics (%)");
	DeltaPhi_TotalGraph->GetXaxis()->SetTitle("#Delta#phi(J/#psi_{1},J/#psi_{2})");
	DeltaPhi_TotalGraph->GetYaxis()->SetRangeUser(0, 45);
	DeltaPhi_TotalGraph->GetXaxis()->SetRangeUser(0, TMath::Pi());
	DeltaPhi_TotalGraph->SetLineWidth(6);
	DeltaPhi_TotalGraph->SetLineColor(kRed);
	DeltaPhi_TotalGraph->SetMarkerStyle(20);
	DeltaPhi_TotalGraph->SetMarkerSize(5);
	DeltaPhi_TotalGraph->SetMarkerColor(kRed);
	DeltaPhi_TotalGraph->Draw("AP same");
	DeltaPhi_LumiGraph->SetLineWidth(2);
	DeltaPhi_LumiGraph->SetLineColor(kBlack);
	DeltaPhi_LumiGraph->SetMarkerStyle(20);
	DeltaPhi_LumiGraph->SetMarkerSize(5);
	DeltaPhi_LumiGraph->SetMarkerColor(kBlack);
	DeltaPhi_LumiGraph->Draw("P same");
	DeltaPhi_BRGraph->SetLineWidth(4);
	DeltaPhi_BRGraph->SetLineColor(kBlue - 1);
	DeltaPhi_BRGraph->SetMarkerStyle(20);
	DeltaPhi_BRGraph->SetMarkerSize(5);
	DeltaPhi_BRGraph->SetMarkerColor(kBlue - 1);
	DeltaPhi_BRGraph->Draw("P same");
	DeltaPhi_CorrectionGraph->SetLineWidth(4);
	DeltaPhi_CorrectionGraph->SetLineColor(kBlue + 2);
	DeltaPhi_CorrectionGraph->SetMarkerStyle(20);
	DeltaPhi_CorrectionGraph->SetMarkerSize(5);
	DeltaPhi_CorrectionGraph->SetMarkerColor(kBlue + 2);
	DeltaPhi_CorrectionGraph->Draw("P same");
	DeltaPhi_CtauGraph->SetLineWidth(4);
	DeltaPhi_CtauGraph->SetLineColor(kBlue);
	DeltaPhi_CtauGraph->SetMarkerStyle(20);
	DeltaPhi_CtauGraph->SetMarkerSize(5);
	DeltaPhi_CtauGraph->SetMarkerColor(kBlue);
	DeltaPhi_CtauGraph->Draw("P same");
	DeltaPhi_SigLxyGraph->SetLineWidth(4);
	DeltaPhi_SigLxyGraph->SetLineColor(kBlue - 7);
	DeltaPhi_SigLxyGraph->SetMarkerStyle(20);
	DeltaPhi_SigLxyGraph->SetMarkerSize(5);
	DeltaPhi_SigLxyGraph->SetMarkerColor(kBlue - 7);
	DeltaPhi_SigLxyGraph->Draw("P same");
	DeltaPhi_FitterGraph->SetLineWidth(4);
	DeltaPhi_FitterGraph->SetLineColor(kBlue - 10);
	DeltaPhi_FitterGraph->SetMarkerStyle(20);
	DeltaPhi_FitterGraph->SetMarkerSize(5);
	DeltaPhi_FitterGraph->SetMarkerColor(kBlue - 10);
	DeltaPhi_FitterGraph->Draw("P same");
	DeltaPhi_FixGraph->SetLineWidth(4);
	DeltaPhi_FixGraph->SetLineColor(kGray);
	DeltaPhi_FixGraph->SetMarkerStyle(20);
	DeltaPhi_FixGraph->SetMarkerSize(5);
	DeltaPhi_FixGraph->SetMarkerColor(kGray);
	DeltaPhi_FixGraph->Draw("P same");

	TLegend Legend_DeltaPhi(0.50, 0.70, 0.89, 0.88);
	Legend_DeltaPhi.SetNColumns(2);
	Legend_DeltaPhi.AddEntry(DeltaPhi_LumiGraph, "Luminosity", "LP");
	Legend_DeltaPhi.AddEntry(DeltaPhi_BRGraph, "BR", "LP");
	Legend_DeltaPhi.AddEntry(DeltaPhi_CorrectionGraph, "Correction", "LP");
	Legend_DeltaPhi.AddEntry(DeltaPhi_CtauGraph, "c#tau shape", "LP");
	Legend_DeltaPhi.AddEntry(DeltaPhi_SigLxyGraph, "Lifetime variable", "LP");
	Legend_DeltaPhi.AddEntry(DeltaPhi_FitterGraph, "Fitter stability", "LP");
	Legend_DeltaPhi.AddEntry(DeltaPhi_FixGraph, "Parameter fix", "LP");
	Legend_DeltaPhi.AddEntry(DeltaPhi_TotalGraph, "Total", "LP");
	Legend_DeltaPhi.SetLineWidth(0);
	Legend_DeltaPhi.DrawClone();

	latex.SetTextSize(0.05);//0.5
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    // default is helvetica-italics
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.22, 0.80, "(e)");

	DeltaPhiCanvas->SaveAs("./SystematicPlot/DeltaPhi.png");

	//
	TCanvas* JPtCanvas = new TCanvas("JPtCanvas", "JPtCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	JPt_TotalGraph->SetTitle("");
	JPt_TotalGraph->GetYaxis()->SetTitle("Relative systematics (%)");
	JPt_TotalGraph->GetXaxis()->SetTitle("p_{T}(J/#psi) (GeV)");
	JPt_TotalGraph->GetYaxis()->SetRangeUser(0, 60);
	JPt_TotalGraph->GetXaxis()->SetRangeUser(10, 40);
	JPt_TotalGraph->SetLineWidth(6);
	JPt_TotalGraph->SetLineColor(kRed);
	JPt_TotalGraph->SetMarkerStyle(20);
	JPt_TotalGraph->SetMarkerSize(5);
	JPt_TotalGraph->SetMarkerColor(kRed);
	JPt_TotalGraph->Draw("AP same");
	JPt_LumiGraph->SetLineWidth(2);
	JPt_LumiGraph->SetLineColor(kBlack);
	JPt_LumiGraph->SetMarkerStyle(20);
	JPt_LumiGraph->SetMarkerSize(5);
	JPt_LumiGraph->SetMarkerColor(kBlack);
	JPt_LumiGraph->Draw("P same");
	JPt_BRGraph->SetLineWidth(4);
	JPt_BRGraph->SetLineColor(kBlue - 1);
	JPt_BRGraph->SetMarkerStyle(20);
	JPt_BRGraph->SetMarkerSize(5);
	JPt_BRGraph->SetMarkerColor(kBlue - 1);
	JPt_BRGraph->Draw("P same");
	JPt_CorrectionGraph->SetLineWidth(4);
	JPt_CorrectionGraph->SetLineColor(kBlue + 2);
	JPt_CorrectionGraph->SetMarkerStyle(20);
	JPt_CorrectionGraph->SetMarkerSize(5);
	JPt_CorrectionGraph->SetMarkerColor(kBlue + 2);
	JPt_CorrectionGraph->Draw("P same");
	JPt_CtauGraph->SetLineWidth(4);
	JPt_CtauGraph->SetLineColor(kBlue);
	JPt_CtauGraph->SetMarkerStyle(20);
	JPt_CtauGraph->SetMarkerSize(5);
	JPt_CtauGraph->SetMarkerColor(kBlue);
	JPt_CtauGraph->Draw("P same");
	JPt_SigLxyGraph->SetLineWidth(4);
	JPt_SigLxyGraph->SetLineColor(kBlue - 7);
	JPt_SigLxyGraph->SetMarkerStyle(20);
	JPt_SigLxyGraph->SetMarkerSize(5);
	JPt_SigLxyGraph->SetMarkerColor(kBlue - 7);
	JPt_SigLxyGraph->Draw("P same");
	JPt_FitterGraph->SetLineWidth(4);
	JPt_FitterGraph->SetLineColor(kBlue - 10);
	JPt_FitterGraph->SetMarkerStyle(20);
	JPt_FitterGraph->SetMarkerSize(5);
	JPt_FitterGraph->SetMarkerColor(kBlue - 10);
	JPt_FitterGraph->Draw("P same");
	JPt_FixGraph->SetLineWidth(4);
	JPt_FixGraph->SetLineColor(kGray);
	JPt_FixGraph->SetMarkerStyle(20);
	JPt_FixGraph->SetMarkerSize(5);
	JPt_FixGraph->SetMarkerColor(kGray);
	JPt_FixGraph->Draw("P same");

	TLegend Legend_JPt(0.50, 0.70, 0.89, 0.88);
	Legend_JPt.SetNColumns(2);
	Legend_JPt.AddEntry(JPt_LumiGraph, "Luminosity", "LP");
	Legend_JPt.AddEntry(JPt_BRGraph, "BR", "LP");
	Legend_JPt.AddEntry(JPt_CorrectionGraph, "Correction", "LP");
	Legend_JPt.AddEntry(JPt_CtauGraph, "c#tau shape", "LP");
	Legend_JPt.AddEntry(JPt_SigLxyGraph, "Lifetime variable", "LP");
	Legend_JPt.AddEntry(JPt_FitterGraph, "Fitter stability", "LP");
	Legend_JPt.AddEntry(JPt_FixGraph, "Parameter fix", "LP");
	Legend_JPt.AddEntry(JPt_TotalGraph, "Total", "LP");
	Legend_JPt.SetLineWidth(0);
	Legend_JPt.DrawClone();

	latex.SetTextSize(0.05);//0.5
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    // default is helvetica-italics
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.22, 0.80, "(f)");

	JPtCanvas->SaveAs("./SystematicPlot/JPt.png");

}
