//Jinfeng, 2025.09.25
//To acquire the DPS fraction by a simultaneous fit

#include "TCanvas.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TTree.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

#include "RooSimultaneous.h"
#include "RooCategory.h"

using namespace RooFit;

void SimFit()
{
	//Scale factor (16)
	double ScaleFactor = 0.001 / (36.303 * 0.05961 * 0.05961);

	//Bin edge
	float DeltaYBinEdge[] = { 0.0,0.5,1.0,1.5,2.0,2.5 };
	float DeltaPhiBinEdge[] = { 0.0,TMath::Pi() / 8,TMath::Pi() / 4,TMath::Pi() * 3 / 8,TMath::Pi() / 2,TMath::Pi() * 5 / 8,TMath::Pi() * 3 / 4,TMath::Pi() * 7 / 8,TMath::Pi() };

	//Data histogram
	float DataDeltaYBinContent[6] = { 3738,1281,842,662,479,338 };
	float DataDeltaYBinError_Statistic[6] = { 67,41,34,31,26,22 };
	float DataDeltaYBinError_Systematic[6] = { 14.937, 22.277, 25.670, 22.536, 27.699, 38.302 };
	float DataDeltaYBinError_Total[6];
	for (int i = 0; i < 6; i++) DataDeltaYBinError_Total[i] = TMath::Sqrt(TMath::Sq(DataDeltaYBinContent[i] * DataDeltaYBinError_Systematic[i] / 100) + TMath::Sq(DataDeltaYBinError_Statistic[i]));

	auto DataDeltaYHistogram_Total = new TH1F("DataDeltaYHistogram_Total", "", 5, DeltaYBinEdge);
	DataDeltaYHistogram_Total->SetStats(kFALSE);
	DataDeltaYHistogram_Total->Sumw2();
	
	///////////////
	float DataDeltaPhiBinContent[8] = { 2024,1122,557,562,473,604,786,710 };
	float DataDeltaPhiBinError_Statistic[8] = { 49,37,26,28,26,29,33,31 };
	float DataDeltaPhiBinError_Systematic[8] = { 18.856, 13.889, 22.245, 32.611, 21.656, 22.841, 25.113, 24.403 };
	float DataDeltaPhiBinError_Total[8];
	for (int i = 0; i < 8; i++) DataDeltaPhiBinError_Total[i] = TMath::Sqrt(TMath::Sq(DataDeltaPhiBinContent[i] * DataDeltaPhiBinError_Systematic[i] / 100) + TMath::Sq(DataDeltaPhiBinError_Statistic[i]));

	auto DataDeltaPhiHistogram_Total = new TH1F("DataDeltaPhiHistogram_Total", "", 8, DeltaPhiBinEdge);
	DataDeltaPhiHistogram_Total->SetStats(kFALSE);
	DataDeltaPhiHistogram_Total->Sumw2();
	
	//Fill data histograms
	for (int i = 1; i < 7; i++)
	{
		DataDeltaYHistogram_Total->SetBinContent(i, DataDeltaYBinContent[i - 1]);
		DataDeltaYHistogram_Total->SetBinError(i, DataDeltaYBinError_Total[i - 1]);
	}

	for (int i = 1; i < 9; i++)
	{
		DataDeltaPhiHistogram_Total->SetBinContent(i, DataDeltaPhiBinContent[i - 1]);
		DataDeltaPhiHistogram_Total->SetBinError(i, DataDeltaPhiBinError_Total[i - 1]);
	}

	DataDeltaYHistogram_Total->Scale(ScaleFactor * 2);
	DataDeltaPhiHistogram_Total->Scale(ScaleFactor * 8 / TMath::Pi());

	//SPS&DPS histograms
	//Parameters
	float SPS_J1Mass;
	float SPS_J2Mass;
	float SPS_DeltaY;
	float SPS_DeltaPhi;
	float SPS_Weight;

	float DPS_J1Mass;
	float DPS_J2Mass;
	float DPS_DeltaY;
	float DPS_DeltaPhi;
	float DPS_Weight;

	//MC file acquired
	TFile* SPS_File = TFile::Open("../Sample/MC_SPS.root");
	TTree* SPS_Tree = new TTree("SPS_Tree", "SPS_Tree");
	SPS_File->GetObject("WeightedTree", SPS_Tree);
	SPS_Tree->SetBranchAddress("J1Mass_", &SPS_J1Mass);
	SPS_Tree->SetBranchAddress("J2Mass_", &SPS_J2Mass);
	SPS_Tree->SetBranchAddress("DeltaY_", &SPS_DeltaY);
	SPS_Tree->SetBranchAddress("DeltaPhi_", &SPS_DeltaPhi);
	SPS_Tree->SetBranchAddress("Weight_sum", &SPS_Weight);

	TFile* DPS_File = TFile::Open("../Sample/MC_DPS.root");
	TTree* DPS_Tree = new TTree("DPS_Tree", "DPS_Tree");
	DPS_File->GetObject("WeightedTree", DPS_Tree);
	DPS_Tree->SetBranchAddress("J1Mass_", &DPS_J1Mass);
	DPS_Tree->SetBranchAddress("J2Mass_", &DPS_J2Mass);
	DPS_Tree->SetBranchAddress("DeltaY_", &DPS_DeltaY);
	DPS_Tree->SetBranchAddress("DeltaPhi_", &DPS_DeltaPhi);
	DPS_Tree->SetBranchAddress("Weight_sum", &DPS_Weight);

	int SPS_Length = SPS_Tree->GetEntries();
	int DPS_Length = DPS_Tree->GetEntries();

	auto SPSDeltaYHistogram = new TH1F("SPSDeltaYHistogram", "", 5, DeltaYBinEdge);
	auto DPSDeltaYHistogram = new TH1F("DPSDeltaYHistogram", "", 5, DeltaYBinEdge);

	auto SPSDeltaPhiHistogram = new TH1F("SPSDeltaPhiHistogram", "", 8, DeltaPhiBinEdge);
	auto DPSDeltaPhiHistogram = new TH1F("DPSDeltaPhiHistogram", "", 8, DeltaPhiBinEdge);

	SPSDeltaYHistogram->Sumw2();
	DPSDeltaYHistogram->Sumw2();
	SPSDeltaPhiHistogram->Sumw2();
	DPSDeltaPhiHistogram->Sumw2();

	SPSDeltaYHistogram->SetStats(kFALSE);
	DPSDeltaYHistogram->SetStats(kFALSE);
	SPSDeltaPhiHistogram->SetStats(kFALSE);
	DPSDeltaPhiHistogram->SetStats(kFALSE);

	//Fill MC histograms
	for (int i = 0; i < SPS_Length; i++) 
	{
		SPS_Tree->GetEntry(i);

		if ((SPS_J1Mass > 2.95) && (SPS_J1Mass < 3.25) && (SPS_J2Mass > 2.95) && (SPS_J2Mass < 3.25) && (SPS_Weight < 5000)) 
		{
			SPSDeltaYHistogram->Fill(SPS_DeltaY, SPS_Weight);
			SPSDeltaPhiHistogram->Fill(SPS_DeltaPhi, SPS_Weight);
		}
	}

	SPSDeltaYHistogram->Scale(ScaleFactor * 2);
	SPSDeltaPhiHistogram->Scale(ScaleFactor * 8 / TMath::Pi());

	for (int i = 0; i < DPS_Length; i++)
	{
		DPS_Tree->GetEntry(i);

		if ((DPS_J1Mass > 2.95) && (DPS_J1Mass < 3.25) && (DPS_J2Mass > 2.95) && (DPS_J2Mass < 3.25) && (DPS_Weight < 5000))
		{
			DPSDeltaYHistogram->Fill(DPS_DeltaY, DPS_Weight);
			DPSDeltaPhiHistogram->Fill(DPS_DeltaPhi, DPS_Weight);
		}
	}

	DPSDeltaYHistogram->Scale(ScaleFactor * 2);
	DPSDeltaPhiHistogram->Scale(ScaleFactor * 8 / TMath::Pi());

	//Prepare the templates
	RooRealVar DeltaY("DeltaY", "#Deltay(J/#psi_{1},J/#psi_{2})", 0, 2.5);

	RooDataHist SPSDeltaYTemplate("SPSDeltaYTemplate", "SPSDeltaYTemplate", DeltaY, SPSDeltaYHistogram);
	RooDataHist DPSDeltaYTemplate("DPSDeltaTemplate", "DPSDeltaTemplate", DeltaY, DPSDeltaYHistogram);
	RooDataHist DataDeltaYTemplate("DataDeltaTemplate", "DataDeltaTemplate", DeltaY, DataDeltaYHistogram_Total);

	RooHistPdf SPSDeltaYPDF("SPSDeltaYPDF", "SPSDeltaYPDF", RooArgSet(DeltaY), SPSDeltaYTemplate);
	RooHistPdf DPSDeltaYPDF("DPSDeltaYPDF", "DPSDeltaYPDF", RooArgSet(DeltaY), DPSDeltaYTemplate);

	////////////////
	RooRealVar DeltaPhi("DeltaPhi", "#Delta#phi(J/#psi_{1}, J/#psi_{2})", 0, TMath::Pi());

	RooDataHist SPSDeltaPhiTemplate("SPSDeltaPhiTemplate", "SPSDeltaPhiTemplate", DeltaPhi, SPSDeltaPhiHistogram);
	RooDataHist DPSDeltaPhiTemplate("DPSDeltaPhiTemplate", "DPSDeltaPhiTemplate", DeltaPhi, DPSDeltaPhiHistogram);
	RooDataHist DataDeltaPhiTemplate("DataDeltaPhiTemplate", "DataDeltaPhiTemplate", DeltaPhi, DataDeltaPhiHistogram_Total);

	RooHistPdf SPSDeltaPhiPDF("SPSDeltaPhiPDF", "SPSDeltaPhiPDF", RooArgSet(DeltaPhi), SPSDeltaPhiTemplate);
	RooHistPdf DPSDeltaPhiPDF("DPSDeltaPhiPDF", "DPSDeltaPhiPDF", RooArgSet(DeltaPhi), DPSDeltaPhiTemplate);

	////////////////
	RooRealVar DPSFraction("DPSFraction", "DPSFraction", 0.5, 0, 1);

	RooAddPdf DeltaYPDF("DeltaYPDF", "DeltaYPDF", RooArgList(DPSDeltaYPDF, SPSDeltaYPDF), RooArgList(DPSFraction));
	RooAddPdf DeltaPhiPDF("DeltaPhiPDF", "DeltaPhiPDF", RooArgList(DPSDeltaPhiPDF, SPSDeltaPhiPDF), RooArgList(DPSFraction));

	//Prepare the category
	RooCategory DataCategory("DataCategory", "DataCategory");
	DataCategory.defineType("DeltaYCategory");
	DataCategory.defineType("DeltaPhiCategory");

	//Prepare the dataset
	RooDataHist DatasetDeltaYDeltaPhi("DatasetDeltaYDeltaPhi", "DatasetDeltaYDeltaPhi", RooArgList(DeltaY, DeltaPhi),
		Index(DataCategory),
		Import({ {"DeltaYCategory",&DataDeltaYTemplate},{"DeltaPhiCategory",&DataDeltaPhiTemplate} }));

	//Prepare the simultaneous PDF
	RooSimultaneous SimPDF("SimPDF", "SimPDF",
		{ {"DeltaYCategory",&DeltaYPDF},{"DeltaPhiCategory",&DeltaPhiPDF} }, DataCategory);

	//Fit
	RooFitResult* FitResult = SimPDF.fitTo(DatasetDeltaYDeltaPhi, Extended(kTRUE), Hesse(kTRUE), Strategy(1), NumCPU(4), Save(kTRUE), Minos(kFALSE));

	//Plot
	RooPlot* DeltaYFrame = DeltaY.frame(Title(" "));

	DatasetDeltaYDeltaPhi.plotOn(DeltaYFrame, Cut("DataCategory==DataCategory::DeltaYCategory"), LineColor(kBlack), LineWidth(2), MarkerSize(5), MarkerStyle(20), Name("Data"));

	SimPDF.plotOn(DeltaYFrame, Slice(DataCategory, "DeltaYCategory"), ProjWData(DataCategory, DatasetDeltaYDeltaPhi), LineColor(kRed), LineWidth(2), Name("Total"));
	SimPDF.plotOn(DeltaYFrame, Components(SPSDeltaYPDF), Slice(DataCategory, "DeltaYCategory"), ProjWData(DataCategory, DatasetDeltaYDeltaPhi), LineColor(kBlue), LineStyle(kDashed), LineWidth(2), Name("SPS"));
	SimPDF.plotOn(DeltaYFrame, Components(DPSDeltaYPDF), Slice(DataCategory, "DeltaYCategory"), ProjWData(DataCategory, DatasetDeltaYDeltaPhi), LineColor(40), DrawOption("F"), FillColor(40), MoveToBack(), Name("DPS"));

	RooHist* DeltaYPull = DeltaYFrame->pullHist("Data", "Total");
	DeltaYPull->SetMarkerStyle(20);
	DeltaYPull->SetMarkerSize(5);
	RooPlot* DeltaYPullFrame = DeltaY.frame(Title(" "));
	DeltaYPullFrame->addPlotable(DeltaYPull, "P");

	//Legend
	TLegend legend(0.65, 0.6, 0.85, 0.88);
	legend.AddEntry(DeltaYFrame->findObject("Data"), "Data", "LEP");

	legend.AddEntry(DeltaYFrame->findObject("Total"), "Total PDF", "L");
	legend.AddEntry(DeltaYFrame->findObject("SPS"), "SPS", "L");
	legend.AddEntry(DeltaYFrame->findObject("DPS"), "DPS", "F");
	legend.SetLineWidth(0);

	TCanvas* DeltaYCanvas = new TCanvas("DeltaYCanvas", "DeltaYCanvas", 3000, 3000);
	DeltaYCanvas->Divide(1, 2);
	DeltaYCanvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	DeltaYFrame->GetYaxis()->SetTitle("d#sigma/d|#Deltay| (pb)");
	DeltaYFrame->GetYaxis()->SetRangeUser(4, 75);
	DeltaYFrame->Draw("same");

	legend.DrawClone();

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.05);//0.5
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.25, 0.88, "CMS");

	latex.SetTextFont(52);                    
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.42, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");

	DeltaYCanvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	DeltaYPullFrame->GetYaxis()->SetTitle("Pull");
	DeltaYPullFrame->GetYaxis()->SetTitleOffset(0.4);
	DeltaYPullFrame->GetYaxis()->SetTitleSize(0.1);
	DeltaYPullFrame->GetYaxis()->SetLabelSize(0.1);
	DeltaYPullFrame->GetXaxis()->SetTitleSize(0.1);
	DeltaYPullFrame->GetXaxis()->SetLabelSize(0.1);
	DeltaYPullFrame->Draw();

	TLine* l = new TLine(0, 0, 2.5, 0);
	l->SetLineWidth(2);
	l->SetLineColor(kRed);
	l->SetLineStyle(kDashed);
	l->Draw("same");

	DeltaYCanvas->SaveAs("./Plot/Template_DeltaY.pdf");

	////////////////
	RooPlot* DeltaPhiFrame = DeltaPhi.frame(Title(" "));

	DatasetDeltaYDeltaPhi.plotOn(DeltaPhiFrame, Cut("DataCategory==DataCategory::DeltaPhiCategory"), LineColor(kBlack), LineWidth(2), MarkerSize(5), MarkerStyle(20), Name("Data"));

	SimPDF.plotOn(DeltaPhiFrame, Slice(DataCategory, "DeltaPhiCategory"), ProjWData(DataCategory, DatasetDeltaYDeltaPhi), LineColor(kRed), LineWidth(2), Name("Total"));
	SimPDF.plotOn(DeltaPhiFrame, Components(SPSDeltaPhiPDF), Slice(DataCategory, "DeltaPhiCategory"), ProjWData(DataCategory, DatasetDeltaYDeltaPhi), LineColor(kBlue), LineWidth(2), LineStyle(kDashed), Name("SPS"));
	SimPDF.plotOn(DeltaPhiFrame, Components(DPSDeltaPhiPDF), Slice(DataCategory, "DeltaPhiCategory"), ProjWData(DataCategory, DatasetDeltaYDeltaPhi), LineColor(40), DrawOption("F"), FillColor(40), MoveToBack(), Name("DPS"));

	//Pull
	RooHist* DeltaPhiPull = DeltaPhiFrame->pullHist("Data", "Total");
	DeltaPhiPull->SetMarkerStyle(20);
	DeltaPhiPull->SetMarkerSize(5);
	RooPlot* DeltaPhiPullFrame = DeltaPhi.frame(Title(" "));
	DeltaPhiPullFrame->addPlotable(DeltaPhiPull, "P");

	//DeltaPhi
	TCanvas* DeltaPhiCanvas = new TCanvas("DeltaPhiCanvas", "DeltaPhiCanvas", 3000, 3000);
	DeltaPhiCanvas->Divide(1, 2);
	DeltaPhiCanvas->cd(1)->SetPad(0.01, 0.20, 0.95, 0.98);
	gPad->SetLeftMargin(0.15);
	DeltaPhiFrame->GetYaxis()->SetRangeUser(5, 50);
	DeltaPhiFrame->GetYaxis()->SetTitle("d#sigma/d|#Delta#phi| (pb)");
	DeltaPhiFrame->Draw("same");
	legend.DrawClone();

	latex.SetTextSize(0.05);//0.5
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.25, 0.88, "CMS");

	latex.SetTextFont(52);                    
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.42, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); // 0.4
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");

	DeltaPhiCanvas->cd(2)->SetPad(0.01, 0.03, 0.95, 0.25);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.25);
	DeltaPhiPullFrame->GetYaxis()->SetTitle("Pull");
	DeltaPhiPullFrame->GetYaxis()->SetTitleOffset(0.4);
	DeltaPhiPullFrame->GetYaxis()->SetTitleSize(0.1);
	DeltaPhiPullFrame->GetYaxis()->SetLabelSize(0.1);
	DeltaPhiPullFrame->GetXaxis()->SetTitleSize(0.1);
	DeltaPhiPullFrame->GetXaxis()->SetLabelSize(0.1);
	DeltaPhiPullFrame->Draw();

	l = new TLine(0, 0, TMath::Pi(), 0);
	l->SetLineWidth(2);
	l->SetLineColor(kRed);
	l->SetLineStyle(kDashed);
	l->Draw("same");

	DeltaPhiCanvas->SaveAs("./Plot/Template_DeltaPhi.pdf");

	//Print
	cout << "DPS fraction: " << DPSFraction.getValV() << " +/- " << DPSFraction.getError() << endl;

}
