//Jinfeng, 2025.09.25
//To draw the distribution of differential cross section

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


void DifferentialPainter() {

	//Scale factor (16)
	double ScaleFactor = 0.001 / (36.303 * 0.05961 * 0.05961);

	//Bin edge (do not add the "overflow edge" here)
	float JJMassBinEdge[] = { 7.5,17.5,27.5,37.5,47.5,57.5,67.5,77.5 };
	float DeltaYBinEdge[] = { 0.0,0.5,1.0,1.5,2.0,2.5 };
	float JJPtBinEdge[] = { 0,5,10,15,20,25,30,35,40 };
	float JJYBinEdge[] = { 0,0.4,0.8,1.2,1.6,2.0 };
	float DeltaPhiBinEdge[] = { 0.0,TMath::Pi() / 8,TMath::Pi() / 4,TMath::Pi() * 3 / 8,TMath::Pi() / 2,TMath::Pi() * 5 / 8,TMath::Pi() * 3 / 4,TMath::Pi() * 7 / 8,TMath::Pi() };
	float JPtBinEdge[] = { 10,15,20,25,30,35,40 };

	//Histogram
	float JJMassBinContent[8] = { 2858,1684,1438,709,280,178,94,105 }; //Event number from the fitting
	float JJMassBinError_Statistic[8] = { 58,45,45,32,21,17,11,13 }; //(1)
	float JJMassBinError_Systematic[8] = { 15.007, 17.800, 18.686, 14.461, 21.236, 15.553, 28.470, 28.190 }; //(%)
	float JJMassBinError_Total[8];
	for (int i = 0; i < 8; i++) JJMassBinError_Total[i] = TMath::Sqrt(TMath::Sq(JJMassBinContent[i] * JJMassBinError_Systematic[i] / 100) + TMath::Sq(JJMassBinError_Statistic[i]));

	float DeltaYBinContent[6] = { 3738,1281,842,662,479,338 };
	float DeltaYBinError_Statistic[6] = { 67,41,34,31,26,22 };
	float DeltaYBinError_Systematic[6] = { 14.937, 22.277, 25.670, 22.536, 27.699, 38.302 };
	float DeltaYBinError_Total[6];
	for (int i = 0; i < 6; i++) DeltaYBinError_Total[i] = TMath::Sqrt(TMath::Sq(DeltaYBinContent[i] * DeltaYBinError_Systematic[i] / 100) + TMath::Sq(DeltaYBinError_Statistic[i]));

	float JJPtBinContent[9] = { 1043,926,748,576,1805,1131,419,308,387 };
	float JJPtBinError_Statistic[9] = { 37,36,32,29,47,38,23,19,22 };
	float JJPtBinError_Systematic[9] = { 18.062, 27.105, 18.705, 21.361, 16.246, 17.868, 12.782, 21.939, 33.512 };
	float JJPtBinError_Total[9];
	for (int i = 0; i < 9; i++) JJPtBinError_Total[i] = TMath::Sqrt(TMath::Sq(JJPtBinContent[i] * JJPtBinError_Systematic[i] / 100) + TMath::Sq(JJPtBinError_Statistic[i]));

	//No overflow bin for Y(JJ)
	float JJYBinContent[5] = { 2079,2161,1679,1134,353 };
	float JJYBinError_Statistic[5] = { 53,54,47,38,22 };
	float JJYBinError_Systematic[5] = { 16.578, 16.737, 17.288, 14.752, 18.240 };
	float JJYBinError_Total[5];
	for (int i = 0; i < 5; i++) JJYBinError_Total[i] = TMath::Sqrt(TMath::Sq(JJYBinContent[i] * JJYBinError_Systematic[i] / 100) + TMath::Sq(JJYBinError_Statistic[i]));

	//No overflow bin for DeltaPhi
	float DeltaPhiBinContent[8] = { 2024,1122,557,562,473,604,786,710 };
	float DeltaPhiBinError_Statistic[8] = { 49,37,26,28,26,29,33,31 };
	float DeltaPhiBinError_Systematic[8] = { 18.856, 13.889, 22.245, 32.611, 21.656, 22.841, 25.113, 24.403 };
	float DeltaPhiBinError_Total[8];
	for (int i = 0; i < 8; i++) DeltaPhiBinError_Total[i] = TMath::Sqrt(TMath::Sq(DeltaPhiBinContent[i] * DeltaPhiBinError_Systematic[i] / 100) + TMath::Sq(DeltaPhiBinError_Statistic[i]));

	//No overflow bin for Pt(J)
	float JPtBinContent[6] = { 5351,1407,475,163,24,6 };
	float JPtBinError_Statistic[6] = { 83,43,25,15,6,1 };
	float JPtBinError_Systematic[6] = { 18.045, 14.826, 16.910, 23.485, 31.321, 38.208 };
	float JPtBinError_Total[6];
	for (int i = 0; i < 6; i++) JPtBinError_Total[i] = TMath::Sqrt(TMath::Sq(JPtBinContent[i] * JPtBinError_Systematic[i] / 100) + TMath::Sq(JPtBinError_Statistic[i]));

	//Histogram to display statistic uncertainty
	auto JJMassHistogram_Statistic = new TH1F("JJMassHistogram_Statistic", "", 7, JJMassBinEdge);
	JJMassHistogram_Statistic->SetStats(kFALSE);
	JJMassHistogram_Statistic->Sumw2();

	//Histogram to retreve systematic uncertainty
	//No plot will be produced 
	auto JJMassHistogram_Systematic = new TH1F("JJMassHistogram_Systematic", "", 7, JJMassBinEdge);
	JJMassHistogram_Systematic->Sumw2();

	//Histogram to display total uncertainty
	auto JJMassHistogram = new TH1F("JJMassHistogram", "", 7, JJMassBinEdge);
	JJMassHistogram->SetStats(kFALSE);
	JJMassHistogram->Sumw2();

	auto DeltaYHistogram_Statistic = new TH1F("DeltaYHistogram_Statistic", "", 5, DeltaYBinEdge);
	DeltaYHistogram_Statistic->SetStats(kFALSE);
	DeltaYHistogram_Statistic->Sumw2();
	auto DeltaYHistogram_Systematic = new TH1F("DeltaYHistogram_Systematic", "", 5, DeltaYBinEdge);
	DeltaYHistogram_Systematic->Sumw2();
	auto DeltaYHistogram = new TH1F("DeltaYHistogram", "", 5, DeltaYBinEdge);
	DeltaYHistogram->SetStats(kFALSE);
	DeltaYHistogram->Sumw2();

	auto JJPtHistogram_Statistic = new TH1F("JJPtHistogram_Statistic", "", 8, JJPtBinEdge);
	JJPtHistogram_Statistic->SetStats(kFALSE);
	JJPtHistogram_Statistic->Sumw2();
	auto JJPtHistogram_Systematic = new TH1F("JJPtHistogram_Systematic", "", 8, JJPtBinEdge);
	JJPtHistogram_Systematic->Sumw2();
	auto JJPtHistogram = new TH1F("JJPtHistogram", "", 8, JJPtBinEdge);
	JJPtHistogram->SetStats(kFALSE);
	JJPtHistogram->Sumw2();

	auto JJYHistogram_Statistic = new TH1F("JJYHistogram_Statistic", "", 5, JJYBinEdge);
	JJYHistogram_Statistic->SetStats(kFALSE);
	JJYHistogram_Statistic->Sumw2();
	auto JJYHistogram_Systematic = new TH1F("JJYHistogram_Systematic", "", 5, JJYBinEdge);
	JJYHistogram_Systematic->Sumw2();
	auto JJYHistogram = new TH1F("JJYHistogram", "", 5, JJYBinEdge);
	JJYHistogram->SetStats(kFALSE);
	JJYHistogram->Sumw2();

	auto DeltaPhiHistogram_Statistic = new TH1F("DeltaPhiHistogram_Statistic", "", 8, DeltaPhiBinEdge);
	DeltaPhiHistogram_Statistic->SetStats(kFALSE);
	DeltaPhiHistogram_Statistic->Sumw2();
	auto DeltaPhiHistogram_Systematic = new TH1F("DeltaPhiHistogram_Systematic", "", 8, DeltaPhiBinEdge);
	DeltaPhiHistogram_Systematic->Sumw2();
	auto DeltaPhiHistogram = new TH1F("DeltaPhiHistogram", "", 8, DeltaPhiBinEdge);
	DeltaPhiHistogram->SetStats(kFALSE);
	DeltaPhiHistogram->Sumw2();

	auto JPtHistogram_Statistic = new TH1F("JPtHistogram_Statistic", "", 6, JPtBinEdge);
	JPtHistogram_Statistic->SetStats(kFALSE);
	JPtHistogram_Statistic->Sumw2();
	auto JPtHistogram_Systematic = new TH1F("JPtHistogram_Systematic", "", 6, JPtBinEdge);
	JPtHistogram_Systematic->Sumw2();
	auto JPtHistogram = new TH1F("JPtHistogram", "", 6, JPtBinEdge);
	JPtHistogram->SetStats(kFALSE);
	JPtHistogram->Sumw2();

	//Fill histograms
	for (int i = 1; i < 9; i++)
	{
		JJMassHistogram_Statistic->SetBinContent(i, JJMassBinContent[i - 1]);
		JJMassHistogram_Statistic->SetBinError(i, JJMassBinError_Statistic[i - 1]);
		JJMassHistogram_Systematic->SetBinContent(i, JJMassBinContent[i - 1]);
		JJMassHistogram_Systematic->SetBinError(i, JJMassBinError_Systematic[i - 1] * JJMassBinContent[i - 1] / 100);
		JJMassHistogram->SetBinContent(i, JJMassBinContent[i - 1]);
		JJMassHistogram->SetBinError(i, JJMassBinError_Total[i - 1]);
	}

	for (int i = 1; i < 7; i++)
	{
		DeltaYHistogram_Statistic->SetBinContent(i, DeltaYBinContent[i - 1]);
		DeltaYHistogram_Statistic->SetBinError(i, DeltaYBinError_Statistic[i - 1]);
		DeltaYHistogram_Systematic->SetBinContent(i, DeltaYBinContent[i - 1]);
		DeltaYHistogram_Systematic->SetBinError(i, DeltaYBinError_Systematic[i - 1] * DeltaYBinContent[i - 1] / 100);
		DeltaYHistogram->SetBinContent(i, DeltaYBinContent[i - 1]);
		DeltaYHistogram->SetBinError(i, DeltaYBinError_Total[i - 1]);
	}
		
	for (int i = 1; i < 10; i++)
	{
		JJPtHistogram_Statistic->SetBinContent(i, JJPtBinContent[i - 1]);
		JJPtHistogram_Statistic->SetBinError(i, JJPtBinError_Statistic[i - 1]);
		JJPtHistogram_Systematic->SetBinContent(i, JJPtBinContent[i - 1]);
		JJPtHistogram_Systematic->SetBinError(i, JJPtBinError_Systematic[i - 1] * JJPtBinContent[i - 1] / 100);
		JJPtHistogram->SetBinContent(i, JJPtBinContent[i - 1]);
		JJPtHistogram->SetBinError(i, JJPtBinError_Total[i - 1]);
	}

	for (int i = 1; i < 6; i++)
	{
		JJYHistogram_Statistic->SetBinContent(i, JJYBinContent[i - 1]);
		JJYHistogram_Statistic->SetBinError(i, JJYBinError_Statistic[i - 1]);
		JJYHistogram_Systematic->SetBinContent(i, JJYBinContent[i - 1]);
		JJYHistogram_Systematic->SetBinError(i, JJYBinError_Systematic[i - 1] * JJYBinContent[i - 1] / 100);
		JJYHistogram->SetBinContent(i, JJYBinContent[i - 1]);
		JJYHistogram->SetBinError(i, JJYBinError_Total[i - 1]);
	}

	for (int i = 1; i < 9; i++)
	{
		DeltaPhiHistogram_Statistic->SetBinContent(i, DeltaPhiBinContent[i - 1]);
		DeltaPhiHistogram_Statistic->SetBinError(i, DeltaPhiBinError_Statistic[i - 1]);
		DeltaPhiHistogram_Systematic->SetBinContent(i, DeltaPhiBinContent[i - 1]);
		DeltaPhiHistogram_Systematic->SetBinError(i, DeltaPhiBinError_Systematic[i - 1] * DeltaPhiBinContent[i - 1] / 100);
		DeltaPhiHistogram->SetBinContent(i, DeltaPhiBinContent[i - 1]);
		DeltaPhiHistogram->SetBinError(i, DeltaPhiBinError_Total[i - 1]);
	}

	for (int i = 1; i < 7; i++)
	{
		JPtHistogram_Statistic->SetBinContent(i, JPtBinContent[i - 1]);
		JPtHistogram_Statistic->SetBinError(i, JPtBinError_Statistic[i - 1]);
		JPtHistogram_Systematic->SetBinContent(i, JPtBinContent[i - 1]);
		JPtHistogram_Systematic->SetBinError(i, JPtBinError_Systematic[i - 1] * JPtBinContent[i - 1] / 100);
		JPtHistogram->SetBinContent(i, JPtBinContent[i - 1]);
		JPtHistogram->SetBinError(i, JPtBinError_Total[i - 1]);
	}

	//Calculate the differential cross section
	JJMassHistogram_Statistic->Scale(ScaleFactor * 0.1);
	JJMassHistogram_Systematic->Scale(ScaleFactor * 0.1);
	JJMassHistogram->Scale(ScaleFactor * 0.1);

	DeltaYHistogram_Statistic->Scale(ScaleFactor * 2);
	DeltaYHistogram_Systematic->Scale(ScaleFactor * 2);
	DeltaYHistogram->Scale(ScaleFactor * 2);

	JJPtHistogram_Statistic->Scale(ScaleFactor / 5);
	JJPtHistogram_Systematic->Scale(ScaleFactor / 5);
	JJPtHistogram->Scale(ScaleFactor / 5);

	JJYHistogram_Statistic->Scale(ScaleFactor / 0.4);
	JJYHistogram_Systematic->Scale(ScaleFactor / 0.4);
	JJYHistogram->Scale(ScaleFactor / 0.4);

	DeltaPhiHistogram_Statistic->Scale(ScaleFactor * 8 / TMath::Pi());
	DeltaPhiHistogram_Systematic->Scale(ScaleFactor * 8 / TMath::Pi());
	DeltaPhiHistogram->Scale(ScaleFactor * 8 / TMath::Pi());

	JPtHistogram_Statistic->Scale(ScaleFactor / 5);
	JPtHistogram_Systematic->Scale(ScaleFactor / 5);
	JPtHistogram->Scale(ScaleFactor / 5);
	
	//Print the result
	cout << "///////////////////" << endl;

	cout << "M(JJ): " << endl;
	for (int i = 1; i < 9; i++) 
	{
		cout << JJMassBinEdge[i - 1] << "-" << JJMassBinEdge[i] << " & ";
		cout << JJMassHistogram->GetBinContent(i) << " \\pm ";
		cout << JJMassHistogram_Statistic->GetBinError(i) << " \\pm ";
		cout << JJMassHistogram_Systematic->GetBinError(i) << "\\\\" << endl;
	}

	cout << "DeltaY: " << endl;
	for (int i = 1; i < 7; i++)
	{
		cout << DeltaYBinEdge[i - 1] << "-" << DeltaYBinEdge[i] << " & ";
		cout << DeltaYHistogram->GetBinContent(i) << " \\pm ";
		cout << DeltaYHistogram_Statistic->GetBinError(i) << " \\pm ";
		cout << DeltaYHistogram_Systematic->GetBinError(i) << "\\\\" << endl;
	}

	cout << "Pt(JJ): " << endl;
	for (int i = 1; i < 10; i++)
	{
		cout << JJPtBinEdge[i - 1] << "-" << JJPtBinEdge[i] << " & ";
		cout << JJPtHistogram->GetBinContent(i) << " \\pm ";
		cout << JJPtHistogram_Statistic->GetBinError(i) << " \\pm ";
		cout << JJPtHistogram_Systematic->GetBinError(i) << "\\\\" << endl;
	}
		
	cout << "Y(JJ): " << endl;
	for (int i = 1; i < 6; i++)
	{
		cout << JJYBinEdge[i - 1] << "-" << JJYBinEdge[i] << " & ";
		cout << JJYHistogram->GetBinContent(i) << " \\pm ";
		cout << JJYHistogram_Statistic->GetBinError(i) << " \\pm ";
		cout << JJYHistogram_Systematic->GetBinError(i) << "\\\\" << endl;
	}
		
	cout << "DeltaPhi: " << endl;
	for (int i = 1; i < 9; i++)
	{
		cout << DeltaPhiBinEdge[i - 1] << "-" << DeltaPhiBinEdge[i] << " & ";
		cout << DeltaPhiHistogram->GetBinContent(i) << " \\pm ";
		cout << DeltaPhiHistogram_Statistic->GetBinError(i) << " \\pm ";
		cout << DeltaPhiHistogram_Systematic->GetBinError(i) << "\\\\" << endl;
	}

	cout << "Pt(J): " << endl;
	for (int i = 1; i < 7; i++)
	{
		cout << JPtBinEdge[i - 1] << "-" << JPtBinEdge[i] << " & ";
		cout << JPtHistogram->GetBinContent(i) << " \\pm ";
		cout << JPtHistogram_Statistic->GetBinError(i) << " \\pm ";
		cout << JPtHistogram_Systematic->GetBinError(i) << "\\\\" << endl;
	}
		
	//Draw
	//JJMass
	TCanvas* JJMassCanvas = new TCanvas("JJMassCanvas", "JJMassCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	JJMassHistogram_Statistic->GetYaxis()->SetTitle("d#sigma/dM_{JJ} (pb/GeV)");
	JJMassHistogram_Statistic->GetYaxis()->SetRangeUser(0, 2.8);
	JJMassHistogram_Statistic->GetXaxis()->SetTitle("M_{J/#psi_{1}J/#psi_{2}} (GeV)");
	JJMassHistogram_Statistic->GetXaxis()->SetRangeUser(7.5, 80);
	JJMassHistogram_Statistic->SetLineWidth(4);
	JJMassHistogram_Statistic->SetFillColor(40);
	JJMassHistogram_Statistic->Draw("E2 same"); //Error is drawn as a filled area
	JJMassHistogram->GetXaxis()->SetRangeUser(7.5, 80);
	JJMassHistogram->SetLineColor(kBlack);
	JJMassHistogram->SetLineWidth(2);
	JJMassHistogram->SetMarkerStyle(20);
	JJMassHistogram->SetMarkerSize(5);
	JJMassHistogram->Draw("E same");

	TLatex overflow;
	overflow.SetTextSize(0.03);
	overflow.SetTextColor(40);
	overflow.SetTextAlign(22);
	overflow.DrawLatex(80.5, 0.8, "Overflow");

	TArrow arrow;
	arrow.SetLineColor(40);
	arrow.SetFillColor(40);
	arrow.DrawArrow(82.5, 0.7, 82.5, 0.2, 0.03, "|>");

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.05);
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.88, 0.85, "(a)");

	JJMassCanvas->SaveAs("DifferentialPlot/Differential_JJMass.pdf");

	//DeltaY
	TCanvas* DeltaYCanvas = new TCanvas("DeltaYCanvas", "DeltaYCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	DeltaYHistogram_Statistic->GetYaxis()->SetTitle("d#sigma/d|#Delta y| (pb)");
	DeltaYHistogram_Statistic->GetYaxis()->SetRangeUser(4, 75);
	DeltaYHistogram_Statistic->GetXaxis()->SetTitle("#Deltay(J/#psi_{1},J/#psi_{2})");
	DeltaYHistogram_Statistic->GetXaxis()->SetRangeUser(0, 3.0);
	DeltaYHistogram_Statistic->SetLineWidth(4);
	DeltaYHistogram_Statistic->SetFillColor(40);
	DeltaYHistogram_Statistic->Draw("E2 same");
	DeltaYHistogram->GetXaxis()->SetRangeUser(0, 3.0);
	DeltaYHistogram->SetLineColor(kBlack);
	DeltaYHistogram->SetLineWidth(2);
	DeltaYHistogram->SetMarkerStyle(20);
	DeltaYHistogram->SetMarkerSize(5);
	DeltaYHistogram->Draw("E same");

	overflow.DrawLatex(2.75, 20.0, "Overflow");
	arrow.DrawArrow(2.75, 18.0, 2.75, 8.5, 0.03, "|>");

	latex.SetTextSize(0.05);
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.88, 0.85, "(b)");

	DeltaYCanvas->SaveAs("DifferentialPlot/Differential_DeltaY.pdf");

	//JJPt
	TCanvas* JJPtCanvas = new TCanvas("JJPtCanvas", "JJPtCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	JJPtHistogram_Statistic->GetYaxis()->SetTitle("d#sigma/dp_{T}_{JJ} (pb/GeV)");
	JJPtHistogram_Statistic->GetYaxis()->SetRangeUser(0, 3.8);
	JJPtHistogram_Statistic->GetXaxis()->SetTitle("p_{T}(J/#psi_{1}J/#psi_{2}) (GeV)");
	JJPtHistogram_Statistic->GetXaxis()->SetRangeUser(0, 45);
	JJPtHistogram_Statistic->SetLineWidth(4);
	JJPtHistogram_Statistic->SetFillColor(40);
	JJPtHistogram_Statistic->Draw("E2 same");
	JJPtHistogram->GetXaxis()->SetRangeUser(0, 45);
	JJPtHistogram->SetLineColor(kBlack);
	JJPtHistogram->SetMarkerStyle(20);
	JJPtHistogram->SetLineWidth(2);
	JJPtHistogram->SetMarkerSize(5);
	JJPtHistogram->Draw("E same");

	overflow.DrawLatex(41.0, 1.5, "Overflow");
	arrow.DrawArrow(42.5, 1.4, 42.5, 1.0, 0.03, "|>");

	latex.SetTextSize(0.05);
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.88, 0.85, "(c)");

	JJPtCanvas->SaveAs("DifferentialPlot/Differential_JJPt.pdf");

	//JJY
	TCanvas* JJYCanvas = new TCanvas("JJYCanvas", "JJYCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	JJYHistogram_Statistic->GetYaxis()->SetTitle("d#sigma/d|y_{JJ}| (pb)");
	JJYHistogram_Statistic->GetYaxis()->SetRangeUser(0, 60);
	JJYHistogram_Statistic->GetXaxis()->SetTitle("|y(J/#psi_{1}J/#psi_{2})|");
	JJYHistogram_Statistic->SetLineWidth(4);
	JJYHistogram_Statistic->SetFillColor(40);
	JJYHistogram_Statistic->Draw("E2 same");
	JJYHistogram->SetLineColor(kBlack);
	JJYHistogram->SetLineWidth(2);
	JJYHistogram->SetMarkerStyle(20);
	JJYHistogram->SetMarkerSize(5);
	JJYHistogram->Draw("E same");

	latex.SetTextSize(0.05);
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.88, 0.85, "(d)");

	JJYCanvas->SaveAs("DifferentialPlot/Differential_JJY.pdf");

	//DeltaPhi
	TCanvas* DeltaPhiCanvas = new TCanvas("DeltaPhiCanvas", "DeltaPhiCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	DeltaPhiHistogram_Statistic->GetYaxis()->SetTitle("d#sigma/d|#Delta#phi| (pb)");
	DeltaPhiHistogram_Statistic->GetYaxis()->SetRangeUser(5, 55);
	DeltaPhiHistogram_Statistic->GetXaxis()->SetTitle("#Delta#phi(J/#psi_{1},J/#psi_{2})");
	DeltaPhiHistogram_Statistic->SetLineWidth(4);
	DeltaPhiHistogram_Statistic->SetFillColor(40);
	DeltaPhiHistogram_Statistic->Draw("E2 same");
	DeltaPhiHistogram->SetLineColor(kBlack);
	DeltaPhiHistogram->SetLineWidth(2);
	DeltaPhiHistogram->SetMarkerStyle(20);
	DeltaPhiHistogram->SetMarkerSize(5);
	DeltaPhiHistogram->Draw("E same");

	latex.SetTextSize(0.05);
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                    
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.88, 0.85, "(e)");

	DeltaPhiCanvas->SaveAs("DifferentialPlot/Differential_DeltaPhi.pdf");

	//JPt
	TCanvas* JPtCanvas = new TCanvas("JPtCanvas", "JPtCanvas", 3000, 3000);
	gPad->SetLeftMargin(0.15);
	gPad->SetLogy();
	JPtHistogram_Statistic->GetYaxis()->SetTitle("d#sigma/dp_{T}(J) (pb/GeV)");
	JPtHistogram_Statistic->GetYaxis()->SetRangeUser(0.005, 20);
	JPtHistogram_Statistic->GetXaxis()->SetTitle("p_{T}(J/#psi) (GeV)");
	JPtHistogram_Statistic->SetLineWidth(4);
	JPtHistogram_Statistic->SetFillColor(40);
	JPtHistogram_Statistic->Draw("E2 same");
	JPtHistogram->SetLineColor(kBlack);
	JPtHistogram->SetLineWidth(2);
	JPtHistogram->SetMarkerStyle(20);
	JPtHistogram->SetMarkerSize(5);
	JPtHistogram->Draw("E same");

	latex.SetTextSize(0.05);
	latex.SetTextFont(61);
	latex.SetTextAlign(33);
	latex.DrawLatex(0.28, 0.88, "CMS");

	latex.SetTextFont(52);                   
	latex.SetTextAlign(33);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.48, 0.87, "Preliminary");

	latex.SetTextFont(42);
	latex.SetTextAlign(31);
	latex.SetTextSize(0.04); 
	latex.DrawLatex(0.88, 0.92, "36.3 fb^{-1} (13 TeV)");
	latex.DrawLatex(0.88, 0.85, "(f)");

	JPtCanvas->SaveAs("DifferentialPlot/Differential_JPt.pdf");

}
