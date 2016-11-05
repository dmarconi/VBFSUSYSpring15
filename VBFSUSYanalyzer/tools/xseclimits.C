#include <TROOT.h>
#include <TPad.h>
#include <TFrame.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <cmath>

void setTDRStyle_mod() {
	  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

	  // For the canvas:
	  //tdrStyle->SetCanvasBorderMode(0);
	  //tdrStyle->SetCanvasColor(kWhite);
	  //tdrStyle->SetCanvasDefH(600); //Height of canvas
	  //tdrStyle->SetCanvasDefW(600); //Width of canvas
	  //tdrStyle->SetCanvasDefX(0);   //POsition on screen
	  //tdrStyle->SetCanvasDefY(0);

	  // For the Pad:
	  //tdrStyle->SetPadBorderMode(0);
	  // tdrStyle->SetPadBorderSize(Width_t size = 1);
	  //tdrStyle->SetPadColor(kWhite);
	  //tdrStyle->SetPadGridX(false);
	  //tdrStyle->SetPadGridY(false);
	  //tdrStyle->SetGridColor(0);
	  //tdrStyle->SetGridStyle(3);
	  //tdrStyle->SetGridWidth(1);

	  // For the frame:
	  //tdrStyle->SetFrameBorderMode(0);
	  //tdrStyle->SetFrameBorderSize(1);
	  //tdrStyle->SetFrameFillColor(0);
	  //tdrStyle->SetFrameFillStyle(0);
	  //tdrStyle->SetFrameLineColor(1);
	  //tdrStyle->SetFrameLineStyle(1);
	  //tdrStyle->SetFrameLineWidth(1);

	  // For the histo:
	  // tdrStyle->SetHistFillColor(1);
	  // tdrStyle->SetHistFillStyle(0);
	  //tdrStyle->SetHistLineColor(1);
	  //tdrStyle->SetHistLineStyle(0);
	  //tdrStyle->SetHistLineWidth(1);
	  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
	  // tdrStyle->SetNumberContours(Int_t number = 20);

	  //tdrStyle->SetEndErrorSize(2);
	  // tdrStyle->SetErrorMarker(20);
	  //tdrStyle->SetErrorX(0.);

	  //tdrStyle->SetMarkerStyle(20);

	  //For the fit/function:
	  //tdrStyle->SetOptFit(1);
	  //tdrStyle->SetFitFormat("5.4g");
	  //tdrStyle->SetFuncColor(2);
	  //tdrStyle->SetFuncStyle(1);
	  //tdrStyle->SetFuncWidth(1);

	  //For the date:
	  //tdrStyle->SetOptDate(0);
	  // tdrStyle->SetDateX(Float_t x = 0.01);
	  // tdrStyle->SetDateY(Float_t y = 0.01);

	  // For the statistics box:
	  //tdrStyle->SetOptFile(0);
	  //tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
	  //tdrStyle->SetStatColor(kWhite);
	  //tdrStyle->SetStatFont(42);
	  //tdrStyle->SetStatFontSize(0.025);
	  //tdrStyle->SetStatTextColor(1);
	  //tdrStyle->SetStatFormat("6.4g");
	  //tdrStyle->SetStatBorderSize(1);
	  //tdrStyle->SetStatH(0.1);
	  //tdrStyle->SetStatW(0.15);
	  // tdrStyle->SetStatStyle(Style_t style = 1001);
	  // tdrStyle->SetStatX(Float_t x = 0);
	  // tdrStyle->SetStatY(Float_t y = 0);

	  // Margins:
	  //tdrStyle->SetPadTopMargin(0.05);
	  //tdrStyle->SetPadBottomMargin(0.13);
	  //tdrStyle->SetPadLeftMargin(0.16);
	  //tdrStyle->SetPadRightMargin(0.02);
	  tdrStyle->SetPadRightMargin(0.12);

	  // For the Global title:

	  //tdrStyle->SetOptTitle(0);
	  //tdrStyle->SetTitleFont(42);
	  //tdrStyle->SetTitleColor(1);
	  //tdrStyle->SetTitleTextColor(1);
	  //tdrStyle->SetTitleFillColor(10);
	  //tdrStyle->SetTitleFontSize(0.05);
	  // tdrStyle->SetTitleH(0); // Set the height of the title box
	  // tdrStyle->SetTitleW(0); // Set the width of the title box
	  // tdrStyle->SetTitleX(0); // Set the position of the title box
	  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
	  // tdrStyle->SetTitleStyle(Style_t style = 1001);
	  // tdrStyle->SetTitleBorderSize(2);

	  // For the axis titles:

	  //tdrStyle->SetTitleColor(1, "XYZ");
	  //tdrStyle->SetTitleFont(42, "XYZ");
	  //tdrStyle->SetTitleSize(0.06, "XYZ");
	  //tdrStyle->SetTitleSize(0.03, "XYZ");
	  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
	  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
	  //tdrStyle->SetTitleXOffset(0.9);
	  //tdrStyle->SetTitleYOffset(1.25);
	  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

	  // For the axis labels:

	  //tdrStyle->SetLabelColor(1, "XYZ");
	  //tdrStyle->SetLabelFont(42, "XYZ");
	  //tdrStyle->SetLabelOffset(0.007, "XYZ");
	  //tdrStyle->SetLabelSize(0.05, "XYZ");

	  // For the axis:

	  //tdrStyle->SetAxisColor(1, "XYZ");
	 //tdrStyle->SetStripDecimals(kTRUE);
	  //tdrStyle->SetTickLength(0.03, "XYZ");
	  //tdrStyle->SetNdivisions(510, "XYZ");
	  //tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	  //tdrStyle->SetPadTickY(1);

	  // Change for log plots:
	  //tdrStyle->SetOptLogx(0);
	  //tdrStyle->SetOptLogy(0);
	  //tdrStyle->SetOptLogz(0);

	  // Postscript options:
	  //tdrStyle->SetPaperSize(20.,20.);
	  // tdrStyle->SetLineScalePS(Float_t scale = 3);
	  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
	  // tdrStyle->SetHeaderPS(const char* header);
	  // tdrStyle->SetTitlePS(const char* pstitle);

	  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
	  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
	  // tdrStyle->SetPaintTextFormat(const char* format = "g");
	  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
	  // tdrStyle->SetTimeOffset(Double_t toffset);
	  // tdrStyle->SetHistMinimumZero(kTRUE);

	  //tdrStyle->SetHatchesLineWidth(5);
	  //tdrStyle->SetHatchesSpacing(0.05);

	  tdrStyle->cd();

}

void setTDRStyle() {
	  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

	  // For the canvas:
	  tdrStyle->SetCanvasBorderMode(0);
	  tdrStyle->SetCanvasColor(kWhite);
	  tdrStyle->SetCanvasDefH(600); //Height of canvas
	  tdrStyle->SetCanvasDefW(600); //Width of canvas
	  tdrStyle->SetCanvasDefX(0);   //POsition on screen
	  tdrStyle->SetCanvasDefY(0);

	  // For the Pad:
	  tdrStyle->SetPadBorderMode(0);
	  // tdrStyle->SetPadBorderSize(Width_t size = 1);
	  tdrStyle->SetPadColor(kWhite);
	  tdrStyle->SetPadGridX(false);
	  tdrStyle->SetPadGridY(false);
	  tdrStyle->SetGridColor(0);
	  tdrStyle->SetGridStyle(3);
	  tdrStyle->SetGridWidth(1);

	  // For the frame:
	  tdrStyle->SetFrameBorderMode(0);
	  tdrStyle->SetFrameBorderSize(1);
	  tdrStyle->SetFrameFillColor(0);
	  tdrStyle->SetFrameFillStyle(0);
	  tdrStyle->SetFrameLineColor(1);
	  tdrStyle->SetFrameLineStyle(1);
	  tdrStyle->SetFrameLineWidth(1);

	  // For the histo:
	  // tdrStyle->SetHistFillColor(1);
	  // tdrStyle->SetHistFillStyle(0);
	  tdrStyle->SetHistLineColor(1);
	  tdrStyle->SetHistLineStyle(0);
	  tdrStyle->SetHistLineWidth(1);
	  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
	  // tdrStyle->SetNumberContours(Int_t number = 20);

	  tdrStyle->SetEndErrorSize(2);
	  // tdrStyle->SetErrorMarker(20);
	  //tdrStyle->SetErrorX(0.);

	  tdrStyle->SetMarkerStyle(20);

	  //For the fit/function:
	  tdrStyle->SetOptFit(1);
	  tdrStyle->SetFitFormat("5.4g");
	  tdrStyle->SetFuncColor(2);
	  tdrStyle->SetFuncStyle(1);
	  tdrStyle->SetFuncWidth(1);

	  //For the date:
	  tdrStyle->SetOptDate(0);
	  // tdrStyle->SetDateX(Float_t x = 0.01);
	  // tdrStyle->SetDateY(Float_t y = 0.01);

	  // For the statistics box:
	  tdrStyle->SetOptFile(0);
	  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
	  tdrStyle->SetStatColor(kWhite);
	  tdrStyle->SetStatFont(42);
	  tdrStyle->SetStatFontSize(0.025);
	  tdrStyle->SetStatTextColor(1);
	  tdrStyle->SetStatFormat("6.4g");
	  tdrStyle->SetStatBorderSize(1);
	  tdrStyle->SetStatH(0.1);
	  tdrStyle->SetStatW(0.15);
	  // tdrStyle->SetStatStyle(Style_t style = 1001);
	  // tdrStyle->SetStatX(Float_t x = 0);
	  // tdrStyle->SetStatY(Float_t y = 0);

	  // Margins:
	  tdrStyle->SetPadTopMargin(0.05);
	  tdrStyle->SetPadBottomMargin(0.13);
	  tdrStyle->SetPadLeftMargin(0.16);
	  //tdrStyle->SetPadRightMargin(0.02);
	  tdrStyle->SetPadRightMargin(0.12);

	  // For the Global title:

	  tdrStyle->SetOptTitle(0);
	  tdrStyle->SetTitleFont(42);
	  tdrStyle->SetTitleColor(1);
	  tdrStyle->SetTitleTextColor(1);
	  tdrStyle->SetTitleFillColor(10);
	  tdrStyle->SetTitleFontSize(0.05);
	  // tdrStyle->SetTitleH(0); // Set the height of the title box
	  // tdrStyle->SetTitleW(0); // Set the width of the title box
	  // tdrStyle->SetTitleX(0); // Set the position of the title box
	  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
	  // tdrStyle->SetTitleStyle(Style_t style = 1001);
	  // tdrStyle->SetTitleBorderSize(2);

	  // For the axis titles:

	  tdrStyle->SetTitleColor(1, "XYZ");
	  tdrStyle->SetTitleFont(42, "XYZ");
	  //tdrStyle->SetTitleSize(0.06, "XYZ");
	  tdrStyle->SetTitleSize(0.03, "XYZ");
	  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
	  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
	  tdrStyle->SetTitleXOffset(0.9);
	  tdrStyle->SetTitleYOffset(1.25);
	  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

	  // For the axis labels:

	  tdrStyle->SetLabelColor(1, "XYZ");
	  tdrStyle->SetLabelFont(42, "XYZ");
	  tdrStyle->SetLabelOffset(0.007, "XYZ");
	  tdrStyle->SetLabelSize(0.05, "XYZ");

	  // For the axis:

	  tdrStyle->SetAxisColor(1, "XYZ");
	  tdrStyle->SetStripDecimals(kTRUE);
	  tdrStyle->SetTickLength(0.03, "XYZ");
	  tdrStyle->SetNdivisions(510, "XYZ");
	  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
	  tdrStyle->SetPadTickY(1);

	  // Change for log plots:
	  tdrStyle->SetOptLogx(0);
	  tdrStyle->SetOptLogy(0);
	  tdrStyle->SetOptLogz(0);

	  // Postscript options:
	  tdrStyle->SetPaperSize(20.,20.);
	  // tdrStyle->SetLineScalePS(Float_t scale = 3);
	  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
	  // tdrStyle->SetHeaderPS(const char* header);
	  // tdrStyle->SetTitlePS(const char* pstitle);

	  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
	  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
	  // tdrStyle->SetPaintTextFormat(const char* format = "g");
	  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
	  // tdrStyle->SetTimeOffset(Double_t toffset);
	  // tdrStyle->SetHistMinimumZero(kTRUE);

	  tdrStyle->SetHatchesLineWidth(5);
	  tdrStyle->SetHatchesSpacing(0.05);

	  tdrStyle->cd();

}

double vbfefficiency(double evenCRcounts, double oddCRcounts){

	return ( (oddCRcounts) / ( (oddCRcounts) + (evenCRcounts) ));
}

double vbfConversionFactor(string taupt, string isoregion) {

	TFile *inputfile = TFile::Open((taupt + "/allQCD_"+ taupt +".root").c_str());
	TH1F* h_ditaucharge;
	TH1F* h_ditauchargeVBFinverted;

	h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	//h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h_ditaucharge").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

	double counts = h_ditaucharge->GetBinContent(3); 
	double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3);


	double vbfeff = vbfefficiency(counts_inverted, counts);

	delete inputfile;
	return (vbfeff / (1. - vbfeff));


}

double LtoTfactor(string taupt) {

	TFile *inputfile = TFile::Open(( taupt + "/allQCD_"+ taupt +"_LtoT.root").c_str());

	TH1F* h1_counts;
	h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

	double oneJetcounts = h1_counts->GetBinContent(2);
	double oneTauMatchcounts = h1_counts->GetBinContent(3);
	double jetAlsoTcounts = h1_counts->GetBinContent(4);
	double fourJetscounts = h1_counts->GetBinContent(5);
	double oneLooseTau = h1_counts->GetBinContent(6);
	double oneLcounts = h1_counts->GetBinContent(7);
	double alsoTcounts = h1_counts->GetBinContent(8);

	double jetToTightProb = jetAlsoTcounts / oneTauMatchcounts;
	double looseToTightProb = alsoTcounts / oneLcounts;
	double twoLooseTo2Tight = jetToTightProb * jetToTightProb;

	delete inputfile;
	return (twoLooseTo2Tight);	

}


TH2F* makeEffPlot(string taupt, string isoregion, string chi, string lsp) {
	//TFile *inputfile = TFile::Open("VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M.root");
	//TFile *inputfile = TFile::Open("VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau195_Chargino200_1M.root");
	TFile *inputfile;
	if ((chi == "chi100") && (lsp == "lsp000")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau095_Chargino100_1M.root").c_str());
	if ((chi == "chi200") && (lsp == "lsp000")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau195_Chargino200_1M.root").c_str());
	if ((chi == "chi300") && (lsp == "lsp000")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M.root").c_str());
	if ((chi == "chi100") && (lsp == "lsp050")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau095_Chargino100_1M.root").c_str());
	if ((chi == "chi200") && (lsp == "lsp050")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau195_Chargino200_1M.root").c_str());
	if ((chi == "chi300") && (lsp == "lsp050")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau295_Chargino300_1M.root").c_str());
	TH2F* h2_DiJetInvMass_vs_MET;
	TH1F* h_count; 
	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h2_DiJetInvMass_vs_MET").c_str())));
	//h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	h_count = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX(); 
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1);
	//cout << "N total events: " << ntotalevents << endl;
	//double ntotalevents = h_ditaucharge->GetBinContent(3);

	TH2F* h2_DiJetInvMass_vs_MET_eff;
	h2_DiJetInvMass_vs_MET_eff = new TH2F ("h2_DiJetInvMass_vs_MET_eff","h2_DiJetInvMass_vs_MET_eff", nbinsx, 0., 240., nbinsy , 0., 2500.);		
	h2_DiJetInvMass_vs_MET_eff->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_eff->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_eff->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double integral = h2_DiJetInvMass_vs_MET->Integral( i, nbinsx, j, nbinsy );
			h2_DiJetInvMass_vs_MET_eff->SetBinContent(i,j, (integral/ntotalevents));

		}
	}
	TCanvas *my_canvas = new TCanvas;
	my_canvas->cd();
	gPad->SetLogz();
	h2_DiJetInvMass_vs_MET_eff->Draw("colz");
	my_canvas->Print(("results/JetInvMass_vs_MET_eff_" + isoregion + "_" + taupt + ".pdf").c_str());
	my_canvas->Close();
	//h2_DiJetInvMass_vs_MET_eff->Clear();
	//delete inputfile;
	return h2_DiJetInvMass_vs_MET_eff;
}

TH2F* makeBackgroundPlot_LtoT(string taupt, string isoregion){


	TFile *inputfile = TFile::Open((taupt + "/allQCD_"+ taupt +".root").c_str());
	//TFile *inputfile = TFile::Open("allQCD.root");
	//TFile *inputfile = TFile::Open("allQCD_taupt20.root");
	TH2F* h2_DiJetInvMass_vs_MET;
	TH2F* h2_DiJetInvMass_vs_MET_LtoT;
	TH1F* h_ditauchargeVBFinverted;
	TH1F* h_count;
	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h2_DiJetInvMass_vs_MET").c_str())));
	//h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h2_DiJetInvMass_vs_MET").c_str())));
	//h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	h_count = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	//h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h_ditaucharge").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));
	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX(); 
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1);
	//double ntotalnumevents = h2_DiJetInvMass_vs_MET->Integral( 0. , nbinsx, 0. , nbinsy );
	//double ntotalevents = h_ditaucharge->GetBinContent(3);

	//TODO convert to the tight plot
	h2_DiJetInvMass_vs_MET_LtoT = new TH2F ("h2_DiJetInvMass_vs_MET_LtoT","h2_DiJetInvMass_vs_MET_LtoT", nbinsx, 0., 240., nbinsy , 0., 2500.);		
	h2_DiJetInvMass_vs_MET_LtoT->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->SetStats(0);

	TH2F* h2_DiJetInvMass_vs_MET_eff;
	h2_DiJetInvMass_vs_MET_eff = new TH2F ("h2_DiJetInvMass_vs_MET_eff","h2_DiJetInvMass_vs_MET_eff", nbinsx, 0., 240., nbinsy , 0., 2500.);		
	h2_DiJetInvMass_vs_MET_eff->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_eff->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_eff->SetStats(0);

	//TODO fill the 2 tight map applying the 2fold method for every bin

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double bincontent = h2_DiJetInvMass_vs_MET->GetBinContent (i, j);
			double vbffactor = vbfConversionFactor(taupt, isoregion);
			double ltotfactor = LtoTfactor(taupt);
			h2_DiJetInvMass_vs_MET_LtoT->SetBinContent(i,j, (bincontent*vbffactor*ltotfactor));

		}
	}
	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double integral = h2_DiJetInvMass_vs_MET_LtoT->Integral( i, nbinsx, j, nbinsy );
			h2_DiJetInvMass_vs_MET_eff->SetBinContent(i,j, (integral/ntotalevents));

		}
	}


	TCanvas *my_canvas = new TCanvas;
	my_canvas->cd();
	gPad->SetLogz();
	double maximum = h2_DiJetInvMass_vs_MET->GetMaximum();
	double minimum = h2_DiJetInvMass_vs_MET->GetMinimum();
	h2_DiJetInvMass_vs_MET_LtoT->GetZaxis()->SetRangeUser(0.000001 , 10000);
	h2_DiJetInvMass_vs_MET_LtoT->SetStats(0);
	h2_DiJetInvMass_vs_MET_LtoT->Draw("colz");
	my_canvas->Print(("JetInvMass_vs_MET_" + isoregion + "_" + taupt + "_LtoT.pdf").c_str());
	gPad->SetLogy();
	TH1D *p_met = h2_DiJetInvMass_vs_MET->ProjectionX();
	p_met->Draw();
	my_canvas->Print("p_met_original.root");
	my_canvas->Close();
	h2_DiJetInvMass_vs_MET_eff->Clear();

	return h2_DiJetInvMass_vs_MET_LtoT;
}

double getSignificance(double signal, double background, double sigma) { 
	return (   (signal) / ( (sigma * 0.5) + sqrt( background + (0.5 * background * background)))  );
}

double getXSection(double efficiency, double luminosity, double background, double sigma, double significance) { 
	return ( (significance * ( (0.5 * sigma) + sqrt( background + (0.5 * background * background))  ) ) / (efficiency * luminosity)   );
}

double getSignalEvents(int xbin, int ybin, TH2F* signal_map, double xsec, double lumi){
	double efficiency = signal_map->GetBinContent(xbin,ybin);
	return (efficiency * xsec * lumi);
}

double getSignalEfficiency(int xbin, int ybin, TH2F* signal_map){
	double efficiency = signal_map->GetBinContent(xbin,ybin);
	return (efficiency);
}


double getBackgroudEvents(int xbin, int ybin, TH2F* background_map){
	int nbinsx = background_map->GetNbinsX(); 
	int nbinsy = background_map->GetNbinsY();
	double events = background_map->Integral( xbin, nbinsx, ybin, nbinsy );
	return (events); 
}


void makeSignificance(string chi, string lsp) {

	double xsec = 0.003106; //pb only for test purposes
	double lumi = 85000.; 
	TH2F* h2_DiJetInvMass_vs_MET_eff_signal;
	TH2F* h2_DiJetInvMass_vs_MET_background;
	h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot("taupt20", "Taui2TightIso", "chi100","lsp000");
	h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT("taupt20", "Tau2LooseIsoInclusive");



	int nbinsx = h2_DiJetInvMass_vs_MET_background->GetNbinsX(); 
	int nbinsy = h2_DiJetInvMass_vs_MET_background->GetNbinsY();

	TH2F* h2_DiJetInvMass_vs_MET_significance;
	h2_DiJetInvMass_vs_MET_significance = new TH2F ("h2_DiJetInvMass_vs_MET_significance","h2_DiJetInvMass_vs_MET_significance", nbinsx, 0., 240., nbinsy , 0., 2500.);		
	h2_DiJetInvMass_vs_MET_significance->GetZaxis()->SetTitle("Significance");
	h2_DiJetInvMass_vs_MET_significance->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_significance->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_significance->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double sensitivity = getSignificance(
					getSignalEvents( i , j, h2_DiJetInvMass_vs_MET_eff_signal, xsec, lumi), 
					getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background),
					2.5 //sigma
					);
			h2_DiJetInvMass_vs_MET_significance->SetBinContent(i, j, sensitivity);
		}
	}

	//paleta(h2_DiJetInvMass_vs_MET_significance);

	TCanvas *my_canvas = new TCanvas;
	my_canvas->cd();
	gPad->SetLogz();

	h2_DiJetInvMass_vs_MET_significance->GetZaxis()->SetRangeUser(0.0000001,10);
	h2_DiJetInvMass_vs_MET_significance->Draw("colz");
	//my_canvas->Print(("JetInvMass_vs_MET_significance_" + isoregion + "_" + taupt + "_LtoT.pdf").c_str());
	//my_canvas->Print("JetInvMass_vs_MET_significance_Tau2TightIso_taupt20.pdf");
	my_canvas->Print("results/JetInvMass_vs_MET_significance_Tau2TightIso_taupt20.root");
	my_canvas->Close();
}

void makeXSection(string taupt,string chi, string lsp) {

	double lumi = 85000.; 
	TH2F* h2_DiJetInvMass_vs_MET_eff_signal;
	TH2F* h2_DiJetInvMass_vs_MET_background;
	h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot(taupt, "Taui2TightIso", chi, lsp);
	h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT(taupt, "baseline");

	int nbinsx = h2_DiJetInvMass_vs_MET_background->GetNbinsX(); 
	int nbinsy = h2_DiJetInvMass_vs_MET_background->GetNbinsY();

	TH1D *p_met = h2_DiJetInvMass_vs_MET_background->ProjectionX();


	TH2F* h2_DiJetInvMass_vs_MET_xsec;
	TH2F* h2_DiJetInvMass_vs_MET_xsec_zoom;
	h2_DiJetInvMass_vs_MET_xsec = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(),("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(), nbinsx, 0., 240., nbinsy , 0., 2500.);		
	h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetTitle("#sigma pb^{-1} [GeV]");
	h2_DiJetInvMass_vs_MET_xsec->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_xsec->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double xsec = getXSection(
					getSignalEfficiency(i,j, h2_DiJetInvMass_vs_MET_eff_signal), 
					lumi, 
					getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background), 
					2.5, //sigma
					2.0 //significance
					); 


			h2_DiJetInvMass_vs_MET_xsec->SetBinContent(i, j, xsec);
		}
	}



	TCanvas *my_canvas_zoom = new TCanvas("mycanvas_zoom","mycanvas_zoom",1024.,768.);
	my_canvas_zoom->cd();
	gPad->SetLogz();
	h2_DiJetInvMass_vs_MET_xsec_zoom = (TH2F*) h2_DiJetInvMass_vs_MET_xsec->Clone();
	h2_DiJetInvMass_vs_MET_xsec_zoom->GetXaxis()->SetRangeUser(70.,200.);		
	h2_DiJetInvMass_vs_MET_xsec_zoom->GetYaxis()->SetRangeUser(0.,1000.);		
	double maximum_zoom = h2_DiJetInvMass_vs_MET_xsec_zoom->GetMaximum();
	double minimum_zoom = h2_DiJetInvMass_vs_MET_xsec_zoom->GetMinimum();
	h2_DiJetInvMass_vs_MET_xsec_zoom ->GetZaxis()->SetRangeUser(minimum_zoom,maximum_zoom);		
	h2_DiJetInvMass_vs_MET_xsec_zoom->Draw("colz");
	gPad->Modified();
	gPad->Update();
	my_canvas_zoom->Print(("results/JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt + "_zoom.pdf").c_str());

	TCanvas *my_canvas = new TCanvas("mycanvas","mycanvas",1024.,768.);
	my_canvas->cd();
	gPad->SetLogz();

	int binmin = h2_DiJetInvMass_vs_MET_xsec->GetMinimumBin();
	double minimum = 100000000000000000.;
	double x_min = 0.;
	double y_min = 0.;
	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			if ((h2_DiJetInvMass_vs_MET_xsec->GetBinContent(i,j) < minimum) && (h2_DiJetInvMass_vs_MET_xsec->GetBinContent(i,j) > 0)){
				minimum = h2_DiJetInvMass_vs_MET_xsec->GetBinContent(i,j);	
				x_min = h2_DiJetInvMass_vs_MET_xsec->GetXaxis()->GetBinCenter(i);
				y_min = h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->GetBinCenter(j);
			}
		}
	}

	cout << "Minimum at " << minimum << " for " << chi << ", " << lsp << " and " << taupt << ": Jetinvmass " << y_min << " MET " << x_min << endl; 

	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetRangeUser(0.001,1000);
	h2_DiJetInvMass_vs_MET_xsec->Draw("colz");
	gPad->Modified();
	gPad->Update();
	gPad->ResizePad();
	//setTDRStyle_mod();
	//gPad->SetRightMargin(0.04);
	gStyle->SetPadRightMargin(0.12);
	
	my_canvas->SetRightMargin(0.04);
	//my_canvas->Print(("results/JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt + ".pdf").c_str());
	my_canvas->Print(("results/JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt + ".root").c_str());
	gPad->SetLogy();
	p_met->Draw();
	my_canvas->Print(("results/p_met" + taupt + ".root").c_str());
	my_canvas->Close();
}

void makeXsecLimPlots(string chi, string lsp){


	makeXSection("taupt20",chi,lsp);
	makeXSection("taupt25",chi,lsp);
	makeXSection("taupt30",chi,lsp);
	makeXSection("taupt35",chi,lsp);
	makeXSection("taupt40",chi,lsp);
	makeXSection("taupt45",chi,lsp);

}

void fullXsecLimScan(){

	makeXsecLimPlots("chi100", "lsp000");
	makeXsecLimPlots("chi200", "lsp000");
	makeXsecLimPlots("chi300", "lsp000");
	makeXsecLimPlots("chi100", "lsp050");
	makeXsecLimPlots("chi200", "lsp050");
	makeXsecLimPlots("chi300", "lsp050");

}
