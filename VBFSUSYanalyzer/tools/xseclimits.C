#include <TROOT.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2D.h>
#include <iostream>
#include <fstream>
#include <cmath>


double vbfefficiency(double evenCRcounts, double oddCRcounts){

	        return ( (oddCRcounts) / ( (oddCRcounts) + (evenCRcounts) ));
}

double vbfConversionFactor(string taupt, string isoregion) {

	//TFile *inputfile = TFile::Open(("allQCD_"+ taupt +".root").c_str());
	TFile *inputfile = TFile::Open("allQCD.root");
	TH1F* h_ditaucharge;
	TH1F* h_ditauchargeVBFinverted;

	h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h_ditaucharge").c_str())));

	double counts = h_ditaucharge->GetBinContent(3); 
	double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3);


	double vbfeff = vbfefficiency(counts_inverted, counts);

	return (vbfeff / (1. - vbfeff));

}

double LtoTfactor(string taupt) {

	TFile *inputfile = TFile::Open(("allQCD_"+ taupt +"_LtoT.root").c_str());

	TH1F* h1_counts;
	h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

	double fourJetscounts = h1_counts->GetBinContent(2);
	double twoLooseTaus = h1_counts->GetBinContent(3);
	double oneLcounts = h1_counts->GetBinContent(4);
	double alsoTcounts = h1_counts->GetBinContent(5);

	//double fourJetscounts_err = binerrors(h1_counts, 2, 2);
	//double twoLooseTaus_err = binerrors(h1_counts, 3, 3);
	//double oneLcounts_err = binerrors(h1_counts, 4, 4);
	//double alsoTcounts_err = binerrors(h1_counts, 5, 5);
	
	double looseToTightProb = alsoTcounts / oneLcounts;
	//double looseToTightProb_err = sqrt( pow( (alsoTcounts_err / oneLcounts )   , 2.) + pow( ( (alsoTcounts * oneLcounts_err) / (oneLcounts * oneLcounts)     )    , 2.) );

	double twoLooseTo2Tight = looseToTightProb * looseToTightProb;
	//double twoLooseTo2Tight_err = 2. * looseToTightProb * looseToTightProb_err;

	return (twoLooseTo2Tight);	

}


TH2F* makeEffPlot(string taupt, string isoregion) {
	TFile *inputfile = TFile::Open("VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M.root");
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
	my_canvas->Print(("JetInvMass_vs_MET_eff_" + isoregion + "_" + taupt + ".pdf").c_str());
	my_canvas->Close();
	//h2_DiJetInvMass_vs_MET_eff->Clear();
	return h2_DiJetInvMass_vs_MET_eff;
}

TH2F* makeBackgroundPlot_LtoT(string taupt, string isoregion){


	//TFile *inputfile = TFile::Open(("allQCD_"+ taupt +".root").c_str());
	TFile *inputfile = TFile::Open("allQCD.root");
	TH2F* h2_DiJetInvMass_vs_MET;
	TH2F* h2_DiJetInvMass_vs_MET_LtoT;
	TH1F* h_ditauchargeVBFinverted;
	TH1F* h_count; 
	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h2_DiJetInvMass_vs_MET").c_str())));
	//h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	h_count = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h_ditaucharge").c_str())));
	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX(); 
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1);
	double ntotalnumevents = h2_DiJetInvMass_vs_MET->Integral( 0. , nbinsx, 0. , nbinsy );
	//cout << "N total events: " << ntotalevents << endl;
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
	//gPad->SetLogy();
	h2_DiJetInvMass_vs_MET_eff->Draw("colz");
	my_canvas->Print(("JetInvMass_vs_MET_eff_" + isoregion + "_" + taupt + "_LtoT.pdf").c_str());
	my_canvas->Close();
	h2_DiJetInvMass_vs_MET_eff->Clear();
	///cout << "Total denominator events for " << isoregion << " / "  <<taupt << " : " << ntotalevents << endl;
	//cout << "Total numerator events for " << isoregion << " / "  <<taupt << " : " << ntotalnumevents << endl;
	return h2_DiJetInvMass_vs_MET_LtoT;
}

double getSensitivity(double signal, double background) { 

	return ( (signal)/ (sqrt(background + pow((0.5 * background),2.)))  );
}

double getSignalEvents(int xbin, int ybin, TH2F* signal_map, double xsec, double lumi){
//TODO implement the proper efficiency input
double efficiency = signal_map->GetBinContent(xbin,ybin);
	return (efficiency * xsec * lumi);
}

double getBackgroudEvents(int xbin, int ybin, TH2F* background_map){
	//TODO implement the proper events input
	int nbinsx = background_map->GetNbinsX(); 
	int nbinsy = background_map->GetNbinsY();
	double events = background_map->Integral( xbin, nbinsx, ybin, nbinsy );
	if (events > 0.) {
		return (events); 
	} else {
		return 0.0000001;
	}
}

void paleta(TH2F* h) {
	//palette settings - completely independent
	const Int_t NRGBs = 6;
	const Int_t NCont = 999;

	Double_t stops[NRGBs] = { 0.00, 0.1, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs]   = { 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.0, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };


	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);

	//gStyle->SetOptStat(0);

	//here the actually interesting code starts
	const Double_t min = 0.000001;
	const Double_t max = 100000;

	//const Int_t nLevels = 999;
	//Double_t levels[nLevels];
	const Int_t nLevels = 15;
	Double_t levels[nLevels];

	  levels[0] = 0.000000001; 
	  levels[1] = 0.00000001; 
	  levels[2] = 0.0000001; 
	  levels[3] = 0.000001; 
	  levels[4] = 0.00001;
	  levels[5] = 0.0001;
	  levels[6] = 0.001;
	  levels[7] = 0.01; 
	  levels[8] = 0.1; 
	  levels[9] = 1.;
	  levels[10] = 10.; 
	  levels[11] = 100.; 
	  levels[12] = 1000.;
	  levels[13] = 10000.; 
	  levels[14] = 100000.; 

	//for(int i = 1; i < nLevels; i++) {
	//	levels[i] = min + (max - min) / (nLevels - 1) * (i);
	//}
	//levels[0] = 0.01;
	//  levels[0] = -1; //Interesting, but this also works as I want!

	TCanvas * c = new TCanvas();
	//TH2D *h  = new TH2D("h", "", 10, 0, 10, 10, 0, 10);
	h->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
	gPad->SetLogz();
	h->DrawClone("colz");// draw "axes", "contents", "statistics box"

	h->GetZaxis()->SetRangeUser(min, max); // ... set the range ...

	h->Draw("z same"); // draw the "color palette"
	c->SaveAs("c.png");
}

void xSecLimits() {

	double xsec = 0.003106; //pb only for test purposes
	double lumi = 85000.; 
	TH2F* h2_DiJetInvMass_vs_MET_eff_signal;
	TH2F* h2_DiJetInvMass_vs_MET_background;
	h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot("taupt20", "Taui2TightIso");
	h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT("taupt20", "Tau2LooseIsoInclusive");

	int nbinsx = h2_DiJetInvMass_vs_MET_background->GetNbinsX(); 
	int nbinsy = h2_DiJetInvMass_vs_MET_background->GetNbinsY();
	
	TH2F* h2_DiJetInvMass_vs_MET_sensitivity;
	h2_DiJetInvMass_vs_MET_sensitivity = new TH2F ("h2_DiJetInvMass_vs_MET_sensitivity","h2_DiJetInvMass_vs_MET_sensitivity", nbinsx, 0., 240., nbinsy , 0., 2500.);		
	h2_DiJetInvMass_vs_MET_sensitivity->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_sensitivity->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_sensitivity->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double sensitivity = getSensitivity(
					getSignalEvents( i , j, h2_DiJetInvMass_vs_MET_eff_signal, xsec, lumi), 
					getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background)
					);
			cout << sensitivity << endl;
			h2_DiJetInvMass_vs_MET_sensitivity->SetBinContent(i, j, sensitivity);
		}
	}

	//paleta(h2_DiJetInvMass_vs_MET_sensitivity);

	TCanvas *my_canvas = new TCanvas;
	my_canvas->cd();
	gPad->SetLogz();
	
	h2_DiJetInvMass_vs_MET_sensitivity->GetZaxis()->SetRangeUser(0.0000001,10000);
	h2_DiJetInvMass_vs_MET_sensitivity->Draw("colz");
	//my_canvas->Print(("JetInvMass_vs_MET_sensitivity_" + isoregion + "_" + taupt + "_LtoT.pdf").c_str());
	my_canvas->Print("JetInvMass_vs_MET_sensitivity_Tau2TightIso_taupt20.pdf");
	my_canvas->Close();
}

