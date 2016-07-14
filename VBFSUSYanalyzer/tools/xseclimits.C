#include <TROOT.h>
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


double vbfefficiency(double evenCRcounts, double oddCRcounts){

	        return ( (oddCRcounts) / ( (oddCRcounts) + (evenCRcounts) ));
}

double vbfConversionFactor(string taupt, string isoregion) {

	TFile *inputfile = TFile::Open(("allQCD_"+ taupt +".root").c_str());
	//TFile *inputfile = TFile::Open("allQCD.root");
	//TFile *inputfile = TFile::Open("allQCD_taupt20.root");
	TH1F* h_ditaucharge;
	TH1F* h_ditauchargeVBFinverted;

	h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	//h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h_ditaucharge").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

	double counts = h_ditaucharge->GetBinContent(3); 
	double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3);


	double vbfeff = vbfefficiency(counts_inverted, counts);

	return (vbfeff / (1. - vbfeff));

}

double LtoTfactor(string taupt) {

	TFile *inputfile = TFile::Open(("allQCD_"+ taupt +"_LtoT.root").c_str());

	TH1F* h1_counts;
	h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

	double oneJetcounts = h1_counts->GetBinContent(2);
	double oneTauMatchcounts = h1_counts->GetBinContent(3);
	double jetAlsoTcounts = h1_counts->GetBinContent(4);
	double fourJetscounts = h1_counts->GetBinContent(5);
	double oneLooseTau = h1_counts->GetBinContent(6);
	double oneLcounts = h1_counts->GetBinContent(7);
	double alsoTcounts = h1_counts->GetBinContent(8);

	//double fourJetscounts_err = binerrors(h1_counts, 2, 2);
	//double twoLooseTaus_err = binerrors(h1_counts, 3, 3);
	//double oneLcounts_err = binerrors(h1_counts, 4, 4);
	//double alsoTcounts_err = binerrors(h1_counts, 5, 5);
	//double jetToTightProb = jetAlsoTcounts / oneJetcounts;
	//double jetToTightProb = jetAlsoTcounts / fourJetscounts;
	double jetToTightProb = jetAlsoTcounts / oneTauMatchcounts;
	double looseToTightProb = alsoTcounts / oneLcounts;
	//double looseToTightProb_err = sqrt( pow( (alsoTcounts_err / oneLcounts )   , 2.) + pow( ( (alsoTcounts * oneLcounts_err) / (oneLcounts * oneLcounts)     )    , 2.) );

	//double twoLooseTo2Tight = looseToTightProb * looseToTightProb;
	//double twoLooseTo2Tight = looseToTightProb * jetToTightProb;
	double twoLooseTo2Tight = jetToTightProb * jetToTightProb;
	//double twoLooseTo2Tight_err = 2. * looseToTightProb * looseToTightProb_err;

	return (twoLooseTo2Tight);	

}


TH2F* makeEffPlot(string taupt, string isoregion) {
	//TFile *inputfile = TFile::Open("VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M.root");
	//TFile *inputfile = TFile::Open("VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau195_Chargino200_1M.root");
	TFile *inputfile = TFile::Open("VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau095_Chargino100_1M.root");
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


	TFile *inputfile = TFile::Open(("allQCD_"+ taupt +".root").c_str());
	//TFile *inputfile = TFile::Open("allQCD.root");
	//TFile *inputfile = TFile::Open("allQCD_taupt20.root");
	TH2F* h2_DiJetInvMass_vs_MET;
	TH2F* h2_DiJetInvMass_vs_MET_LtoT;
	TH1F* h_ditauchargeVBFinverted;
	TH1F* h_count;
	cout << "bla bla" << endl;	
	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h2_DiJetInvMass_vs_MET").c_str())));
	//h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h2_DiJetInvMass_vs_MET").c_str())));
	//h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	h_count = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	//h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h_ditaucharge").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));
	cout << "bla bla bla" << endl;	
	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX(); 
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1);
	//double ntotalnumevents = h2_DiJetInvMass_vs_MET->Integral( 0. , nbinsx, 0. , nbinsy );
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
			//cout << "LtoT: " << ltotfactor << endl;
			//cout << "VBFfactor: " << vbffactor <<endl;
			h2_DiJetInvMass_vs_MET_LtoT->SetBinContent(i,j, (bincontent*vbffactor*ltotfactor));
			//h2_DiJetInvMass_vs_MET_LtoT->SetBinContent(i,j, (bincontent*ltotfactor));

		}
	}

	cout << "cazzo" << endl;
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
	//cout << "minimum: " << minimum << endl;
	//h2_DiJetInvMass_vs_MET->GetZaxis()->SetRangeUser((minimum * 0.00000001),maximum);
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
	///cout << "Total denominator events for " << isoregion << " / "  <<taupt << " : " << ntotalevents << endl;
	//cout << "Total numerator events for " << isoregion << " / "  <<taupt << " : " << ntotalnumevents << endl;
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


void makeSignificance() {

	double xsec = 0.003106; //pb only for test purposes
	double lumi = 85000.; 
	TH2F* h2_DiJetInvMass_vs_MET_eff_signal;
	TH2F* h2_DiJetInvMass_vs_MET_background;
	h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot("taupt20", "Taui2TightIso");
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
			cout << sensitivity << endl;
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
	my_canvas->Print("JetInvMass_vs_MET_significance_Tau2TightIso_taupt20.root");
	my_canvas->Close();
}

void makeXSection() {

	double lumi = 85000.; 
	TH2F* h2_DiJetInvMass_vs_MET_eff_signal;
	TH2F* h2_DiJetInvMass_vs_MET_background;
	//h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot("taupt20", "Taui2TightIso");
	//h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT("taupt20", "Tau2LooseIsoInclusive");
	h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot("taupt20", "Taui2TightIso");
	//h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT("taupt20", "TauAnyIso");
	//h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT("taupt20", "TauAnyIsoPlusNones");
	h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT("taupt20", "baseline");

	cout << "bla bla bla bla" << endl;	
	
	int nbinsx = h2_DiJetInvMass_vs_MET_background->GetNbinsX(); 
	int nbinsy = h2_DiJetInvMass_vs_MET_background->GetNbinsY();
	
	TH1D *p_met = h2_DiJetInvMass_vs_MET_background->ProjectionX();

	
	TH2F* h2_DiJetInvMass_vs_MET_xsec;
	TH2F* h2_DiJetInvMass_vs_MET_xsec_zoom;
	h2_DiJetInvMass_vs_MET_xsec = new TH2F ("h2_DiJetInvMass_vs_MET_xsec","h2_DiJetInvMass_vs_MET_xsec", nbinsx, 0., 240., nbinsy , 0., 2500.);		
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


			cout << xsec << endl;
			h2_DiJetInvMass_vs_MET_xsec->SetBinContent(i, j, xsec);
		}
	}


	cout << "bla bla bla bla bla" << endl;	

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
	my_canvas_zoom->Print("JetInvMass_vs_MET_xsec_Tau2TightIso_taupt20_zoom.root");
	
	TCanvas *my_canvas = new TCanvas("mycanvas","mycanvas",1024.,768.);
	my_canvas->cd();
	gPad->SetLogz();
	
	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetRangeUser(0.001,1000);
	//h2_DiJetInvMass_vs_MET_xsec->Draw("colz text");
	h2_DiJetInvMass_vs_MET_xsec->Draw("colz");
	my_canvas->Print("JetInvMass_vs_MET_xsec_Tau2TightIso_taupt20.root");
	gPad->SetLogy();
	p_met->Draw();
	my_canvas->Print("p_met.root");
	//my_canvas->Print("JetInvMass_vs_MET_xsec_Tau2TightIso_taupt20.pdf");
	my_canvas->Close();
}
