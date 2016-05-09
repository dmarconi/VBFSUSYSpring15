#include <TROOT.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <iostream>
#include <fstream>
#include <cmath>

double vbfefficiency(double evenCRcounts, double oddCRcounts){

	        return ( (oddCRcounts) / ( (oddCRcounts) + (evenCRcounts) ));
}

double vbfConversionFactor(string taupt, string isoregion) {

	TFile *inputfile = TFile::Open(("allQCD_"+ taupt +".root").c_str());
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

void makeEffPlot(string taupt, string isoregion) {
	TFile *inputfile = TFile::Open(("allQCD_"+ taupt +".root").c_str());
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
	//gPad->SetLogy();
	h2_DiJetInvMass_vs_MET_eff->Draw("colz");
	my_canvas->Print(("JetInvMass_vs_MET_eff_" + isoregion + "_" + taupt + ".pdf").c_str());
	my_canvas->Close();
	h2_DiJetInvMass_vs_MET_eff->Clear();
}

void makeEffPlot_LtoT(string taupt, string isoregion){


	TFile *inputfile = TFile::Open(("allQCD_"+ taupt +".root").c_str());
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
	cout << "Total denominator events for " << isoregion << " / "  <<taupt << " : " << ntotalevents << endl;
	cout << "Total numerator events for " << isoregion << " / "  <<taupt << " : " << ntotalnumevents << endl;
}

void mjjvsmeteff() {

	makeEffPlot("taupt20", "Tau2LooseIsoInclusive");	
	makeEffPlot_LtoT("taupt20", "Tau2LooseIsoInclusive");
	makeEffPlot("taupt25", "Tau2LooseIsoInclusive");	
	makeEffPlot_LtoT("taupt25", "Tau2LooseIsoInclusive");
	makeEffPlot("taupt30", "Tau2LooseIsoInclusive");	
	makeEffPlot_LtoT("taupt30", "Tau2LooseIsoInclusive");
	makeEffPlot("taupt35", "Tau2LooseIsoInclusive");	
	makeEffPlot_LtoT("taupt35", "Tau2LooseIsoInclusive");
	makeEffPlot("taupt40", "Tau2LooseIsoInclusive");	
	makeEffPlot_LtoT("taupt40", "Tau2LooseIsoInclusive");
}	
