#include <TROOT.h>
#include <TPad.h>
#include <TFrame.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TColor.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <array>


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
	if ((chi == "chi100") && (lsp == "lsp000")){ inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau095_Chargino100_1M.root").c_str());}
	else if ((chi == "chi200") && (lsp == "lsp000")){ inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau195_Chargino200_1M.root").c_str());}
	else if ((chi == "chi300") && (lsp == "lsp000")){ inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M.root").c_str());}
	else if ((chi == "chi400") && (lsp == "lsp000")){ inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau395_Chargino400_1M.root").c_str());}
	else if ((chi == "chi500") && (lsp == "lsp000")){ inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau495_Chargino500_1M.root").c_str());}
	else if ((chi == "chi100") && (lsp == "lsp050")){ inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau095_Chargino100_1M.root").c_str());}
	else if ((chi == "chi200") && (lsp == "lsp050")){ inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau195_Chargino200_1M.root").c_str());}
	else if ((chi == "chi300") && (lsp == "lsp050")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau295_Chargino300_1M.root").c_str());
	else if ((chi == "chi400") && (lsp == "lsp050")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau395_Chargino400_1M.root").c_str());
	else if ((chi == "chi500") && (lsp == "lsp050")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau495_Chargino500_1M.root").c_str());

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

	TH2F* h2_DiJetInvMass_vs_MET;
	TH2F* h2_DiJetInvMass_vs_MET_LtoT;
	TH1F* h_ditauchargeVBFinverted;
	TH1F* h_count;

	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h2_DiJetInvMass_vs_MET").c_str())));
	h_count = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1);

	h2_DiJetInvMass_vs_MET_LtoT = new TH2F ("h2_DiJetInvMass_vs_MET_LtoT","h2_DiJetInvMass_vs_MET_LtoT", nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_LtoT->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->SetStats(0);

	TH2F* h2_DiJetInvMass_vs_MET_eff;
	h2_DiJetInvMass_vs_MET_eff = new TH2F ("h2_DiJetInvMass_vs_MET_eff","h2_DiJetInvMass_vs_MET_eff", nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_eff->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_eff->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_eff->SetStats(0);

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
	h2_DiJetInvMass_vs_MET_xsec->SetTitle("CMS Work");
	h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->SetTitle("M_{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetTitle("#sigma_{lim} [pb]");
	h2_DiJetInvMass_vs_MET_xsec->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->SetTitleOffset(1.40);
	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetTitleOffset(1.13);
	h2_DiJetInvMass_vs_MET_xsec->SetStats(0);

	//Converting model point into doubles
	string chi_value;
	string lsp_value;
	string taupt_label;

	if (chi == "chi100") {chi_value = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 100 GeV";}
	else if (chi == "chi200") {chi_value = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 200 GeV";}
	else if (chi == "chi300") {chi_value = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 300 GeV";}
	else if (chi == "chi400") {chi_value = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 400 GeV";}
	else if (chi == "chi500") {chi_value = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 500 GeV";}

	if (lsp == "lsp000") {lsp_value = "m(#tilde{#chi}^{0}_{1}) = 0 GeV";}
	else if (lsp == "lsp050") {lsp_value = "m(#tilde{#chi}^{0}_{1}) = 50 GeV";}

	if (taupt == "taupt20") {taupt_label = "#tau_{Pt} = 20 GeV";}
	else if (taupt == "taupt25") {taupt_label = "#tau_{pt} = 25 GeV";}
	else if (taupt == "taupt30") {taupt_label = "#tau_{pt} = 30 GeV";}
	else if (taupt == "taupt35") {taupt_label = "#tau_{pt} = 35 GeV";}
	else if (taupt == "taupt40") {taupt_label = "#tau_{pt} = 40 GeV";}
	else if (taupt == "taupt45") {taupt_label = "#tau_{pt} = 45 GeV";}

	//if ((chi == "chi100") && (lsp == "lsp000")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau095_Chargino100_1M.root").c_str());

	//TH2F* h2_DiJetInvMass_vs_MET;


	//defining legend
	TLegend* leg = new TLegend(0.63,0.68,0.86,0.9);
	leg->SetTextSize(0.02);
	leg->AddEntry((TObject*)0, (chi_value).c_str(), "");
	leg->AddEntry((TObject*)0, "m(#tilde{#chi}^{#pm}_{1}) - m(#tilde{#tau}_{1}) = 5 GeV", "");
	leg->AddEntry((TObject*)0, (lsp_value).c_str(), "");
	leg->AddEntry((TObject*)0, (taupt_label).c_str(), "");

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
	h2_DiJetInvMass_vs_MET_xsec_zoom->GetXaxis()->SetRangeUser(30.,200.);
	h2_DiJetInvMass_vs_MET_xsec_zoom->GetYaxis()->SetRangeUser(0.,1000.);
	double maximum_zoom = h2_DiJetInvMass_vs_MET_xsec_zoom->GetMaximum();
	double minimum_zoom = h2_DiJetInvMass_vs_MET_xsec_zoom->GetMinimum();
	h2_DiJetInvMass_vs_MET_xsec_zoom ->GetZaxis()->SetRangeUser(minimum_zoom,maximum_zoom);
	h2_DiJetInvMass_vs_MET_xsec_zoom->Draw("colz");
	leg->Draw();
	gPad->SetRightMargin(0.14);
	gPad->SetLogz();
	gPad->Update();
	my_canvas_zoom->Print(("results/JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt + "_zoom.pdf").c_str());

	TCanvas *my_canvas = new TCanvas("mycanvas","mycanvas",1024.,768.);
	my_canvas->cd();

	int binmin = h2_DiJetInvMass_vs_MET_xsec->GetMinimumBin();
	double minimum = FLT_MAX;
	double x_min = 0.;
	double y_min = 0.;
	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			if ((h2_DiJetInvMass_vs_MET_xsec->GetBinContent(i,j) < minimum) && (h2_DiJetInvMass_vs_MET_xsec->GetBinContent(i,j) > 0)){
				minimum = h2_DiJetInvMass_vs_MET_xsec->GetBinContent(i,j);
				x_min = h2_DiJetInvMass_vs_MET_xsec->GetXaxis()->GetBinCenter(i);
				x_min = x_min - 0.5 * h2_DiJetInvMass_vs_MET_xsec->GetXaxis()->GetBinWidth(i);
				y_min = h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->GetBinCenter(j);
				y_min = y_min - 0.5 * h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->GetBinWidth(i);
			}
		}
	}

	double taupt_value = 0.;
	if (taupt == "taupt20") {taupt_value = 20.;} else
	if (taupt == "taupt25") {taupt_value = 25.;} else
	if (taupt == "taupt30") {taupt_value = 30.;} else
	if (taupt == "taupt35") {taupt_value = 35.;} else
	if (taupt == "taupt40") {taupt_value = 40.;} else
	if (taupt == "taupt45") taupt_value = 45.;

	//cout << "Minimum at " << minimum << " for " << chi << ", " << lsp << " and " << taupt << ": Jetinvmass " << y_min << " MET " << x_min << endl;
	cout << minimum << " & $<$ " << taupt_value << " & $<$ " << y_min << "  & $<$ " << x_min << " \\\\ " << endl;

	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetRangeUser(0.001,1000);
	h2_DiJetInvMass_vs_MET_xsec->Draw("colz");
	leg->Draw();

	gPad->SetRightMargin(0.14);
	gPad->SetLogz();
	gPad->Update();

	my_canvas->SetRightMargin(0.14);
	my_canvas->Print(("results/JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt + ".pdf").c_str());
	//my_canvas->Print(("results/JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt + ".root").c_str());
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
	makeXsecLimPlots("chi400", "lsp000");
	makeXsecLimPlots("chi500", "lsp000");
	makeXsecLimPlots("chi100", "lsp050");
	makeXsecLimPlots("chi200", "lsp050");
	makeXsecLimPlots("chi300", "lsp050");
	makeXsecLimPlots("chi400", "lsp050");
	makeXsecLimPlots("chi500", "lsp050");

}

void makeComparisonPlots() {

	//Arrays declarations

	Int_t n1 = 13;
	Int_t n2 = 5;
	Double_t cmsxsec_x[n1], cmsxsec_y[n1], cmsxsec_x_err[n1] ,cmsxsec_y_err[n1];
	Double_t vbfxsec_x[n2], vbfxsec_y[n2], vbfxsec_x_err[n2] ,vbfxsec_y_err[n2];

	//Filling arrays with VBF analysis and CMS cross section limits

	vbfxsec_x[0] = 100;
	vbfxsec_x[1] = 200;
	vbfxsec_x[2] = 300;
	vbfxsec_x[3] = 400;
	vbfxsec_x[4] = 500;

	vbfxsec_x_err[0] = 0.;
	vbfxsec_x_err[1] = 0.;
	vbfxsec_x_err[2] = 0.;
	vbfxsec_x_err[3] = 0.;
	vbfxsec_x_err[4] = 0.;

	vbfxsec_y[0] = 0.0327396;
	vbfxsec_y[1] = 0.0333074;
	vbfxsec_y[2] = 0.0345751;
	vbfxsec_y[3] = 0.0302742;
	vbfxsec_y[4] = 0.0345751;

	vbfxsec_y_err[0] = 0.00680598 + sqrt(pow( 0.00368, 2.) + pow( 0.001193, 2.));
	vbfxsec_y_err[1] = 0.0069229+ sqrt(pow(0.00375, 2.) + pow( 0.00121, 2.));
	vbfxsec_y_err[2] =  0.00701739 + sqrt(pow(0.00382, 2.) + pow( 0.00123, 2.));
	vbfxsec_y_err[3] =  0.00458659 + sqrt(pow(0.00271169, 2.) + pow( 0.000385888, 2.));
	vbfxsec_y_err[4] =  0.00631726 + sqrt(pow(0.00342821, 2.) + pow( 0.000495253, 2.));


	cmsxsec_x[0] = 100;
	cmsxsec_x[1] = 125;
	cmsxsec_x[2] = 150;
	cmsxsec_x[3] = 175;
	cmsxsec_x[4] = 200;
	cmsxsec_x[5] = 225;
	cmsxsec_x[6] = 250;
	cmsxsec_x[7] = 275;
	cmsxsec_x[8] = 300;
	cmsxsec_x[9] = 325;
	cmsxsec_x[10] = 350;
	cmsxsec_x[11] = 375;
	cmsxsec_x[12] = 400;

	cmsxsec_x_err[0] = 0.;
	cmsxsec_x_err[1] = 0.;
	cmsxsec_x_err[2] = 0.;
	cmsxsec_x_err[3] = 0.;
	cmsxsec_x_err[4] = 0.;
	cmsxsec_x_err[5] = 0.;
	cmsxsec_x_err[6] = 0.;
	cmsxsec_x_err[7] = 0.;
	cmsxsec_x_err[8] = 0.;
	cmsxsec_x_err[9] = 0.;
	cmsxsec_x_err[10] = 0.;
	cmsxsec_x_err[11] = 0.;
	cmsxsec_x_err[12] = 0.;


	cmsxsec_y[0] = 22670.1;
	cmsxsec_y[1] = 10034.8;
	cmsxsec_y[2] = 5180.86;
	cmsxsec_y[3] = 2953.28;
	cmsxsec_y[4] = 1807.39;
	cmsxsec_y[5] = 1165.09;
	cmsxsec_y[6] = 782.487;
	cmsxsec_y[7] = 543.03;
	cmsxsec_y[8] = 386.936;
	cmsxsec_y[9] = 281.911;
	cmsxsec_y[10] = 209.439;
	cmsxsec_y[11] = 158.06;
	cmsxsec_y[12] = 121.013;

	cmsxsec_y_err[0] = 973.967;
	cmsxsec_y_err[1] = 457.604;
	cmsxsec_y_err[2] = 253.223;
	cmsxsec_y_err[3] = 154.386;
	cmsxsec_y_err[4] = 101.316;
	cmsxsec_y_err[5] = 68.8042;
	cmsxsec_y_err[6] = 48.7463;
	cmsxsec_y_err[7] = 35.4083;
	cmsxsec_y_err[8] = 26.3602;
	cmsxsec_y_err[9] = 20.0201;
	cmsxsec_y_err[10] = 15.4539;
	cmsxsec_y_err[11] = 12.0956;
	cmsxsec_y_err[12] = 9.61659;

	for(Int_t i = 0; i < n1; i++) cmsxsec_y[i] *= 0.001;
	for(Int_t i = 0; i < n1; i++) cmsxsec_y_err[i] *= 0.001;


	//Set plot minimum
	double y_max = cmsxsec_y[0] * 1.5;
	double y_min = vbfxsec_y[2] * 0.5;
	double x_max = 550.;
	double x_min = 0.;

	//TGraph definition
	TGraph *gr1 = new TGraphErrors (n1, cmsxsec_x, cmsxsec_y, cmsxsec_x_err, cmsxsec_y_err);
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(21);

	TGraph *gr2 = new TGraphErrors (n2, vbfxsec_x, vbfxsec_y, vbfxsec_x_err, vbfxsec_y_err);
	gr2->SetMarkerColor(6);
	gr2->SetMarkerStyle(21);

	//TH1F definition
	TH1F* lim_comparison_bkg = new TH1F ("lim_comparison","lim_comparison", 13, 0. , 550.);
	lim_comparison_bkg->SetTitle("CMS Work");
	lim_comparison_bkg->GetYaxis()->SetTitle("#sigma [pb]");
	lim_comparison_bkg->GetXaxis()->SetTitle("m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) [GeV]");
	lim_comparison_bkg->GetYaxis()->SetTitleOffset(1.40);
	lim_comparison_bkg->GetYaxis()->SetRangeUser(y_min,y_max);
	lim_comparison_bkg->GetXaxis()->SetRangeUser(0.,450.);
	lim_comparison_bkg->SetStats(0);

	TH1F* lim_comparison = new TH1F ("lim_comparison","lim_comparison", 13, x_min - 25. , x_max - 25.);
	lim_comparison->SetStats(0);
	lim_comparison->SetMarkerStyle(21);
	lim_comparison->SetMarkerColor(6);
	lim_comparison->SetBinContent( 3, vbfxsec_y[0]);
	lim_comparison->SetBinContent( 5, vbfxsec_y[1]);
	lim_comparison->SetBinContent( 7, vbfxsec_y[2]);
	lim_comparison->SetBinContent( 9, vbfxsec_y[3]);
	lim_comparison->SetBinContent( 11, vbfxsec_y[4]);
	lim_comparison->SetBinError( 3, vbfxsec_y_err[0]);
	lim_comparison->SetBinError( 5, vbfxsec_y_err[1]);
	lim_comparison->SetBinError( 7, vbfxsec_y_err[2]);
	lim_comparison->SetBinError( 9, vbfxsec_y_err[3]);
	lim_comparison->SetBinError( 11, vbfxsec_y_err[4]);


	// create a multigraph and draw it
	TMultiGraph  *mg  = new TMultiGraph();
	mg->Add(gr1);
  //mg->Add(gr2);
	//mg->GetXaxis()->SetTitle("E_{#gamma} (GeV)");
  //mg->GetYaxis()->SetTitle("Coefficients");

	//defining legend
	TLegend* leg = new TLegend(0.63,0.68,0.86,0.9);
	leg->SetTextSize(0.04);
	leg->AddEntry(gr1, "#sigma^{CMS} [pb]", "LP");
	leg->AddEntry(lim_comparison, "#sigma_{lim}^{VBF} [pb]", "Pe");

	//Canvas definition
	TCanvas *my_canvas = new TCanvas("mycanvas","mycanvas",600.,600.);
	my_canvas->cd();
	gPad->SetLogy();

	//Drawing workaround
	lim_comparison_bkg->Draw();
	lim_comparison->Draw("Pe+same");
  mg->Draw("LP");
	leg->Draw("same");
	my_canvas->Update();
	my_canvas->Print("results/xsec_confront.pdf");

}
