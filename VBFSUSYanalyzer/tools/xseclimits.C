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
	else if ((chi == "chi100") && (lsp == "lsp050")){ inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau095_Chargino100_1M.root").c_str());}
	else if ((chi == "chi200") && (lsp == "lsp050")){ inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau195_Chargino200_1M.root").c_str());}
	else if ((chi == "chi300") && (lsp == "lsp050")) inputfile = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau295_Chargino300_1M.root").c_str());

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
	h2_DiJetInvMass_vs_MET_xsec->SetTitle("CMS Work in Progress");
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
	makeXsecLimPlots("chi100", "lsp050");
	makeXsecLimPlots("chi200", "lsp050");
	makeXsecLimPlots("chi300", "lsp050");

}

void makeComparisonPlots() {

	Int_t n1 = 13;
	Int_t n2 = 3;
	Double_t x1[n1], y1[n1], x1_err[n1] ,y1_err[n1];
	Double_t x2[n2], y2[n2], x2_err[n2] ,y2_err[n2];

	x2[0] = 100;
	x2[1] = 200;
	x2[2] = 300;

	x2_err[0] = 0.;
	x2_err[1] = 0.;
	x2_err[2] = 0.;

	y2[0] = 0.0327396;
	y2[1] = 0.0333074;
	y2[2] = 0.0345751;

	y2_err[0] = 0.00680598 + sqrt(pow( 0.00368, 2.) + pow( 0.001193, 2.));
	y2_err[1] = 0.0069229+ sqrt(pow(0.00375, 2.) + pow( 0.00121, 2.));
	y2_err[2] =  0.00701739 + sqrt(pow(0.00382, 2.) + pow( 0.00123, 2.));

	x1[0] = 100;
	x1[1] = 125;
	x1[2] = 150;
	x1[3] = 175;
	x1[4] = 200;
	x1[5] = 225;
	x1[6] = 250;
	x1[7] = 275;
	x1[8] = 300;
	x1[9] = 325;
	x1[10] = 350;
	x1[11] = 375;
	x1[12] = 400;

	x1_err[0] = 0.;
	x1_err[1] = 0.;
	x1_err[2] = 0.;
	x1_err[3] = 0.;
	x1_err[4] = 0.;
	x1_err[5] = 0.;
	x1_err[6] = 0.;
	x1_err[7] = 0.;
	x1_err[8] = 0.;
	x1_err[9] = 0.;
	x1_err[10] = 0.;
	x1_err[11] = 0.;
	x1_err[12] = 0.;


	y1[0] = 22670.1;
	y1[1] = 10034.8;
	y1[2] = 5180.86;
	y1[3] = 2953.28;
	y1[4] = 1807.39;
	y1[5] = 1165.09;
	y1[6] = 782.487;
	y1[7] = 543.03;
	y1[8] = 386.936;
	y1[9] = 281.911;
	y1[10] = 209.439;
	y1[11] = 158.06;
	y1[12] = 121.013;

	y1_err[0] = 973.967;
	y1_err[1] = 457.604;
	y1_err[2] = 253.223;
	y1_err[3] = 154.386;
	y1_err[4] = 101.316;
	y1_err[5] = 68.8042;
	y1_err[6] = 48.7463;
	y1_err[7] = 35.4083;
	y1_err[8] = 26.3602;
	y1_err[9] = 20.0201;
	y1_err[10] = 15.4539;
	y1_err[11] = 12.0956;
	y1_err[12] = 9.61659;

	for(Int_t i = 0; i < n1; i++) y1[i] *= 0.001;
	for(Int_t i = 0; i < n1; i++) y1_err[i] *= 0.001;

	double y_max = y1[0] * 1.5;
	double y_min = y2[2] * 0.5;
	double x_max = 450.;
	double x_min = 0.;

	TGraph *gr1 = new TGraphErrors (n1, x1, y1, x1_err, y1_err);
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(21);

	TGraph *gr2 = new TGraphErrors (n2, x2, y2, x2_err, y2_err);
	gr2->SetMarkerColor(6);
	gr2->SetMarkerStyle(21);

	TH1F* lim_comparison_bkg = new TH1F ("lim_comparison","lim_comparison", 9, 0. , 450.);
	lim_comparison_bkg->SetTitle("CMS Work in Progress");
	lim_comparison_bkg->GetYaxis()->SetTitle("#sigma [pb]");
	lim_comparison_bkg->GetXaxis()->SetTitle("m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) [GeV]");
	lim_comparison_bkg->GetYaxis()->SetTitleOffset(1.40);
	lim_comparison_bkg->GetYaxis()->SetRangeUser(y_min,y_max);
	lim_comparison_bkg->GetXaxis()->SetRangeUser(0.,450.);
	lim_comparison_bkg->SetStats(0);

	TH1F* lim_comparison = new TH1F ("lim_comparison","lim_comparison", 9, x_min - 25. , x_max - 25.);
	lim_comparison->SetStats(0);
	lim_comparison->SetMarkerStyle(21);
	lim_comparison->SetMarkerColor(6);
	lim_comparison->SetBinContent( 3, y2[0]);
	lim_comparison->SetBinContent( 5, y2[1]);
	lim_comparison->SetBinContent( 7, y2[2]);
	lim_comparison->SetBinError( 3, y2_err[0]);
	lim_comparison->SetBinError( 5, y2_err[1]);
	lim_comparison->SetBinError( 7, y2_err[2]);


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

	TCanvas *my_canvas = new TCanvas("mycanvas","mycanvas",600.,600.);
	my_canvas->cd();
	gPad->SetLogy();

	lim_comparison_bkg->Draw();
	lim_comparison->Draw("Pe+same");
  mg->Draw("LP");
	leg->Draw("same");
	my_canvas->Update();
	my_canvas->Print("results/xsec_confront.pdf");

}
