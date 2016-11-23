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
#include <TH2F.h>
#include <TH2D.h>
#include <TProfile.h>
#include <iostream>
#include <fstream>
#include <cmath>

double vbfefficiency(double evenCRcounts, double oddCRcounts){

	return ( (oddCRcounts) / ( (oddCRcounts) + (evenCRcounts) ));
}

double vbfefficiency_statunc(double evenCRcounts,double evenCRcounts_statunc,
                              double oddCRcounts, double oddCRcounts_statunc){

  return (  (  (evenCRcounts * oddCRcounts_statunc)  +  (oddCRcounts * evenCRcounts_statunc)  )  /
            (  (evenCRcounts + oddCRcounts ) * (evenCRcounts + oddCRcounts )  )  );
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

double vbfConversionFactor_statunc(string taupt, string isoregion) {

	TFile *inputfile = TFile::Open((taupt + "/allQCD_"+ taupt +".root").c_str());
	TH1F* h_ditaucharge;
	TH1F* h_ditauchargeVBFinverted;

	h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	//h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h_ditaucharge").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

  double counts = h_ditaucharge->GetBinContent(3);
  double counts_statunc = h_ditaucharge->GetBinError(3);
  double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3);
  double counts_inverted_statunc = h_ditauchargeVBFinverted->GetBinError(3);


	double vbfeff = vbfefficiency(counts_inverted, counts);
  double vbfeff_statunc = vbfefficiency_statunc(counts_inverted,counts_inverted_statunc, counts, counts_statunc);

	delete inputfile;
	return ((2. * vbfeff * vbfeff_statunc) / ( (1. - vbfeff) * (1. - vbfeff) ) );


}


double LtoTfactor(string taupt) {

	TFile *inputfile = TFile::Open(( taupt + "/allQCD_"+ taupt +"_LtoT.root").c_str());

	TH1F* h1_counts;
	h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

	double oneTauMatchcounts = h1_counts->GetBinContent(3);
	double jetAlsoTcounts = h1_counts->GetBinContent(4);

	double jetToTightProb = jetAlsoTcounts / oneTauMatchcounts;
	double twoLooseTo2Tight = jetToTightProb * jetToTightProb;

	delete inputfile;
	return (twoLooseTo2Tight);

}

double LtoTfactor_statunc(string taupt) {

	TFile *inputfile = TFile::Open(( taupt + "/allQCD_"+ taupt +"_LtoT.root").c_str());

	TH1F* h1_counts;
	h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

  double oneTauMatchcounts = h1_counts->GetBinContent(3);
  double oneTauMatchcounts_statunc = h1_counts->GetBinError(3);
  double jetAlsoTcounts = h1_counts->GetBinContent(4);
  double jetAlsoTcounts_statunc = h1_counts->GetBinError(4);

  double jetToTightProb = jetAlsoTcounts / oneTauMatchcounts;
	double jetToTightProb_statunc = (jetAlsoTcounts / oneTauMatchcounts) *
                                  sqrt( ( (oneTauMatchcounts_statunc * oneTauMatchcounts_statunc) / (oneTauMatchcounts * oneTauMatchcounts) ) +
                                  ( (jetAlsoTcounts_statunc * jetAlsoTcounts_statunc) / (jetAlsoTcounts * jetAlsoTcounts) ) );
  double twoLooseTo2Tight = jetToTightProb * jetToTightProb;
  double twoLooseTo2Tight_statunc = 2. * jetToTightProb * jetToTightProb_statunc;

	delete inputfile;
	return (twoLooseTo2Tight_statunc);

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

TH2F* makeEffPlotStatUnc(string taupt, string isoregion, string chi, string lsp) {
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
  double ntotalevents_statunc = h_count->GetBinError(1);
	//cout << "N total events: " << ntotalevents << endl;
	//double ntotalevents = h_ditaucharge->GetBinContent(3);

	TH2F* h2_DiJetInvMass_vs_MET_eff_statunc;
	h2_DiJetInvMass_vs_MET_eff_statunc = new TH2F ("h2_DiJetInvMass_vs_MET_eff_statunc","h2_DiJetInvMass_vs_MET_eff_statunc", nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_eff_statunc->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_eff_statunc->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_eff_statunc->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
      double integral_statunc = 0.;
      double integral = h2_DiJetInvMass_vs_MET->IntegralAndError( i, nbinsx, j, nbinsy, integral_statunc );
       h2_DiJetInvMass_vs_MET_eff_statunc->SetBinContent(i,j,
        //error propagation
        ( (integral/ntotalevents) *
        sqrt( (integral_statunc / integral) * (integral_statunc / integral)
            + (ntotalevents_statunc / ntotalevents) * (ntotalevents_statunc / ntotalevents) ) )
      );
		}
	}
	return h2_DiJetInvMass_vs_MET_eff_statunc;
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

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double bincontent = h2_DiJetInvMass_vs_MET->GetBinContent (i, j);
			double vbffactor = vbfConversionFactor(taupt, isoregion);
			double ltotfactor = LtoTfactor(taupt);
			h2_DiJetInvMass_vs_MET_LtoT->SetBinContent(i,j, (bincontent*vbffactor*ltotfactor));

		}
	}

	return h2_DiJetInvMass_vs_MET_LtoT;
}

TH2F* makeBackgroundPlot_LtoT_StatUnc(string taupt, string isoregion){


	TFile *inputfile = TFile::Open((taupt + "/allQCD_"+ taupt +".root").c_str());

	TH2F* h2_DiJetInvMass_vs_MET;
	TH2F* h2_DiJetInvMass_vs_MET_LtoT_statunc;
	TH1F* h_ditauchargeVBFinverted;
	TH1F* h_count;

	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h2_DiJetInvMass_vs_MET").c_str())));
	h_count = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();

	h2_DiJetInvMass_vs_MET_LtoT_statunc = new TH2F ("h2_DiJetInvMass_vs_MET_LtoT_statunc","h2_DiJetInvMass_vs_MET_LtoT_statunc", nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_LtoT_statunc->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT_statunc->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT_statunc->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
      double bincontent = h2_DiJetInvMass_vs_MET->GetBinContent (i, j);
      double bincontent_statunc = h2_DiJetInvMass_vs_MET->GetBinError (i, j);
      double vbffactor = vbfConversionFactor(taupt, isoregion);
      double vbffactor_statunc = vbfConversionFactor_statunc(taupt, isoregion);
      double ltotfactor = LtoTfactor(taupt);
      double ltotfactor_statunc = LtoTfactor(taupt);

			h2_DiJetInvMass_vs_MET_LtoT_statunc->SetBinContent(i,j,
        sqrt(pow( (bincontent_statunc*vbffactor*ltotfactor) , 2.) + pow( (bincontent*vbffactor_statunc*ltotfactor) , 2.) + pow( (bincontent*vbffactor*ltotfactor_statunc) , 2.) )
      );
		}
	}

	return h2_DiJetInvMass_vs_MET_LtoT_statunc;
}

double getSignalEfficiency(int xbin, int ybin, TH2F* signal_map){
	double efficiency = signal_map->GetBinContent(xbin,ybin);
	return (efficiency);
}

double getSignalEfficiencyStatUnc(int xbin, int ybin, TH2F* signal_map){
	double efficiency = signal_map->GetBinError(xbin,ybin);
	return (efficiency);
}

double getBackgroudEvents(int xbin, int ybin, TH2F* background_map){
	int nbinsx = background_map->GetNbinsX();
	int nbinsy = background_map->GetNbinsY();
	double events = background_map->Integral( xbin, nbinsx, ybin, nbinsy );
	return (events);
}

double getBackgroudEventsStatUnc(int xbin, int ybin, TH2F* background_map_unc){
	int nbinsx = background_map_unc->GetNbinsX();
	int nbinsy = background_map_unc->GetNbinsY();
  double squaredsum = 0.;

  for (int i = xbin; i < nbinsx; i++) {
		for (int j = ybin; j < nbinsy; j++) {
      squaredsum += (background_map_unc->GetBinContent(i,j)) * (background_map_unc->GetBinContent(i,j));
    }
  }

	return (sqrt(squaredsum));
}

double getXSection(double efficiency, double luminosity, double background, double sigma, double significance) {
	return ( (significance * ( (0.5 * sigma) + sqrt( background + (0.5 * background * background))  ) ) / (efficiency * luminosity)   );
}

double getXSectionStatUnc(double efficiency, double efficiency_statunc, double luminosity, double background, double background_statunc, double sigma, double significance) {

  double var_a = sqrt(background);
  double var_b = sqrt( 0.5 * background + 1.);

  double bkg_part = sqrt( background + (0.5 * background * background));
  double bkg_part_statunc = sqrt(0.5 * ( (var_b / var_a) + 0.5 * ( var_a / var_b ) ) ) * background_statunc;

  double var_c = significance / luminosity;
  double var_d = 0.5 * sigma;

  double xsec = var_c * (var_d + bkg_part) / efficiency;
  double xsec_statunc = xsec * sqrt( pow((bkg_part_statunc/bkg_part),2.) + pow((efficiency_statunc/efficiency),2.) );

  return (  xsec_statunc  );
}

void makeXSection(string taupt,string chi, string lsp) {

  //create eff plot

  //create bg plot

  double lumi = 85000.;
  TH2F* h2_DiJetInvMass_vs_MET_eff_signal;
  TH2F* h2_DiJetInvMass_vs_MET_eff_signal_stat;
  TH2F* h2_DiJetInvMass_vs_MET_background;
  TH2F* h2_DiJetInvMass_vs_MET_background_stat;

  h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot(taupt, "Taui2TightIso", chi, lsp);
  h2_DiJetInvMass_vs_MET_eff_signal_stat = makeEffPlotStatUnc(taupt, "Taui2TightIso", chi, lsp);
  h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT(taupt, "baseline");
  h2_DiJetInvMass_vs_MET_background_stat = makeBackgroundPlot_LtoT_StatUnc(taupt, "baseline");

  int nbinsx = h2_DiJetInvMass_vs_MET_background->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET_background->GetNbinsY();

  TH2F* h2_DiJetInvMass_vs_MET_xsec;
  h2_DiJetInvMass_vs_MET_xsec = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(),("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(), nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec->SetTitle("CMS Work in Progress");
	h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->SetTitle("M_{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetTitle("#sigma_{lim} pb");
	h2_DiJetInvMass_vs_MET_xsec->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->SetTitleOffset(1.40);
	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetTitleOffset(1.13);
	h2_DiJetInvMass_vs_MET_xsec->SetStats(0);

  TH2F* h2_DiJetInvMass_vs_MET_xsec_stat;
  h2_DiJetInvMass_vs_MET_xsec_stat = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(),("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(), nbinsx, 0., 240., nbinsy , 0., 2500.);


  for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double xsec = getXSection(
					getSignalEfficiency(i,j, h2_DiJetInvMass_vs_MET_eff_signal),
					lumi,
					getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background),
					2.5, //sigma
					2.0 //significance
					);

      double xsec_statunc = getXSectionStatUnc(
        getSignalEfficiency(i,j, h2_DiJetInvMass_vs_MET_eff_signal),
        getSignalEfficiencyStatUnc(i,j, h2_DiJetInvMass_vs_MET_eff_signal_stat),
        lumi,
        getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background),
        getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background_stat),
        2.5, //sigma
        2.0 //significance
        );

        h2_DiJetInvMass_vs_MET_xsec->SetBinContent(i, j, xsec);
        h2_DiJetInvMass_vs_MET_xsec_stat->SetBinContent(i, j, xsec_statunc);

        cout << "xsec: " << xsec << " unc: " << xsec_statunc << endl;

		}
	}


}

void makeXsecLimPlots(string chi, string lsp){

	makeXSection("taupt20",chi,lsp);

}

void fullXsecLimScan(){

	makeXsecLimPlots("chi100", "lsp000");

}
