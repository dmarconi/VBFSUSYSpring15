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

double vbfConversionFactorMCSyst(string taupt, string isoregion, double variation) {

	TFile *inputfile = TFile::Open((taupt + "/allQCD_"+ taupt +".root").c_str());
	TH1F* h_ditaucharge;
	TH1F* h_ditauchargeVBFinverted;

	h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	//h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h_ditaucharge").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

	double counts = h_ditaucharge->GetBinContent(3) + (variation * h_ditaucharge->GetBinContent(3));
	double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3) + (variation * h_ditauchargeVBFinverted->GetBinContent(3));


	double vbfeff = vbfefficiency(counts_inverted, counts);

	delete inputfile;
	return (vbfeff / (1. - vbfeff));


}

double vbfConversionFactoVBFSyst(string taupt, string isoregion, double variation) {

	TFile *inputfile = TFile::Open((taupt + "/allQCD_"+ taupt +".root").c_str());
	TH1F* h_ditaucharge;
	TH1F* h_ditauchargeVBFinverted;

	h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
	//h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedObjectSelection/h_ditaucharge").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

	double counts = h_ditaucharge->GetBinContent(3);
	double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3);


	double vbfeff = (1. + variation) * vbfefficiency(counts_inverted, counts);

	delete inputfile;
	return (vbfeff / (1. - vbfeff));


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

double LtoTfactorSyst(string taupt, double variation) {

	TFile *inputfile = TFile::Open(( taupt + "/allQCD_"+ taupt +"_LtoT.root").c_str());

	TH1F* h1_counts;
	h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

	double oneTauMatchcounts = h1_counts->GetBinContent(3) + (variation * h1_counts->GetBinContent(3));
	double jetAlsoTcounts = h1_counts->GetBinContent(4) + (variation * h1_counts->GetBinContent(4));

	double jetToTightProb = jetAlsoTcounts / oneTauMatchcounts;
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
			double eff = integral / ntotalevents;
			double eff_statunc = eff * sqrt (pow((ntotalevents_statunc/ntotalevents),2.) + pow((integral_statunc/integral),2.));
      h2_DiJetInvMass_vs_MET_eff_statunc->SetBinContent(i,j,eff_statunc);
		}
	}
	return h2_DiJetInvMass_vs_MET_eff_statunc;
}

TH2F* makeEffPlotSyst(string taupt, string isoregion, string chi, string lsp, double variation) {
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
	double ntotalevents = h_count->GetBinContent(1) + (variation * h_count->GetBinContent(1));
	//double ntotalevents = h_ditaucharge->GetBinContent(3);

	TH2F* h2_DiJetInvMass_vs_MET_eff;
	h2_DiJetInvMass_vs_MET_eff = new TH2F ("h2_DiJetInvMass_vs_MET_eff","h2_DiJetInvMass_vs_MET_eff", nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_eff->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_eff->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_eff->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double integral = h2_DiJetInvMass_vs_MET->Integral( i, nbinsx, j, nbinsy ) + (variation * h2_DiJetInvMass_vs_MET->Integral( i, nbinsx, j, nbinsy ));
			h2_DiJetInvMass_vs_MET_eff->SetBinContent(i,j, (integral/ntotalevents));
		}
	}
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

TH2F* makeBackgroundPlot_LtoT_MCSyst(string taupt, string isoregion, double variation){


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
	double ntotalevents = h_count->GetBinContent(1) + (variation * h_count->GetBinContent(1));

	h2_DiJetInvMass_vs_MET_LtoT = new TH2F ("h2_DiJetInvMass_vs_MET_LtoT","h2_DiJetInvMass_vs_MET_LtoT", nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_LtoT->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double bincontent = h2_DiJetInvMass_vs_MET->GetBinContent(i,j) + (variation * h2_DiJetInvMass_vs_MET->GetBinContent(i,j));
			double vbffactor = vbfConversionFactorMCSyst(taupt, isoregion, variation);
			double ltotfactor = LtoTfactorSyst(taupt,variation);
			h2_DiJetInvMass_vs_MET_LtoT->SetBinContent(i,j, (bincontent*vbffactor*ltotfactor));

		}
	}

	return h2_DiJetInvMass_vs_MET_LtoT;
}

TH2F* makeBackgroundPlot_LtoT_VBFSyst(string taupt, string isoregion, double variation){


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
	double ntotalevents = h_count->GetBinContent(1) + (variation * h_count->GetBinContent(1));

	h2_DiJetInvMass_vs_MET_LtoT = new TH2F ("h2_DiJetInvMass_vs_MET_LtoT","h2_DiJetInvMass_vs_MET_LtoT", nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_LtoT->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double bincontent = h2_DiJetInvMass_vs_MET->GetBinContent(i,j);
			double vbffactor = vbfConversionFactoVBFSyst(taupt, isoregion, variation);
			double ltotfactor = LtoTfactor(taupt);
			h2_DiJetInvMass_vs_MET_LtoT->SetBinContent(i,j, (bincontent*vbffactor*ltotfactor));

		}
	}

	return h2_DiJetInvMass_vs_MET_LtoT;
}


double getSignalEfficiency(int xbin, int ybin, TH2F* signal_map){
	double efficiency = signal_map->GetBinContent(xbin,ybin);
	return (efficiency);
}

double getSignalEfficiencyStatUnc(int xbin, int ybin, TH2F* signal_map){
	double efficiency = signal_map->GetBinContent(xbin,ybin);
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

  double var_c = significance * ((0.5 * sigma) + bkg_part );
	double var_c_statunc = significance * bkg_part_statunc;
  double var_d = efficiency * luminosity;
	double var_d_statunc = luminosity * efficiency_statunc;

  double xsec = var_c / var_d;
	//double xsec_statunc = xsec * sqrt( pow((bkg_part_statunc/bkg_part),2.) + pow((efficiency_statunc/efficiency),2.) );
	double xsec_statunc = xsec * sqrt( pow((var_c_statunc/var_c),2.) + pow((var_d_statunc/var_d),2.) );

  return (  xsec_statunc  );
}

void makeXSection(string taupt,string chi, string lsp) {

  double lumi = 85000.;
  TH2F* h2_DiJetInvMass_vs_MET_eff_signal;
	TH2F* h2_DiJetInvMass_vs_MET_eff_signal_stat;
	TH2F* h2_DiJetInvMass_vs_MET_eff_signal_mcsystup;
	TH2F* h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown;
  TH2F* h2_DiJetInvMass_vs_MET_background;
	TH2F* h2_DiJetInvMass_vs_MET_background_stat;
	TH2F* h2_DiJetInvMass_vs_MET_background_mcsystup;
	TH2F* h2_DiJetInvMass_vs_MET_background_vbfsystup;
	TH2F* h2_DiJetInvMass_vs_MET_background_mcsystdown;
	TH2F* h2_DiJetInvMass_vs_MET_background_vbfsystdown;


  h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot(taupt, "Taui2TightIso", chi, lsp);
	h2_DiJetInvMass_vs_MET_eff_signal_stat = makeEffPlotStatUnc(taupt, "Taui2TightIso", chi, lsp);
	h2_DiJetInvMass_vs_MET_eff_signal_mcsystup = makeEffPlotSyst(taupt, "Taui2TightIso", chi, lsp, 0.5);
	h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown = makeEffPlotSyst(taupt, "Taui2TightIso", chi, lsp, -0.5);
	h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT(taupt, "baseline");
	h2_DiJetInvMass_vs_MET_background_stat = makeBackgroundPlot_LtoT_StatUnc(taupt, "baseline");
	h2_DiJetInvMass_vs_MET_background_mcsystup = makeBackgroundPlot_LtoT_MCSyst(taupt, "baseline", 0.5);
	h2_DiJetInvMass_vs_MET_background_mcsystdown = makeBackgroundPlot_LtoT_MCSyst(taupt, "baseline", -0.5);
	h2_DiJetInvMass_vs_MET_background_vbfsystup = makeBackgroundPlot_LtoT_VBFSyst(taupt, "baseline", 0.175);
	h2_DiJetInvMass_vs_MET_background_vbfsystdown = makeBackgroundPlot_LtoT_VBFSyst(taupt, "baseline", -0.076);

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
	TH2F* h2_DiJetInvMass_vs_MET_xsec_mcsystup;
	TH2F* h2_DiJetInvMass_vs_MET_xsec_mcsystdown;
	TH2F* h2_DiJetInvMass_vs_MET_xsec_vbfsystup;
	TH2F* h2_DiJetInvMass_vs_MET_xsec_vbfsystdown;
	h2_DiJetInvMass_vs_MET_xsec_stat = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(),("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(), nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_mcsystup = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(),("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(), nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_mcsystdown = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(),("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(), nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_vbfsystup = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(),("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(), nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_vbfsystdown = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(),("JetInvMass_vs_MET_xsec_" + chi + "_" + lsp + "_" + taupt).c_str(), nbinsx, 0., 240., nbinsy , 0., 2500.);


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

				double xsec_mcsystup = getXSection(
						getSignalEfficiency(i,j, h2_DiJetInvMass_vs_MET_eff_signal_mcsystup),
						lumi,
						getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background_mcsystup),
						2.5, //sigma
						2.0 //significance
						);
				double xsec_mcsystdown = getXSection(
						getSignalEfficiency(i,j, h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown),
						lumi,
						getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background_mcsystdown),
						2.5, //sigma
						2.0 //significance
						);
				double xsec_vbfsystup = getXSection(
						getSignalEfficiency(i,j, h2_DiJetInvMass_vs_MET_eff_signal),
						lumi,
						getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background_vbfsystup),
						2.5, //sigma
						2.0 //significance
						);
				double xsec_vbfsystdown = getXSection(
						getSignalEfficiency(i,j, h2_DiJetInvMass_vs_MET_eff_signal),
						lumi,
						getBackgroudEvents( i, j, h2_DiJetInvMass_vs_MET_background_vbfsystdown),
						2.5, //sigma
						2.0 //significance
						);
						h2_DiJetInvMass_vs_MET_xsec->SetBinContent(i, j, xsec);
						h2_DiJetInvMass_vs_MET_xsec_stat->SetBinContent(i, j, xsec_statunc);
						h2_DiJetInvMass_vs_MET_xsec_mcsystup->SetBinContent(i, j, xsec_mcsystup);
						h2_DiJetInvMass_vs_MET_xsec_vbfsystup->SetBinContent(i, j, xsec_vbfsystup);
						h2_DiJetInvMass_vs_MET_xsec_mcsystdown->SetBinContent(i, j, xsec_mcsystdown);
						h2_DiJetInvMass_vs_MET_xsec_vbfsystdown->SetBinContent(i, j, xsec_vbfsystdown);
		}
	}

	//Cross section minimum
	int binmin = h2_DiJetInvMass_vs_MET_xsec->GetMinimumBin();
	double minimum = FLT_MAX;
	double x_min = 0.;
	double y_min = 0.;
	int i_min = -1;
	int j_min = -1;
	double taupt_value = 0.;
	if (taupt == "taupt20") {taupt_value = 20.;} else
	if (taupt == "taupt25") {taupt_value = 25.;} else
	if (taupt == "taupt30") {taupt_value = 30.;} else
	if (taupt == "taupt35") {taupt_value = 35.;} else
	if (taupt == "taupt40") {taupt_value = 40.;} else
	if (taupt == "taupt45") taupt_value = 45.;

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			if ((h2_DiJetInvMass_vs_MET_xsec->GetBinContent(i,j) < minimum) && (h2_DiJetInvMass_vs_MET_xsec->GetBinContent(i,j) > 0)){
				minimum = h2_DiJetInvMass_vs_MET_xsec->GetBinContent(i,j);
				x_min = h2_DiJetInvMass_vs_MET_xsec->GetXaxis()->GetBinCenter(i);
				x_min = x_min - 0.5 * h2_DiJetInvMass_vs_MET_xsec->GetXaxis()->GetBinWidth(i);
				y_min = h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->GetBinCenter(j);
				y_min = y_min - 0.5 * h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->GetBinWidth(i);
				i_min = i;
				j_min = j;
			}
		}
	}
	double minimum_statunc = h2_DiJetInvMass_vs_MET_xsec_stat->GetBinContent(i_min,j_min);
	double minimum_mcsystup_var = h2_DiJetInvMass_vs_MET_xsec_mcsystup->GetBinContent(i_min,j_min) - minimum;
	double minimum_mcsystdown_var = minimum - h2_DiJetInvMass_vs_MET_xsec_mcsystdown->GetBinContent(i_min,j_min);
	double minimum_vbfsystup_var = h2_DiJetInvMass_vs_MET_xsec_vbfsystup->GetBinContent(i_min,j_min) - minimum;
	double minimum_vbfsystdown_var = minimum - h2_DiJetInvMass_vs_MET_xsec_vbfsystdown->GetBinContent(i_min,j_min);

	cout << " for " << chi << ", " << lsp << " and " << taupt << endl;
	cout << "$"<<minimum << "\\pm"<< minimum_statunc << "^{+" << minimum_mcsystup_var << " + " << minimum_vbfsystup_var << "}_{-" << minimum_mcsystdown_var << "-" << minimum_vbfsystdown_var << "}$ & $<$ " << taupt_value << " & $<$ " << y_min << "  & $<$ " << x_min << " \\\\ " << endl;


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
