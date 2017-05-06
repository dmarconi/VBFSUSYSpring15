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

void th2fDebugReader(TH2F* &inputplot) {

	cout << "-------------------(TH2FDEBUGREADER) Fetching values from " << inputplot->GetName() << "-------------------" << endl;

	int nbinsx = inputplot->GetNbinsX();
	int nbinsy = inputplot->GetNbinsY();

	for (int j = 0; j < nbinsy; j++) {

		for (int i = 0; i < nbinsx; i++) {
			double value = inputplot->GetBinContent( i, j);
			cout << value << " ";
		}
		cout << endl;
	}

	cout << "--------------------------------------------(TH2FDEBUGREADER) END-------------------------------------------" << endl;

}

//calculates the VBF efficiency
double vbfefficiency(double evenCRcounts, double oddCRcounts){

	return ( (oddCRcounts) / ( (oddCRcounts) + (evenCRcounts) ));
}

//calculates the VBF efficiency stat uncertainty
double vbfefficiency_statunc(double evenCRcounts,double evenCRcounts_statunc,
	double oddCRcounts, double oddCRcounts_statunc){

		return (  (  (evenCRcounts * oddCRcounts_statunc)  +  (oddCRcounts * evenCRcounts_statunc)  )  /
		(  (evenCRcounts + oddCRcounts ) * (evenCRcounts + oddCRcounts )  )  );
	}

	//calculates the VBF conversion factor
	double vbfConversionFactor(TFile* &inputfile, string isoregion) {

		TH1F* h_ditaucharge;
		TH1F* h_ditauchargeVBFinverted;

		h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
		h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

		double counts = h_ditaucharge->GetBinContent(3);
		double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3);


		double vbfeff = vbfefficiency(counts_inverted, counts);

		return (vbfeff / (1. - vbfeff));

		//h_ditaucharge->Delete();
		//h_ditauchargeVBFinverted->Delete();

	}

	//calculates the VBF conversion factor
	//statistical uncertainty
	double vbfConversionFactor_statunc(TFile* &inputfile, string isoregion) {

		TH1F* h_ditaucharge;
		TH1F* h_ditauchargeVBFinverted;

		h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
		h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

		double counts = h_ditaucharge->GetBinContent(3);
		double counts_statunc = h_ditaucharge->GetBinError(3);
		double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3);
		double counts_inverted_statunc = h_ditauchargeVBFinverted->GetBinError(3);


		double vbfeff = vbfefficiency(counts_inverted, counts);
		double vbfeff_statunc = vbfefficiency_statunc(counts_inverted,counts_inverted_statunc, counts, counts_statunc);

		return ((2. * vbfeff * vbfeff_statunc) / ( (1. - vbfeff) * (1. - vbfeff) ) );

		h_ditaucharge->Delete();
		h_ditauchargeVBFinverted->Delete();
	}

	//calculates the VBF conversion
	//Monte Carlo systematic uncertainty
	double vbfConversionFactorMCSyst(TFile* &inputfile, string isoregion, double variation) {

		TH1F* h_ditaucharge;
		TH1F* h_ditauchargeVBFinverted;

		h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
		h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

		double counts = h_ditaucharge->GetBinContent(3) + (variation * h_ditaucharge->GetBinContent(3));
		double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3) + (variation * h_ditauchargeVBFinverted->GetBinContent(3));


		double vbfeff = vbfefficiency(counts_inverted, counts);

		return (vbfeff / (1. - vbfeff));

		h_ditaucharge->Delete();
		h_ditauchargeVBFinverted->Delete();

	}

	//calculates the VBF conversion
	//VBF conversion systematic uncertainty
	double vbfConversionFactoVBFSyst(TFile* &inputfile, string isoregion, double variation) {

		TH1F* h_ditaucharge;
		TH1F* h_ditauchargeVBFinverted;

		h_ditaucharge = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h_ditaucharge").c_str())));
		h_ditauchargeVBFinverted = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

		double counts = h_ditaucharge->GetBinContent(3);
		double counts_inverted = h_ditauchargeVBFinverted->GetBinContent(3);


		double vbfeff = (1. + variation) * vbfefficiency(counts_inverted, counts);

		return (vbfeff / (1. - vbfeff));

		h_ditaucharge->Delete();
		h_ditauchargeVBFinverted->Delete();

	}

	//calculates the LtoT conversion factor
	double LtoTfactor(TFile* &inputfile) {

		TH1F* h1_counts;
		h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

		double oneTauMatchcounts = h1_counts->GetBinContent(3);
		double jetAlsoTcounts = h1_counts->GetBinContent(4);

		double jetToTightProb = jetAlsoTcounts / oneTauMatchcounts;
		double twoLooseTo2Tight = jetToTightProb * jetToTightProb;

		return (twoLooseTo2Tight);

		h1_counts->Delete();

	}

	//calculates the LtoT conversion factor
	//statistical uncertainty
	double LtoTfactor_statunc(TFile* &inputfile) {

		TH1F* h1_counts;
		h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

		double oneTauMatchcounts = h1_counts->GetBinContent(3);
		double oneTauMatchcounts_statunc = h1_counts->GetBinError(3);
		double jetAlsoTcounts = h1_counts->GetBinContent(4);
		double jetAlsoTcounts_statunc = h1_counts->GetBinError(4);

		double jetToTightProb = jetAlsoTcounts / oneTauMatchcounts;
		double jetToTightProb_statunc = (jetAlsoTcounts / oneTauMatchcounts) *
																	sqrt( ( (oneTauMatchcounts_statunc * oneTauMatchcounts_statunc)
																	/ (oneTauMatchcounts * oneTauMatchcounts) ) +
																	( (jetAlsoTcounts_statunc * jetAlsoTcounts_statunc) / (jetAlsoTcounts * jetAlsoTcounts) ) );
		double twoLooseTo2Tight = jetToTightProb * jetToTightProb;
		double twoLooseTo2Tight_statunc = 2. * jetToTightProb * jetToTightProb_statunc;

		return (twoLooseTo2Tight_statunc);

		h1_counts->Delete();

	}

	//calculates the LtoT conversion factor
	//systematic uncertainty
	double LtoTfactorSyst(TFile* &inputfile, double variation) {

		TH1F* h1_counts;
		h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

		double oneTauMatchcounts = h1_counts->GetBinContent(3) +
															(variation * h1_counts->GetBinContent(3));
		double jetAlsoTcounts = h1_counts->GetBinContent(4) +
														(variation * h1_counts->GetBinContent(4));

		double jetToTightProb = jetAlsoTcounts / oneTauMatchcounts;
		double twoLooseTo2Tight = jetToTightProb * jetToTightProb;

		return (twoLooseTo2Tight);

		h1_counts->Delete();

	}

	//calculates the signal events for a given
	//histogram bin
	double getSignalEfficiency(int xbin, int ybin, TH2F* signal_map){
		double efficiency = signal_map->GetBinContent(xbin,ybin);
		return (efficiency);
	}

	//calculates the signal events stat uncertainty
	//for a given histogram bin
	double getSignalEfficiencyStatUnc(int xbin, int ybin, TH2F* signal_map){
		double efficiency = signal_map->GetBinContent(xbin,ybin);
		return (efficiency);
	}

	//calculates the background events for a given
	//histogram bin
	double getBackgroudEvents(int xbin, int ybin, TH2F* background_map){
		int nbinsx = background_map->GetNbinsX();
		int nbinsy = background_map->GetNbinsY();
		double events = background_map->Integral( xbin, nbinsx, ybin, nbinsy );
		return (events);
	}

	//calculates the background events stat uncertainty
	//for a given histogram bin
	double getBackgroudEventsStatUnc(int xbin, int ybin, TH2F* background_map_unc){
		int nbinsx = background_map_unc->GetNbinsX();
		int nbinsy = background_map_unc->GetNbinsY();
		double squaredsum = 0.;

		for (int i = xbin; i < nbinsx; i++) {
			for (int j = ybin; j < nbinsy; j++) {
				squaredsum += (background_map_unc->GetBinContent(i,j)) *
				(background_map_unc->GetBinContent(i,j));
			}
		}

		return (sqrt(squaredsum));
	}


	//calculates the cross section limit
	double getXSection(double efficiency, double luminosity, double background, double sigma, double significance) {
		return ( (significance * ( (0.5 * sigma) +
		sqrt( background + (0.5 * background * background))  ) ) /
		(efficiency * luminosity)
	);
}

//calculates the cross section limit
//stat uncertainty
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

//creates the signal efficiency plot for
//a given signal sample and selection
TH2F* makeEffPlot(string taupt, string isoregion, string chi, string stau, string lsp, TFile* &inputfile) {

	TH2F* h2_DiJetInvMass_vs_MET;
	TH1F* h_count;

	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h2_DiJetInvMass_vs_MET").c_str())));
	h_count = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1);

	TH2F* h2_DiJetInvMass_vs_MET_eff;
	h2_DiJetInvMass_vs_MET_eff = new TH2F (("h2_DiJetInvMass_vs_MET_eff" + chi + "_" + stau + "_" + lsp + "_" + taupt + "_stat").c_str(),
	("h2_DiJetInvMass_vs_MET_eff" + chi + "_" + stau + "_" + lsp + "_" + taupt + "_stat").c_str(),
	nbinsx, 0., 240., nbinsy , 0., 2500.);
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

	h2_DiJetInvMass_vs_MET->Delete();
	h_count->Delete();
	h2_DiJetInvMass_vs_MET_eff->Delete();

}

//creates the signal efficiency statistical
//uncertainty plot for a given signal sample and selection
TH2F* makeEffPlotStatUnc(string taupt, string isoregion, string chi, string stau, string lsp, TFile* &inputfile) {

	TH2F* h2_DiJetInvMass_vs_MET;
	TH1F* h_count;

	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h2_DiJetInvMass_vs_MET").c_str())));
	h_count = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1);
	double ntotalevents_statunc = h_count->GetBinError(1);

	TH2F* h2_DiJetInvMass_vs_MET_eff_statunc;
	h2_DiJetInvMass_vs_MET_eff_statunc = new TH2F (("h2_DiJetInvMass_vs_MET_eff_statunc_" + chi + "_" + stau + "_" + lsp + "_" + taupt + "_stat").c_str(),
																			("h2_DiJetInvMass_vs_MET_eff_statunc_" + chi + "_" + stau + "_" + lsp + "_" + taupt + "_stat").c_str(),
																			nbinsx, 0., 240., nbinsy , 0., 2500.);
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

	h2_DiJetInvMass_vs_MET->Delete();
	h_count->Delete();
	h2_DiJetInvMass_vs_MET_eff_statunc->Delete();

}

//creates the signal efficiency systematic
//uncertainty plot for a given signal sample and selection
TH2F* makeEffPlotSyst(string taupt, string isoregion, string chi, string stau, string lsp, double variation, TFile* &inputfile) {

	string var_label;
	if (variation > 00) var_label = "mcsystup";
	else if (variation < 00) var_label = "mcsystdown";
	else if (variation == 0.) cout << "ERROR: something's wrong with the syst uncertainty variation" << endl;

	TH2F* h2_DiJetInvMass_vs_MET;
	TH1F* h_count;

	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/h2_DiJetInvMass_vs_MET").c_str())));
	h_count = ((TH1F*)(inputfile->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1) + (variation * h_count->GetBinContent(1));

	TH2F* h2_DiJetInvMass_vs_MET_eff_syst;
	h2_DiJetInvMass_vs_MET_eff_syst = new TH2F (("h2_DiJetInvMass_vs_MET_eff_" + var_label + "_" + chi + "_" + stau + "_" + lsp + "_" + taupt + "_stat").c_str(),
																							("h2_DiJetInvMass_vs_MET_eff_" + var_label + "_" + chi + "_" + stau + "_" + lsp + "_" + taupt + "_stat").c_str(),
																							nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_eff_syst->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_eff_syst->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_eff_syst->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double integral = h2_DiJetInvMass_vs_MET->Integral( i, nbinsx, j, nbinsy ) + (variation * h2_DiJetInvMass_vs_MET->Integral( i, nbinsx, j, nbinsy ));
			h2_DiJetInvMass_vs_MET_eff_syst->SetBinContent(i,j, (integral/ntotalevents));
		}
	}
	return h2_DiJetInvMass_vs_MET_eff_syst;

	h2_DiJetInvMass_vs_MET->Delete();
	h_count->Delete();
	h2_DiJetInvMass_vs_MET_eff_syst->Delete();

}

//creates a histogram storing the LtoT
//conversion factor for a given selection
TH2F* makeBackgroundPlot_LtoT(string isoregion, TFile* &inputfile_qcd, TFile* &inputfile_ltot){

	TH2F* h2_DiJetInvMass_vs_MET;
	TH2F* h2_DiJetInvMass_vs_MET_LtoT;
	TH1F* h_ditauchargeVBFinverted;
	TH1F* h_count;


	if (inputfile_qcd->IsZombie()) cout << "ERROR: couldn't open inputfile_qcd file for taupt" << endl;
	if (inputfile_ltot->IsZombie()) cout << "ERROR: couldn't open inputfile_ltot file for taupt" << endl;

	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile_qcd->Get(("demo/" + isoregion +"VBFInvertedSelection/h2_DiJetInvMass_vs_MET").c_str())));
	h_count = ((TH1F*)(inputfile_qcd->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile_qcd->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

	th2fDebugReader(h2_DiJetInvMass_vs_MET);

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
			double vbffactor = vbfConversionFactor(inputfile_qcd, isoregion);
			double ltotfactor = LtoTfactor(inputfile_ltot);
			h2_DiJetInvMass_vs_MET_LtoT->SetBinContent(i,j, (bincontent*vbffactor*ltotfactor));

		}
	}


	return h2_DiJetInvMass_vs_MET_LtoT;

	h2_DiJetInvMass_vs_MET->Delete();
	h2_DiJetInvMass_vs_MET_LtoT->Delete();
	h_ditauchargeVBFinverted->Delete();
	h_count->Delete();

}

//creates a histogram storing the LtoT factor
//statistical uncertainty for a given selection
TH2F* makeBackgroundPlot_LtoT_StatUnc(string isoregion, TFile* &inputfile_qcd, TFile* &inputfile_ltot){

	TH2F* h2_DiJetInvMass_vs_MET;
	TH2F* h2_DiJetInvMass_vs_MET_LtoT_statunc;
	TH1F* h_ditauchargeVBFinverted;
	TH1F* h_count;

	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile_qcd->Get(("demo/" + isoregion +"VBFInvertedSelection/h2_DiJetInvMass_vs_MET").c_str())));
	h_count = ((TH1F*)(inputfile_qcd->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile_qcd->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

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
			double vbffactor = vbfConversionFactor(inputfile_qcd, isoregion);
			double vbffactor_statunc = vbfConversionFactor_statunc(inputfile_qcd, isoregion);
			double ltotfactor = LtoTfactor(inputfile_ltot);
			double ltotfactor_statunc = LtoTfactor_statunc(inputfile_ltot);

			h2_DiJetInvMass_vs_MET_LtoT_statunc->SetBinContent(i,j,
				sqrt(pow( (bincontent_statunc*vbffactor*ltotfactor) , 2.) + pow( (bincontent*vbffactor_statunc*ltotfactor) , 2.) + pow( (bincontent*vbffactor*ltotfactor_statunc) , 2.) )
			);
		}
	}

	return h2_DiJetInvMass_vs_MET_LtoT_statunc;

	h2_DiJetInvMass_vs_MET->Delete();
	h2_DiJetInvMass_vs_MET_LtoT_statunc->Delete();
	h_ditauchargeVBFinverted->Delete();
	h_count->Delete();

}


//creates a histogram storing the LtoT factor
//Monte Carlo systematic uncertainty for a given selection
TH2F* makeBackgroundPlot_LtoT_MCSyst(string isoregion, TFile* &inputfile_qcd, TFile* &inputfile_ltot, double variation){

	TH2F* h2_DiJetInvMass_vs_MET;
	TH2F* h2_DiJetInvMass_vs_MET_LtoT;
	TH1F* h_ditauchargeVBFinverted;
	TH1F* h_count;

	string var_label;
	if (variation > 00) var_label = "mcsystup";
	else if (variation < 00) var_label = "mcsystdown";
	else if (variation == 0.) cout << "ERROR: something's wrong with the syst uncertainty variation" << endl;

	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile_qcd->Get(("demo/" + isoregion +"VBFInvertedSelection/h2_DiJetInvMass_vs_MET").c_str())));
	h_count = ((TH1F*)(inputfile_qcd->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile_qcd->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1) + (variation * h_count->GetBinContent(1));

	h2_DiJetInvMass_vs_MET_LtoT = new TH2F (
									("h2_DiJetInvMass_vs_MET_LtoT_" + var_label).c_str(),
									("h2_DiJetInvMass_vs_MET_LtoT_" + var_label).c_str(),
									nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_LtoT->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double bincontent = h2_DiJetInvMass_vs_MET->GetBinContent(i,j) + (variation * h2_DiJetInvMass_vs_MET->GetBinContent(i,j));
			double vbffactor = vbfConversionFactorMCSyst(inputfile_qcd, isoregion, variation);
			double ltotfactor = LtoTfactorSyst(inputfile_ltot,variation);
			h2_DiJetInvMass_vs_MET_LtoT->SetBinContent(i,j, (bincontent*vbffactor*ltotfactor));

		}
	}

	return h2_DiJetInvMass_vs_MET_LtoT;

	h2_DiJetInvMass_vs_MET->Delete();
	h2_DiJetInvMass_vs_MET_LtoT->Delete();
	h_ditauchargeVBFinverted->Delete();
	h_count->Delete();

}

//creates a histogram storing the LtoT factor
//VBF systematic uncertainty for a given selection
TH2F* makeBackgroundPlot_LtoT_VBFSyst(string isoregion, TFile* &inputfile_qcd, TFile* &inputfile_ltot, double variation){

	TH2F* h2_DiJetInvMass_vs_MET;
	TH2F* h2_DiJetInvMass_vs_MET_LtoT;
	TH1F* h_ditauchargeVBFinverted;
	TH1F* h_count;
	string var_label;
	if (variation > 00) var_label = "vbfsystup";
	else if (variation < 00) var_label = "vbfsystdown";
	else if (variation == 0.) cout << "ERROR: something's wrong with the syst uncertainty variation" << endl;

	h2_DiJetInvMass_vs_MET = ((TH2F*)(inputfile_qcd->Get(("demo/" + isoregion +"VBFInvertedSelection/h2_DiJetInvMass_vs_MET").c_str())));
	h_count = ((TH1F*)(inputfile_qcd->Get(("demo/" + isoregion +"ObjectSelection/counts").c_str())));
	h_ditauchargeVBFinverted = ((TH1F*)(inputfile_qcd->Get(("demo/" + isoregion +"VBFInvertedSelection/h_ditaucharge").c_str())));

	int nbinsx = h2_DiJetInvMass_vs_MET->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET->GetNbinsY();
	double ntotalevents = h_count->GetBinContent(1) + (variation * h_count->GetBinContent(1));

	h2_DiJetInvMass_vs_MET_LtoT = new TH2F (
									("h2_DiJetInvMass_vs_MET_LtoT_" + var_label).c_str(),
									("h2_DiJetInvMass_vs_MET_LtoT_" + var_label).c_str(),
									nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_LtoT->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_LtoT->SetStats(0);

	for (int i = 0; i < nbinsx; i++) {

		for (int j = 0; j < nbinsy; j++) {
			double bincontent = h2_DiJetInvMass_vs_MET->GetBinContent(i,j);
			double vbffactor = vbfConversionFactoVBFSyst(inputfile_qcd, isoregion, variation);
			double ltotfactor = LtoTfactor(inputfile_ltot);
			h2_DiJetInvMass_vs_MET_LtoT->SetBinContent(i,j, (bincontent*vbffactor*ltotfactor));

		}
	}

	return h2_DiJetInvMass_vs_MET_LtoT;

	h2_DiJetInvMass_vs_MET->Delete();
	h2_DiJetInvMass_vs_MET_LtoT->Delete();
	h_ditauchargeVBFinverted->Delete();
	h_count->Delete();

}

//calculates the cross section limit minimum for a given
//signal sample and selection and creates a LaTex output
void makeXSectionUnc(string taupt, string chi, string stau, string lsp, double lumi, bool debug) {

	if (debug == true) cout << "INPUT(makeXSectionUnc)"
	<< " taupt = " << taupt
	<< " chi = " << chi
	<< " stau = " << stau
	<< " lsp = " << lsp
	<< " lumi = " << lumi
	<< " debug = " << debug
	<< endl;

	//open input file
	if (debug == true) cout << "Opening signal input file..." << endl;

	TFile *inputfile_signal;
	if ((chi == "chi100") && (stau == "stau095") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau095_Chargino100_1M.root").c_str());
	else if ((chi == "chi200") && (stau == "stau195") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau195_Chargino200_1M.root").c_str());
	else if ((chi == "chi300") && (stau == "stau295") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M.root").c_str());
	else if ((chi == "chi400") && (stau == "stau395") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau395_Chargino400_1M.root").c_str());
	else if ((chi == "chi500") && (stau == "stau495") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau495_Chargino500_1M.root").c_str());

	else if ((chi == "chi100") && (stau == "stau195") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau095_Chargino100_1M.root").c_str());
	else if ((chi == "chi200") && (stau == "stau195") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau195_Chargino200_1M.root").c_str());
	else if ((chi == "chi300") && (stau == "stau295") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau295_Chargino300_1M.root").c_str());
	else if ((chi == "chi400") && (stau == "stau395") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau395_Chargino400_1M.root").c_str());
	else if ((chi == "chi500") && (stau == "stau495") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau495_Chargino500_1M.root").c_str());

	else if ((chi == "chi100") && (stau == "stau050") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau050_Chargino100_1M.root").c_str());
	else if ((chi == "chi200") && (stau == "stau100") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau100_Chargino200_1M.root").c_str());
	else if ((chi == "chi300") && (stau == "stau150") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau150_Chargino300_1M.root").c_str());
	else if ((chi == "chi400") && (stau == "stau200") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau200_Chargino400_1M.root").c_str());
	else if ((chi == "chi500") && (stau == "stau250") && (lsp == "lsp000")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau250_Chargino500_1M.root").c_str());

	else if ((chi == "chi100") && (stau == "stau075") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau075_Chargino100_1M.root").c_str());
	else if ((chi == "chi200") && (stau == "stau125") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau125_Chargino200_1M.root").c_str());
	else if ((chi == "chi300") && (stau == "stau175") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau175_Chargino300_1M.root").c_str());
	else if ((chi == "chi400") && (stau == "stau225") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau225_Chargino400_1M.root").c_str());
	else if ((chi == "chi500") && (stau == "stau275") && (lsp == "lsp050")) inputfile_signal = TFile::Open((taupt + "/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP050_Stau275_Chargino500_1M.root").c_str());
	else if (debug == true) cout << "ERROR: bechmark point not found" << endl;

	if (inputfile_signal->IsZombie()) cout << "ERROR: couldn't open file for chi" << chi << " stau" << stau << " lsp" << lsp << endl;

	//definition and filling of necessary input plots

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

	if (debug == true) cout << "Creating signal efficiency related plots..." <<endl;
	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_eff_signal" <<endl;
	h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot(taupt, "Taui2TightIso", chi, stau, lsp, inputfile_signal);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_eff_signal);
	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_eff_signal_stat" <<endl;
	h2_DiJetInvMass_vs_MET_eff_signal_stat = makeEffPlotStatUnc(taupt, "Taui2TightIso", chi, stau, lsp, inputfile_signal);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_eff_signal_stat);
	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_eff_signal_mcsystup" <<endl;
	h2_DiJetInvMass_vs_MET_eff_signal_mcsystup = makeEffPlotSyst(taupt, "Taui2TightIso", chi, stau, lsp, 0.5, inputfile_signal);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_eff_signal_mcsystup);
	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown" <<endl;
	h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown = makeEffPlotSyst(taupt, "Taui2TightIso", chi, stau, lsp, -0.5, inputfile_signal);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown);

	if (debug == true) cout << "Creating background related plots..." <<endl;

	if (debug == true) cout << "Opening signal input files..." << endl;
	TFile* inputfile_qcd = TFile::Open((taupt + "/allQCD_"+ taupt +".root").c_str());
	if (inputfile_qcd->IsZombie()) cout << "ERROR: couldn't open allQCD file for taupt" << taupt << endl;

	TFile* inputfile_ltot = TFile::Open(( taupt + "/allQCD_"+ taupt +"_LtoT.root").c_str());
	if (inputfile_ltot->IsZombie()) cout << "ERROR: couldn't open allQCD_LtoT file for taupt" << taupt << endl;

	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_background" <<endl;
	h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT("baseline", inputfile_qcd, inputfile_ltot);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_background);
	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_background_stat" <<endl;
	h2_DiJetInvMass_vs_MET_background_stat = makeBackgroundPlot_LtoT_StatUnc("baseline", inputfile_qcd, inputfile_ltot);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_stat);
	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_background_mcsystup" <<endl;
	h2_DiJetInvMass_vs_MET_background_mcsystup = makeBackgroundPlot_LtoT_MCSyst("baseline", inputfile_qcd, inputfile_ltot, 0.5);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_mcsystup);
	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_background_mcsystdown" <<endl;
	h2_DiJetInvMass_vs_MET_background_mcsystdown = makeBackgroundPlot_LtoT_MCSyst("baseline", inputfile_qcd, inputfile_ltot, -0.5);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_mcsystdown);
	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_background_vbfsystup" <<endl;
	h2_DiJetInvMass_vs_MET_background_vbfsystup = makeBackgroundPlot_LtoT_VBFSyst("baseline", inputfile_qcd, inputfile_ltot, 0.175);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_vbfsystup);
	if (debug == true) cout << "...h2_DiJetInvMass_vs_MET_background_vbfsystdown" <<endl;
	h2_DiJetInvMass_vs_MET_background_vbfsystdown = makeBackgroundPlot_LtoT_VBFSyst("baseline", inputfile_qcd, inputfile_ltot, -0.076);
	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_vbfsystdown);

	if (debug == true) cout << "Creating xsec related plots..." <<endl;

	if (debug == true) th2fDebugReader(h2_DiJetInvMass_vs_MET_background);

	int nbinsx = h2_DiJetInvMass_vs_MET_background->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET_background->GetNbinsY();

	TH2F* h2_DiJetInvMass_vs_MET_xsec;
	h2_DiJetInvMass_vs_MET_xsec = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
	("JetInvMass_vs_MET_xsec_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
	nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec->SetTitle("CMS Work");
	h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->SetTitle("M_{(jet,jet)} [GeV]");
	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetTitle("#sigma_{lim} pb");
	h2_DiJetInvMass_vs_MET_xsec->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
	h2_DiJetInvMass_vs_MET_xsec->GetYaxis()->SetTitleOffset(1.40);
	h2_DiJetInvMass_vs_MET_xsec->GetZaxis()->SetTitleOffset(1.13);
	h2_DiJetInvMass_vs_MET_xsec->SetStats(0);

	//declaration of the cross section and
	// its uncertainties plots

	TH2F* h2_DiJetInvMass_vs_MET_xsec_stat;
	TH2F* h2_DiJetInvMass_vs_MET_xsec_mcsystup;
	TH2F* h2_DiJetInvMass_vs_MET_xsec_mcsystdown;
	TH2F* h2_DiJetInvMass_vs_MET_xsec_vbfsystup;
	TH2F* h2_DiJetInvMass_vs_MET_xsec_vbfsystdown;

	h2_DiJetInvMass_vs_MET_xsec_stat = new TH2F (
			("JetInvMass_vs_MET_xsec_stat_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			("JetInvMass_vs_MET_xsec_stat_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_mcsystup = new TH2F (
			("JetInvMass_vs_MET_xsec_mcsystup_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			("JetInvMass_vs_MET_xsec_mcsystup_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_mcsystdown = new TH2F (
			("JetInvMass_vs_MET_xsec_mcsystdown_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			("JetInvMass_vs_MET_xsec_mcsystdown_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_vbfsystup = new TH2F (
			("JetInvMass_vs_MET_xsec_vbfsystup_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			("JetInvMass_vs_MET_xsec_vbfsystup_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_vbfsystdown = new TH2F (
			("JetInvMass_vs_MET_xsec_vbfsystdown_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			("JetInvMass_vs_MET_xsec_vbfsystdown_" + chi + "_" + stau + "_" + lsp + "_" + taupt).c_str(),
			nbinsx, 0., 240., nbinsy , 0., 2500.);

	//this loops cycles trough every xsec plots bin
	// and sets it to the proper value
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

	if (debug == true) cout << "Searching for the Xsec minimum..." <<endl;

	//Searchin the cross section minimum
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

	//minimum output in LaTex format
	cout << " for " << chi << ", " << lsp << " and " << taupt << endl;
	cout << "$"<<minimum << "\\pm"<< minimum_statunc << "^{+" << minimum_mcsystup_var << " + " << minimum_vbfsystup_var << "}_{-" << minimum_mcsystdown_var << "-" << minimum_vbfsystdown_var << "}$ & $<$ " << taupt_value << " & $<$ " << y_min << "  & $<$ " << x_min << " \\\\ " << endl;

	//clearing memory

	// inputfile_signal->Close();
	// inputfile_qcd->Close();
	// inputfile_ltot->Close();

	// inputfile_signal->Delete();
	// inputfile_qcd->Delete();
	// inputfile_ltot->Delete();

	h2_DiJetInvMass_vs_MET_eff_signal->Delete();
	h2_DiJetInvMass_vs_MET_eff_signal_stat->Delete();
	h2_DiJetInvMass_vs_MET_eff_signal_mcsystup->Delete();
	h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown->Delete();

	h2_DiJetInvMass_vs_MET_background->Delete();
	h2_DiJetInvMass_vs_MET_background_stat->Delete();
	h2_DiJetInvMass_vs_MET_background_mcsystup->Delete();
	h2_DiJetInvMass_vs_MET_background_vbfsystup->Delete();
	h2_DiJetInvMass_vs_MET_background_mcsystdown->Delete();
	h2_DiJetInvMass_vs_MET_background_vbfsystdown->Delete();

	h2_DiJetInvMass_vs_MET_xsec->Delete();
	h2_DiJetInvMass_vs_MET_xsec_stat->Delete();
	h2_DiJetInvMass_vs_MET_xsec_mcsystup->Delete();
	h2_DiJetInvMass_vs_MET_xsec_mcsystdown->Delete();
	h2_DiJetInvMass_vs_MET_xsec_vbfsystup->Delete();
	h2_DiJetInvMass_vs_MET_xsec_vbfsystdown->Delete();
}

//repeats the cross section study for
//a given signal sample
void makeXsecLimPlots(string chi, string stau, string lsp, double lumi, bool debug){

	if (debug == true) cout << "taupt20" << endl;
	makeXSectionUnc("taupt20",chi,stau,lsp,lumi,debug);
	if (debug == true) cout << "taupt25" << endl;
	makeXSectionUnc("taupt25",chi,stau,lsp,lumi,debug);
	if (debug == true) cout << "taupt30" << endl;
	makeXSectionUnc("taupt30",chi,stau,lsp,lumi,debug);
	if (debug == true) cout << "taupt35" << endl;
	makeXSectionUnc("taupt35",chi,stau,lsp,lumi,debug);
	if (debug == true) cout << "taupt40" << endl;
	makeXSectionUnc("taupt40",chi,stau,lsp,lumi,debug);
	if (debug == true) cout << "taupt45" << endl;
	makeXSectionUnc("taupt45",chi,stau,lsp,lumi,debug);


}

//main function. Runs the cross section
//uncertainty study over all the signal samples available
void fullXsecLimScan(){

	double luminosity =  85000.;
	bool debug = false;

	if (debug == true) cout << "Processing xsec limit for luminosity =" << luminosity <<endl;

	if (debug == true) cout << "Signal scenario: chi100 stau095 lsp000" <<endl;
	makeXsecLimPlots("chi100", "stau095", "lsp000", luminosity, debug);
	if (debug == true) cout << "Signal scenario: chi200 stau195 lsp000" <<endl;
	makeXsecLimPlots("chi200", "stau195", "lsp000", luminosity, debug);
	if (debug == true) cout << "Signal scenario: chi300 stau295 lsp000" <<endl;
	makeXsecLimPlots("chi300", "stau295", "lsp000", luminosity, debug);
	if (debug == true) cout << "Signal scenario: chi500 stau395 lsp000" <<endl;
	makeXsecLimPlots("chi400", "stau395", "lsp000", luminosity, debug);
	if (debug == true) cout << "Signal scenario: chi500 stau495 lsp000" <<endl;
	makeXsecLimPlots("chi500", "stau495", "lsp000", luminosity, debug);

	if (debug == true) cout << "Signal scenario: chi100 stau095 lsp050" <<endl;
	makeXsecLimPlots("chi100", "stau095", "lsp050", luminosity, debug);
	if (debug == true) cout << "Signal scenario: chi200 stau195 lsp050" <<endl;
	makeXsecLimPlots("chi200", "stau195", "lsp050", luminosity, debug);
	if (debug == true) cout << "Signal scenario: chi300 stau295 lsp050" <<endl;
	makeXsecLimPlots("chi300", "stau295", "lsp050", luminosity, debug);
	if (debug == true) cout << "Signal scenario: chi500 stau395 lsp050" <<endl;
	makeXsecLimPlots("chi400", "stau395", "lsp050", luminosity, debug);
	if (debug == true) cout << "Signal scenario: chi500 stau495 lsp050" <<endl;
	makeXsecLimPlots("chi500", "stau495", "lsp050", luminosity, debug);

}
