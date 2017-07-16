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
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <fstream>
#include <iostream>
#include <fstream>
#include <cmath>

//useful to store all
//xsec limit information
struct xsecLim {
	double value;
  double err_stat;
	double err_systVBFup;
	double err_systVBFdown;
	double err_systMCup;
	double err_systMCdown;
} ;

//transates the xsec limit result into string
string xsecLim_to_string(xsecLim xsecLim_min){

	string filestream;

	double err_systVBF;
	if (fabs(xsecLim_min.err_systVBFup) > fabs(xsecLim_min.err_systVBFdown)) err_systVBF = fabs(xsecLim_min.err_systVBFup);
	else err_systVBF = fabs(xsecLim_min.err_systVBFdown);

	double err_systMC;
	if (fabs(xsecLim_min.err_systMCup) > fabs(xsecLim_min.err_systMCdown)) err_systMC = fabs(xsecLim_min.err_systMCup);
	else err_systMC = fabs(xsecLim_min.err_systMCdown);

	double err_tot = xsecLim_min.err_stat + sqrt(pow( err_systVBF, 2.) + pow( err_systMC, 2.));

	filestream = std::to_string(xsecLim_min.value) + " " + std::to_string(err_tot);

	return filestream;
}

//read cross section values and error
//from txt file
void xsec_reader_values(string ifile, Int_t lines, vector<double> & vbfxsec_y, vector<double> & vbfxsec_y_err) {

	ifstream infile;
	infile.open(("results/txt/out_xsecmin_" + ifile + ".txt").c_str());

	double value;
	double error;

	for (int i = 1; i <= lines; i++){
		infile >> value >> error;
		vbfxsec_y.push_back(value);
		vbfxsec_y_err.push_back(error);
	}
}

//handmade TH2F reader
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

void plotXsecLim(string taupt, string chi, string stau, string lsp, TH2F* & h2_xsec_lim, bool debug) {

	TH2F* h2_xsec_lim_zoom;

	//Converting model point into doubles

	string chi_label;
	string stau_label;
	string lsp_label;
	string taupt_label;

	if (chi == "Chargino100") {chi_label = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 100 GeV";}
	else if (chi == "Chargino200") {chi_label = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 200 GeV";}
	else if (chi == "Chargino300") {chi_label = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 300 GeV";}
	else if (chi == "Chargino400") {chi_label = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 400 GeV";}
	else if (chi == "Chargino500") {chi_label = "m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) = 500 GeV";}

	if ((stau == "Stau095") || (stau == "Stau195") || (stau == "Stau295") ||
			(stau == "Stau395") || (stau == "Stau495")
		) {stau_label = "m(#tilde{#chi}^{#pm}_{1}) - m(#tilde{#tau}^{#pm}) = 5 GeV";}
	else {stau_label = "m(#tilde{#tau}^{#pm}) = m(#tilde{#chi}^{#pm}_{1})/2 + m(#tilde{#chi}^{0}_{1}/2";}

	if (lsp == "LSP000") {lsp_label = "m(#tilde{#chi}^{0}_{1}) = 0 GeV";}
	else if (lsp == "LSP050") {lsp_label = "m(#tilde{#chi}^{0}_{1}) = 50 GeV";}
	else if (lsp == "LSP150") {lsp_label = "m(#tilde{#chi}^{0}_{1}) = 150 GeV";}
	else if (lsp == "LSP250") {lsp_label = "m(#tilde{#chi}^{0}_{1}) = 250 GeV";}
	else if (lsp == "LSP350") {lsp_label = "m(#tilde{#chi}^{0}_{1}) = 350 GeV";}
	else if (lsp == "LSP450") {lsp_label = "m(#tilde{#chi}^{0}_{1}) = 450 GeV";}

	if (taupt == "taupt20") {taupt_label = "#tau_{Pt} = 20 GeV";}
	else if (taupt == "taupt25") {taupt_label = "#tau_{pt} = 25 GeV";}
	else if (taupt == "taupt30") {taupt_label = "#tau_{pt} = 30 GeV";}
	else if (taupt == "taupt35") {taupt_label = "#tau_{pt} = 35 GeV";}
	else if (taupt == "taupt40") {taupt_label = "#tau_{pt} = 40 GeV";}
	else if (taupt == "taupt45") {taupt_label = "#tau_{pt} = 45 GeV";}

	//defining legend
	TLegend* leg = new TLegend(0.63,0.68,0.86,0.9);
	leg->SetTextSize(0.02);
	leg->AddEntry((TObject*)0, (chi_label).c_str(), "");
	leg->AddEntry((TObject*)0, (stau_label).c_str(), "");
	leg->AddEntry((TObject*)0, (lsp_label).c_str(), "");
	leg->AddEntry((TObject*)0, (taupt_label).c_str(), "");

	TCanvas *my_canvas_zoom = new TCanvas("mycanvas_zoom","mycanvas_zoom",1024.,768.);
	my_canvas_zoom->cd();
	gPad->SetLogz();
	h2_xsec_lim_zoom = (TH2F*) h2_xsec_lim->Clone();
	h2_xsec_lim_zoom->GetXaxis()->SetRangeUser(30.,200.);
	h2_xsec_lim_zoom->GetYaxis()->SetRangeUser(0.,1000.);
	double maximum_zoom = h2_xsec_lim_zoom->GetMaximum();
	double minimum_zoom = h2_xsec_lim_zoom->GetMinimum();
	h2_xsec_lim_zoom ->GetZaxis()->SetRangeUser(minimum_zoom,maximum_zoom);
	h2_xsec_lim_zoom->Draw("colz");
	leg->Draw();
	gPad->SetRightMargin(0.14);
	gPad->SetLogz();
	gPad->Update();
	my_canvas_zoom->Print(("results/JetInvMass_vs_MET_xseclim_" + chi + "_" + stau + "_" + lsp + "_" + taupt + "_zoom.pdf").c_str());

	TCanvas *my_canvas = new TCanvas("mycanvas","mycanvas",1024.,768.);
	my_canvas->cd();

	h2_xsec_lim->GetZaxis()->SetRangeUser(0.001,1000);
	h2_xsec_lim->Draw("colz");
	leg->Draw();

	gPad->SetRightMargin(0.14);
	gPad->SetLogz();
	gPad->Update();

	my_canvas->SetRightMargin(0.14);
	my_canvas->Print(("results/JetInvMass_vs_MET_xseclim_" + chi + "_" + stau + "_" + lsp + "_" + taupt + ".pdf").c_str());

	my_canvas_zoom->Close();
	my_canvas->Close();
	h2_xsec_lim_zoom->Delete();

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

		h_ditaucharge->Delete();
		h_ditauchargeVBFinverted->Delete();

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
xsecLim findXsecLim(string cutcfg, string chi, string stau, string lsp, std::ofstream& out_latex, double lumi, bool debug) {

	if (debug == kTRUE) cout << "INPUT(findXsecLimMin)"
	<< " cutcfg = " << cutcfg
	<< " chi = " << chi
	<< " stau = " << stau
	<< " lsp = " << lsp
	<< " lumi = " << lumi
	<< " debug = " << debug
	<< endl;

	//xsecLim struct
	xsecLim xsecLim_min;

	//open input files
	if (debug == kTRUE) cout << "Opening inputs file..." << endl;

	TFile* inputfile_qcd = TFile::Open((cutcfg + "/allQCD.root").c_str());
	if (inputfile_qcd->IsZombie()) cout << "ERROR: couldn't open allQCD file for taupt" << cutcfg << endl;

	TFile* inputfile_ltot = TFile::Open(( cutcfg + "/allQCD_LtoT.root").c_str());
	if (inputfile_ltot->IsZombie()) cout << "ERROR: couldn't open allQCD_LtoT file for taupt" << cutcfg << endl;

	TFile *inputfile_signal = TFile::Open((cutcfg + "/VBFC1pmN2_C1ToTau_N2ToTauTau_" + lsp +"_" + stau + "_" + chi + "_1M.root").c_str());
	if (inputfile_signal->IsZombie()) cout << "ERROR: couldn't open file for " << chi << " " << stau << " " << lsp << endl;

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

	if (debug == kTRUE) cout << "Creating signal efficiency related plots..." <<endl;
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_eff_signal" <<endl;
	h2_DiJetInvMass_vs_MET_eff_signal = makeEffPlot(cutcfg, "Taui2TightIso", chi, stau, lsp, inputfile_signal);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_eff_signal);
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_eff_signal_stat" <<endl;
	h2_DiJetInvMass_vs_MET_eff_signal_stat = makeEffPlotStatUnc(cutcfg, "Taui2TightIso", chi, stau, lsp, inputfile_signal);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_eff_signal_stat);
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_eff_signal_mcsystup" <<endl;
	h2_DiJetInvMass_vs_MET_eff_signal_mcsystup = makeEffPlotSyst(cutcfg, "Taui2TightIso", chi, stau, lsp, 0.5, inputfile_signal);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_eff_signal_mcsystup);
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown" <<endl;
	h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown = makeEffPlotSyst(cutcfg, "Taui2TightIso", chi, stau, lsp, -0.5, inputfile_signal);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_eff_signal_mcsystdown);

	if (debug == kTRUE) cout << "Creating background related plots..." <<endl;
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_background" <<endl;
	h2_DiJetInvMass_vs_MET_background = makeBackgroundPlot_LtoT("baseline", inputfile_qcd, inputfile_ltot);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_background);
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_background_stat" <<endl;
	h2_DiJetInvMass_vs_MET_background_stat = makeBackgroundPlot_LtoT_StatUnc("baseline", inputfile_qcd, inputfile_ltot);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_stat);
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_background_mcsystup" <<endl;
	h2_DiJetInvMass_vs_MET_background_mcsystup = makeBackgroundPlot_LtoT_MCSyst("baseline", inputfile_qcd, inputfile_ltot, 0.5);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_mcsystup);
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_background_mcsystdown" <<endl;
	h2_DiJetInvMass_vs_MET_background_mcsystdown = makeBackgroundPlot_LtoT_MCSyst("baseline", inputfile_qcd, inputfile_ltot, -0.5);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_mcsystdown);
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_background_vbfsystup" <<endl;
	h2_DiJetInvMass_vs_MET_background_vbfsystup = makeBackgroundPlot_LtoT_VBFSyst("baseline", inputfile_qcd, inputfile_ltot, 0.175);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_vbfsystup);
	if (debug == kTRUE) cout << "...h2_DiJetInvMass_vs_MET_background_vbfsystdown" <<endl;
	h2_DiJetInvMass_vs_MET_background_vbfsystdown = makeBackgroundPlot_LtoT_VBFSyst("baseline", inputfile_qcd, inputfile_ltot, -0.076);
	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_background_vbfsystdown);

	if (debug == kTRUE) cout << "Creating xsec related plots..." <<endl;

	if (debug == kTRUE) th2fDebugReader(h2_DiJetInvMass_vs_MET_background);

	int nbinsx = h2_DiJetInvMass_vs_MET_background->GetNbinsX();
	int nbinsy = h2_DiJetInvMass_vs_MET_background->GetNbinsY();

	TH2F* h2_DiJetInvMass_vs_MET_xsec;
	h2_DiJetInvMass_vs_MET_xsec = new TH2F (("JetInvMass_vs_MET_xsec_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
	("JetInvMass_vs_MET_xsec_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
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
			("JetInvMass_vs_MET_xsec_stat_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
			("JetInvMass_vs_MET_xsec_stat_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
			nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_mcsystup = new TH2F (
			("JetInvMass_vs_MET_xsec_mcsystup_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
			("JetInvMass_vs_MET_xsec_mcsystup_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
			nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_mcsystdown = new TH2F (
			("JetInvMass_vs_MET_xsec_mcsystdown_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
			("JetInvMass_vs_MET_xsec_mcsystdown_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
			nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_vbfsystup = new TH2F (
			("JetInvMass_vs_MET_xsec_vbfsystup_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
			("JetInvMass_vs_MET_xsec_vbfsystup_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
			nbinsx, 0., 240., nbinsy , 0., 2500.);
	h2_DiJetInvMass_vs_MET_xsec_vbfsystdown = new TH2F (
			("JetInvMass_vs_MET_xsec_vbfsystdown_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
			("JetInvMass_vs_MET_xsec_vbfsystdown_" + chi + "_" + stau + "_" + lsp + "_" + cutcfg).c_str(),
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

	if (debug == kTRUE) cout << "Searching for the Xsec minimum..." <<endl;

	//Searchin the cross section minimum
	int binmin = h2_DiJetInvMass_vs_MET_xsec->GetMinimumBin();
	double minimum = FLT_MAX;
	double x_min = 0.;
	double y_min = 0.;
	int i_min = -1;
	int j_min = -1;

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
	if (debug == kTRUE) cout << "Writing LaTex output to file..." <<endl;
	out_latex << " for " << chi << ", " << stau << ", " << lsp << " and " << cutcfg << endl;
	out_latex << "$"<<minimum << "\\pm"<< minimum_statunc << "^{+"
						<< minimum_mcsystup_var << " + " << minimum_vbfsystup_var << "}_{-"
						<< minimum_mcsystdown_var << "-" << minimum_vbfsystdown_var << "}$ & $<$ "
						<< y_min << "  & $<$ " << x_min << " \\\\ "
						<< endl;

	xsecLim_min.value = minimum;
	xsecLim_min.err_stat = minimum_statunc;
	xsecLim_min.err_systVBFup = minimum_vbfsystup_var;
	xsecLim_min.err_systVBFdown = minimum_vbfsystdown_var;
	xsecLim_min.err_systMCup = minimum_mcsystup_var;
	xsecLim_min.err_systMCdown = minimum_mcsystdown_var;

	return (xsecLim_min);

	//clearing memory

	inputfile_signal->Close();
	inputfile_qcd->Close();
	inputfile_ltot->Close();

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

void fullXsecLimScan(){

  double luminosity =  85000.;
	Bool_t debug = kFALSE;

  xsecLim v_xsecLim;

  ofstream out_latex("results/txt/out_latex.txt");
  ofstream out_xsecmin_lspcomp_stauclose_cutcfg1("results/txt/out_xsecmin_lspcomp_stauclose_cutcfg1.txt");
  ofstream out_xsecmin_lspcomp_stauclose_cutcfg2("results/txt/out_xsecmin_lspcomp_stauclose_cutcfg2.txt");
  ofstream out_xsecmin_lspcomp_stauclose_cutcfg3("results/txt/out_xsecmin_lspcomp_stauclose_cutcfg3.txt");
  ofstream out_xsecmin_lspcomp_stauclose_cutcfg4("results/txt/out_xsecmin_lspcomp_stauclose_cutcfg4.txt");

  if (debug == kTRUE) cout << "Processing xsec limit for luminosity =" << luminosity <<endl;

  if (debug == kTRUE) cout << "Signal scenario: Chargino100 Stau075 LSP050" <<endl;
	out_xsecmin_lspcomp_stauclose_cutcfg1 << xsecLim_to_string(findXsecLim("cutcfg1","Chargino100","Stau075","LSP050",out_latex,luminosity,debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino200 Stau175 LSP150" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg1 << xsecLim_to_string(findXsecLim("cutcfg1","Chargino200","Stau175","LSP150",out_latex,luminosity,debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino300 Stau275 LSP250" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg1 << xsecLim_to_string(findXsecLim("cutcfg1","Chargino300", "Stau275", "LSP250", out_latex, luminosity, debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino500 Stau375 LSP350" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg1 << xsecLim_to_string(findXsecLim("cutcfg1","Chargino400", "Stau375", "LSP350", out_latex, luminosity, debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino500 Stau475 LSP450" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg1 << xsecLim_to_string(findXsecLim("cutcfg1","Chargino500", "Stau475", "LSP450", out_latex, luminosity, debug)) <<endl;

  if (debug == kTRUE) cout << "Signal scenario: Chargino100 Stau075 LSP050" <<endl;
	out_xsecmin_lspcomp_stauclose_cutcfg2 << xsecLim_to_string(findXsecLim("cutcfg2","Chargino100","Stau075","LSP050",out_latex,luminosity,debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino200 Stau175 LSP150" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg2 << xsecLim_to_string(findXsecLim("cutcfg2","Chargino200","Stau175","LSP150",out_latex,luminosity,debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino300 Stau275 LSP250" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg2 << xsecLim_to_string(findXsecLim("cutcfg2","Chargino300", "Stau275", "LSP250", out_latex, luminosity, debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino500 Stau375 LSP350" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg2 << xsecLim_to_string(findXsecLim("cutcfg2","Chargino400", "Stau375", "LSP350", out_latex, luminosity, debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino500 Stau475 LSP450" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg2 << xsecLim_to_string(findXsecLim("cutcfg2","Chargino500", "Stau475", "LSP450", out_latex, luminosity, debug)) <<endl;

  if (debug == kTRUE) cout << "Signal scenario: Chargino100 Stau075 LSP050" <<endl;
	out_xsecmin_lspcomp_stauclose_cutcfg3 << xsecLim_to_string(findXsecLim("cutcfg3","Chargino100","Stau075","LSP050",out_latex,luminosity,debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino200 Stau175 LSP150" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg3 << xsecLim_to_string(findXsecLim("cutcfg3","Chargino200","Stau175","LSP150",out_latex,luminosity,debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino300 Stau275 LSP250" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg3 << xsecLim_to_string(findXsecLim("cutcfg3","Chargino300", "Stau275", "LSP250", out_latex, luminosity, debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino500 Stau375 LSP350" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg3 << xsecLim_to_string(findXsecLim("cutcfg3","Chargino400", "Stau375", "LSP350", out_latex, luminosity, debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino500 Stau475 LSP450" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg3 << xsecLim_to_string(findXsecLim("cutcfg3","Chargino500", "Stau475", "LSP450", out_latex, luminosity, debug)) <<endl;

  if (debug == kTRUE) cout << "Signal scenario: Chargino100 Stau075 LSP050" <<endl;
	out_xsecmin_lspcomp_stauclose_cutcfg4 << xsecLim_to_string(findXsecLim("cutcfg4","Chargino100","Stau075","LSP050",out_latex,luminosity,debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino200 Stau175 LSP150" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg4 << xsecLim_to_string(findXsecLim("cutcfg4","Chargino200","Stau175","LSP150",out_latex,luminosity,debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino300 Stau275 LSP250" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg4 << xsecLim_to_string(findXsecLim("cutcfg4","Chargino300", "Stau275", "LSP250", out_latex, luminosity, debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino500 Stau375 LSP350" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg4 << xsecLim_to_string(findXsecLim("cutcfg4","Chargino400", "Stau375", "LSP350", out_latex, luminosity, debug)) <<endl;
  if (debug == kTRUE) cout << "Signal scenario: Chargino500 Stau475 LSP450" <<endl;
  out_xsecmin_lspcomp_stauclose_cutcfg4 << xsecLim_to_string(findXsecLim("cutcfg4","Chargino500", "Stau475", "LSP450", out_latex, luminosity, debug)) <<endl;

}

void makeComparisonPlot() {

	//Arrays declarations

	Int_t n1 = 17;
	Int_t n2 = 5;
	Double_t cmsxsec_x[n1], cmsxsec_y[n1], cmsxsec_x_err[n1] ,cmsxsec_y_err[n1];
	Double_t vbfxsec_x[n2], vbfxsec_y[n2], vbfxsec_x_err[n2] ,vbfxsec_y_err[n2];
  Double_t vbfxsec_cutcfg1_y[n2], vbfxsec_cutcfg2_y[n2], vbfxsec_cutcfg3_y[n2], vbfxsec_cutcfg4_y[n2];
  Double_t vbfxsec_cutcfg1_y_err[n2], vbfxsec_cutcfg2_y_err[n2], vbfxsec_cutcfg3_y_err[n2], vbfxsec_cutcfg4_y_err[n2];

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

	vector<double> v_value;
	vector<double> v_error;

  xsec_reader_values("lspcomp_stauclose", n2, v_value, v_error);

	//converter from vector<double> to Double_t array
	for (int i = 0; i < n2; i++) {
		vbfxsec_y[i] = v_value[i];
		vbfxsec_y_err[i] = v_error[i];
	}

  v_value.clear();
  v_error.clear();

	xsec_reader_values("lspcomp_stauclose_cutcfg1", n2, v_value, v_error);

	//converter from vector<double> to Double_t array
	for (int i = 0; i < n2; i++) {
		vbfxsec_cutcfg1_y[i] = v_value[i];
		vbfxsec_cutcfg1_y_err[i] = v_error[i];
	}

  v_value.clear();
  v_error.clear();

  xsec_reader_values("lspcomp_stauclose_cutcfg2", n2, v_value, v_error);

  //converter from vector<double> to Double_t array
  for (int i = 0; i < n2; i++) {
    vbfxsec_cutcfg2_y[i] = v_value[i];
    vbfxsec_cutcfg2_y_err[i] = v_error[i];
  }

  v_value.clear();
  v_error.clear();

  xsec_reader_values("lspcomp_stauclose_cutcfg3", n2, v_value, v_error);

	//converter from vector<double> to Double_t array
	for (int i = 0; i < n2; i++) {
		vbfxsec_cutcfg3_y[i] = v_value[i];
		vbfxsec_cutcfg3_y_err[i] = v_error[i];
	}

  v_value.clear();
  v_error.clear();

  xsec_reader_values("lspcomp_stauclose_cutcfg4", n2, v_value, v_error);

  //converter from vector<double> to Double_t array
  for (int i = 0; i < n2; i++) {
    vbfxsec_cutcfg4_y[i] = v_value[i];
    vbfxsec_cutcfg4_y_err[i] = v_error[i];
  }

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
	cmsxsec_x[13] = 425;
	cmsxsec_x[14] = 450;
	cmsxsec_x[15] = 475;
	cmsxsec_x[16] = 500;
	// cmsxsec_x[17] = 525;
	// cmsxsec_x[18] = 550;

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
	cmsxsec_x_err[13] = 0.;
	cmsxsec_x_err[14] = 0.;
	cmsxsec_x_err[15] = 0.;
	cmsxsec_x_err[16] = 0.;
	// cmsxsec_x_err[17] = 0.;
	// cmsxsec_x_err[18] = 0.;

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
	cmsxsec_y[13] = 93.771;
	cmsxsec_y[14] = 73.4361;
	cmsxsec_y[15] = 58.0811;
	cmsxsec_y[16] = 46.3533;
	// cmsxsec_y[17] = 37.2636;
	// cmsxsec_y[18] = 30.1656;

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
	cmsxsec_y_err[13] = 7.73547;
	cmsxsec_y_err[14] = 6.2389;
	cmsxsec_y_err[15] = 5.05005;
	cmsxsec_y_err[16] = 4.16461;
	// cmsxsec_y_err[17] = 3.46763;
	// cmsxsec_y_err[18] = 2.88065;

	//converting from
	for(Int_t i = 0; i < n1; i++) cmsxsec_y[i] *= 0.001;
	for(Int_t i = 0; i < n1; i++) cmsxsec_y_err[i] *= 0.001;


	//Set plot boundaries

	double y_max_temp = -FLT_MAX;
	double y_min_temp = FLT_MAX;

	for(Int_t i = 0; i < n1; i++) {
		if (cmsxsec_y[i] > y_max_temp) y_max_temp = cmsxsec_y[i];
		if (cmsxsec_y[i] < y_min_temp) y_min_temp = cmsxsec_y[i];
	}

	for(Int_t i = 0; i < n2; i++) {
		if (vbfxsec_y[i] > y_max_temp) y_max_temp = vbfxsec_y[i];
		if (vbfxsec_y[i] < y_min_temp) y_min_temp = vbfxsec_y[i];
    if (vbfxsec_cutcfg1_y[i] > y_max_temp) y_max_temp = vbfxsec_cutcfg1_y[i];
    if (vbfxsec_cutcfg1_y[i] < y_min_temp) y_min_temp = vbfxsec_cutcfg1_y[i];
    if (vbfxsec_cutcfg2_y[i] > y_max_temp) y_max_temp = vbfxsec_cutcfg2_y[i];
    if (vbfxsec_cutcfg2_y[i] < y_min_temp) y_min_temp = vbfxsec_cutcfg2_y[i];
    if (vbfxsec_cutcfg3_y[i] > y_max_temp) y_max_temp = vbfxsec_cutcfg3_y[i];
    if (vbfxsec_cutcfg3_y[i] < y_min_temp) y_min_temp = vbfxsec_cutcfg3_y[i];
    if (vbfxsec_cutcfg4_y[i] > y_max_temp) y_max_temp = vbfxsec_cutcfg4_y[i];
    if (vbfxsec_cutcfg4_y[i] < y_min_temp) y_min_temp = vbfxsec_cutcfg4_y[i];
	}

	double y_max = y_max_temp * 1.5;
	double y_min =  y_min_temp * 0.5;
	double x_max = 550.;
	double x_min = 0.;

	//TGraph definition
	TGraph *gr1 = new TGraphErrors (n1, cmsxsec_x, cmsxsec_y, cmsxsec_x_err, cmsxsec_y_err);
	gr1->SetMarkerColor(4);
	gr1->SetMarkerStyle(21);

	TGraph *gr2 = new TGraphErrors (n2, vbfxsec_x, vbfxsec_y, vbfxsec_x_err, vbfxsec_y_err);
	gr2->SetMarkerColor(6);
	gr2->SetMarkerStyle(21);

  TGraph *gr3 = new TGraphErrors (n2, vbfxsec_x, vbfxsec_cutcfg1_y, vbfxsec_x_err, vbfxsec_cutcfg1_y_err);
  gr3->SetMarkerColor(7);
  gr3->SetMarkerStyle(21);

  TGraph *gr4 = new TGraphErrors (n2, vbfxsec_x, vbfxsec_cutcfg2_y, vbfxsec_x_err, vbfxsec_cutcfg2_y_err);
  gr4->SetMarkerColor(8);
  gr4->SetMarkerStyle(21);

  TGraph *gr5 = new TGraphErrors (n2, vbfxsec_x, vbfxsec_cutcfg3_y, vbfxsec_x_err, vbfxsec_cutcfg3_y_err);
  gr5->SetMarkerColor(11);
  gr5->SetMarkerStyle(21);

  TGraph *gr6 = new TGraphErrors (n2, vbfxsec_x, vbfxsec_cutcfg4_y, vbfxsec_x_err, vbfxsec_cutcfg4_y_err);
  gr6->SetMarkerColor(12);
  gr6->SetMarkerStyle(21);


	//TH1F definition
	TH1F* lim_comparison_bkg = new TH1F ("lim_comparison","lim_comparison", 22, 0. , 550.);
	lim_comparison_bkg->SetTitle("CMS Work");
	lim_comparison_bkg->GetYaxis()->SetTitle("#sigma [pb]");
	lim_comparison_bkg->GetXaxis()->SetTitle("m(#tilde{#chi}^{#pm}_{1}) = m(#tilde{#chi}^{0}_{2}) [GeV]");
	lim_comparison_bkg->GetYaxis()->SetTitleOffset(1.40);
	lim_comparison_bkg->GetYaxis()->SetRangeUser(y_min,y_max);
	lim_comparison_bkg->GetXaxis()->SetRangeUser(0.,550.);
	lim_comparison_bkg->SetStats(0);

	TH1F* lim_comparison = new TH1F ("lim_comparison","lim_comparison", 22, x_min , x_max );
	lim_comparison->SetStats(0);
	lim_comparison->SetMarkerStyle(21);
	lim_comparison->SetMarkerColor(6);


	// create a multigraph and draw it
	TMultiGraph  *mg  = new TMultiGraph();
	mg->Add(gr1);
  mg->Add(gr2);
  mg->Add(gr3);
  mg->Add(gr4);
  mg->Add(gr5);
  mg->Add(gr6);


	//Converting scenario into legend
	string stau_label;
	string lsp_label;

	lsp_label = "m(#tilde{#chi}^{0}_{1}) = m(#tilde{#chi}^{#pm}_{1}) - 50 GeV";
	stau_label = "m(#tilde{#tau}^{#pm}) = m(#tilde{#chi}^{#pm}_{1})/2 + m(#tilde{#chi}^{0}_{1})/2";

	//defining legend
	TLegend* leg_1 = new TLegend(0.53,0.63,0.86,0.9);
	leg_1->SetTextSize(0.03);
	leg_1->AddEntry(gr1, "#sigma^{CMS}", "LP");
  leg_1->AddEntry(gr2, "#sigma_{lim}^{VBF}", "Pe");
  leg_1->AddEntry(gr3, "#sigma_{lim}^{VBF}(cut set 1)", "Pe");
  leg_1->AddEntry(gr4, "#sigma_{lim}^{VBF}(cut set 2)", "Pe");
  leg_1->AddEntry(gr5, "#sigma_{lim}^{VBF}(cut set 3)", "Pe");
  leg_1->AddEntry(gr6, "#sigma_{lim}^{VBF}(cut set 4)", "Pe");


	//defining legend
  TLegend* leg_2 = new TLegend(0.23,0.13,0.56,0.23);
	leg_2->SetTextSize(0.015);
	leg_2->AddEntry((TObject*)0, (stau_label).c_str(), "");
	leg_2->AddEntry((TObject*)0, (lsp_label).c_str(), "");

	//Canvas definition
	TCanvas *my_canvas = new TCanvas("mycanvas","mycanvas",600.,600.);
	my_canvas->cd();
	gPad->SetLogy();

	//Drawing workaround
	lim_comparison_bkg->Draw();
	lim_comparison->Draw("Pe+same");
  mg->Draw("LP");
	leg_1->Draw("same");
	leg_2->Draw("same");
	my_canvas->Update();
	my_canvas->Print("results/out_xsecmin_cutvar.pdf");

}
