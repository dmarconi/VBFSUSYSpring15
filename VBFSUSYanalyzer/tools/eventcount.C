#include <TROOT.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double vbfefficiency(double evenCRcounts, double oddCRcounts, double evenCRnonQCDbg, double oddCRnonQCDbg){

	        return ( (oddCRcounts - oddCRnonQCDbg) / ( (oddCRcounts - oddCRnonQCDbg) + (evenCRcounts - evenCRnonQCDbg) ));
}

double vbfefficiency_staterr(double evenCRcounts, double evenCRcounts_err,double oddCRcounts, double oddCRcounts_err, double evenCRnonQCDbg, double evenCRnonQCDbg_err, double oddCRnonQCDbg, double oddCRnonQCDbg_err) {
	double temp_evenCRcounts = evenCRcounts - evenCRnonQCDbg;
	double temp_evenCRcounts_err = sqrt( pow(evenCRcounts_err,2.) + pow(evenCRnonQCDbg_err,2.)  );
	double temp_oddCRcounts = oddCRcounts - oddCRnonQCDbg;
	double temp_oddCRcounts_err = sqrt( pow(oddCRcounts_err,2.) + pow(oddCRnonQCDbg_err,2.)  );
	return sqrt( pow(  (( temp_evenCRcounts * temp_oddCRcounts_err)/(pow((temp_oddCRcounts + temp_evenCRcounts),2.))), 2.)  + pow(  ((temp_oddCRcounts * temp_evenCRcounts_err)/(pow((temp_oddCRcounts + evenCRcounts),2.))), 2.)    );
}

double binerrors(TH1F* h1, int minMETbin, int maxMETbin){
	double temperror = 0;
	for (int i = minMETbin; i < maxMETbin + 1; i++){
		temperror = temperror + pow(h1->GetBinError(i),2.);
	}
	return sqrt(temperror);
}

void eventcount(string plotname) {

	TFile *inputfile = TFile::Open("allQCD.root");

	TH1F* h1_TauAnyIso;
	TH1F* h1_TauAnyIsoVBFInverted;
	TH1F* h1_TauAnyIsoPlusNones;
	TH1F* h1_TauAnyIsoPlusNonesVBFInverted;
	TH1F* h1_TauAntiMediumIso;
	TH1F* h1_Tau2LooseIso;
	TH1F* h1_TauAntiMediumIsoVBFInverted;
	TH1F* h1_Tau2LooseIsoVBFInverted;
	TH1F* h1_TauAntiTightIso;
	TH1F* h1_Tau2MediumIso;
	TH1F* h1_TauAntiTightIsoVBFInverted;
	TH1F* h1_Tau2MediumIsoVBFInverted;
	TH1F* h1_Tau1TightIso;
	TH1F* h1_Tau1TightIsoVBFInverted;
	TH1F* h1_Tau2TightIso;
	TH1F* h1_Tau2TightIsoVBFInverted;

	cout << "FLAG1" << endl;
	h1_TauAnyIso = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoObjectSelection/" + plotname).c_str())));
	h1_TauAnyIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	cout << "FLAG2" << endl;
	h1_TauAnyIsoPlusNones = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoPlusNonesObjectSelection/" + plotname).c_str())));
	h1_TauAnyIsoPlusNonesVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoPlusNonesVBFInvertedObjectSelection/" + plotname).c_str())));
	cout << "FLAG3" << endl;
	h1_Tau2LooseIso = ((TH1F*)(inputfile->Get(("demo/Tau2LooseIsoObjectSelection/" + plotname).c_str())));
	h1_TauAntiMediumIso = ((TH1F*)(inputfile->Get(("demo/TauAntiMediumIsoObjectSelection/" + plotname).c_str())));
	cout << "FLAG4" << endl;
	h1_TauAntiMediumIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAntiMediumIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_Tau2LooseIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau2LooseIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_TauAntiTightIso = ((TH1F*)(inputfile->Get(("demo/TauAntiTightIsoObjectSelection/" + plotname).c_str())));
	cout << "FLAG5" << endl;
	h1_Tau2MediumIso = ((TH1F*)(inputfile->Get(("demo/Tau2MediumIsoObjectSelection/" + plotname).c_str())));
	h1_TauAntiTightIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAntiTightIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	cout << "FLAG6" << endl;
	h1_Tau2MediumIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau2MediumIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	cout << "FLAG7" << endl;
	h1_Tau1TightIso = ((TH1F*)(inputfile->Get(("demo/Tau1TightIsoObjectSelection/" + plotname).c_str())));
	h1_Tau1TightIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau1TightIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	cout << "FLAG8" << endl;
	h1_Tau2TightIso = ((TH1F*)(inputfile->Get(("demo/Taui2TightIsoObjectSelection/" + plotname).c_str())));
	cout << "FLAG9" << endl;
	h1_Tau2TightIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau2TightIsoVBFInvertedObjectSelection/" + plotname).c_str())));


	cout << "FLAG1" << endl;
	TH1F* h1_AnyIsocounts;
	TH1F* h1_AnyIsocountsVBFInverted;

	double TauAnyIsocounts = h1_TauAnyIso->GetBinContent(3);
	double TauAnyIsoVBFInvertedcounts = h1_TauAnyIsoVBFInverted->GetBinContent(3);
	double TauAnyIsoPlusNonescounts = h1_TauAnyIsoPlusNones->GetBinContent(3);
	cout << "FLAG2" << endl;
	double TauAnyIsoPlusNonesVBFInvertedcounts = h1_TauAnyIsoPlusNonesVBFInverted->GetBinContent(3);
	double TauAntiMediumIsocounts = h1_TauAntiMediumIso->GetBinContent(3);
	double Tau2LooseIsocounts = h1_Tau2LooseIso->GetBinContent(3);
	double TauAntiMediumIsoVBFInvertedcounts = h1_TauAntiMediumIsoVBFInverted->GetBinContent(3);
	cout << "FLAG3" << endl;
	double Tau2LooseIsoVBFInvertedcounts = h1_Tau2LooseIsoVBFInverted->GetBinContent(3);
	double TauAntiTightIsocounts = h1_TauAntiTightIso->GetBinContent(3);
	double Tau2MediumIsocounts = h1_Tau2MediumIso->GetBinContent(3);
	cout << "FLAG4" << endl;
	double TauAntiTightIsoVBFInvertedcounts = h1_TauAntiTightIsoVBFInverted->GetBinContent(3);
	double Tau2MediumIsoVBFInvertedcounts = h1_Tau2MediumIsoVBFInverted->GetBinContent(3);
	double Tau1TightIsocounts = h1_Tau1TightIso->GetBinContent(3);
	cout << "FLAG5" << endl;
	double Tau1TightIsoVBFInvertedcounts = h1_Tau1TightIsoVBFInverted->GetBinContent(3);
	cout << "FLAG5" << endl;
	double Tau2TightIsocounts = h1_Tau2TightIso->GetBinContent(3);
	cout << "FLAG5" << endl;
	double Tau2TightIsoVBFInvertedcounts = h1_Tau2TightIsoVBFInverted->GetBinContent(3);

	cout << "FLAG6" << endl;
	double TauAnyIsocounts_err = binerrors(h1_TauAnyIso, 2, 4);
	double TauAnyIsoVBFInvertedcounts_err = binerrors(h1_TauAnyIsoVBFInverted, 2, 4);
	double TauAnyIsoPlusNonescounts_err = binerrors(h1_TauAnyIsoPlusNones, 2, 4);
	cout << "FLAG7" << endl;
	double TauAnyIsoPlusNonesVBFInvertedcounts_err = binerrors(h1_TauAnyIsoPlusNonesVBFInverted, 2, 4);
	double TauAntiMediumIsocounts_err = 	binerrors(h1_TauAntiMediumIso, 2, 4);
	double TauAntiMediumIsoVBFInvertedcounts_err = binerrors(h1_TauAntiMediumIsoVBFInverted, 2, 4);
	double Tau2LooseIsocounts_err = binerrors(h1_Tau2LooseIso, 2, 4);
	cout << "FLAG8" << endl;
	double Tau2LooseIsoVBFInvertedcounts_err = binerrors(h1_Tau2LooseIsoVBFInverted, 2, 4);
	double TauAntiTightIsocounts_err =	binerrors(h1_TauAntiTightIso, 2, 4);
	double TauAntiTightIsoVBFInvertedcounts_err = binerrors(h1_TauAntiTightIsoVBFInverted, 2, 4);
	cout << "FLAG9" << endl;
	double Tau2MediumIsocounts_err = binerrors(h1_Tau2MediumIso, 2, 4);
	double Tau2MediumIsoVBFInvertedcounts_err = binerrors(h1_Tau2MediumIsoVBFInverted, 2, 4);
	double Tau1TightIsocounts_err = binerrors(h1_Tau1TightIso, 2, 4);
	double Tau1TightIsoVBFInvertedcounts_err = binerrors(h1_Tau1TightIsoVBFInverted, 2, 4);
	cout << "FLAG0" << endl;
	double Tau2TightIsocounts_err = binerrors(h1_Tau2TightIso, 2, 4);
	double Tau2TightIsoVBFInvertedcounts_err = binerrors(h1_Tau2TightIsoVBFInverted, 2, 4);
	
	cout << "FLAG1" << endl;
	cout << endl;
	cout << endl;
	cout << "-------EVENT COUNT--------" << endl;
	cout << "Tau2TightIso: " << Tau2TightIsocounts << " +/- " << Tau2TightIsocounts_err << " ( % " << (Tau2TightIsocounts_err/Tau2TightIsocounts)*100. << " )"<< endl;
	cout << "Tau2TightIsoVBFInverted: " << Tau2TightIsoVBFInvertedcounts << " +/- " << Tau2TightIsoVBFInvertedcounts_err << " ( % " << (Tau2TightIsoVBFInvertedcounts_err/Tau2TightIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "Tau2MediumIso: " << Tau2MediumIsocounts << " +/- " << Tau2MediumIsocounts_err << " ( % " << (Tau2MediumIsocounts_err/Tau2MediumIsocounts)*100. << " )"<< endl;
	cout << "Tau2MediumIsoVBFInverted: " << Tau2MediumIsoVBFInvertedcounts << " +/- " << Tau2MediumIsoVBFInvertedcounts_err << " ( % " << (Tau2MediumIsoVBFInvertedcounts_err/Tau2MediumIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "TauAntiTightIso: " << TauAntiTightIsocounts << " +/- " << TauAntiTightIsocounts_err << " ( % " << (TauAntiTightIsocounts_err/TauAntiTightIsocounts)*100. << " )"<< endl;
	cout << "TauAntiTightIsoVBFInverted: " << TauAntiTightIsoVBFInvertedcounts << " +/- " << TauAntiTightIsoVBFInvertedcounts_err << " ( % " << (TauAntiTightIsoVBFInvertedcounts_err/TauAntiTightIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "Tau2LooseIso: " << Tau2LooseIsocounts << " +/- " << Tau2LooseIsocounts_err << " ( % " << (Tau2LooseIsocounts_err/Tau2LooseIsocounts)*100. << " )"<< endl;
	cout << "Tau2LooseIsoVBFInverted: " << Tau2LooseIsoVBFInvertedcounts << " +/- " << Tau2LooseIsoVBFInvertedcounts_err << " ( % " << (Tau2LooseIsoVBFInvertedcounts_err/Tau2LooseIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "TauAntiMediumIso: " << TauAntiMediumIsocounts << " +/- " << TauAntiMediumIsocounts_err << " ( % " << (TauAntiMediumIsocounts_err/TauAntiMediumIsocounts)*100. << " )"<< endl;
	cout << "TauAntiMediumIsoVBFInverted: " << TauAntiMediumIsoVBFInvertedcounts << " +/- " << TauAntiMediumIsoVBFInvertedcounts_err << " ( % " << (TauAntiMediumIsoVBFInvertedcounts_err/TauAntiMediumIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "Tau1TightIso: " << Tau1TightIsocounts << " +/- " << Tau1TightIsocounts_err << " ( % " << (Tau1TightIsocounts_err/Tau1TightIsocounts)*100. << " )"<< endl;
	cout << "Tau1TightIsoVBFInverted: " << Tau1TightIsoVBFInvertedcounts << " +/- " << Tau1TightIsoVBFInvertedcounts_err << " ( % " << (Tau1TightIsoVBFInvertedcounts_err/Tau1TightIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "TauAnyIso: " << TauAnyIsocounts << " +/- " <<  TauAnyIsocounts_err << " ( % " << (TauAnyIsocounts_err/TauAnyIsocounts)*100. << " )"<< endl;
	cout << "TauAnyIsoVBFInverted: " << TauAnyIsoVBFInvertedcounts << " +/- " << TauAnyIsoVBFInvertedcounts_err << " ( % " << (TauAnyIsoVBFInvertedcounts_err/TauAnyIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "TauAnyIsoPlusNones: " << TauAnyIsoPlusNonescounts << " +/- " << TauAnyIsoPlusNonescounts_err << " ( % " << (TauAnyIsoPlusNonescounts_err/TauAnyIsoPlusNonescounts)*100. << " )"<< endl;
	cout << "TauAnyIsoPlusNonesVBFInverted: " << TauAnyIsoPlusNonesVBFInvertedcounts << " +/- " << TauAnyIsoPlusNonesVBFInvertedcounts_err << " ( % " << (TauAnyIsoPlusNonesVBFInvertedcounts_err/TauAnyIsoPlusNonesVBFInvertedcounts)*100. << " )"<< endl;
	cout << endl;
	cout << endl;

	double Tau2TightIso_eff = vbfefficiency(Tau2TightIsoVBFInvertedcounts, Tau2TightIsocounts , 0., 0.);
	double Tau2MediumIso_eff = vbfefficiency(Tau2MediumIsoVBFInvertedcounts, Tau2MediumIsocounts , 0., 0.); 
	double TauAntiTightIso_eff = vbfefficiency(TauAntiTightIsoVBFInvertedcounts, TauAntiTightIsocounts , 0., 0.); 
	double TauAntiMediumIso_eff = vbfefficiency(TauAntiMediumIsoVBFInvertedcounts, TauAntiMediumIsocounts , 0., 0.);
	double Tau2LooseIso_eff = vbfefficiency(Tau2LooseIsoVBFInvertedcounts, Tau2LooseIsocounts , 0., 0.);
	double Tau1TightIso_eff = vbfefficiency(Tau1TightIsoVBFInvertedcounts, Tau1TightIsocounts , 0., 0.);
	double TauAnyIso_eff = vbfefficiency(TauAnyIsoVBFInvertedcounts, TauAnyIsocounts , 0., 0.);
	double TauAnyIsoPlusNones_eff = vbfefficiency(TauAnyIsoPlusNonesVBFInvertedcounts, TauAnyIsoPlusNonescounts , 0., 0.);

	double Tau2TightIso_eff_err = vbfefficiency_staterr(Tau2TightIsoVBFInvertedcounts, Tau2TightIsoVBFInvertedcounts_err,Tau2TightIsocounts, Tau2TightIsocounts_err, 0., 0., 0., 0.);
	double Tau2MediumIso_eff_err = vbfefficiency_staterr( Tau2MediumIsoVBFInvertedcounts, Tau2MediumIsoVBFInvertedcounts_err, Tau2MediumIsocounts, Tau2MediumIsocounts_err, 0., 0., 0., 0.);
	double TauAntiTightIso_eff_err = vbfefficiency_staterr( TauAntiTightIsoVBFInvertedcounts, TauAntiTightIsoVBFInvertedcounts_err, TauAntiTightIsocounts, TauAntiTightIsocounts_err, 0., 0., 0., 0.);
	double Tau2LooseIso_eff_err = vbfefficiency_staterr( Tau2LooseIsoVBFInvertedcounts, Tau2LooseIsoVBFInvertedcounts_err, Tau2LooseIsocounts, Tau2LooseIsocounts_err, 0., 0., 0., 0.);
	double TauAntiMediumIso_eff_err = vbfefficiency_staterr( TauAntiMediumIsoVBFInvertedcounts, TauAntiMediumIsoVBFInvertedcounts_err, TauAntiMediumIsocounts, TauAntiMediumIsocounts_err, 0., 0., 0., 0.);
	double Tau1TightIso_eff_err = vbfefficiency_staterr( Tau1TightIsoVBFInvertedcounts, Tau1TightIsoVBFInvertedcounts_err, Tau1TightIsocounts, Tau1TightIsocounts_err, 0., 0., 0., 0.);
	double TauAnyIso_eff_err = vbfefficiency_staterr( TauAnyIsoVBFInvertedcounts, TauAnyIsoVBFInvertedcounts_err, TauAnyIsocounts, TauAnyIsocounts_err, 0., 0., 0., 0.);
	double TauAnyIsoPlusNones_eff_err = vbfefficiency_staterr( TauAnyIsoPlusNonesVBFInvertedcounts, TauAnyIsoPlusNonesVBFInvertedcounts_err, TauAnyIsoPlusNonescounts, TauAnyIsoPlusNonescounts_err, 0., 0., 0., 0.);


	cout << "-------VBF EFFICIENCY--------" << endl;
	cout << "Tau2TightIso_eff: " << Tau2TightIso_eff << " +/- " << Tau2TightIso_eff_err << " ( % " << (Tau2TightIso_eff_err/Tau2TightIso_eff)*100. << " )"<< endl;
	cout << "Tau2MediumIso_eff: " << Tau2MediumIso_eff << " +/- " << Tau2MediumIso_eff_err << " ( % " << (Tau2MediumIso_eff_err/Tau2MediumIso_eff)*100. << " )"<< endl;
	cout << "TauAntiTightIso_eff: " << TauAntiTightIso_eff << " +/- " << TauAntiTightIso_eff_err << " ( % " << (TauAntiTightIso_eff_err/TauAntiTightIso_eff)*100. << " )"<< endl;
	cout << "Tau2LooseIso_eff: " << Tau2LooseIso_eff << " +/- " << Tau2LooseIso_eff_err << " ( % " << (Tau2LooseIso_eff_err/Tau2LooseIso_eff)*100. << " )"<< endl;
	cout << "TauAntiMediumIso_eff: " << TauAntiMediumIso_eff << " +/- " << TauAntiMediumIso_eff_err << " ( % " << (TauAntiMediumIso_eff_err/TauAntiMediumIso_eff)*100. << " )"<< endl;
	cout << "Tau1TightIso_eff: " << Tau1TightIso_eff << " +/- " << Tau1TightIso_eff_err << " ( % " << (Tau1TightIso_eff_err/Tau1TightIso_eff)*100. << " )"<< endl;
	cout << "TauAnyIso_eff: " << TauAnyIso_eff << " +/- " << TauAnyIso_eff_err << " ( % " << (TauAnyIso_eff_err/TauAnyIso_eff)*100. << " )"<< endl;
	cout << "TauAnyIsoPlusNones_eff: " << TauAnyIsoPlusNones_eff << " +/- " << TauAnyIsoPlusNones_eff_err << " ( % " << (TauAnyIsoPlusNones_eff_err/TauAnyIsoPlusNones_eff)*100. << " )"<< endl;

	cout << endl;

	cout << "//-------------------LATEX OUTPUT------------------------//" << endl;
	cout << "\\begin{table}[ht]" << endl;
	cout << " \\centering{ " << endl;
	cout << "\\tabcolsep=0.05cm " << endl;
	cout << "\\begin{tabular}{| l | c | c |}  " << endl;
	cout << "\\hline\\hline " << endl;
	cout << " $ \\tau isolation region $     & VBF Selection    & Inverted VBF Selection \\\\ [0.5ex] \\hline " << endl;
	cout << "$ Only 2 Tight $    &$ "<< Tau2TightIsocounts << " \\pm  "<< " \\% " << (Tau2TightIsocounts_err/Tau2TightIsocounts)*100.  << " $  &$  "
		<< Tau2TightIsoVBFInvertedcounts << "\\pm "<< "  \\% " << (Tau2TightIsoVBFInvertedcounts_err/Tau2TightIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ 1Tight (+M+L+N) $   &$ "<< Tau1TightIsocounts << " \\pm  "<< " \\% " << (Tau1TightIsocounts_err/Tau1TightIsocounts)*100. << " $  &$  "
		<< Tau1TightIsoVBFInvertedcounts << "\\pm "<< "  \\% " << (Tau1TightIsoVBFInvertedcounts_err/Tau1TightIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ Only 2Medium $    &$ "<< Tau2MediumIsocounts << " \\pm  "<< " \\% " << (Tau2MediumIsocounts_err/Tau2MediumIsocounts)*100. << " $  &$  "
		<< Tau2MediumIsoVBFInvertedcounts << "\\pm "<< " \\% " << (Tau2MediumIsoVBFInvertedcounts_err/Tau2MediumIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ AntiTight (M+L+N) $    &$ "<< TauAntiTightIsocounts << " \\pm  "<< " \\% " << (TauAntiTightIsocounts_err/TauAntiTightIsocounts)*100. << " $  &$  "
		<< TauAntiTightIsoVBFInvertedcounts << "\\pm "<< " \\% " << (TauAntiTightIsoVBFInvertedcounts_err/TauAntiTightIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ Only 2 Loose $    &$ "<< Tau2LooseIsocounts << " \\pm  "<< " \\% " << (Tau2LooseIsocounts_err/Tau2LooseIsocounts)*100. << " $  &$  "
		<< Tau2LooseIsoVBFInvertedcounts << "\\pm "<< " \\% " << (Tau2LooseIsoVBFInvertedcounts_err/Tau2LooseIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ AntiMedium (L+N) $    &$ "<< TauAntiMediumIsocounts << " \\pm  "<< " \\% " << (TauAntiMediumIsocounts_err/TauAntiMediumIsocounts)*100. << " $  &$  "
		<< TauAntiMediumIsoVBFInvertedcounts << "\\pm "<< " \\% " << (TauAntiMediumIsoVBFInvertedcounts_err/TauAntiMediumIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ Any iso (+N) $    &$ "<< TauAnyIsoPlusNonescounts << " \\pm  "<< " \\% " << (TauAnyIsoPlusNonescounts_err/TauAnyIsoPlusNonescounts)*100. << " $  &$  "
		<< TauAnyIsoPlusNonesVBFInvertedcounts << "\\pm "<< " \\% " << (TauAnyIsoPlusNonesVBFInvertedcounts_err/TauAnyIsoPlusNonesVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "\\hline\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "} " << endl;
	cout << "\\end{table}" << endl;
	cout << endl;
	cout << endl;
	cout << endl;
	cout << "\\begin{table}[ht]" << endl;
	cout << " \\centering{ " << endl;
	cout << "\\tabcolsep=0.05cm " << endl;
	cout << "\\begin{tabular}{| l | c | }  " << endl;
	cout << "\\hline\\hline " << endl;
//	cout << " Region   & $\\epsilon_{VBF} $" << endl;
//TODO end the latex script 
	cout << " Variable     & Only 2 Tight region     & 1Tight region (+M+L+N)     & Medium (+L+N)       & Loose (+N)   & AnyIso (+N) \\\\ [0.5ex] \\hline " << endl;
	cout << "$\\epsilon^{QCD}_{VBF}$    &$ "<< Tau2TightIso_eff << " \\pm  "<< " \\% " << (Tau2TightIso_eff_err/Tau2TightIso_eff)*100. << " $  &$  "
		<< Tau1TightIso_eff << "\\pm "<< " \\% " << (Tau1TightIso_eff_err/Tau1TightIso_eff)*100. << " $  &$  "
		<< TauAntiTightIso_eff << "\\pm "<< " \\% " << (TauAntiTightIso_eff_err/TauAntiTightIso_eff)*100. << " $  &$  "
		<< TauAntiMediumIso_eff << "\\pm "<< " \\% " << (TauAntiMediumIso_eff_err/TauAntiMediumIso_eff)*100. << " $  &$  "
		<< TauAnyIsoPlusNones_eff << "\\pm "<< " \\% " << (TauAnyIsoPlusNones_eff_err/TauAnyIsoPlusNones_eff)*100. << " $ \\\\ " << endl;
	cout << "\\hline\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "} " << endl;
	cout << "\\end{table}" << endl;

}
