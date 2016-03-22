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
	TH1F* h1_TauLooseIso;
	TH1F* h1_TauLooseIsoVBFInverted;
	TH1F* h1_TauMediumIso;
	TH1F* h1_TauMediumIsoVBFInverted;
	TH1F* h1_Tau1TightIso;
	TH1F* h1_Tau1TightIsoVBFInverted;
	TH1F* h1_TauTightIso;
	TH1F* h1_TauTightIsoVBFInverted;

	h1_TauAnyIso = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoObjectSelection/" + plotname).c_str())));
	h1_TauAnyIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_TauAnyIsoPlusNones = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoPlusNonesObjectSelection/" + plotname).c_str())));
	h1_TauAnyIsoPlusNonesVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoPlusNonesVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_TauLooseIso = ((TH1F*)(inputfile->Get(("demo/TauLooseIsoObjectSelection/" + plotname).c_str())));
	h1_TauLooseIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauLooseIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_TauMediumIso = ((TH1F*)(inputfile->Get(("demo/TauMediumIsoObjectSelection/" + plotname).c_str())));
	h1_TauMediumIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauMediumIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_Tau1TightIso = ((TH1F*)(inputfile->Get(("demo/Tau1TightIsoObjectSelection/" + plotname).c_str())));
	h1_Tau1TightIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau1TightIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_TauTightIso = ((TH1F*)(inputfile->Get(("demo/TauTightIsoObjectSelection/" + plotname).c_str())));
	h1_TauTightIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauTightIsoVBFInvertedObjectSelection/" + plotname).c_str())));


	TH1F* h1_AnyIsocounts;
	TH1F* h1_AnyIsocountsVBFInverted;

	double TauAnyIsocounts = h1_TauAnyIso->GetBinContent(3);
	double TauAnyIsoVBFInvertedcounts = h1_TauAnyIsoVBFInverted->GetBinContent(3);
	double TauAnyIsoPlusNonescounts = h1_TauAnyIsoPlusNones->GetBinContent(3);
	double TauAnyIsoPlusNonesVBFInvertedcounts = h1_TauAnyIsoPlusNonesVBFInverted->GetBinContent(3);
	double TauLooseIsocounts = h1_TauLooseIso->GetBinContent(3);
	double TauLooseIsoVBFInvertedcounts = h1_TauLooseIsoVBFInverted->GetBinContent(3);
	double TauMediumIsocounts = h1_TauMediumIso->GetBinContent(3);
	double TauMediumIsoVBFInvertedcounts = h1_TauMediumIsoVBFInverted->GetBinContent(3);
	double Tau1TightIsocounts = h1_Tau1TightIso->GetBinContent(3);
	double Tau1TightIsoVBFInvertedcounts = h1_Tau1TightIsoVBFInverted->GetBinContent(3);
	double TauTightIsocounts = h1_TauTightIso->GetBinContent(3);
	double TauTightIsoVBFInvertedcounts = h1_TauTightIsoVBFInverted->GetBinContent(3);

	double TauAnyIsocounts_err = binerrors(h1_TauAnyIso, 2, 4);
	double TauAnyIsoVBFInvertedcounts_err = binerrors(h1_TauAnyIsoVBFInverted, 2, 4);
	double TauAnyIsoPlusNonescounts_err = binerrors(h1_TauAnyIsoPlusNones, 2, 4);
	double TauAnyIsoPlusNonesVBFInvertedcounts_err = binerrors(h1_TauAnyIsoPlusNonesVBFInverted, 2, 4);
	double TauLooseIsocounts_err = 	binerrors(h1_TauLooseIso, 2, 4);
	double TauLooseIsoVBFInvertedcounts_err = binerrors(h1_TauLooseIsoVBFInverted, 2, 4);
	double TauMediumIsocounts_err =	binerrors(h1_TauMediumIso, 2, 4);
	double TauMediumIsoVBFInvertedcounts_err = binerrors(h1_TauMediumIsoVBFInverted, 2, 4);
	double Tau1TightIsocounts_err = binerrors(h1_Tau1TightIso, 2, 4);
	double Tau1TightIsoVBFInvertedcounts_err = binerrors(h1_Tau1TightIsoVBFInverted, 2, 4);
	double TauTightIsocounts_err = binerrors(h1_TauTightIso, 2, 4);
	double TauTightIsoVBFInvertedcounts_err = binerrors(h1_TauTightIsoVBFInverted, 2, 4);
	
	cout << endl;
	cout << endl;
	cout << "-------EVENT COUNT--------" << endl;
	cout << "TauTightIso: " << TauTightIsocounts << " +/- " << TauTightIsocounts_err << " ( % " << (TauTightIsocounts_err/TauTightIsocounts)*100. << " )"<< endl;
	cout << "TauTightIsoVBFInverted: " << TauTightIsoVBFInvertedcounts << " +/- " << TauTightIsoVBFInvertedcounts_err << " ( % " << (TauTightIsoVBFInvertedcounts_err/TauTightIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "TauMediumIso: " << TauMediumIsocounts << " +/- " << TauMediumIsocounts_err << " ( % " << (TauMediumIsocounts_err/TauMediumIsocounts)*100. << " )"<< endl;
	cout << "TauMediumIsoVBFInverted: " << TauMediumIsoVBFInvertedcounts << " +/- " << TauMediumIsoVBFInvertedcounts_err << " ( % " << (TauMediumIsoVBFInvertedcounts_err/TauMediumIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "TauLooseIso: " << TauLooseIsocounts << " +/- " << TauLooseIsocounts_err << " ( % " << (TauLooseIsocounts_err/TauLooseIsocounts)*100. << " )"<< endl;
	cout << "TauLooseIsoVBFInverted: " << TauLooseIsoVBFInvertedcounts << " +/- " << TauLooseIsoVBFInvertedcounts_err << " ( % " << (TauLooseIsoVBFInvertedcounts_err/TauLooseIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "Tau1TightIso: " << Tau1TightIsocounts << " +/- " << Tau1TightIsocounts_err << " ( % " << (Tau1TightIsocounts_err/Tau1TightIsocounts)*100. << " )"<< endl;
	cout << "Tau1TightIsoVBFInverted: " << Tau1TightIsoVBFInvertedcounts << " +/- " << Tau1TightIsoVBFInvertedcounts_err << " ( % " << (Tau1TightIsoVBFInvertedcounts_err/Tau1TightIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "TauAnyIso: " << TauAnyIsocounts << " +/- " <<  TauAnyIsocounts_err << " ( % " << (TauAnyIsocounts_err/TauAnyIsocounts)*100. << " )"<< endl;
	cout << "TauAnyIsoVBFInverted: " << TauAnyIsoVBFInvertedcounts << " +/- " << TauAnyIsoVBFInvertedcounts_err << " ( % " << (TauAnyIsoVBFInvertedcounts_err/TauAnyIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "TauAnyIsoPlusNones: " << TauAnyIsoPlusNonescounts << " +/- " << TauAnyIsoPlusNonescounts_err << " ( % " << (TauAnyIsoPlusNonescounts_err/TauAnyIsoPlusNonescounts)*100. << " )"<< endl;
	cout << "TauAnyIsoPlusNonesVBFInverted: " << TauAnyIsoPlusNonesVBFInvertedcounts << " +/- " << TauAnyIsoPlusNonesVBFInvertedcounts_err << " ( % " << (TauAnyIsoPlusNonesVBFInvertedcounts_err/TauAnyIsoPlusNonesVBFInvertedcounts)*100. << " )"<< endl;
	cout << endl;
	cout << endl;

	double TauTightIso_eff = vbfefficiency(TauTightIsoVBFInvertedcounts, TauTightIsocounts , 0., 0.);
	double TauMediumIso_eff = vbfefficiency(TauMediumIsoVBFInvertedcounts, TauMediumIsocounts , 0., 0.); 
	double TauLooseIso_eff = vbfefficiency(TauLooseIsoVBFInvertedcounts, TauLooseIsocounts , 0., 0.);
	double Tau1TightIso_eff = vbfefficiency(Tau1TightIsoVBFInvertedcounts, Tau1TightIsocounts , 0., 0.);
	double TauAnyIso_eff = vbfefficiency(TauAnyIsoVBFInvertedcounts, TauAnyIsocounts , 0., 0.);
	double TauAnyIsoPlusNones_eff = vbfefficiency(TauAnyIsoPlusNonesVBFInvertedcounts, TauAnyIsoPlusNonescounts , 0., 0.);

	double TauTightIso_eff_err = vbfefficiency_staterr(TauTightIsoVBFInvertedcounts, TauTightIsoVBFInvertedcounts_err,TauTightIsocounts, TauTightIsocounts_err, 0., 0., 0., 0.);
	double TauMediumIso_eff_err = vbfefficiency_staterr( TauMediumIsoVBFInvertedcounts, TauMediumIsoVBFInvertedcounts_err, TauMediumIsocounts, TauMediumIsocounts_err, 0., 0., 0., 0.);
	double TauLooseIso_eff_err = vbfefficiency_staterr( TauLooseIsoVBFInvertedcounts, TauLooseIsoVBFInvertedcounts_err, TauLooseIsocounts, TauLooseIsocounts_err, 0., 0., 0., 0.);
	double Tau1TightIso_eff_err = vbfefficiency_staterr( Tau1TightIsoVBFInvertedcounts, Tau1TightIsoVBFInvertedcounts_err, Tau1TightIsocounts, Tau1TightIsocounts_err, 0., 0., 0., 0.);
	double TauAnyIso_eff_err = vbfefficiency_staterr( TauAnyIsoVBFInvertedcounts, TauAnyIsoVBFInvertedcounts_err, TauAnyIsocounts, TauAnyIsocounts_err, 0., 0., 0., 0.);
	double TauAnyIsoPlusNones_eff_err = vbfefficiency_staterr( TauAnyIsoPlusNonesVBFInvertedcounts, TauAnyIsoPlusNonesVBFInvertedcounts_err, TauAnyIsoPlusNonescounts, TauAnyIsoPlusNonescounts_err, 0., 0., 0., 0.);


	cout << "-------VBF EFFICIENCY--------" << endl;
	cout << "TauTightIso_eff: " << TauTightIso_eff << " +/- " << TauTightIso_eff_err << " ( % " << (TauTightIso_eff_err/TauTightIso_eff)*100. << " )"<< endl;
	cout << "TauMediumIso_eff: " << TauMediumIso_eff << " +/- " << TauMediumIso_eff_err << " ( % " << (TauMediumIso_eff_err/TauMediumIso_eff)*100. << " )"<< endl;
	cout << "TauLooseIso_eff: " << TauLooseIso_eff << " +/- " << TauLooseIso_eff_err << " ( % " << (TauLooseIso_eff_err/TauLooseIso_eff)*100. << " )"<< endl;
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
	cout << "$ Only 2 Tight $    &$ "<< TauTightIsocounts << " \\pm  "<< " \\% " << (TauTightIsocounts_err/TauTightIsocounts)*100.  << " $  &$  "
		<< TauTightIsoVBFInvertedcounts << "\\pm "<< "  \\% " << (TauTightIsoVBFInvertedcounts_err/TauTightIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ 1Tight (+M+L+N) $   &$ "<< Tau1TightIsocounts << " \\pm  "<< " \\% " << (Tau1TightIsocounts_err/Tau1TightIsocounts)*100. << " $  &$  "
		<< Tau1TightIsoVBFInvertedcounts << "\\pm "<< "  \\% " << (Tau1TightIsoVBFInvertedcounts_err/Tau1TightIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ Medium (+L+N) $    &$ "<< TauMediumIsocounts << " \\pm  "<< " \\% " << (TauMediumIsocounts_err/TauMediumIsocounts)*100. << " $  &$  "
		<< TauMediumIsoVBFInvertedcounts << "\\pm "<< " \\% " << (TauMediumIsoVBFInvertedcounts_err/TauMediumIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ Loose (+N) $    &$ "<< TauLooseIsocounts << " \\pm  "<< " \\% " << (TauLooseIsocounts_err/TauLooseIsocounts)*100. << " $  &$  "
		<< TauLooseIsoVBFInvertedcounts << "\\pm "<< " \\% " << (TauLooseIsoVBFInvertedcounts_err/TauLooseIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
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
	cout << "\\begin{tabular}{| l | c | c | c | c | c |}  " << endl;
	cout << "\\hline\\hline " << endl;
	cout << " Variable     & Only 2 Tight region     & 1Tight region (+M+L+N)     & Medium (+L+N)       & Loose (+N)   & AnyIso (+N) \\\\ [0.5ex] \\hline " << endl;
	cout << "$\\epsilon^{QCD}_{VBF}$    &$ "<< TauTightIso_eff << " \\pm  "<< " \\% " << (TauTightIso_eff_err/TauTightIso_eff)*100. << " $  &$  "
		<< Tau1TightIso_eff << "\\pm "<< " \\% " << (Tau1TightIso_eff_err/Tau1TightIso_eff)*100. << " $  &$  "
		<< TauMediumIso_eff << "\\pm "<< " \\% " << (TauMediumIso_eff_err/TauMediumIso_eff)*100. << " $  &$  "
		<< TauLooseIso_eff << "\\pm "<< " \\% " << (TauLooseIso_eff_err/TauLooseIso_eff)*100. << " $  &$  "
		<< TauAnyIsoPlusNones_eff << "\\pm "<< " \\% " << (TauAnyIsoPlusNones_eff_err/TauAnyIsoPlusNones_eff)*100. << " $ \\\\ " << endl;
	cout << "\\hline\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "} " << endl;
	cout << "\\end{table}" << endl;

}
