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
	TH1F* h1_Tau2LooseIsoExclusive;
	TH1F* h1_Tau2LooseIsoInclusive;
	TH1F* h1_TauAntiMediumIsoVBFInverted;
	TH1F* h1_Tau2LooseIsoExclusiveVBFInverted;
	TH1F* h1_Tau2LooseIsoInclusiveVBFInverted;
	TH1F* h1_TauAntiTightIso;
	TH1F* h1_Tau2MediumIsoExclusive;
	TH1F* h1_Tau2MediumIsoInclusive;
	TH1F* h1_TauAntiTightIsoVBFInverted;
	TH1F* h1_Tau2MediumIsoExclusiveVBFInverted;
	TH1F* h1_Tau2MediumIsoInclusiveVBFInverted;
	TH1F* h1_Tau1TightIso;
	TH1F* h1_Tau1TightIsoVBFInverted;
	TH1F* h1_Tau2TightIso;
	TH1F* h1_Tau2TightIsoVBFInverted;

	h1_TauAnyIso = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoObjectSelection/" + plotname).c_str())));
	h1_TauAnyIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_TauAnyIsoPlusNones = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoPlusNonesObjectSelection/" + plotname).c_str())));
	h1_TauAnyIsoPlusNonesVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAnyIsoPlusNonesVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_Tau2LooseIsoExclusive = ((TH1F*)(inputfile->Get(("demo/Tau2LooseIsoExclusiveObjectSelection/" + plotname).c_str())));
	h1_Tau2LooseIsoInclusive = ((TH1F*)(inputfile->Get(("demo/Tau2LooseIsoInclusiveObjectSelection/" + plotname).c_str())));
	h1_TauAntiMediumIso = ((TH1F*)(inputfile->Get(("demo/TauAntiMediumIsoObjectSelection/" + plotname).c_str())));
	h1_TauAntiMediumIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAntiMediumIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_Tau2LooseIsoExclusiveVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau2LooseIsoExclusiveVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_Tau2LooseIsoInclusiveVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau2LooseIsoInclusiveVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_TauAntiTightIso = ((TH1F*)(inputfile->Get(("demo/TauAntiTightIsoObjectSelection/" + plotname).c_str())));
	h1_Tau2MediumIsoExclusive = ((TH1F*)(inputfile->Get(("demo/Tau2MediumIsoExclusiveObjectSelection/" + plotname).c_str())));
	h1_Tau2MediumIsoInclusive = ((TH1F*)(inputfile->Get(("demo/Tau2MediumIsoInclusiveObjectSelection/" + plotname).c_str())));
	h1_TauAntiTightIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/TauAntiTightIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_Tau2MediumIsoExclusiveVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau2MediumIsoExclusiveVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_Tau2MediumIsoInclusiveVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau2MediumIsoInclusiveVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_Tau1TightIso = ((TH1F*)(inputfile->Get(("demo/Tau1TightIsoObjectSelection/" + plotname).c_str())));
	h1_Tau1TightIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau1TightIsoVBFInvertedObjectSelection/" + plotname).c_str())));
	h1_Tau2TightIso = ((TH1F*)(inputfile->Get(("demo/Taui2TightIsoObjectSelection/" + plotname).c_str())));
	h1_Tau2TightIsoVBFInverted = ((TH1F*)(inputfile->Get(("demo/Tau2TightIsoVBFInvertedObjectSelection/" + plotname).c_str())));


	TH1F* h1_AnyIsocounts;
	TH1F* h1_AnyIsocountsVBFInverted;

	double TauAnyIsocounts = h1_TauAnyIso->GetBinContent(3);
	double TauAnyIsoVBFInvertedcounts = h1_TauAnyIsoVBFInverted->GetBinContent(3);
	double TauAnyIsoPlusNonescounts = h1_TauAnyIsoPlusNones->GetBinContent(3);
	double TauAnyIsoPlusNonesVBFInvertedcounts = h1_TauAnyIsoPlusNonesVBFInverted->GetBinContent(3);
	double TauAntiMediumIsocounts = h1_TauAntiMediumIso->GetBinContent(3);
	double Tau2LooseIsoExclusivecounts = h1_Tau2LooseIsoExclusive->GetBinContent(3);
	double Tau2LooseIsoInclusivecounts = h1_Tau2LooseIsoInclusive->GetBinContent(3);
	double TauAntiMediumIsoVBFInvertedcounts = h1_TauAntiMediumIsoVBFInverted->GetBinContent(3);
	double Tau2LooseIsoExclusiveVBFInvertedcounts = h1_Tau2LooseIsoExclusiveVBFInverted->GetBinContent(3);
	double Tau2LooseIsoInclusiveVBFInvertedcounts = h1_Tau2LooseIsoInclusiveVBFInverted->GetBinContent(3);
	double TauAntiTightIsocounts = h1_TauAntiTightIso->GetBinContent(3);
	double Tau2MediumIsoExclusivecounts = h1_Tau2MediumIsoExclusive->GetBinContent(3);
	double Tau2MediumIsoInclusivecounts = h1_Tau2MediumIsoInclusive->GetBinContent(3);
	double TauAntiTightIsoVBFInvertedcounts = h1_TauAntiTightIsoVBFInverted->GetBinContent(3);
	double Tau2MediumIsoExclusiveVBFInvertedcounts = h1_Tau2MediumIsoExclusiveVBFInverted->GetBinContent(3);
	double Tau2MediumIsoInclusiveVBFInvertedcounts = h1_Tau2MediumIsoInclusiveVBFInverted->GetBinContent(3);
	double Tau1TightIsocounts = h1_Tau1TightIso->GetBinContent(3);
	double Tau1TightIsoVBFInvertedcounts = h1_Tau1TightIsoVBFInverted->GetBinContent(3);
	double Tau2TightIsocounts = h1_Tau2TightIso->GetBinContent(3);
	double Tau2TightIsoVBFInvertedcounts = h1_Tau2TightIsoVBFInverted->GetBinContent(3);

	double TauAnyIsocounts_err = binerrors(h1_TauAnyIso, 2, 4);
	double TauAnyIsoVBFInvertedcounts_err = binerrors(h1_TauAnyIsoVBFInverted, 2, 4);
	double TauAnyIsoPlusNonescounts_err = binerrors(h1_TauAnyIsoPlusNones, 2, 4);
	double TauAnyIsoPlusNonesVBFInvertedcounts_err = binerrors(h1_TauAnyIsoPlusNonesVBFInverted, 2, 4);
	double TauAntiMediumIsocounts_err = 	binerrors(h1_TauAntiMediumIso, 2, 4);
	double TauAntiMediumIsoVBFInvertedcounts_err = binerrors(h1_TauAntiMediumIsoVBFInverted, 2, 4);
	double Tau2LooseIsoExclusivecounts_err = binerrors(h1_Tau2LooseIsoExclusive, 2, 4);
	double Tau2LooseIsoInclusivecounts_err = binerrors(h1_Tau2LooseIsoInclusive, 2, 4);
	double Tau2LooseIsoExclusiveVBFInvertedcounts_err = binerrors(h1_Tau2LooseIsoExclusiveVBFInverted, 2, 4);
	double Tau2LooseIsoInclusiveVBFInvertedcounts_err = binerrors(h1_Tau2LooseIsoInclusiveVBFInverted, 2, 4);
	double TauAntiTightIsocounts_err =	binerrors(h1_TauAntiTightIso, 2, 4);
	double TauAntiTightIsoVBFInvertedcounts_err = binerrors(h1_TauAntiTightIsoVBFInverted, 2, 4);
	double Tau2MediumIsoExclusivecounts_err = binerrors(h1_Tau2MediumIsoExclusive, 2, 4);
	double Tau2MediumIsoInclusivecounts_err = binerrors(h1_Tau2MediumIsoInclusive, 2, 4);
	double Tau2MediumIsoExclusiveVBFInvertedcounts_err = binerrors(h1_Tau2MediumIsoExclusiveVBFInverted, 2, 4);
	double Tau2MediumIsoInclusiveVBFInvertedcounts_err = binerrors(h1_Tau2MediumIsoInclusiveVBFInverted, 2, 4);
	double Tau1TightIsocounts_err = binerrors(h1_Tau1TightIso, 2, 4);
	double Tau1TightIsoVBFInvertedcounts_err = binerrors(h1_Tau1TightIsoVBFInverted, 2, 4);
	double Tau2TightIsocounts_err = binerrors(h1_Tau2TightIso, 2, 4);
	double Tau2TightIsoVBFInvertedcounts_err = binerrors(h1_Tau2TightIsoVBFInverted, 2, 4);
	
	cout << endl;
	cout << endl;
	cout << "-------EVENT COUNT--------" << endl;
	cout << "Tau2TightIso: " << Tau2TightIsocounts << " +/- " << Tau2TightIsocounts_err << " ( % " << (Tau2TightIsocounts_err/Tau2TightIsocounts)*100. << " )"<< endl;
	cout << "Tau2TightIsoVBFInverted: " << Tau2TightIsoVBFInvertedcounts << " +/- " << Tau2TightIsoVBFInvertedcounts_err << " ( % " << (Tau2TightIsoVBFInvertedcounts_err/Tau2TightIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "Tau2MediumIsoExclusive: " << Tau2MediumIsoExclusivecounts << " +/- " << Tau2MediumIsoExclusivecounts_err << " ( % " << (Tau2MediumIsoExclusivecounts_err/Tau2MediumIsoExclusivecounts)*100. << " )"<< endl;
	cout << "Tau2MediumIsoInclusive: " << Tau2MediumIsoInclusivecounts << " +/- " << Tau2MediumIsoInclusivecounts_err << " ( % " << (Tau2MediumIsoInclusivecounts_err/Tau2MediumIsoInclusivecounts)*100. << " )"<< endl;
	cout << "Tau2MediumIsoExclusiveVBFInverted: " << Tau2MediumIsoExclusiveVBFInvertedcounts << " +/- " << Tau2MediumIsoExclusiveVBFInvertedcounts_err << " ( % " << (Tau2MediumIsoExclusiveVBFInvertedcounts_err/Tau2MediumIsoExclusiveVBFInvertedcounts)*100. << " )"<< endl;
	cout << "Tau2MediumIsoInclusiveVBFInverted: " << Tau2MediumIsoInclusiveVBFInvertedcounts << " +/- " << Tau2MediumIsoInclusiveVBFInvertedcounts_err << " ( % " << (Tau2MediumIsoInclusiveVBFInvertedcounts_err/Tau2MediumIsoInclusiveVBFInvertedcounts)*100. << " )"<< endl;
	cout << "TauAntiTightIso: " << TauAntiTightIsocounts << " +/- " << TauAntiTightIsocounts_err << " ( % " << (TauAntiTightIsocounts_err/TauAntiTightIsocounts)*100. << " )"<< endl;
	cout << "TauAntiTightIsoVBFInverted: " << TauAntiTightIsoVBFInvertedcounts << " +/- " << TauAntiTightIsoVBFInvertedcounts_err << " ( % " << (TauAntiTightIsoVBFInvertedcounts_err/TauAntiTightIsoVBFInvertedcounts)*100. << " )"<< endl;
	cout << "Tau2LooseIsoExclusive: " << Tau2LooseIsoExclusivecounts << " +/- " << Tau2LooseIsoExclusivecounts_err << " ( % " << (Tau2LooseIsoExclusivecounts_err/Tau2LooseIsoExclusivecounts)*100. << " )"<< endl;
	cout << "Tau2LooseIsoInclusive: " << Tau2LooseIsoInclusivecounts << " +/- " << Tau2LooseIsoInclusivecounts_err << " ( % " << (Tau2LooseIsoInclusivecounts_err/Tau2LooseIsoInclusivecounts)*100. << " )"<< endl;
	cout << "Tau2LooseIsoExclusiveVBFInverted: " << Tau2LooseIsoExclusiveVBFInvertedcounts << " +/- " << Tau2LooseIsoExclusiveVBFInvertedcounts_err << " ( % " << (Tau2LooseIsoExclusiveVBFInvertedcounts_err/Tau2LooseIsoExclusiveVBFInvertedcounts)*100. << " )"<< endl;
	cout << "Tau2LooseIsoInclusiveVBFInverted: " << Tau2LooseIsoInclusiveVBFInvertedcounts << " +/- " << Tau2LooseIsoInclusiveVBFInvertedcounts_err << " ( % " << (Tau2LooseIsoInclusiveVBFInvertedcounts_err/Tau2LooseIsoInclusiveVBFInvertedcounts)*100. << " )"<< endl;
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
	double Tau2MediumIsoExclusive_eff = vbfefficiency(Tau2MediumIsoExclusiveVBFInvertedcounts, Tau2MediumIsoExclusivecounts , 0., 0.); 
	double Tau2MediumIsoInclusive_eff = vbfefficiency(Tau2MediumIsoInclusiveVBFInvertedcounts, Tau2MediumIsoInclusivecounts , 0., 0.); 
	double TauAntiTightIso_eff = vbfefficiency(TauAntiTightIsoVBFInvertedcounts, TauAntiTightIsocounts , 0., 0.); 
	double TauAntiMediumIso_eff = vbfefficiency(TauAntiMediumIsoVBFInvertedcounts, TauAntiMediumIsocounts , 0., 0.);
	double Tau2LooseIsoExclusive_eff = vbfefficiency(Tau2LooseIsoExclusiveVBFInvertedcounts, Tau2LooseIsoExclusivecounts , 0., 0.);
	double Tau2LooseIsoInclusive_eff = vbfefficiency(Tau2LooseIsoInclusiveVBFInvertedcounts, Tau2LooseIsoInclusivecounts , 0., 0.);
	double Tau1TightIso_eff = vbfefficiency(Tau1TightIsoVBFInvertedcounts, Tau1TightIsocounts , 0., 0.);
	double TauAnyIso_eff = vbfefficiency(TauAnyIsoVBFInvertedcounts, TauAnyIsocounts , 0., 0.);
	double TauAnyIsoPlusNones_eff = vbfefficiency(TauAnyIsoPlusNonesVBFInvertedcounts, TauAnyIsoPlusNonescounts , 0., 0.);

	double Tau2TightIso_eff_err = vbfefficiency_staterr(Tau2TightIsoVBFInvertedcounts, Tau2TightIsoVBFInvertedcounts_err,Tau2TightIsocounts, Tau2TightIsocounts_err, 0., 0., 0., 0.);
	double Tau2MediumIsoExclusive_eff_err = vbfefficiency_staterr( Tau2MediumIsoExclusiveVBFInvertedcounts, Tau2MediumIsoExclusiveVBFInvertedcounts_err, Tau2MediumIsoExclusivecounts, Tau2MediumIsoExclusivecounts_err, 0., 0., 0., 0.);
	double Tau2MediumIsoInclusive_eff_err = vbfefficiency_staterr( Tau2MediumIsoInclusiveVBFInvertedcounts, Tau2MediumIsoInclusiveVBFInvertedcounts_err, Tau2MediumIsoInclusivecounts, Tau2MediumIsoInclusivecounts_err, 0., 0., 0., 0.);
	double TauAntiTightIso_eff_err = vbfefficiency_staterr( TauAntiTightIsoVBFInvertedcounts, TauAntiTightIsoVBFInvertedcounts_err, TauAntiTightIsocounts, TauAntiTightIsocounts_err, 0., 0., 0., 0.);
	double Tau2LooseIsoExclusive_eff_err = vbfefficiency_staterr( Tau2LooseIsoExclusiveVBFInvertedcounts, Tau2LooseIsoExclusiveVBFInvertedcounts_err, Tau2LooseIsoExclusivecounts, Tau2LooseIsoExclusivecounts_err, 0., 0., 0., 0.);
	double Tau2LooseIsoInclusive_eff_err = vbfefficiency_staterr( Tau2LooseIsoInclusiveVBFInvertedcounts, Tau2LooseIsoInclusiveVBFInvertedcounts_err, Tau2LooseIsoInclusivecounts, Tau2LooseIsoInclusivecounts_err, 0., 0., 0., 0.);
	double TauAntiMediumIso_eff_err = vbfefficiency_staterr( TauAntiMediumIsoVBFInvertedcounts, TauAntiMediumIsoVBFInvertedcounts_err, TauAntiMediumIsocounts, TauAntiMediumIsocounts_err, 0., 0., 0., 0.);
	double Tau1TightIso_eff_err = vbfefficiency_staterr( Tau1TightIsoVBFInvertedcounts, Tau1TightIsoVBFInvertedcounts_err, Tau1TightIsocounts, Tau1TightIsocounts_err, 0., 0., 0., 0.);
	double TauAnyIso_eff_err = vbfefficiency_staterr( TauAnyIsoVBFInvertedcounts, TauAnyIsoVBFInvertedcounts_err, TauAnyIsocounts, TauAnyIsocounts_err, 0., 0., 0., 0.);
	double TauAnyIsoPlusNones_eff_err = vbfefficiency_staterr( TauAnyIsoPlusNonesVBFInvertedcounts, TauAnyIsoPlusNonesVBFInvertedcounts_err, TauAnyIsoPlusNonescounts, TauAnyIsoPlusNonescounts_err, 0., 0., 0., 0.);


	cout << "-------VBF EFFICIENCY--------" << endl;
	cout << "Tau2TightIso_eff: " << Tau2TightIso_eff << " +/- " << Tau2TightIso_eff_err << " ( % " << (Tau2TightIso_eff_err/Tau2TightIso_eff)*100. << " )"<< endl;
	cout << "Tau2MediumIsoExclusive_eff: " << Tau2MediumIsoExclusive_eff << " +/- " << Tau2MediumIsoExclusive_eff_err << " ( % " << (Tau2MediumIsoExclusive_eff_err/Tau2MediumIsoExclusive_eff)*100. << " )"<< endl;
	cout << "Tau2MediumIsoInclusive_eff: " << Tau2MediumIsoInclusive_eff << " +/- " << Tau2MediumIsoInclusive_eff_err << " ( % " << (Tau2MediumIsoInclusive_eff_err/Tau2MediumIsoInclusive_eff)*100. << " )"<< endl;
	cout << "1M (+L+N)_eff: " << TauAntiTightIso_eff << " +/- " << TauAntiTightIso_eff_err << " ( % " << (TauAntiTightIso_eff_err/TauAntiTightIso_eff)*100. << " )"<< endl;
	cout << "Tau2LooseIsoExclusive_eff: " << Tau2LooseIsoExclusive_eff << " +/- " << Tau2LooseIsoExclusive_eff_err << " ( % " << (Tau2LooseIsoExclusive_eff_err/Tau2LooseIsoExclusive_eff)*100. << " )"<< endl;
	cout << "Tau2LooseIsoInclusive_eff: " << Tau2LooseIsoInclusive_eff << " +/- " << Tau2LooseIsoInclusive_eff_err << " ( % " << (Tau2LooseIsoInclusive_eff_err/Tau2LooseIsoInclusive_eff)*100. << " )"<< endl;
	cout << "1L (+N)_eff: " << TauAntiMediumIso_eff << " +/- " << TauAntiMediumIso_eff_err << " ( % " << (TauAntiMediumIso_eff_err/TauAntiMediumIso_eff)*100. << " )"<< endl;
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
	cout << "$ 2 Medium exclusive $    &$ "<< Tau2MediumIsoExclusivecounts << " \\pm  "<< " \\% " << (Tau2MediumIsoExclusivecounts_err/Tau2MediumIsoExclusivecounts)*100. << " $  &$  "
		<< Tau2MediumIsoExclusiveVBFInvertedcounts << "\\pm "<< " \\% " << (Tau2MediumIsoExclusiveVBFInvertedcounts_err/Tau2MediumIsoExclusiveVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ 2 Medium Inclusive $    &$ "<< Tau2MediumIsoInclusivecounts << " \\pm  "<< " \\% " << (Tau2MediumIsoInclusivecounts_err/Tau2MediumIsoInclusivecounts)*100. << " $  &$  "
		<< Tau2MediumIsoInclusiveVBFInvertedcounts << "\\pm "<< " \\% " << (Tau2MediumIsoInclusiveVBFInvertedcounts_err/Tau2MediumIsoInclusiveVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ 1M (+L+N) $    &$ "<< TauAntiTightIsocounts << " \\pm  "<< " \\% " << (TauAntiTightIsocounts_err/TauAntiTightIsocounts)*100. << " $  &$  "
		<< TauAntiTightIsoVBFInvertedcounts << "\\pm "<< " \\% " << (TauAntiTightIsoVBFInvertedcounts_err/TauAntiTightIsoVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ 2 Loose Exclusive$    &$ "<< Tau2LooseIsoExclusivecounts << " \\pm  "<< " \\% " << (Tau2LooseIsoExclusivecounts_err/Tau2LooseIsoExclusivecounts)*100. << " $  &$  "
		<< Tau2LooseIsoExclusiveVBFInvertedcounts << "\\pm "<< " \\% " << (Tau2LooseIsoExclusiveVBFInvertedcounts_err/Tau2LooseIsoExclusiveVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ 2 Loose Inclusive $    &$ "<< Tau2LooseIsoInclusivecounts << " \\pm  "<< " \\% " << (Tau2LooseIsoInclusivecounts_err/Tau2LooseIsoInclusivecounts)*100. << " $  &$  "
		<< Tau2LooseIsoInclusiveVBFInvertedcounts << "\\pm "<< " \\% " << (Tau2LooseIsoInclusiveVBFInvertedcounts_err/Tau2LooseIsoInclusiveVBFInvertedcounts)*100. << " $ \\\\ " << endl;
	cout << "$ 1L (+N) $    &$ "<< TauAntiMediumIsocounts << " \\pm  "<< " \\% " << (TauAntiMediumIsocounts_err/TauAntiMediumIsocounts)*100. << " $  &$  "
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
	cout << " CR     &$ \\epsilon^{QCD}_{VBF}$ \\\\ [0.5ex] \\hline " << endl;
	cout << " Only 2 Tight     &$ "<< Tau2TightIso_eff << " \\pm  "<< " \\% " << (Tau2TightIso_eff_err/Tau2TightIso_eff)*100.  <<  "$  \\\\ " << endl;
	cout << " 1Tight (+M+L+N)     &$ "<< Tau1TightIso_eff << " \\pm  "<< " \\% " << (Tau1TightIso_eff_err/Tau1TightIso_eff)*100.  <<  "$  \\\\ " << endl;
	cout << " 2 Medium Exclusive     &$ "<< Tau2MediumIsoExclusive_eff << " \\pm  "<< " \\% " << (Tau2MediumIsoExclusive_eff_err/Tau2MediumIsoExclusive_eff)*100.  <<  "$  \\\\ " << endl;
	cout << " 2 Medium Inclusive    &$ "<< Tau2MediumIsoInclusive_eff << " \\pm  "<< " \\% " << (Tau2MediumIsoInclusive_eff_err/Tau2MediumIsoInclusive_eff)*100.  <<  "$  \\\\ " << endl;
	cout << " 1M (+L+N)     &$ "<< TauAntiTightIso_eff << " \\pm  "<< " \\% " << (TauAntiTightIso_eff_err/TauAntiTightIso_eff)*100.  <<  "$  \\\\ " << endl;
	cout << " 2 Loose Exclusive     &$ "<< Tau2LooseIsoExclusive_eff << " \\pm  "<< " \\% " << (Tau2LooseIsoExclusive_eff_err/Tau2LooseIsoExclusive_eff)*100.  <<  "$  \\\\ " << endl;
	cout << " 2 Loose Inclusive     &$ "<< Tau2LooseIsoInclusive_eff << " \\pm  "<< " \\% " << (Tau2LooseIsoInclusive_eff_err/Tau2LooseIsoInclusive_eff)*100.  <<  "$  \\\\ " << endl;
	cout << " 1L (+N)     &$ "<< TauAntiMediumIso_eff << " \\pm  "<< " \\% " << (TauAntiMediumIso_eff_err/TauAntiMediumIso_eff)*100.  <<  "$  \\\\ " << endl;
	cout << " Any iso (+N)     &$ "<< TauAnyIsoPlusNones_eff << " \\pm  "<< " \\% " << (TauAnyIsoPlusNones_eff_err)*100.  <<  "$  \\\\ " << endl;
	cout << "\\hline\\hline" << endl;
	cout << "\\end{tabular}" << endl;
	cout << "} " << endl;
	cout << "\\end{table}" << endl;


}
