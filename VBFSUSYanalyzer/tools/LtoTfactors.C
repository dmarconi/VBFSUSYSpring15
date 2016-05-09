#include <TROOT.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double binerrors(TH1F* h1, int minMETbin, int maxMETbin){
	double temperror = 0;
	for (int i = minMETbin; i < maxMETbin + 1; i++){
		temperror = temperror + pow(h1->GetBinError(i),2.);
	}
	return sqrt(temperror);
}

void LtoTfactors() {

	TFile *inputfile = TFile::Open("allQCD.root");

	TH1F* h1_counts;
	h1_counts = ((TH1F*)(inputfile->Get("demo/baselineObjectSelection/counts")));

	double fourJetscounts = h1_counts->GetBinContent(2);
	double twoLooseTaus = h1_counts->GetBinContent(3);
	double oneLcounts = h1_counts->GetBinContent(4);
	double alsoTcounts = h1_counts->GetBinContent(5);

	double fourJetscounts_err = binerrors(h1_counts, 2, 2);
	double twoLooseTaus_err = binerrors(h1_counts, 3, 3);
	double oneLcounts_err = binerrors(h1_counts, 4, 4);
	double alsoTcounts_err = binerrors(h1_counts, 5, 5);
	
	double looseToTightProb = alsoTcounts / oneLcounts;
	double looseToTightProb_err = sqrt( pow( (alsoTcounts_err / oneLcounts )   , 2.) + pow( ( (alsoTcounts * oneLcounts_err) / (oneLcounts * oneLcounts)     )    , 2.) );

	double twoLooseTo2Tight = looseToTightProb * looseToTightProb;
	double twoLooseTo2Tight_err = 2. * looseToTightProb * looseToTightProb_err;

	cout << endl;
	cout << endl;
	cout << "-------EVENT COUNT--------" << endl;
	cout << "4Jetsevents: " << fourJetscounts << " +/- % " << (fourJetscounts_err/fourJetscounts)*100. << " )"<< endl;
	cout << "at least 2 isolates Taus: " << twoLooseTaus << " +/- % " << (twoLooseTaus_err/twoLooseTaus)*100. << " )"<< endl;
	cout << "1 fake Loose Tau found: " << oneLcounts << " +/- % " << (oneLcounts_err/oneLcounts)*100. << " )"<< endl;
	cout << "is also a Tight Tau: " << alsoTcounts << " +/- % " << (alsoTcounts_err/alsoTcounts)*100. << " )"<< endl;
	cout << endl;
	cout << endl;
	cout << "1 Loose to 1 Tight probability: " << looseToTightProb << " +/- % " << (looseToTightProb_err/looseToTightProb)*100. << " )"<< endl;
	cout << endl;
	cout << endl;
	cout << "2 Loose to 2 Tight probability: " << twoLooseTo2Tight << " +/- % " << (twoLooseTo2Tight_err/twoLooseTo2Tight)*100. << " )"<< endl;
	

}
