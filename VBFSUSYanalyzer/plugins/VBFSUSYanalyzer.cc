// -*- C++ -*-
//
// Package:    VBFSUSYSpring15/VBFSUSYanalyzer
// Class:      VBFSUSYanalyzer
// 
/**\class VBFSUSYanalyzer VBFSUSYanalyzer.cc VBFSUSYSpring15/VBFSUSYanalyzer/plugins/VBFSUSYanalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Daniele Marconi
//         Created:  Sun, 17 Jan 2016 14:11:37 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <utility>

// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// -- Root

#include "TROOT.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

// useful structs

struct MassAndIndex {
	std::string label;
	unsigned int first;
	unsigned int second;
	double Mass;
	double dR;
	int signEta;
	double dEta;

	MassAndIndex (const std::string & inputlabel){
		label = inputlabel;
		first = 99999;
		second = 99999;
		Mass = -1.;
		dR = 0.;
		signEta = 1;
		dEta = 0;
	}
};

struct TauProperties {
	std::string label;
	unsigned int first;
	unsigned int second;
	double Mass;
	double dR;
	double dEta;
	double cosDphi;
	int charge;

	TauProperties (const std::string & inputlabel){
		label = inputlabel;
		first = 99999;
		second = 99999;
		Mass = -1.;
		dR = 0.;
		dEta=0.;
		cosDphi = 1;
		charge = -1;
	}
};


struct MyEventCollection {

	std::string label;
	bool goodVertex;
	std::vector <const pat::Tau*> tau;
	std::vector <const pat::Jet*> jet;
	std::vector <const pat::Jet*> bjet;
	std::vector <const pat::MET*> met;
	int NVtx;
	int PUinteractions;

	void init(const std::string & inputlabel) {
		label = inputlabel;
		goodVertex = false;
		tau.clear();
		jet.clear();
		bjet.clear();
		met.clear();
	}
	void clear() {
		goodVertex = false;
		tau.clear();
		jet.clear();
		bjet.clear();
		met.clear();
	}
};

struct MyHistoCollection {

	std::string label;
	edm::Service<TFileService> fs;
	TH1F* h_count;
	TH1F* h_njet;
	TH1F* h_jetpt;
	TH1F* h_jeteta;
	TH1F* h_jet1pt;
	TH1F* h_jet1eta;
	TH1F* h_jet2pt;
	TH1F* h_jet2eta;
	TH1F* h_dijetinvariantmass;
	TH1F* h_dijetdeltaeta;

	TH1F* h_tau1pt;
	TH1F* h_tau1eta;
	TH1F* h_tau2pt;
	TH1F* h_tau2eta;
	TH1F* h_ditauinvariantmass;
	TH1F* h_ditaucharge;
	TH1F* h_ditaucosdeltaphi;
	TH1F* h_ditaudeltaeta;

	TH1F* h_met;

	TH1F* h_ht;
	TH1F* h_ht_withtau;

	TH1F* h_jetTauDistanceFirst;
	TH1F* h_jetTauDistanceSecond;

	TH2F* h2_DiJetInvMass_vs_DiJetDEta;
	TH2F* h2_tau1pt_vs_tau2pt;

	void init(const std::string & inputlabel) {

		label = inputlabel;
		TFileDirectory subDir = fs->mkdir( inputlabel );
		//f->mkdir(inputlabel.c_str());
		//f->cd(inputlabel.c_str());

		h_count = subDir.make<TH1F>("counts", "", 1,0,1);
		h_count->SetBit(TH1::kCanRebin);
		h_count->SetStats(0);
		h_count->Fill("NoCuts",0);

		h_njet = subDir.make<TH1F>("h_njet", "h_njet", 21, -0.5, 20.5);
		h_njet->GetXaxis()->SetTitle("number of jets not matched to #tau");
		h_jetpt = subDir.make<TH1F>("h_jetpt", "h_jetpt", 50, 0., 500.);
		h_jetpt->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
		h_jeteta = subDir.make<TH1F>("h_jeteta", "h_jeteta", 30 , -5., 5.);
		h_jeteta->GetXaxis()->SetTitle("#eta^{jet}");
		h_jet1pt = subDir.make<TH1F>("h_jet1pt", "h_jet1pt", 50, 0., 500.);
		h_jet1pt->GetXaxis()->SetTitle("p_{T}^{jet 1} [GeV]");
		h_jet1eta = subDir.make<TH1F>("h_jet1eta", "h_jet1eta", 50 , -5., 5.);
		h_jet1eta->GetXaxis()->SetTitle("#eta^{jet 1}");
		h_jet2pt = subDir.make<TH1F>("h_jet2pt", "h_jet2pt", 50, 0., 500.);
		h_jet2pt->GetXaxis()->SetTitle("p_{T}^{jet 2} [GeV]");
		h_jet2eta = subDir.make<TH1F>("h_jet2eta", "h_jet2eta", 50 , -5., 5.);
		h_jet2eta->GetXaxis()->SetTitle("#eta^{jet 2}");
		h_dijetinvariantmass = subDir.make<TH1F>("h_dijetinvariantmass","h_dijetinvariantmass", 10, 0., 2500.);
		h_dijetinvariantmass->GetXaxis()->SetTitle("M^{(jet,jet)} [GeV]");
		h_dijetdeltaeta = subDir.make<TH1F> ("h_dijetdeltaeta", "h_dijetdeltaeta", 20, 0., 10.);
		h_dijetdeltaeta->GetXaxis()->SetTitle("#Delta#eta^{jj}");

		h_tau1pt = subDir.make<TH1F>("h_tau1pt", "h_tau1pt", 50, 0., 500.);
		h_tau1pt->GetXaxis()->SetTitle("p_{T}^{#tau 1} [GeV]");
		h_tau1eta = subDir.make<TH1F>("h_tau1eta", "h_tau1eta", 30 , -3., 3.);
		h_tau1eta->GetXaxis()->SetTitle("#eta^{#tau 1}");
		h_tau2pt = subDir.make<TH1F>("h_tau2pt", "h_tau2pt", 50, 0., 500.);
		h_tau2pt->GetXaxis()->SetTitle("p_{T}^{#tau 2} [GeV]");
		h_tau2eta = subDir.make<TH1F>("h_tau2eta", "h_tau2eta", 30 , -3., 3.);
		h_tau2eta->GetXaxis()->SetTitle("#eta^{#tau 2}");
		h_ditauinvariantmass = subDir.make<TH1F>("h_ditauinvariantmass", "h_ditauinvariantmass", 50, 0., 500.);
		h_ditauinvariantmass->GetXaxis()->SetTitle("M^{(#tau,#tau)} [GeV]");
		h_ditaucharge = subDir.make<TH1F>("h_ditaucharge", "h_ditaucharge", 5, -4., 6.);
		h_ditaucharge->GetXaxis()->SetTitle("sign(#tau^{1}) #upoint sign(#tau^{2})");
		h_ditaucosdeltaphi = subDir.make<TH1F>("h_ditaucosdeltaphi", "h_ditaucosdeltaphi", 50, -1.1, 1.1);
		h_ditaucosdeltaphi->GetXaxis()->SetTitle("cos(#Delta#phi(#tau,#tau))");
		h_ditaudeltaeta = subDir.make<TH1F> ("h_ditaudeltaeta", "h_ditaudeltaeta", 20, 0., 10.);
		h_ditaudeltaeta->GetXaxis()->SetTitle("#Delta#eta(#tau,#tau)");

		h_met = subDir.make<TH1F>("h_met", "h_met", 24, 0., 240.);
		h_met->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");

		h_ht = subDir.make<TH1F>("h_ht", "h_ht", 50, 0., 1300.);
		h_ht->GetXaxis()->SetTitle("H_{T} [GeV]");
		h_ht_withtau = subDir.make<TH1F>("h_ht_withtau", "h_ht_withtau", 50, 0., 1300.);
		h_ht_withtau->GetXaxis()->SetTitle("H_{T}+#Sigma p_{T}^{#tau} [GeV]");

		h_jetTauDistanceFirst = subDir.make<TH1F>("h_jetTauDistanceFirst", "h_jetTauDistanceFirst", 25, 0., 0.5);
		h_jetTauDistanceFirst->GetXaxis()->SetTitle("#DeltaR(jet,#tau^{1})");
		h_jetTauDistanceSecond = subDir.make<TH1F>("h_jetTauDistanceSecond", "h_jetTauDistanceSecond", 25, 0., 0.5);
		h_jetTauDistanceSecond->GetXaxis()->SetTitle("#DeltaR(jet,#tau^{2})");

		h2_DiJetInvMass_vs_DiJetDEta = new TH2F("h2_DiJetInvMass_vs_DiJetDEta","h2_DiJetInvMass_vs_DiJetDEta", 20, 0., 10., 10, 0., 2500.);
		h2_DiJetInvMass_vs_DiJetDEta->GetXaxis()->SetTitle("#Delta#eta^{jj}");
		h2_DiJetInvMass_vs_DiJetDEta->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
		h2_tau1pt_vs_tau2pt = new TH2F("h2_tau1pt_vs_tau2pt","correlation of first and second p_{T}^{#tau}", 50, 0., 500., 50, 0., 500.);
		h2_tau1pt_vs_tau2pt->GetXaxis()->SetTitle("p_{T}^{#tau 1}");
		h2_tau1pt_vs_tau2pt->GetYaxis()->SetTitle("p_{T}^{#tau 2}");

	}
};


//useful functions
//-----------------
//leading jet finder
//_________________

std::pair<unsigned int,unsigned int> LeadingJets(MyEventCollection collection)
{
	unsigned int     	j1=99999;
	unsigned int     	j2=99999;
	double 		pt1=-1.;
	double 		pt2=-1.;

	for(unsigned int j=0; j<collection.jet.size();j++)
	{
		if(pt1 < collection.jet[j]->pt()){

			if(pt2 < pt1){ pt2=pt1; j2=j1; }

			j1=j;
			pt1=collection.jet[j]->pt();

		}

		if((pt2 < collection.jet[j]->pt()) && (pt1 > collection.jet[j]->pt())){ j2=j; pt2=collection.jet[j]->pt(); }
	}

	return std::make_pair(j1,j2);
}

//-----------------
//leading tau finder
//_________________

std::pair<unsigned int,unsigned int> LeadingTaus(MyEventCollection collection)
{
	unsigned int     	t1=99999;
	unsigned int     	t2=99999;
	double 		pt1=-1.;
	double 		pt2=-1.;

	for(unsigned int t=0; t<collection.tau.size();t++)
	{
		if(pt1 < collection.tau[t]->pt()){

			if(pt2 < pt1){ pt2=pt1; t2=t1; }

			t1=t;
			pt1=collection.tau[t]->pt();

		}

		if((pt2 < collection.tau[t]->pt()) && (pt1 > collection.tau[t]->pt())){ t2=t; pt2=collection.tau[t]->pt(); }
	}

	return std::make_pair(t1,t2);
}

double TauJetMinDistance(MyEventCollection collection, const pat::Jet &jet)
{
	double minDeltaRtauJet = 99999.;
	for(unsigned int t =0;t<collection.tau.size();++t){
		TLorentzVector v_tau;
		TLorentzVector v_jet;
		v_tau.SetPtEtaPhiM(collection.tau[t]->pt(),collection.tau[t]->eta(),collection.tau[t]->phi(),collection.tau[t]->mass());
		v_jet.SetPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.mass());
		double temp_mindeltaRtaujet = v_tau.DeltaR(v_jet); 
		if (temp_mindeltaRtaujet < minDeltaRtauJet) minDeltaRtauJet = temp_mindeltaRtaujet;
	} 
	return minDeltaRtauJet;
}

//--------------------------
//2-jet max inv. mass finder and properties of dijet system
//__________________________

MassAndIndex Inv2jMassIndex(MyEventCollection collection)
{	
	struct MassAndIndex Inv2jMass("Inv2jMass");
	TLorentzVector jet1_4v;
	TLorentzVector jet2_4v; 

	double Mass=0.;
	double dEtaCheck=0.;
	unsigned int first=99999;
	unsigned int second=99999;

	bool pass=false; 	//passing all selection cuts for the dijet system
	bool failOne=false;	//failing dijet eta cut
	bool failTwo=false;	//failing dijet eta sign cut

	for(unsigned int j1=1; j1<collection.jet.size(); j1++)
	{
		for(unsigned int j2=0; j2<j1; j2++)
		{
			jet1_4v.SetPtEtaPhiE(collection.jet[j1]->pt(), collection.jet[j1]->eta(), collection.jet[j1]->phi(), collection.jet[j1]->energy());
			jet2_4v.SetPtEtaPhiE(collection.jet[j2]->pt(), collection.jet[j2]->eta(), collection.jet[j2]->phi(), collection.jet[j2]->energy());

			int sign=collection.jet[j1]->eta() * collection.jet[j2]->eta();
			double dEta=fabs(collection.jet[j1]->eta() - collection.jet[j2]->eta());

			TLorentzVector dijet_4v = jet1_4v + jet2_4v;

			//select the pair of dijets with highest invariant mass passing all cuts
			if(dEta>=4.2 && sign<0 && dijet_4v.M()>=250+Mass) {Mass = dijet_4v.M(); first = j1; second = j2; pass=true;}
			else if(!pass){
				if(dEta>dEtaCheck && sign<0 && dijet_4v.M()>250) {Mass = dijet_4v.M(); first = j1; second = j2; dEtaCheck=dEta; failOne=true;} //search for pair failing only the dijet delta eta cut
				else if(!failOne)
				{
					if(dEta>dEtaCheck && dijet_4v.M()>250) {Mass = dijet_4v.M(); first = j1; second = j2; dEtaCheck=dEta; failTwo=true;} //search for a pair failing dijet delta eta and eta sign cuts
					else if(!failTwo && dijet_4v.M()>Mass && dEta>dEtaCheck) {Mass = dijet_4v.M(); first = j1; second = j2; dEtaCheck=dEta;}  //failing all cuts
				}
			}
		}
	}
	if(first < 99999 && second < 99999)
	{
		double dR=jet1_4v.DeltaR(jet2_4v);
		int sign=(collection.jet[first]->eta() * collection.jet[second]->eta())/fabs(collection.jet[first]->eta() * collection.jet[second]->eta());
		double dEta=fabs(collection.jet[first]->eta() - collection.jet[second]->eta());
		Inv2jMass.Mass=Mass;
		Inv2jMass.first=first;
		Inv2jMass.second=second;
		Inv2jMass.dR=dR;
		Inv2jMass.signEta=sign;
		Inv2jMass.dEta=dEta;
		//std::cout<<"chose jets "<<first<<" and "<<second<<std::endl;
	}    
	return Inv2jMass;
}

//--------------------------
//2-tau system properties
//__________________________  

TauProperties Inv2tMassIndex(MyEventCollection collection)
{
	struct TauProperties Inv2tMass("Inv2tMass");
	TLorentzVector tau1_4v;
	TLorentzVector tau2_4v; 

	std::pair<unsigned int,unsigned int> tauIndex=LeadingTaus(collection);
	Inv2tMass.first = tauIndex.first;
	Inv2tMass.second = tauIndex.second;

	if(tauIndex.first < 99999 && tauIndex.second < 99999)
	{
		tau1_4v.SetPtEtaPhiM(collection.tau[tauIndex.first]->pt(), collection.tau[tauIndex.first]->eta(), collection.tau[tauIndex.first]->phi(), collection.tau[tauIndex.first]->mass());
		tau2_4v.SetPtEtaPhiM(collection.tau[tauIndex.second]->pt(), collection.tau[tauIndex.second]->eta(), collection.tau[tauIndex.second]->phi(), collection.tau[tauIndex.second]->mass());

		TLorentzVector ditau_4v = tau1_4v + tau2_4v;

		double dR=tau1_4v.DeltaR(tau2_4v);
		int charge = collection.tau[tauIndex.first]->charge() * collection.tau[tauIndex.second]->charge();
		double cosdeltaphiDiTau = cos(tau1_4v.DeltaPhi(tau2_4v));
		double dEta=fabs(collection.tau[tauIndex.first]->eta() - collection.tau[tauIndex.second]->eta());

		Inv2tMass.Mass = ditau_4v.M();
		Inv2tMass.dR = dR;
		Inv2tMass.charge = charge;
		Inv2tMass.cosDphi = cosdeltaphiDiTau;
		Inv2tMass.dEta=dEta;
	}
	return Inv2tMass;
}

// Fill histograms


void fillHistoCollection (MyHistoCollection &inputHistoCollection, MyEventCollection inputEventCollection, double weight) {

	// ---------------------
	// -- fill histograms --
	// ---------------------	  
	
	using namespace std;


	bool verbose=false;  
	//JETS	  

	//find index of leading jets
	std::pair<unsigned int,unsigned int> jetIndex=LeadingJets(inputEventCollection);  

	//find max value for 2-jet-mass
	MassAndIndex Inv2j = Inv2jMassIndex(inputEventCollection);

	//define ht
	double ht_jets=0.;

	//JET SEL
	for (unsigned int j = 0;j<inputEventCollection.jet.size();++j){
		inputHistoCollection.h_jetpt->Fill(inputEventCollection.jet[j]->pt(),weight); //fill jet-pt-histogram
		inputHistoCollection.h_jeteta->Fill(inputEventCollection.jet[j]->eta(),weight); //fill jet-eta-histogram
		ht_jets+=inputEventCollection.jet[j]->pt(); //add up scalar pt to ht
	}

	//fill jet count
	inputHistoCollection.h_njet->Fill( (int)inputEventCollection.jet.size(),weight );
	if(verbose)std::cout<<"Pass selection -> Fill njet="<<inputEventCollection.jet.size()<<", weight="<<weight<<std::endl;

	//fill jet pt indizes
	if (jetIndex.first < 99999)  {
		inputHistoCollection.h_jet1pt->Fill(inputEventCollection.jet[jetIndex.first]->pt(),weight);
		inputHistoCollection.h_jet1eta->Fill(inputEventCollection.jet[jetIndex.first]->eta(),weight);
	}
	if (jetIndex.second < 99999) {
		inputHistoCollection.h_jet2pt->Fill(inputEventCollection.jet[jetIndex.second]->pt(),weight);
		inputHistoCollection.h_jet2eta->Fill(inputEventCollection.jet[jetIndex.second]->eta(),weight);
	}

	//fill 2-jet-event inv. mass and eta-difference
	if ( Inv2j.Mass >= 0 ) {
		inputHistoCollection.h_dijetinvariantmass ->Fill(Inv2j.Mass,weight);
		inputHistoCollection.h_dijetdeltaeta ->Fill(Inv2j.dEta,weight);
	}

	//fill ht distribution
	inputHistoCollection.h_ht -> Fill(ht_jets,weight);

	//fill jet tau distance distribution
	//if(jetIndex.first<99999 && jetIndex.second<99999){
	//	inputHistoCollection.h_jetTauDistanceFirst->Fill(TauJetMinDistance(inputEventCollection,inputEventCollection.jet[jetIndex.first]->eta, inputEventCollection.jet[jetIndex.first]->phi));
	//	inputHistoCollection.h_jetTauDistanceSecond->Fill(TauJetMinDistance(inputEventCollection,inputEventCollection.jet[jetIndex.second]->eta, inputEventCollection.jet[jetIndex.second]->phi));	 
	//}
	//____________________________________________________________________________________________

	//TAUS

	//find tau properties
	TauProperties Inv2t = Inv2tMassIndex(inputEventCollection);
	TLorentzVector tau1_4v;
	TLorentzVector tau2_4v;

	//set ht of taus to default
	double ht_jetsPtau=ht_jets;

	for(unsigned int t =0;t<inputEventCollection.tau.size();++t){
		//add up scalar sum of tau pt to ht
		ht_jetsPtau+=inputEventCollection.tau[t]->pt();
	}

	if ( (Inv2t.first < 99999) && (Inv2t.second < 99999) ) {
		//determine leading two tau invariant mass
		inputHistoCollection.h_ditauinvariantmass ->Fill(Inv2t.Mass,weight);

		//fill tau charge and  cosdeltaphi and deltaeta and 2Dpt-plot
		inputHistoCollection.h_ditaucharge ->Fill(Inv2t.charge,weight);
		inputHistoCollection.h_ditaucosdeltaphi ->Fill(Inv2t.cosDphi,weight);	
		inputHistoCollection.h_ditaudeltaeta->Fill(Inv2t.dEta,weight);
		inputHistoCollection.h2_tau1pt_vs_tau2pt->Fill(inputEventCollection.tau[Inv2t.first]->pt(),inputEventCollection.tau[Inv2t.second]->pt(),weight);   
	}

	//fill tau pt and eta
	if (Inv2t.first < 99999) {
		inputHistoCollection.h_tau1pt->Fill(inputEventCollection.tau[Inv2t.first]->pt(),weight);
		inputHistoCollection.h_tau1eta->Fill(inputEventCollection.tau[Inv2t.first]->eta(),weight);
	}
	if (Inv2t.second < 99999) {
		inputHistoCollection.h_tau2pt->Fill(inputEventCollection.tau[Inv2t.second]->pt(),weight);
		inputHistoCollection.h_tau2eta->Fill(inputEventCollection.tau[Inv2t.second]->eta(),weight);
	}

	//fill ht with taus included
	inputHistoCollection.h_ht_withtau -> Fill(ht_jetsPtau,weight);

	// MET
	inputHistoCollection.h_met -> Fill(inputEventCollection.met[0]->pt(),weight);

	//fill DiJetInvMass_vs_DiJetDEta
	inputHistoCollection.h2_DiJetInvMass_vs_DiJetDEta -> Fill(Inv2j.dEta, Inv2j.Mass,weight);

	//________________________________________

} 

// class declaration
//

class VBFSUSYanalyzer : public edm::EDAnalyzer {
	public:
		explicit VBFSUSYanalyzer(const edm::ParameterSet&);
		~VBFSUSYanalyzer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
		//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
		//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;


		// ---------event collections-----------------------------

		MyEventCollection baselineObjectSelectionCollection;
		MyEventCollection TauLooseIsoObjectSelectionCollection;
		MyEventCollection TauMediumIsoObjectSelectionCollection;
		MyEventCollection Tau1TightIsoObjectSelectionCollection;
		MyEventCollection TauTightIsoObjectSelectionCollection;

		// ---------histograms-----------------------------
		edm::Service<TFileService> fs;	
		TH1F* count;
		MyHistoCollection myHistoColl_baselineSelection;
		MyHistoCollection myHistoColl_TauLooseIsoObjectSelection;
		MyHistoCollection myHistoColl_TauMediumIsoObjectSelection;
		MyHistoCollection myHistoColl_Tau1TightIsoObjectSelection;
		MyHistoCollection myHistoColl_TauTightIsoObjectSelection;

		// ----------member data ---------------------------

		edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
		edm::EDGetTokenT<pat::MuonCollection> muonToken_;
		edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
		edm::EDGetTokenT<pat::TauCollection> tauToken_;
		edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
		edm::EDGetTokenT<pat::JetCollection> jetToken_;
		edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
		edm::EDGetTokenT<pat::METCollection> metToken_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
VBFSUSYanalyzer::VBFSUSYanalyzer(const edm::ParameterSet& iConfig):

	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
	electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
	tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
	fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
	metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))

{

	//Event collection init

	baselineObjectSelectionCollection.init("baselineObjectSelection");
	TauLooseIsoObjectSelectionCollection.init("TauLooseIsoObjectSelection");
	TauMediumIsoObjectSelectionCollection.init("TauMediumIsoObjectSelection");
	TauTightIsoObjectSelectionCollection.init("Tau1TightIsoObjectSelection");
	TauTightIsoObjectSelectionCollection.init("TauTightIsoObjectSelection");


	//histogram initialization	
	count = fs->make<TH1F>("counts", "", 1,0,1);
	count->SetBit(TH1::kCanRebin);
	count->SetStats(0);
	count->Fill("NoCuts",0);
	myHistoColl_baselineSelection.init("baselineObjectSelection");
	myHistoColl_TauLooseIsoObjectSelection.init("TauLooseIsoObjectSelection");
	myHistoColl_TauLooseIsoObjectSelection.h_count->Fill("AtLeast2taus",0.);	
	myHistoColl_TauMediumIsoObjectSelection.init("TauMediumIsoObjectSelection");
	myHistoColl_TauMediumIsoObjectSelection.h_count->Fill("AtLeast2taus",0.);	
	myHistoColl_Tau1TightIsoObjectSelection.init("Tau1TightIsoObjectSelection");
	myHistoColl_Tau1TightIsoObjectSelection.h_count->Fill("AtLeast2taus",0.);	
	myHistoColl_TauTightIsoObjectSelection.init("TauTightIsoObjectSelection");
	myHistoColl_TauTightIsoObjectSelection.h_count->Fill("AtLeast2taus",0.);	


}


VBFSUSYanalyzer::~VBFSUSYanalyzer()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
VBFSUSYanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace std;


	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vtxToken_, vertices);
	if (vertices->empty()) return; // skip the event if no PV found
	const reco::Vertex &PV = vertices->front();

	edm::Handle<pat::MuonCollection> muons;
	iEvent.getByToken(muonToken_, muons);
	//for (const pat::Muon &mu : *muons) {
	//	if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
	//	printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
	//			mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
	//}

	edm::Handle<pat::ElectronCollection> electrons;
	iEvent.getByToken(electronToken_, electrons);
	//for (const pat::Electron &el : *electrons) {
	//	if (el.pt() < 5) continue;
	//printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d\n",
	//		el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(), el.passConversionVeto());
	//}

	edm::Handle<pat::PhotonCollection> photons;
	iEvent.getByToken(photonToken_, photons);
	//for (const pat::Photon &pho : *photons) {
	//	if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
	//	printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
	//			pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
	//}


	edm::Handle<pat::TauCollection> taus;
	iEvent.getByToken(tauToken_, taus);
	//for (const pat::Tau &tau : *taus) {
	//	if (tau.pt() < 20) continue;
	//	printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
	//			tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
	//}


	edm::Handle<pat::JetCollection> jets;
	iEvent.getByToken(jetToken_, jets);
	//int ijet = 0;
	//for (const pat::Jet &j : *jets) {
	//	if (j.pt() < 20) continue;
	//	printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n",
	//			j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")), j.userFloat("pileupJetId:fullDiscriminant"));
	//	if ((++ijet) == 1) { // for the first jet, let's print the leading constituents
	//		std::vector<reco::CandidatePtr> daus(j.daughterPtrVector());
	//		std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
	//		for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
	//			const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
	//			printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i2,cand.pt(),cand.dz(PV.position()),cand.pdgId());
	//		}
	//	}
	//}


	edm::Handle<pat::METCollection> mets;
	iEvent.getByToken(metToken_, mets);
	const pat::MET &met = mets->front();
	//printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
	//		met.pt(), met.phi(), met.sumEt(),
	//		met.genMET()->pt(),
	//		met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));

	//printf("\n");

	//std::vector<int> tights;
	//std::vector<int> mediums;
	//std::vector<int> looses;
	//std::vector<int> nones;

	std::vector<const pat::Tau*> tights;
	std::vector<const pat::Tau*> mediums;
	std::vector<const pat::Tau*> looses;
	std::vector<const pat::Tau*> nones;

	//smart tau selection
	//cout << "DEBUG: Tau vector size: " << tau.size() << endl;
	for (const pat::Tau &tau : *taus) {
		//cout << "DEBUG: Tau counter: " << t << endl;
		if(!(	fabs(tau.eta()) <= 2.1                              					)) continue;
		//cout << "DEBUG: Eta cut passed for tau : " << t << endl;
		//cout << "DEBUG: Tau pt : " << tau[t].pt << endl;
		//cout << "DEBUG: Tau eta : " << tau[t].eta << endl;
		//cout << "DEBUG: Tau phi : " << tau[t].phi << endl;
		//OLDID    //if(!(       tau[t].pt >= 45.                                            				)) continue;
		if(!(       tau.pt() >= 20.                                            				)) continue;
		//cout << "DEBUG: Pt cut passed for tau : " << t << endl;
		//OLDID  //if(!(       tau[t].tauID_againstElectronMediumMVA5 > 0.5                				)) continue;
		if(!(       tau.tauID("againstElectronVLooseMVA5") > 0.5                				)) continue;
		//cout << "DEBUG: Electron veto cut passed for tau : " << t << endl;
		if(!(       tau.tauID("againstMuonLoose3") > 0.5                        				)) continue;
		//cout << "DEBUG: Muon veto cut passed for tau : " << t << endl;
		if(!(       tau.leadChargedHadrCand()->pt() >= 5.0                      				)) continue;
		//cout << "DEBUG: leadChargedHadrCand_pt cut passed for tau : " << t << endl;
		if(!(       (tau.tauID("decayModeFindingNewDMs") > 0.5) && (tau.signalChargedHadrCands().size() == 1)	)) continue;
		//cout << "DEBUG: decayModeFindingNewDMs cut passed for tau : " << t << endl;

		baselineObjectSelectionCollection.tau.push_back(&tau);

		if(!(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")  <= 0.5)) tights.push_back(&tau);
		else if(!(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits")  <= 0.5)) mediums.push_back(&tau);
		else if(!(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")  <= 0.5)) looses.push_back(&tau);
		else nones.push_back(&tau);

	}

	if(tights.size()==2) for(unsigned int t =0;t<tights.size();++t) {TauTightIsoObjectSelectionCollection.tau.push_back(tights[t]);}
	else if(tights.size()==1 && (mediums.size()+looses.size()+nones.size())==1) {tights.insert(tights.end(),mediums.begin(), mediums.end()); tights.insert(tights.end(),looses.begin(), looses.end()); tights.insert(tights.end(),nones.begin(), nones.end()); for(unsigned int t =0;t<tights.size();++t) {Tau1TightIsoObjectSelectionCollection.tau.push_back(tights[t]);}}
	else if(mediums.size()>=1 && (mediums.size()+looses.size()+nones.size())==2) {mediums.insert(mediums.end(), looses.begin(), looses.end()); mediums.insert(mediums.end(), nones.begin(), nones.end()); for(unsigned int t =0;t<mediums.size();++t) {TauMediumIsoObjectSelectionCollection.tau.push_back(mediums[t]);}}
	else if(looses.size()>=1 && (looses.size()+nones.size())==2) {looses.insert(looses.end(), nones.begin(), nones.end()); for(unsigned int t =0;t<looses.size();++t) {TauLooseIsoObjectSelectionCollection.tau.push_back(looses[t]);}}
	//else if(nones.size()==2) for(unsigned int t =0;t<nones.size();++t) {int i=nones[t]; TauNoIsoObjectSelectionCollection.tau.push_back(&tau[i]);}


	// jet && bjet selection
	// ? id ?
	//cout << "jet.size(): " << jet.size() << endl;
	for(const pat::Jet &jet : *jets){
		if(!(      jet.pt() >= 30.                                                	)) continue;  // Original value 20
		if(!(      fabs(jet.eta()) <= 5.0                                          )) continue;

		double baseDistance = TauJetMinDistance(baselineObjectSelectionCollection, jet);
		double mainDistance = TauJetMinDistance(TauTightIsoObjectSelectionCollection, jet);
		double T1Distance = TauJetMinDistance(Tau1TightIsoObjectSelectionCollection, jet);
		double mediumDistance = TauJetMinDistance(TauMediumIsoObjectSelectionCollection, jet);
		double looseDistance = TauJetMinDistance(TauLooseIsoObjectSelectionCollection, jet);
		//double NoDistance = TauJetMinDistance(TauNoIsoObjectSelectionCollection, jet[j].eta, jet[j].phi);

		bool jetid=true;
		if(!(      jet.neutralHadronEnergyFraction() < 0.99                                        )) jetid=false;
		if(!(      jet.neutralEmEnergyFraction() < 0.99                                            )) jetid=false;
		if(!(      jet.numberOfDaughters() > 1                                                     )) jetid=false;
		if(fabs(jet.eta()) < 2.4) {
			if(!(      jet.chargedHadronEnergyFraction() > 0                        			)) jetid=false;
			if(!(      jet.chargedEmEnergyFraction() < 0.99                            		)) jetid=false;
			if(!(      jet.chargedHadronMultiplicity() > 0					)) jetid=false;
		}
		
		//Filling Jet collection
		if(      /*jet[j].pt >= 50.  &&*/ jetid		){
			if(	baseDistance >= 0.3	) baselineObjectSelectionCollection.jet.push_back(&jet);	
				  if(	mainDistance >= 0.3	) TauTightIsoObjectSelectionCollection.jet.push_back(&jet);
				  if(	T1Distance >= 0.3	) Tau1TightIsoObjectSelectionCollection.jet.push_back(&jet);
				  if(	mediumDistance >= 0.3	) TauMediumIsoObjectSelectionCollection.jet.push_back(&jet);
				  if(	looseDistance >= 0.3	) TauLooseIsoObjectSelectionCollection.jet.push_back(&jet);
			//	  if(	NoDistance  >= 0.3	) TauNoIsoObjectSelectionCollection.jet.push_back(&jet[j]);
		}
		
		//Filling bJet collection
		if(fabs(jet.eta()) <= 2.4 && jet.bDiscriminator("combinedSecondaryVertexBJetTags") /*jet[j].bDiscriminator_combinedSecondaryVertexBJetTags*/ > 0.244    ){
			if(	baseDistance >= 0.3	) baselineObjectSelectionCollection.bjet.push_back(&jet);	
						  if(	mainDistance >= 0.3	) TauTightIsoObjectSelectionCollection.bjet.push_back(&jet);
						  if(	T1Distance >= 0.3	) Tau1TightIsoObjectSelectionCollection.bjet.push_back(&jet);
						  if(	mediumDistance >= 0.3	) TauMediumIsoObjectSelectionCollection.bjet.push_back(&jet);
						  if(	looseDistance >= 0.3	) TauLooseIsoObjectSelectionCollection.bjet.push_back(&jet);
			//			  if(	NoDistance  >= 0.3	) TauNoIsoObjectSelectionCollection.bjet.push_back(&jet[j]);
		}
	}	

	//MET selection
	baselineObjectSelectionCollection.met.push_back(&met);
	TauLooseIsoObjectSelectionCollection.met.push_back(&met);
	TauMediumIsoObjectSelectionCollection.met.push_back(&met);
	Tau1TightIsoObjectSelectionCollection.met.push_back(&met);
	TauTightIsoObjectSelectionCollection.met.push_back(&met);

	
	//Filling count plot
	count->Fill("NoCuts",1.);

	//------------------------------------------//
	//------- EVENT SELECTION START ------------//
	//------------------------------------------//
	
	myHistoColl_baselineSelection.h_count->Fill("NoCuts",1.);	
	myHistoColl_TauLooseIsoObjectSelection.h_count->Fill("NoCuts",1.);
	myHistoColl_TauMediumIsoObjectSelection.h_count->Fill("NoCuts",1.);
	myHistoColl_Tau1TightIsoObjectSelection.h_count->Fill("NoCuts",1.);
	myHistoColl_TauTightIsoObjectSelection.h_count->Fill("NoCuts",1.);
	
	//Filling Histograms for baseline selection
	fillHistoCollection ( myHistoColl_baselineSelection, baselineObjectSelectionCollection,1.);

	
	  //Tau requirements and plot fillings
	if(!((int)(TauLooseIsoObjectSelectionCollection.tau.size() >= 2))){			//check if there is at least min taus in the event
			myHistoColl_TauLooseIsoObjectSelection.h_count->Fill("AtLeast2taus",1.);	
			fillHistoCollection ( myHistoColl_TauLooseIsoObjectSelection, TauLooseIsoObjectSelectionCollection,1.);
	    }
	if(!((int)(TauMediumIsoObjectSelectionCollection.tau.size() >= 2))){			//check if there is at least min taus in the event
			myHistoColl_TauMediumIsoObjectSelection.h_count->Fill("AtLeast2taus",1.);	
			fillHistoCollection ( myHistoColl_TauMediumIsoObjectSelection, TauMediumIsoObjectSelectionCollection,1.);
	    }
	if(!((int)(Tau1TightIsoObjectSelectionCollection.tau.size() >= 2))){			//check if there is at least min taus in the event
			myHistoColl_Tau1TightIsoObjectSelection.h_count->Fill("AtLeast2taus",1.);	
			fillHistoCollection ( myHistoColl_Tau1TightIsoObjectSelection, Tau1TightIsoObjectSelectionCollection,1.);
	    }
	if(!((int)(TauTightIsoObjectSelectionCollection.tau.size() >= 2))){			//check if there is at least min taus in the event
			myHistoColl_TauTightIsoObjectSelection.h_count->Fill("AtLeast2taus",1.);	
			fillHistoCollection ( myHistoColl_TauTightIsoObjectSelection, TauTightIsoObjectSelectionCollection,1.);
	    }


	//clearing event collections
	baselineObjectSelectionCollection.clear();
	TauLooseIsoObjectSelectionCollection.clear();
	TauMediumIsoObjectSelectionCollection.clear();
	Tau1TightIsoObjectSelectionCollection.clear();
	TauTightIsoObjectSelectionCollection.clear();

}



// ------------ method called once each job just before starting event loop  ------------
	void 
VBFSUSYanalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
VBFSUSYanalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   VBFSUSYanalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   VBFSUSYanalyzer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   VBFSUSYanalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   VBFSUSYanalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VBFSUSYanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VBFSUSYanalyzer);
