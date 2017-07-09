// -*- C++ -*-
//
// Package:    VBFSUSYSpring15/VBFSUSYLtoTfactors
// Class:      VBFSUSYLtoTfactors
//
/**\class VBFSUSYLtoTfactors VBFSUSYLtoTfactors.cc VBFSUSYSpring15/VBFSUSYLtoTfactors/plugins/VBFSUSYLtoTfactors.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Daniele Marconi
//         Created:  Sun, 17 Jan 2016 14:11:37 GMT
//
//


//---------------------------------
// system include files
//---------------------------------

#include <memory>
#include <string>
#include <vector>
#include <utility>

//---------------------------------
// user include files
//---------------------------------

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

//---------------------------------
// -- Root
//---------------------------------

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

//---------------------------------
// useful structs
//---------------------------------

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

//---------------------------------
struct MyEventCollection_LtoT {
//---------------------------------

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

struct MyHistoCollection_LtoT {

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
	TH2F* h2_DiJetInvMass_vs_MET;

	//---------------------------------
	//Plot init
	//---------------------------------

	void init(const std::string & inputlabel) {

		label = inputlabel;
		TFileDirectory subDir = fs->mkdir( inputlabel );
		TH1::SetDefaultSumw2();

		h_count = subDir.make<TH1F>("counts", "", 1,0,1);
		h_count->SetBit(TH1::kCanRebin);
		h_count->SetStats(0);
		h_count->Fill("NoCuts", 0);
		h_count->Fill("4 Jet",0.);
		h_count->Fill("1 Tau Match", 0.);
		h_count->Fill("Tau Match also 1Tight", 0.);
		h_count->Fill("4 Jets", 0.);
		h_count->Fill("1 Tau", 0.);
		h_count->Fill("1 Fake Loose Tau",0.);
		h_count->Fill("1 Fake also Tight",0.);

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

		h2_DiJetInvMass_vs_DiJetDEta = new TH2F("h2_DiJetInvMass_vs_DiJetDEta","h2_DiJetInvMass_vs_DiJetDEta", 20, 0., 10., 40, 0., 2500.);
		h2_DiJetInvMass_vs_DiJetDEta->GetXaxis()->SetTitle("#Delta#eta^{jj}");
		h2_DiJetInvMass_vs_DiJetDEta->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
		h2_tau1pt_vs_tau2pt = new TH2F("h2_tau1pt_vs_tau2pt","correlation of first and second p_{T}^{#tau}", 50, 0., 500., 50, 0., 500.);
		h2_tau1pt_vs_tau2pt->GetXaxis()->SetTitle("p_{T}^{#tau 1}");
		h2_tau1pt_vs_tau2pt->GetYaxis()->SetTitle("p_{T}^{#tau 2}");
		h2_DiJetInvMass_vs_MET = subDir.make<TH2F>("h2_DiJetInvMass_vs_MET","h2_DiJetInvMass_vs_MET", 24, 0., 240., 40, 0., 2500.);
		h2_DiJetInvMass_vs_MET->GetYaxis()->SetTitle("M^{(jet,jet)} [GeV]");
		h2_DiJetInvMass_vs_MET->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");

	}
};

//---------------------------------
// class declaration
//---------------------------------

class VBFSUSYLtoTfactors : public edm::EDAnalyzer {
	public:
		explicit VBFSUSYLtoTfactors(const edm::ParameterSet&);
		~VBFSUSYLtoTfactors();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		void fillHistoCollection (MyHistoCollection_LtoT &, MyEventCollection_LtoT, double, bool verbose_);
		std::pair<unsigned int,unsigned int> LeadingJets(MyEventCollection_LtoT);
		std::pair<unsigned int,unsigned int> LeadingTaus(MyEventCollection_LtoT);
		double TauJetMinDistance(MyEventCollection_LtoT , const pat::Jet &);
		TauProperties Inv2tMassIndex(MyEventCollection_LtoT);
		MassAndIndex Inv2jMassIndex(MyEventCollection_LtoT);
		void makeSelection (MyHistoCollection_LtoT &, MyEventCollection_LtoT, double , bool verbose_);

		//---------------------------------
		// ---------event collections-----------------------------
		//---------------------------------

		MyEventCollection_LtoT baselineObjectSelectionCollection;
		MyEventCollection_LtoT TauAnyIsoObjectSelectionCollection;
		MyEventCollection_LtoT TauAnyIsoPlusNonesObjectSelectionCollection;

		//---------------------------------
		// ---------histograms-----------------------------
		//---------------------------------

		edm::Service<TFileService> fs;
		TH1F* count;
		MyHistoCollection_LtoT myHistoColl_baselineSelection;
		MyHistoCollection_LtoT myHistoColl_TauAnyIsoObjectSelection;
		MyHistoCollection_LtoT myHistoColl_TauAnyIsoPlusNonesObjectSelection;

		//---------------------------------
		// ----------member data ---------------------------
		//---------------------------------

		edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
		edm::EDGetTokenT<pat::MuonCollection> muonToken_;
		edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
		edm::EDGetTokenT<pat::TauCollection> tauToken_;
		edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
		edm::EDGetTokenT<pat::JetCollection> jetToken_;
		edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
		edm::EDGetTokenT<pat::METCollection> metToken_;

		bool verbose_;
		double weight_;
		double taupt_;

};

//---------------------------------
// constants, enums and typedefs
//---------------------------------

//---------------------------------
// static data member definitions
//---------------------------------

//---------------------------------
// constructors and destructor
//---------------------------------

VBFSUSYLtoTfactors::VBFSUSYLtoTfactors(const edm::ParameterSet& iConfig):

	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
	electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
	tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
	fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
	metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),

	verbose_(iConfig.getParameter<bool>("verbose")),
	weight_(iConfig.getParameter<double>("eventweight")),
	taupt_(iConfig.getParameter<double>("taupt")),
	met_(iConfig.getParameter<double>("met")),
	mjj_(iConfig.getParameter<double>("mjj")),
	jetpt_(iConfig.getParameter<double>("jetpt"))

{
	//---------------------------------
	//Event collection init
	//---------------------------------

	baselineObjectSelectionCollection.init("baselineObjectSelection");
	TauAnyIsoObjectSelectionCollection.init("TauAnyIsoObjectSelection");
	TauAnyIsoObjectSelectionCollection.init("TauAnyIsoPlusNonesObjectSelection");
	//---------------------------------
	//histogram initialization
	//---------------------------------

	count = fs->make<TH1F>("counts", "", 1,0,1);
	count->SetBit(TH1::kCanRebin);
	count->SetStats(0);
	count->Fill("NoCuts",0);
	myHistoColl_baselineSelection.init("baselineObjectSelection");
	myHistoColl_TauAnyIsoObjectSelection.init("TauAnyIsoObjectSelection");
	myHistoColl_TauAnyIsoPlusNonesObjectSelection.init("TauAnyIsoPlusNonesObjectSelection");

}


VBFSUSYLtoTfactors::~VBFSUSYLtoTfactors()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//---------------------------------
// member functions
//---------------------------------

//---------------------------------
// ------------ method called for each event  ------------
//---------------------------------

	void
VBFSUSYLtoTfactors::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	using namespace std;

	//CMS Handles initialization

	edm::Handle<reco::VertexCollection> vertices;
	iEvent.getByToken(vtxToken_, vertices);
	if (vertices->empty()) return; // skip the event if no PV found
	const reco::Vertex &PV = vertices->front();

	edm::Handle<pat::MuonCollection> muons;
	iEvent.getByToken(muonToken_, muons);
	if (verbose_){
		for (const pat::Muon &mu : *muons) {
			if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
			printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
					mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
			}
	}

	edm::Handle<pat::ElectronCollection> electrons;
	iEvent.getByToken(electronToken_, electrons);
	if (verbose_){
		for (const pat::Electron &el : *electrons) {
			if (el.pt() < 5) continue;
			printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), pass conv veto %d\n",
				el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.passConversionVeto());
			}
	}

	edm::Handle<pat::PhotonCollection> photons;
	iEvent.getByToken(photonToken_, photons);
	if (verbose_){
		for (const pat::Photon &pho : *photons) {
			if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
			printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
					pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
				}
	}

	edm::Handle<pat::TauCollection> taus;
	iEvent.getByToken(tauToken_, taus);
	if (verbose_){
		for (const pat::Tau &tau : *taus) {
			if (tau.pt() < 20) continue;
			printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
				tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
		}
	}

	edm::Handle<pat::JetCollection> jets;
	iEvent.getByToken(jetToken_, jets);
	if (verbose_){
		int ijet = 0;
		for (const pat::Jet &j : *jets) {
			if (j.pt() < 20) continue;
			printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n",
					j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")), j.userFloat("pileupJetId:fullDiscriminant"));
					if ((++ijet) == 1) { // for the first jet, let's print the leading constituents
					std::vector<reco::CandidatePtr> daus(j.daughterPtrVector());
					std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
					for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
						const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
						printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i2,cand.pt(),cand.dz(PV.position()),cand.pdgId());
					}
				}
			}
		}

	edm::Handle<pat::METCollection> mets;
	iEvent.getByToken(metToken_, mets);
	const pat::MET &met = mets->front();
	if (verbose_){
		printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
				met.pt(), met.phi(), met.sumEt(),
				met.genMET()->pt(),
				met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));
}


	//smart tau selection
	std::vector<const pat::Tau*> anyiso;
	std::vector<const pat::Tau*> anyisoplusnones;
	std::vector<const pat::Tau*> nones;

	for (const pat::Tau &tau : *taus) {
		if (verbose_) cout << "Tau object"<< endl;
		if(!(	fabs(tau.eta()) <= 2.1                              					)) continue;
		if (verbose_) cout << "Eta cut passed for tau" << endl;
		if (verbose_) cout << "Tau pt : " << tau.pt() << endl;
		if (verbose_) cout << "Tau eta : " << tau.eta() << endl;
		if (verbose_) cout << "Tau phi : " << tau.phi() << endl;
		if(!(       tau.pt() >= taupt_                                            				)) continue;
		if (verbose_) cout << "Pt cut passed for tau" << endl;
		if(!(       tau.tauID("againstElectronVLooseMVA5") > 0.5                				)) continue;
		if (verbose_) cout << "Electron veto cut passed for tau" << endl;
		if(!(       tau.tauID("againstMuonLoose3") > 0.5                        				)) continue;
		if (verbose_) cout << "Muon veto cut passed for tau" << endl;
		if(!(       tau.leadChargedHadrCand()->pt() >= 5.0                      				)) continue;
		if (verbose_) cout << "leadChargedHadrCand_pt cut passed for tau" << endl;
		if(!(       (tau.tauID("decayModeFindingNewDMs") > 0.5) && (tau.signalChargedHadrCands().size() < 4)	)) continue;

		baselineObjectSelectionCollection.tau.push_back(&tau);

		if( tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")  > 0.5) {
			anyiso.push_back(&tau);
		}
		else nones.push_back(&tau);

	}


	if(anyiso.size()>=1) {
		anyisoplusnones.insert(anyisoplusnones.end(), anyiso.begin(), anyiso.end());
		anyisoplusnones.insert(anyisoplusnones.end(), nones.begin(), nones.end());
		for(unsigned int t =0;t<anyiso.size();++t) {
			TauAnyIsoObjectSelectionCollection.tau.push_back(anyiso[t]);
		}
		for(unsigned int t =0;t<anyisoplusnones.size();++t){
			TauAnyIsoPlusNonesObjectSelectionCollection.tau.push_back(anyisoplusnones[t]);
		}
	}


	for(const pat::Jet &jet : *jets){
		if(!(      jet.pt() >= 30.                                                	)) continue;  // Original value 20
		if(!(      fabs(jet.eta()) <= 5.0                                          )) continue;

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
			baselineObjectSelectionCollection.jet.push_back(&jet);
			TauAnyIsoObjectSelectionCollection.jet.push_back(&jet);
			TauAnyIsoPlusNonesObjectSelectionCollection.jet.push_back(&jet);
		}

		//Filling bJet collection
		if(fabs(jet.eta()) <= 2.4 && jet.bDiscriminator("combinedSecondaryVertexBJetTags") /*jet[j].bDiscriminator_combinedSecondaryVertexBJetTags*/ > 0.244    ){
			baselineObjectSelectionCollection.bjet.push_back(&jet);
			TauAnyIsoObjectSelectionCollection.bjet.push_back(&jet);
			TauAnyIsoPlusNonesObjectSelectionCollection.bjet.push_back(&jet);
		}
	}

	//MET selection
	baselineObjectSelectionCollection.met.push_back(&met);
	TauAnyIsoObjectSelectionCollection.met.push_back(&met);
	TauAnyIsoPlusNonesObjectSelectionCollection.met.push_back(&met);


	//Filling count plot
	count->Fill("NoCuts",weight_);

	//------------------------------------------//
	//------- EVENT SELECTION START ------------//
	//------------------------------------------//

	//Filling Histograms for baseline selection

	makeSelection ( myHistoColl_baselineSelection,
									baselineObjectSelectionCollection,
									weight_,
									verbose_);
	makeSelection ( myHistoColl_TauAnyIsoObjectSelection,
									TauAnyIsoObjectSelectionCollection,
									weight_,
									verbose_);
	makeSelection ( myHistoColl_TauAnyIsoPlusNonesObjectSelection,
									TauAnyIsoPlusNonesObjectSelectionCollection,
									weight_,
									verbose_);

	//clearing event collections
	baselineObjectSelectionCollection.clear();
	TauAnyIsoObjectSelectionCollection.clear();
	TauAnyIsoPlusNonesObjectSelectionCollection.clear();

}



// ------------ method called once each job just before starting event loop  ------------
	void
VBFSUSYLtoTfactors::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void
VBFSUSYLtoTfactors::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   VBFSUSYLtoTfactors::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   VBFSUSYLtoTfactors::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   VBFSUSYLtoTfactors::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   VBFSUSYLtoTfactors::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

void
VBFSUSYLtoTfactors::fillHistoCollection (MyHistoCollection_LtoT &inputHistoCollection, MyEventCollection_LtoT inputEventCollection, double weight_, bool _verbose) {

	// ---------------------
	// -- fill histograms --
	// ---------------------

	using namespace std;


	//JETS

	//find index of leading jets
	std::pair<unsigned int,unsigned int> jetIndex=LeadingJets(inputEventCollection);

	//find max value for 2-jet-mass
	MassAndIndex Inv2j = Inv2jMassIndex(inputEventCollection);

	//define ht
	double ht_jets=0.;

	//JET SEL
	for (unsigned int j = 0;j<inputEventCollection.jet.size();++j){
		inputHistoCollection.h_jetpt->Fill(inputEventCollection.jet[j]->pt(),weight_); //fill jet-pt-histogram
		inputHistoCollection.h_jeteta->Fill(inputEventCollection.jet[j]->eta(),weight_); //fill jet-eta-histogram
		ht_jets+=inputEventCollection.jet[j]->pt(); //add up scalar pt to ht
	}

	//fill jet count
	inputHistoCollection.h_njet->Fill( (int)inputEventCollection.jet.size(),weight_ );
	if(verbose_)std::cout<<"Pass selection -> Fill njet="<<inputEventCollection.jet.size()<<", weight="<<weight_<<std::endl;

	//fill jet pt indizes
	if (jetIndex.first < 99999)  {
		inputHistoCollection.h_jet1pt->Fill(inputEventCollection.jet[jetIndex.first]->pt(),weight_);
		inputHistoCollection.h_jet1eta->Fill(inputEventCollection.jet[jetIndex.first]->eta(),weight_);
		if(verbose_)std::cout<<"Pass selection -> Fill jet1pt="<<inputEventCollection.jet[jetIndex.first]->pt()<<", weight="<<weight_<<std::endl;
		if(verbose_)std::cout<<"Pass selection -> Fill jet1eta="<<inputEventCollection.jet[jetIndex.first]->eta()<<", weight="<<weight_<<std::endl;
	}
	if (jetIndex.second < 99999) {
		inputHistoCollection.h_jet2pt->Fill(inputEventCollection.jet[jetIndex.second]->pt(),weight_);
		inputHistoCollection.h_jet2eta->Fill(inputEventCollection.jet[jetIndex.second]->eta(),weight_);
		if(verbose_)std::cout<<"Pass selection -> Fill jet2pt="<<inputEventCollection.jet[jetIndex.second]->pt()<<", weight="<<weight_<<std::endl;
		if(verbose_)std::cout<<"Pass selection -> Fill jet2eta="<<inputEventCollection.jet[jetIndex.second]->eta()<<", weight="<<weight_<<std::endl;
	}

	//fill 2-jet-event inv. mass and eta-difference
	if ( Inv2j.Mass >= 0 ) {
		inputHistoCollection.h_dijetinvariantmass ->Fill(Inv2j.Mass,weight_);
		inputHistoCollection.h_dijetdeltaeta ->Fill(Inv2j.dEta,weight_);
		if(verbose_)std::cout<<"Pass selection -> Fill dijetinvariantmass="<< Inv2j.Mass <<", weight="<<weight_<<std::endl;
		if(verbose_)std::cout<<"Pass selection -> Fill dijetdeltaeta="<< Inv2j.dEta <<", weight="<<weight_<<std::endl;
	}

	//fill ht distribution
	inputHistoCollection.h_ht -> Fill(ht_jets,weight_);
	if(verbose_)std::cout<<"Pass selection -> Fill ht_jets="<< ht_jets <<", weight="<<weight_<<std::endl;

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
		inputHistoCollection.h_ditauinvariantmass ->Fill(Inv2t.Mass,weight_);
		if(verbose_)std::cout<<"Pass selection -> Fill ditauinvariantmass="<< Inv2t.Mass <<", weight="<<weight_<<std::endl;

		//fill tau charge and  cosdeltaphi and deltaeta and 2Dpt-plot
		inputHistoCollection.h_ditaucharge ->Fill(Inv2t.charge,weight_);
		if(verbose_)std::cout<<"Pass selection -> Fill ditaucharge="<< Inv2t.charge <<", weight="<<weight_<<std::endl;
		inputHistoCollection.h_ditaucosdeltaphi ->Fill(Inv2t.cosDphi,weight_);
		if(verbose_)std::cout<<"Pass selection -> Fill ditaucosdeltaphi="<< Inv2t.cosDphi <<", weight="<<weight_<<std::endl;
		inputHistoCollection.h_ditaudeltaeta->Fill(Inv2t.dEta,weight_);
		if(verbose_)std::cout<<"Pass selection -> Fill ditaudeltaeta="<< Inv2t.dEta <<", weight="<<weight_<<std::endl;
		inputHistoCollection.h2_tau1pt_vs_tau2pt->Fill(inputEventCollection.tau[Inv2t.first]->pt(),inputEventCollection.tau[Inv2t.second]->pt(),weight_);
	}

	//fill tau pt and eta
	if (Inv2t.first < 99999) {
		inputHistoCollection.h_tau1pt->Fill(inputEventCollection.tau[Inv2t.first]->pt(),weight_);
		inputHistoCollection.h_tau1eta->Fill(inputEventCollection.tau[Inv2t.first]->eta(),weight_);
		if(verbose_)std::cout<<"Pass selection -> Fill tau1pt="<< inputEventCollection.tau[Inv2t.first]->pt() <<", weight="<<weight_<<std::endl;
		if(verbose_)std::cout<<"Pass selection -> Fill tau1eta="<< inputEventCollection.tau[Inv2t.first]->eta() <<", weight="<<weight_<<std::endl;
	}
	if (Inv2t.second < 99999) {
		inputHistoCollection.h_tau2pt->Fill(inputEventCollection.tau[Inv2t.second]->pt(),weight_);
		inputHistoCollection.h_tau2eta->Fill(inputEventCollection.tau[Inv2t.second]->eta(),weight_);
		if(verbose_)std::cout<<"Pass selection -> Fill tau2pt="<< inputEventCollection.tau[Inv2t.second]->pt() <<", weight="<<weight_<<std::endl;
		if(verbose_)std::cout<<"Pass selection -> Fill tau2eta="<< inputEventCollection.tau[Inv2t.second]->eta() <<", weight="<<weight_<<std::endl;
	}

	//fill ht with taus included
	inputHistoCollection.h_ht_withtau -> Fill(ht_jetsPtau,weight_);
	if(verbose_)std::cout<<"Pass selection -> Fill ht_withtau="<< ht_jetsPtau <<", weight="<<weight_<<std::endl;

	// MET
	inputHistoCollection.h_met -> Fill(inputEventCollection.met[0]->pt(),weight_);
	if(verbose_)std::cout<<"Pass selection -> Fill met="<< inputEventCollection.met[0]->pt() <<", weight="<<weight_<<std::endl;

	//fill DiJetInvMass_vs_DiJetDEta
	inputHistoCollection.h2_DiJetInvMass_vs_DiJetDEta -> Fill(Inv2j.dEta, Inv2j.Mass,weight_);

	//fill DiJetInvMass_vs_MET
	inputHistoCollection.h2_DiJetInvMass_vs_MET -> Fill(inputEventCollection.met[0]->pt(), Inv2j.Mass,weight_);
	//________________________________________

}

//-----------------
//leading jet finder
//_________________

std::pair<unsigned int,unsigned int> VBFSUSYLtoTfactors::LeadingJets(MyEventCollection_LtoT collection)
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

std::pair<unsigned int,unsigned int> VBFSUSYLtoTfactors::LeadingTaus(MyEventCollection_LtoT collection)
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

double TauJetMinDistance(MyEventCollection_LtoT collection, const pat::Jet &jet)
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

MassAndIndex VBFSUSYLtoTfactors::Inv2jMassIndex(MyEventCollection_LtoT collection)
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

TauProperties VBFSUSYLtoTfactors::Inv2tMassIndex(MyEventCollection_LtoT collection)
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

void VBFSUSYLtoTfactors::makeSelection (MyHistoCollection_LtoT &inputHistoCollection, MyEventCollection_LtoT inputEventCollection, double weight_, bool verbose_) {

	inputHistoCollection.h_count->Fill("NoCuts",weight_);

	bool oneJet = false;
	bool foundTauMatch = false;
	bool isAlsoTightTau = false;

	if (inputEventCollection.jet.size() >= 4) oneJet = true;

	for(unsigned int j=1; j<inputEventCollection.jet.size(); j++){

		TLorentzVector v_tau;
		TLorentzVector v_jet;

		for(unsigned int t=1; t<inputEventCollection.tau.size(); t++){

			v_tau.SetPtEtaPhiM(inputEventCollection.tau[t]->pt(),inputEventCollection.tau[t]->eta(),inputEventCollection.tau[t]->phi(),inputEventCollection.tau[t]->mass());
			v_jet.SetPtEtaPhiM(inputEventCollection.jet[j]->pt(),inputEventCollection.jet[j]->eta(),inputEventCollection.jet[j]->phi(),inputEventCollection.jet[j]->mass());
			double temp_deltaRtaujet = v_tau.DeltaR(v_jet);
			if ( temp_deltaRtaujet < 0.3 ) {
				foundTauMatch = true;
				if( inputEventCollection.tau[t]->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")  > 0.5) isAlsoTightTau = true;
			}
		}
	}


	if(oneJet) inputHistoCollection.h_count->Fill("4 Jet", weight_);
	if(foundTauMatch) inputHistoCollection.h_count->Fill("1 Tau Match", weight_);
	if(isAlsoTightTau) inputHistoCollection.h_count->Fill("Tau Match also 1Tight", weight_);

	//check if there are 4 jets in the event
	if((int)(inputEventCollection.jet.size() >= 4)){
		inputHistoCollection.h_count->Fill("4 Jets",weight_);
	} else return;

	//check if there are 4 jets in the event
	if((int)(inputEventCollection.tau.size() >= 1)){
		inputHistoCollection.h_count->Fill("1 Tau",weight_);
	} else return;

	bool foundLooseMatch = false;
	bool isAlsoTight = false;

	for(unsigned int t=1; t<inputEventCollection.tau.size(); t++){

		TLorentzVector v_tau;
		TLorentzVector v_jet;
		for(unsigned int j=1; j<inputEventCollection.jet.size(); j++){
			v_tau.SetPtEtaPhiM(inputEventCollection.tau[t]->pt(),inputEventCollection.tau[t]->eta(),inputEventCollection.tau[t]->phi(),inputEventCollection.tau[t]->mass());
			v_jet.SetPtEtaPhiM(inputEventCollection.jet[j]->pt(),inputEventCollection.jet[j]->eta(),inputEventCollection.jet[j]->phi(),inputEventCollection.jet[j]->mass());
			double temp_deltaRtaujet = v_tau.DeltaR(v_jet);
			if ( temp_deltaRtaujet < 0.3 ) {
				if( inputEventCollection.tau[t]->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5) foundLooseMatch = true;
				if( inputEventCollection.tau[t]->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")  > 0.5) isAlsoTight = true;
			}
		}

	}

	if(foundLooseMatch){
		inputHistoCollection.h_count->Fill("1 Fake Loose Tau",weight_);
	} else return;

	if(isAlsoTight){
		inputHistoCollection.h_count->Fill("1 Fake also Tight",weight_);
	} else return;

	fillHistoCollection ( inputHistoCollection,
												inputEventCollection,
												weight_,
												verbose_
											);

	return;

}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VBFSUSYLtoTfactors::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VBFSUSYLtoTfactors);
