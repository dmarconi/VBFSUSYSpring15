import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/0CD21733-437A-E511-97BF-44A8423C4026.root'
    )
)

process.demo = cms.EDAnalyzer('VBFSUSYanalyzer',
		vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
		muons = cms.InputTag("slimmedMuons"),
		electrons = cms.InputTag("slimmedElectrons"),
		taus = cms.InputTag("slimmedTaus"),
		photons = cms.InputTag("slimmedPhotons"),
		jets = cms.InputTag("slimmedJets"),
		fatjets = cms.InputTag("slimmedJetsAK8"),
		mets = cms.InputTag("slimmedMETs"),
)

process.TFileService = cms.Service("TFileService", 
					fileName = cms.string('histodemo.root')
				)

process.p = cms.Path(process.demo)
