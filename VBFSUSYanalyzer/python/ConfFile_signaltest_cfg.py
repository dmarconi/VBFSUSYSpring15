import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'file:/nfs/dust/cms/user/dmarconi/workdir/VBFsignalSamples_Madgraph/testcmssw/CMSSW_7_4_14/src/genfragment_GEN_SIM_RECOBEFMIX_DIGI_L1_L1Reco_RECO.root'
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
