import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import sys

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

outfile = sys.argv[3]

filelist = open(sys.argv[2]).readlines()

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring( filelist )
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
		taupt = cms.double(50.),
		met = cms.double(0.),
		mjj = cms.double(0.),
		jetpt = cms.double(30.),
		eventweight = cms.double(1.),
		verbose = cms.bool(False),
)

process.TFileService = cms.Service("TFileService",
					fileName = cms.string(outfile)
				)

process.p = cms.Path(process.demo)
