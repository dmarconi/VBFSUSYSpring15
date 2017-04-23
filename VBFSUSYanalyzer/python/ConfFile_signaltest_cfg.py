import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'file:///nfs/dust/cms/user/dmarconi/workdir/gridcontrol_makeSamples/miniaodv2.root'
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_5.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_50.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_51.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_52.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_53.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_54.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_55.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_56.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_57.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_58.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_59.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_6.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_60.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_61.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_62.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_63.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_64.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_65.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_66.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_67.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_68.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_69.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_7.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_70.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_71.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_72.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_73.root'
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_51.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_52.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_53.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_54.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_55.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_56.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_57.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_58.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_59.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_6.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_60.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_61.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_62.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_63.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_64.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_65.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_66.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_67.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_68.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_69.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_7.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_70.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_71.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_72.root',
#'file:///nfs/dust/cms/user/dmarconi/workdir/samples/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_1M/MINIAODSIM/VBFC1pmN2_C1ToTau_N2ToTauTau_LSP000_Stau295_Chargino300_MINIAODSIM_73.root'
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
        taupt = cms.double(20.),
        eventweight = cms.double(1.),
		verbose = cms.bool(True),
		)

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('histodemo.root')
		)

process.p = cms.Path(process.demo)
