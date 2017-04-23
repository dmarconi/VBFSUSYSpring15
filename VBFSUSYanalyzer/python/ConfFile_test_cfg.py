import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/006E8D1C-7D6F-E511-B8DA-02163E012ED1.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/02346408-706F-E511-8C19-0025905C2C86.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/04F73055-7D6F-E511-BBB3-02163E00B0E6.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/08421225-7D6F-E511-B333-0CC47A04CFF6.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/10D0D759-7E6F-E511-AB80-02163E00C49A.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/10F0DC4D-5B6F-E511-A964-7845C4FC363E.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/1CDED031-706F-E511-AD59-0025905C95F8.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/2040F543-5B6F-E511-BAC5-00266CFAE074.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/280D9A48-7D6F-E511-B0A6-02163E016963.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/2867B93B-706F-E511-B0B0-0025905C4300.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/2882200F-706F-E511-BD06-0025905C3D96.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/32592BE1-706F-E511-A51B-0025904C6626.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/32AFA923-7D6F-E511-A438-02163E00B7E1.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/344F00F1-706F-E511-BD3A-0025905C3E66.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/3813187C-7E6F-E511-A669-02163E015DD4.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/3E7452F8-6F6F-E511-AFD8-0025905C4264.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/46F3353A-7D6F-E511-A20B-001E67C9AF38.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/4848CE07-706F-E511-84A7-0025905C42F2.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/48560C2F-706F-E511-9E40-0025905C2CEA.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/4896AFE6-706F-E511-B258-0025904C66F4.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/4A901A63-7E6F-E511-96BA-02163E00F513.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/4E444202-7D6F-E511-A886-003048F0E2C2.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/52428632-706F-E511-AAD4-0025905C42A8.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/56B7D4A2-7E6F-E511-B55D-02163E016747.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/5A67CB3F-EC6E-E511-9DE4-0025905A48D0.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/5E5E933F-5B6F-E511-BD61-00266CF9B318.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/6075686F-7E6F-E511-A5B6-02163E013ED1.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/60B3ED32-706F-E511-A6CA-0025904CF758.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/66906B6A-7D6F-E511-9DBC-02163E00BDE7.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/68CA39EF-706F-E511-8DF9-0025905C3D98.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/6A13DBA5-7D6F-E511-B206-02163E00E5C7.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/6A4B00A2-7E6F-E511-9859-02163E00E72D.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/6CF4EE13-7D6F-E511-9542-02163E010D48.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/6EB822FF-6F6F-E511-8DDE-0025904C66EC.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/729E593A-EC6E-E511-8874-0025905B85D0.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/7A9E1C6E-5B6F-E511-9FE6-00266CFAEA48.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/7E258335-706F-E511-A6F2-0025905C975E.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/7ED3F816-706F-E511-AFA0-0025905C4300.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/8A8A7452-7D6F-E511-AC9F-02163E00CD84.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/8C670326-7D6F-E511-9349-02163E016AD7.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/903FFD99-7E6F-E511-85EA-02163E016A8E.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/924DAE88-7E6F-E511-A0B0-02163E013BEA.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/9C410E35-5B6F-E511-B5F5-3417EBE644C2.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/A668F990-7E6F-E511-AFFA-02163E0166D5.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/AA879BE3-706F-E511-886E-0025904C678A.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/BC23E90F-5D6F-E511-93B1-008CFA007B98.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/C4E7D36F-7D6F-E511-9C2D-02163E01678E.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/D03ABB1F-7D6F-E511-B83A-02163E013D28.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/D09EA6A6-7E6F-E511-9019-02163E012DD3.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/D2D78013-7F6F-E511-BF93-02163E01679C.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/D4898277-7E6F-E511-82F5-001E675793A8.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/D4DB358B-7E6F-E511-BABE-003048FEAFB4.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/DA2F8DBB-7E6F-E511-9542-02163E0152D0.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/DA5216AB-7E6F-E511-85AA-02163E016952.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/E6F81A3C-EC6E-E511-916B-0025905A6084.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/EEE8C895-7D6F-E511-A7A8-0025904B2C68.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/F054743B-7D6F-E511-A4A2-02163E0169AB.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/F0F48A74-7E6F-E511-963C-002590494D1A.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/F8189895-7E6F-E511-8040-02163E016745.root',
'file:/pnfs/desy.de/cms/tier2/store/mc/RunIISpring15MiniAODv2/DYBBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/F85D88B0-7E6F-E511-8231-02163E00CE00.root'
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
        taupt = cms.double(20),
        verbose = cms.bool(False),

		)

process.TFileService = cms.Service("TFileService",
		fileName = cms.string('histodemo.root')
		)

process.p = cms.Path(process.demo)
