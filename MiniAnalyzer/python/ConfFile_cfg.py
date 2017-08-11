import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#Choose how many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#Select the MiniAOD file to process
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring(
	'/store/mc/RunIISummer16MiniAODv2/AToZhToLLTauTau_M-340_13TeV_madgraph_4f_LO/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/80000/1C5025FA-9ABE-E611-B244-B083FECF83AB.root'
	)
)

process.demo = cms.EDAnalyzer('MiniAnalyzer',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
)


process.p = cms.Path(process.demo)
