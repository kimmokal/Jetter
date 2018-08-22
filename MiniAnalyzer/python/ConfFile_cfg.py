import FWCore.ParameterSet.Config as cms
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.Utilities.FileUtils as FileUtils
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

#process.source.lumisToProcess = LumiList.LumiList(filename = 'goodList.json').getVLuminosityBlockRange()

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

#QG likelihood
process.load("Jetter.MiniAnalyzer.QGLikelihood_cfi")
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets=cms.InputTag("slimmedJets")
process.QGTagger.jetsLabel = cms.string("QGL_AK4PFchs")

#File service
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('Jetter.root')

#Choose how many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Select the MiniAOD file to process
# With pileup:
# fileList = FileUtils.loadListFromFile('RunIISummer16MiniAODv2_QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_filelist.txt')
# Without pileup:
# fileList = FileUtils.loadListFromFile('RunIISummer16MiniAODv2_QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_filelist.txt')

#'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/08278E4E-E4EF-E611-8BD7-FA163E3ABA64.root'
#'/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/NoPU_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/DCD55D97-61F0-E611-9C1B-FA163EAD13C1.root'

process.source = cms.Source("PoolSource",
	# fileNames = cms.untracked.vstring(*fileList)
	fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/08278E4E-E4EF-E611-8BD7-FA163E3ABA64.root')
)

process.demo = cms.EDAnalyzer('MiniAnalyzer',
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    mets = cms.InputTag("slimmedMETs"),
    pfCands = cms.InputTag("packedPFCandidates"),
    packed = cms.InputTag("packedGenParticles"),
    genJets = cms.InputTag("slimmedGenJets"),
	genEventInfo = cms.InputTag("generator")

)


process.p = cms.Path(process.QGTagger + process.demo)

process.MessageLogger.cerr.FwkReport.reportEvery = 10
