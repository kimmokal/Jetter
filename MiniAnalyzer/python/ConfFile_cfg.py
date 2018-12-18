import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("AK4jets")

# QG likelihood
process.load("Jetter.MiniAnalyzer.QGLikelihood_cfi")
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag("slimmedJets")
process.QGTagger.jetsLabel = cms.string("QGL_AK4PFchs")

# File service
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('Jetter.root')

# Load up the filelist
fileList = FileUtils.loadListFromFile('RunIISummer16MiniAODv2_QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8_PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1_filelist.txt')

process.source = cms.Source("PoolSource",
	## Use whole data set
	# fileNames = cms.untracked.vstring(*fileList)
	## Use just one file
	fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/08278E4E-E4EF-E611-8BD7-FA163E3ABA64.root')
)

process.AK4jets = cms.EDAnalyzer('MiniAnalyzer',
	## jet, PF and generator level collections ##
	jets = cms.InputTag("slimmedJets"),
	pfCands = cms.InputTag("packedPFCandidates"),
	genJets = cms.InputTag("slimmedGenJets"),
	genEventInfo = cms.InputTag("generator"),
	## good primary vertices ##
	vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
	confGoodVtxNdof = cms.double(4),
	confGoodVtxZ = cms.double(24),
	confGoodVtxRho = cms.double(2),
	## pileup and rhos ##
	pileupInfo = cms.InputTag('slimmedAddPileupInfo'),
	pfRhoAll = cms.InputTag('fixedGridRhoFastjetAll'),
	pfRhoCentral = cms.InputTag('fixedGridRhoFastjetCentral'),
	pfRhoCentralNeutral = cms.InputTag('fixedGridRhoFastjetCentralNeutral'),
	pfRhoCentralChargedPileUp = cms.InputTag('fixedGridRhoFastjetCentralChargedPileUp'),
)

# Choose how many events to process (-1 = all)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Report execution progress
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.p = cms.Path(process.QGTagger + process.AK4jets)
