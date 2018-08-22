from CRABClient.UserUtilities import config, getUsernameFromSiteDB
import os

config = config()

config.General.requestName = 'QCD_jetTuples_pythia8_PUMoriond17_v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_cfg.py'
config.JobType.inputFiles = ['QGL_cmssw8020_v2.db']

config.Data.inputDataset = '/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIISummer16MiniAODv2-PUMoriond17_magnetOn_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag = 'QCD_jetTuples_pythia8_PUMoriond17_v2'

config.Site.storageSite = 'T2_FI_HIP'
