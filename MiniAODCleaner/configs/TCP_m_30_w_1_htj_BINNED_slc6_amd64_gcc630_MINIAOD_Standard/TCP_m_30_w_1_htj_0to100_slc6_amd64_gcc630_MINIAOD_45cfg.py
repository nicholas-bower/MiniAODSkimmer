
import FWCore.ParameterSet.Config as cms
######
# Configuration to run tau ReReco+PAT at MiniAOD samples
# M. Bluj, NCBJ Warsaw
# based on work of J. Steggemann, CERN
# Created: 9 Nov. 2017
#With additional implementation for Muon/Electron Cleaning from Jets
# Redwan Md Habibullah 8 July 2021
######

######

import PhysicsTools.PatAlgos.tools.helpers as configtools
from PhysicsTools.PatAlgos.tools.helpers import cloneProcessingSnippet
from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag
from FWCore.ParameterSet.MassReplace import massSearchReplaceParam


#runSignal = True
runSignal=False
#maxEvents = 1000
maxEvents=-1
appendOutput = True
isMC = True
# If 'reclusterJets' set true a new collection of uncorrected ak4PFJets is
# built to seed taus (as at RECO), otherwise standard slimmedJets are used
reclusterJets = True
# reclusterJets = False

# set true for upgrade studies
phase2 = False
# phase2 = True

# Output mode
outMode = 2  # store original MiniAOD and new selectedPatTaus
# outMode = 1 #store original MiniAOD, new selectedPatTaus, and all PFtau products as in AOD (except of unsuported ones)


print('Running Tau reco&id with MiniAOD inputs:')
print('	 Run on signal:', runSignal)
print('	 Recluster jets:', reclusterJets)
print('	 Use Phase2 settings:', phase2)
print('	 Output mode:', outMode)

#####
#from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Era_Run2_2017_cff import Run2_2017
from Configuration.Eras.Modifier_run2_miniAOD_devel_cff import run2_miniAOD_devel
from Configuration.ProcessModifiers.run2_miniAOD_UL_cff import run2_miniAOD_UL  
#era = Run2_2018
era = Run2_2017,run2_miniAOD_devel,run2_miniAOD_UL
if phase2:
    from Configuration.Eras.Era_Phase2_timing_cff import Phase2_timing
    era = Phase2_timing
process = cms.Process("TAURECO",  Run2_2017, run2_miniAOD_devel,run2_miniAOD_UL)
# for CH reco
process.load("Configuration.StandardSequences.MagneticField_cff")
if not phase2:
    process.load("Configuration.Geometry.GeometryRecoDB_cff")
else:
    process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')

#####
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource", fileNames=readFiles, secondaryFileNames=secFiles)

process.maxEvents = cms.untracked.PSet(
    input=cms.untracked.int32(maxEvents)
)
print('	 Max events:', process.maxEvents.input.value())

if runSignal:
    readFiles.extend([
        #'file:patMiniAOD_standard.root'
       'root://cmseos.fnal.gov//store/user/nbower/Events/TCP_m_30_w_1_htj_0to100_slc6_amd64_gcc630_MINIAOD/TCP_m_30_w_1_htj_0to100_slc6_amd64_gcc630_MINIAOD_45.root'
    ])
else:
    readFiles.extend([
        #'file:patMiniAOD_standard.root'
        'root://cmseos.fnal.gov//store/user/nbower/Events/TCP_m_30_w_1_htj_0to100_slc6_amd64_gcc630_MINIAOD/TCP_m_30_w_1_htj_0to100_slc6_amd64_gcc630_MINIAOD_45.root'
        #'/store/relval/CMSSW_10_5_0_pre1/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_103X_mcRun2_asymptotic_v3-v1/20000/A5CBC261-E3AB-C842-896F-E6AFB38DD22F.root'
    ])

#####
#import RecoTauTag.Configuration.tools.adaptToRunAtMiniAOD as tauAtMiniTools
import MiniAODSkimmer.MiniAODCleaner.adaptToRunAtMiniAODCustom as tauAtMiniToolsCustom

#####
print "Step : 1 - Added Paths for RecoCleaned "

tauAtMiniToolsCustom.addTauReRecoCustom(process)



#####
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if not phase2:
    process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')
else:
    process.GlobalTag = GlobalTag(
        process.GlobalTag, 'auto:phase2_realistic', '')

#####
# mode = 0: store original MiniAOD and new selectedPatTaus
# mode = 1: store original MiniAOD, new selectedPatTaus, and all PFtau products as in AOD (except of unsuported ones)
print "Step : 2 - Declare Outputs"

process.output = tauAtMiniToolsCustom.setOutputModule(mode=outMode)
if runSignal:
    process.output.fileName = 'root://cmseos.fnal.gov//store/user/nbower/Events/TCP_m_30_w_1_htj_BINNED_slc6_amd64_gcc630_MINIAOD_Skimmed_Standard/0to100//TCP_m_30_w_1_htj_0to100_slc6_amd64_gcc630_MINIAOD_45_SKIMMEDLowPT.root'
    if reclusterJets:
        process.output.fileName = 'root://cmseos.fnal.gov//store/user/nbower/Events/TCP_m_30_w_1_htj_BINNED_slc6_amd64_gcc630_MINIAOD_Skimmed_Standard/0to100//TCP_m_30_w_1_htj_0to100_slc6_amd64_gcc630_MINIAOD_45_SKIMMEDLowPT.root'
else:
    process.output.fileName = 'root://cmseos.fnal.gov//store/user/nbower/Events/TCP_m_30_w_1_htj_BINNED_slc6_amd64_gcc630_MINIAOD_Skimmed_Standard/0to100//TCP_m_30_w_1_htj_0to100_slc6_amd64_gcc630_MINIAOD_45_SKIMMEDLowPT.root'
    if reclusterJets:
        process.output.fileName = 'root://cmseos.fnal.gov//store/user/nbower/Events/TCP_m_30_w_1_htj_BINNED_slc6_amd64_gcc630_MINIAOD_Skimmed_Standard/0to100//TCP_m_30_w_1_htj_0to100_slc6_amd64_gcc630_MINIAOD_45_SKIMMEDLowPT.root'
##### Modify ouput by Hand#####

if appendOutput:
    #process.output.outputCommands.append('keep *_selectedPatTaus_*_*')
    #process.output.outputCommands.append('keep *_selectedPatTausElectronCleaned_*_*')
    #process.output.outputCommands.append('keep *_selectedPatTausMuonCleaned_*_*')
    process.output.outputCommands.append('keep *_slimmedTausUnCleaned_*_*')
    process.output.outputCommands.append('keep *_slimmedTausElectronCleaned_*_*')
    process.output.outputCommands.append('keep *_slimmedTausMuonCleaned_*_*')
    process.output.outputCommands.append('keep *_lumiSummary_*_*')
    
 
#####

tauAtMiniToolsCustom.addTauReRecoCustom(process)



#process.out = cms.EndPath(process.output)
#process.schedule.append(process.out)
#####
print "Step : 3 - Adapt Tau Reco to MiniAOD inputs"

tauAtMiniToolsCustom.adaptTauToMiniAODReReco(process, reclusterJets)
print "Step : 4 - Lower Pt Standard Taus"

###### lowering Pt of Standard Taus ######
minJetPt = 5

process.ak4PFJetsLegacyHPSPiZeros.minJetPt = minJetPt
process.combinatoricRecoTaus.minJetPt = minJetPt
process.recoTauAK4Jets08RegionPAT.minJetPt = minJetPt
process.ak4PFJetsRecoTauChargedHadrons.minJetPt = minJetPt
process.selectedPatTaus.cut = cms.string('pt > 8.0 && abs(eta)<2.3 && tauID(\'decayModeFinding\')> 0.5')
##########################################
#### Lower Tau Pt ElectronCleaned Taus###########

print "Step : 5 - Lower Pt ElectronCleaned Taus"

jetPt=5
tauPt=8

getattr(process,'selectedPatTausElectronCleaned').cut = cms.string("pt > {} && abs(eta) < 2.3 && tauID('decayModeFinding')> 0.5".format(tauPt))
process.ak4PFJetsLegacyHPSPiZerosElectronCleaned.minJetPt = jetPt
process.combinatoricRecoTausElectronCleaned.minJetPt = jetPt
process.recoTauAK4Jets08RegionPATElectronCleaned.minJetPt = jetPt
process.ak4PFJetsRecoTauChargedHadronsElectronCleaned.minJetPt = jetPt

##########################################
#### Lower Tau Pt MuonCleaned Taus###########

print "Step : 6 - Lower Pt MuonCleaned Taus"

getattr(process,'selectedPatTausMuonCleaned').cut = cms.string("pt > {} && abs(eta) < 2.3 && tauID('decayModeFinding')> 0.5".format(tauPt))
process.ak4PFJetsLegacyHPSPiZerosMuonCleaned.minJetPt = jetPt
process.combinatoricRecoTausMuonCleaned.minJetPt = jetPt
process.recoTauAK4Jets08RegionPATMuonCleaned.minJetPt = jetPt
process.ak4PFJetsRecoTauChargedHadronsMuonCleaned.minJetPt = jetPt


tauAtMiniToolsCustom.addFurtherSkimming(process)

process.out = cms.EndPath(process.output)
process.schedule.append(process.out)

###########################################
process.load('FWCore.MessageService.MessageLogger_cfi')
if process.maxEvents.input.value() > 10:
    process.MessageLogger.cerr.FwkReport.reportEvery = process.maxEvents.input.value()//10
if process.maxEvents.input.value() > 10000 or process.maxEvents.input.value() < 0:
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#####
process.options = cms.untracked.PSet(
)
#process.options.numberOfThreads = cms.untracked.uint32(4)
process.options.numberOfThreads=cms.untracked.uint32(1)
process.options.numberOfStreams = cms.untracked.uint32(0)
print('	 No. of threads:', process.options.numberOfThreads.value(), ', no. of streams:', process.options.numberOfStreams.value())

process.options = cms.untracked.PSet(
    process.options,
    wantSummary=cms.untracked.bool(True)
)

dump_file = open('dump_rerunMiniAODClean.py','w')
dump_file.write(process.dumpPython())


