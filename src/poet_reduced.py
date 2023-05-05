import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
import sys

#---- sys.argv takes the parameters given as input cmsRun PhysObjectExtractor/python/poet_cfg.py <isData (default=False)>
#----  e.g: cmsRun PhysObjectExtractor/python/poet_cfg.py True
#---- NB the first two parameters are always "cmsRun" and the config file name
#---- Work with data (if False, assumed MC simulations)
#---- This needs to be in agreement with the input files/datasets below.
if len(sys.argv) > 2:
    isData = eval(sys.argv[2])
else:
    isData = False
isMC = True
if isData: isMC = False

process = cms.Process("POET")

#---- Configure the framework messaging system
#---- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"
process.MessageLogger.categories.append("POET")
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit=cms.untracked.int32(-1))
process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True))

#---- Select the maximum number of events to process (if -1, run over all events)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

#---- Needed configuration for dealing with transient tracks if required
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#---- Define the test source files to be read using the xrootd protocol (root://), or local files (file:)
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
        #'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/00DF0A73-17C2-E511-B086-E41D2D08DE30.root'   
        'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000/003C3C2D-06EA-E511-8381-0023AEEEB799.root'
	#'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000/00940DCD-0DEC-E511-BF97-34E6D7E38781.root'
	#'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/00000/00C5EB13-87EA-E511-B076-BC305B390A32.root'
	#'root://eospublic.cern.ch//eos/opendata/cms/mc/RunIIFall15MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext3-v1/00000/020D0AF1-4BC2-E511-BDFC-0026B95ADB18.root'


        )
)
if isData:
    process.source.fileNames = cms.untracked.vstring(
        'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/SingleMuon/MINIAOD/16Dec2015-v1/10000/00006301-CAA8-E511-AD39-549F35AD8BC9.root'
#	'root://eospublic.cern.ch//eos/opendata/cms/Run2015D/SingleElectron/MINIAOD/08Jun2016-v1/10000/001A703B-B52E-E611-BA13-0025905A60B6.root'
        )

    #---- Apply the data quality JSON file filter. This example is for 2015 data
    #---- It needs to be done after the process.source definition
    #---- Make sure the location of the file agrees with your setup
    goodJSON = "data/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt"
    myLumis = LumiList.LumiList(filename=goodJSON).getCMSSWString().split(",")
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)


#---- These two lines are needed if you require access to the conditions database. E.g., to get jet energy corrections, trigger prescales, etc.
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#---- Uncomment and arrange a line like this if you are getting access to the conditions database through CVMFS snapshot files (requires installing CVMFS client)
#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/76X_dataRun2_16Dec2015_v0.db')
#---- If the container has local DB files available, uncomment lines like the ones below instead of the corresponding lines above
if isData: process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/76X_dataRun2_16Dec2015_v0.db')
else: process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/76X_mcRun2_asymptotic_RunIIFall15DR76_v1.db')
#---- The global tag must correspond to the needed epoch (comment out if no conditions needed)
if isData: process.GlobalTag.globaltag = '76X_dataRun2_16Dec2015_v0'
else: process.GlobalTag.globaltag = "76X_mcRun2_asymptotic_RunIIFall15DR76_v1"

#----- Configure POET analyzers -----#

process.mymuons = cms.EDAnalyzer('MuonAnalyzer', 
                                 muons = cms.InputTag("slimmedMuons"), 
                                 vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
				 pruned=cms.InputTag("prunedGenParticles"),
				 input_particle = cms.vstring("1:13"))

#---- Module to store trigger objects (functional but not fully developed yet) -------#
#process.mytrigobjs = cms.EDAnalyzer('TriggObjectAnalyzer', objects = cms.InputTag("selectedPatTrigger"))

#---- Example on how to add trigger information
#---- To include it, uncomment the lines below and include the
#---- module in the final path
#process.mytriggers = cms.EDAnalyzer('TriggerAnalyzer',
#                              processName = cms.string("HLT"),
#                              #---- These are example of OR of triggers for 2015
#                              #---- Wildcards * and ? are accepted (with usual meanings)
#                              #---- If left empty, all triggers will run              
#                              triggerPatterns = cms.vstring("HLT_IsoMu20_v*","HLT_IsoTkMu20_v*"), 
#                              triggerResults = cms.InputTag("TriggerResults","","HLT")
#                              )


#------------Example of simple trigger module with parameters by hand-------------------#


#----------- Turn on a trigger filter by adding this module to the the final path below -------#


#process.mygenparticle = cms.EDAnalyzer('GenParticleAnalyzer', 
#                                       pruned=cms.InputTag("prunedGenParticles"),
#                                       muons = cms.InputTag("slimmedMuons"),
#                                       vertices=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                       #---- Collect particles with specific "status:pdgid"
                                       #---- if 0:0, collect them all 
#                                       input_particle = cms.vstring("1:13"))

#----- Begin Jet correction setup -----#
JecString = 'MC'
if isData: JecString = 'DATA'

from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import updatedPatJets

#----- Apply the noise jet ID filter -----#
process.looseAK4Jets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                    filterParams = pfJetIDSelector.clone(),
                                    src = cms.InputTag("slimmedJets"))
process.looseAK8Jets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                                    filterParams = pfJetIDSelector.clone(),
                                    src = cms.InputTag("slimmedJetsAK8"))

#----- Apply the final jet energy corrections for 2015 -----#
process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(src = cms.InputTag("looseAK4Jets"))
if isData: process.patJetCorrFactorsReapplyJEC.levels.append('L2L3Residual')
process.slimmedJetsNewJEC = updatedPatJets.clone(
    jetSource = cms.InputTag("looseAK4Jets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC")),
)
process.patJetCorrFactorsReapplyJECAK8 = updatedPatJetCorrFactors.clone(
        src = cms.InputTag("looseAK8Jets"),
        levels = ['L1FastJet', 'L2Relative', 'L3Absolute'],
        payload = 'AK8PFchs'
        )
if isData: process.patJetCorrFactorsReapplyJECAK8.levels.append('L2L3Residual')
process.slimmedJetsAK8NewJEC = updatedPatJets.clone(
    jetSource = cms.InputTag("looseAK8Jets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK8")),
)

#----- Configure the POET jet analyzers -----#
process.myjets = cms.EDAnalyzer('JetAnalyzer', 
				jets = cms.InputTag("slimmedJetsNewJEC"),
				isData = cms.bool(isData),
				jetJECUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_Uncertainty_AK4PFchs.txt'), 
                                jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_PtResolution_AK4PFchs.txt'),
                                jerSFName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_SF_AK4PFchs.txt'),
				)
#process.myfatjets = cms.EDAnalyzer('FatjetAnalyzer', 
#				fatjets = cms.InputTag("slimmedJetsAK8NewJEC"),
#				isData = cms.bool(isData),
#                                jecL2Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L2Relative_AK8PFchs.txt'), 
#                                jecL3Name = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_'+JecString+'_L3Absolute_AK8PFchs.txt'), 
#                                jecResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_DATA_L2L3Residual_AK8PFchs.txt'), 
#				jetJECUncName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_Uncertainty_AK8PFchs.txt'), 
#                                jerResName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_PtResolution_AK8PFchs.txt'),
#                                jerSFName = cms.FileInPath('PhysObjectExtractorTool/PhysObjectExtractor/JEC/Fall15_25nsV2_MC_SF_AK4PFchs.txt'), # AK8 == AK4
#				)

#----- Propagate the jet energy corrections to the MET -----#
from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
from PhysicsTools.PatUtils.patPFMETCorrections_cff import patPFMetT1T2Corr

process.uncorrectedMet = cms.EDProducer("RecoMETExtractor",
        correctionLevel = cms.string('raw'),
        metSource = cms.InputTag("slimmedMETs", "", "@skipCurrentProcess")
        )

addMETCollection(process, labelName="uncorrectedPatMet", metSource="uncorrectedMet")
process.uncorrectedPatMet.addGenMET = False

#----- Evaluate the Type-1 correction -----#

process.Type1CorrForNewJEC = patPFMetT1T2Corr.clone(
        isMC = cms.bool(isMC),
        src = cms.InputTag("slimmedJetsNewJEC"),
        )
process.slimmedMETsNewJEC = cms.EDProducer('CorrectedPATMETProducer',
        src = cms.InputTag('uncorrectedPatMet'),
        srcCorrections = cms.VInputTag(cms.InputTag('Type1CorrForNewJEC', 'type1'))
        )

#----- Configure the POET MET analyzer -----#
process.mymets = cms.EDAnalyzer('MetAnalyzer',mets=cms.InputTag("slimmedMETsNewJEC"),rawmets=cms.InputTag("uncorrectedPatMet"))


#----- RUN THE JOB! -----#
process.TFileService = cms.Service("TFileService", fileName=cms.string("wjets_example.root"))

if isData:
	process.p = cms.Path(
                             )
else:
	process.p = cms.Path(process.mymuons+process.looseAK4Jets+process.patJetCorrFactorsReapplyJEC+
                             process.slimmedJetsNewJEC+process.myjets+process.looseAK8Jets+process.patJetCorrFactorsReapplyJECAK8+
                             process.slimmedJetsAK8NewJEC+process.uncorrectedMet+process.uncorrectedPatMet+
                             process.Type1CorrForNewJEC+process.slimmedMETsNewJEC+process.mymets
                             )
