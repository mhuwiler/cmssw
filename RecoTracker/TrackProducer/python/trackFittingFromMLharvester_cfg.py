# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step4 --conditions auto:phase2_realistic_T15 -s HARVESTING:@trackingOnlyValidation+@trackingOnlyDQM --scenario pp --filetype DQM --geometry Extended2026D49 --era Phase2C9 --mc -n 100 --filein file:step3_inDQM.root --fileout file:step4.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9

process = cms.Process('HARVESTING',Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.DQMSaverAtRunEnd_cff')
process.load('Configuration.StandardSequences.Harvesting_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring('file:step3_RAW2DIGI_RECO_VALIDATION_DQM_inDQM.root')
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step4 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
process.mlTrack2generalTracksEfficiencies = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring(
        "Tracking/GNNTracksWRTgeneralTracks",
        "Tracking/GNNTracksWRTTrackSelector",
    ),
    verbose        = cms.untracked.uint32(0),
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "eff_pt              'Relative Efficiency vs Pt;#P_T;relative efficiency'               ref_matched_pt          ref_pt          eff",
        "eff_eta             'Relative Efficiency vs Eta;#eta;relative efficiency'              ref_matched_eta         ref_eta         eff",
        "eff_phi             'Relative Efficiency vs Phi;#phi;relative efficiency'              ref_matched_phi         ref_phi         eff",
        "eff_dxy             'Relative Efficiency vs dxy;d_{xy};relative efficiency'            ref_matched_dxy         ref_dxy         eff",
        "eff_dz              'Relative Efficiency vs dz;d_{z};relative efficiency'              ref_matched_dz          ref_dz          eff",
        "eff_dxyWRTpv        'Relative Efficiency vs dxyWRTpv;d_{xy};relative efficiency'       ref_matched_dxyWRTpv    ref_dxyWRTpv    eff",
        "eff_dzWRTpv         'Relative Efficiency vs dzWRTpv;d_{z};relative efficiency'         ref_matched_dzWRTpv     ref_dzWRTpv     eff",
        "eff_charge          'Relative Efficiency vs charge;charge;relative efficiency'         ref_matched_charge      ref_charge      eff",
        "eff_hits            'Relative Efficiency vs hits;number of hits;relative efficiency'   ref_matched_hits        ref_hits        eff",

        "fakeRate_pt         'Relative Fake Rate vs Pt;#P_T;relative fake rate'                 mon_unMatched_pt        mon_pt          eff",
        "fakeRate_eta        'Relative Fake Rate vs Eta;#eta;relative fake rate'                mon_unMatched_eta       mon_eta         eff",
        "fakeRate_phi        'Relative Fake Rate vs Phi;#phi;relative fake rate'                mon_unMatched_phi       mon_phi         eff",
        "fakeRate_dxy        'Relative Fake Rate vs dxy;d_{xy};relative fake rate'              mon_unMatched_dxy       mon_dxy         eff",
        "fakeRate_dz         'Relative Fake Rate vs dz;d_{z};relative fake rate'                mon_unMatched_dz        mon_dz          eff",
        "fakeRate_dxyWRTpv   'Relative Fake Rate vs dxyWRTpv;d_{xy};relative fake rate'         mon_unMatched_dxyWRTpv  mon_dxyWRTpv    eff",
        "fakeRate_dzWRTpv    'Relative Fake Rate vs dzWRTpv;d_{z};relative fake rate'           mon_unMatched_dzWRTpv   mon_dzWRTpv     eff",
        "fakeRate_charge     'Relative Fake Rate vs charge;charge;relative fake rate'           mon_unMatched_charge    mon_charge      eff",
        "fakeRate_hits       'Relative Fake Rate vs hits;number of hits;relative fake rate'     mon_unMatched_hits      mon_hits        eff",
    ),
)

# Output definition

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

# Path and EndPath definitions
#process.genHarvesting = cms.Path(process.postValidation_gen)
#process.validationprodHarvesting = cms.Path(process.hltpostvalidation_prod+process.postValidation_gen)
#process.validationHarvestingHI = cms.Path(process.postValidationHI)
#process.dqmHarvestingExtraHLT = cms.Path(process.DQMOffline_SecondStep_ExtraHLT+process.DQMOffline_Certification)
#process.alcaHarvesting = cms.Path()
#process.validationHarvestingFS = cms.Path(process.recoMuonPostProcessors+process.postValidationTracking+process.MuIsoValPostProcessor+process.calotowersPostProcessor+process.hcalSimHitsPostProcessor+process.hcaldigisPostProcessor+process.hcalrechitsPostProcessor+process.electronPostValidationSequence+process.photonPostProcessor+process.pfJetClient+process.pfMETClient+process.pfJetResClient+process.pfElectronClient+process.rpcRecHitPostValidation_step+process.makeBetterPlots+process.bTagCollectorSequenceMCbcl+process.METPostProcessor+process.L1GenPostProcessor+process.bdHadronTrackPostProcessor+process.MuonGEMHitsPostProcessors+process.MuonGEMDigisPostProcessors+process.MuonGEMRecHitsPostProcessors+process.hgcalPostProcessor+process.MuonME0DigisPostProcessors+process.MuonME0SegPostProcessors+process.trackerphase2ValidationHarvesting+process.postValidation_gen)
#process.validationpreprodHarvesting = cms.Path(process.postValidation_preprod+process.hltpostvalidation_preprod+process.postValidation_gen)
#process.validationHarvestingNoHLT = cms.Path(process.postValidation+process.postValidation_gen)
#process.validationHarvesting = cms.Path(process.postValidation+process.hltpostvalidation+process.postValidation_gen)
#process.validationpreprodHarvestingNoHLT = cms.Path(process.postValidation_preprod+process.postValidation_gen)
#process.dqmHarvestingPOGMC = cms.Path(process.DQMOffline_SecondStep_PrePOGMC)
#process.dqmHarvestingFakeHLT = cms.Path(process.DQMOffline_SecondStep_FakeHLT+process.DQMOffline_Certification)
#process.validationHarvestingMiniAOD = cms.Path(process.JetPostProcessor+process.METPostProcessorHarvesting+process.postValidationMiniAOD)
#process.dqmHarvesting = cms.Path(process.DQMOffline_SecondStep+process.DQMOffline_Certification)
#process.postValidation_trackingOnly_step = cms.Path(process.postValidation_trackingOnly)
#process.DQMHarvestTracking_step = cms.Path(process.DQMHarvestTracking)
process.trackToTrackCompaarison_step = cms.Path(process.mlTrack2generalTracksEfficiencies)
process.dqmsave_step = cms.Path(process.DQMSaver)

# Schedule definition
process.schedule = cms.Schedule(process.trackToTrackCompaarison_step,process.dqmsave_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)



# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
