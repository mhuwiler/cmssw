import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer

pixelverticesvalid = DQMEDAnalyzer("SiPixelValidateVerticesGPU",
    trackCollectionsrc = cms.InputTag("pixelTrackSoA"),
    beamSpotsrc = cms.InputTag("offlineBeamSpot"),
    mightGet = cms.optional.untracked.vstring,
    src = cms.InputTag("pixelVertexCUDA")
)
