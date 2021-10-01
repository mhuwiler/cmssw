import FWCore.ParameterSet.Config as cms

from DQM.TrackingMonitorSource.trackToTrackComparisonHists_cfi import trackToTrackComparisonHists
mlTrack2generalTracks = trackToTrackComparisonHists.clone()
mlTrack2generalTracks.monitoredTrack          = cms.InputTag("trackCollectionKFfromML") # your tracks collection
mlTrack2generalTracks.referenceTrack            = cms.InputTag("generalTracks")
mlTrack2generalTracks.monitoredBeamSpot  = cms.InputTag("offlineBeamSpot")
mlTrack2generalTracks.referenceBeamSpot    = cms.InputTag("offlineBeamSpot")
mlTrack2generalTracks.topDirName                = cms.string("Tracking/GNNTracksWRTgeneralTracks/")
mlTrack2generalTracks.monitoredPrimaryVertices = cms.InputTag("offlinePrimaryVertices")
mlTrack2generalTracks.referencePrimaryVertices   = cms.InputTag("offlinePrimaryVertices")


mlTrack2generalTracks.histoPSet = trackToTrackComparisonHists.histoPSet .clone(
    Eta_rangeMin = cms.double(-4.0),
    Eta_rangeMax = cms.double(4.0),
    Eta_nbin = cms.uint32(160),
    Dxy_rangeMin = cms.double(-5),
    Dxy_rangeMax = cms.double(5),
    Dxy_nbin = cms.uint32(1500)
)

mlTrackToTrackSelector = mlTrack2generalTracks.clone()
mlTrackToTrackSelector.referenceTrack= cms.InputTag("trackSelector")
mlTrackToTrackSelector.topDirName = cms.string("Tracking/GNNTracksWRTTrackSelector/")


mlToTrackMonitoring = cms.Sequence(mlTrack2generalTracks+mlTrackToTrackSelector)

