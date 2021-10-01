
import FWCore.ParameterSet.Config as cms

#fragment = cms.ProcessFragment( "TrackFittingFromML" )
# import propagator cff 


from TrackingTools.MaterialEffects.MaterialPropagator_cfi import *

from TrackingTools.MaterialEffects.OppositeMaterialPropagator_cfi import *

from RecoTracker.TkTrackingRegions.globalTrackingRegionFromBeamSpot_cfi import globalTrackingRegionFromBeamSpot as _globalTrackingRegionFromBeamSpot
trackingRegion = _globalTrackingRegionFromBeamSpot.clone(
    RegionPSet = dict(
        nSigmaZ = cms.double( 4.0 ), #TODO: tune 
        ptMin = cms.double( 1.0 ), # match threshold from ML 
        originRadius = cms.double( 1.0 ) #TODO: tune
    )
)

trackSelector = cms.EDFilter('TrackSelector',
    src = cms.InputTag('generalTracks'),
    cut = cms.string("abs(eta)<=2.4&pt>=1.0")
)

trackCollectionKFfromML = cms.EDProducer ("TrackFitterFromML", 
    propagator = cms.string("PropagatorWithMaterial"), 
    oppositePropagator = cms.string("PropagatorWithMaterialOpposite"), 
    ttRecHitBuilder = cms.string("PixelTTRHBuilderWithoutAngle"), 
    #trackerGeometry = cms.InputTag(""),   
    #magneticField = cms.InputTag(""), 
    beamSpot = cms.InputTag("offlineBeamSpot"), 
    trackingRegion = cms.InputTag("trackingRegion"), 
    doTest = cms.bool(True), 
    tracks = cms.InputTag("generalTracks") 
)

trackCollectionFromMLSelector = trackCollectionKFfromML.clone(
    tracks = cms.InputTag("trackSelector")
)


trackFittingKFFromRecHit = cms.Sequence(cms.ignore(trackSelector)+trackingRegion+trackCollectionKFfromML+trackCollectionFromMLSelector)

