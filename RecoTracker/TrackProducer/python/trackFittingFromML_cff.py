
import FWCore.ParameterSet.Config as cms

#fragment = cms.ProcessFragment( "TrackFittingFromML" )
# import propagator cff 


from TrackingTools.MaterialEffects.MaterialPropagator_cfi import *

from TrackingTools.MaterialEffects.OppositeMaterialPropagator_cfi import *

from RecoTracker.TkTrackingRegions.globalTrackingRegionFromBeamSpot_cfi import globalTrackingRegionFromBeamSpot as _globalTrackingRegionFromBeamSpot
trackingRegion = _globalTrackingRegionFromBeamSpot.clone(
    RegionPSet = dict(
        nSigmaZ = cms.double( 4.0 ), #TODO: tune 
        ptMin = cms.double( 0.1 ), # match threshold from ML 
        originRadius = cms.double( 0.02 ) #TODO: tune
    )
)

trackCollectionKFfromML = cms.EDProducer ("TrackFitterFromML", 
    propagator = cms.string("PropagatorWithMaterial"), 
    oppositePropagator = cms.string("PropagatorWithMaterialOpposite"), 
    ttRecHitBuilder = cms.string("PixelTTRHBuilderWithoutAngle"), 
    #trackerGeometry = cms.InputTag(""),   
    #magneticField = cms.InputTag(""), 
    beamSpot = cms.InputTag("offlineBeamSpot"), 
    trackingRegion = cms.InputTag("trackingRegion"), 
    doTest = cms.bool(False), 
    tracks = cms.InputTag("generalTracks")
)


trackFittingKFFromRecHit = cms.Sequence(trackingRegion+trackCollectionKFfromML)

