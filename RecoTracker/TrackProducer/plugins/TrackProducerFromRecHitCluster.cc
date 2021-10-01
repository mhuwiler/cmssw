// -*- C++ -*-
//
// Package:    TrackingML/TrackFitterFromML
// Class:      TrackFitterFromML
//
/**\class TrackFitterFromML TrackFitterFromML.cc TrackingML/TrackFitterFromML/plugins/TrackFitterFromML.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marc Huwiler
//         Created:  Mon, 27 Sep 2021 12:15:10 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "RecoPixelVertexing/PixelTrackFitting/interface/KFBasedPixelFitter.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "RecoPixelVertexing/PixelTrackFitting/plugins/storeTracks.h"
#include "RecoPixelVertexing/PixelTrackFitting/interface/TracksWithHits.h"
//#include "RecoTracker/TkTrackingRegions/interface/GlobalTrackingRegion.h"



//
// class declaration
//

class TrackFitterFromML : public edm::stream::EDProducer<> {
public:
    explicit TrackFitterFromML(const edm::ParameterSet&);
    ~TrackFitterFromML();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    void beginStream(edm::StreamID) override;
    void produce(edm::Event&, const edm::EventSetup&) override;
    void endStream() override;

    std::vector<std::vector<const TrackingRecHit *> > getRechitsFromTracks(const edm::Event&); 

    edm::Handle<reco::BeamSpot> hBeamSpot; 

    //edm::EDPutTokenT<reco::TrackCollection> trackPutToken; 

    const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken;
    const edm::EDGetTokenT<edm::OwnVector<TrackingRegion> > trackingRegionToken; 
    const edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> ttTrackBuilderToken;
    const edm::ESGetToken<Propagator, TrackingComponentsRecord> trackPropagator;
    const edm::ESGetToken<Propagator, TrackingComponentsRecord> trackPropagatorOppositeToken;
    const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken;
    const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> fieldToken; 
    const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoToken;
    bool doTest = false; 
    edm::EDGetTokenT<reco::TrackCollection> tracksToken; 

    PixelFitterBase *fitter = nullptr; 

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
TrackFitterFromML::TrackFitterFromML(const edm::ParameterSet& iConfig)
    : beamSpotToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
      trackingRegionToken(consumes<edm::OwnVector<TrackingRegion> >(iConfig.getParameter<edm::InputTag>("trackingRegion"))), 
      ttTrackBuilderToken(esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("ttRecHitBuilder")))),
      trackPropagator(esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("propagator")))),
      trackPropagatorOppositeToken(esConsumes(edm::ESInputTag("", iConfig.getParameter<std::string>("oppositePropagator")))),
      trackerGeometryToken(esConsumes()),
      fieldToken(esConsumes()), 
      tTopoToken(esConsumes<TrackerTopology, TrackerTopologyRcd>()), 
      doTest(iConfig.getParameter<bool>("doTest")),
      tracksToken(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))) 
{
    produces<reco::TrackCollection>();
    produces<TrackingRecHitCollection>(); 
    produces<reco::TrackExtraCollection>(); 

    
    if (doTest) 
    {
        
    }
}

TrackFitterFromML::~TrackFitterFromML() 
{
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void TrackFitterFromML::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
    using namespace edm;

    std::unique_ptr<reco::TrackCollection> tracks = std::make_unique<reco::TrackCollection>(); 

    iEvent.getByToken(beamSpotToken, hBeamSpot); 


    std::vector<std::vector<const TrackingRecHit *> > mlProtoTracks; // TODO: get this from the DNN somehow 

    if (doTest) 
    {
        // Get the RecHits from tracks 
        mlProtoTracks = getRechitsFromTracks(iEvent); 
    }

    tracks->reserve(mlProtoTracks.size()); // TODO: do we really want to do this? TrackingRegion might reduce the phase space 

    fitter = new KFBasedPixelFitter(&iSetup.getData(trackPropagator),
                                    &iSetup.getData(trackPropagatorOppositeToken),
                                    &iSetup.getData(ttTrackBuilderToken),
                                    &iSetup.getData(trackerGeometryToken),
                                    &iSetup.getData(fieldToken),
                                    hBeamSpot.product());  

    const TrackingRegion *trackingRegion = &(*iEvent.get(trackingRegionToken).begin()); // TODO: loop over all tracking region 

    //Retrieve tracker topology from geometry
    //edm::ESHandle<TrackerTopology> tTopoHandle = iSetup.getHandle(tTopoToken);
    auto topology= iSetup.getData(tTopoToken); 

    pixeltrackfitting::TracksWithRecHits trackHitMap; 

    for (auto hits : mlProtoTracks) 
    {
        assert(fitter); 
        /*std::cout << "Before fitting step" << std::endl; 
        // For each collection of RecHits provided by the DNN, aka protoTrack, we perform a Kalman Fit
        std::unique_ptr<reco::Track> tempTrack = std::move(fitter->run(recHits, *trackingRegion, iSetup)); //TODO: in version 12, the setup is no longer required 
        std::cout << "After fitting step" << std::endl; 
        for (auto hit : recHits) 
        {
            tempTrack->appendHitPattern(&(*hit), topology); 
        }
        //tracks->emplace_back(tempTrack.release()); 
        trackHitMap.push_back(std::make_pair<reco::Track, std::vector)
        std::cout << "After pushing into collection" << std::endl; */

        // fitting
        std::unique_ptr<reco::Track> track = fitter->run(hits, *trackingRegion, iSetup); 
        if (!track) continue; 

        // add tracks
        trackHitMap.emplace_back(track.release(), hits);
    }

    std::cout << "Number of (protoTracks, fittedTracks): " << mlProtoTracks.size() << ", " << trackHitMap.size() << std::endl; 
    
    storeTracks(iEvent, trackHitMap, topology); 

    //iEvent.put(trackPutToken, std::move(tracks)); 

    delete fitter; 
}

std::vector<std::vector<const TrackingRecHit *> > TrackFitterFromML::getRechitsFromTracks(const edm::Event& iEvent) 
{
    edm::Handle<reco::TrackCollection> hTracks; 
    iEvent.getByToken(tracksToken, hTracks); 

    std::vector<std::vector<const TrackingRecHit *> > recHitCollection; 

    for (auto track : *hTracks.product()) 
    {
        // Access the RecHits
        std::vector<const TrackingRecHit*> recHitTrack; 
        for (auto recHit : track.recHits()) 
        {
            // Fill the new datafromat with the recHits 
            if (!recHit->isValid()) continue; 
            recHitTrack.push_back(recHit); 
        }
        if (recHitTrack.size()<3) continue; 
        recHitCollection.push_back(recHitTrack); 
    }

    return recHitCollection; 

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void TrackFitterFromML::beginStream(edm::StreamID) 
{
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void TrackFitterFromML::endStream() 
{
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
TrackFitterFromML::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
TrackFitterFromML::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TrackFitterFromML::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TrackFitterFromML::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackFitterFromML::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.add<std::string>("propagator", "PropagatorWithMaterial");
    desc.add<std::string>("oppositePropagator", "PropagatorWithMaterialOpposite"); 
    desc.add<std::string>("ttRecHitBuilder", "PixelTTRHBuilderWithoutAngle"); 
    //desc.add<edm::InputTag>("trackerGeometry", edm::InputTag("pixelVertexSoA"));  
    //desc.add<edm::InputTag>("magneticField", edm::InputTag("GlobalTrackingRegionFromBeamSpotEDProducer")); 
    desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));  
    desc.add<edm::InputTag>("trackingRegion", edm::InputTag("GlobalTrackingRegionFromBeamSpotEDProducer")); 
    desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks")); 
    desc.add<bool>("doTest", false); 
    auto label = "trackProducerFromRecHitCluster"; // Is this the name to access the collection? 
    descriptions.add(label, desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackFitterFromML);
