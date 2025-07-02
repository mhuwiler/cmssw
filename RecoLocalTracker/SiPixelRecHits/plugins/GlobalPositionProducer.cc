// -*- C++ -*-
//
// Package:    RecoLocalTracker/SiPixelRecHits
// Class:      GlobalPositionProducer
//
/**\class GlobalPositionProducer GlobalPositionProducer.cc RecoLocalTracker/SiPixelRecHits/plugins/GlobalPositionProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marc Huwiler
//         Created:  Mon, 23 Jun 2025 14:15:09 GMT
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

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/OwnVector.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"


typedef edm::OwnVector<TrackingRecHit,edm::ClonePolicy<TrackingRecHit> > TrackingRecHitCollection; 


//
// class declaration
//

class GlobalPositionProducer : public edm::stream::EDProducer<> {
public:
  explicit GlobalPositionProducer(const edm::ParameterSet&);
  ~GlobalPositionProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<TrackingRecHitCollection> token_;

  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;

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
GlobalPositionProducer::GlobalPositionProducer(const edm::ParameterSet& iConfig) {
  //register your products
  produces<edm::ValueMap<GlobalPoint>>();

  geomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();

  /* Examples
  

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
  */
  //now do what ever other initialization is needed
}

GlobalPositionProducer::~GlobalPositionProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void GlobalPositionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
  using namespace edm;
  /* This is an event example
  //Read 'ExampleData' from the Event
  ExampleData const& in = iEvent.get(inToken_);

  //Use the ExampleData to create an ExampleData2 which 
  // is put into the Event
  iEvent.put(std::make_unique<ExampleData2>(in));
  */

  /* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  SetupData& setup = iSetup.getData(setupToken_);
  */

  edm::Handle<TrackingRecHitCollection> hits;
  iEvent.getByToken(token_, hits);

  //edm::ESHandle<TrackerGeometry> geom;
  //iSetup.get<TrackerGeometry>().get(geom);

  const TrackerGeometry& geometry = iSetup.getData(geomToken_);

  std::unique_ptr<ValueMap<GlobalPoint> > refs; 
  //refs.reserve(hits.size()); 


  std::vector<GlobalPoint> values; 
  //values.reserve(hits.size()); 

  for (size_t i = 0; i < hits->size(); ++i) 
  {
  	const auto& hit = (*hits)[i];
  	auto det = geometry.idToDet(hit.geographicalId());
  	GlobalPoint gp = det->surface().toGlobal(hit.localPosition());
  	values.push_back(gp);
  	//refs->push_back(edm::Ref<TrackingRecHitCollection>(hits, i));
  }

  edm::ValueMap<GlobalPoint>::Filler filler(*refs);

  filler.insert(hits, values.begin(), values.end()); 

  filler.fill(); 


  iEvent.put(std::move(refs)); 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void GlobalPositionProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void GlobalPositionProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
GlobalPositionProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
GlobalPositionProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GlobalPositionProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GlobalPositionProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GlobalPositionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GlobalPositionProducer);
