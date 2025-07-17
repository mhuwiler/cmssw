// -*- C++ -*-
//
// Package:    L1TriggerScouting/TauTagging
// Class:      L1TScPhase2CLUEJets
//
/**\class L1TScPhase2CLUEJets L1TScPhase2CLUEJets.cc L1TriggerScouting/TauTagging/plugins/L1TScPhase2CLUEJets.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Marc Huwiler
//         Created:  Wed, 16 Jul 2025 15:11:12 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1ScoutingSoA/interface/alpaka/CLUEsteringCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/TauClusterCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/TauTagging/interface/L1TScPhase2Common.h"

//
// class declaration
//

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

class L1TScPhase2CLUEJets : public edm::stream::EDProducer<> {
public:
  explicit L1TScPhase2CLUEJets(const edm::ParameterSet&);
  ~L1TScPhase2CLUEJets() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  const device::EDGetToken<PFCandidateCollection> pf_token_;
  const device::EDGetToken<CLUEsteringCollection> clusters_token_;
  const device::EDPutToken<CLUEJetsCollection> jets_token_;

  const bool debug; 

  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

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
L1TScPhase2CLUEJets::L1TScPhase2CLUEJets(const edm::ParameterSet& iConfig) 
  : EDProducer<>(params),
        pf_token_{consumes(params.getParameter<edm::InputTag>("pf"))},
        cluestes_token_{consumes(params.getParameter<edm::InputTag>("clusters"))},
        jets_token_{produces()},
        debug_(params.getUntrackedParameter<bool>("debug")) 
{

  //register your products
  /* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
  */
  //now do what ever other initialization is needed
}

L1TScPhase2CLUEJets::~L1TScPhase2CLUEJets() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void L1TScPhase2CLUEJets::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  const auto& pf = event.get(pf_token_);
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
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void L1TScPhase2CLUEJets::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void L1TScPhase2CLUEJets::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TScPhase2CLUEJets::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
L1TScPhase2CLUEJets::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TScPhase2CLUEJets::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TScPhase2CLUEJets::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1TScPhase2CLUEJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

//define this as a plug-in
DEFINE_FWK_MODULE(l1sc::L1TScPhase2CLUEJets);
