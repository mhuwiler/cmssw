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
//#include <memory>

// user include files
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Utilities/interface/StreamID.h"

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
//#include "L1TriggerScouting/TauTagging/plugins/alpaka/JETKernel.cc"
#include "L1TriggerScouting/TauTagging/interface/alpaka/JETConcatenation.h"

//
// class declaration
//
constexpr int batchsize = 16; 


namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  class L1TScPhase2CLUEJets : public stream::EDProducer<> {
  public:
    explicit L1TScPhase2CLUEJets(const edm::ParameterSet&);
    ~L1TScPhase2CLUEJets() override;

    void produce(device::Event &event, const device::EventSetup &event_setup) override;

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    

    const device::EDGetToken<PFCandidateCollection> pf_token_;
    const device::EDGetToken<CLUEsteringCollection> clusters_token_;
    const device::EDPutToken<TauClusterCollection> jets_token_;

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
    : EDProducer<>(iConfig),
      pf_token_{consumes(iConfig.getParameter<edm::InputTag>("pf"))},
      clusters_token_{consumes(iConfig.getParameter<edm::InputTag>("clusters"))},
      jets_token_{produces()},
      debug(iConfig.getUntrackedParameter<bool>("debug")) 
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



// ------------ method called to produce the data  ------------
  void L1TScPhase2CLUEJets::produce(device::Event &event, const device::EventSetup &evtsetup) {
    //using namespace edm;

    const auto& pf = event.get(pf_token_);

    const auto& clusters = event.get(clusters_token_);

    const int32_t n = clusters.view().metadata().size(); 

    uint32_t numJets = 0; 

    for (int32_t i=0; i<n; i++) 
    {

      // Accessing the column "cluster" by its name as a finctional
      auto a = clusters.view().cluster()[i]; 

      std::cout << a << std::endl; 

      if (a > numJets) numJets = a; 
    }

    numJets +=1; // numbering of clusters starts at 0


    //JETConcatenationKernel concatenate; 


    /*alpaka::exec<Tag>(
        event.queue(),
        workDiv,
        concatenate, 
        pf, 
        clusters, 
        nJets);
    alpaka::wait(event.queue()); */
    //JETConcatenation concat; 

    Concatenate<Acc1D>(event.queue(), pf, clusters, numJets);



    int nConst = batchsize*numJets; 

    // creating the output collection
    auto jets = TauClusterCollection(nConst, event.queue());
    jets.zeroInitialise(event.queue());
  
  }

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
  void L1TScPhase2CLUEJets::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("pf");
    desc.add<edm::InputTag>("clusters");
    //desc.add<edm::InputTag>("NNTaus");
    desc.addUntracked<bool>("debug", false);
    descriptions.addWithDefaultLabel(desc);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

//define this as a plug-in
DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2CLUEJets);
