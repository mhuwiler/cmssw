#ifndef L1TriggerScouting_TauTagging_plugins_alpaka_JETConcatenation_h
#define L1TriggerScouting_TauTagging_plugins_alpaka_JETConcatenation_h

// libs
#include <alpaka/alpaka.hpp>
// typedefs
#include "DataFormats/L1ScoutingSoA/interface/alpaka/CLUEsteringCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateCollection.h"
// heterogeneous
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {


extern void  Concatenate(Queue& queue, const PFCandidateCollection& pf, const CLUEsteringCollection& clusters, const uint32_t clusters_num); 


}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif  // L1TriggerScouting_TauTagging_plugins_alpaka_JETConcatenation_h