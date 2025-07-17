#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_TauClusterCollection_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_TauClusterCollection_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TauClusterHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TauClustersSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/CopyToHost.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  /**
   * make the names from the top-level `l1sc` namespace visible for unqualified lookup
   * inside the `ALPAKA_ACCELERATOR_NAMESPACE::l1sc` namespace
   */
  using namespace ::l1sc;

  using TauClusterCollection = PortableCollection<TauClusterSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::TauClusterCollection, l1sc::TauClusterHostCollection);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_TauClusterCollection_h