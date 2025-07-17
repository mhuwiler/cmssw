#ifndef DataFormats_L1ScoutingSoA_interface_TauClusterHostCollection_h
#define DataFormats_L1ScoutingSoA_interface_TauClusterHostCollection_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/TauClustersSoA.h"

namespace l1sc {

  using TauClusterHostCollection = PortableHostCollection<TauClusterSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_TauClusterHostCollection_h