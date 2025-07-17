#ifndef DataFormats_L1ScoutingSoA_interface_JETClusterSoA_h
#define DataFormats_L1ScoutingSoA_interface_JETClusterSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(JETClusterLayout,
                      SOA_COLUMN(float, eta),
                      SOA_COLUMN(float, phi),
                      SOA_COLUMN(float, pt),
                      SOA_COLUMN(float, z0),
                      SOA_COLUMN(float, dxy),
                      SOA_COLUMN(int16_t, pdgid))

  using JETClusterSoA = JETClusterLayout<>;
  using JETClusterSoAView = JETClusterSoA::View;
  using JETClusterSoAConstView = JETClusterSoA::ConstView;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_JETClusterSoA_h