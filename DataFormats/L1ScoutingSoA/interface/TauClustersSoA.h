#ifndef DataFormats_L1ScoutingSoA_interface_TauClusterSoA_h
#define DataFormats_L1ScoutingSoA_interface_TauClusterSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/SoATemplate/interface/SoAView.h"

namespace l1sc {

  GENERATE_SOA_LAYOUT(TauClusterLayout,
                      SOA_COLUMN(float, pt),
                      SOA_COLUMN(float, deltaeta),
                      SOA_COLUMN(float, deltaphi),
                      SOA_COLUMN(float, vz),
                      SOA_COLUMN(int8_t, charge),
                      SOA_COLUMN(bool, ishargedhad), 
                      SOA_COLUMN(bool, isneutralhad),
                      SOA_COLUMN(bool, iselectron),
                      SOA_COLUMN(bool, ismuon),
                      SOA_COLUMN(bool, isphoton))

  using TauClusterSoA = TauClusterLayout<>;
  using TauClusterSoAView = TauClusterSoA::View;
  using TauClusterSoAConstView = TauClusterSoA::ConstView;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_TauClustersSoA_h