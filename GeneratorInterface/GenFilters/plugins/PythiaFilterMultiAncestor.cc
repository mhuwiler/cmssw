/*
This filter allows to select events where a given particle originates *either directly 
(daughter) or indirectly ((grand)^n-daughter)*  from a list of possible ancestors.

It also allows to filter on the immediate daughters of the said particle.

Kinematic selections can also be applied on the particle and its daughters.

The example below shows how to select a jpsi->mumu coming from *any* b hadron (notice
that MotherIDs is just 5, i.e. b-quark) however long the decay chain is.
This includes, for example, B->Jpsi + X as well as B->Psi(2S)(->JPsi) X

process.jpsi_from_bhadron_filter = cms.EDFilter("PythiaFilterMultiAncestor",
    DaughterIDs = cms.untracked.vint32(-13, 13),
    DaughterMaxEtas = cms.untracked.vdouble(3., 3.),
    DaughterMaxPts = cms.untracked.vdouble(100000.0, 100000.0),
    DaughterMinEtas = cms.untracked.vdouble(-2.6, -2.6),
    DaughterMinPts = cms.untracked.vdouble(2.5, 2.5),
    MaxEta = cms.untracked.double(3.0),
    MinEta = cms.untracked.double(-3.0),
    MinPt = cms.untracked.double(6.0),
    MotherIDs = cms.untracked.vint32(5),
    ParticleID = cms.untracked.int32(443)
)
*/

// system include files
#include <memory>
#include <iostream>

// user include files
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// #include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/global/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "GeneratorInterface/GenFilters/plugins/MCFilterZboostHelper.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

using namespace edm;
using namespace std;

namespace edm {
  class HepMCProduct;
}

class PythiaFilterMultiAncestor : public edm::global::EDFilter<> {
public:
  explicit PythiaFilterMultiAncestor(const edm::ParameterSet&);

  bool filter(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

private:
  bool isAncestor(HepMC::GenParticle* particle, int IDtoMatch, bool chargeConj = false) const;

  bool hasDaughters(const std::vector<int>& daughters, const HepMC::GenParticle* particle, bool chargeConj = false) const; 

  const edm::EDGetTokenT<edm::HepMCProduct> token_;
  const int particleID;
  const double minpcut;
  const double maxpcut;
  const double minptcut;
  const double maxptcut;
  const double minetacut;
  const double maxetacut;
  const double minrapcut;
  const double maxrapcut;
  const double minphicut;
  const double maxphicut;

  const int status;
  const std::vector<int> motherIDs;
  const std::vector<int> daughterIDs;
  const std::vector<double> daughterMinPts;
  const std::vector<double> daughterMaxPts;
  const std::vector<double> daughterMinEtas;
  const std::vector<double> daughterMaxEtas;

  const int processID;

  const double betaBoost;
  const bool considerCC = true; // TODO: change to an option 

};

PythiaFilterMultiAncestor::PythiaFilterMultiAncestor(const edm::ParameterSet& iConfig)
    : token_(consumes<edm::HepMCProduct>(
          edm::InputTag(iConfig.getUntrackedParameter("moduleLabel", std::string("generator")), "unsmeared"))),
      particleID(iConfig.getUntrackedParameter("ParticleID", 0)),
      minpcut(iConfig.getUntrackedParameter("MinP", 0.)),
      maxpcut(iConfig.getUntrackedParameter("MaxP", 10000.)),
      minptcut(iConfig.getUntrackedParameter("MinPt", 0.)),
      maxptcut(iConfig.getUntrackedParameter("MaxPt", 10000.)),
      minetacut(iConfig.getUntrackedParameter("MinEta", -10.)),
      maxetacut(iConfig.getUntrackedParameter("MaxEta", 10.)),
      minrapcut(iConfig.getUntrackedParameter("MinRapidity", -20.)),
      maxrapcut(iConfig.getUntrackedParameter("MaxRapidity", 20.)),
      minphicut(iConfig.getUntrackedParameter("MinPhi", -3.5)),
      maxphicut(iConfig.getUntrackedParameter("MaxPhi", 3.5)),
      status(iConfig.getUntrackedParameter("Status", 0)),
      motherIDs(iConfig.getUntrackedParameter("MotherIDs", std::vector<int>{0})),
      daughterIDs(iConfig.getUntrackedParameter("DaughterIDs", std::vector<int>{0})),
      daughterMinPts(iConfig.getUntrackedParameter("DaughterMinPts", std::vector<double>{0.})),
      daughterMaxPts(iConfig.getUntrackedParameter("DaughterMaxPts", std::vector<double>{10000.})),
      daughterMinEtas(iConfig.getUntrackedParameter("DaughterMinEtas", std::vector<double>{-10.})),
      daughterMaxEtas(iConfig.getUntrackedParameter("DaughterMaxEtas", std::vector<double>{10.})),
      processID(iConfig.getUntrackedParameter("ProcessID", 0)),
      betaBoost(iConfig.getUntrackedParameter("BetaBoost", 0.)) {
  //now do what ever initialization is needed
}

// ------------ access the full genealogy ---------
bool PythiaFilterMultiAncestor::isAncestor(HepMC::GenParticle* particle, int IDtoMatch, bool chargeConj) const 
{
  bool result = false; 

  for (HepMC::GenVertex::particle_iterator ancestor = particle->production_vertex()->particles_begin(HepMC::ancestors);
       ancestor != particle->production_vertex()->particles_end(HepMC::ancestors); // If multiple mothers with same ID are required, will return true possibly on the same particle 
       ++ancestor) 
  {
    // std::cout << __LINE__ << "]\t particle's PDG ID " << particle->pdg_id()
    //                       << " \t particle's ancestor's PDG ID " << (*ancestor)->pdg_id()
    //                       << " \t ID to match " << IDtoMatch << std::endl;

    if (((*ancestor)->pdg_id() == IDtoMatch) && !chargeConj) 
    {
      //  std::cout << __LINE__ << "]\t found!" << std::endl;
      result = true; 
    }
    else if (((*ancestor)->pdg_id() == -IDtoMatch) && chargeConj)
    {
      result = true;
    }
  }

  // std::cout << __LINE__ << "]\t nope, no luck" << std::endl;
  return result;
}

bool PythiaFilterMultiAncestor::hasDaughters(const std::vector<int>& daughters, const HepMC::GenParticle* particle, const bool chargeConj) const 
{
    bool result = false; 

    int coeff = 1; 

    if (chargeConj) coeff = -1; 

    uint good_dau = 0;
    int idx = -1; 
    bool matchingTable[particle->end_vertex()->particles_out_size()][daughters.size()]; 
    for (HepMC::GenVertex::particle_iterator dau = particle->end_vertex()->particles_begin(HepMC::children);
         dau != particle->end_vertex()->particles_end(HepMC::children); ++dau) 
    {
      idx++; 
      for (unsigned int i = 0; i < daughters.size(); ++i) 
      {
        // if a daughter has its pdgID among the desired ones, apply kin cuts on it
        // if it survives, add a notch to the counter
        if ((*dau)->pdg_id() == coeff*daughterIDs[i]) 
        {
          //std::cout << "Particle matching " << std::endl; 
          if (((*dau)->momentum().perp() > daughterMinPts[i]) && ((*dau)->momentum().perp() < daughterMaxPts[i]) && 
              ((*dau)->momentum().eta() > daughterMinEtas[i]) && ((*dau)->momentum().eta() < daughterMaxEtas[i])) 
          {
              ++good_dau; 
          }
        }
      }
    }
    return (good_dau >= daughterIDs.size()); 
}

// ------------ method called to produce the data  ------------
bool PythiaFilterMultiAncestor::filter(edm::StreamID, edm::Event& iEvent, const edm::EventSetup&) const
{
  using namespace edm;
  bool accepted = false;
  Handle<HepMCProduct> evt;
  iEvent.getByToken(token_, evt);
  bool isCC = false; 


  const HepMC::GenEvent* myGenEvent = evt->GetEvent();

  if (processID == 0 || processID == myGenEvent->signal_process_id()) 
  {
    for (HepMC::GenEvent::particle_const_iterator p = myGenEvent->particles_begin(); p != myGenEvent->particles_end(); ++p) 
    {
      HepMC::FourVector momentum = MCFilterZboostHelper::zboost((*p)->momentum(), betaBoost);
      double rapidity = 0.5 * log((momentum.e() + momentum.pz()) / (momentum.e() - momentum.pz()));

      int pid = (*p)->pdg_id(); 
      if(pid == particleID) 
      {
          isCC = false; // Reset it in each loop 
      }
      else if (pid == -particleID) 
      {
          isCC = true; 
          if (!considerCC) continue; 
      }
      else 
      {
        continue; 
      }

      if (momentum.rho() > minpcut && momentum.rho() < maxpcut &&
          (*p)->momentum().perp() > minptcut && (*p)->momentum().perp() < maxptcut && momentum.eta() > minetacut &&
          momentum.eta() < maxetacut && rapidity > minrapcut && rapidity < maxrapcut && (*p)->momentum().phi() > minphicut &&
          (*p)->momentum().phi() < maxphicut) 
      {
        // Check the status of the particle 
        bool statusPass = ((status == 0) || ((*p)->status() == status)); 

        // find the mother
        bool momFound = false; 
        for (std::vector<int>::const_iterator motherID = motherIDs.begin(); motherID != motherIDs.end(); ++motherID) 
        {
          if ((*motherID == 0) || isAncestor(*p, *motherID, isCC)) momFound = true; // If one of moms is found, set to true 

        }

        // find the daughters
        if (statusPass && momFound && (!daughterIDs.empty())) 
        {
          // if you got this far it means that the mother was found
          // now let's check the daughters
          // use a counter, if there's enough daughters that match the pdg and kinematic
          // criteria accept the event
          if (hasDaughters(daughterIDs, *p, isCC)) accepted = true; 
        }
        std::cout << "Has daughters: " << accepted << ", " << hasDaughters(daughterIDs, *p, isCC) << std::endl; 
      }
      // only need to satisfy the conditions _once_
      if (accepted)
        break;
    }

  } else 
  {
    accepted = true;
  }

  if (accepted) 
  {
    return true;
  } 
  else 
  {
    return false;
  }
}

DEFINE_FWK_MODULE(PythiaFilterMultiAncestor);
