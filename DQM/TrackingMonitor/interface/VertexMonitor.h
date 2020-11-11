#ifndef VertexMonitor_H
#define VertexMonitor_H
// -*- C++ -*-
//
//
/**\class VertexMonitor VertexMonitor.cc
Monitoring source for general quantities related to vertex
*/

// system include files
#include <memory>

// user include files
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

class GetLumi;

class VertexMonitor {
public:
  typedef dqm::legacy::DQMStore DQMStore;
  typedef dqm::legacy::MonitorElement MonitorElement;

  VertexMonitor(const edm::ParameterSet&, const edm::InputTag&, const edm::InputTag&, std::string pvLabel);
  VertexMonitor(const edm::ParameterSet&,
                const edm::InputTag&,
                const edm::InputTag&,
                std::string pvLabel,
                edm::ConsumesCollector& iC);

  virtual ~VertexMonitor();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  virtual void initHisto(DQMStore::IBooker& ibooker);
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &); // From PrimaryVertexMonitor

  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
private:
  void pvTracksPlots(const reco::Vertex &v); // From PrimaryVertexMonitor
  void vertexPlots(const reco::Vertex &v, const reco::BeamSpot &beamSpot, int i); // From PrimaryVertexMonitor

  edm::ParameterSet conf_;

  edm::InputTag primaryVertexInputTag_;
  edm::InputTag selectedPrimaryVertexInputTag_;
  std::string label_;

  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDGetTokenT<reco::VertexCollection> selpvToken_;

  edm::EDGetTokenT<reco::VertexCollection> vertexToken_; // Duplicate from PrimaryVertexMonitor
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_; // From PrimaryVertexMonitor
  using VertexScore = edm::ValueMap<float>; // From PrimaryVertexMonitor
  edm::EDGetTokenT<VertexScore> scoreToken_; // From PrimaryVertexMonitor

  edm::InputTag vertexInputTag_, beamSpotInputTag_; // Duplicate from PrimaryVertexMonitor

  std::string dqmLabel; // Duplicate from PrimaryVertexMonitor

  std::string TopFolderName_; // From PrimaryVertexMonitor
  std::string AlignmentLabel_; // From PrimaryVertexMonitor
  int ndof_; // From PrimaryVertexMonitor
  bool errorPrinted_; // From PrimaryVertexMonitor

  GetLumi* lumiDetails_;

  MonitorElement* NumberOfPVtx;
  MonitorElement* NumberOfPVtxVsBXlumi;
  MonitorElement* NumberOfPVtxVsGoodPVtx;
  MonitorElement* NumberOfGoodPVtx;
  MonitorElement* NumberOfGoodPVtxVsBXlumi;
  MonitorElement* FractionOfGoodPVtx;
  MonitorElement* FractionOfGoodPVtxVsBXlumi;
  MonitorElement* FractionOfGoodPVtxVsGoodPVtx;
  MonitorElement* FractionOfGoodPVtxVsPVtx;
  MonitorElement* NumberOfFakePVtx;
  MonitorElement* NumberOfFakePVtxVsBXlumi;
  MonitorElement* NumberOfFakePVtxVsGoodPVtx;
  MonitorElement* NumberOfBADndofPVtx;
  MonitorElement* NumberOfBADndofPVtxVsBXlumi;
  MonitorElement* NumberOfBADndofPVtxVsGoodPVtx;

  MonitorElement* Chi2oNDFVsGoodPVtx;
  MonitorElement* Chi2oNDFVsBXlumi;
  MonitorElement* Chi2ProbVsGoodPVtx;
  MonitorElement* Chi2ProbVsBXlumi;

  MonitorElement* GoodPVtxSumPt;
  MonitorElement* GoodPVtxSumPtVsBXlumi;
  MonitorElement* GoodPVtxSumPtVsGoodPVtx;

  MonitorElement* GoodPVtxNumberOfTracks;
  MonitorElement* GoodPVtxNumberOfTracksVsBXlumi;
  MonitorElement* GoodPVtxNumberOfTracksVsGoodPVtx;
  MonitorElement* GoodPVtxNumberOfTracksVsGoodPVtxNdof;

  MonitorElement* GoodPVtxChi2oNDFVsGoodPVtx;
  MonitorElement* GoodPVtxChi2oNDFVsBXlumi;
  MonitorElement* GoodPVtxChi2ProbVsGoodPVtx;
  MonitorElement* GoodPVtxChi2ProbVsBXlumi;

  // From PrimaryVertexMonitor
  MonitorElement *nbvtx, *nbgvtx, *nbtksinvtx[2], *trksWeight[2], *score[2];
  MonitorElement *tt[2];
  MonitorElement *xrec[2], *yrec[2], *zrec[2], *xDiff[2], *yDiff[2], *xerr[2], *yerr[2], *zerr[2];
  MonitorElement *xerrVsTrks[2], *yerrVsTrks[2], *zerrVsTrks[2];
  MonitorElement *ntracksVsZ[2];
  MonitorElement *vtxchi2[2], *vtxndf[2], *vtxprob[2], *nans[2];
  MonitorElement *type[2];
  MonitorElement *bsX, *bsY, *bsZ, *bsSigmaZ, *bsDxdz, *bsDydz, *bsBeamWidthX, *bsBeamWidthY, *bsType;

  MonitorElement *sumpt, *ntracks, *weight, *chi2ndf, *chi2prob;
  MonitorElement *dxy, *dxy2, *dz, *dxyErr, *dzErr;
  MonitorElement *dxyVsPhi_pt1, *dzVsPhi_pt1;
  MonitorElement *dxyVsEta_pt1, *dzVsEta_pt1;
  MonitorElement *dxyVsPhi_pt10, *dzVsPhi_pt10;
  MonitorElement *dxyVsEta_pt10, *dzVsEta_pt10;


  bool doAllPlots_;
  bool doPlotsVsBXlumi_;
  bool doPlotsVsGoodPVtx_;

  std::string histname;  //for naming the histograms according to algorithm used
};
#endif
