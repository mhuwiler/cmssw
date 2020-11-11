
/*
 *  See header file for a description of this class.
 *
 *  \author:  Mia Tosi,40 3-B32,+41227671609 
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DQM/TrackingMonitor/interface/VertexMonitor.h"

#include "DQM/TrackingMonitor/interface/GetLumi.h"

// From PrimaryVertexMonitor.cc
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "TMath.h"


using namespace reco;
using namespace edm;


VertexMonitor::VertexMonitor(const edm::ParameterSet& iConfig,
                             const edm::InputTag& primaryVertexInputTag,
                             const edm::InputTag& selectedPrimaryVertexInputTag,
                             std::string pvLabel)
    : conf_(iConfig),
      primaryVertexInputTag_(primaryVertexInputTag),
      selectedPrimaryVertexInputTag_(selectedPrimaryVertexInputTag),
      label_(pvLabel),
      NumberOfPVtx(nullptr),
      NumberOfPVtxVsBXlumi(nullptr),
      NumberOfPVtxVsGoodPVtx(nullptr),
      NumberOfGoodPVtx(nullptr),
      NumberOfGoodPVtxVsBXlumi(nullptr),
      FractionOfGoodPVtx(nullptr),
      FractionOfGoodPVtxVsBXlumi(nullptr),
      FractionOfGoodPVtxVsGoodPVtx(nullptr),
      FractionOfGoodPVtxVsPVtx(nullptr),
      NumberOfBADndofPVtx(nullptr),
      NumberOfBADndofPVtxVsBXlumi(nullptr),
      NumberOfBADndofPVtxVsGoodPVtx(nullptr),
      GoodPVtxSumPt(nullptr),
      GoodPVtxSumPtVsBXlumi(nullptr),
      GoodPVtxSumPtVsGoodPVtx(nullptr),
      GoodPVtxNumberOfTracks(nullptr),
      GoodPVtxNumberOfTracksVsBXlumi(nullptr),
      GoodPVtxNumberOfTracksVsGoodPVtx(nullptr),
      GoodPVtxNumberOfTracksVsGoodPVtxNdof(nullptr),
      GoodPVtxChi2oNDFVsGoodPVtx(nullptr),
      GoodPVtxChi2oNDFVsBXlumi(nullptr),
      GoodPVtxChi2ProbVsGoodPVtx(nullptr),
      GoodPVtxChi2ProbVsBXlumi(nullptr),
      doAllPlots_(conf_.getParameter<bool>("doAllPlots")),
      doPlotsVsBXlumi_(conf_.getParameter<bool>("doPlotsVsBXlumi")),
      doPlotsVsGoodPVtx_(conf_.getParameter<bool>("doPlotsVsGoodPVtx"))

{
  //now do what ever initialization is needed
  if (doPlotsVsBXlumi_)
    lumiDetails_ = new GetLumi(iConfig.getParameter<edm::ParameterSet>("BXlumiSetup"));
}

VertexMonitor::VertexMonitor(const edm::ParameterSet& iConfig,
                             const edm::InputTag& primaryVertexInputTag,
                             const edm::InputTag& selectedPrimaryVertexInputTag,
                             std::string pvLabel,
                             edm::ConsumesCollector& iC)
    : VertexMonitor(iConfig, primaryVertexInputTag, selectedPrimaryVertexInputTag, pvLabel) {
  if (doPlotsVsBXlumi_)
    lumiDetails_ = new GetLumi(iConfig.getParameter<edm::ParameterSet>("BXlumiSetup"), iC);

  pvToken_ = iC.consumes<reco::VertexCollection>(primaryVertexInputTag_);
  selpvToken_ = iC.consumes<reco::VertexCollection>(selectedPrimaryVertexInputTag_);
}

VertexMonitor::~VertexMonitor() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //  if (lumiDetails_) delete lumiDetails_;
}

//
// member functions
//

// -- Analyse
// ------------ method called for each event  ------------
// ------------------------------------------------------- //
void VertexMonitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  double bxlumi = -999.;
  if (doPlotsVsBXlumi_) 
  {
    bxlumi = lumiDetails_->getValue(iEvent);
  }
  //  std::cout << "bxlumi : " << bxlumi << std::endl; // add this as debug output

  size_t totalNumPV = 0;
  size_t totalNumBADndofPV = 0;
  size_t totalNumGoodPV = 0;

  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<reco::VertexCollection> selpvHandle;
  iEvent.getByToken(selpvToken_, selpvHandle);

  Handle<reco::VertexCollection> recVtxs; // Probably to be deleted, as identical to one of the above 
  iEvent.getByToken(vertexToken_, recVtxs);

  Handle<VertexScore> scores;
  iEvent.getByToken(scoreToken_, scores);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  if (not pvHandle.isValid()) return;

  if (recVtxs.isValid() == false || beamSpotHandle.isValid() == false) 
  {
    edm::LogWarning("PrimaryVertexMonitor")
        << " Some products not available in the event: VertexCollection " << vertexInputTag_ << " " << recVtxs.isValid()
        << " BeamSpot " << beamSpotInputTag_ << " " << beamSpotHandle.isValid() << ". Skipping plots for this event";
    return;
  }

  totalNumPV = pvHandle->size();
  //      std::cout << "totalNumPV : " << totalNumPV << std::endl;
  for (reco::VertexCollection::const_iterator pv = pvHandle->begin(); pv != pvHandle->end(); ++pv) 
  {
    //--- count pv w/ ndof < 4
    if (pv->ndof() < 4.)
    {
      totalNumBADndofPV++;
    }
  }
    
  NumberOfPVtx->Fill(totalNumPV);
  NumberOfBADndofPVtx->Fill(totalNumBADndofPV);
  if (doPlotsVsBXlumi_) 
  {
    NumberOfPVtxVsBXlumi->Fill(bxlumi, totalNumPV);
    NumberOfBADndofPVtxVsBXlumi->Fill(bxlumi, totalNumBADndofPV);
  }

  
  if (not selpvHandle.isValid()) return; 
    
  totalNumGoodPV = selpvHandle->size();
  
  //  std::cout << "totalNumGoodPV: " << totalNumGoodPV << std::endl;
  if (doPlotsVsGoodPVtx_) 
  {
    NumberOfPVtxVsGoodPVtx->Fill(totalNumGoodPV, totalNumPV);
    NumberOfBADndofPVtxVsGoodPVtx->Fill(totalNumGoodPV, totalNumBADndofPV);
  }

  double fracGoodPV = static_cast<double>(totalNumGoodPV)/static_cast<double>(totalNumPV);
  //  std::cout << "fracGoodPV: " << fracGoodPV << std::endl;

  NumberOfGoodPVtx->Fill(totalNumGoodPV);
  FractionOfGoodPVtx->Fill(fracGoodPV);
  if (doPlotsVsBXlumi_) {
    NumberOfGoodPVtxVsBXlumi->Fill(bxlumi, totalNumGoodPV);
    FractionOfGoodPVtxVsBXlumi->Fill(bxlumi, fracGoodPV);
  }
  if (doPlotsVsGoodPVtx_) {
    FractionOfGoodPVtxVsGoodPVtx->Fill(totalNumGoodPV, fracGoodPV);
    FractionOfGoodPVtxVsPVtx->Fill(totalNumPV, fracGoodPV);
  }

  if (!selpvHandle->empty()) {
    double sumpt = 0;
    size_t ntracks = 0;
    double chi2ndf = 0.;
    double chi2prob = 0.;

    if (!selpvHandle->at(0).isFake()) {
      reco::Vertex pv = selpvHandle->at(0);

      ntracks = pv.tracksSize();
      chi2ndf = pv.normalizedChi2();
      chi2prob = TMath::Prob(pv.chi2(), (int)pv.ndof());

      for (reco::Vertex::trackRef_iterator itrk = pv.tracks_begin(); itrk != pv.tracks_end(); ++itrk) {
        double pt = (**itrk).pt();
        sumpt += pt * pt;
      }
      GoodPVtxSumPt->Fill(sumpt);
      GoodPVtxNumberOfTracks->Fill(ntracks);

      if (doPlotsVsBXlumi_) {
        GoodPVtxSumPtVsBXlumi->Fill(bxlumi, sumpt);
        GoodPVtxNumberOfTracksVsBXlumi->Fill(bxlumi, ntracks);
        GoodPVtxChi2oNDFVsBXlumi->Fill(bxlumi, chi2ndf);
        GoodPVtxChi2ProbVsBXlumi->Fill(bxlumi, chi2prob);
      }
      if (doPlotsVsGoodPVtx_) {
        GoodPVtxSumPtVsGoodPVtx->Fill(totalNumGoodPV, sumpt);
        GoodPVtxNumberOfTracksVsGoodPVtx->Fill(totalNumGoodPV, ntracks);
        GoodPVtxChi2oNDFVsGoodPVtx->Fill(totalNumGoodPV, chi2ndf);
        GoodPVtxChi2ProbVsGoodPVtx->Fill(totalNumGoodPV, chi2prob);
      }
    }
  }
}

// From PrimaryVertexMonitor
void VertexMonitor::pvTracksPlots(const Vertex& v) 
{
  if (!v.isValid())
    return;
  if (v.isFake())
    return;

  if (v.tracksSize() == 0) {
    ntracks->Fill(0);
    return;
  }

  const math::XYZPoint myVertex(v.position().x(), v.position().y(), v.position().z());

  size_t nTracks = 0;
  float sumPT = 0.;
  const int cmToUm = 10000;

  for (reco::Vertex::trackRef_iterator t = v.tracks_begin(); t != v.tracks_end(); t++) {
    bool isHighPurity = (**t).quality(reco::TrackBase::highPurity);
    if (!isHighPurity)
      continue;

    float pt = (**t).pt();
    if (pt < 1.)
      continue;

    nTracks++;

    float eta = (**t).eta();
    float phi = (**t).phi();

    float w = v.trackWeight(*t);
    float chi2NDF = (**t).normalizedChi2();
    float chi2Prob = TMath::Prob((**t).chi2(), (int)(**t).ndof());
    float Dxy = (**t).dxy(myVertex) * cmToUm;  // is it needed ?
    float Dz = (**t).dz(myVertex) * cmToUm;    // is it needed ?
    float DxyErr = (**t).dxyError() * cmToUm;
    float DzErr = (**t).dzError() * cmToUm;

    sumPT += pt * pt;

    // fill MEs
    weight->Fill(w);
    chi2ndf->Fill(chi2NDF);
    chi2prob->Fill(chi2Prob);
    dxy->Fill(Dxy);
    dxy2->Fill(Dxy);
    dz->Fill(Dz);
    dxyErr->Fill(DxyErr);
    dzErr->Fill(DzErr);

    dxyVsPhi_pt1->Fill(phi, Dxy);
    dzVsPhi_pt1->Fill(phi, Dz);
    dxyVsEta_pt1->Fill(eta, Dxy);
    dzVsEta_pt1->Fill(eta, Dz);

    if (pt < 10.)
      continue;
    dxyVsPhi_pt10->Fill(phi, Dxy);
    dzVsPhi_pt10->Fill(phi, Dz);
    dxyVsEta_pt10->Fill(eta, Dxy);
    dzVsEta_pt10->Fill(eta, Dz);
  }
  ntracks->Fill(float(nTracks));
  sumpt->Fill(sumPT);
}

// From PrimaryVertexMonitor
void VertexMonitor::vertexPlots(const Vertex& v, const BeamSpot& beamSpot, int i) 
{
  if (i < 0 || i > 1)
    return;
  if (!v.isValid())
    type[i]->Fill(2.);
  else if (v.isFake())
    type[i]->Fill(1.);
  else
    type[i]->Fill(0.);

  if (v.isValid() && !v.isFake()) {
    float weight = 0;
    for (reco::Vertex::trackRef_iterator t = v.tracks_begin(); t != v.tracks_end(); t++)
      weight += v.trackWeight(*t);
    trksWeight[i]->Fill(weight);
    nbtksinvtx[i]->Fill(v.tracksSize());
    ntracksVsZ[i]->Fill(v.position().z() - beamSpot.z0(), v.tracksSize());

    vtxchi2[i]->Fill(v.chi2());
    vtxndf[i]->Fill(v.ndof());
    vtxprob[i]->Fill(ChiSquaredProbability(v.chi2(), v.ndof()));

    xrec[i]->Fill(v.position().x());
    yrec[i]->Fill(v.position().y());
    zrec[i]->Fill(v.position().z());

    float xb = beamSpot.x0() + beamSpot.dxdz() * (v.position().z() - beamSpot.z0());
    float yb = beamSpot.y0() + beamSpot.dydz() * (v.position().z() - beamSpot.z0());
    xDiff[i]->Fill((v.position().x() - xb) * 10000);
    yDiff[i]->Fill((v.position().y() - yb) * 10000);

    xerr[i]->Fill(v.xError() * 10000);
    yerr[i]->Fill(v.yError() * 10000);
    zerr[i]->Fill(v.zError() * 10000);
    xerrVsTrks[i]->Fill(weight, v.xError() * 10000);
    yerrVsTrks[i]->Fill(weight, v.yError() * 10000);
    zerrVsTrks[i]->Fill(weight, v.zError() * 10000);

    nans[i]->Fill(1., edm::isNotFinite(v.position().x()) * 1.);
    nans[i]->Fill(2., edm::isNotFinite(v.position().y()) * 1.);
    nans[i]->Fill(3., edm::isNotFinite(v.position().z()) * 1.);

    int index = 3;
    for (int k = 0; k != 3; k++) {
      for (int j = k; j != 3; j++) {
        index++;
        nans[i]->Fill(index * 1., edm::isNotFinite(v.covariance(k, j)) * 1.);
        // in addition, diagonal element must be positive
        if (j == k && v.covariance(k, j) < 0) {
          nans[i]->Fill(index * 1., 1.);
        }
      }
    }
  }
}

// ------------ method called once each job just before starting event loop  ------------
void VertexMonitor::initHisto(DQMStore::IBooker& ibooker) {
  // parameters from the configuration
  std::string MEFolderName = conf_.getParameter<std::string>("PVFolderName");

  // get binning from the configuration
  edm::ParameterSet ParametersGoodPVtx = conf_.getParameter<edm::ParameterSet>("GoodPVtx");
  int GoodPVtxBin = ParametersGoodPVtx.getParameter<int>("GoodPVtxBin");
  double GoodPVtxMin = ParametersGoodPVtx.getParameter<double>("GoodPVtxMin");
  double GoodPVtxMax = ParametersGoodPVtx.getParameter<double>("GoodPVtxMax");

  edm::ParameterSet ParametersNTrkPVtx = conf_.getParameter<edm::ParameterSet>("NTrkPVtx");
  int NTrkPVtxBin = ParametersNTrkPVtx.getParameter<int>("NTrkPVtxBin");
  double NTrkPVtxMin = ParametersNTrkPVtx.getParameter<double>("NTrkPVtxMin");
  double NTrkPVtxMax = ParametersNTrkPVtx.getParameter<double>("NTrkPVtxMax");

  edm::ParameterSet ParametersSumPtPVtx = conf_.getParameter<edm::ParameterSet>("SumPtPVtx");
  int SumPtPVtxBin = ParametersSumPtPVtx.getParameter<int>("SumPtPVtxBin");
  double SumPtPVtxMin = ParametersSumPtPVtx.getParameter<double>("SumPtPVtxMin");
  double SumPtPVtxMax = ParametersSumPtPVtx.getParameter<double>("SumPtPVtxMax");

  // book histo
  // ----------------------//
  ibooker.setCurrentFolder(MEFolderName + "/" + label_);

  histname = "NumberOfPVtx_" + label_;
  NumberOfPVtx = ibooker.book1D(histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax);
  NumberOfPVtx->setAxisTitle("Number of PV", 1);
  NumberOfPVtx->setAxisTitle("Number of Events", 2);

  histname = "NumberOfGoodPVtx_" + label_;
  NumberOfGoodPVtx = ibooker.book1D(histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax);
  NumberOfGoodPVtx->setAxisTitle("Number of Good PV", 1);
  NumberOfGoodPVtx->setAxisTitle("Number of Events", 2);

  histname = "FractionOfGoodPVtx_" + label_;
  FractionOfGoodPVtx = ibooker.book1D(histname, histname, 100, 0., 1.);
  FractionOfGoodPVtx->setAxisTitle("fraction of Good PV", 1);
  FractionOfGoodPVtx->setAxisTitle("Number of Events", 2);

  histname = "NumberOfBADndofPVtx_" + label_;
  NumberOfBADndofPVtx = ibooker.book1D(histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax);
  NumberOfBADndofPVtx->setAxisTitle("Number of BADndof #PV", 1);
  NumberOfBADndofPVtx->setAxisTitle("Number of Events", 2);

  histname = "GoodPVtxSumPt_" + label_;
  GoodPVtxSumPt = ibooker.book1D(histname, histname, SumPtPVtxBin, SumPtPVtxMin, SumPtPVtxMax);
  GoodPVtxSumPt->setAxisTitle("primary vertex #Sum p_{T}^{2} [GeV^{2}/c^{2}]", 1);
  GoodPVtxSumPt->setAxisTitle("Number of events", 2);

  histname = "GoodPVtxNumberOfTracks_" + label_;
  GoodPVtxNumberOfTracks = ibooker.book1D(histname, histname, NTrkPVtxBin, NTrkPVtxMin, NTrkPVtxMax);
  GoodPVtxNumberOfTracks->setAxisTitle("primary vertex number of tracks", 1);
  GoodPVtxNumberOfTracks->setAxisTitle("Number of events", 2);

  //histname = 

  if (doPlotsVsBXlumi_) {
    // get binning from the configuration
    edm::ParameterSet BXlumiParameters = conf_.getParameter<edm::ParameterSet>("BXlumiSetup");
    int BXlumiBin = BXlumiParameters.getParameter<int>("BXlumiBin");
    double BXlumiMin = BXlumiParameters.getParameter<double>("BXlumiMin");
    double BXlumiMax = BXlumiParameters.getParameter<double>("BXlumiMax");

    ibooker.setCurrentFolder(MEFolderName + "/" + label_ + "/PUmonitoring/");

    histname = "NumberOfPVtxVsBXlumi_" + label_;
    NumberOfPVtxVsBXlumi =
        ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, GoodPVtxMin, GoodPVtxMax * 3, "");
    NumberOfPVtxVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    NumberOfPVtxVsBXlumi->setAxisTitle("lumi BX [10^{30}Hzcm^{-2}]", 1);
    NumberOfPVtxVsBXlumi->setAxisTitle("Mean number of PV", 2);

    histname = "NumberOfGoodPVtxVsBXlumi_" + label_;
    NumberOfGoodPVtxVsBXlumi =
        ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, GoodPVtxMin, GoodPVtxMax * 3, "");
    NumberOfGoodPVtxVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    NumberOfGoodPVtxVsBXlumi->setAxisTitle("lumi BX [10^{30}Hzcm^{-2}]", 1);
    NumberOfGoodPVtxVsBXlumi->setAxisTitle("Mean number of PV", 2);

    histname = "FractionOfGoodPVtxVsBXlumi_" + label_;
    FractionOfGoodPVtxVsBXlumi = ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, 0., 1.5, "");
    FractionOfGoodPVtxVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    FractionOfGoodPVtxVsBXlumi->setAxisTitle("lumi BX [10^{30}Hzcm^{-2}]", 1);
    FractionOfGoodPVtxVsBXlumi->setAxisTitle("Mean number of PV", 2);

    histname = "NumberOfBADndofPVtxVsBXlumi_" + label_;
    NumberOfBADndofPVtxVsBXlumi =
        ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, GoodPVtxMin, GoodPVtxMax * 3, "");
    NumberOfBADndofPVtxVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    NumberOfBADndofPVtxVsBXlumi->setAxisTitle("BADndof #PV", 1);
    NumberOfBADndofPVtxVsBXlumi->setAxisTitle("Number of Events", 2);

    histname = "GoodPVtxSumPtVsBXlumi_" + label_;
    GoodPVtxSumPtVsBXlumi =
        ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, SumPtPVtxMin, SumPtPVtxMax * 3, "");
    GoodPVtxSumPtVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    GoodPVtxSumPtVsBXlumi->setAxisTitle("lumi BX [10^{30}Hzcm^{-2}]", 1);
    GoodPVtxSumPtVsBXlumi->setAxisTitle("Mean pv #Sum p_{T}^{2} [GeV^{2}/c]^{2}", 2);

    histname = "GoodPVtxNumberOfTracksVsBXlumi_" + label_;
    GoodPVtxNumberOfTracksVsBXlumi =
        ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, NTrkPVtxMin, NTrkPVtxMax * 3, "");
    GoodPVtxNumberOfTracksVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    GoodPVtxNumberOfTracksVsBXlumi->setAxisTitle("lumi BX [10^{30}Hzcm^{-2}]", 1);
    GoodPVtxNumberOfTracksVsBXlumi->setAxisTitle("Mean pv number of tracks", 2);

    // get binning from the configuration
    double Chi2NDFMin = conf_.getParameter<double>("Chi2NDFMin");
    double Chi2NDFMax = conf_.getParameter<double>("Chi2NDFMax");

    double Chi2ProbMin = conf_.getParameter<double>("Chi2ProbMin");
    double Chi2ProbMax = conf_.getParameter<double>("Chi2ProbMax");

    histname = "Chi2oNDFVsBXlumi_" + label_;
    Chi2oNDFVsBXlumi =
        ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, Chi2NDFMin, Chi2NDFMax * 3, "");
    Chi2oNDFVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    Chi2oNDFVsBXlumi->setAxisTitle("lumi BX [10^{30}Hzcm^{-2}]", 1);
    Chi2oNDFVsBXlumi->setAxisTitle("Mean #chi^{2}/ndof", 2);

    histname = "Chi2ProbVsBXlumi_" + label_;
    Chi2ProbVsBXlumi =
        ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, Chi2ProbMin, Chi2ProbMax * 3, "");
    Chi2ProbVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    Chi2ProbVsBXlumi->setAxisTitle("lumi BX [10^{30}Hzcm^{-2}]", 1);
    Chi2ProbVsBXlumi->setAxisTitle("Mean #chi^{2}/prob", 2);

    histname = "GoodPVtxChi2oNDFVsBXlumi_" + label_;
    GoodPVtxChi2oNDFVsBXlumi =
        ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, Chi2NDFMin, Chi2NDFMax, "");
    GoodPVtxChi2oNDFVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    GoodPVtxChi2oNDFVsBXlumi->setAxisTitle("lumi BX [10^{30}Hzcm^{-2}]", 1);
    GoodPVtxChi2oNDFVsBXlumi->setAxisTitle("Mean PV #chi^{2}/ndof", 2);

    histname = "GoodPVtxChi2ProbVsBXlumi_" + label_;
    GoodPVtxChi2ProbVsBXlumi =
        ibooker.bookProfile(histname, histname, BXlumiBin, BXlumiMin, BXlumiMax, Chi2ProbMin, Chi2ProbMax * 3, "");
    GoodPVtxChi2ProbVsBXlumi->getTH1()->SetCanExtend(TH1::kAllAxes);
    GoodPVtxChi2ProbVsBXlumi->setAxisTitle("lumi BX [10^{30}Hzcm^{-2}]", 1);
    GoodPVtxChi2ProbVsBXlumi->setAxisTitle("Mean PV #chi^{2}/prob", 2);
  }

  if (doPlotsVsGoodPVtx_) {
    ibooker.setCurrentFolder(MEFolderName + "/" + label_ + "/PUmonitoring/VsGoodPVtx");

    histname = "NumberOfPVtxVsGoodPVtx_" + label_;
    NumberOfPVtxVsGoodPVtx = ibooker.bookProfile(
        histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax, GoodPVtxMin, GoodPVtxMax * 3, "");
    NumberOfPVtxVsGoodPVtx->getTH1()->SetCanExtend(TH1::kAllAxes);
    NumberOfPVtxVsGoodPVtx->setAxisTitle("Number of Good PV", 1);
    NumberOfPVtxVsGoodPVtx->setAxisTitle("Mean number of PV", 2);

    histname = "FractionOfGoodPVtxVsGoodPVtx_" + label_;
    FractionOfGoodPVtxVsGoodPVtx = ibooker.bookProfile(
        histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax, GoodPVtxMin, GoodPVtxMax * 3, "");
    FractionOfGoodPVtxVsGoodPVtx->getTH1()->SetCanExtend(TH1::kAllAxes);
    FractionOfGoodPVtxVsGoodPVtx->setAxisTitle("Number of Good PV", 1);
    FractionOfGoodPVtxVsGoodPVtx->setAxisTitle("Mean fraction of Good PV", 2);

    histname = "FractionOfGoodPVtxVsPVtx_" + label_;
    FractionOfGoodPVtxVsPVtx = ibooker.bookProfile(
        histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax, GoodPVtxMin, GoodPVtxMax * 3, "");
    FractionOfGoodPVtxVsPVtx->getTH1()->SetCanExtend(TH1::kAllAxes);
    FractionOfGoodPVtxVsPVtx->setAxisTitle("Number of Good PV", 1);
    FractionOfGoodPVtxVsPVtx->setAxisTitle("Mean number of Good PV", 2);

    histname = "NumberOfBADndofPVtxVsGoodPVtx_" + label_;
    NumberOfBADndofPVtxVsGoodPVtx = ibooker.bookProfile(
        histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax, GoodPVtxMin, GoodPVtxMax * 3, "");
    NumberOfBADndofPVtxVsGoodPVtx->getTH1()->SetCanExtend(TH1::kAllAxes);
    NumberOfBADndofPVtxVsGoodPVtx->setAxisTitle("Number of Good PV", 1);
    NumberOfBADndofPVtxVsGoodPVtx->setAxisTitle("Mean Number of BAD PV", 2);

    histname = "GoodPVtxSumPtVsGoodPVtx_" + label_;
    GoodPVtxSumPtVsGoodPVtx = ibooker.bookProfile(
        histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax, SumPtPVtxMin, SumPtPVtxMax * 3, "");
    GoodPVtxSumPtVsGoodPVtx->getTH1()->SetCanExtend(TH1::kAllAxes);
    GoodPVtxSumPtVsGoodPVtx->setAxisTitle("Number of Good PV", 1);
    GoodPVtxSumPtVsGoodPVtx->setAxisTitle("Mean pv #Sum p_{T}^{2} [GeV^{2}/c]^{2}", 2);

    histname = "GoodPVtxNumberOfTracksVsGoodPVtx_" + label_;
    GoodPVtxNumberOfTracksVsGoodPVtx = ibooker.bookProfile(
        histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax, NTrkPVtxMin, NTrkPVtxMax * 3, "");
    GoodPVtxNumberOfTracksVsGoodPVtx->getTH1()->SetCanExtend(TH1::kAllAxes);
    GoodPVtxNumberOfTracksVsGoodPVtx->setAxisTitle("Number of Good PV", 1);
    GoodPVtxNumberOfTracksVsGoodPVtx->setAxisTitle("Mean pv number of tracks", 2);

    // get binning from the configuration
    double Chi2NDFMin = conf_.getParameter<double>("Chi2NDFMin");
    double Chi2NDFMax = conf_.getParameter<double>("Chi2NDFMax");

    double Chi2ProbMin = conf_.getParameter<double>("Chi2ProbMin");
    double Chi2ProbMax = conf_.getParameter<double>("Chi2ProbMax");

    histname = "GoodPVtxChi2oNDFVsGoodPVtx_" + label_;
    GoodPVtxChi2oNDFVsGoodPVtx =
        ibooker.bookProfile(histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax, Chi2NDFMin, Chi2NDFMax * 3, "");
    GoodPVtxChi2oNDFVsGoodPVtx->getTH1()->SetCanExtend(TH1::kAllAxes);
    GoodPVtxChi2oNDFVsGoodPVtx->setAxisTitle("Number of Good PV", 1);
    GoodPVtxChi2oNDFVsGoodPVtx->setAxisTitle("Mean PV #chi^{2}/ndof", 2);

    histname = "GoodPVtxChi2ProbVsGoodPVtx_" + label_;
    GoodPVtxChi2ProbVsGoodPVtx = ibooker.bookProfile(
        histname, histname, GoodPVtxBin, GoodPVtxMin, GoodPVtxMax, Chi2ProbMin, Chi2ProbMax * 3, "");
    GoodPVtxChi2ProbVsGoodPVtx->getTH1()->SetCanExtend(TH1::kAllAxes);
    GoodPVtxChi2ProbVsGoodPVtx->setAxisTitle("Number of Good PV", 1);
    GoodPVtxChi2ProbVsGoodPVtx->setAxisTitle("Mean PV #chi^{2}/prob", 2);
  }
}

void VertexMonitor::bookHistograms(DQMStore::IBooker& iBooker, edm::Run const&, edm::EventSetup const&) { // Why those additional arguments? 
  std::string dqmLabel = "";

  //
  // Book all histograms.
  //

  //  get the store
  dqmLabel = TopFolderName_ + "/" + vertexInputTag_.label();
  iBooker.setCurrentFolder(dqmLabel);

  //   xPos = iBooker.book1D ("xPos","x Coordinate" ,100, -0.1, 0.1);

  nbvtx = iBooker.book1D("vtxNbr", "Reconstructed Vertices in Event", 80, -0.5, 79.5);
  nbgvtx = iBooker.book1D("goodvtxNbr", "Reconstructed Good Vertices in Event", 80, -0.5, 79.5);

  // to be configured each year...
  auto vposx = conf_.getParameter<double>("Xpos");
  auto vposy = conf_.getParameter<double>("Ypos");

  nbtksinvtx[0] = iBooker.book1D("otherVtxTrksNbr", "Reconstructed Tracks in Vertex (other Vtx)", 40, -0.5, 99.5);
  ntracksVsZ[0] = iBooker.bookProfile(
      "otherVtxTrksVsZ", "Reconstructed Tracks in Vertex (other Vtx) vs Z", 80, -20., 20., 50, 0, 100, "");
  ntracksVsZ[0]->setAxisTitle("z-bs", 1);
  ntracksVsZ[0]->setAxisTitle("#tracks", 2);

  score[0] = iBooker.book1D("otherVtxScore", "sqrt(score) (other Vtx)", 100, 0., 400.);
  trksWeight[0] = iBooker.book1D("otherVtxTrksWeight", "Total weight of Tracks in Vertex (other Vtx)", 40, 0, 100.);
  vtxchi2[0] = iBooker.book1D("otherVtxChi2", "#chi^{2} (other Vtx)", 100, 0., 200.);
  vtxndf[0] = iBooker.book1D("otherVtxNdf", "ndof (other Vtx)", 100, 0., 200.);
  vtxprob[0] = iBooker.book1D("otherVtxProb", "#chi^{2} probability (other Vtx)", 100, 0., 1.);
  nans[0] = iBooker.book1D("otherVtxNans", "Illegal values for x,y,z,xx,xy,xz,yy,yz,zz (other Vtx)", 9, 0.5, 9.5);

  nbtksinvtx[1] = iBooker.book1D("tagVtxTrksNbr", "Reconstructed Tracks in Vertex (tagged Vtx)", 100, -0.5, 99.5);
  ntracksVsZ[1] = iBooker.bookProfile(
      "tagVtxTrksVsZ", "Reconstructed Tracks in Vertex (tagged Vtx) vs Z", 80, -20., 20., 50, 0, 100, "");
  ntracksVsZ[1]->setAxisTitle("z-bs", 1);
  ntracksVsZ[1]->setAxisTitle("#tracks", 2);

  score[1] = iBooker.book1D("tagVtxScore", "sqrt(score) (tagged Vtx)", 100, 0., 400.);
  trksWeight[1] = iBooker.book1D("tagVtxTrksWeight", "Total weight of Tracks in Vertex (tagged Vtx)", 100, 0, 100.);
  vtxchi2[1] = iBooker.book1D("tagVtxChi2", "#chi^{2} (tagged Vtx)", 100, 0., 200.);
  vtxndf[1] = iBooker.book1D("tagVtxNdf", "ndof (tagged Vtx)", 100, 0., 200.);
  vtxprob[1] = iBooker.book1D("tagVtxProb", "#chi^{2} probability (tagged Vtx)", 100, 0., 1.);
  nans[1] = iBooker.book1D("tagVtxNans", "Illegal values for x,y,z,xx,xy,xz,yy,yz,zz (tagged Vtx)", 9, 0.5, 9.5);

  xrec[0] = iBooker.book1D("otherPosX", "Position x Coordinate (other Vtx)", 100, vposx - 0.1, vposx + 0.1);
  yrec[0] = iBooker.book1D("otherPosY", "Position y Coordinate (other Vtx)", 100, vposy - 0.1, vposy + 0.1);
  zrec[0] = iBooker.book1D("otherPosZ", "Position z Coordinate (other Vtx)", 100, -20., 20.);
  xDiff[0] = iBooker.book1D("otherDiffX", "X distance from BeamSpot (other Vtx)", 100, -500, 500);
  yDiff[0] = iBooker.book1D("otherDiffY", "Y distance from BeamSpot (other Vtx)", 100, -500, 500);
  xerr[0] = iBooker.book1D("otherErrX", "Uncertainty x Coordinate (other Vtx)", 100, 0., 100);
  yerr[0] = iBooker.book1D("otherErrY", "Uncertainty y Coordinate (other Vtx)", 100, 0., 100);
  zerr[0] = iBooker.book1D("otherErrZ", "Uncertainty z Coordinate (other Vtx)", 100, 0., 100);
  xerrVsTrks[0] = iBooker.book2D(
      "otherErrVsWeightX", "Uncertainty x Coordinate vs. track weight (other Vtx)", 100, 0, 100., 100, 0., 100);
  yerrVsTrks[0] = iBooker.book2D(
      "otherErrVsWeightY", "Uncertainty y Coordinate vs. track weight (other Vtx)", 100, 0, 100., 100, 0., 100);
  zerrVsTrks[0] = iBooker.book2D(
      "otherErrVsWeightZ", "Uncertainty z Coordinate vs. track weight (other Vtx)", 100, 0, 100., 100, 0., 100);

  xrec[1] = iBooker.book1D("tagPosX", "Position x Coordinate (tagged Vtx)", 100, vposx - 0.1, vposx + 0.1);
  yrec[1] = iBooker.book1D("tagPosY", "Position y Coordinate (tagged Vtx)", 100, vposy - 0.1, vposy + 0.1);
  zrec[1] = iBooker.book1D("tagPosZ", "Position z Coordinate (tagged Vtx)", 100, -20., 20.);
  xDiff[1] = iBooker.book1D("tagDiffX", "X distance from BeamSpot (tagged Vtx)", 100, -500, 500);
  yDiff[1] = iBooker.book1D("tagDiffY", "Y distance from BeamSpot (tagged Vtx)", 100, -500, 500);
  xerr[1] = iBooker.book1D("tagErrX", "Uncertainty x Coordinate (tagged Vtx)", 100, 0., 100);
  yerr[1] = iBooker.book1D("tagErrY", "Uncertainty y Coordinate (tagged Vtx)", 100, 0., 100);
  zerr[1] = iBooker.book1D("tagErrZ", "Uncertainty z Coordinate (tagged Vtx)", 100, 0., 100);
  xerrVsTrks[1] = iBooker.book2D(
      "tagErrVsWeightX", "Uncertainty x Coordinate vs. track weight (tagged Vtx)", 100, 0, 100., 100, 0., 100);
  yerrVsTrks[1] = iBooker.book2D(
      "tagErrVsWeightY", "Uncertainty y Coordinate vs. track weight (tagged Vtx)", 100, 0, 100., 100, 0., 100);
  zerrVsTrks[1] = iBooker.book2D(
      "tagErrVsWeightZ", "Uncertainty z Coordinate vs. track weight (tagged Vtx)", 100, 0, 100., 100, 0., 100);

  type[0] = iBooker.book1D("otherType", "Vertex type (other Vtx)", 3, -0.5, 2.5);
  type[1] = iBooker.book1D("tagType", "Vertex type (tagged Vtx)", 3, -0.5, 2.5);
  for (int i = 0; i < 2; ++i) {
    type[i]->setBinLabel(1, "Valid, real");
    type[i]->setBinLabel(2, "Valid, fake");
    type[i]->setBinLabel(3, "Invalid");
  }

  //  get the store
  dqmLabel = TopFolderName_ + "/" + beamSpotInputTag_.label();
  iBooker.setCurrentFolder(dqmLabel);

  bsX = iBooker.book1D("bsX", "BeamSpot x0", 100, -0.1, 0.1);
  bsY = iBooker.book1D("bsY", "BeamSpot y0", 100, -0.1, 0.1);
  bsZ = iBooker.book1D("bsZ", "BeamSpot z0", 100, -2., 2.);
  bsSigmaZ = iBooker.book1D("bsSigmaZ", "BeamSpot sigmaZ", 100, 0., 10.);
  bsDxdz = iBooker.book1D("bsDxdz", "BeamSpot dxdz", 100, -0.0003, 0.0003);
  bsDydz = iBooker.book1D("bsDydz", "BeamSpot dydz", 100, -0.0003, 0.0003);
  bsBeamWidthX = iBooker.book1D("bsBeamWidthX", "BeamSpot BeamWidthX", 100, 0., 100.);
  bsBeamWidthY = iBooker.book1D("bsBeamWidthY", "BeamSpot BeamWidthY", 100, 0., 100.);
  bsType = iBooker.book1D("bsType", "BeamSpot type", 4, -1.5, 2.5);
  bsType->setBinLabel(1, "Unknown");
  bsType->setBinLabel(2, "Fake");
  bsType->setBinLabel(3, "LHC");
  bsType->setBinLabel(4, "Tracker");

  //  get the store
  dqmLabel = TopFolderName_ + "/" + AlignmentLabel_;
  iBooker.setCurrentFolder(dqmLabel);

  int TKNoBin = conf_.getParameter<int>("TkSizeBin");
  double TKNoMin = conf_.getParameter<double>("TkSizeMin");
  double TKNoMax = conf_.getParameter<double>("TkSizeMax");

  int DxyBin = conf_.getParameter<int>("DxyBin");
  double DxyMin = conf_.getParameter<double>("DxyMin");
  double DxyMax = conf_.getParameter<double>("DxyMax");

  int DzBin = conf_.getParameter<int>("DzBin");
  double DzMin = conf_.getParameter<double>("DzMin");
  double DzMax = conf_.getParameter<double>("DzMax");

  int PhiBin = conf_.getParameter<int>("PhiBin");
  double PhiMin = conf_.getParameter<double>("PhiMin");
  double PhiMax = conf_.getParameter<double>("PhiMax");

  int EtaBin = conf_.getParameter<int>("EtaBin");
  double EtaMin = conf_.getParameter<double>("EtaMin");
  double EtaMax = conf_.getParameter<double>("EtaMax");

  ntracks = iBooker.book1D("ntracks", "number of PV tracks (p_{T} > 1 GeV)", TKNoBin, TKNoMin, TKNoMax);
  ntracks->setAxisTitle("Number of PV Tracks (p_{T} > 1 GeV) per Event", 1);
  ntracks->setAxisTitle("Number of Event", 2);

  weight = iBooker.book1D("weight", "weight of PV tracks (p_{T} > 1 GeV)", 100, 0., 1.);
  weight->setAxisTitle("weight of PV Tracks (p_{T} > 1 GeV) per Event", 1);
  weight->setAxisTitle("Number of Event", 2);

  sumpt = iBooker.book1D("sumpt", "#Sum p_{T} of PV tracks (p_{T} > 1 GeV)", 100, -0.5, 249.5);
  chi2ndf = iBooker.book1D("chi2ndf", "PV tracks (p_{T} > 1 GeV) #chi^{2}/ndof", 100, 0., 20.);
  chi2prob = iBooker.book1D("chi2prob", "PV tracks (p_{T} > 1 GeV) #chi^{2} probability", 100, 0., 1.);

  dxy = iBooker.book1D("dxy", "PV tracks (p_{T} > 1 GeV) d_{xy} (#mum)", DxyBin, DxyMin, DxyMax);
  dxy2 = iBooker.book1D("dxyzoom", "PV tracks (p_{T} > 1 GeV) d_{xy} (#mum)", DxyBin, DxyMin / 5., DxyMax / 5.);
  dxyErr = iBooker.book1D("dxyErr", "PV tracks (p_{T} > 1 GeV) d_{xy} error (#mum)", 100, 0., 2000.);
  dz = iBooker.book1D("dz", "PV tracks (p_{T} > 1 GeV) d_{z} (#mum)", DzBin, DzMin, DzMax);
  dzErr = iBooker.book1D("dzErr", "PV tracks (p_{T} > 1 GeV) d_{z} error(#mum)", 100, 0., 10000.);

  dxyVsPhi_pt1 = iBooker.bookProfile("dxyVsPhi_pt1",
                                     "PV tracks (p_{T} > 1 GeV) d_{xy} (#mum) VS track #phi",
                                     PhiBin,
                                     PhiMin,
                                     PhiMax,
                                     DxyBin,
                                     DxyMin,
                                     DxyMax,
                                     "");
  dxyVsPhi_pt1->setAxisTitle("PV track (p_{T} > 1 GeV) #phi", 1);
  dxyVsPhi_pt1->setAxisTitle("PV track (p_{T} > 1 GeV) d_{xy} (#mum)", 2);

  dzVsPhi_pt1 = iBooker.bookProfile("dzVsPhi_pt1",
                                    "PV tracks (p_{T} > 1 GeV) d_{z} (#mum) VS track #phi",
                                    PhiBin,
                                    PhiMin,
                                    PhiMax,
                                    DzBin,
                                    DzMin,
                                    DzMax,
                                    "");
  dzVsPhi_pt1->setAxisTitle("PV track (p_{T} > 1 GeV) #phi", 1);
  dzVsPhi_pt1->setAxisTitle("PV track (p_{T} > 1 GeV) d_{z} (#mum)", 2);

  dxyVsEta_pt1 = iBooker.bookProfile("dxyVsEta_pt1",
                                     "PV tracks (p_{T} > 1 GeV) d_{xy} (#mum) VS track #eta",
                                     EtaBin,
                                     EtaMin,
                                     EtaMax,
                                     DxyBin,
                                     DxyMin,
                                     DxyMax,
                                     "");
  dxyVsEta_pt1->setAxisTitle("PV track (p_{T} > 1 GeV) #eta", 1);
  dxyVsEta_pt1->setAxisTitle("PV track (p_{T} > 1 GeV) d_{xy} (#mum)", 2);

  dzVsEta_pt1 = iBooker.bookProfile("dzVsEta_pt1",
                                    "PV tracks (p_{T} > 1 GeV) d_{z} (#mum) VS track #eta",
                                    EtaBin,
                                    EtaMin,
                                    EtaMax,
                                    DzBin,
                                    DzMin,
                                    DzMax,
                                    "");
  dzVsEta_pt1->setAxisTitle("PV track (p_{T} > 1 GeV) #eta", 1);
  dzVsEta_pt1->setAxisTitle("PV track (p_{T} > 1 GeV) d_{z} (#mum)", 2);

  dxyVsPhi_pt10 = iBooker.bookProfile("dxyVsPhi_pt10",
                                      "PV tracks (p_{T} > 10 GeV) d_{xy} (#mum) VS track #phi",
                                      PhiBin,
                                      PhiMin,
                                      PhiMax,
                                      DxyBin,
                                      DxyMin,
                                      DxyMax,
                                      "");
  dxyVsPhi_pt10->setAxisTitle("PV track (p_{T} > 10 GeV) #phi", 1);
  dxyVsPhi_pt10->setAxisTitle("PV track (p_{T} > 10 GeV) d_{xy} (#mum)", 2);

  dzVsPhi_pt10 = iBooker.bookProfile("dzVsPhi_pt10",
                                     "PV tracks (p_{T} > 10 GeV) d_{z} (#mum) VS track #phi",
                                     PhiBin,
                                     PhiMin,
                                     PhiMax,
                                     DzBin,
                                     DzMin,
                                     DzMax,
                                     "");
  dzVsPhi_pt10->setAxisTitle("PV track (p_{T} > 10 GeV) #phi", 1);
  dzVsPhi_pt10->setAxisTitle("PV track (p_{T} > 10 GeV) d_{z} (#mum)", 2);

  dxyVsEta_pt10 = iBooker.bookProfile("dxyVsEta_pt10",
                                      "PV tracks (p_{T} > 10 GeV) d_{xy} (#mum) VS track #eta",
                                      EtaBin,
                                      EtaMin,
                                      EtaMax,
                                      DxyBin,
                                      DxyMin,
                                      DxyMax,
                                      "");
  dxyVsEta_pt10->setAxisTitle("PV track (p_{T} > 10 GeV) #eta", 1);
  dxyVsEta_pt10->setAxisTitle("PV track (p_{T} > 10 GeV) d_{xy} (#mum)", 2);

  dzVsEta_pt10 = iBooker.bookProfile("dzVsEta_pt10",
                                     "PV tracks (p_{T} > 10 GeV) d_{z} (#mum) VS track #eta",
                                     EtaBin,
                                     EtaMin,
                                     EtaMax,
                                     DzBin,
                                     DzMin,
                                     DzMax,
                                     "");
  dzVsEta_pt10->setAxisTitle("PV track (p_{T} > 10 GeV) #eta", 1);
  dzVsEta_pt10->setAxisTitle("PV track (p_{T} > 10 GeV) d_{z} (#mum)", 2);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void VertexMonitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
