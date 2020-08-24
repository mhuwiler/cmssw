#include <cuda_runtime.h>
#include "CUDADataFormats/Common/interface/Product.h"
#include "CUDADataFormats/Common/interface/HostProduct.h"
#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "HeterogeneousCore/CUDACore/interface/ScopedContext.h"
#include "RecoLocalTracker/SiPixelRecHits/interface/pixelCPEforGPU.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

class SiPixelValidateVerticesFromSoA : public DQMEDAnalyzer {

public:
  explicit SiPixelValidateVerticesFromSoA(const edm::ParameterSet& iConfig);
  ~SiPixelValidateVerticesFromSoA() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void analyze(const edm::Event& iEvent, edm::EventSetup const& iSetup) override;
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;
  
  edm::EDGetTokenT<ZVertexHeterogeneous> tokenvertexsoa_;
  edm::EDGetTokenT<reco::BeamSpot> tokenbeamSpot_;
  edm::EDGetTokenT<PixelTrackHeterogeneous> tokenTracks_;  

  std::string topFolderName_;
  
  MonitorElement* hnVertices;
  MonitorElement* hvx;
  MonitorElement* hvy;
  MonitorElement* hvz;
  MonitorElement* ndof; 
  MonitorElement* chi2; 
  MonitorElement* pt2; 
  MonitorElement* wv; 
  MonitorElement* chi2ndof; 
  MonitorElement* tracksumPt2; 
  MonitorElement* ptvsNtracks; 
  MonitorElement* ndofvsNtracks; 

  MonitorElement* ntracks; 
  MonitorElement* trackQuality; 
  MonitorElement* associatedTrack; 
  MonitorElement* trackQualityPass; 
  MonitorElement* trackPt; 
  MonitorElement* trackEta; 
  MonitorElement* trackPhi; 
  MonitorElement* trackCharge; 
  MonitorElement* trackZip; 
  MonitorElement* trackTip; 


  // LWF like plots 
  MonitorElement* hvr; 
  MonitorElement* ntrackslwf; 
  MonitorElement* nverticeslwf; 
  MonitorElement* zposlwf; 
  MonitorElement* xposlwf; 
  MonitorElement* yposlwf; 
  MonitorElement* chi2lwf; 
  MonitorElement* ndoflwf; 
 
};

SiPixelValidateVerticesFromSoA::SiPixelValidateVerticesFromSoA(const edm::ParameterSet& iConfig)
    : tokenvertexsoa_(consumes<ZVertexHeterogeneous>(iConfig.getParameter<edm::InputTag>("src"))),
      tokenbeamSpot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotsrc"))),
      //tokenIndToEdm_(consumes<IndToEdm>(conf.getParameter<edm::InputTag>("TrackCollectionsrc"))), 
      tokenTracks_(consumes<PixelTrackHeterogeneous>(iConfig.getParameter<edm::InputTag>("trackCollectionsrc")))
{
  topFolderName_ = "SiPixelHeterogeneousV/PixelVerticesSoA";
}

void SiPixelValidateVerticesFromSoA::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& run , edm::EventSetup const& es) {
  ibooker.cd();
  ibooker.setCurrentFolder(topFolderName_);
  hnVertices = ibooker.book1D("nVerticesgpu", ";# of vertices per event;", 150, -0.5, 149.5);
  hvx = ibooker.book1D("vertex_xpos", ";vertex x position;", 10, -0.5, 0.5);
  hvy = ibooker.book1D("vertex_ypos", ";vertex y position;", 10, -0.5, 0.5);
  hvz = ibooker.book1D("vertex_zpos", ";vertex z position;", 10, -20., 20.);
  ndof = ibooker.book1D("ndof", ";number of degrees of freedom per vertex (ndof);", 171, -0.5, 170.5);
  chi2 = ibooker.book1D("chi2", ";chi^{2} of vertex;", 151, -0.5, 150.5); 
  std::vector<Float_t> pt2bins = {0., 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200., 500., 1000., 2000., 5000., 10000.};  // Set to match the LWF binning
  pt2 = ibooker.book1D("pt2", ";p_{T}^{2} of tracks per vertex;", pt2bins.size()-1, pt2bins.data());   // 450, -0.5, 4499.5
  wv = ibooker.book1D("wv", ";Error^{2} on z position of vertex (1/weight);", 100, -0.5, 1.5); 
  ntracks = ibooker.book1D("ntracks", ";# of tracks per vertex;", 101, -0.5, 199.5); 
  std::vector<Float_t> ntracksbins = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 22.0, 26.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0, 200.0}; 
  //ntrackslwf = ibooker.book1D("ntrackslwf", ";# of tracks per vertex (legacy binning);", ntracksbins.size()-1, ntracksbins.data()); 
  trackQuality = ibooker.book1D("trackQuality", ";Quality of tracks associated with vertices;", 6, -0.5, 5.5); 
  associatedTrack = ibooker.book1D("associated tracks", ";track associated;", 100, -1.5, 98.5); 
  trackQualityPass = ibooker.book1D("trackQualityPass", ";Quality > loose of tracks forming a vertex;", 6, -0.5, 5.5); 
  trackPt = ibooker.book1D("trackPt", ";p_{T} of tracks associated with vertices;", 700, -0.5, 700.5); 
  trackCharge = ibooker.book1D("trackCharge", ";Charge of the tracks associated with vertices;", 5, -2.5, 2.5); 
  trackZip = ibooker.book1D("trackZip", ";z IP of tracks associated with vertices;", 19, -9.5, 9.5); 
  trackTip = ibooker.book1D("trackTip", ";T IP of tracks associated with vertices;", 21, -10.5, 10.5); 
  trackEta = ibooker.book1D("trackEta", ";#eta of tracks associated with vertices;", 11, -5.5, 5.5);
  trackPhi = ibooker.book1D("trackPhi", ";#phi of tracks associated with vertices;", 7, -3.5, 3.5);  
  chi2ndof = ibooker.book1D("chi2ndof", ";#chi^{2}/ndof per vertex;", 100, -0.5, 99.5); 
  tracksumPt2 = ibooker.book1D("tracksumPt2", ";Sum of p_{T} of the tracks per vertex;", 450, -0.5, 4499.5); 
  ptvsNtracks = ibooker.book2D("ptvsNtracks", "p_{T}^{2} as function of the number of tracks;number of tracks per vertex;p_{T}^{2} of vertex", 101, -0.5, 199.5, 1500, -0.5, 1499.5); 
  ndofvsNtracks = ibooker.book2D("ndofvsNtracks", "ndof as function of the number of tracks;number of tracks;number of degrees of freedom (ndof)", 101, -0.5, 199.5, 171, -0.5, 170.5); 

  hvr = ibooker.book1D("vertex_rpos", ";vertex R position;", 400, -0.5, 1.5); 

  // LWF histos
  nverticeslwf = ibooker.book1D("nverticeslwf", ";# of vertices per event;", 80, -0.5, 79.5); //nverticeslwf = ibooker.book1D("nverticeslwf", ";# of vertices per event;", 100, 0., 200.); 
  zposlwf = ibooker.book1D("zposlwf", ";vertex z position;", 100, -20., 20.); //  zpos = ibooker.book1D("zpos", ";vertex z position;", 120, -60., 60.); 
  xposlwf = ibooker.book1D("xposlwf", ";vertex x position;", 100, 0., 0.2); //  xpos = ibooker.book1D("xpos", ";vertex x position;", 120, -0.6, 0.6); 
  yposlwf = ibooker.book1D("yposlwf", ";vertex y position;", 120, -0.1, 0.1); //  ypos = ibooker.book1D("ypos", ";vertex y position;", 120, -0.6, 0.6); 
  ntrackslwf = ibooker.book1D("ntrackslwf", ";# number of tracks per vertex;", 40, -0.5, 99.5); 
  chi2lwf = ibooker.book1D("chi2lwf", ";chi^{2} of vertex;", 100, 0., 200.); 
  ndoflwf = ibooker.book1D("ndoflwf", ";number of degrees of freedom per vertex (ndof)", 100, 0., 200.); 

  


}

void SiPixelValidateVerticesFromSoA::analyze(const edm::Event& iEvent, edm::EventSetup const& es) {
  
  const auto& vertexsoa = *(iEvent.get(tokenvertexsoa_).get());

  assert(&vertexsoa); 

  edm::ESHandle<GeometricDet> geomDetHandle; //taken from: https://github.com/cms-sw/cmssw/blob/master/DPGAnalysis/SiStripTools/plugins/DetIdSelectorTest.cc
  iSetup.get<IdealGeometryRecord>().get(geomDetHandle);
  const auto detids = TrackerGeometryUtils::getSiStripDetIds(*geomDetHandle);

  auto nv = vertexsoa.nvFinal;
  if(nv > vertexsoa.MAXVTX)   return;
  hnVertices->Fill(nv);
  nverticeslwf->Fill(nv); 

  const auto& tracksoa = *(iEvent.get(tokenTracks_)); 
  const auto *quality = tracksoa.qualityData(); 

  assert(&tracksoa); 

  uint32_t maxNumTracks = tracksoa.stride(); // This is the dimension of the soas, not the actual number of tracks. We iterate over them and break the loop if nHits is 0 (index no longer used). 
  

  //beamSpot
  edm::Handle<reco::BeamSpot> beampspotHandle;
  iEvent.getByToken(tokenbeamSpot_, beampspotHandle);
  float x0 = 0., y0 = 0., z0 = 0., dxdz = 0., dydz = 0.;
  if (!beampspotHandle.isValid()) {
    edm::LogWarning("SiPixelValidateVerticesFromSoA") << "No beamspot found. returning vertexes with (0,0,Z) ";
  } else {
    const reco::BeamSpot &bs = *beampspotHandle;
    x0 = bs.x0();
    y0 = bs.y0();
    z0 = bs.z0();
    dxdz = bs.dxdz();
    dydz = bs.dydz();
  }

  for (unsigned int i = 0; i < nv; i++) {
    auto z = vertexsoa.zv[i];
    auto x = x0 + dxdz * z;
    auto y = y0 + dydz * z;
    z += z0;
    auto r = sqrt(x*x+y*y); 
    hvx->Fill(x);
    hvy->Fill(y);
    hvz->Fill(z);
    hvr->Fill(r); 
    ndof->Fill(vertexsoa.ndof[i]); 
    chi2->Fill(vertexsoa.chi2[i]); 
    pt2->Fill(vertexsoa.ptv2[i]); 
    wv->Fill(1./static_cast<Double_t>(vertexsoa.wv[i])); 
    zposlwf->Fill(z); 
    xposlwf->Fill(x); 
    yposlwf->Fill(y); 
    chi2lwf->Fill(vertexsoa.chi2[i]); 
    ndoflwf->Fill(vertexsoa.ndof[i]); 
    if (vertexsoa.ndof[i]>0) 
    {
      chi2ndof->Fill(static_cast<Double_t>(vertexsoa.chi2[i])/static_cast<Double_t>(vertexsoa.ndof[i])); 
    }


    Double_t sumPt2 = 0; 

    // Start accesing the tracks 
    uint32_t nTk = 0; // Total number of tracks assiciated with the vertex. 
    uint32_t nTkLoose = 0; // Number of tracks with quality 'loose'. 
    uint32_t nTkGood = 0; // Number of good tracks (with quality higher than qualityThreshold)
    auto tkQualityThres = trackQuality::loose; // Quality above which (included) the track is considered of good quality.   

    for (uint32_t iTk=0; iTk<maxNumTracks; iTk++) // put here vertexsoa.MAXTRACKS ? 
    {
      auto nHits = tracksoa.nHits(iTk); 
      if (nHits==0) break; // Since we are looping over the size of the soa, we need to escape at the point where the elements are no longer used. 
      auto qual = quality[iTk];   
      

      if (vertexsoa.idv[iTk] == vertexsoa.sortInd[i]) // The track is associated with this vertex 
      {
        
        
        // Filling histograms
        associatedTrack->Fill(vertexsoa.idv[iTk]); 
        trackQuality->Fill(qual);   
        trackPt->Fill(tracksoa.pt(iTk)); 
        trackEta->Fill(tracksoa.eta(iTk)); 
        trackPhi->Fill(tracksoa.phi(iTk)); 
        trackCharge->Fill(tracksoa.charge(iTk)); 
        trackZip->Fill(tracksoa.zip(iTk));
        trackTip->Fill(tracksoa.tip(iTk)); 

        // Making stats 
        nTk++;    

        sumPt2 += tracksoa.pt(iTk)*tracksoa.pt(iTk); 

        if (qual >= tkQualityThres) 
        {
          nTkGood++; 
          trackQualityPass->Fill(qual); 
        }        

        if (qual == trackQuality::loose)  
        {
          nTkLoose++;
        }  
        
        std::cout << "Track: " << iTk  << " is associated to vetex: " << i << ", out of: " << nv << " matching" << std::endl;  
      }
      else 
      {
        //std::cout << "Track: " << iTk  << " is not associated to vetex: " << i << ", out of: " << nv << "not matching" << std::endl; 
      }

    }  

    std::cout << "Num tracks: " << nTk << ", (out of max number: " << maxNumTracks << "), thereof good tracks: " << nTkGood 
      << ", N loose tracks: " << nTkLoose 
      << std::endl;   

    ntracks->Fill(nTk);  
    ntrackslwf->Fill(nTk); 

    tracksumPt2->Fill(sumPt2);  

    ptvsNtracks->Fill(nTk, vertexsoa.ptv2[i]); 

    ndofvsNtracks->Fill(nTk, vertexsoa.ndof[i]); 


  }

  //if (maxNumTracks > vertexsoa.MAXTRACKS) 
  //{
  //  std::cerr << "The maximum number of tracks exceeds the maximum number of tracks in the vertex. " << std::endl; 
  //  return;
  //} 
}
  

void SiPixelValidateVerticesFromSoA::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("trackCollectionsrc", edm::InputTag("pixelTracks"));
  desc.add<edm::InputTag>("beamSpotsrc", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("src", edm::InputTag("pixelVertexSoA"));
  descriptions.add("SiPixelValidateVerticesFromSoA", desc);
}
DEFINE_FWK_MODULE(SiPixelValidateVerticesFromSoA);
