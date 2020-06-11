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
  hnVertices = ibooker.book1D("nVerticesgpu", ";#of vertices;", 150, -0.5, 149.5);
  hvx = ibooker.book1D("vertex_xpos", ";xPos;", 100, -0.5, 0.5);
  hvy = ibooker.book1D("vertex_ypos", ";yPos;", 10, -0.5, 0.5);
  hvz = ibooker.book1D("vertex_zpos", ";zPos;", 10, -20., 20.);
  ndof = ibooker.book1D("ndof", ";ndof;", 171, -0.5, 170.5);
  chi2 = ibooker.book1D("chi2", ";chi2;", 100, -0.5, 1.5); 
  pt2 = ibooker.book1D("pt2", ";p_{T}^{2} of tracks forming vertex;", 450, -0.5, 4499.5); 
  wv = ibooker.book1D("wv", ";Output weight (1/error^2 on z position) of vertex;", 100, -0.5, 1.5); 
  ntracks = ibooker.book1D("ntracks", ";#of tracks;", 101, -0.5, 199.5); 
  trackQuality = ibooker.book1D("trackQuality", ";Track quality;", 6, -0.5, 5.5); 
  associatedTrack = ibooker.book1D("associated tracks", ";track associated;", 100, -1.5, 98.5); 
  trackQualityPass = ibooker.book1D("trackQualityPass", ";Quality of tracks passing the threshold;", 6, -0.5, 5.5); 
  trackPt = ibooker.book1D("trackPt", ";p_T of tracks associated with vertices;", 700, -0.5, 700.5); 
  trackCharge = ibooker.book1D("trackCharge", ";Charge of the track;", 5, -2.5, 2.5); 
  trackZip = ibooker.book1D("trackZip", ";ZIP of tracks;", 19, -9.5, 9.5); 
  trackTip = ibooker.book1D("trackTip", ";TIP of tracks;", 21, -10.5, 10.5); 
  trackEta = ibooker.book1D("trackEta", ";#eta of the tracks;", 11, -5.5, 5.5);
  trackPhi = ibooker.book1D("trackPhi", ";#phi of the tracks;", 7, -3.5, 3.5);  
  chi2ndof = ibooker.book1D("chi2ndof", ";#chi2/ndof per vertex;", 100, -0.5, 99.5); 
  tracksumPt2 = ibooker.book1D("tracksumPt2", ";Sum of p_{T} per vertex;", 450, -0.5, 4499.5); 
  ptvsNtracks = ibooker.book2D("ptvsNtracks", "Pt^2 as function of the number of tracks;# of tracks per vertex;p_{T}^{2}", 101, -0.5, 199.5, 1500, -0.5, 1499.5); 


}

void SiPixelValidateVerticesFromSoA::analyze(const edm::Event& iEvent, edm::EventSetup const& es) {
  
  const auto& vertexsoa = *(iEvent.get(tokenvertexsoa_).get());

  assert(&vertexsoa); 

  auto nv = vertexsoa.nvFinal;
  if(nv > vertexsoa.MAXVTX)   return;
  hnVertices->Fill(nv);

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
    hvx->Fill(x);
    hvy->Fill(y);
    hvz->Fill(z);
    ndof->Fill(vertexsoa.ndof[i]); 
    chi2->Fill(vertexsoa.chi2[i]); 
    pt2->Fill(vertexsoa.ptv2[i]); 
    wv->Fill(1./static_cast<Double_t>(vertexsoa.wv[i])); 
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

    tracksumPt2->Fill(sumPt2);  

    ptvsNtracks->Fill(nTk, vertexsoa.ptv2[i]); 


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
