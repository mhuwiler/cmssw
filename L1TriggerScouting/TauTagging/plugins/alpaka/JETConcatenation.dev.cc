// Check that ALPAKA_HOST_ONLY is not defined during device compilation:
#ifdef ALPAKA_HOST_ONLY
#error ALPAKA_HOST_ONLY defined in device compilation
#endif

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "L1TriggerScouting/TauTagging/interface/alpaka/JETConcatenation.h"
//#include "L1TriggerScouting/JetClusteringTagging/interface/alpaka/Utils.h"
//#include "L1TriggerScouting/JetClusteringTagging/interface/alpaka/Clustering.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/TauClusterCollection.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

using namespace cms::alpakatools;


// Insertion sorting
    ALPAKA_FN_ACC void insertionSort(float* data, int* spectator, int N)
    {
      for (int i = 1; i < N; ++i) 
      {
        float key = data[i];
        int temp = spectator[i];
        int j = i - 1;
        while (j >= 0 && data[j] < key) 
        {
          data[j + 1] = data[j];
          spectator[j + 1] = spectator[j]; 
          --j;
        }
        data[j + 1] = key;
        spectator[j + 1] = temp; 
      }
    }


    ALPAKA_FN_ACC void to4vector(float pt, float eta, float phi, float fmass) 
    {
     
    }

    ALPAKA_FN_ACC void jet4vector(float* pt, float* eta, float* phi, float* mass, int N, float& jetPt, float& jetEta, float& jetPhi, float& jetE) 
    {
      float x = 0; 
      float y = 0; 
      float z = 0; 
      float E = 0; 
      
      for (int i=0; i<N; i++) 
      {
        float currentPt = pt[i]; 

        float px = currentPt * cos(phi[i]); 

        float py = currentPt * sin(phi[i]); 

        float pz = currentPt * sinh(eta[i]); 

        x += px; 
        y += py; 
        z += pz; 

        E += sqrt(px*px + py*py + pz*pz + mass[i]*mass[i]); 

      }

      jetPt = sqrt(x*x + y*y); 
      jetEta = 0.; 
      if (jetPt > 0.) jetEta = asinh(z/jetPt); // 0.5 * log((p + pz_sum) / (p - pz_sum))
      jetPhi = atan2(y, x); 


    }


class JETConcatenationKernel 
{
public: 
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc> > >
    ALPAKA_FN_ACC void operator()(const TAcc& acc, PFCandidateCollection::ConstView pf, CLUEsteringCollection::ConstView clusters, TauClusterCollection::View jets) const 
    {
      printf("Starting kernel\n"); 

      //using Dim = alpaka::Dim<TAcc>;
        //using Idx = alpaka::Idx<TAcc>;
      using Vec = alpaka::Vec<alpaka::Dim<TAcc>, alpaka::Idx<TAcc> >;
      using Vec1D = alpaka::Vec<alpaka::DimInt<1u>, alpaka::Idx<TAcc> >;

      Vec const globalThreadIdx = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc);
        Vec const globalThreadExtent = alpaka::getWorkDiv<alpaka::Grid, alpaka::Threads>(acc);

        // Map the three dimensional thread index into a
        // one dimensional thread index space. We call it
        // linearize the thread index.
        uint32_t const idx = alpaka::mapIdx<1u>(globalThreadIdx, globalThreadExtent).front();


        // First get the number of clusters
        // TODO: have it produced as metadata by the clustering
        uint32_t numJets = 0; 

        for (int32_t i=0; i<pf.metadata().size(); i++) 
        {

          // Accessing the column "cluster" by its name as a finctional
          auto a = clusters.cluster()[i]; 

          printf("Cluster number %u\n", a); 

          if (a > numJets) numJets = a; 
        }

        numJets +=1; // numbering of clusters starts at 0



        // Loop over the PF collection, and extract the candidates with cluster number matching the thread number
        // TODO: make sure only nJets threads operate

        if (idx > numJets) return; 


        const int N = 128; 
        const int P = 16; 

        // Extracting the indices and the corresponding pt values from the PF candidates matching the cluster number
        int indices[N] = {-1}; 
        int *ind = indices; 

        float pt[N] = {-999}; 
        float *currentpt = pt; 

        for (int i = 0; i < pf.metadata().size(); i++) // make this copy in parallel for GPU
        {
          printf("i value %d, idx: %d\n", i, idx); 
          //if (clusters.cluster()[i] == idx) //copy the PF locally
          {
            *ind = i; 
            ind++; 
            *currentpt = pf.pt()[i]; 
            currentpt++; 
          }
        }


        for (uint32_t i = 0; i < N; i++) 
        {
          printf("Pt value %f, index: %d\n", pt[i], indices[i]); 
        }

        // Sorting the arrays according to pT 
        insertionSort(pt, indices, N); 
        printf("Sorted:\n"); 
        
        for (uint32_t i = 0; i < N; i++) 
        {
          printf("Pt value %f, index: %d\n", pt[i], indices[i]); 
        }

        jets.pt()[0] = 1.; 

        std::copy(pt, pt+P, &jets.pt()[0]); 

        printf("Output collection:\n"); 
        for (uint32_t i = 0; i < P; i++) 
        {
          printf("Pt value %f\n", jets.pt()[i]); 
        }


        // Filling the other data columns
        float masses[P]; 
        for (int i=0; i<P; i++) 
        {
          jets.deltaeta()[i] = pf.eta()[indices[i]]; 
          jets.deltaphi()[i] = pf.phi()[indices[i]]; 
          jets.vz()[i] = pf.z0()[indices[i]]; 
          //jets.dxy()[i] = pf.dxy()[indices[i]]; 
          //jets.pdgid()[i] = pf.pdgid()[indices[i]]; 
          masses[i] = 0.139; //pf.mass()[i]; 
        }

        // Compute the jet 4 momentum
        float jetpt, jeteta, jetphi, jetm = 0; 
        jet4vector(jets.pt(), jets.deltaeta(), jets.deltaphi(), masses, P, jetpt, jeteta, jetphi, jetm); 

        printf("Pt value %f, eta: %f, phi: %f, m: %f\n", jetpt, jeteta, jetphi, jetm); 


        for (int i=0; i<P; i++) {
          jets.deltaeta()[i] -= jeteta; 
        }




    }

};


// Function to launch the kernel
//template <typename TAcc>
void Concatenate(Queue& queue, const PFCandidateCollection& pf, const CLUEsteringCollection& clusters, TauClusterCollection& jets, const int Nclusters) 
{
  uint32_t threads_per_block = Nclusters;
  uint32_t blocks_per_grid = 1;        
  auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);
  //alpaka::exec<Acc1D>(queue, grid, JETConcatenationKernel{}, pf.const_view(), clusters.const_view()); //, Nclusters
  alpaka::exec<Acc1D>(queue, grid, JETConcatenationKernel{}, pf.const_view(), clusters.const_view(), jets.view()); //, Nclusters
  //alpaka::wait(queue);
}


}  // namespace ALPAKA_ACCELERATOR_NAMESPACE