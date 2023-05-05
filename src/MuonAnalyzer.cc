// -*- C++ -*-
//
// Package:    GenParticle/GenParticleAnalyzer
// Class:      GenParticleAnalyzer
//
 
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//class to extract MC information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//classes to save data
#include "TTree.h"
#include "TFile.h"
#include<vector>


//class to extract Muon information
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

//TransientTrack and IPTools for impact parameter
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/Math/interface/deltaR.h"


// All of these libraries above are the standars for muons and gen particles. The unique exception is deltaR.


class MuonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonAnalyzer(const edm::ParameterSet&);
      ~MuonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;  
      std::vector<std::string>  particle;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;

      // ----------member data ---------------------------
      TTree *mtree;      
      //muons 
      int nummuon; 
      std::vector<float> muon_e;
      std::vector<float> muon_pt;
      std::vector<float> muon_px;
      std::vector<float> muon_py;
      std::vector<float> muon_pz;
      std::vector<float> muon_eta;
      std::vector<float> muon_phi;
      std::vector<float> muon_ch;
      std::vector<float> muon_genpartidx;

      
     //gens
      
      int numGenPart;
      std::vector<int> GenPart_status;
      std::vector<float> GenPart_pt;
      std::vector<float> GenPart_eta;
      std::vector<float> GenPart_mass;
      std::vector<int> GenPart_pdgId;
      std::vector<float> GenPart_phi;
      std::vector<float> GenPart_px;
      std::vector<float> GenPart_py;
      std::vector<float> GenPart_pz;
      std::vector<int> GenPart_motherId;
      
};


MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig):
 muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
 vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
 particle(iConfig.getParameter<std::vector<std::string> >("input_particle")), 
 prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned")))
{
   edm::Service<TFileService> fs;
   mtree = fs->make<TTree>("Events", "Events");
   //gens 
   
   mtree->Branch("numGenPart",&numGenPart);
   mtree->GetBranch("numGenPart")->SetTitle("number of generator particles");
   mtree->Branch("GenPart_pt",&GenPart_pt);
   mtree->GetBranch("GenPart_pt")->SetTitle("generator particle transverse momentum");
   mtree->Branch("GenPart_eta",&GenPart_eta);
   mtree->GetBranch("GenPart_eta")->SetTitle("generator particle pseudorapidity");
   mtree->Branch("GenPart_mass",&GenPart_mass);
   mtree->GetBranch("GenPart_mass")->SetTitle("generator particle mass");
   mtree->Branch("GenPart_pdgId",&GenPart_pdgId);
   mtree->GetBranch("GenPart_pdgId")->SetTitle("generator particle PDG id");
   mtree->Branch("GenPart_phi",&GenPart_phi);
   mtree->GetBranch("GenPart_phi")->SetTitle("generator particle azimuthal angle of momentum vector");
   mtree->Branch("GenPart_px",&GenPart_px);
   mtree->GetBranch("GenPart_px")->SetTitle("generator particle x coordinate of momentum vector");
   mtree->Branch("GenPart_py",&GenPart_py);
   mtree->GetBranch("GenPart_py")->SetTitle("generator particle y coordinate of momentum vector");
   mtree->Branch("GenPart_pz",&GenPart_pz);
   mtree->GetBranch("GenPart_pz")->SetTitle("generator particle z coordinate of momentum vector");
   mtree->Branch("GenPart_status",&GenPart_status);
   mtree->GetBranch("GenPart_status")->SetTitle("Particle status. 1=stable");
   mtree->Branch("GenPart_motherId", &GenPart_motherId);
   mtree->GetBranch("GenPart_motherId")->SetTitle("MotherID");
   
   
   //muons
   
   mtree->Branch("numbermuon",&nummuon);
   mtree->GetBranch("numbermuon")->SetTitle("number of muons");
   mtree->Branch("muon_e",&muon_e);
   mtree->GetBranch("muon_e")->SetTitle("muon energy");
   mtree->Branch("muon_pt",&muon_pt);
   mtree->GetBranch("muon_pt")->SetTitle("muon transverse momentum");
   mtree->Branch("muon_px",&muon_px);
   mtree->GetBranch("muon_px")->SetTitle("muon momentum x-component");
   mtree->Branch("muon_py",&muon_py);
   mtree->GetBranch("muon_py")->SetTitle("muon momentum y-component");
   mtree->Branch("muon_pz",&muon_pz);
   mtree->GetBranch("muon_pz")->SetTitle("muon momentum z-component");
   mtree->Branch("muon_eta",&muon_eta);
   mtree->GetBranch("muon_eta")->SetTitle("muon pseudorapidity");
   mtree->Branch("muon_phi",&muon_phi);
   mtree->GetBranch("muon_phi")->SetTitle("muon polar angle");
   mtree->Branch("muon_ch",&muon_ch);
   mtree->GetBranch("muon_ch")->SetTitle("muon charge");
   mtree->Branch("muon_genpartidx",&muon_genpartidx);
   mtree->GetBranch("muon_genpartidx")->SetTitle("index into genParticle list for MC matching to status==1 muons");
   
}


MuonAnalyzer::~MuonAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void
MuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   
   //muons
   Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);
   Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   

   nummuon = 0;
   muon_e.clear();
   muon_pt.clear();
   muon_px.clear();
   muon_py.clear();
   muon_pz.clear();
   muon_eta.clear();
   muon_phi.clear();
   muon_ch.clear();
  // muon_genpartidx.clear(); 
  
    for (const pat::Muon &mu : *muons)
    {
      muon_e.push_back(mu.energy());
      muon_pt.push_back(mu.pt());
      muon_px.push_back(mu.px());
      muon_py.push_back(mu.py());
      muon_pz.push_back(mu.pz());
      muon_eta.push_back(mu.eta());
      muon_phi.push_back(mu.phi());
      muon_ch.push_back(mu.charge());
      //muon_genpartidx.push_back(-1);
      
      nummuon++;
    }   
    

   //gens
   Handle<edm::View<reco::GenParticle> > pruned;
   iEvent.getByToken(prunedGenToken_,pruned);

   numGenPart=0;
   GenPart_pt.clear();
   GenPart_eta.clear();
   GenPart_mass.clear();
   GenPart_pdgId.clear();
   GenPart_phi.clear();
   GenPart_px.clear();
   GenPart_py.clear();
   GenPart_pz.clear();
   GenPart_status.clear();
   GenPart_motherId.clear();
    
 std::vector<reco::GenParticle> interestingGenParticles;
 for (auto it = pruned->begin(); it != pruned->end(); it++) {
   const auto status = it->status();
   const auto pdgId = std::abs(it->pdgId());
   if (status == 1 && pdgId == 13) { 
     interestingGenParticles.emplace_back(*it);
   }
 }
 
   for (auto mu = muons->begin(); mu!= muons->end(); mu++) {
    float minDeltaR = 999.0;
    int idx = -1;
    for (auto igp = interestingGenParticles.begin(); igp!= interestingGenParticles.end();igp++) {
	const auto tmp = deltaR(igp->p4(), mu->p4());
	if (tmp < minDeltaR) {
	    minDeltaR = tmp;
	    idx= igp - interestingGenParticles.begin();
	}
    }
    if (idx!=-1){
	    auto g = interestingGenParticles.begin() + idx ;
	    GenPart_motherId.push_back(g->mother()->pdgId());
	    GenPart_pt.push_back(g->pt());
	    GenPart_eta.push_back(g->eta());
	    GenPart_mass.push_back(g->mass());
	    GenPart_pdgId.push_back(g->pdgId());
	    GenPart_phi.push_back(g->phi());
	    GenPart_status.push_back(g->status());
	    GenPart_px.push_back(g->px());
	    GenPart_py.push_back(g->py());
	    GenPart_pz.push_back(g->pz());
    
    	numGenPart++;
    }
  }   
  

  mtree->Fill();
  return; 
         
}



// ------------ method called once each job just before starting event loop  ------------
void
MuonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MuonAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonAnalyzer);

