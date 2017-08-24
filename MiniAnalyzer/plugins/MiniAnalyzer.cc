// -*- C++ -*-
//
// Package:    Jetter/MiniAnalyzer
// Class:      MiniAnalyzer
// 
//**\class MiniAnalyzer MiniAnalyzer.cc Jetter/MiniAnalyzer/plugins/MiniAnalyzer.cc
//
// Original Author:  Petra-Maria Ekroos
//         Created:  Wed, 02 Aug 2017 14:50:21 GMT
//
//
// system include files
#include <memory>

#include <vector>
#include <cmath>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"


//MiniAOD PAT libraries
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


//ROOT libraries
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"


//#include "Jetter/MiniAnalyzer/plugins/MiniAnalyzer.h"

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
    public:
        explicit MiniAnalyzer(const edm::ParameterSet&);
        ~MiniAnalyzer();

        struct PFV {float pT,dR,dTheta, mass;};
        static PFV pfv;


    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        edm::EDGetTokenT<pat::MuonCollection> muonToken_;
        edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
        edm::EDGetTokenT<pat::TauCollection> tauToken_;
        edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
        edm::EDGetTokenT<pat::JetCollection> jetToken_;
        edm::EDGetTokenT<pat::METCollection> metToken_;
        edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
       
        TFile* outputFile;
        TTree* jetTree;
        
        //defining the jet specific parameters
        float jetPt;
        float jetEta;
        float jetPhi;
        float mass;
        float et;

        unsigned int event;
        unsigned int run;
        unsigned int lumi;
        //unsigned int bx;
        
        //typedef struct {Float_t pT,deltaR,deltaTheta,mass,type;} PF;

};

MiniAnalyzer::PFV MiniAnalyzer::pfv;

MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
    
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))

{
    //sets the output file, tree and the parameters to be saved to the tree

    outputFile = new TFile("nanotuple_particles.root","recreate");
    jetTree = new TTree("jetTree", "A simplified jet tree");

    jetTree->Branch("jetPt", &jetPt, "jetPt/F");
    jetTree->Branch("jetEta", &jetEta, "jetEta/F");
    jetTree->Branch("jetPhi", &jetPhi, "jetPhi/F");
    jetTree->Branch("mass", &mass, "mass/F");
    jetTree->Branch("et", &et, "et/F");
    	
    jetTree->Branch("event", &event, "event/l");
    jetTree->Branch("run", &run, "run/l");
    jetTree->Branch("lumi", &lumi, "lumi/l");
    //jetTree->Branch("bx", &bx, "bx/l");


    jetTree->Branch("pfv", &pfv, "pT:dR:dTheta:mass");
    


}

MiniAnalyzer::~MiniAnalyzer()
{

//at the end, write data into tree
	jetTree->Write();
	outputFile->Close();

}


void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    //kaikki partikkelit läpi

    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    edm::Handle<pat::PhotonCollection> photons;
    iEvent.getByToken(photonToken_, photons);
    edm::Handle<pat::TauCollection> taus;
    iEvent.getByToken(tauToken_, taus);
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    edm::Handle<pat::PackedCandidateCollection> pfs;
    iEvent.getByToken(pfToken_, pfs);  


    for (const pat::Jet &j : *jets) {
        if (j.pt() < 20) continue;

        //adding jet parameters to jet-based tree
        jetPt = j.pt();
        jetEta = j.eta();
        jetPhi = j.phi();
        mass = j.mass();
        et = j.et();
        

        //adding event information to jet-based tree
        event = iEvent.id().event();
        run = iEvent.id().run();
        lumi = iEvent.id().luminosityBlock();
        //bx = iEvent.id().bx();

//        for (const pat::PackedCandidate &pf : *pfs) {

        for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
            const pat::PackedCandidate &pf = (*pfs)[i];

            float deltaEta = (j.eta()-pf.eta());
            float deltaPhi = std::abs(j.phi()-pf.phi()); if (deltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
        
            if ( (deltaEta < 0.5) && (deltaPhi < 0.5) ) continue;
                
            pfv.pT = pf.pt();
            pfv.dR = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
            //sqrt((deltaEta*deltaEta)+(deltaPhi*deltaPhi));
            pfv.dTheta = std::abs(j.theta() - pf.theta());
            pfv.mass = pf.mass();
            //type = pf.

            /*
                
                pt = pf.pt();
                eta = pf.eta();
                phi = pf.phi();
                mass = pf.mass();
                et = pf.et();
            */
                        
        }

        jetTree->Fill();
    }




}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
