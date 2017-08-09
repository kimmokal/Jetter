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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

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

//#include "Jetter/MiniAnalyzer/plugins/MiniAnalyzer.h"

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
    public:
        explicit MiniAnalyzer(const edm::ParameterSet&);
        ~MiniAnalyzer();

    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<pat::JetCollection> jetToken_;
        
        TFile* outputFile;
        TTree* tree;
        
        //defining the jet specific parameters
        float pT;
        float eta;
        float phi;
        float E;
        unsigned int event;
        unsigned int run;
        unsigned int lumi;

};


MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
    
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
{
    //sets the output file, tree and the parameters to be saved to the tree

    outputFile = new TFile("nanotuple.root","recreate");
    jetTree = new TTree("jetTree", "A simplified jet tree");

    jetTree->Branch("pT", &pT, "pT/F");
    jetTree->Branch("eta", &eta, "eta/F");
    jetTree->Branch("phi", &phi, "phi/F");
    jetTree->Branch("mass", &mass, "mass/F");
    jetTree->Branch("et", &et, "et/F");
	
    jetTree->Branch("event", &event, "event/l");
    jetTree->Branch("run", &run, "run/l");
    jetTree->Branch("lumi", &lumi, "lumi/l");

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
    
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    for (const pat::Jet &j : *jets) {
        if (j.pt() < 20) continue;
        
        //adding jet parameters to jet-based tree
        pT = j.pt();
        eta = j.eta();
        phi = j.phi();
        mass = j.mass();
        et = j.et();

        //adding event information to jet-based tree
        event = iEvent.id().event();
        run = iEvent.id().run();
        lumi = iEvent.id().lumi();


        jetTree->Fill();
    }


}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
