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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

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
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        edm::EDGetTokenT<pat::JetCollection> jetToken_;
        
        TFile* outputFile;
        TTree* tree;

        float pT;
        float eta;
        float phi;
        float E;

        unsigned int event;
        unsigned int run;
        unsigned int lumi;

};


MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets")))
{

//tanne branch maarittelyt
    outputFile = new TFile("nanotuple.root","recreate");
    tree = new TTree("jet_tree", "A simplified jet tree");

    tree->Branch("pT", &pT, "pT/F");
    tree->Branch("eta", &eta, "eta/F");
    tree->Branch("phi", &phi, "phi/F");
    tree->Branch("E", &E, "E/F");
	
    tree->Branch("event", &event, "event/l");
    tree->Branch("run", &run, "run/l");
    tree->Branch("lumi", &lumi, "lumi/l");

}

MiniAnalyzer::~MiniAnalyzer()
{

//lopettaessa
	tree->Write();
	outputFile->Close();

}


void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    for (const pat::Jet &j : *jets) {
        if (j.pt() < 20) continue;
        pT = j.pt();

        tree->Fill();
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
