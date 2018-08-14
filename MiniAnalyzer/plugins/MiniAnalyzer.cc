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
// Modified by: Kimmo Kallonen
//
// system include files
#include <memory>

#include <vector>
#include <string>
#include <cmath>

#include <iostream>
#include <fstream>

//user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"

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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//parton
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

//hadron-level definition
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

//generator-level event information
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//ROOT libraries
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"



using namespace std;

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
    public:
        explicit MiniAnalyzer(const edm::ParameterSet&);
        ~MiniAnalyzer();

        int ngenv, nPF;
        struct PFV {float pT,dR,dTheta, mass;};
        static const int kMaxPF = 5000;

	// Jet constituent variables
	Float_t jetPF_pT[kMaxPF];
	Float_t jetPF_pTrel[kMaxPF];
	Float_t jetPF_dR[kMaxPF];
	Float_t jetPF_dTheta[kMaxPF];
	Float_t jetPF_mass[kMaxPF];
	int jetPF_id[kMaxPF];
        int jetPF_fromPV[kMaxPF];

	// Gen particle variables
        Float_t genPF_pT[kMaxPF];
        Float_t genPF_dR[kMaxPF];
        Float_t genPF_dTheta[kMaxPF];
        Float_t genPF_mass[kMaxPF];
	int genPF_id[kMaxPF];

	//Jet image variables
	Float_t PF_pT[kMaxPF];
	Float_t PF_dR[kMaxPF];
	Float_t PF_dTheta[kMaxPF];
        Float_t PF_dPhi[kMaxPF];
        Float_t PF_dEta[kMaxPF];
	Float_t PF_mass[kMaxPF];
        int PF_id[kMaxPF];
        int PF_fromPV[kMaxPF];


    private:
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

        // ----------member data ---------------------------
        edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
        edm::EDGetTokenT<pat::MuonCollection> muonToken_;
        edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
        edm::EDGetTokenT<pat::TauCollection> tauToken_;
        edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
        edm::EDGetTokenT<pat::JetCollection> jetToken_;
        edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
        edm::EDGetTokenT<pat::METCollection> metToken_;
        edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
      	edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
	edm::EDGetTokenT<reco::GenJetCollection> EDMGenJetsToken_;
        edm::EDGetTokenT<GenEventInfoProduct> genEventInfoToken_;

        edm::EDGetTokenT<edm::ValueMap<float> > qglToken_;
        edm::EDGetTokenT<edm::ValueMap<float> > ptDToken_;
        edm::EDGetTokenT<edm::ValueMap<float> > axis2Token_;
        edm::EDGetTokenT<edm::ValueMap<int>   > multToken_;

	std::string t_qgtagger;

        TFile* outputFile;
        TTree* jetTree;

        //defining the jet specific parameters
        Float_t jetPt;
        Float_t jetEta;
        Float_t jetPhi;
        Float_t jetMass;
        Float_t jetGirth;

	Float_t jetRawPt;
	Float_t jetRawEta;
	Float_t jetRawPhi;
	Float_t jetRawMass;

	unsigned int jetChargedHadronMult;
	unsigned int jetNeutralHadronMult;
	unsigned int jetChargedMult;
	unsigned int jetNeutralMult;
	unsigned int jetMult;

	unsigned int jetLooseID;
	unsigned int jetTightID;

        unsigned int event;
        unsigned int run;
        unsigned int lumi;

	unsigned int eventJetMult;
	unsigned int jetPtOrder;
        
        Float_t dPhiJetsLO;
        Float_t dEtaJetsLO;
        Float_t alpha;

	unsigned int partonFlav;
	unsigned int hadronFlav;
	unsigned int physFlav;

	unsigned int isPartonUDS;
	unsigned int isPartonG;
	unsigned int isPartonOther;
	unsigned int isPhysUDS;
	unsigned int isPhysG;
	unsigned int isPhysOther;

        Float_t genJetPt;
        Float_t genJetEta;
        Float_t genJetPhi;
        Float_t genJetMass;

        Float_t pthat;
        Float_t eventWeight;

	Float_t mMaxY;
	Float_t mMinGenPt;

	//quark-gluon stuff
	Float_t jetQGl;
	Float_t QG_ptD;
	Float_t QG_axis2;
	unsigned int QG_mult;

};

//MiniAnalyzer::PFV MiniAnalyzer::pfv[kMaxPF];
//MiniAnalyzer::PFV MiniAnalyzer::genv[kMaxPF];

MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):

    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
    pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
    packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
    EDMGenJetsToken_(consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
    genEventInfoToken_(consumes <GenEventInfoProduct> (edm::InputTag(std::string("generator")))),
    qglToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
    ptDToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"))),
    axis2Token_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"))),
    multToken_(consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult")))

{

//    mMaxY = iConfig.getParameter<double>("maxY");
    mMinGenPt = iConfig.getUntrackedParameter<double>("minGenPt", 30);

    //sets the output file, tree and the parameters to be saved to the tree

    outputFile = new TFile("QCD_jettuples_pythia8_PU.root","recreate");
    jetTree = new TTree("jetTree", "Jet tree");
    jetTree->SetMaxTreeSize(1000000000); 	// Set max file size to around 1 gigabyte

    jetTree->Branch("jetPt", &jetPt, "jetPt/F");
    jetTree->Branch("jetEta", &jetEta, "jetEta/F");
    jetTree->Branch("jetPhi", &jetPhi, "jetPhi/F");
    jetTree->Branch("jetMass", &jetMass, "jetMass/F");
    jetTree->Branch("jetGirth", &jetGirth, "jetGirth/F");

    jetTree->Branch("jetRawPt", &jetRawPt, "jetRawPt/F");
    jetTree->Branch("jetRawEta", &jetRawEta, "jetRawEta/F");
    jetTree->Branch("jetRawPhi", &jetRawPhi, "jetRawPhi/F");
    jetTree->Branch("jetRawMass", &jetRawMass, "jetRawMass/F");

    jetTree->Branch("jetChargedHadronMult", &jetChargedHadronMult, "jetChargedHadronMult/I");
    jetTree->Branch("jetNeutralHadronMult", &jetNeutralHadronMult, "jetNeutralHadronMult/I");
    jetTree->Branch("jetChargedMult", &jetChargedMult, "jetChargedMult/I");
    jetTree->Branch("jetNeutralMult", &jetNeutralMult, "jetNeutralMult/I");
    jetTree->Branch("jetMult", &jetMult, "jetMult/I");

    jetTree->Branch("jetPF_pT", &jetPF_pT, "jetPF_pT[jetMult]/F");
    jetTree->Branch("jetPF_pTrel", &jetPF_pTrel, "jetPF_pTrel[jetMult]/F");
    jetTree->Branch("jetPF_dR", &jetPF_dR, "jetPF_dR[jetMult]/F");
    jetTree->Branch("jetPF_dTheta", &jetPF_dTheta, "jetPF_dTheta[jetMult]/F");
    jetTree->Branch("jetPF_mass", &jetPF_mass, "jetPF_mass[jetMult]/F");
    jetTree->Branch("jetPF_id", &jetPF_id, "jetPF_id[jetMult]/I");
    jetTree->Branch("jetPF_fromPV", &jetPF_fromPV, "jetPF_fromPV[jetMult]/I");

    jetTree->Branch("ng",&ngenv,"ng/I");

    jetTree->Branch("genPF_pT", &genPF_pT, "genPF_pT[ng]/F");
    jetTree->Branch("genPF_dR", &genPF_dR, "genPF_dR[ng]/F");
    jetTree->Branch("genPF_dTheta", &genPF_dTheta, "genPF_dTheta[ng]/F");
    jetTree->Branch("genPF_mass", &genPF_mass, "genPF_mass[ng]/F");
    jetTree->Branch("genPF_id", &genPF_id, "genPF_id[ng]/I");

    jetTree->Branch("jetLooseID", &jetLooseID, "jetLooseID/I");
    jetTree->Branch("jetTightID", &jetTightID, "jetTightID/I");

    jetTree->Branch("genJetPt", &genJetPt, "genJetPt/F");
    jetTree->Branch("genJetEta", &genJetEta, "genJetEta/F");
    jetTree->Branch("genJetPhi", &genJetPhi, "genJetPhi/F");
    jetTree->Branch("genJetMass", &genJetMass, "genJetMass/F");

    jetTree->Branch("pthat", &pthat, "pthat/F");
    jetTree->Branch("eventWeight", &eventWeight, "eventWeight/F");

    jetTree->Branch("event", &event, "event/l");
    jetTree->Branch("run", &run, "run/I");
    jetTree->Branch("lumi", &lumi, "lumi/I");

    jetTree->Branch("eventJetMult", &eventJetMult, "eventJetMult/I");
    jetTree->Branch("jetPtOrder", &jetPtOrder, "jetPtOrder/I");

    jetTree->Branch("dPhiJetsLO", &dPhiJetsLO, "dPhiJetsLO/F");
    jetTree->Branch("dEtaJetsLO", &dEtaJetsLO, "dEtaJetsLO/F");
    jetTree->Branch("alpha", &alpha, "alpha/F");

    jetTree->Branch("partonFlav", &partonFlav, "partonFlav/I");
    jetTree->Branch("hadronFlav", &hadronFlav, "hadronFlav/I");
    jetTree->Branch("physFlav", &physFlav, "physFlav/I");

    jetTree->Branch("jetQGl", &jetQGl, "jetQGl/F");
    jetTree->Branch("QG_ptD", &QG_ptD, "QG_ptD/F");
    jetTree->Branch("QG_axis2", &QG_axis2, "QG_axis2/F");
    jetTree->Branch("QG_mult", &QG_mult, "QG_mult/I");

    jetTree->Branch("isPartonUDS", &isPartonUDS, "isPartonUDS/I");
    jetTree->Branch("isPartonG", &isPartonG, "isPartonG/I");
    jetTree->Branch("isPartonOther", &isPartonOther, "isPartonOther/I");
    jetTree->Branch("isPhysUDS", &isPhysUDS, "isPhysUDS/I");
    jetTree->Branch("isPhysG", &isPhysG, "isPhysG/I");
    jetTree->Branch("isPhysOther", &isPhysOther, "isPhysOther/I");

    // Jet image variables
    jetTree->Branch("nPF", &nPF, "nPF/I");
    jetTree->Branch("PF_pT", &PF_pT, "PF_pT[nPF]/F");
    jetTree->Branch("PF_dR", &PF_dR, "PF_dR[nPF]/F");
    jetTree->Branch("PF_dTheta", &PF_dTheta, "PF_dTheta[nPF]/F");
    jetTree->Branch("PF_dPhi", &PF_dPhi, "PF_dPhi[nPF]/F");
    jetTree->Branch("PF_dEta", &PF_dEta, "PF_dEta[nPF]/F");
    jetTree->Branch("PF_mass", &PF_mass, "cPF_mass[nPF]/F");
    jetTree->Branch("PF_id", &PF_id, "PF_id[nPF]/I");
    jetTree->Branch("PF_fromPV", &PF_fromPV, "PF_fromPV[nPF]/I");

}

MiniAnalyzer::~MiniAnalyzer()
{

//at the end, write data into tree (needed only when running locally on lxplus
	outputFile = jetTree->GetCurrentFile();
	outputFile->Write();
	outputFile->Close();

}

// Create jet struct for storing the jet index within the event
struct JetIndexed {
	pat::Jet jet;
	unsigned int eventIndex;
	JetIndexed(pat::Jet j, unsigned int eIdx) : jet(j), eventIndex(eIdx) {}
};

// Create a sort function to compare the jet pTs for later pT-ordering
struct higher_pT_sort
{
	inline bool operator() (const JetIndexed& jet1, const JetIndexed& jet2)
	{
		return ( jet1.jet.pt() > jet2.jet.pt() );
	}
};


void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace reco;
    using namespace pat;

    //vector<LorentzVector> mGenJets;

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
    edm::Handle<pat::JetCollection> fatjets;
    iEvent.getByToken(fatjetToken_, fatjets);
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    edm::Handle<pat::PackedCandidateCollection> pfs;
    iEvent.getByToken(pfToken_, pfs);
        // Packed particles are all the status 1, so usable to remake jets
        // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
    edm::Handle<edm::View<pat::PackedGenParticle> > packed;
    iEvent.getByToken(packedGenToken_, packed);
    edm::Handle<GenEventInfoProduct> genEventInfo;
    iEvent.getByToken(genEventInfoToken_, genEventInfo);

    edm::Handle<edm::ValueMap<float>> qglHandle;
    iEvent.getByToken(qglToken_, qglHandle);
    edm::Handle<edm::ValueMap<float>> ptDHandle;
    iEvent.getByToken(ptDToken_, ptDHandle);
    edm::Handle<edm::ValueMap<float>> axis2Handle;
    iEvent.getByToken(axis2Token_, axis2Handle);
    edm::Handle<edm::ValueMap<int>> multHandle;
    iEvent.getByToken(multToken_, multHandle);

    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByToken(EDMGenJetsToken_, genJets);

    // Create a vector to add the jets to
    vector<JetIndexed> sortedJets;
    vector<JetIndexed> selectedJets;

    // Loop over the jets for pT-ordering within the event
    // Decide on using AK4 or AK8 jet algorithm. If AK8 -> replce jets with fatjets
    int iJetR = -1;
    for(pat::JetCollection::const_iterator jetIt = jets->begin(); jetIt!=jets->end(); ++jetIt) {
        const pat::Jet &jet = *jetIt;
	++iJetR;

	sortedJets.push_back( JetIndexed( jet, iJetR) );

        // Select
        if ( (jet.pt() > 30) && (fabs(jet.eta()) < 2.5) ) {
                selectedJets.push_back( JetIndexed( jet, iJetR) );
        }
    }

    // Sort the jets in pT-ordering
    std::sort(sortedJets.begin(), sortedJets.end(), higher_pT_sort());
    std::sort(selectedJets.begin(), selectedJets.end(), higher_pT_sort());

    // Loop over the selected jets in pT order
    for (unsigned int ptIdx = 0; ptIdx < selectedJets.size(); ++ptIdx) {
        // Make cuts on the event level
        if (sortedJets.size() < 2) continue;
        if (fabs(sortedJets[0].jet.eta()) > 2.5 || fabs(sortedJets[1].jet.eta()) > 2.5) continue;
        if (fabs(sortedJets[0].jet.pt()) < 30 || fabs(sortedJets[1].jet.pt()) < 30) continue;

	JetIndexed idxJet = selectedJets[ptIdx];
	const pat::Jet j = idxJet.jet;
	int iJetRef = idxJet.eventIndex;

        //adding jet parameters to jet-based tree
        jetPt = j.pt();
        jetEta = j.eta();
        jetPhi = j.phi();
        jetMass = j.mass();

	jetRawPt = j.correctedJet("Uncorrected").pt();
	jetRawEta = j.correctedJet("Uncorrected").eta();
	jetRawPhi = j.correctedJet("Uncorrected").phi();
	jetRawMass = j.correctedJet("Uncorrected").mass();

	jetChargedHadronMult = j.chargedHadronMultiplicity();
	jetNeutralHadronMult = j.neutralHadronMultiplicity();
	jetChargedMult = j.chargedMultiplicity();
	jetNeutralMult = j.neutralMultiplicity();

	jetPtOrder = ptIdx;

	//jetID
	jetLooseID = 0;
	jetTightID = 0;

	Float_t nhf = j.neutralHadronEnergyFraction();
	Float_t nemf = j.neutralEmEnergyFraction();
	Float_t chf = j.chargedHadronEnergyFraction();
	Float_t cemf = j.chargedEmEnergyFraction();
	unsigned int numconst = j.chargedMultiplicity() + j.neutralMultiplicity();
	unsigned int chm = j.chargedMultiplicity();

	if (abs(j.eta())<=2.7 && (numconst>1 && nhf<0.99 && nemf<0.99) && ((abs(j.eta())<=2.4 && chf>0 && chm>0 && cemf<0.99) || abs(j.eta())>2.4)) {
		jetLooseID = 1;
		if (nhf<0.90 && nemf<0.90) {
			jetTightID = 1;
		}
	}

        //add variables for deltaPhi and deltaEta for the two leading jets of the event
        dPhiJetsLO = deltaPhi(sortedJets[0].jet.phi(), sortedJets[1].jet.phi());
        dEtaJetsLO = sortedJets[0].jet.eta() - sortedJets[1].jet.eta();

        //the alpha variable is the third jet's pT divided by the average of the two leading jets' pT
        alpha = 0;
        //first make sure that there are at least 3 jets in the event
        if(sortedJets.size() > 2) {
                Float_t leadingPtAvg = (sortedJets[0].jet.pt() + sortedJets[1].jet.pt()) * 0.5;
                alpha = sortedJets[2].jet.pt() / leadingPtAvg;
        }

	//assign flavours for each jet
	partonFlav = abs(j.partonFlavour());
	hadronFlav = abs(j.hadronFlavour());
        physFlav = 0;
	if (j.genParton()) physFlav = j.genParton()->pdgId();

	isPartonUDS = 0;
	isPartonG = 0;
	isPartonOther = 0;
	isPhysUDS = 0;
	isPhysG = 0;
	isPhysOther = 0;

	//parton definition for flavours
	if(partonFlav == 1 || partonFlav == 2 || partonFlav == 3) {
		isPartonUDS = 1;
	} else if(partonFlav == 21) {
		isPartonG = 1;
	} else {
		isPartonOther = 1;
	}

	//physics definition for flavours
	if(abs(physFlav) == 1 || abs(physFlav) == 2 || abs(physFlav) == 3) {
		isPhysUDS = 1;
	} else if(abs(physFlav) == 21) {
		isPhysG = 1;
	} else {
		isPhysOther = 1;
	}

	edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(jets, iJetRef));
        jetQGl = (*qglHandle)[jetRef];
	QG_ptD = (*ptDHandle)[jetRef];
	QG_axis2 = (*axis2Handle)[jetRef];
	QG_mult = (*multHandle)[jetRef];

        //adding event information to jet-based tree
        event = iEvent.id().event();
        run = iEvent.id().run();
        lumi = iEvent.id().luminosityBlock();

        pthat = 0;
        if (genEventInfo->hasBinningValues()) {
		pthat = genEventInfo->binningValues()[0];
        }
        eventWeight = genEventInfo->weight();

	eventJetMult = selectedJets.size();

	// Loop over the pf candidates contained inside the jet (first sorting them in pT-order)
        // Here the jet girth is also calculated
        jetGirth = 0;

	std::vector<reco::CandidatePtr> pfCands = j.daughterPtrVector();
	std::sort(pfCands.begin(), pfCands.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt(); });
	int njetpf(0);

        unsigned int pfCandsSize = pfCands.size();
	for (unsigned int i = 0; i < pfCandsSize; ++i) {
		const pat::PackedCandidate &pf = dynamic_cast<const pat::PackedCandidate &>(*pfCands[i]);
                float dEta = pf.eta()-j.eta();
                float dPhi = deltaPhi(pf.phi(), j.phi());
                float dY = pf.rapidity() - j.rapidity();

		jetPF_pT[njetpf] = pf.pt();
                jetPF_pTrel[njetpf] = pf.pt() / j.pt();
                jetPF_mass[njetpf] = pf.mass();
                jetPF_dR[njetpf] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
                jetPF_dTheta[njetpf] = std::atan2(dPhi, dEta);
		jetPF_id[njetpf] = pf.pdgId();
                jetPF_fromPV[njetpf] = pf.fromPV();

                jetGirth += sqrt(dY*dY + dPhi*dPhi) * pf.pt()/j.pt();
		++njetpf;
	}
	jetMult = njetpf;


        // MC genJet & genPF variables
        genJetPt = 0;
        genJetEta = 0;
        genJetPhi = 0;
        genJetMass = 0;
        int ng(0);

        if(j.genJet()) {
                const reco::GenJet* gj = j.genJet();
                genJetPt = gj->pt();
                genJetEta = gj->eta();
                genJetPhi = gj->phi();
                genJetMass = gj->mass();


                // TLorentzVector g(0,0,0,0);
                // Loop over the genJet constituents after sorting them in pT-order
                std::vector<const pat::PackedGenParticle*> genParticles;
                for (unsigned int i = 0; i < gj->numberOfDaughters(); ++i) {
                        const pat::PackedGenParticle* genParticle = dynamic_cast<const pat::PackedGenParticle*>(gj->daughter(i));
                        genParticles.push_back( genParticle );
                }
                std::sort(genParticles.begin(), genParticles.end(), [](const pat::PackedGenParticle* p1, const pat::PackedGenParticle* p2) {return p1->pt() > p2->pt(); });

                unsigned int genParticlesSize = genParticles.size();
                for (unsigned int i = 0; i != genParticlesSize; ++i) {
                        const pat::PackedGenParticle* genParticle = dynamic_cast<const pat::PackedGenParticle*>(genParticles[i]);

                        genPF_pT[ng] = genParticle->pt();

                        float deltaEta = (genParticle->eta()-gj->eta());
                        //float deltaPhi = std::fabs((*packed)[i].phi()-j.phi());
                        //if (deltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
                        float DeltaPhi = deltaPhi(genParticle->phi(), gj->phi());
                        // later: TLorentzVector::DeltaPhi() => dPhi = pf.DeltaPhi(j);

                        genPF_pT[ng] = genParticle->pt();
                        genPF_dR[ng] = deltaR(gj->eta(), gj->phi(), genParticle->eta(), genParticle->phi());
                        genPF_dTheta[ng] = std::atan2(DeltaPhi, deltaEta);
                        genPF_mass[ng] = genParticle->mass();
                        genPF_id[ng] = genParticle->pdgId();
                        ++ng;

                        // if ( genPF_dR[i] < 0.4 )
                        // g += TLorentzVector((*packed)[i].px(), (*packed)[i].py(), (*packed)[i].pz(), (*packed)[i].energy());
                }
                ngenv = ng;
        }


	// Create the jet images that include particles also outside the jet
        // PF Particle loop
                if (!(kMaxPF < pfs->size()))
                assert(kMaxPF > pfs->size());
                int npfs(0);

        unsigned int pfsSize = pfs->size();
        for (unsigned int i = 0; i != pfsSize; ++i) {
            const pat::PackedCandidate &pf = (*pfs)[i];

            float deltaEta = (pf.eta()-j.eta());
            float DeltaPhi = deltaPhi(pf.phi(),j.phi());
	    //if (DeltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
            // later: TLorentzVector::DeltaPhi() => dPhi = pf.DeltaPhi(j);

            if ( (fabs(deltaEta) > 1.0) || (fabs(DeltaPhi) > 1.0) ) continue;

                PF_pT[npfs] = pf.pt();
                PF_dR[npfs] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
                PF_dTheta[npfs] = std::atan2(DeltaPhi, deltaEta);
                PF_dPhi[npfs] = DeltaPhi;
                PF_dEta[npfs] = deltaEta;
                PF_mass[npfs] = pf.mass();
                PF_id[npfs] = pf.pdgId();
                PF_fromPV[npfs] =  pf.fromPV();
                ++npfs;

        } // for jetImage pfs
        nPF = npfs;

	jetTree->Fill();

    } // for jets


} // analyze


//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
