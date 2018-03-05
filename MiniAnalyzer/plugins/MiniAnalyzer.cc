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

#include <iostream>

// user include files
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

        int npfv, ngenv, ncPF, nnPF;
        struct PFV {float pT,dR,dTheta, mass;};
        static const int kMaxPF = 5000;

	// Jet constituent variables
	Float_t jetcPF_pT[kMaxPF];
	Float_t jetcPF_pTrel[kMaxPF];
	Float_t jetcPF_dR[kMaxPF];
	Float_t jetcPF_dTheta[kMaxPF];
	Float_t jetcPF_mass[kMaxPF];

	Float_t jetnPF_pT[kMaxPF];
	Float_t jetnPF_pTrel[kMaxPF];
	Float_t jetnPF_dR[kMaxPF];
	Float_t jetnPF_dTheta[kMaxPF];
	Float_t jetnPF_mass[kMaxPF];
/*
        Float_t pf_pT[kMaxPF];
        Float_t pf_dR[kMaxPF];
        Float_t pf_dTheta[kMaxPF];
        Float_t pf_mass[kMaxPF];
*/
	// Gen particle variables
        Float_t gen_pT[kMaxPF];
        Float_t gen_dR[kMaxPF];
        Float_t gen_dTheta[kMaxPF];
        Float_t gen_mass[kMaxPF];

	//Jet image variables
	Float_t cPF_pT[kMaxPF];
	Float_t cPF_dR[kMaxPF];
	Float_t cPF_dTheta[kMaxPF];
	Float_t cPF_mass[kMaxPF];

	Float_t nPF_pT[kMaxPF];
	Float_t nPF_dR[kMaxPF];
	Float_t nPF_dTheta[kMaxPF];
	Float_t nPF_mass[kMaxPF];

	//int getMatchedPartonGen(edm::Event const&, std::vector<reco::GenJet>::const_iterator&);
	

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

	Float_t jetRawPt;
	Float_t jetRawEta;
	Float_t jetRawPhi;
	Float_t jetRawMass;

	unsigned int jetChargedHadronMult;
	unsigned int jetNeutralHadronMult;
	unsigned int jetChargedMult;
	unsigned int jetNeutralMult;

	unsigned int jetLooseID;
	unsigned int jetTightID;
	
        unsigned int event;
        unsigned int run;
        unsigned int lumi;

	unsigned int eventJetMult;
	unsigned int jetPtOrder;
        
	unsigned int partonFlav;
	unsigned int hadronFlav;
	unsigned int physFlav;

	unsigned int isPartonUDS;
	unsigned int isPartonG;
	unsigned int isPartonOther;
	unsigned int isPhysUDS;
	unsigned int isPhysG;
	unsigned int isPhysOther;

        Float_t genPt;
        Float_t genEta;
        Float_t genPhi;
        Float_t genMass;
        
        Float_t pthat;

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
    qglToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "qgLikelihood"))),
    ptDToken_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "ptD"))),
    axis2Token_(consumes<edm::ValueMap<float>>(edm::InputTag("QGTagger", "axis2"))),
    multToken_(consumes<edm::ValueMap<int>>(edm::InputTag("QGTagger", "mult")))

{

//    mMaxY = iConfig.getParameter<double>("maxY");
    mMinGenPt = iConfig.getUntrackedParameter<double>("minGenPt", 30);

    //sets the output file, tree and the parameters to be saved to the tree

    outputFile = new TFile("nanotuples2_PU.root","recreate");
    jetTree = new TTree("jetTree", "Jet tree");

    jetTree->Branch("jetPt", &jetPt, "jetPt/F");
    jetTree->Branch("jetEta", &jetEta, "jetEta/F");
    jetTree->Branch("jetPhi", &jetPhi, "jetPhi/F");
    jetTree->Branch("jetMass", &jetMass, "jetMass/F");

    jetTree->Branch("jetRawPt", &jetRawPt, "jetRawPt/F");
    jetTree->Branch("jetRawEta", &jetRawEta, "jetRawEta/F");
    jetTree->Branch("jetRawPhi", &jetRawPhi, "jetRawPhi/F");
    jetTree->Branch("jetRawMass", &jetRawMass, "jetRawMass/F");

    jetTree->Branch("jetChargedHadronMult", &jetChargedHadronMult, "jetChargedHadronMult/I");
    jetTree->Branch("jetNeutralHadronMult", &jetNeutralHadronMult, "jetNeutralHadronMult/I");
    jetTree->Branch("jetChargedMult", &jetChargedMult, "jetChargedMult/I");    
    jetTree->Branch("jetNeutralMult", &jetNeutralMult, "jetNeutralMult/I");

    jetTree->Branch("jetcPF_pT", &jetcPF_pT, "jetcPF_pT[jetChargedMult]/F");
    jetTree->Branch("jetcPF_pTrel", &jetcPF_pTrel, "jetcPF_pTrel[jetChargedMult]/F");
    jetTree->Branch("jetcPF_dR", &jetcPF_dR, "jetcPF_dR[jetChargedMult]/F");
    jetTree->Branch("jetcPF_dTheta", &jetcPF_dTheta, "jetcPF_dTheta[jetChargedMult]/F");
    jetTree->Branch("jetcPF_mass", &jetcPF_mass, "jetcPF_mass[jetChargedMult]/F");

    jetTree->Branch("jetnPF_pT", &jetnPF_pT, "jetnPF_pT[jetNeutralMult]/F");
    jetTree->Branch("jetnPF_pTrel", &jetnPF_pTrel, "jetnPF_pTrel[jetNeutralMult]/F");
    jetTree->Branch("jetnPF_dR", &jetnPF_dR, "jetnPF_dR[jetNeutralMult]/F");
    jetTree->Branch("jetnPF_dTheta", &jetnPF_dTheta, "jetnPF_dTheta[jetNeutralMult]/F");
    jetTree->Branch("jetnPF_mass", &jetnPF_mass, "jetnPF_mass[jetNeutralMult]/F");

    jetTree->Branch("jetLooseID", &jetLooseID, "jetLooseID/I");
    jetTree->Branch("jetTightID", &jetTightID, "jetTightID/I");

    jetTree->Branch("genPt", &genPt, "genPt/F");
    jetTree->Branch("genEta", &genEta, "genEta/F");
    jetTree->Branch("genPhi", &genPhi, "genPhi/F");
    jetTree->Branch("genMass", &genMass, "genMass/F");
    
    jetTree->Branch("pthat", &pthat, "pthat/F");

    jetTree->Branch("jetPtOrder", &jetPtOrder, "jetPtOrder/I");
    
    jetTree->Branch("event", &event, "event/l");
    jetTree->Branch("run", &run, "run/I");
    jetTree->Branch("lumi", &lumi, "lumi/I");

    jetTree->Branch("eventJetMult", &eventJetMult, "eventJetMult/I");
    jetTree->Branch("jetPtOrder", &jetPtOrder, "jetPtOrder/I");

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

/*
    jetTree->Branch("np",&npfv,"np/I");
    //jetTree->Branch("pf", pfv, "pT[np]/F:dR[np]/F:dTheta[np]/F:"
    //	    "mass[np]/F");
    jetTree->Branch("pf_pT", &pf_pT, "pT[np]/F");
    jetTree->Branch("pf_dR", &pf_dR, "dR[np]/F");
    jetTree->Branch("pf_dTheta", &pf_dTheta, "dTheta[np]/F");
    jetTree->Branch("pf_mass", &pf_mass, "mass[np]/F");
*/
    // Jet image variables
    jetTree->Branch("ncPF", &ncPF, "ncPF/I");
    jetTree->Branch("cPF_pT", &cPF_pT, "cPF_pT[ncPF]/F");
    jetTree->Branch("cPF_dR", &cPF_dR, "cPF_dR[ncPF]/F");
    jetTree->Branch("cPF_dTheta", &cPF_dTheta, "cPF_dTheta[ncPF]/F");
    jetTree->Branch("cPF_mass", &cPF_mass, "cPF_mass[ncPF]/F");

    jetTree->Branch("nnPF", &nnPF, "nnPF/I");
    jetTree->Branch("nPF_pT", &nPF_pT, "nPF_pT[nnPF]/F");
    jetTree->Branch("nPF_dR", &nPF_dR, "nPF_dR[nnPF]/F");
    jetTree->Branch("nPF_dTheta", &nPF_dTheta, "nPF_dTheta[nnPF]/F");
    jetTree->Branch("nPF_mass", &nPF_mass, "nPF_mass[nnPF]/F");

    jetTree->Branch("ng",&ngenv,"ng/I");
    //jetTree->Branch("gen", genv, "pT[ng]/F:dR[ng]/F:dTheta[ng]/F:"
    //	    "mass[ng]/F");
    
    jetTree->Branch("gen_pT", &gen_pT, "gen_pT[ng]/F");
    jetTree->Branch("gen_dR", &gen_dR, "gen_dR[ng]/F");
    jetTree->Branch("gen_dTheta", &gen_dTheta, "gen_dTheta[ng]/F");
    jetTree->Branch("gen_mass", &gen_mass, "gen_mass[ng]/F");
    
}

MiniAnalyzer::~MiniAnalyzer()
{

//at the end, write data into tree
	jetTree->Write();
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
    //edm::Handle<reco::JetFlavourInfoMatchingCollection> theJetFlavourInfos;
    //iEvent.getByToken(jetFlavourInfosToken_, theJetFlavourInfos);

/*
	//Genjet parton flavor loop
	for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
    	if (i_gen->pt() > 30 && (fabs(i_gen->y()) < mMaxY)) {
    		//int FlavourGen = getMatchedPartonGen(event, i_gen); TO DO!!
    		//if(FlavourGen < -100) cout<<FlavourGen<<" "<<i_gen->pt()<<" "<<i_gen->eta()<<" "<<i_gen->phi()<<endl; TO DO!!
    		//GenFlavour.push_back(FlavourGen); TO DO!!
    	}


	//Genjet hadron flavor loop
        for ( reco::JetFlavourInfoMatchingCollection::const_iterator j = theJetFlavourInfos->begin(); j != theJetFlavourInfos->end(); ++j ) {
    	if (i_gen->pt() > mMinGenPt && fabs(i_gen->y()) < mMaxY) {
	    	reco::JetFlavourInfo aInfo = (*j).second;
	    	int FlavourGenHadron = aInfo.getHadronFlavour();
	      	if(FlavourGenHadron==5) cout<<FlavourGenHadron<<" "<<aJet->pt()<<" "<<aJet->eta()<<" "<<aJet->phi()<<" HADRONFLAV"<<endl;
	      	GenHadronFlavour.push_back(FlavourGenHadron);
         	}
        } //for genjet hadron flavor
*/
//	}


    // Create a vector to add the jets to
    vector<JetIndexed> selectedJets;
    
    // Loop over the jets for pT-ordering within the event
    int iJetR = -1;
    for(pat::JetCollection::const_iterator jetIt = jets->begin(); jetIt!=jets->end(); ++jetIt) {
        const pat::Jet &jet = *jetIt;
	++iJetR;

	// Select
	if ( (jet.pt() < 20) || (fabs(jet.eta()) > 1.3) ) continue;

	selectedJets.push_back( JetIndexed(jet, iJetR) );
    }

    // Sort the jets in pT-ordering
    std::sort(selectedJets.begin(), selectedJets.end(), higher_pT_sort());


    // Loop over the jets
//    int iJetRef = -1;
    // Decide on using AK4 or AK8 jet algorithm. If AK8 -> replce jets with fatjets
//    for(pat::JetCollection::const_iterator jetIt = jets->begin(); jetIt!=jets->end(); ++jetIt) {
//	const pat::Jet &j = *jetIt;
//	++iJetRef;

      	// Select
//      	if ( (j.pt() < 20) || (fabs(j.eta()) > 1.3) ) continue;


// Loop over the selected jets in pT order
    for (unsigned int ptIdx = 0; ptIdx < selectedJets.size(); ++ptIdx) {

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

	//assign flavours for each jet
	partonFlav = abs(j.partonFlavour());
	hadronFlav = abs(j.hadronFlavour());
        physFlav = 0;
	if (j.genParton()) physFlav = abs(j.genParton()->pdgId());

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
	if(physFlav == 1 || physFlav == 2 || physFlav == 3) {
		isPhysUDS = 1;
	} else if(physFlav == 21) {
		isPhysG = 1;
	} else {
		isPhysOther = 1;
	}

	edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(jets, iJetRef));
        jetQGl = (*qglHandle)[jetRef];
	QG_ptD = (*ptDHandle)[jetRef];
	QG_axis2 = (*axis2Handle)[jetRef];
	QG_mult = (*multHandle)[jetRef];

	genPt = 0;
	genEta = 0;
	genPhi = 0;
	genMass = 0;
        //adding MCjet parameters
	if(j.genJet()) {
		genPt = j.genJet()->pt();
		genEta = j.genJet()->eta();
		genPhi = j.genJet()->phi();
		genMass = j.genJet()->mass();
	}

        //adding event information to jet-based tree
        event = iEvent.id().event();
        run = iEvent.id().run();
        lumi = iEvent.id().luminosityBlock();

	eventJetMult = selectedJets.size();
       
	// Loop over the pf candidates contained inside the jet (first sorting them in pT-order)
	std::vector<reco::CandidatePtr> pfCands = j.daughterPtrVector();
	std::sort(pfCands.begin(), pfCands.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt(); });
	int nnpf(0);
	int ncpf(0);
	for (unsigned int i = 0; i < pfCands.size(); ++i) {
		const pat::PackedCandidate &pf = dynamic_cast<const pat::PackedCandidate &>(*pfCands[i]);
                float dEta = pf.eta()-j.eta();
                float dPhi = deltaPhi(pf.phi(), j.phi());		
		if (pf.charge() == 0) {
                	jetnPF_pT[nnpf] = pf.pt();
                	jetnPF_pTrel[nnpf] = pf.pt() / j.pt();
                	jetnPF_mass[nnpf] = pf.mass();
                	jetnPF_dR[nnpf] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
                	jetnPF_dTheta[nnpf] = std::atan2(dPhi, dEta);
			++nnpf;
		} else {
                	jetcPF_pT[ncpf] = pf.pt();
                	jetcPF_pTrel[ncpf] = pf.pt() / j.pt();
                	jetcPF_mass[ncpf] = pf.mass();
                	jetcPF_dR[ncpf] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
                	jetcPF_dTheta[ncpf] = std::atan2(dPhi, dEta);
			++ncpf;
		}
	}	
	jetNeutralMult = nnpf;
	jetChargedMult = ncpf;
	

	// Create the jet images that include particles also outside the jet
	// particle loop starts here
/*		if (!(kMaxPF < pfs->size()))
		assert(kMaxPF > pfs->size());
		int np(0);

        for (unsigned int i = 0; i != pfs->size(); ++i) {
            const pat::PackedCandidate &pf = (*pfs)[i];

            float deltaEta = (pf.eta()-j.eta());
            float DeltaPhi = deltaPhi(pf.phi(),j.phi());
	    	//if (DeltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
	    	// later: TLorentzVector::DeltaPhi() => dPhi = pf.DeltaPhi(j);

            if ( (fabs(deltaEta) > 1.0) || (fabs(DeltaPhi) > 1.0) ) continue;
                
            pf_pT[np] = pf.pt();
            pf_dR[np] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
            pf_dTheta[np] = std::atan2(DeltaPhi, deltaEta);
            pf_mass[np] = pf.mass();                        
	    ++np;

        } // for pfs
	
		npfv = np;
*/
        // PF Particle loop
                if (!(kMaxPF < pfs->size()))
                assert(kMaxPF > pfs->size());
                int nc(0); // charged
		int nn(0);  // neutral

        for (unsigned int i = 0; i != pfs->size(); ++i) {
            const pat::PackedCandidate &pf = (*pfs)[i];

            float deltaEta = (pf.eta()-j.eta());
            float DeltaPhi = deltaPhi(pf.phi(),j.phi());
	    //if (DeltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
            // later: TLorentzVector::DeltaPhi() => dPhi = pf.DeltaPhi(j);

            if ( (fabs(deltaEta) > 1.0) || (fabs(DeltaPhi) > 1.0) ) continue;

	    if (pf.charge() == 0) {
                nPF_pT[nn] = pf.pt();
                nPF_dR[nn] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
                nPF_dTheta[nn] = std::atan2(DeltaPhi, deltaEta);
                nPF_mass[nn] = pf.mass();
                ++nn;
	    } else {
                cPF_pT[nc] = pf.pt();
                cPF_dR[nc] = deltaR(j.eta(), j.phi(), pf.eta(), pf.phi());
                cPF_dTheta[nc] = std::atan2(DeltaPhi, deltaEta);
                cPF_mass[nc] = pf.mass();
                ++nc;
	    }

        } // for pfs

                ncPF = nc;
		nnPF = nn;

	// Gen particle loop
		TLorentzVector g(0,0,0,0);
		int ng(0);

        for (unsigned int i = 0; i != packed->size(); ++i) {

            float deltaEta = ((*packed)[i].eta()-j.eta());
            //float deltaPhi = std::fabs((*packed)[i].phi()-j.phi());
	    	//if (deltaPhi>(M_PI)) deltaPhi-=(2*M_PI);
	    	float DeltaPhi = deltaPhi((*packed)[i].phi(),j.phi());
	    	// later: TLorentzVector::DeltaPhi() => dPhi = pf.DeltaPhi(j);

            if ( (fabs(deltaEta) > 1.0) || (fabs(DeltaPhi) > 1.0) ) continue;
                
            gen_pT[ng] = (*packed)[i].pt();
            gen_dR[ng] = deltaR(j.eta(), j.phi(), (*packed)[i].eta(), (*packed)[i].phi());
	    	gen_dTheta[ng] = std::atan2(DeltaPhi, deltaEta);
            gen_mass[ng] = (*packed)[i].mass();
	    	++ng;
	    
	    	if ( gen_dR[i] < 0.4 )
	    	g += TLorentzVector((*packed)[i].px(), (*packed)[i].py(), (*packed)[i].pz(), (*packed)[i].energy());
        
        } // for packed gen pfs
	
		ngenv = ng;

	    jetTree->Fill();

    } // for jets


} // analyze


///// Matching Flavour //////
/*int MiniAnalyzer::getMatchedPartonGen(edm::Event const& event, GenJetCollection::const_iterator i_gen)
{

  int jetFlavour = -100;
  bool switchB = 0;
  bool switchC = 0;

  edm::Handle<reco::GenParticleCollection> genParticles;
  event.getByToken(mgenParticles, genParticles);

  for (size_t i = 0; i < genParticles->size(); ++i) {
      const GenParticle & genIt = (*genParticles)[i];
      int pdgId = genIt.pdgId();
      double DeltaR = deltaR(genIt.p4().eta(), genIt.p4().phi(), i_gen->eta(), i_gen->phi());
      double DeltaRmin = 0.3;
      if (DeltaR < DeltaRmin ){

	DeltaRmin=DeltaR;
	if(abs(pdgId)==5) { jetFlavour=5; switchB=true;}
	if(abs(pdgId)==4) { jetFlavour=4; switchC=true;}
	if(abs(pdgId)<=3 && abs(pdgId)>=1) { jetFlavour=1; }
	if(abs(pdgId)==21) { jetFlavour=21; }
      }
      
      if (switchB) {jetFlavour=5;}
      if (switchC && !switchB) {jetFlavour=4;}

  }

  return jetFlavour;
}
*/


//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
