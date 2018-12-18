# Jetter - a Jet nanotuple producer from CMS MiniAOD

This is a CMSSW module for producing mostly flat tuples from 13 TeV Run2 MC samples.

The code works on lxplus and CMSSW_8_0_26, and requires grid permissions to access data:
```voms-proxy-init --voms cms```


## Setting up

First setup your GitHub:
```
git config --global user.name 'Your Name'
git config --global user.email 'your@ema.il'
git config --global user.github 'username'
```
Make sure you register your SSH key: https://help.github.com/articles/generating-ssh-keys


Creating the work area
```
ssh -Y CERNusername@lxplus.cern.ch

mkdir WorkingArea
cd WorkingArea
cmsrel CMSSW_8_0_26
cd CMSSW_8_0_26/src
cmsenv
git clone https://github.com/kimmokal/Jetter/
scram b
cd Jetter/MiniAnalyzer
```

## Running the code

If desired, change the input file and the output name in ```python/ConfFile_cfg.py```.

To generate the tuples, execute
```
    cmsRun python/ConfFile_cfg.py
```

## Content of the tuples

Variables are saved to ROOT trees jet by jet.

In this version, the reconstructed jets are AK4 jets clustered from Particle Flow candidates. The standard L1+L2+L3+residual energy corrections are applied to the jets and pileup is reduced using the CHS algorithm.

| Data variable | Type | Description |
| :------------ | :---------- | :---------- |
| jetPt | Float_t | Transverse momentum of the jet |
| jetEta | Float_t | Pseudorapidity (η) of the jet |
| jetPhi | Float_t | Azimuthal angle (ϕ) of the jet |
| jetMass | Float_t | Mass of the jet |
| jetGirth | Float_t | Girth of the jet (as defined in arXiv:1106.3076 [hep-ph]) |
| jetArea | Float_t | Catchment area of the jet; used for jet energy corrections |
| jetRawPt | Float_t | Transverse momentum of the jet before the energy corrections |
| jetRawMass | Float_t | Mass of the jet before the energy corrections |
| jetLooseID | UInt_t | Binary variable indicating whether the jet passes Loose JetID criteria (https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016) |
| jetTightID | UInt_t | Binary variable indicating whether the jet passes Tight JetID criteria (see above link) |
| jetGenMatch | UInt_t | 1: if a matched generator level jet exists; 0: if no match was found |
| jetQGl | Float_t | Quark-Gluon likelihood as determined by the current BDT discriminator (https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood) |
| QG_ptD | Float_t | Jet energy variable (see above link) |
| QG_axis2 | Float_t | Minor axis of the jet (see above link) |
| QG_mult | UInt_t | Jet constituent multiplicity with additional cuts (see above link) |
| partonFlav | Int_t | Flavour of the jet, as defined by the parton-based definition (https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools) |
| hadronFlav | Int_t | Flavour of the jet, as defined by the hadron-based definition (see above link) |
| physFlav | Int_t | Flavour of the jet, as defined by the 'physics' definition (see above link) |
| isPartonUDS | UInt_t | Indicates light quark (Up, Down, Strange) jets: partonFlav = 1, 2, 3 |
| isPartonG | UInt_t | Indicates gluon jets: partonFlav = 21 |
| isPartonOther | UInt_t | Indicates any other kind of jet: partonflav != 1, 2, 3, 21 |
| isPhysUDS | UInt_t | Indicates light quark (Up, Down, Strange) jets: physFlav = 1, 2, 3 |
| isPhysG | UInt_t | Indicates gluon jets: physFlav = 21 |
| isPhysOther | UInt_t | Indicates any other kind of jet: physFlav != 1, 2, 3, 21 |
| jetChargedHadronMult | UInt_t | Multiplicity of charged hadron jet constituents |
| jetNeutralHadronMult | UInt_t | Multiplicity of neutral hadron jet constituents |
| jetChargedMult | UInt_t | Multiplicity of charged jet constituents |
| jetNeutralMult | UInt_t | Multiplicity of neutral jet constituents |
| jetMult | UInt_t | Multiplicity of jet constituents |
| nPF | UInt_t | Number of particle flow (PF) candidates (jet constituent particles reconstructed by the particle flow algorithm); contains all particles within |Δϕ<1| < 1 && |Δη<1| < 1 from the center of the jet |
| PF_pT[nPF] | Float_t | Transverse momentum of a PF candidate |
| PF_dR[nPF] | Float_t | Distance of a PF candidate to the center of the jet |
| PF_dTheta[nPF] | Float_t | Polar angle (θ) of a PF candidate |
| PF_dPhi[nPF] | Float_t | Azimuthal angle (ϕ) of a PF candidate |
| PF_dEta[nPF] | Float_t | Pseudorapidity (η) of a PF candidate |
| PF_mass[nPF] | Float_t | Mass of a PF candidate |
| PF_id[nPF] | Int_t | Generator level particle identifier for the particle flow candidates, as defined in the PDG particle numbering scheme |
| PF_fromPV[nPF] | UInt_t | A number indicating how tightly a particle is associated with the primary vertex (ranges from 3 to 0) |
| PF_fromAK4Jet | UInt_t | 1: if the particle flow candidate is a constituent of the reconstructed AK4 jet; 0: if it is not a constituent of the jet |
| genJetPt | Float_t | Transverse momentum of the matched generator level jet |
| genJetEta | Float_t | Pseudorapidity (η) of the matched generator level jet |
| genJetPhi | Float_t | Azimuthal angle (ϕ) of the matched generator level jet |
| genJetMass | Float_t | Mass of the matched generator level jet |
| nGenJetPF | UInt_t | Number of particles in the matched generator level jet |
| genPF_pT[nGenJetPF] | Float_t | Transverse momentum of a particle in the matched generator level jet |
| genPF_dR[nGenJetPF] | Float_t | Distance of a particle to the center of the matched generator level jet |
| genPF_dTheta[nGenJetPF] | Float_t | Polar angle (θ) of a particle in the matched generator level jet |
| genPF_mass[nGenJetPF] | Float_t | Mass of a particle in the matched generator level jet |
| genPF_id[nGenJetPF] | Int_t | Generator level particle identifier for the particles in the matched generator level jet, as defined in the PDG particle numbering scheme |
| eventJetMult | UInt_t | Multiplicity of jets in the event |
| jetPtOrder | UInt_t | Indicates the ranking number of the jet, as the jets are ordered by their transverse momenta within a single event |
| dPhiJetsLO | Float_t | The phi difference of the two leading jets |
| dEtaJetsLO | Float_t | The eta difference of the two leading jets |
| alpha | Float_t | If there are at least 3 jets in the event, alpha is the third jet's transverse momentum divided by the average transverse momentum of the two leading jets |
| event | ULong64_t | Event number |
| run | UInt_t | Run number |
| lumi | UInt_t | Luminosity block |
| pthat | Float_t | Transverse momentum of the generated hard process |
| eventWeight | Float_t | Weight assigned to the generated event |
| rhoAll | Float_t | The median density (in GeV/A) of pile-up contamination per event; computed from all PF candidates of the event |
| rhoCentral | Float_t | Same as above, computed from all PF candidates with |eta| < 2.5 |
| rhoCentralNeutral | Float_t | Same as above, computed from all neutral PF candidates with |eta| < 2.5 |
| rhoCentralChargedPileUp | Float_t | Same as above, computed from all PF charged hadrons associated to pileup vertices and with |eta| < 2.5 |
| PV_npvsGood | UInt_t | The number of good reconstructed primary vertices |
| Pileup_nPU | UInt_t | The number of pileup interactions that have been added to the event in the current bunch crossing |
| Pileup_nTrueInt | Float_t | The true mean number of the poisson distribution for this event from which the number of interactions in each bunch crossing has been sampled |

There are also variables for 'jet images' (arXiv:1612.01551 [hep-ph]), for which every Particle Flow (PF) candidate within the area Δη<1 && Δϕ<1 around the jet center is saved. This includes particles which are not part of the jet.
