# Jetter - a Jet nanotuple producer from CMS MiniAOD

This project is a CMSSW module producing mostly flat tuples from run 2 data/MC.

So far the code works on lxplus and CMSSW_8_0_26, and requires grid permissions to access data:
```voms-proxy-init -voms cms```


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

Change the input file in python/ConfFile_cfg.py and if desired, the output filename in plugins/MiniAnalyzer.cc .

To generate the tuples, run
```
    cmsRun python/ConfFile_cfg.py
```
And to view the tuples in ROOT
```
    root -l nanotuple_*
    TBrowser a
```
from there use the graphical interface.


## Content of the tuples

Data variables are saved to ROOT trees jet by jet.

In this version, the reconstructed jets are AK4 jets clustered from Particle Flow candidates. The standard L1+L2+L3+residual energy corrections are applied to the jets and pileup is reduced using the CHS algorithm.

Some of the variables (such as jetRawPhi) are admittedly redundant.

| Data variable | Description |
| :------------ | :---------- |
| jetPt | Transverse momentum of the jet |
| jetEta | Pseudorapidity (η) of the jet |
| jetPhi | Azimuthal angle (ϕ) of the jet |
| jetMass | Mass of the jet |
| jetGirth | Girth of the jet (as defined in arXiv:1106.3076 [hep-ph]) |
| jetRawPt | Transverse momentum of the jet before the energy corrections|
| jetRawEta | Pseudorapidity (η) of the jet before the energy corrections |
| jetRawPhi | The azimuthal angle (ϕ) of the jet before the energy corrections |
| jetRawMass | Mass of the jet before the energy corrections |
| jetLooseID | Binary variable indicating whether the jet passes Loose JetID criteria (https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016) |
| jetTightID | Binary variable indicating whether the jet passes Tight JetID criteria (see above link) |
| jetChargedHadronMult | Multiplicity of charged hadron jet constituents |
| jetNeutralHadronMult | Multiplicity of neutral hadron jet constituents |
| jetChargedMult | Multiplicity of charged jet constituents |
| jetNeutralMult | Multiplicity of neutral jet constituents |
| jetMult | Multiplicity of jet constituents |
| genJetPt | Transverse momentum of the matching generator level jet (i.e. MC truth) |
| genJetEta | Pseudorapidity (η) of the matching generator level jet |
| genJetPhi | Azimuthal angle (ϕ) of the matching generator level jet |
| genJetMass | Mass of the matching generator level jet |
| pthat | Transverse momentum of the generated hard process |
| eventWeight | Weight assigned to the generated event |
| jetPtOrder | Indicates the ranking number of the jet, as the jets are ordered by their transverse momenta within a single event |
| event | Event number |
| run | Run number |
| lumi | Luminosity block |
| eventJetMult | Multiplicity of jets in the event |
| dPhiJetsLO | The phi difference of the two leading jets |
| dEtaJetsLO | The eta difference of the two leading jets |
| alpha | If there are at least 3 jets in the event, alpha is the third jet's pT divided by the average pT of the two leading jets |
| partonFlav | Flavour of the jet, as defined by the parton-based definition (https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools) |
| hadronFlav | Flavour of the jet, as defined by the hadron-based definition (see above link) |
| physFlav | Flavour of the jet, as defined by the 'physics' definition (see above link) |
| jetQGl | Quark-Gluon likelihood as determined by the current BDT discriminator (https://twiki.cern.ch/twiki/bin/viewauth/CMS/QuarkGluonLikelihood) |
| QG_ptD | Jet energy variable (see above link) |
| QG_axis2 | Minor axis of the jet (see above link) |
| QG_mult | Jet constituent multiplicity (see above link) |
| isPartonUDS | Indicates light quark (Up, Down, Strange) jets: partonFlav = 1, 2, 3 |
| isPartonG | Indicates gluon jets: partonFlav = 21 |
| isPartonOther | Indicates any other kind of jet: partonflav != 1, 2, 3, 21 |
| isPhysUDS | Indicates light quark (Up, Down, Strange) jets: physFlav = 1, 2, 3 |
| isPhysG | Indicates gluon jets: physFlav = 21 |
| isPhysOther | Indicates any other kind of jet: physFlav != 1, 2, 3, 21|
| jetPF_pT[jetMult] | Transverse momentum of a jet constituent (Particle Flow (PF) candidate) |
| jetPF_pTrel[jetMult] | Relative transverse momentum of a jet constituent with respect to the jet transverse momentum (jetPF_pT / jetPt)|
| jetPF_dR[jetMult] | Distance of a jet constituent to the center of the jet |
| jetPF_dTheta[jetMult] | Polar angle (θ) of a jet constituent |
| jetPF_mass[jetMult] | Mass of a jet constituent |
| jetPF_id[jetMult] | Particle identifier of a jet constituent, as defined by the PDG numbering scheme (For PF candidates: charged hadron = 211/-211; neutral hadron = 130; photon = 22) |
| jetPF_fromPV[jetMult] | A number indicating how tightly a particle is associated with the primary vertex (ranges from 3 to 0) |
| ng | Number of generator level particles in the corresponding generator level jet |
| genPF_pT[ng] | Transverse momentum of a gen level jet constituent |
| genPF_dR[ng] | Distance of a gen level jet constituent to the center of the gen jet |
| genPF_dTheta[ng] | Polar angle (θ) of a gen level jet constituent |
| genPF_mass[ng] | Mass of a gen level jet constituent |
| genPF_id[ng] | Particle identifier of a gen level jet constituent, as defined by the PDG numbering scheme |

There are also variables for 'jet images' (arXiv:1612.01551 [hep-ph]), for which every Particle Flow (PF) candidate within the area Δη<1 && Δϕ<1 around the jet center is saved. This includes particles which are not part of the jet.

| Data variable | Description |
| :------------ | :---------- |
| nPF | Number of PF candidates in the jet image |
| PF_pT[nPF] | Transverse momentum of a PF candidate |
| PF_dR[nPF] | Distance of a PF candidate to the center of the jet |
| PF_dTheta[nPF] | Polar angle (θ) of a PF candidate |
| PF_mass[nPF] | Mass of a PF candidate |
| PF_id[nPF] | Particle identifier as defined by the PDG convention (For PF candidates: charged hadron = 211/-211; neutral hadron = 130; photon = 22) |
| PF_fromPV[nPF] | A number indicating how tightly a particle is associated with the primary vertex (ranges from 3 to 0) |



## Unfinished features

- [ ] Needed Variables
	- Phase space cut variables
