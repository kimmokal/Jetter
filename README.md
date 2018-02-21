# Jetter - a Jet nanotuple producer from CMS MiniAOD

This project is a CMSSW module producing mostly flat tuples from run 2 data/MC

(A little bit under construction)

So far the code works on lxplus and CMSSW_8_0_6, and requires grid permissions to access data:
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
cmsrel CMSSW_8_0_6
cd CMSSW_8_0_6/src
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

There are five different "groups" of data variables that are stored in the jet-based tuples. 

| Jet | GenJet | Particle Flow candidate (PF) | GenParticle | Event |
| :------------- | :------------- | :------------- | :------------- | :------------- |
| jetPt | genPt | pf_pT[np] | gen_pT[ng] | event |
| jetEta | getEta | pf_dR[np] | gen_dR[ng] | run |
| jetPhi | genPhi| pf_dTheta[np] | gen_dTheta[ng] | lumi |
| jetMass | genMass| pf_mass[np] | gen_mass[ng] |
| | | npfv | ngenv |

As a general rule, variable with underscore refers to particles (to avoid confusion between the two get variables)



## Unfinished features

- [ ] Needed Variables
	- Charged Hadron Multiplicity
	- Charged Particle Transverse Momenta
	- Neutral Hadron Multiplicity
	- Neutral Particle Transverse Momenta
	- Relative pT of a charged jet constituent with respect to the jet pT
	- Relative pT of a neutral jet constituent with respect to the jet pT
	- ∆R(cPF) : ∆R to jet axis of a charged candidate
	- ∆R(nPF) : ∆R to jet axis of a neutral candidate
	- JetID
	- Phase space cut variables

- [ ] ```pdg_id``` for gen particles and condensed id for 
```
charged hadron = pi+/pi- 211/-211
neutral hadron -> 130
photon -> 111
```
- [ ] Vertex index (pf_ivtx)
	- Comparison of tracks, save the ivtx if a match for particle is found. Default value -1 (no vertex)
```
0 = 1st primary
>=1 : pu_vertex
-1 : no vertex
```

- [ ] Accounting for weight/pt-hat

- [ ] Not properly tested: Jet flavors for genjets
