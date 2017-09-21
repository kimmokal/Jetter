# Jetter - a Jet nanotuple producer from CMS MiniAOD

This project is a CMSSW module producing mostly flat tuples from run 2 data/MC

(A little bit under construction)

So far the code works on lxplus and CMSSW_8_0_6, and requires grid permissions to access data ```(voms-proxy-init -voms cms)```


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
git clone https://github.com/ekuruusu/Jetter/
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
