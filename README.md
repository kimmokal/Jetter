# Jetter - A jet nanotuple producer from CMS MiniAOD data

This project is a CMSSW module producing mostly flat tuples from run 2 data/MC, to be used in DNN training.

## Setting up

First setup your GitHub or use these dummy values:
```
git config --global user.name 'Your Name'
git config --global user.email 'your@ema.il'
git config --global user.github 'username'
```
Make sure you register your SSH key in github: https://help.github.com/articles/generating-ssh-keys

As the code works on lxplus and CMSSW_8_0_6.

Creating the work area
```
ssh -Y CERNusername@lxplus006.cern.ch

cmsrel CMSSW_8_0_6
cd CMSSW_8_0_6/src
cmsenv
git clone https://github.com/ekuruusu/Jetter/
scram b
cd Jetter/MiniAnalyzer/python
```

## Running the code

Change the source file in python/ConfFile_cfg.py 

To run the code and generate the tuples:
```
    cmsRun ConfFile_cfg.py
```

To view the tuples in ROOT:
```
    root -l nanotuple_*
    TBrowser a
```
