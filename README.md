# HadronicSMEFT
```
cmsrel CMSSW_10_6_28
cd CMSSW_10_6_28/src
cmsenv
git cms-init
git clone https://github.com/HephyAnalysisSW/HadronicSMEFT
git clone --branch UL https://github.com/HephyAnalysisSW/Analysis
git clone https://github.com/HephyAnalysisSW/Samples
git clone https://github.com/HephyAnalysisSW/RootTools
git clone https://github.com/HephyAnalysisSW/NanoAODJMARTools.git PhysicsTools/NanoAODJMARTools
cp $CMSSW_BASE/src/PhysicsTools/NanoAODJMARTools/xmlfiles/* $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/
scram setup fastjet
scram setup fastjet-contrib
scram b -j40
```

## Delphes

```
cd $CMSSW_BASE/..
git clone https://github.com/TTXPheno/delphes.git
#patch $CMSSW_BASE/../delphes/cards/delphes_card_CMS.tcl < $CMSSW_BASE/src/TTXPheno/patches/slim_delphes.diff # Reduce Delphes output
cd delphes
./configure
sed -i -e 's/c++0x/c++17/g' Makefile
make -j 4 
```
