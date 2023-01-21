#!/bin/sh

#python genPostProcessingFatJets.py --overwrite --addReweights --sample tt1LepHad #SPLIT100
python genPostProcessingFatJets.py --overwrite --removeDelphesFiles --addReweights --targetDir v6 --delphesEra ATLAS --sample tschRefPoint #SPLIT100
python genPostProcessingFatJets.py --overwrite --removeDelphesFiles --addReweights --targetDir v6 --delphesEra ATLAS --sample tschRefPointNoWidthRW #SPLIT100
