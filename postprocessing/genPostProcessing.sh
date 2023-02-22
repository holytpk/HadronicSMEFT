
python genPostProcessing.py --addReweights --targetDir v8 --trainingCoefficients ctGRe ctGIm ctWRe ctWIm ctBRe ctBIm cHtbRe cHtbIm cHt --sample TT01jDebug --process ttbar --delphesEra ATLAS --removeDelphesFiles #SPLIT199 
python genPostProcessing.py --addReweights --targetDir v8 --trainingCoefficients ctGRe ctGIm ctWRe ctWIm ctBRe ctBIm cHtbRe cHtbIm cHt --sample TT01j1l    --process ttbar --delphesEra ATLAS --removeDelphesFiles #SPLIT200 

python genPostProcessing.py --addReweights --targetDir v8 --trainingCoefficients cHDD cHbox cW cWtil cHW cHWtil cHWB cHB cHWBtil --sample WhadZlepJJ --process WhadZlep --delphesEra ATLAS --removeDelphesFiles #SPLIT200 

#python genPostProcessing.py --addReweights --targetDir v8 --sample TT01j2l --delphesEra ATLAS --removeDelphesFiles #SPLIT200 
#python genPostProcessing.py --addReweights --targetDir v8 --sample TT01j2l_HT800 --delphesEra ATLAS --removeDelphesFiles #SPLIT200 
