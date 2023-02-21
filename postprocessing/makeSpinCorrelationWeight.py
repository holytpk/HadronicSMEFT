from HadronicSMEFT.Samples.genTopJets_v1_pp import TT01j1l
from HadronicSMEFT.Tools.user                 import plot_directory
from HadronicSMEFT.Tools.delphesCutInterpreter import cutInterpreter

#from Samples.nanoAOD.
ttInc = FWLiteSample.fromDas/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM

selection = cutInterpreter.cutString("singlelep-AK8pt500-AK8merged-njet4p-btag1p") 
# Measure efficiency of HT cut

