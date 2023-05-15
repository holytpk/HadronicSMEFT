'''
Compute efficiency of boosted selection in central powheg semilep ttbar.
'''

import os
import ROOT
from HadronicSMEFT.Samples.genTopJets_v1_pp import TT01j1l
from HadronicSMEFT.Tools.user                 import plot_directory
from HadronicSMEFT.Tools.delphesCutInterpreter import cutInterpreter
from RootTools.core.standard import *
import HadronicSMEFT.Tools.user as user
import Analysis.Tools.syncer as syncer

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--selection', action='store', default='singlelep-AK8pt500-AK8merged-njet4p-btag1p', help="Selection?")
args = argParser.parse_args()

# Logger
import HadronicSMEFT.Tools.logger as _logger
import RootTools.core.logger as _logger_rt
logger    = _logger.get_logger(   "INFO", logFile = None)
logger_rt = _logger_rt.get_logger("INFO", logFile = None)

directory = "/scratch-cbe/users/robert.schoefbeck/HadronicSMEFT/postprocessed/gen/v8/"

ttSemiLepInc = Sample.fromDirectory( "ttSemiLepInc", [os.path.join(directory, subdir) for subdir in ["TTToSemiLeptonic_UL16", "TTToSemiLeptonic_UL16APV", "TTToSemiLeptonic_UL17", "TTToSemiLeptonic_UL18"]]) 
#ttSemiLepInc = Sample.fromDirectory( "ttSemiLepInc", [os.path.join(directory, subdir) for subdir in ["TTToSemiLeptonic_UL16", "TTToSemiLeptonic_UL17", "TTToSemiLeptonic_UL18"]]) 
#ttSemiLepInc.reduceFiles(to=10)

#selection = "singlelep-AK8pt500-AK8merged-njet4p-btag1p" 
#selectionString = cutInterpreter.cutString("singlelep-AK8pt500-njet4p-btag1p") 
selectionString = cutInterpreter.cutString(args.selection) 


for var, name, texX in [ 
        ("parton_xi_nn+parton_xi_rr+parton_xi_kk", "D", "#xi_{nn}+#xi_{rr}+#xi_{kk}"), 
        ("parton_cos_phi", "cosPhi", "cos(#phi(l,l))"), 
        ]:

    h_inclusive = ttSemiLepInc.get1DHistoFromDraw( var, [20,-1,1], weightString = "(1.)")
    plot = Plot.fromHisto(args.selection+"_inclusive_"+name, histos = [[h_inclusive]], texX = texX)
    plotting.draw( plot, logY=False, plot_directory = os.path.join( user.plot_directory, "ttSemiLepInc" ) ) 

    eff_inclusive = ttSemiLepInc.get1DHistoFromDraw( var, [20,-1,1], selectionString = selectionString, weightString = "(1.)")
    eff_inclusive.style = styles.lineStyle( ROOT.kBlue, errors=True)
    eff_inclusive.Sumw2()

    plot = Plot.fromHisto(args.selection+"_"+name, histos = [[eff_inclusive]], texX = texX)
    plotting.draw( plot, logY=False, plot_directory = os.path.join( user.plot_directory, "ttSemiLepInc" ), copyIndexPHP = True) 

    eff_inclusive.Divide(h_inclusive) 

    plot = Plot.fromHisto(args.selection+"_efficiency_"+name, histos = [[eff_inclusive]], texX = texX)
    plotting.draw( plot, logY=False, plot_directory = os.path.join( user.plot_directory, "ttSemiLepInc" ), copyIndexPHP = True) 

    syncer.sync()
