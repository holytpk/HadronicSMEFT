''' Benchmark samples for TopEFT (EDM)'''

# standard imports
import os
import ROOT

# RootTools
from RootTools.core.standard import *

# Logging
import logging
logger = logging.getLogger(__name__)

#TMB
from TMB.Samples.color import color

# Logging
if __name__ == "__main__":
    import TMB.Tools.logger as logger
    logger = logger.get_logger('INFO')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('INFO')

pp_dir       = "/scratch-cbe/users/robert.schoefbeck/HadronicSMEFT/predict_output/"
gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/"

TT01j1l_HT800_ctGIm = Sample.fromDirectory("TT01j1l_HT800_ctGIm", texName = "t#bar{t} 1l H_{T}>800", directory = [os.path.join( pp_dir, "ctGIm" )]) 
TT01j1l_HT800_ctGIm.reweight_pkl = os.path.join(gridpack_directory, "TT01jDebug_reweight_card.pkl")
