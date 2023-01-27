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
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/"
pp_dir       = "/scratch-cbe/users/robert.schoefbeck/HadronicSMEFT/postprocessed/gen/v8/"

TT01j = Sample.fromDirectory("TT01jDebug", texName = "TT01jDebug", directory = [os.path.join( pp_dir, "TT01jDebug" )]) 
TT01j.reweight_pkl = os.path.join(gridpack_directory, "TT01jDebug_reweight_card.pkl")
