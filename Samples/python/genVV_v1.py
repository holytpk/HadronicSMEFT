''' GEN samples for HadronicSMEFT'''

# standard imports
import os

# RootTools
from RootTools.core.standard import *

# Logging
if __name__ == "__main__":
    import HadronicSMEFT.Tools.logger as logger
    logger = logger.get_logger('DEBUG')
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger('DEBUG')

# HadronicSMEFT
from HadronicSMEFT.Tools.user import cache_directory, gridpack_directory 

# sqlite3 sample cache file
dbFile = os.path.join( cache_directory, 'sample_cache', 'genVV_v1.db')
overwrite = False

# for flavor analysis 

gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/"

WhadZlepJJ = FWLiteSample.fromDAS("WhadZlepJJ", "/PNet_WhadZlepJJ/schoef-PNet_WhadZlepJJ-d567ccc03c3872ad4f96d95d59dd48f0/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/")
WhadZlepJJ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WhadZlepJJ_EWK_LO_SM_mjj100_pTj10_reweight_card.pkl"

WZto1L1Nu_HT800 = FWLiteSample.fromDAS("WZto1L1Nu_HT800", "/PNet_WZto1L1Nu_HT800/schoef-PNet_WZto1L1Nu_HT800-6f9b96e49dc64ec77fff158a39815098/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/")
WZto1L1Nu_HT800.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto1L1Nu_reweight_card.pkl"
