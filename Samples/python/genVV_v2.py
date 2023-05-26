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
dbFile = os.path.join( cache_directory, 'sample_cache', 'genVV_v2.db')
overwrite = False

# for flavor analysis 

gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/"

WZto1L1Nu_HT300 = FWLiteSample.fromDAS("WZto1L1Nu_HT300", "/SMEFTNet_WZto1L1Nu_HT300/schoef-SMEFTNet_WZto1L1Nu_HT300-56aff94bea8f52b00c95191f2684039c/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WZto1L1Nu_HT300.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto1L1Nu_HT800_reweight_card.pkl"

WZto2L_HT300 = FWLiteSample.fromDAS("WZto2L_HT300", "/SMEFTNet_WZto2L_HT300/schoef-SMEFTNet_WZto2L_HT300-b975dbcbae79b8cb5652f82f2f9f15c4/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WZto2L_HT300.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto2L_HT300_reweight_card.pkl"

WhadZlepJJ = FWLiteSample.fromDAS("WhadZlepJJ", "/SMEFTNet_WhadZlepJJ/schoef-SMEFTNet_WhadZlepJJ-d567ccc03c3872ad4f96d95d59dd48f0/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WhadZlepJJ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WhadZlepJJ_EWK_LO_SM_mjj100_pTj10_reweight_card.pkl"

WlepZhadJJ = FWLiteSample.fromDAS("WlepZhadJJ", "/SMEFTNet_WlepZhadJJ/schoef-SMEFTNet_WlepZhadJJ-78d452a89a799a4e9828293ffc399337/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WlepZhadJJ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WlepZhadJJ_EWK_LO_SM_mjj100_pTj10_reweight_card.pkl"
