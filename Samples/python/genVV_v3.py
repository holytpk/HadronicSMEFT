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
dbFile = os.path.join( cache_directory, 'sample_cache', 'genVV_v3.db')
overwrite = False

gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/"

WZto1L1Nu_HT300 = FWLiteSample.fromDAS("WZto1L1Nu_HT300", "/SMEFTNet_v3_WZto1L1Nu_HT300/schoef-SMEFTNet_v3_WZto1L1Nu_HT300-c02f0719055f37601ffced1868606ecc/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WZto1L1Nu_HT300.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto1L1NuNoRef_HT300_reweight_card.pkl"

WZto2L_HT300 = FWLiteSample.fromDAS("WZto2L_HT300", "/SMEFTNet_v3_WZto2L_HT300/schoef-SMEFTNet_v3_WZto2L_HT300-30e3aa1ed71123ab50b02a1bc9720d07/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WZto2L_HT300.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto2LNoRef_HT300_reweight_card.pkl"

WhadZlepJJEWK = FWLiteSample.fromDAS("WhadZlepJJEWK", "/SMEFTNet_v3_WhadZlepJJ/schoef-SMEFTNet_v3_WhadZlepJJ-5ce880e6626b2530231b9a4b214156df/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WhadZlepJJEWK.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WhadZlepJJEWKNoRef_reweight_card.pkl"
#
WlepZhadJJEWK = FWLiteSample.fromDAS("WlepZhadJJEWK", "/SMEFTNet_v3_WlepZhadJJ/schoef-SMEFTNet_v3_WlepZhadJJ-9808a07c6bf8f6faf90ec7b78fedd3eb/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WlepZhadJJEWK.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WlepZhadJJEWKNoRef_reweight_card.pkl"

WhadZlepJJ = FWLiteSample.fromDAS("WhadZlepJJ", "/SMEFTNet_v4_WhadZlepJJ/schoef-SMEFTNet_v4_WhadZlepJJ-5ce880e6626b2530231b9a4b214156df/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WhadZlepJJ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WhadZlepJJNoRef_reweight_card.pkl"
#
WlepZhadJJ = FWLiteSample.fromDAS("WlepZhadJJ", "/SMEFTNet_v4_WlepZhadJJ/schoef-SMEFTNet_v4_WlepZhadJJ-9808a07c6bf8f6faf90ec7b78fedd3eb/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WlepZhadJJ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WlepZhadJJNoRef_reweight_card.pkl"
