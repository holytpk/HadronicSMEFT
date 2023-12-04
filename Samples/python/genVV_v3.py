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

gridpack_directory = "/groups/hephy/cms/robert.schoefbeck/gridpacks/ParticleNet/"

#WZto1L1Nu_HT300 = FWLiteSample.fromDAS("WZto1L1Nu_HT300", "/SMEFTNet_v3_WZto1L1Nu_HT300/schoef-SMEFTNet_v3_WZto1L1Nu_HT300-c02f0719055f37601ffced1868606ecc/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
#WZto1L1Nu_HT300.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto1L1NuNoRef_HT300_reweight_card.pkl"
#
#WZto2L_HT300_Ref = FWLiteSample.fromDAS("WZto2L_HT300_Ref", "/SMEFTNet_v4_WZto2L_HT300/schoef-SMEFTNet_v4_WZto2L_HT300-b975dbcbae79b8cb5652f82f2f9f15c4/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
#WZto2L_HT300_Ref.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto2L_HT300_reweight_card.pkl"
#
#WZto2L_HT300_Ref_ext = FWLiteSample.fromDAS("WZto2L_HT300_Ref_ext", "/SMEFTNet_v4_ext_WZto2L_HT300/schoef-SMEFTNet_v4_ext_WZto2L_HT300-b975dbcbae79b8cb5652f82f2f9f15c4/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
#WZto2L_HT300_Ref_ext.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto2L_HT300_reweight_card.pkl"

WZto2L_HT300_v5 = FWLiteSample.fromDAS("WZto2L_HT300_v5", "/SMEFTNet_v5_WZto2L_HT300_Ref/schoef-SMEFTNet_v5_WZto2L_HT300_Ref-842e2448f1cf8999d7ad665e208a2e9f/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WZto2L_HT300_v5.reweight_pkl = "/groups/hephy/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto2L_HT300_Ref_v5_reweight_card.pkl"

WZto2L_HT300_NoRef_v5_v2 = FWLiteSample.fromDAS("WZto2L_HT300_NoRef_v5_v2", "/SMEFTNet_v5_v2_WZto2LNoRef_HT300/schoef-SMEFTNet_v5_v2_WZto2LNoRef_HT300-bc19e596fe8b10b51c22bf6b67fbb6cc/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
WZto2L_HT300_NoRef_v5_v2.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WZto2LNoRef_HT300_reweight_card.pkl"

#WhadZlepJJEWK = FWLiteSample.fromDAS("WhadZlepJJEWK", "/SMEFTNet_v3_WhadZlepJJ/schoef-SMEFTNet_v3_WhadZlepJJ-5ce880e6626b2530231b9a4b214156df/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
#WhadZlepJJEWK.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WhadZlepJJEWKNoRef_reweight_card.pkl"
##
#WlepZhadJJEWK = FWLiteSample.fromDAS("WlepZhadJJEWK", "/SMEFTNet_v3_WlepZhadJJ/schoef-SMEFTNet_v3_WlepZhadJJ-9808a07c6bf8f6faf90ec7b78fedd3eb/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
#WlepZhadJJEWK.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WlepZhadJJEWKNoRef_reweight_card.pkl"
#
#WhadZlepJJ = FWLiteSample.fromDAS("WhadZlepJJ", "/SMEFTNet_v4_WhadZlepJJ/schoef-SMEFTNet_v4_WhadZlepJJ-5ce880e6626b2530231b9a4b214156df/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
#WhadZlepJJ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WhadZlepJJNoRef_reweight_card.pkl"
##
#WlepZhadJJ = FWLiteSample.fromDAS("WlepZhadJJ", "/SMEFTNet_v4_WlepZhadJJ/schoef-SMEFTNet_v4_WlepZhadJJ-9808a07c6bf8f6faf90ec7b78fedd3eb/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
#WlepZhadJJ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WlepZhadJJNoRef_reweight_card.pkl"
#
#WG_HT300 = FWLiteSample.fromDAS("WG_HT300", "/SMEFTNet_v4_WG_HT300/schoef-SMEFTNet_v4_WG_HT300-6a2c92f8213d1f2251274a9b114b51d3/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
#WG_HT300.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WG_HT300_reweight_card.pkl"
#
#WJets_HT300 = FWLiteSample.fromDAS("WJets_HT300", "/SMEFTNet_v4_WJets_HT300/schoef-SMEFTNet_v4_WJets_HT300-e979169dcc7eede63393cf75a7ef663a/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
#DY_HT300 = FWLiteSample.fromDAS("DY_HT300", "/SMEFTNet_v4_DY_HT300/schoef-SMEFTNet_v4_DY_HT300-3d08d5678ff01e219241ce9532334619/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/", overwrite=overwrite)
