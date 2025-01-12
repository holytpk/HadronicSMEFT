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
dbFile = os.path.join( cache_directory, 'sample_cache', 'genTopJets_v1.db')
overwrite = False

# for flavor analysis 

gridpack_directory = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/"

tt1LepHad_test = FWLiteSample.fromFiles("tt1LepHad_test", ["/eos/vbc/group/cms/robert.schoefbeck/ParticleNet/test/GEN_LO_0j_102X.root"])
tt1LepHad_test.reweight_pkl = os.path.join(gridpack_directory, "tt01j-1l-NPtHad_HT800_slc7_amd64_gcc700_CMSSW_10_6_19_tarball.pkl")

TT01jDebug = FWLiteSample.fromDAS("TT01jDebug", "/PNet/schoef-PNet-af3a260270385fd3f7df95e2de8dce51/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/")
TT01jDebug.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/TT01jDebug_reweight_card.pkl" 

TT01j_HT800_ext = FWLiteSample.fromDAS("TT01j_HT800_ext", "/PNet_TT01j1l_HT800-ext/schoef-PNet_TT01j1l_HT800-ext-511efaae0884df23dd8cb8291ee53fc9/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/")
TT01j_HT800_ext.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/TT01jDebug_reweight_card.pkl" 

TT01j1l = FWLiteSample.fromDAS("TT01j1l", "/PNet_TT01j1l/schoef-PNet_TT01j1l-c64478b79b2fe2c66cf5b1ff13d188f0/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/")
TT01j1l.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/TT01jDebug_reweight_card.pkl" 
#
TT01j2l       = FWLiteSample.fromDAS("TT01j2l", "/PNet_TT01j2l/schoef-PNet_TT01j2l-a0241147cf6c0844d7923af96307ec5a/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/")
TT01j2l.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/TT01jDebug_reweight_card.pkl" 

TT01j2l_HT800 = FWLiteSample.fromDAS("TT01j2l_HT800", "/PNet_TT01j2l_HT800/schoef-PNet_TT01j2l_HT800-1e05253e09d02ca0e67caa8c74dd84b2/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/")
TT01j2l_HT800.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/TT01jDebug_reweight_card.pkl" 

WhadZlepJJ = FWLiteSample.fromDAS("WhadZlepJJ", "/PNet_WhadZlepJJ/schoef-PNet_WhadZlepJJ-d567ccc03c3872ad4f96d95d59dd48f0/USER", dbFile=dbFile, instance="phys03", prefix="root://eos.grid.vbc.ac.at/")
WhadZlepJJ.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/WhadZlepJJ_EWK_LO_SM_mjj100_pTj10_reweight_card.pkl" 

#tmp = FWLiteSample.fromFiles("GEN", ["/users/robert.schoefbeck/CMS/CMSSW_10_6_27/src/Samples/crab/gen/GEN_LO_0j_102X.root"]) 
#tmp.reweight_pkl = "/eos/vbc/group/cms/robert.schoefbeck/gridpacks/ParticleNet/TT01jDebug_reweight_card.pkl" 

tsch_test = FWLiteSample.fromFiles("tsch_test", ["/eos/vbc/group/cms/robert.schoefbeck/ParticleNet/test-tsch/GEN_LO_0j_102X.root"])
tsch_test.reweight_pkl = os.path.join(gridpack_directory, "t-sch_reweight_card.pkl")

tschRefPoint_test = FWLiteSample.fromFiles("tschRefPoint_test", ["/eos/vbc/group/cms/robert.schoefbeck/ParticleNet/test-tsch-RefPoint/GEN_LO_0j_102X.root"])
tschRefPoint_test.reweight_pkl = os.path.join(gridpack_directory, "t-sch-RefPoint_reweight_card.pkl")

#tschRefPointNoWidthRW = FWLiteSample.fromDAS("tschRefPointNoWidthRW", "/t-sch-RefPoint-noWidthRW/schoef-t-sch-RefPoint-noWidthRW-ad4a531db5c6d25218664ba7bdd18ceb/USER", "phys03", dbFile = dbFile, overwrite=overwrite, skipCheck=True)
#tschRefPointNoWidthRW.reweight_pkl = os.path.join(gridpack_directory, "t-sch-RefPoint-noWidthRW_reweight_card.pkl")

tschRefPoint = FWLiteSample.fromDirectory("tschRefPoint", "/groups/hephy/cms/robert.schoefbeck/ParticleNet/GEN/t-sch-RefPoint/")
tschRefPoint.reweight_pkl = "/groups/hephy/cms/robert.schoefbeck/ParticleNet/GEN/t-sch-RefPoint/t-sch-RefPoint_reweight_card.pkl"

tschRefPointNoWidthRW = FWLiteSample.fromDirectory("tschRefPointNoWidthRW", "/groups/hephy/cms/robert.schoefbeck/ParticleNet/GEN/t-sch-RefPoint-noWidthRW/")
tschRefPointNoWidthRW.reweight_pkl = "/groups/hephy/cms/robert.schoefbeck/ParticleNet/GEN/t-sch-RefPoint-noWidthRW/t-sch-RefPoint-noWidthRW_reweight_card.pkl"
