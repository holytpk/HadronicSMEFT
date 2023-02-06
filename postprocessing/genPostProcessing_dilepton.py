#!/usr/bin/env python
'''  produce fat jet ntuple 
'''
#
# Standard imports and batch mode
#
import os, sys, imp, uuid
import copy, shutil
import ROOT
ROOT.gROOT.SetBatch(True)
from math                             import sqrt, cos, sin, pi, acos, cosh, sinh
from numpy import sign
#RootTools
from RootTools.core.standard          import *

#Analysis
from Analysis.Tools.WeightInfo        import WeightInfo
from Analysis.Tools.HyperPoly         import HyperPoly
from Analysis.Tools.GenSearch         import GenSearch
from Analysis.Tools.helpers           import deltaPhi, deltaR, deltaR2, checkRootFile
from Analysis.Tools.DelphesProducer   import DelphesProducer

# HadronicSMEFT
import HadronicSMEFT.Tools.user                 as user
from HadronicSMEFT.Tools.genObjectSelection     import genJetId
from HadronicSMEFT.Tools.DelphesObjectSelection import isGoodRecoLepton, isGoodRecoJet, isGoodRecoPhoton
#import HadronicSMEFT.Tools.fixTVecMul

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--overwrite',          action='store',      nargs='?', choices = ['none', 'all', 'target'], default = 'none', help='Overwrite?')#, default = True)
argParser.add_argument('--targetDir',          action='store',      default='v1')
argParser.add_argument('--sample',             action='store',      default='tt1LepHad', help="Name of the sample loaded from fwlite_benchmarks. Only if no inputFiles are specified")
argParser.add_argument('--inputFiles',         action='store',      nargs = '*', default=[])
argParser.add_argument('--delphesEra',         action='store',      default = "RunII", choices = ["RunII", "ATLAS", "RunIICentral", "RunIInoDelphesIso", "RunIIPileUp", "PhaseII"], help="specify delphes era")
argParser.add_argument('--addReweights',       action='store_true',   help="Add reweights?")
argParser.add_argument('--nJobs',              action='store',      nargs='?', type=int, default=1,  help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',                action='store',      nargs='?', type=int, default=0,  help="Run only job i")
argParser.add_argument('--removeDelphesFiles', action='store_true',   help="remove Delphes file after postprocessing?")
argParser.add_argument('--interpolationOrder', action='store',      nargs='?', type=int, default=2,  help="Interpolation order for EFT weights.")
argParser.add_argument('--trainingCoefficients', action='store',    nargs='?', default=['ctWRe', 'ctWIm'],  help="Training vectors for particle net")
args = argParser.parse_args()

#
# Logger
#
import HadronicSMEFT.Tools.logger as _logger
import RootTools.core.logger as _logger_rt
logger    = _logger.get_logger(   args.logLevel, logFile = None)
logger_rt = _logger_rt.get_logger(args.logLevel, logFile = None)

# Load sample either from 
if len(args.inputFiles)>0:
    logger.info( "Input files found. Ignoring 'sample' argument. Files: %r", args.inputFiles)
    sample = FWLiteSample( args.targetSampleName, args.inputFiles)
else:
    sample_file = "$CMSSW_BASE/python/HadronicSMEFT/Samples/genTopJets_v1.py"
    samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
    sample = getattr( samples, args.sample )
    logger.debug( 'Loaded sample %s with %i files.', sample.name, len(sample.files) )

maxEvents = -1
if args.small: 
    args.targetDir += "_small"
    maxEvents       = 100 
    sample.files=sample.files[:1]

# CMSSW FastJet & CMSSW wrappers
# The FastJet-contrib wrappers (ECF and Nsubjettiness) are in https://github.com/HephyAnalysisSW/NanoAODJMARTools.git
import fastjet
ak8 = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.8, fastjet.E_scheme)
softDrop       = ROOT.SoftDropWrapper(0.0, 0.1, 0.8, 200)
nSubjettiness  = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6 )

# Energy correlators
ecf  = ROOT.ECFWrapper()
ecfs = [
 ('ecf1',       ( 1, 1., 1., "ECF" )),
 ('ecf2',       ( 2, 1., 1., "ECF" )),
 ('ecf3',       ( 3, 1., 1., "ECF" )),
 ('ecfC1',      ( 1, 1., 1., "C" )),
 ('ecfC2',      ( 2, 1., 1., "C" )),
 ('ecfC3',      ( 3, 1., 1., "C" )),
 ('ecfD',       ( 2, 1., 1., "D" )),
 ('ecfDbeta2',  ( 2, 2., 2., "D" )),
 ('ecfM1',      ( 1, 1., 1., "M" )),
 ('ecfM2',      ( 2, 1., 1., "M" )),
 ('ecfM3',      ( 3, 1., 1., "M" )),
 ('ecfM1beta2', ( 1, 2., 2., "M" )),
 ('ecfM2beta2', ( 2, 2., 2., "M" )),
 ('ecfM3beta2', ( 3, 2., 2., "M" )),
 ('ecfN1',      ( 1, 1., 1., "N" )),
 ('ecfN2',      ( 2, 1., 1., "N" )),
 ('ecfN3',      ( 3, 1., 1., "N" )),
 ('ecfN1beta2', ( 1, 2., 2., "N" )),
 ('ecfN2beta2', ( 2, 2., 2., "N" )),
 ('ecfN3beta2', ( 3, 2., 2., "N" )),
 ('ecfU1',      ( 1, 1., 1., "U" )),
 ('ecfU2',      ( 2, 1., 1., "U" )),
 ('ecfU3',      ( 3, 1., 1., "U" )),
 ('ecfU1beta2', ( 1, 2., 2., "U" )),
 ('ecfU2beta2', ( 2, 2., 2., "U" )),
 ('ecfU3beta2', ( 3, 2., 2., "U" )),
]

for _, args_ in ecfs:
    ecf.addECF( *args_ )

# variables (to be stored)
variables = []

# Load reweight pickle file if supposed to keep weights. 
if args.addReweights:
    weightInfo = WeightInfo( sample.reweight_pkl )
    weightInfo.set_order( args.interpolationOrder ) 
    # Determine coefficients for storing in vector
    # Sort Ids wrt to their position in the card file

    # weights from base base points 
    weight_base      = TreeVariable.fromString( "weight[base/F]")
    weight_base.nMax = weightInfo.nid
    variables.append(weight_base)

    # coefficients for the weight parametrization
    param_vector      = TreeVariable.fromString( "p[C/F]" )
    param_vector.nMax = HyperPoly.get_ndof(weightInfo.nvar, args.interpolationOrder)
    hyperPoly         = HyperPoly( args.interpolationOrder )
    variables.append(param_vector)
    variables.append(TreeVariable.fromString( "chi2_ndof/F"))
    def interpret_weight(weight_id):
        str_s = weight_id.rstrip('_nlo').split('_')
        res={}
        for i in range(len(str_s)/2):
            res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
        return res

    # Suddenly only lower case weight.id ... who on earth does such things?
    weightInfo_data_lower = {k.lower():val for k, val in weightInfo.data.iteritems()}
    weightInfo_data_lower.update(weightInfo.data)

    target_coeff = {}
    for i_comb, comb in enumerate(weightInfo.combinations[1:]):
        variables.append( "target_%s/F"%("_".join(comb)) )
        target_coeff[i_comb] = "_".join(comb)

# output directory
output_directory = os.path.join(user.skim_output_directory, 'gen', args.targetDir, sample.name) 

if not os.path.exists( output_directory ): 
    try:
        os.makedirs( output_directory )
    except OSError:
        pass
    logger.info( "Created output directory %s", output_directory )

# Run only job number "args.job" from total of "args.nJobs"
if args.nJobs>1:
    n_files_before = len(sample.files)
    sample = sample.split(args.nJobs)[args.job]
    n_files_after  = len(sample.files)
    logger.info( "Running job %i/%i over %i files from a total of %i.", args.job, args.nJobs, n_files_after, n_files_before)

# tmp_output_directory
tmp_output_directory  = os.path.join( user.tmp_output_directory, str(uuid.uuid3(uuid.NAMESPACE_OID, sample.name)) )
try:    #Avoid trouble with race conditions in multithreading
    os.makedirs(tmp_output_directory)
    logger.info( "Created output directory %s.", tmp_output_directory )
except:
    pass

# output file & log files
output_filename =  os.path.join(output_directory, sample.name + '.root')
_logger.   add_fileHandler( output_filename.replace('.root', '.log'), args.logLevel )
_logger_rt.add_fileHandler( output_filename.replace('.root', '_rt.log'), args.logLevel )

# small helpers
def makeP4(cand):
     v = ROOT.TLorentzVector()
     v.SetPtEtaPhiM(cand.pt(),cand.eta(),cand.phi(),cand.mass())

     return v

def varnames( vec_vars ):
    return [v.split('/')[0] for v in vec_vars.split(',')]

def addIndex( collection ):
    for i  in range(len(collection)):
        collection[i]['index'] = i

def vecSumPt(*args):
    return sqrt( sum([o['pt']*cos(o['phi']) for o in args],0.)**2 + sum([o['pt']*sin(o['phi']) for o in args],0.)**2 )

def fill_vector_collection( event, collection_name, collection_varnames, objects):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects):
        for var in collection_varnames:
            getattr(event, collection_name+"_"+var)[i_obj] = obj[var]

def fill_vector( event, collection_name, collection_varnames, obj):
    for var in collection_varnames:
        try:
            setattr(event, collection_name+"_"+var, obj[var] )
        except TypeError as e:
            logger.error( "collection_name %s var %s obj[var] %r", collection_name, var,  obj[var] )
            raise e
        except KeyError as e:
            logger.error( "collection_name %s var %s obj[var] %r", collection_name, var,  obj[var] )
            raise e


def get_spincorrelation_angles(antilepton_4v, lepton_4v, topQ_4v, atopQ_4v):

    # anti-lepton is +ve charged lepton coming from top quark decay chain 
    # lepton is -ve charged lepton coming from anti-top quark decay chain

    p4_top = copy.deepcopy(topQ_4v)
    p4_antitop = copy.deepcopy(atopQ_4v)
    p4_l_plus = copy.deepcopy(antilepton_4v)
    p4_l_minus = copy.deepcopy(lepton_4v)

    # Boosting into tt rest frame

    boost_tt = (p4_top + p4_antitop).BoostVector()

    p4_top.Boost(-boost_tt)
    p4_antitop.Boost(-boost_tt)

    p4_l_plus.Boost(-boost_tt)
    p4_l_minus.Boost(-boost_tt)

    p4_l_plus.Boost(-p4_top.BoostVector())
    p4_l_minus.Boost(-p4_antitop.BoostVector())

    # Unit vectors for leptons 

    l_plus = p4_l_plus.Vect().Unit()
    l_minus = p4_l_minus.Vect().Unit()

    # Create k,r,n axes

    k_hat = p4_top.Vect().Unit()
    p_hat = ROOT.TVector3(0,0,1) # Eq. 4.13 Bernreuther!
    
    y = k_hat.Dot(p_hat)
    r = sqrt( 1. - (y*y) )
    sign_ =  float(sign(y)) 
        
    r_hat = sign_*(1./r)*(p_hat - (y*k_hat))
    n_hat = sign_*(1./r)*(p_hat.Cross(k_hat))

    # Bernreuther Table 5  
    n_a = sign_*n_hat
    n_b = -sign_*n_hat
    r_a = sign_*r_hat
    r_b = -sign_*r_hat
    k_a =  k_hat
    k_b = -k_hat

    # Eq 4.21 Bernreuther (k* and r*)
    sign_star = float(sign(abs(p4_top.Rapidity()) - abs(p4_antitop.Rapidity())))
    k_a_star  = sign_star*k_hat
    k_b_star  = -sign_star*k_hat
    r_a_star  = sign_star*sign_*r_hat
    r_b_star  = -sign_star*sign_*r_hat

    # Determine cos(angles)
    # Bernreuther Eq. 4.7

    cosThetaPlus_n  = n_a.Dot(l_plus)
    cosThetaMinus_n = n_b.Dot(l_minus)
    cosThetaPlus_r  = r_a.Dot(l_plus)
    cosThetaMinus_r = r_b.Dot(l_minus)
    cosThetaPlus_k  = k_a.Dot(l_plus)
    cosThetaMinus_k = k_b.Dot(l_minus)

    cosThetaPlus_r_star  = r_a_star.Dot(l_plus)
    cosThetaMinus_r_star = r_b_star.Dot(l_minus)
    cosThetaPlus_k_star  = k_a_star.Dot(l_plus)
    cosThetaMinus_k_star = k_b_star.Dot(l_minus)

    cos_phi     = l_plus.Dot(l_minus)

    return (cosThetaPlus_n, cosThetaMinus_n, cosThetaPlus_r, cosThetaMinus_r, cosThetaPlus_k, cosThetaMinus_k, cosThetaPlus_r_star, cosThetaMinus_r_star, cosThetaPlus_k_star, cosThetaMinus_k_star, cos_phi)

# EDM standard variables
variables  += ["run/I", "lumi/I", "evt/l"]

# lepton vector 
lep_vars       =  "pt/F,eta/F,phi/F,pdgId/I,status/I"
lep_extra_vars =  "mother_pdgId/I,grandmother_pdgId/I"
lep_varnames   =  varnames( lep_vars )
lep_all_varnames = lep_varnames + varnames(lep_extra_vars)
variables     += ["genLep[%s]"%(','.join([lep_vars, lep_extra_vars]))]

# parton top quarks
#variables += ["parton_hadTop_pt/F", "parton_hadTop_eta/F", "parton_hadTop_phi/F", "parton_hadTop_mass/F", "parton_hadTop_pdgId/I", "parton_hadTop_decayAngle_theta/F", "parton_hadTop_decayAngle_phi/F"]
#variables += ["parton_hadTop_q1_pt/F",  "parton_hadTop_q1_eta/F",  "parton_hadTop_q1_phi/F",  "parton_hadTop_q1_mass/F",  "parton_hadTop_q1_pdgId/I"]
#variables += ["parton_hadTop_q2_pt/F",  "parton_hadTop_q2_eta/F",  "parton_hadTop_q2_phi/F",  "parton_hadTop_q2_mass/F",  "parton_hadTop_q2_pdgId/I"]
#variables += ["parton_hadTop_b_pt/F",   "parton_hadTop_b_eta/F",   "parton_hadTop_b_phi/F",   "parton_hadTop_b_mass/F",   "parton_hadTop_b_pdgId/I"]
#variables += ["parton_hadTop_W_pt/F",   "parton_hadTop_W_eta/F",   "parton_hadTop_W_phi/F",   "parton_hadTop_W_mass/F", "parton_hadTop_W_pdgId/I"]

variables += ["parton_lepTop1_pt/F",   "parton_lepTop1_eta/F",   "parton_lepTop1_phi/F",   "parton_lepTop1_mass/F",   "parton_lepTop1_pdgId/I"]
variables += ["parton_lepTop1_lep_pt/F",      "parton_lepTop1_lep_eta/F",      "parton_lepTop1_lep_phi/F",      "parton_lepTop1_lep_mass/F",      "parton_lepTop1_lep_pdgId/I"]
variables += ["parton_lepTop1_nu_pt/F",       "parton_lepTop1_nu_eta/F",       "parton_lepTop1_nu_phi/F",       "parton_lepTop1_nu_mass/F",       "parton_lepTop1_nu_pdgId/I"]
variables += ["parton_lepTop1_b_pt/F", "parton_lepTop1_b_eta/F", "parton_lepTop1_b_phi/F", "parton_lepTop1_b_mass/F", "parton_lepTop1_b_pdgId/I"]
variables += ["parton_lepTop1_W_pt/F", "parton_lepTop1_W_eta/F", "parton_lepTop1_W_phi/F", "parton_lepTop1_W_mass/F", "parton_lepTop1_W_pdgId/I"]

variables += ["parton_lepTop2_pt/F",   "parton_lepTop2_eta/F",   "parton_lepTop2_phi/F",   "parton_lepTop2_mass/F",   "parton_lepTop2_pdgId/I"]
variables += ["parton_lepTop2_lep_pt/F",      "parton_lepTop2_lep_eta/F",      "parton_lepTop2_lep_phi/F",      "parton_lepTop2_lep_mass/F",      "parton_lepTop2_lep_pdgId/I"]
variables += ["parton_lepTop2_nu_pt/F",       "parton_lepTop2_nu_eta/F",       "parton_lepTop2_nu_phi/F",       "parton_lepTop2_nu_mass/F",       "parton_lepTop2_nu_pdgId/I"]
variables += ["parton_lepTop2_b_pt/F", "parton_lepTop2_b_eta/F", "parton_lepTop2_b_phi/F", "parton_lepTop2_b_mass/F", "parton_lepTop2_b_pdgId/I"]
variables += ["parton_lepTop2_W_pt/F", "parton_lepTop2_W_eta/F", "parton_lepTop2_W_phi/F", "parton_lepTop2_W_mass/F", "parton_lepTop2_W_pdgId/I"]

variables += ["parton_cosThetaPlus_n/F", "parton_cosThetaMinus_n/F", "parton_cosThetaPlus_r/F", "parton_cosThetaMinus_r/F", "parton_cosThetaPlus_k/F", "parton_cosThetaMinus_k/F", 
              "parton_cosThetaPlus_r_star/F", "parton_cosThetaMinus_r_star/F", "parton_cosThetaPlus_k_star/F", "parton_cosThetaMinus_k_star/F",
              "parton_xi_nn/F", "parton_xi_rr/F", "parton_xi_kk/F", 
              "parton_xi_nr_plus/F", "parton_xi_nr_minus/F", "parton_xi_rk_plus/F", "parton_xi_rk_minus/F", "parton_xi_nk_plus/F", "parton_xi_nk_minus/F", 
              "parton_cos_phi/F", "parton_cos_phi_lab/F", "parton_abs_delta_phi_ll_lab/F",
             ] 

variables += ["delphesJet_pt/F", "delphesJet_eta/F", "delphesJet_phi/F", "delphesJet_mass/F", "delphesJet_nConstituents/I"] 
variables += ['delphesJet_SDmass/F', 
              'delphesJet_SDsubjet0_eta/F', 'delphesJet_SDsubjet0_deltaEta/F', 'delphesJet_SDsubjet0_phi/F', 'delphesJet_SDsubjet0_deltaPhi/F', 'delphesJet_SDsubjet0_deltaR/F', 'delphesJet_SDsubjet0_mass/F', 
              'delphesJet_SDsubjet1_eta/F', 'delphesJet_SDsubjet1_deltaEta/F', 'delphesJet_SDsubjet1_phi/F', 'delphesJet_SDsubjet1_deltaPhi/F', 'delphesJet_SDsubjet1_deltaR/F', 'delphesJet_SDsubjet1_mass/F', 
              'delphesJet_tau1/F', 'delphesJet_tau2/F', 'delphesJet_tau3/F', 'delphesJet_tau4/F', 'delphesJet_tau21/F', 'delphesJet_tau32/F']
'''
for i_ecf, (name, _) in enumerate( ecfs ):
    variables.append( "delphesJet_%s/F"%name )

variables += ["delphesJet_dR_matched_hadTop_parton/F", "delphesJet_dR_lepTop_parton/F", "delphesJet_dR_hadTop_q1/F", "delphesJet_dR_hadTop_q2/F", 
              "delphesJet_dR_hadTop_W/F",              "delphesJet_dR_hadTop_b/F", "delphesJet_dR_hadTop_maxq1q2b/F"]
'''
if args.addReweights:
    # for each Wilson coefficient listed in args.trainingCoefficients, store a separate length-3 ntuple of ('w0'*10**6, 'w1', 'w2') to facilitate particle-net training 
    for coefficient in args.trainingCoefficients:    
        variables += [VectorTreeVariable.fromString("%s[coeff/F]"%coefficient, nMax=3 )]

# enumeration according to PF type: https://github.com/cms-sw/cmssw/blob/master/DataFormats/ParticleFlowCandidate/interface/PFCandidate.h#L44-L52
categories = [
    {'name':'e',   'type':2, 'eflow_func':lambda p:abs(p['pdgId'])==11, 'gen_func':lambda p:abs(p['pdgId'])==11}, #electrons 
    {'name':'mu',  'type':3, 'eflow_func':lambda p:abs(p['pdgId'])==13, 'gen_func':lambda p:abs(p['pdgId'])==13}, #muons 
    {'name':'ph',  'type':4, 'eflow_func':lambda p:abs(p['pdgId'])==22, 'gen_func':lambda p:p['pdgId']==22}, #photons 
    {'name':'chh', 'type':1, 'eflow_func':lambda p:abs(p['pdgId'])>100 and p['charge']!=0, 'gen_func':lambda p:abs(p['pdgId'])>100 and p['charge']!=0 }, #charged hadrons 
    {'name':'neh', 'type':5, 'eflow_func':lambda p:abs(p['pdgId'])>100 and p['charge']==0, 'gen_func':lambda p:abs(p['pdgId'])>100 and p['charge']==0 }, # neutral hadrons
]

#cand_vars           =  "pt/F,etarel/F,phirel/F,eta/F,phi/F,pdgId/I"
#cand_varnames       =  varnames( cand_vars )
cand_vars        =  "pt/F,etarel/F,phirel/F,eta/F,phi/F,pdgId/I,charge/I,type/I"
cand_varnames    =  varnames( cand_vars )

nCandMax = 200
variables.append(VectorTreeVariable.fromString("eflow[%s]"%(cand_vars), nMax=nCandMax ))

# reconstructed leptons
recoLep_vars       = "pt/F,eta/F,phi/F,pdgId/I,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F,genMatched/I,MT/F"
variables         += ["recoLep[%s]"%recoLep_vars]
recoLep_varnames   = varnames( recoLep_vars )

# reconstructed jets
#btagWPs = ["loose"]#, "medium", "tight"] #, "looswMTD", "mediumMTD", "tightMTD"]
#default_btagWP = "loose"
variables.append( "nBTag/I" )
recoJet_vars    = 'pt/F,eta/F,phi/F,bTag/I,nCharged/I,nNeutrals/I'#,matchGenBJet/I'#,pt_JEC_up/F,pt_JEC_up/F'

variables += ["recoJet[%s]"%recoJet_vars]
recoJet_varnames = varnames( recoJet_vars )
variables += ["recoBj0_%s"%var for var in recoJet_vars.split(',')]
variables += ["recoBj1_%s"%var for var in recoJet_vars.split(',')]

variables += ["recoMet_pt/F", "recoMet_phi/F"]
variables += ["delphesGenMet_pt/F", "delphesGenMet_phi/F"]

logger.info( "Running over files: %s", ", ".join(sample.files ) )

readers = []

# FWLite reader 
products = {
    'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
    'gen':{'type':'GenEventInfoProduct', 'label':'generator'},
    'gp':{'type':'vector<reco::GenParticle>', 'label':("genParticles")},
    'ak8GenJets':{'type':'vector<reco::GenJet>', 'label':("ak8GenJetsNoNu")},
}

fwliteReader = sample.fwliteReader( products = products )
readers.append( fwliteReader )

# Delphes reader if we run Delphes
from HadronicSMEFT.Tools.DelphesReaderEFlow     import DelphesReader
if args.delphesEra == 'RunII':
    delphesCard = 'delphes_card_CMS'
elif args.delphesEra == 'ATLAS':
    delphesCard = 'delphes_card_ATLAS'
elif args.delphesEra == 'RunIICentral':
    delphesCard = 'delphes_card_CMS_Central'
elif args.delphesEra == 'RunIInoDelphesIso':
    delphesCard = 'delphes_card_CMS_noLepIso'
elif args.delphesEra == 'RunIIPileUp':
    delphesCard = 'delphes_card_CMS_PileUp'
elif args.delphesEra == 'PhaseII':
    delphesCard = 'CMS_PhaseII/CMS_PhaseII_200PU_v03'

# some ad-hoc DELPHES jet selection
if "ATLAS" in args.delphesEra:
    _isGoodRecoJet    = lambda j:isGoodRecoJet(j, minNCharged=0, minNNeutrals=5)
    _isGoodRecoLepton = lambda l:isGoodRecoLepton(l, maxIso = 0.4)
else:
    _isGoodRecoJet    = isGoodRecoJet 
    _isGoodRecoLepton = isGoodRecoLepton

# Check whether we have to do anything
if os.path.exists( output_filename ) and checkRootFile( output_filename, checkForObjects=["Events"]) and args.overwrite =='none' :
    logger.info( "File %s found. Quit.", output_filename )
    sys.exit(0)

# relocate original
sample.copy_files( os.path.join(tmp_output_directory, "input") )

# run Delphes
delphes_file = os.path.join( output_directory, 'delphes', sample.name+'.root' )
if      ( not os.path.exists( delphes_file )) or \
        ( os.path.exists( delphes_file ) and not checkRootFile( delphes_file, checkForObjects=["Delphes"])) or \
        args.overwrite in ['all']:
    logger.debug( "Reproducing delphes file %s", delphes_file)
    delphesProducer = DelphesProducer( card = delphesCard )
    delphesProducer.produce( sample.files, delphes_file)
delphesReader = DelphesReader( Sample.fromFiles( delphes_file, delphes_file, treeName = "Delphes" ) ) # RootTools version
readers.append( delphesReader )

# TreeMaker initialisation
tmp_dir     = ROOT.gDirectory
output_file = ROOT.TFile( output_filename, 'recreate')
output_file.cd()
maker = TreeMaker(
    #sequence  = [ filler ],
    variables = [ (TreeVariable.fromString(x) if type(x)==str else x) for x in variables ],
    treeName = "Events"
    )
tmp_dir.cd()

def filler( event ):

    event.run, event.lumi, event.evt = fwliteReader.evt
    if fwliteReader.position % 100==0: logger.info("At event %i/%i", fwliteReader.position, fwliteReader.nEvents)

    # Weight based 
    if args.addReweights:
        event.nweight = weightInfo.nid
        lhe_weights = fwliteReader.products['lhe'].weights()
        weights      = []
        param_points = []
        for weight in lhe_weights:
            # Store nominal weight (First position!)
            weight_id = weight.id.rstrip('_nlo')
            if weight_id in ['rwgt_1','dummy']: 
                event.rw_nominal = weight.wgt
            #print "Hello weight", weight_id, ( weight_id.lower() in weightInfo_data_lower.keys()) 
            if not weight_id.lower() in weightInfo_data_lower.keys(): 
                continue
            pos = weightInfo_data_lower[weight_id]
            #print "pos", weight.wgt, event.weight_base[pos]
            event.weight_base[pos] = weight.wgt
            weights.append( weight.wgt )
            interpreted_weight = interpret_weight(weight_id.lower()) 
            #for var in weightInfo.variables:
            #    getattr( event, "rw_"+var )[pos] = interpreted_weight[var]
            # weight data for interpolation
            if not hyperPoly.initialized: param_points.append( tuple(interpreted_weight[var.lower()] for var in weightInfo.variables) )

        # get list of values of ref point in specific order
        ref_point_coordinates = [weightInfo.ref_point_coordinates[var] for var in weightInfo.variables]

        # Initialize with Reference Point
        if not hyperPoly.initialized: 
            #print "evt,run,lumi", event.run, event.lumi, event.evt
            #print "ref point", ref_point_coordinates, "param_points", param_points
            #for i_p, p in enumerate(param_points):
                #print "weight", i_p, weights[i_p], " ".join([ "%s=%3.2f"%( weightInfo.variables[i], p[i]) for i in range(len(p)) if p[i]!=0])
            hyperPoly.initialize( param_points, ref_point_coordinates )

        coeff = hyperPoly.get_parametrization( weights )
        event.np = hyperPoly.ndof
        event.chi2_ndof = hyperPoly.chi2_ndof(coeff, weights)
        #logger.debug( "chi2_ndof %f coeff %r", event.chi2_ndof, coeff )
        if event.chi2_ndof>10**-6: logger.warning( "chi2_ndof is large: %f", event.chi2_ndof )
        for n in xrange(hyperPoly.ndof):
            event.p_C[n] = coeff[n]
            if n>0:
                setattr( event, "target_"+target_coeff[n-1], coeff[n]/coeff[0] )

        # convinience coefficient vectors for particlenet training
        for coefficient in args.trainingCoefficients:
            setattr(event, "n"+coefficient, 3)
            getattr(event, coefficient+"_coeff")[0] = event.p_C[0]*10**6
            index_lin  = weightInfo.combinations.index((coefficient,))
            index_quad = weightInfo.combinations.index((coefficient, coefficient))
            getattr(event, coefficient+"_coeff")[1] = event.p_C[index_lin]/event.p_C[0] 
            getattr(event, coefficient+"_coeff")[2] = event.p_C[index_quad]/event.p_C[0] 

    # genJets
    ak8GenJets = fwliteReader.products['ak8GenJets']
    genJets    = filter( genJetId, ak8GenJets )

    # All gen particles
    gp        = fwliteReader.products['gp']
    #gpPacked  = fwliteReader.products['gpPacked']

    # for searching
    search  = GenSearch( gp )

    # generated leptons from SM bosons
    genLeps    = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId()) in [11, 12, 13, 14, 15, 16]  and abs(p.mother(0).pdgId()) in [22, 23, 24, 25], gp)]
    genLeps.sort( key = lambda p: -p[1].pt() )
    genLeps_from_bosons =  [first for first, last in genLeps]
    genLeps_dict = [ {var: getattr(last, var)() for var in lep_varnames} for first, last in genLeps ]
    addIndex( genLeps_dict )
    for i_genLep, (first, last) in enumerate(genLeps):
        mother = first.mother(0) if first.numberOfMothers()>0 else None
        if mother is not None:
            mother_pdgId      = mother.pdgId()
            mother_ascend     = search.ascend(mother)
            grandmother       = mother_ascend.mother(0) if mother_ascend.numberOfMothers()>0 else None
            grandmother_pdgId = grandmother.pdgId() if grandmother is not None else 0
        else:
            mother_pdgId = 0
            grandmother_pdgId = 0
        genLeps_dict[i_genLep]['mother_pdgId']      = mother_pdgId
        genLeps_dict[i_genLep]['grandmother_pdgId'] = grandmother_pdgId
    fill_vector_collection( event, "genLep", lep_all_varnames, genLeps_dict )

    # all gen-tops
    #top_partons = filter( lambda p:abs(p.pdgId())==6 and search.isLast(p), gp)
    #top_partons.sort( key=lambda p:-p.pt() )

    W_partons = filter( lambda p:abs(p.pdgId())==24 and search.isFirst(p), gp)
    #W_partons.sort( key=lambda p:-p.pt() )
    b_partons = filter( lambda p:abs(p.pdgId())==5 and search.isFirst(p) and abs(p.mother(0).pdgId())==6, gp)
    #b_partons.sort( key=lambda p:-p.pt() )

    #for i_W, W in enumerate(W_partons):
    #    print "W",i_W,W_partons[i_W].pdgId(), W_partons[i_W].pt(), W_partons[i_W].eta(),W_partons[i_W].phi()
    #for i_b, b in enumerate(b_partons):
    #    print "b",i_b,b_partons[i_b].pdgId(), b_partons[i_b].pt(), b_partons[i_b].eta(),b_partons[i_b].phi(), "mother",b_partons[i_b].mother(0).pdgId()

    # require at least two W close to the resonance
    if len(W_partons)<2: return

    # sanity
    assert len(W_partons)==len(b_partons)==2 , "Not a ttbar candidate event!"
    assert W_partons[0].pdgId()+W_partons[1].pdgId()==0, "W bosons not OS!" 
    assert b_partons[0].pdgId()+b_partons[1].pdgId()==0, "b quarks not OS!" 

    # two lists of Ws and bs -> resolve sorting ambiguity
    # flip sequence of b in case W[0] and b[0] can not come from the same top/anti-top
    if not ( (W_partons[0].pdgId() == 24 and b_partons[0].pdgId() == 5) or (W_partons[0].pdgId() == -24 and b_partons[0].pdgId() == -5) ):
        b_partons.reverse()
        #print "Reversed" # never happens because of how MG fills LHE

    # make top quark partons from the two lists of W_partons = [W1, W2] and b_partons = [b1, b2] -> t1=W1+b1, t2=W2+b2. The ambiguity is already resolved.
    top_partons = []
    for i_W, W in enumerate(W_partons):
        
        W_p4 = makeP4(W)
        b_p4 = makeP4(b_partons[i_W])
 
        top_partons.append( {'p4': W_p4+b_p4, 'pdgId':6 if W.pdgId()==24 else -6, 
            'W_p4':W_p4, 'b_p4':b_p4, 'W':W, 'b':b_partons[i_W], 
            'W_pt':W_p4.Pt(), 'b_pt':b_p4.Pt(),
        } )

    hadTop_partons = []
    lepTop_partons = []
    for t in top_partons:
        t['pt']   = t['p4'].Pt()
        t['eta']  = t['p4'].Eta()
        t['phi']  = t['p4'].Phi()
        t['mass'] = t['p4'].M()
        W_last = search.descend(t['W'])
        t['isHadronic'] = ( not (abs(W_last.daughter(0).pdgId()) in [11,12,13,14,15,16] ))
        if t['isHadronic']:
            if W_last.daughter(0).pdgId()%2==1: # down-type quarks have odd pdgId 
                t['q1'], t['q2'] = W_last.daughter(0), W_last.daughter(1)
            else:
                t['q1'], t['q2'] = W_last.daughter(1), W_last.daughter(0)
            t['q1_p4'] = makeP4(t['q1'].p4())
            t['q2_p4'] = makeP4(t['q2'].p4())
            hadTop_partons.append(t)
        else:
            if abs(W_last.daughter(0).pdgId()) in [12,14,16]:
                t['nu'], t['lep'] = W_last.daughter(0), W_last.daughter(1) 
            else: 
                t['lep'], t['nu'] = W_last.daughter(0), W_last.daughter(1) 
            lepTop_partons.append(t)

    hadTop_partons.sort( key=lambda t:-t['pt'] )
    lepTop_partons.sort( key=lambda t:-t['pt'] )
    #lepTop_parton = lepTop_partons[0] if len(lepTop_partons)>0 else None
    #lepTop1_parton, lepTop2_parton = lepTop_partons[0], lepTop_partons[1] if len(lepTop_partons)>1 else None, None
    lepTop1_parton, lepTop2_parton = None, None

    if len(lepTop_partons)>1:

	lepTop1_parton = lepTop_partons[0]
	lepTop2_parton = lepTop_partons[1]

        if(lepTop2_parton['pdgId']==6):
	    kt = lepTop2_parton
	    lepTop2_parton = lepTop1_parton
	    lepTop1_parton = kt

    #if lepTop1_parton is not None and lepTop1_parton is not None:
    #	    print lepTop1_parton['pdgId'], lepTop2_parton['pdgId']
    #print "had",len(hadTop_partons), "lep", len(lepTop_partons)

    # make AK8 jets from Delphes EFlow
    '''
    eflow_event = delphesReader.EFlowTrack()+delphesReader.EFlowNeutralHadron()+delphesReader.EFlowPhoton()
    eflow_event.sort( key=lambda p:-p['pt'] )
    eflow_pseudoJets = map( lambda p: fastjet.PseudoJet( p['pt']*cos(p['phi']), p['pt']*sin(p['phi']), p['pt']*sinh(p['eta']), p['pt']*cosh(p['eta']) ), eflow_event ) 
    for i_p, p in enumerate(eflow_pseudoJets):
        p.set_user_index( i_p )
    clustSeq      = fastjet.ClusterSequence( eflow_pseudoJets, ak8 )
    delphesJets   = fastjet.sorted_by_pt(clustSeq.inclusive_jets()) 
    #print delphesJets

    try:
        delphesJet = max( filter( lambda j: abs(j.m() - 172.5) < 50, delphesJets), key = lambda j:j.pt() )
    except ValueError:
        delphesJet = None

    #for i_delphesJet, delphesJet in enumerate(delphesJets):
    #    if delphesJet.pt()>400:
    #        print "AK8", len(delphesJet.constituents()), "mass", delphesJet.m(), delphesJet.pt(), delphesJet.eta(), delphesJet.phi()

    # AK8 Delphes RECO jet candidate
    matched_hadTop_parton = None
    if delphesJet is not None:

        event.delphesJet_mass = delphesJet.m()
        event.delphesJet_pt   = delphesJet.pt()
        event.delphesJet_eta  = delphesJet.eta()
        event.delphesJet_phi  = delphesJet.phi()
        event.delphesJet_nConstituents  = len(delphesJet.constituents())

        delphesJet_constituents = [ eflow_event[p.user_index()] for p in delphesJet.constituents() ]
        delphesJet_constituents.sort(key = lambda p:-p['pt'] )
        delphesJet_constituents = delphesJet_constituents[:nCandMax]

        # softdrop
        eflowCandsVec = ROOT.vector("TLorentzVector")()
        for p in delphesJet_constituents:
            eflowCandsVec.push_back( ROOT.TLorentzVector( p['pt']*cos(p['phi']), p['pt']*sin(p['phi']), p['pt']*sinh(p['eta']), p['pt']*cosh(p['eta'])) )
        eflowSDJets = softDrop.result( eflowCandsVec )
        if eflowSDJets.size()>=1:
            eflowSDJet = eflowSDJets[0]
            event.delphesJet_SDmass = eflowSDJet.m() # softdrop mass

            eflowSDSubJets = eflowSDJet.pieces()
            if len(eflowSDSubJets)>0:
                event.delphesJet_SDsubjet0_eta      = eflowSDSubJets[0].eta()
                event.delphesJet_SDsubjet0_deltaEta = eflowSDSubJets[0].eta() - delphesJet.eta()
                event.delphesJet_SDsubjet0_phi      = eflowSDSubJets[0].phi()
                event.delphesJet_SDsubjet0_deltaPhi = deltaPhi(delphesJet.phi(), eflowSDSubJets[0].phi(), returnAbs=False)
                event.delphesJet_SDsubjet0_deltaR   = sqrt( event.delphesJet_SDsubjet0_deltaEta**2 + event.delphesJet_SDsubjet0_deltaPhi**2 )
                event.delphesJet_SDsubjet0_mass     = eflowSDSubJets[0].m()
            if len(eflowSDSubJets)>1:
                event.delphesJet_SDsubjet1_eta      = eflowSDSubJets[1].eta()
                event.delphesJet_SDsubjet1_deltaEta = eflowSDSubJets[1].eta() - delphesJet.eta()
                event.delphesJet_SDsubjet1_phi      = eflowSDSubJets[1].phi()
                event.delphesJet_SDsubjet1_deltaPhi = deltaPhi(delphesJet.phi(), eflowSDSubJets[1].phi(), returnAbs=False)
                event.delphesJet_SDsubjet1_deltaR   = sqrt( event.delphesJet_SDsubjet1_deltaEta**2 + event.delphesJet_SDsubjet1_deltaPhi**2 )
                event.delphesJet_SDsubjet1_mass     = eflowSDSubJets[1].m()

        ns_tau = nSubjettiness.getTau( 4, eflowCandsVec )
        event.delphesJet_tau1 = ns_tau[0]
        event.delphesJet_tau2 = ns_tau[1]
        event.delphesJet_tau3 = ns_tau[2]
        event.delphesJet_tau4 = ns_tau[3]
        event.delphesJet_tau21 = ns_tau[1]/ns_tau[0] if ns_tau[0]>0 else 0
        event.delphesJet_tau32 = ns_tau[2]/ns_tau[1] if ns_tau[1]>0 else 0

        ecf.setParticles( eflowCandsVec )
        result = ecf.result()
        for i_ecf, (name, _) in enumerate( ecfs ):
             setattr(event, "delphesJet_"+name, result[i_ecf] )

        delphesJet_dict       = {'eta':event.delphesJet_eta, 'phi':event.delphesJet_phi}
        try:
            closest_hadTop_parton = min( hadTop_partons, key=lambda p: deltaR( delphesJet_dict, p) )
            matched_hadTop_parton = closest_hadTop_parton if deltaR( closest_hadTop_parton, delphesJet_dict )<0.6 else None
        except ValueError:
            matched_hadTop_parton, closest_hadTop_parton = None, None
            
        if matched_hadTop_parton is not None:
            event.delphesJet_dR_matched_hadTop_parton = deltaR( delphesJet_dict, matched_hadTop_parton) 
            event.delphesJet_dR_lepTop_parton         = deltaR( delphesJet_dict, lepTop_parton) 

            event.delphesJet_dR_hadTop_q1  = deltaR( {'phi':matched_hadTop_parton['q1_p4'].Phi(),'eta':matched_hadTop_parton['q1_p4'].Eta()}, delphesJet_dict )
            event.delphesJet_dR_hadTop_q2  = deltaR( {'phi':matched_hadTop_parton['q2_p4'].Phi(),'eta':matched_hadTop_parton['q2_p4'].Eta()}, delphesJet_dict )
            event.delphesJet_dR_hadTop_W   = deltaR( {'phi':matched_hadTop_parton['W_p4'] .Phi(),'eta':matched_hadTop_parton['W_p4'].Eta()}, delphesJet_dict )
            event.delphesJet_dR_hadTop_b   = deltaR( {'phi':matched_hadTop_parton['b_p4'] .Phi(),'eta':matched_hadTop_parton['b_p4'].Eta()}, delphesJet_dict )

            event.delphesJet_dR_hadTop_maxq1q2b = max( [ event.delphesJet_dR_hadTop_q1, event.delphesJet_dR_hadTop_q2, event.delphesJet_dR_hadTop_b] )
            #if event.delphesJet_pt>500:
            #    print "Matched!", "mass", matched_hadTop_parton['mass'], event.delphesJet_mass, 'pt/eta/phi', matched_hadTop_parton['pt'],matched_hadTop_parton['eta'],matched_hadTop_parton['phi']

        for p in delphesJet_constituents:
            p['phirel'] = deltaPhi(delphesJet.phi(), p['phi'], returnAbs=False)
            p['etarel'] = p['eta'] - delphesJet.eta() 
            if not p.has_key('charge'):p['charge']=0
        for cat in categories:
            for cand in filter( cat['eflow_func'], delphesJet_constituents):
                cand['type'] = cat['type']
                if cat['name'] in ['ph', 'neh']:
                    cand['charge'] = 0

        fill_vector_collection( event, "eflow", cand_varnames, delphesJet_constituents)
    '''
    '''
    # Top quark parton: The delphesJet-matched hadronically decaying top quark parton and, if the delphesJet has no match, the leading hadronic top quark parton
    if matched_hadTop_parton is not None:
        hadTop_parton = matched_hadTop_parton
    elif len(hadTop_partons)>0:
        hadTop_parton = hadTop_partons[0]
    else:
        hadTop_parton = None

    if hadTop_parton is not None:        

        event.parton_hadTop_pt      = hadTop_parton['pt']
        event.parton_hadTop_eta     = hadTop_parton['eta']
        event.parton_hadTop_phi     = hadTop_parton['phi']
        event.parton_hadTop_mass    = hadTop_parton['mass']
        event.parton_hadTop_pdgId   = hadTop_parton['pdgId']

        event.parton_hadTop_q1_pt   = hadTop_parton['q1'].pt()
        event.parton_hadTop_q1_eta  = hadTop_parton['q1'].eta()
        event.parton_hadTop_q1_phi  = hadTop_parton['q1'].phi()
        event.parton_hadTop_q1_mass = hadTop_parton['q1'].mass()
        event.parton_hadTop_q1_pdgId= hadTop_parton['q1'].pdgId()

        event.parton_hadTop_q2_pt   = hadTop_parton['q2'].pt()
        event.parton_hadTop_q2_eta  = hadTop_parton['q2'].eta()
        event.parton_hadTop_q2_phi  = hadTop_parton['q2'].phi()
        event.parton_hadTop_q2_mass = hadTop_parton['q2'].mass()
        event.parton_hadTop_q2_pdgId= hadTop_parton['q2'].pdgId()

        event.parton_hadTop_b_pt    = hadTop_parton['b'].pt()
        event.parton_hadTop_b_eta   = hadTop_parton['b'].eta()
        event.parton_hadTop_b_phi   = hadTop_parton['b'].phi()
        event.parton_hadTop_b_mass  = hadTop_parton['b'].mass()
        event.parton_hadTop_b_pdgId = hadTop_parton['b'].pdgId()

        event.parton_hadTop_W_pt    = hadTop_parton['W'].pt()
        event.parton_hadTop_W_eta   = hadTop_parton['W'].eta()
        event.parton_hadTop_W_phi   = hadTop_parton['W'].phi()
        event.parton_hadTop_W_mass  = hadTop_parton['W'].mass()
        event.parton_hadTop_W_pdgId = hadTop_parton['W'].pdgId()
        
        # compute theta and phi (Suman approved)
        beam = ROOT.TLorentzVector()
        beam.SetPxPyPzE(0,0,6500,6500)

        boost_t = hadTop_parton['p4'].BoostVector()

        # copy the vectors, originals will still be needed
        W_p4  = copy.deepcopy(hadTop_parton['W_p4'])
        q1_p4 = copy.deepcopy(hadTop_parton['q1_p4'])
        q2_p4 = copy.deepcopy(hadTop_parton['q2_p4'])
        W_p4  .Boost(-boost_t)
        q1_p4 .Boost(-boost_t)
        q2_p4 .Boost(-boost_t)

        n_scatter = ((beam.Vect().Unit()).Cross(W_p4.Vect())).Unit()
        n_decay   = (q1_p4.Vect().Cross(q2_p4.Vect())).Unit()
        sign_flip =  1 if ( ((n_scatter.Cross(n_decay))*(W_p4.Vect())) > 0 ) else -1

        try:
            event.parton_hadTop_decayAngle_phi = sign_flip*acos(n_scatter.Dot(n_decay))
        except ValueError:
            event.parton_hadTop_decayAngle_phi = -100

        boost_W = W_p4.BoostVector()
        q1_p4.Boost(-boost_W)

        try:
            event.parton_hadTop_decayAngle_theta = (W_p4).Angle(q1_p4.Vect())
        except ValueError:
            event.parton_hadTop_decayAngle_theta = -100

        # let's not confuse ourselves later on
        del W_p4, q1_p4, q2_p4
    '''
    # Leptonic top quark partons, for basic ttbar parton-level kinematics 
    if lepTop1_parton is not None:        

        event.parton_lepTop1_pt      = lepTop1_parton['pt']
        event.parton_lepTop1_eta     = lepTop1_parton['eta']
        event.parton_lepTop1_phi     = lepTop1_parton['phi']
        event.parton_lepTop1_mass    = lepTop1_parton['mass']
        event.parton_lepTop1_pdgId   = lepTop1_parton['pdgId']

        event.parton_lepTop1_lep_pt         = lepTop1_parton['lep'].pt()
        event.parton_lepTop1_lep_eta        = lepTop1_parton['lep'].eta()
        event.parton_lepTop1_lep_phi        = lepTop1_parton['lep'].phi()
        event.parton_lepTop1_lep_mass       = lepTop1_parton['lep'].mass()
        event.parton_lepTop1_lep_pdgId      = lepTop1_parton['lep'].pdgId()
        event.parton_lepTop1_nu_pt          = lepTop1_parton['nu'].pt()
        event.parton_lepTop1_nu_eta         = lepTop1_parton['nu'].eta()
        event.parton_lepTop1_nu_phi         = lepTop1_parton['nu'].phi()
        event.parton_lepTop1_nu_mass        = lepTop1_parton['nu'].mass()
        event.parton_lepTop1_nu_pdgId       = lepTop1_parton['nu'].pdgId()

        event.parton_lepTop1_b_pt    = lepTop1_parton['b'].pt()
        event.parton_lepTop1_b_eta   = lepTop1_parton['b'].eta()
        event.parton_lepTop1_b_phi   = lepTop1_parton['b'].phi()
        event.parton_lepTop1_b_mass  = lepTop1_parton['b'].mass()
        event.parton_lepTop1_b_pdgId = lepTop1_parton['b'].pdgId()

        event.parton_lepTop1_W_pt    = lepTop1_parton['W'].pt()
        event.parton_lepTop1_W_eta   = lepTop1_parton['W'].eta()
        event.parton_lepTop1_W_phi   = lepTop1_parton['W'].phi()
        event.parton_lepTop1_W_mass  = lepTop1_parton['W'].mass()
        event.parton_lepTop1_W_pdgId = lepTop1_parton['W'].pdgId()

    if lepTop2_parton is not None:

        event.parton_lepTop2_pt      = lepTop2_parton['pt']
        event.parton_lepTop2_eta     = lepTop2_parton['eta']
        event.parton_lepTop2_phi     = lepTop2_parton['phi']
        event.parton_lepTop2_mass    = lepTop2_parton['mass']
        event.parton_lepTop2_pdgId   = lepTop2_parton['pdgId']

        event.parton_lepTop2_lep_pt         = lepTop2_parton['lep'].pt()
        event.parton_lepTop2_lep_eta        = lepTop2_parton['lep'].eta()
        event.parton_lepTop2_lep_phi        = lepTop2_parton['lep'].phi()
        event.parton_lepTop2_lep_mass       = lepTop2_parton['lep'].mass()
        event.parton_lepTop2_lep_pdgId      = lepTop2_parton['lep'].pdgId()
        event.parton_lepTop2_nu_pt          = lepTop2_parton['nu'].pt()
        event.parton_lepTop2_nu_eta         = lepTop2_parton['nu'].eta()
        event.parton_lepTop2_nu_phi         = lepTop2_parton['nu'].phi()
        event.parton_lepTop2_nu_mass        = lepTop2_parton['nu'].mass()
        event.parton_lepTop2_nu_pdgId       = lepTop2_parton['nu'].pdgId()

        event.parton_lepTop2_b_pt    = lepTop2_parton['b'].pt()
        event.parton_lepTop2_b_eta   = lepTop2_parton['b'].eta()
        event.parton_lepTop2_b_phi   = lepTop2_parton['b'].phi()
        event.parton_lepTop2_b_mass  = lepTop2_parton['b'].mass()
        event.parton_lepTop2_b_pdgId = lepTop2_parton['b'].pdgId()

        event.parton_lepTop2_W_pt    = lepTop2_parton['W'].pt()
        event.parton_lepTop2_W_eta   = lepTop2_parton['W'].eta()
        event.parton_lepTop2_W_phi   = lepTop2_parton['W'].phi()
        event.parton_lepTop2_W_mass  = lepTop2_parton['W'].mass()
        event.parton_lepTop2_W_pdgId = lepTop2_parton['W'].pdgId()



    # full ttbar gen-information available
    #if hadTop_parton and lepTop_parton:
    if lepTop1_parton and lepTop2_parton:
	
	assert lepTop1_parton['pdgId']+lepTop2_parton['pdgId']==0, "Tops not OS!"
	top, antitop = (lepTop1_parton, lepTop2_parton) if lepTop1_parton['pdgId']==6 else (lepTop2_parton, lepTop1_parton)

        top_4v, antitop_4v = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        top_4v.SetPtEtaPhiM(top['pt'],top['eta'],top['phi'],top['mass'])
        antitop_4v.SetPtEtaPhiM(antitop['pt'],antitop['eta'],antitop['phi'],antitop['mass'])

        lep_4v, alep_4v = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        lep_4v.SetPtEtaPhiM((antitop['lep'].p4()).Pt(), (antitop['lep'].p4()).Eta(), (antitop['lep'].p4()).Phi(), (antitop['lep'].p4()).M())
        alep_4v.SetPtEtaPhiM((top['lep'].p4()).Pt(), (top['lep'].p4()).Eta(), (top['lep'].p4()).Phi(), (top['lep'].p4()).M())

        cos_phi_lab           = lep_4v.Vect().Unit().Dot(alep_4v.Vect().Unit())
        abs_delta_phi_ll_lab  = abs( deltaPhi( lep_4v.Phi(), alep_4v.Phi() ) )

        event.parton_cos_phi_lab = cos_phi_lab
        event.parton_abs_delta_phi_ll_lab   = abs_delta_phi_ll_lab

        event.parton_cosThetaPlus_n, event.parton_cosThetaMinus_n, event.parton_cosThetaPlus_r, event.parton_cosThetaMinus_r, event.parton_cosThetaPlus_k, event.parton_cosThetaMinus_k, event.parton_cosThetaPlus_r_star, event.parton_cosThetaMinus_r_star, event.parton_cosThetaPlus_k_star, event.parton_cosThetaMinus_k_star, event.parton_cos_phi = get_spincorrelation_angles(alep_4v,lep_4v,top_4v,antitop_4v)

        # TOP-18-006 table 1 http://cds.cern.ch/record/2649926/files/TOP-18-006-pas.pdf
        
        event.parton_xi_nn = event.parton_cosThetaPlus_n*event.parton_cosThetaMinus_n 
        event.parton_xi_rr = event.parton_cosThetaPlus_r*event.parton_cosThetaMinus_r 
        event.parton_xi_kk = event.parton_cosThetaPlus_k*event.parton_cosThetaMinus_k

        event.parton_xi_nr_plus = event.parton_cosThetaPlus_n*event.parton_cosThetaMinus_r + event.parton_cosThetaPlus_r*event.parton_cosThetaMinus_n
        event.parton_xi_nr_minus= event.parton_cosThetaPlus_n*event.parton_cosThetaMinus_r - event.parton_cosThetaPlus_r*event.parton_cosThetaMinus_n
        event.parton_xi_rk_plus = event.parton_cosThetaPlus_r*event.parton_cosThetaMinus_k + event.parton_cosThetaPlus_k*event.parton_cosThetaMinus_r
        event.parton_xi_rk_minus= event.parton_cosThetaPlus_r*event.parton_cosThetaMinus_k - event.parton_cosThetaPlus_k*event.parton_cosThetaMinus_r
        event.parton_xi_nk_plus = event.parton_cosThetaPlus_n*event.parton_cosThetaMinus_k + event.parton_cosThetaPlus_k*event.parton_cosThetaMinus_n
        event.parton_xi_nk_minus= event.parton_cosThetaPlus_n*event.parton_cosThetaMinus_k - event.parton_cosThetaPlus_k*event.parton_cosThetaMinus_n

    # Delphes AK4 jets 
    allRecoJets = delphesReader.jets()

    ## add btag info
    #for i_btagWP, btagWP in enumerate(btagWPs):
    #    count = 0
    #    for jet in allRecoJets:
    #        btag = ( jet["bTag"] & (2**i_btagWP) > 0 ) # Read b-tag bitmap
    #        jet["bTag_"+btagWP] = btag

    # read jets
    recoJets =  filter( lambda j: _isGoodRecoJet(j), allRecoJets)
    recoJets.sort( key = lambda p:-p['pt'] )
    addIndex( recoJets )

    #print len(allRecoJets), len(recoJets)

    ## upgrade JEC are flavor dependent
    #for jet in allRecoJets:
    #    #btag_ = jet ["bTag_"+default_btagWP]
    #    jet["matchGenBJet"] = min( [999]+[ deltaR( jet, trueBJet ) for trueBJet in trueBjets ] )<0.4

    # make reco b jets
    recoBJets    = filter( lambda j:     j['bTag'] and abs(j['eta'])<2.4 , recoJets )
    recoNonBJets = filter( lambda j:not (j['bTag'] and abs(j['eta'])<2.4), recoJets )

    event.nBTag = len( recoBJets )

    # select AK4 bjets
    recoBj0, recoBj1 = ( recoBJets + recoNonBJets + [None, None] )[:2]

    if recoBj0: fill_vector( event, "recoBj0", recoJet_varnames, recoBj0)
    if recoBj1: fill_vector( event, "recoBj1", recoJet_varnames, recoBj1)

    # read leptons
    allRecoLeps = delphesReader.muons() + delphesReader.electrons()
    allRecoLeps.sort( key = lambda p:-p['pt'] )
    recoLeps =  filter( _isGoodRecoLepton, allRecoLeps )

    #delphesGenLeptons = filter( lambda p: abs(p['pdgId']) in [11,13] and p['status']==1, delphesReader.genParticles() )
    # gen-match leptons with delphes particles
    for recoLep in allRecoLeps:
        #recoLep['genMatched'] = any( deltaR( recoLep, genLep )<0.1 for genLep in delphesGenLeptons )
        recoLep['genMatched'] = any( deltaR( recoLep, genLep )<0.1 for genLep in genLeps_dict )
        #print recoLep, recoLep['genMatched'], recoLep['genMatched2']

        #print recoLep['genMatched'], [deltaR( recoLep, genLep ) for genLep in delphesGenLeptons], recoLep

    # cross-cleaning of reco-objects
    nrecoLeps_uncleaned = len( recoLeps )
    recoLeps = filter( lambda l: (min([999]+[deltaR2(l, j) for j in recoJets if j['pt']>30]) > 0.4**2 ), recoLeps )
    #logger.info( "Before photon cleaning: %i after: %i allRecoLeps: %i, recoLeps %i", nrecoLeps_uncleaned, len(recoLeps), len( allRecoLeps ), len( recoLeps ) )

    # give index to leptons
    addIndex( recoLeps )

    # MET
    recoMet = delphesReader.met()[0]
    event.recoMet_pt  = recoMet['pt']
    event.recoMet_phi = recoMet['phi']

    delphesGenMet = delphesReader.genMet()[0]
    event.delphesGenMet_pt  = delphesGenMet['pt']
    event.delphesGenMet_phi = delphesGenMet['phi']

    # transverse mass for each lepton
    for recoLep in allRecoLeps:
        recoLep['MT'] = sqrt(2*event.recoMet_pt*recoLep['pt']*(1-cos(deltaPhi(recoLep['phi'],event.recoMet_phi))))

    # Store leptons and jets
    fill_vector_collection( event, "recoLep",    recoLep_varnames, recoLeps )
    fill_vector_collection( event, "recoJet",    recoJet_varnames, recoJets )
            
    maker.fill()
    maker.event.init()

counter = 0
for reader in readers:
    reader.start()
maker.start()

while readers[0].run( ):
    for reader in readers[1:]:
        reader.run()

    filler( maker.event )
         
    counter += 1
    if counter == maxEvents:  break

logger.info( "Done with running over %i events.", readers[0].nEvents )

output_file.cd()
maker.tree.Write()
output_file.Close()

logger.info( "Written output file %s", output_filename )

##cleanup delphes file:
if os.path.exists( output_filename ) and args.removeDelphesFiles:
    os.remove( delphes_file )
    logger.info( "Removing Delphes file %s", delphes_file )

if os.path.exists(tmp_output_directory):
    shutil.rmtree(tmp_output_directory)
    logger.info( "Cleaned tmp directory %s", tmp_output_directory )
