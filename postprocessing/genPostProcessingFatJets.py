#!/usr/bin/env python
'''  produce fat jet ntuple 
'''
#
# Standard imports and batch mode
#
import os, sys, imp
import ROOT
ROOT.gROOT.SetBatch(True)
from math                             import sqrt, cos, sin, pi, acos, cosh, sinh

#RootTools
from RootTools.core.standard          import *

#Analysis
from Analysis.Tools.WeightInfo        import WeightInfo
from Analysis.Tools.HyperPoly         import HyperPoly
from Analysis.Tools.GenSearch         import GenSearch
from Analysis.Tools.helpers           import deltaPhi, deltaR, checkRootFile
from Analysis.Tools.DelphesProducer   import DelphesProducer

# HadronicSMEFT
from HadronicSMEFT.Tools.user         import skim_output_directory
from HadronicSMEFT.Tools.genObjectSelection import genJetId

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--overwrite',          action='store',      nargs='?', choices = ['none', 'all', 'target'], default = 'none', help='Overwrite?')#, default = True)
argParser.add_argument('--targetDir',          action='store',      default='v1')
argParser.add_argument('--sample',             action='store',      default='tt1LepHad', help="Name of the sample loaded from fwlite_benchmarks. Only if no inputFiles are specified")
argParser.add_argument('--inputFiles',         action='store',      nargs = '*', default=[])
argParser.add_argument('--delphesEra',         action='store',      default = None, choices = ["RunII", "ATLAS", "RunIICentral", "RunIInoDelphesIso", "RunIIPileUp", "PhaseII"], help="specify delphes era")
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
output_directory = os.path.join(skim_output_directory, 'gen', args.targetDir, sample.name) 

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

# output file & log files
output_filename =  os.path.join(output_directory, sample.name + '.root')
_logger.   add_fileHandler( output_filename.replace('.root', '.log'), args.logLevel )
_logger_rt.add_fileHandler( output_filename.replace('.root', '_rt.log'), args.logLevel )


# small helpers
def varnames( vec_vars ):
    return [v.split('/')[0] for v in vec_vars.split(',')]

def addIndex( collection ):
    for i  in range(len(collection)):
        collection[i]['index'] = i

def fill_vector_collection( event, collection_name, collection_varnames, objects):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects):
        for var in collection_varnames:
            getattr(event, collection_name+"_"+var)[i_obj] = obj[var]

# EDM standard variables
variables  += ["run/I", "lumi/I", "evt/l"]


#if args.addReweights:
#    variables.append('rw_nominal/F')
#    # Lumi weight 1fb / w_0
#    variables.append("ref_lumiweight1fb/F")

variables += ["parton_top_pt/F", "parton_top_eta/F", "parton_top_phi/F", "parton_top_mass/F", "parton_top_pdgId/I", "parton_top_decayAngle_theta/F", "parton_top_decayAngle_phi/F"]
variables += ["parton_q1_pt/F",  "parton_q1_eta/F",  "parton_q1_phi/F",  "parton_q1_mass/F",  "parton_q1_pdgId/I"]
variables += ["parton_q2_pt/F",  "parton_q2_eta/F",  "parton_q2_phi/F",  "parton_q2_mass/F",  "parton_q2_pdgId/I"]
variables += ["parton_b_pt/F",   "parton_b_eta/F",   "parton_b_phi/F",   "parton_b_mass/F",   "parton_b_pdgId/I"]
variables += ["parton_W_pt/F",   "parton_W_eta/F",   "parton_W_phi/F",   "parton_W_mass/F"]

if args.delphesEra is not None:
    variables += ["delphesJet_pt/F", "delphesJet_eta/F", "delphesJet_phi/F", "delphesJet_mass/F", "delphesJet_nConstituents/I"] 
    variables += ['delphesJet_SDmass/F', 
                  'delphesJet_SDsubjet0_eta/F', 'delphesJet_SDsubjet0_deltaEta/F', 'delphesJet_SDsubjet0_phi/F', 'delphesJet_SDsubjet0_deltaPhi/F', 'delphesJet_SDsubjet0_deltaR/F', 'delphesJet_SDsubjet0_mass/F', 
                  'delphesJet_SDsubjet1_eta/F', 'delphesJet_SDsubjet1_deltaEta/F', 'delphesJet_SDsubjet1_phi/F', 'delphesJet_SDsubjet1_deltaPhi/F', 'delphesJet_SDsubjet1_deltaR/F', 'delphesJet_SDsubjet1_mass/F', 
                  'delphesJet_tau1/F', 'delphesJet_tau2/F', 'delphesJet_tau3/F', 'delphesJet_tau4/F', 'delphesJet_tau21/F', 'delphesJet_tau32/F']
    for i_ecf, (name, _) in enumerate( ecfs ):
        variables.append( "delphesJet_%s/F"%name )
    variables += ["dR_delphesJet_q1/F", "dR_delphesJet_q2/F", "dR_delphesJet_W/F", "dR_delphesJet_b/F", "dR_delphesJet_top/F", "dR_delphesJet_maxq1q2b/F"]

variables += ["genJet_pt/F", "genJet_eta/F", "genJet_phi/F", "genJet_mass/F", "genJet_nConstituents/I", "genJet_isMuon/I", "genJet_isElectron/I", "genJet_isPhoton/I"]
variables += ['genJet_SDmass/F', 
              'genJet_SDsubjet0_eta/F', 'genJet_SDsubjet0_deltaEta/F', 'genJet_SDsubjet0_phi/F', 'genJet_SDsubjet0_deltaPhi/F', 'genJet_SDsubjet0_deltaR/F', 'genJet_SDsubjet0_mass/F', 
              'genJet_SDsubjet1_eta/F', 'genJet_SDsubjet1_deltaEta/F', 'genJet_SDsubjet1_phi/F', 'genJet_SDsubjet1_deltaPhi/F', 'genJet_SDsubjet1_deltaR/F', 'genJet_SDsubjet1_mass/F', 
              'genJet_tau1/F', 'genJet_tau2/F', 'genJet_tau3/F', 'genJet_tau4/F', 'genJet_tau21/F', 'genJet_tau32/F']
for i_ecf, (name, _) in enumerate( ecfs ):
    variables.append( "genJet_%s/F"%name )

variables += ["dR_genJet_q1/F", "dR_genJet_q2/F", "dR_genJet_W/F", "dR_genJet_b/F", "dR_genJet_top/F", "dR_genJet_maxq1q2b/F"]

variables += ["gen_cand_sum_pt/F"]

if args.addReweights:
    # for each Wilson coefficient listed in args.trainingCoefficients, store a separate length-3 ntuple of ('w0'*10**6, 'w1', 'w2') to facilitate particle-net training 
    for coefficient in args.trainingCoefficients:    
        variables += [VectorTreeVariable.fromString("%s[coeff/F]"%coefficient, nMax=3 )]

categories = [
    {'name':'e',   'eflow_func':lambda p:abs(p['pdgId'])==11, 'gen_func':lambda p:abs(p.pdgId())==11}, #electrons 
    {'name':'mu',  'eflow_func':lambda p:abs(p['pdgId'])==13, 'gen_func':lambda p:abs(p.pdgId())==13}, #muons 
    {'name':'ph',  'eflow_func':lambda p:abs(p['pdgId'])==22, 'gen_func':lambda p:p.pdgId()==22}, #photons 
    {'name':'chh', 'eflow_func':lambda p:abs(p['pdgId'])>100 and p['charge']!=0, 'gen_func':lambda p:abs(p.pdgId())>100 and p.charge()!=0 }, #charged hadrons 
    {'name':'neh', 'eflow_func':lambda p:abs(p['pdgId'])>100 and p['charge']==0, 'gen_func':lambda p:abs(p.pdgId())>100 and p.charge()==0 }, # neutral hadrons
]

cand_vars           =  "pt/F,etarel/F,phirel/F,eta/F,phi/F,pdgId/I"
cand_varnames       =  varnames( cand_vars )
cand_ch_vars        =  "pt/F,etarel/F,phirel/F,eta/F,phi/F,pdgId/I,charge/I"
cand_ch_varnames    =  varnames( cand_ch_vars )
for cat in categories:
    if cat['name'] in ['ph', 'neh']:
        variables.append( VectorTreeVariable.fromString("gen_%s[%s]"%(cat['name'], cand_vars), nMax=1000 ) )
        variables.append( VectorTreeVariable.fromString("eflow_%s[%s]"%(cat['name'], cand_vars), nMax=100 ) )
    else:
        variables.append( VectorTreeVariable.fromString("gen_%s[%s]"%(cat['name'], cand_ch_vars), nMax=1000 ) )
        variables.append( VectorTreeVariable.fromString("eflow_%s[%s]"%(cat['name'], cand_ch_vars), nMax=100 ) )

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
if args.delphesEra is not None:
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
if os.path.exists( output_filename ) and checkRootFile( output_filename, checkForObjects=["Events"]) and args.overwrite =='none' :
    logger.info( "File %s found. Quit.", output_filename )
    sys.exit(0)
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

    if args.delphesEra is not None:

        eflow_event = delphesReader.EFlowTrack()+delphesReader.EFlowNeutralHadron()+delphesReader.EFlowPhoton()
        eflow_event.sort( key=lambda p:-p['pt'] )
        eflow_pseudoJets = map( lambda p: fastjet.PseudoJet( p['pt']*cos(p['phi']), p['pt']*sin(p['phi']), p['pt']*sinh(p['eta']), p['pt']*cosh(p['eta'])), eflow_event ) 
        for i_p, p in enumerate(eflow_pseudoJets):
            p.set_user_index( i_p )
        clustSeq      = fastjet.ClusterSequence( eflow_pseudoJets, ak8 )
        delphesJets   = fastjet.sorted_by_pt(clustSeq.inclusive_jets()) 

    # genJets
    ak8GenJets = fwliteReader.products['ak8GenJets']
    genJets = filter( genJetId, ak8GenJets )

    # All gen particles
    gp        = fwliteReader.products['gp']
    #gpPacked  = fwliteReader.products['gpPacked']

    # for searching
    search  = GenSearch( gp )

    # all gen-tops
    top_partons = filter( lambda p:abs(p.pdgId())==6 and search.isLast(p), gp)

    # loop over hadronic parton_tops
    for i_top_parton, top_parton in enumerate(top_partons):

        W, is_hadronic_top_parton = None, False

        t_daughters = [ top_parton.daughter(i) for i in range(top_parton.numberOfDaughters())]
        logger.debug( 'top',i_top_parton,"pdgId", top_parton.pdgId(),"d_pdgid", [d.pdgId() for d in t_daughters] )

        # we must always have a b
        b = next( (b_ for b_ in t_daughters if abs(b_.pdgId()) == 5), None)
        W = next( (W_ for W_ in t_daughters if abs(W_.pdgId()) == 24), None)

        # W resonance
        if W is not None: 
            W_last = search.descend(W) 
            if abs(W_last.daughter(0).pdgId()) in [1,2,3,4]:
                quark1 = W_last.daughter(0)
                quark2 = W_last.daughter(1)
                is_hadronic_top_parton = True
                is_hadronic_top_3body  = False
        # 3-body
        elif len(t_daughters)==3:
            quark1 = next( (q for q in t_daughters if abs(q.pdgId()) in [1,2,3,4]), None)
            quark2 = next( (q for q in t_daughters if abs(q.pdgId()) in [1,2,3,4] and not q.pdgId()==quark1.pdgId()), None)
            if quark1 is not None and quark2 is not None:
                is_hadronic_top_parton = True
                is_hadronic_top_3body  = True
        # leptonic
        else:
            quark1, quark2 = None, None

        if not is_hadronic_top_parton or b is None:
            continue

        # quark2 is the anti-quark
        if quark1.pdgId()<0:
            quark1, quark2 = quark2, quark1

        event.parton_top_pt   = top_parton.pt()
        event.parton_top_eta  = top_parton.eta()
        event.parton_top_phi  = top_parton.phi()
        event.parton_top_mass = top_parton.mass()
        event.parton_top_pdgId= top_parton.pdgId()

        event.parton_q1_pt   = quark1.pt()
        event.parton_q1_eta  = quark1.eta()
        event.parton_q1_phi  = quark1.phi()
        event.parton_q1_mass = quark1.mass()
        event.parton_q1_pdgId= quark1.pdgId()

        event.parton_q2_pt   = quark2.pt()
        event.parton_q2_eta  = quark2.eta()
        event.parton_q2_phi  = quark2.phi()
        event.parton_q2_mass = quark2.mass()
        event.parton_q2_pdgId= quark2.pdgId()

        event.parton_b_pt    = b.pt()
        event.parton_b_eta   = b.eta()
        event.parton_b_phi   = b.phi()
        event.parton_b_mass  = b.mass()
        event.parton_b_pdgId = b.pdgId()

        # 4-vectors
        vec_quark1, vec_quark2, vec_b, vec_W, vec_top = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()

        vec_quark1.SetPtEtaPhiM(quark1.pt(),quark1.eta(),quark1.phi(),quark1.mass())
        vec_quark2.SetPtEtaPhiM(quark2.pt(),quark2.eta(),quark2.phi(),quark2.mass())
        vec_b.SetPtEtaPhiM(b.pt(),b.eta(),b.phi(),b.mass())
        vec_t = vec_quark1 + vec_quark2 + vec_b
        vec_W = vec_quark1 + vec_quark2

        event.parton_W_pt   = vec_W.Pt()
        event.parton_W_eta  = vec_W.Eta()
        event.parton_W_phi  = vec_W.Phi()
        event.parton_W_mass = vec_W.M()
        
        #if W is not None: print "W   radiation", W.pt() - search.descend( W ).pt() 
        #current_b = search.ascend( b )
        #
        #while True:
        #    print "b-pt", current_b.pt(), "nD", current_b.numberOfDaughters()
        #    if current_b.numberOfDaughters()>1:
        #        for i_d in range(current_b.numberOfDaughters()):
        #            d = current_b.daughter(i_d)
        #            if abs(d.pdgId())==5:
        #                current_b_ = d
        #                print "Found daughter b", d.pdgId(), d.pt()
        #            else:
        #                print "Radiation",d.pdgId(),d.pt()
        #        current_b = current_b_
        #    else: 
        #        print "Less than 2 daughters"
        #        break
    
        #    
        #print "b   radiation", search.ascend( b ).pt() - b.pt() 
        #print "Q1  radiation", search.ascend( quark1 ).pt() - quark1.pt() 
        #print "Q2  radiation", search.ascend( quark2 ).pt() - quark2.pt() 
        #print ("is_hadronic_top_3body",is_hadronic_top_3body, "partonTop_pt", event.partonTop_pt, "3body: parton_top_pt", event.parton_top_pt)
        #print

        # compute theta and phi
        beam = ROOT.TLorentzVector()
        beam.SetPxPyPzE(0,0,6500,6500)

        boost_t = vec_t.BoostVector()

        vec_W.Boost(-boost_t)
        vec_quark1.Boost(-boost_t)
        vec_quark2.Boost(-boost_t)

        n_scatter = ((beam.Vect().Unit()).Cross(vec_W.Vect())).Unit()
        n_decay   = (vec_quark1.Vect().Cross(vec_quark2.Vect())).Unit()
        sign_flip =  1 if ( ((n_scatter.Cross(n_decay))*(vec_W.Vect())) > 0 ) else -1

        try:
            event.parton_top_decayAngle_phi = sign_flip*acos(n_scatter.Dot(n_decay))
        except ValueError:
            event.parton_top_decayAngle_phi = -100

        boost_W = vec_W.BoostVector()
        vec_quark1.Boost(-boost_W)

        try:
            event.parton_top_decayAngle_theta = (vec_W).Angle(vec_quark1.Vect())
        except ValueError:
            event.parton_top_decayAngle_theta = -100

        # Delphes jets
        if args.delphesEra is not None:
            matched_delphesJet = next( (j for j in delphesJets if deltaR( {'eta':j.eta(),'phi':j.phi()}, {'eta':top_parton.eta(), 'phi':top_parton.phi()}) < 0.6), None)
            if matched_delphesJet:

                event.delphesJet_pt  = matched_delphesJet.pt()
                event.delphesJet_eta = matched_delphesJet.eta()
                event.delphesJet_phi = matched_delphesJet.phi()
                event.delphesJet_mass= matched_delphesJet.m()
                event.delphesJet_nConstituents  = len( matched_delphesJet.constituents() )

                event.dR_delphesJet_q1  = deltaR( {'phi':event.parton_q1_phi,  'eta':event.parton_q1_eta},  {'phi':event.delphesJet_phi, 'eta':event.delphesJet_eta} )
                event.dR_delphesJet_q2  = deltaR( {'phi':event.parton_q2_phi,  'eta':event.parton_q2_eta},  {'phi':event.delphesJet_phi, 'eta':event.delphesJet_eta} )
                event.dR_delphesJet_W   = deltaR( {'phi':event.parton_W_phi,   'eta':event.parton_W_eta},   {'phi':event.delphesJet_phi, 'eta':event.delphesJet_eta} )
                event.dR_delphesJet_b   = deltaR( {'phi':event.parton_b_phi,   'eta':event.parton_b_eta},   {'phi':event.delphesJet_phi, 'eta':event.delphesJet_eta} )
                event.dR_delphesJet_top = deltaR( {'phi':event.parton_top_phi, 'eta':event.parton_top_eta}, {'phi':event.delphesJet_phi, 'eta':event.delphesJet_eta} )

                event.dR_delphesJet_maxq1q2b = max( [ event.dR_delphesJet_q1, event.dR_delphesJet_q2, event.dR_delphesJet_b ] )

                eflow_jet = [ eflow_event[p.user_index()] for p in matched_delphesJet.constituents()[:100] ]
                eflow_jet.sort(key = lambda p:-p['pt'] )
                for p in eflow_jet:
                    p['phirel'] = deltaPhi(matched_delphesJet.phi(), p['phi'], returnAbs=False)
                    p['etarel'] = p['eta'] - matched_delphesJet.eta() 

                for cat in categories:
                    cands = filter( cat['eflow_func'], eflow_jet )
                    if cat['name'] in ['ph', 'neh']:
                        fill_vector_collection( event, "eflow_"+cat['name'], cand_varnames, cands)
                    else:
                        fill_vector_collection( event, "eflow_"+cat['name'], cand_ch_varnames, cands)

                # softdrop
                eflowCandsVec = ROOT.vector("TLorentzVector")()
                for p in eflow_jet:
                    eflowCandsVec.push_back( ROOT.TLorentzVector( p['pt']*cos(p['phi']), p['pt']*sin(p['phi']), p['pt']*sinh(p['eta']), p['pt']*cosh(p['eta'])) )
                eflowSDJets = softDrop.result( eflowCandsVec )
                if eflowSDJets.size()>=1:
                    eflowSDJet = eflowSDJets[0]
                    event.delphesJet_SDmass = eflowSDJet.m() # softdrop mass

                    eflowSDSubJets = eflowSDJet.pieces()
                    if len(eflowSDSubJets)>0:
                        event.delphesJet_SDsubjet0_eta      = eflowSDSubJets[0].eta()
                        event.delphesJet_SDsubjet0_deltaEta = eflowSDSubJets[0].eta() - matched_delphesJet.eta()
                        event.delphesJet_SDsubjet0_phi      = eflowSDSubJets[0].phi()
                        event.delphesJet_SDsubjet0_deltaPhi = deltaPhi(matched_delphesJet.phi(), eflowSDSubJets[0].phi(), returnAbs=False)
                        event.delphesJet_SDsubjet0_deltaR   = sqrt( event.delphesJet_SDsubjet0_deltaPhi**2 + event.delphesJet_SDsubjet0_deltaEta**2 )
                        event.delphesJet_SDsubjet0_mass     = eflowSDSubJets[0].m()
                    if len(eflowSDSubJets)>1:
                        event.delphesJet_SDsubjet1_eta      = eflowSDSubJets[1].eta()
                        event.delphesJet_SDsubjet1_deltaEta = eflowSDSubJets[1].eta() - matched_delphesJet.eta()
                        event.delphesJet_SDsubjet1_phi      = eflowSDSubJets[1].phi()
                        event.delphesJet_SDsubjet1_deltaPhi = deltaPhi(matched_delphesJet.phi(), eflowSDSubJets[1].phi(), returnAbs=False)
                        event.delphesJet_SDsubjet1_deltaR   = sqrt( event.delphesJet_SDsubjet1_deltaPhi**2 + event.delphesJet_SDsubjet1_deltaEta**2 )
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
                    setattr( event, "delphesJet_%s"%name, result[i_ecf] )

        # match gen jet
        matched_genJet = next( (j for j in genJets if deltaR( {'eta':j.eta(),'phi':j.phi()}, {'eta':top_parton.eta(), 'phi':top_parton.phi()}) < 0.6), None)

        if matched_genJet:
            gen_particles = list(matched_genJet.getJetConstituentsQuick())

            event.genJet_pt  = matched_genJet.pt()
            event.genJet_eta = matched_genJet.eta()
            event.genJet_phi = matched_genJet.phi()
            event.genJet_mass= matched_genJet.mass()
            event.genJet_nConstituents  = len( gen_particles )
            event.genJet_isMuon         = matched_genJet.isMuon()
            event.genJet_isElectron     = matched_genJet.isElectron()
            event.genJet_isPhoton       = matched_genJet.isPhoton()

            event.gen_cand_sum_pt = sum([c.pt() for c in gen_particles],0)

            event.dR_genJet_q1  = deltaR( {'phi':event.parton_q1_phi,  'eta':event.parton_q1_eta},  {'phi':event.genJet_phi, 'eta':event.genJet_eta} )
            event.dR_genJet_q2  = deltaR( {'phi':event.parton_q2_phi,  'eta':event.parton_q2_eta},  {'phi':event.genJet_phi, 'eta':event.genJet_eta} )
            event.dR_genJet_W   = deltaR( {'phi':event.parton_W_phi,   'eta':event.parton_W_eta},   {'phi':event.genJet_phi, 'eta':event.genJet_eta} )
            event.dR_genJet_b   = deltaR( {'phi':event.parton_b_phi,   'eta':event.parton_b_eta},   {'phi':event.genJet_phi, 'eta':event.genJet_eta} )
            event.dR_genJet_top = deltaR( {'phi':event.parton_top_phi, 'eta':event.parton_top_eta}, {'phi':event.genJet_phi, 'eta':event.genJet_eta} )

            event.dR_genJet_maxq1q2b = max( [ event.dR_genJet_q1, event.dR_genJet_q2, event.dR_genJet_b ] )

            count = 0 
            for cat in categories:
                cands = filter( cat['gen_func'], gen_particles )
                cands_list = [ {'pt':c.pt(), 'eta':c.eta(), 'phi':c.phi(), 'charge':c.charge(), 'pdgId':c.pdgId()} for c in cands ]
                for p in cands_list:
                    p['phirel'] = deltaPhi(event.genJet_phi, p['phi'], returnAbs=False)
                    p['etarel'] = p['eta'] - event.genJet_eta

                cands_list.sort( key = lambda p:-p['pt'] )
                addIndex( cands_list )
                if cat['name'] in ['ph', 'neh']:
                    fill_vector_collection( event, "gen_"+cat['name'], cand_varnames, cands_list)
                else:
                    fill_vector_collection( event, "gen_"+cat['name'], cand_ch_varnames, cands_list)
                count+=len(cands)
            assert count == len( gen_particles ), "Missing a gen particle in categorization!!"

            # softdrop
            genCandsVec = ROOT.vector("TLorentzVector")()
            for p in matched_genJet.getJetConstituentsQuick() :
                genCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
            genSDJets = softDrop.result( genCandsVec )
            if genSDJets.size()>=1:
                genSDJet = genSDJets[0]
                event.genJet_SDmass = genSDJet.m() # softdrop mass

                genSDSubJets = genSDJet.pieces()
                if len(genSDSubJets)>0:
                    event.genJet_SDsubjet0_eta      = genSDSubJets[0].eta()
                    event.genJet_SDsubjet0_deltaEta = genSDSubJets[0].eta() - matched_genJet.eta()
                    event.genJet_SDsubjet0_phi      = genSDSubJets[0].phi()
                    event.genJet_SDsubjet0_deltaPhi = deltaPhi(matched_genJet.phi(), genSDSubJets[0].phi(), returnAbs=False)
                    event.genJet_SDsubjet0_deltaR   = sqrt( event.genJet_SDsubjet0_deltaPhi**2 + event.genJet_SDsubjet0_deltaEta**2 ) 
                    event.genJet_SDsubjet0_mass     = genSDSubJets[0].m()
                if len(genSDSubJets)>1:
                    event.genJet_SDsubjet1_eta      = genSDSubJets[1].eta()
                    event.genJet_SDsubjet1_deltaEta = genSDSubJets[1].eta() - matched_genJet.eta()
                    event.genJet_SDsubjet1_phi      = genSDSubJets[1].phi()
                    event.genJet_SDsubjet1_deltaPhi = deltaPhi(matched_genJet.phi(), genSDSubJets[1].phi(), returnAbs=False)
                    event.genJet_SDsubjet1_deltaR   = sqrt( event.genJet_SDsubjet1_deltaPhi**2 + event.genJet_SDsubjet1_deltaEta**2 ) 
                    event.genJet_SDsubjet1_mass     = genSDSubJets[1].m()
            
            ns_tau = nSubjettiness.getTau( 4, genCandsVec )
            event.genJet_tau1 = ns_tau[0]
            event.genJet_tau2 = ns_tau[1]
            event.genJet_tau3 = ns_tau[2]
            event.genJet_tau4 = ns_tau[3]
            event.genJet_tau21 = ns_tau[1]/ns_tau[0] if ns_tau[0]>0 else 0
            event.genJet_tau32 = ns_tau[2]/ns_tau[1] if ns_tau[1]>0 else 0

            ecf.setParticles( genCandsVec )
            result = ecf.result()
            for i_ecf, (name, _) in enumerate( ecfs ):
                setattr( event, "genJet_%s"%name, result[i_ecf] )

            # only fill if we have everything. This is a per-top ntuple, not per-events
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
if os.path.exists( output_filename ) and args.delphesEra is not None and args.removeDelphesFiles:
    os.remove( delphes_file )
    logger.info( "Removing Delphes file %s", delphes_file )
