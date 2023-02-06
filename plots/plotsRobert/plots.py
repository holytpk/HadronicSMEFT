#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('/tmp/delete.png')

import copy
import operator
import random
from math                           import sqrt, cos, sin, pi, isnan, sinh, cosh, log, copysign

# Analysis
import Analysis.Tools.syncer        as syncer
from   Analysis.Tools.WeightInfo    import WeightInfo
from   Analysis.Tools.helpers       import deltaPhi, deltaR, getObjDict, getCollection

# RootTools
from RootTools.core.standard        import *

# HadronicSMEFT
from HadronicSMEFT.Tools.user                 import plot_directory
from HadronicSMEFT.Tools.delphesCutInterpreter import cutInterpreter
#from HadronicSMEFT.Tools.genObjectSelection   import isGoodGenJet

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='delphes-v3')
argParser.add_argument('--selection',          action='store',      default=None)
argParser.add_argument('--sample',             action='store',      default='TT01j1l_HT800',)
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
argParser.add_argument('--reweight',                                action='store_true',     help='Reweight jet pt?')

args = argParser.parse_args()

# Logger'singlelep-WHJet' if sample.name=='WH' else 'dilep-ZHJet-onZ'
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:         plot_directory += "_small"
if args.reweight: args.plot_directory += "-reweight"

plot_directory = os.path.join(plot_directory, args.plot_directory,  args.sample, args.selection if args.selection is not None else "inc")

# Import samples
import HadronicSMEFT.Samples.genTopJets_v1_pp as samples
    
sample = getattr( samples, args.sample)
 
# WeightInfo
sample.weightInfo = WeightInfo(sample.reweight_pkl)
sample.weightInfo.set_order(2)
sample.read_variables = [VectorTreeVariable.fromString( "p[C/F]", nMax=200 )]

eft_configs = [
    {'color':ROOT.kBlack,       'param':{}, 'tex':"SM"},
    {'color':ROOT.kMagenta-4,   'param':{'ctWRe':-1},  'tex':"Re(c_{tW})=-1",  'binning':[20,0,1.5]},
    {'color':ROOT.kMagenta+2,   'param':{'ctWRe':1},   'tex':"Re(c_{tW})=1",   'binning':[20,0,1.5]},
    {'color':ROOT.kGreen-4,     'param':{'ctGRe':-1},  'tex':"Re(c_{tG})=-1",  'binning':[20,0,1.5]},
    {'color':ROOT.kGreen+2,     'param':{'ctGRe':1},   'tex':"Re(c_{tG})=1",   'binning':[20,0,1.5]},
    {'color':ROOT.kBlue-4,      'param':{'ctGIm':-1},  'tex':"Im(c_{tG})=-1",  'binning':[20,0,1.5]},
    {'color':ROOT.kBlue+2,      'param':{'ctGIm':1},   'tex':"Im(c_{tG})=1",   'binning':[20,0,1.5]},
    ]

for eft in eft_configs:
    eft['func'] = sample.weightInfo.get_weight_func(**eft['param']) 
    eft['name'] = "_".join( ["TT01j"] + ( ["SM"] if len(eft['param'])==0 else [ "_".join([key, str(val)]) for key, val in sorted(eft['param'].iteritems())] ) )

# append fake weight
eft_configs.append( {'color':ROOT.kRed, 'param':{'test':0}, 'tex':'exp((p_{T}-500)*log(3.5))', 'binning':[20,0,1.5], 'name':'exp', 'func':lambda event, sample: event.p_C[0]*3.5**((event.delphesJet_pt-500)/1500.)})
#def make_fake_weight( event, sample ):
#    event.eft_weights.append( 1.5**(event.delphesJet_pt/2000) )
#sequence.append( make_fake_weight )

sequence = []
def make_eft_weights( event, sample):
    #if sample.name!=sample.name:
    #    return
    SM_ref_weight         = eft_configs[0]['func'](event, sample)
    event.eft_weights     = [1] + [eft['func'](event, sample)/SM_ref_weight for eft in eft_configs[1:]]

    #print SM_ref_weight, event.delphesJet_pt,  1.5**((event.delphesJet_pt-500)/1500.)

stack       = Stack( )

sequence.append( make_eft_weights )

eft_weights = [] 

for i_eft, eft in enumerate(eft_configs):
    stack.append( [sample] )
    eft_weights.append( [lambda event, sample, i_eft=i_eft: event.eft_weights[i_eft]] )

lumi  = 1
def weight_getter( branches ):
    getters = [ operator.attrgetter(branch) for branch in branches ]
    def getter( event, sample ):
#        for b, g in zip( branches, getters ):
#            print b, g(event)
#        print
        return reduce( operator.mul , [ g(event) for g in getters ], lumi ) 
    return getter

# Read variables and sequences
jetVars          = ['pt/F', 'eta/F', 'phi/F', 'bTag/F']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','isolationVar/F']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
    "nBTag/I", 
    "recoMet_pt/F", "recoMet_phi/F",
    #"genMet_pt/F", "genMet_phi/F",
    "nrecoJet/I",
    "recoJet[%s]"%(",".join(jetVars)),
    "nrecoLep/I",
    "recoLep[%s]"%(",".join(lepVars)),
    #"lumiweight1fb/F",
    "evt/l", "run/I", "lumi/I",

    "parton_hadTop_pt/F", "parton_hadTop_eta/F", "parton_hadTop_phi/F", "parton_hadTop_mass/F", 
    "parton_hadTop_W_pt/F", "parton_hadTop_W_eta/F", "parton_hadTop_W_phi/F", "parton_hadTop_W_mass/F", 
    "parton_hadTop_q1_pt/F", "parton_hadTop_q1_eta/F", "parton_hadTop_q1_phi/F", "parton_hadTop_q1_mass/F", 
    "parton_hadTop_q2_pt/F", "parton_hadTop_q2_eta/F", "parton_hadTop_q2_phi/F", "parton_hadTop_q2_mass/F", 
    "parton_hadTop_b_pt/F", "parton_hadTop_b_eta/F", "parton_hadTop_b_phi/F", "parton_hadTop_b_mass/F",

    "delphesJet_dR_hadTop_q1/F", " delphesJet_dR_hadTop_q2/F", "delphesJet_dR_hadTop_W/F", "delphesJet_dR_hadTop_b/F", 
    "delphesJet_dR_matched_hadTop_parton/F", "delphesJet_dR_hadTop_maxq1q2b/F", 

    "parton_hadTop_decayAngle_phi/F", "parton_hadTop_decayAngle_phi/F",
    "parton_cosThetaPlus_n/F", "parton_cosThetaMinus_n/F", "parton_cosThetaPlus_r/F", "parton_cosThetaMinus_r/F", "parton_cosThetaPlus_k/F", "parton_cosThetaMinus_k/F", 
    "parton_cosThetaPlus_r_star/F", "parton_cosThetaMinus_r_star/F", "parton_cosThetaPlus_k_star/F", "parton_cosThetaMinus_k_star/F",
    "parton_xi_nn/F", "parton_xi_rr/F", "parton_xi_kk/F", 
    "parton_xi_nr_plus/F", "parton_xi_nr_minus/F", "parton_xi_rk_plus/F", "parton_xi_rk_minus/F", "parton_xi_nk_plus/F", "parton_xi_nk_minus/F", 
    "parton_cos_phi/F", "parton_cos_phi_lab/F", "parton_abs_delta_phi_ll_lab/F",
 
    "delphesJet_pt/F", "delphesJet_eta/F", "delphesJet_phi/F", "delphesJet_mass/F", "delphesJet_nConstituents/I",
    'delphesJet_SDmass/F', 
    'delphesJet_SDsubjet0_eta/F', 'delphesJet_SDsubjet0_deltaEta/F', 'delphesJet_SDsubjet0_phi/F', 'delphesJet_SDsubjet0_deltaPhi/F', 'delphesJet_SDsubjet0_deltaR/F', 'delphesJet_SDsubjet0_mass/F', 
    'delphesJet_SDsubjet1_eta/F', 'delphesJet_SDsubjet1_deltaEta/F', 'delphesJet_SDsubjet1_phi/F', 'delphesJet_SDsubjet1_deltaPhi/F', 'delphesJet_SDsubjet1_deltaR/F', 'delphesJet_SDsubjet1_mass/F', 
    'delphesJet_tau1/F', 'delphesJet_tau2/F', 'delphesJet_tau3/F', 'delphesJet_tau4/F', 'delphesJet_tau21/F', 'delphesJet_tau32/F',
    'delphesJet_ecf1/F', 'delphesJet_ecf2/F', 'delphesJet_ecf3/F', 'delphesJet_ecfC1/F', 'delphesJet_ecfC2/F', 'delphesJet_ecfC3/F', 'delphesJet_ecfD/F', 'delphesJet_ecfDbeta2/F', 'delphesJet_ecfM1/F', 'delphesJet_ecfM2/F', 'delphesJet_ecfM3/F', 'delphesJet_ecfM1beta2/F', 'delphesJet_ecfM2beta2/F', 'delphesJet_ecfM3beta2/F', 'delphesJet_ecfN1/F', 'delphesJet_ecfN2/F', 'delphesJet_ecfN3/F', 'delphesJet_ecfN1beta2/F', 'delphesJet_ecfN2beta2/F', 'delphesJet_ecfN3beta2/F', 'delphesJet_ecfU1/F', 'delphesJet_ecfU2/F', 'delphesJet_ecfU3/F', 'delphesJet_ecfU1beta2/F', 'delphesJet_ecfU2beta2/F', 'delphesJet_ecfU3beta2/F', 
]

preselection = [ 
]

selectionString  = "&&".join( [ c[1] for c in preselection] + ([cutInterpreter.cutString(args.selection)] if args.selection is not None else []))

for sample in stack.samples:
    if selectionString != "":
        sample.addSelectionString( selectionString )
    if args.small:
        #sample.reduceFiles( factor = 30 )
        sample.reduceFiles( to = 15 )

#BITs
#import sys, os, time
#sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
#if signal.name.startswith('WH'):
#    import TMB.BIT.configs.WH_delphes_bkgs as config
#    bits        = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/WH_delphes/v2/")
#    bits_bkgs   = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/WH_delphes_bkgs/first_try/")
#elif signal.name.startswith('ZH'):
#    import TMB.BIT.configs.ZH_delphes_bkgs_comb as config
#    bits        = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes/v2/")
#    #bits_bkgs   = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes_bkgs/first_try/")
#    bits_bkgs   = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes_bkgs_comb/v2/")
#
#bits = [
#    ("BIT_bkgs_cHW",             bits_bkgs[('cHW',)],             ([10,-.2,.8] if args.signal=='WH' else [10, -1,9])), 
#    ("BIT_bkgs_cHW_cHW",         bits_bkgs[('cHW','cHW')],        ([10, 0,10] if args.signal=='WH' else [10,0,4])), 
#    ("BIT_bkgs_cHWtil",          bits_bkgs[('cHWtil',)],          [20,-1,1]), 
#    ("BIT_bkgs_cHWtil_cHWtil",   bits_bkgs[('cHWtil','cHWtil')],  [10,0,4]), 
#    ("BIT_cHWtil_cHWtil",        bits[('cHWtil','cHWtil')],  [10,0,4]), 
#    ("BIT_cHWtil_cHWtil_wide",   bits[('cHWtil','cHWtil')],  [30,0,10]), 
##    ("BIT_bkgs_cHj3",            bits_bkgs[('cHj3',)],            [20,-1,1]), 
##    ("BIT_bkgs_cHj3_cHj3",       bits_bkgs[('cHj3','cHj3')],      [30,-5,5]), 
#]
#sequence.extend( config.sequence )
#
#def bit_predict( event, sample ):
#
#    for var, func in config.mva_variables:
#        setattr( event, var, func(event, sample) )
#    
#    # get model inputs assuming lstm
#    event.features = config.predict_inputs( event, sample )
#    for name, model, _ in bits:
#        #print has_lstm, flat_variables, lstm_jets
#        prediction = model.predict( event.features )
#        setattr( event, name, prediction )
##        if not prediction>-float('inf'):
##            print name, prediction, [[getattr( event, mva_variable) for mva_variable, _ in config.mva_variables]]
##            print "mva_m3", event.mva_m3, "m3", event.m3, "event.nJetGood", event.nJetGood
##            raise RuntimeError("Found NAN prediction?")
#
##    # make optimal discriminator for each cfg
##    for eft_config in eft_configs:
##        param = eft_config['param']
##        if len(param)!=1: continue
##        wc_, val_ = list(param.iteritems())[0]
##        setattr( event, "opt_%s_%f"%(wc_, val_), getattr( event, "BIT_bkgs_%s"%wc_ ) + 0.5*val_*getattr( event, "BIT_bkgs_%s_%s"%(wc_, wc_) )) 
#
#sequence.append( bit_predict )
#
## load keras models
#from keras.models import load_model
#
#if signal.name.startswith('ZH'):
#    keras_models = [
#        ("ZH_TT_WJets", load_model("/groups/hephy/cms/robert.schoefbeck/TMB/models/ZH_TT_WJets/ZH_delphes_bkgs/multiclass_model.h5")),
#    ]
#else:
#    keras_models = [
#        ("WH_TT_WJets", load_model("/groups/hephy/cms/robert.schoefbeck/TMB/models/WH_TT_WJets/WH_delphes_bkgs/multiclass_model.h5")),
#    ]
#
#def keras_predict( event, sample ):
#
#    # get model inputs assuming lstm
#    for name, model in keras_models:
#        prediction = model.predict( event.features.reshape(1,-1) )
#
#        #print prediction
#        for i_val, val in enumerate( prediction[0] ):
#            setattr( event, name+'_'+config.training_samples[i_val].name, val)
#
#sequence.append( keras_predict )

plots        = []

#for model_name, _, binning in bits:
#
#    # 1D discriminator
#    plots.append(Plot(
#        name = model_name,
#        texX = model_name, texY = 'Number of Events / 10 GeV',
#        attribute = lambda event, sample, model_name=model_name: getattr(event, model_name),
#        #binning=Binning.fromThresholds([0, 0.5, 1, 2,3,4,10]),
#        binning   = binning,
#        addOverFlowBin = 'upper',
#    ))
#
#for model_name, model in keras_models:
#    for i_tr_s, tr_s in enumerate( config.training_samples ):
#        disc_name = model_name+'_'+config.training_samples[i_tr_s].name
#        plots.append(Plot(
#            texX = disc_name, texY = 'Number of Events',
#            name = "keras_"+disc_name, 
#            attribute = lambda event, sample, disc_name=disc_name: getattr( event, disc_name ),
#            binning=[50, 0, 1],
#        ))

#for eft_config in eft_configs:
#    param = eft_config['param']
#    if len(param)!=1: continue
#    wc_, val_ = list(param.iteritems())[0]
#    name =  "opt_%s_%f"%(wc_, val_)
#    plots.append(Plot(
#        texX = "q(%s=%3.2f)"%(wc_, val_), texY = 'Number of Events',
#        name =  name, 
#        attribute = lambda event, sample, disc_name=name: getattr( event, disc_name ),
#        binning=eft_config['binning'],
#    ))

##features
#for i_key, (key, _) in enumerate( config.mva_variables ):
#    plots.append(Plot( name = key.replace("mva_", ""),
#      texX = config.plot_options[key]['tex'], texY = 'Number of Events',
#      attribute = lambda event, sample, i_key=i_key: event.features[i_key],
#      binning   =  config.plot_options[key]['binning'],
#    ))

# Use some defaults
Plot.setDefaults(stack = stack, weight = eft_weights)#, addOverFlowBin="upper")

if args.reweight:
 
    logger.info("Compute pt reweighting")
    reweight_pt = Plot( name = "delphesJet_pt_reweight",
      attribute      = lambda event, sample: event.delphesJet_pt,
      texX           = "AK8 jet p_{T}", texY = 'Number of Events',
      binning        =  [30, 500, 2000],
      #addOverFlowBin = "upper",
    )

    plotting.fill([ reweight_pt ], read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

    reweight_histos = reweight_pt.histos_added
    ratio = [ reweight_histos[0][0].Clone() for i in range(len(reweight_histos)) ]
    for i_r, r in enumerate(ratio):
        r.Divide(reweight_histos[i_r][0])

    def make_reweight_histo( i_config ):

        def rw_( event, sample ):
            if i_config==0: return 1
            rw_val = ratio[i_config].GetBinContent( ratio[i_config].FindBin( min([event.delphesJet_pt,1999] )))
            return rw_val

        return rw_

    reweight = [ make_reweight_histo( i_eft ) for i_eft, eft in enumerate(eft_configs) ]

    eft_weights = [] 
    for i_eft, eft in enumerate(eft_configs):
        eft_weights.append( [lambda event, sample, i_eft=i_eft: event.eft_weights[i_eft]*reweight[i_eft](event, sample)] )

    # Reset some defaults
    Plot.setDefaults(stack = stack, weight = eft_weights)#, addOverFlowBin="upper")

plots.append(Plot( name = "parton_hadTop_pt",
  texX = "parton p_{T}(t)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_pt,
  binning   =  [20,0,1000],
))

plots.append(Plot( name = "parton_hadTop_eta",
  texX = "parton #eta(t)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_eta,
  binning   =  [20,-4,4],
))

plots.append(Plot( name = "parton_hadTop_phi",
  texX = "parton #phi(t)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_phi,
  binning   =  [20,-pi,pi],
))

plots.append(Plot( name = "parton_hadTop_mass",
  texX = "parton M(t)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_mass,
  binning   =  [50,150,200],
))

plots.append(Plot( name = "parton_hadTop_W_pt",
  texX = "p_{T}(W)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_W_pt,
  binning   =  [20,0,1000],
))

plots.append(Plot( name = "parton_hadTop_W_eta",
  texX = "#eta(W)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_W_eta,
  binning   =  [20,-4,4],
))

plots.append(Plot( name = "parton_hadTop_W_phi",
  texX = "#phi(W)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_W_phi,
  binning   =  [20,-pi,pi],
))

plots.append(Plot( name = "parton_hadTop_W_mass",
  texX = "M(W)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_W_mass,
  binning   =  [50,50,100],
))

plots.append(Plot( name = "parton_hadTop_q1_pt",
  texX = "p_{T}(q1)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_q1_pt,
  binning   =  [20,0,1000],
))

plots.append(Plot( name = "parton_hadTop_q1_eta",
  texX = "#eta(q1)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_q1_eta,
  binning   =  [20,-4,4],
))

plots.append(Plot( name = "parton_hadTop_q1_phi",
  texX = "#phi(q1)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_q1_phi,
  binning   =  [20,-pi,pi],
))

plots.append(Plot( name = "parton_hadTop_q1_mass",
  texX = "M(q1)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_q1_mass,
  binning   =  [50,150,200],
))

plots.append(Plot( name = "parton_hadTop_q2_pt",
  texX = "p_{T}(q2)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_q2_pt,
  binning   =  [20,0,1000],
))

plots.append(Plot( name = "parton_hadTop_q2_eta",
  texX = "#eta(q2)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_q2_eta,
  binning   =  [20,-4,4],
))

plots.append(Plot( name = "parton_hadTop_q2_phi",
  texX = "#phi(q2)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_q2_phi,
  binning   =  [20,-pi,pi],
))

plots.append(Plot( name = "parton_hadTop_q2_mass",
  texX = "M(q2)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_q2_mass,
  binning   =  [50,150,200],
))

plots.append(Plot( name = "parton_hadTop_b_pt",
  texX = "p_{T}(b)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_b_pt,
  binning   =  [20,0,1000],
))

plots.append(Plot( name = "parton_hadTop_b_eta",
  texX = "#eta(b)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_b_eta,
  binning   =  [20,-4,4],
))

plots.append(Plot( name = "parton_hadTop_b_phi",
  texX = "#phi(b)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_b_phi,
  binning   =  [20,-pi,pi],
))

plots.append(Plot( name = "parton_hadTop_b_mass",
  texX = "M(b)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_b_mass,
  binning   =  [20,0,10],
))

plots.append(Plot( name = "parton_cosThetaPlus_n",
  texX = "cos(#theta^{+}_{n})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaPlus_n,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cosThetaPlus_r",
  texX = "cos(#theta^{+}_{r})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaPlus_r,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cosThetaPlus_k",
  texX = "cos(#theta^{+}_{k})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaPlus_k,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cosThetaMinus_n",
  texX = "cos(#theta^{-}_{n})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaMinus_n,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cosThetaMinus_r",
  texX = "cos(#theta^{-}_{r})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaMinus_r,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cosThetaMinus_k",
  texX = "cos(#theta^{-}_{k})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaMinus_k,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cosThetaMinus_r_star",
  texX = "cos(#theta^{-*}_{r})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaMinus_r_star,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cosThetaPlus_r_star",
  texX = "cos(#theta^{+*}_{r})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaPlus_r_star,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cosThetaMinus_k_star",
  texX = "cos(#theta^{-*}_{k})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaMinus_k_star,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cosThetaPlus_k_star",
  texX = "cos(#theta^{+*}_{k})", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cosThetaPlus_k_star,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_xi_nn",
  texX = "#xi_{nn}", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_xi_nn,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_xi_rr",
  texX = "#xi_{rr}", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_xi_rr,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_xi_kk",
  texX = "#xi_{kk}", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_xi_kk,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_xi_D",
  texX = "D=#sum#xi_{ii}", texY = 'Number of Events',
  attribute = lambda event, sample: (event.parton_xi_nn+event.parton_xi_rr+event.parton_xi_kk),
  binning   =  [24,-1.2,1.2],
))

plots.append(Plot( name = "parton_xi_nr_plus",
  texX = "#xi_{nr}^{+}", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_xi_nr_plus,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_xi_nr_minus",
  texX = "#xi_{nr}^{-}", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_xi_nr_minus,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_xi_rk_plus",
  texX = "#xi_{nr}^{+}", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_xi_rk_plus,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_xi_rk_minus",
  texX = "#xi_{nr}^{-}", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_xi_rk_minus,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_xi_nk_plus",
  texX = "#xi_{nr}^{+}", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_xi_nk_plus,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_xi_nk_minus",
  texX = "#xi_{nr}^{-}", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_xi_nk_minus,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cos_phi",
  texX = "parton cos(#phi(l,l))", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cos_phi,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_cos_phi_lab",
  texX = "parton cos(#phi(l,l)) lab", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_cos_phi_lab,
  binning   =  [20,-1,1],
))

plots.append(Plot( name = "parton_abs_delta_phi_ll_lab",
  texX = "parton |#Delta#phi(l,l)| lab", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_abs_delta_phi_ll_lab,
  binning   =  [20,0,pi],
))

plots.append(Plot( name = "deltaR_jet_q1_wide",
  texX = "#DeltaR(gen-jet, q1)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_q1, 
  binning   =  [60,0,6],
))
plots.append(Plot( name = "deltaR_jet_q2_wide",
  texX = "#DeltaR(gen-jet, q2)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_q2,
  binning   =  [60,0,6],
))

plots.append(Plot( name = "deltaR_jet_W_wide",
  texX = "#DeltaR(gen-jet, W)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_W,
  binning   =  [60,0,6],
))

plots.append(Plot( name = "deltaR_jet_b_wide",
  texX = "#DeltaR(gen-jet, b)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_b,
  binning   =  [60,0,6],
))

plots.append(Plot( name = "deltaR_jet_t_wide",
  texX = "#DeltaR(gen-jet, t)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_matched_hadTop_parton,
  binning   =  [60,0,6],
))

plots.append(Plot( name = "deltaR_jet_maxq1q2b_wide",
  texX = "max #DeltaR(gen-jet, [q1, q2, b])", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_maxq1q2b,
  binning   =  [60,0,6],
))

plots.append(Plot( name = "deltaR_jet_q1",
  texX = "#DeltaR(gen-jet, q1)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_q1, 
  binning   =  [60,0,1.2],
))
plots.append(Plot( name = "deltaR_jet_q2",
  texX = "#DeltaR(gen-jet, q2)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_q2,
  binning   =  [60,0,1.2],
))

plots.append(Plot( name = "deltaR_jet_W",
  texX = "#DeltaR(gen-jet, W)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_W,
  binning   =  [60,0,1.2],
))

plots.append(Plot( name = "deltaR_jet_b",
  texX = "#DeltaR(gen-jet, b)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_b,
  binning   =  [60,0,1.2],
))

plots.append(Plot( name = "deltaR_jet_t",
  texX = "#DeltaR(gen-jet, t)", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_matched_hadTop_parton,
  binning   =  [60,0,1.2],
))

plots.append(Plot( name = "deltaR_jet_maxq1q2b",
  texX = "max #DeltaR(gen-jet, [q1, q2, b])", texY = 'Number of Events',
  attribute = lambda event, sample: event.delphesJet_dR_hadTop_maxq1q2b,
  binning   =  [60,0,1.2],
))

plots.append(Plot( name = "parton_hadTop_decayAngle_phi",
  texX = "#phi (gen decay)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_decayAngle_phi,
  binning   =  [20,-pi,pi],
))

plots.append(Plot( name = "parton_hadTop_decayAngle_theta",
  texX = "#theta (gen decay)", texY = 'Number of Events',
  attribute = lambda event, sample: event.parton_hadTop_decayAngle_phi,
  binning   =  [20,0,pi],
))

plots.append(Plot( name = "nrecoJet",
  texX = "N_{jet}", texY = 'Number of Events',
  attribute = lambda event, sample: event.nrecoJet,
  binning   =  [10,0,10],
))

plots.append(Plot( name = "delphesJet_pt",
  attribute = lambda event, sample: event.delphesJet_pt,
  texX = "AK8 jet p_{T}", texY = 'Number of Events',
  binning   =  [50,500,2000],
))

plots.append(Plot( name = "delphesJet_eta",
  attribute = lambda event, sample: event.delphesJet_eta,
  texX = "AK8 jet #eta", texY = 'Number of Events',
  binning   =  [40,-4,4],
))

plots.append(Plot( name = "delphesJet_phi",
  attribute = lambda event, sample: event.delphesJet_phi,
  texX = "AK8 jet #phi", texY = 'Number of Events',
  binning   =  [40,0,2*pi],
))

plots.append(Plot( name = "delphesJet_mass",
  attribute = lambda event, sample: event.delphesJet_mass,
  texX = "AK8 jet mass", texY = 'Number of Events',
  binning   =  [50,0,300],
))

plots.append(Plot( name = "delphesJet_nConstituents",
  attribute = lambda event, sample: event.delphesJet_nConstituents,
  texX = "AK8 jet nConstituents", texY = 'Number of Events',
  binning   =  [50,0,150],
))

plots.append(Plot( name = "delphesJet_SDmass",
  attribute = lambda event, sample: event.delphesJet_SDmass,
  texX = "AK8 jet M_{SD}", texY = 'Number of Events',
  binning   =  [50,0,300],
))

plots.append(Plot( name = "delphesJet_SDsubjet0_deltaEta",
  attribute = lambda event, sample: event.delphesJet_SDsubjet0_deltaEta,
  texX = "AK8 jet SD #eta(subjet_{0})", texY = 'Number of Events',
  binning   =  [50,-1,1],
))

plots.append(Plot( name = "delphesJet_SDsubjet0_deltaPhi",
  attribute = lambda event, sample: event.delphesJet_SDsubjet0_deltaPhi,
  texX = "AK8 jet SD #phi(subjet_{0})", texY = 'Number of Events',
  binning   =  [50,-1,1],
))

plots.append(Plot( name = "delphesJet_SDsubjet0_deltaR",
  attribute = lambda event, sample: event.delphesJet_SDsubjet0_deltaR,
  texX = "AK8 jet SD #Delta R(subjet_{0})", texY = 'Number of Events',
  binning   =  [50,0,1],
))

plots.append(Plot( name = "delphesJet_SDsubjet0_mass",
  attribute = lambda event, sample: event.delphesJet_SDsubjet0_mass,
  texX = "AK8 jet SD M(subjet_{0})", texY = 'Number of Events',
  binning   =  [50,0,1500],
))

plots.append(Plot( name = "delphesJet_SDsubjet1_deltaEta",
  attribute = lambda event, sample: event.delphesJet_SDsubjet1_deltaEta,
  texX = "AK8 jet SD #eta(subjet_{1})", texY = 'Number of Events',
  binning   =  [50,-1,1],
))

plots.append(Plot( name = "delphesJet_SDsubjet1_deltaPhi",
  attribute = lambda event, sample: event.delphesJet_SDsubjet1_deltaPhi,
  texX = "AK8 jet SD #phi(subjet_{1})", texY = 'Number of Events',
  binning   =  [50,-1,1],
))

plots.append(Plot( name = "delphesJet_SDsubjet1_deltaR",
  attribute = lambda event, sample: event.delphesJet_SDsubjet1_deltaR,
  texX = "AK8 jet SD #Delta R(subjet_{1})", texY = 'Number of Events',
  binning   =  [50,0,1],
))

plots.append(Plot( name = "delphesJet_SDsubjet1_mass",
  attribute = lambda event, sample: event.delphesJet_SDsubjet1_mass,
  texX = "AK8 jet SD M(subjet_{1})", texY = 'Number of Events',
  binning   =  [50,0,1500],
))

plots.append(Plot( name = "delphesJet_tau1",
  attribute = lambda event, sample: event.delphesJet_tau1,
  texX = "AK8 jet #tau_{1}", texY = 'Number of Events',
  binning   =  [50,0,.6],
))

plots.append(Plot( name = "delphesJet_tau2",
  attribute = lambda event, sample: event.delphesJet_tau2,
  texX = "AK8 jet #tau_{2}", texY = 'Number of Events',
  binning   =  [50,0,.5],
))

plots.append(Plot( name = "delphesJet_tau3",
  attribute = lambda event, sample: event.delphesJet_tau3,
  texX = "AK8 jet #tau_{3}", texY = 'Number of Events',
  binning   =  [50,0,.4],
))

plots.append(Plot( name = "delphesJet_tau4",
  attribute = lambda event, sample: event.delphesJet_tau4,
  texX = "AK8 jet #tau_{4}", texY = 'Number of Events',
  binning   =  [50,0,.3],
))

plots.append(Plot( name = "delphesJet_tau21",
  attribute = lambda event, sample: event.delphesJet_tau21,
  texX = "AK8 jet #tau_{21}", texY = 'Number of Events',
  binning   =  [50,0,1],
))

plots.append(Plot( name = "delphesJet_tau32",
  attribute = lambda event, sample: event.delphesJet_tau32,
  texX = "AK8 jet #tau_{32}", texY = 'Number of Events',
  binning   =  [50,0,1],
))

plots.append(Plot( name = "delphesJet_ecf1",
  binning=[50,0,2000], texX="AK8 jet ecf1",
  attribute = lambda event, sample: event.delphesJet_ecf1,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecf2",
  binning=[50,0,400000], texX="AK8 jet ecf2",
  attribute = lambda event, sample: event.delphesJet_ecf2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecf3",
  binning=[50,0,4000000], texX="AK8 jet ecf3",
  attribute = lambda event, sample: event.delphesJet_ecf3,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfC1",
  binning=[50,0,.5], texX="AK8 jet ecfC1",
  attribute = lambda event, sample: event.delphesJet_ecfC1,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfC2",
  binning=[50,0,.5], texX="AK8 jet ecfC2",
  attribute = lambda event, sample: event.delphesJet_ecfC2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfC3",
  binning=[50,0,.5], texX="AK8 jet ecfC3",
  attribute = lambda event, sample: event.delphesJet_ecfC3,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfD",
  binning=[50,0,8], texX="AK8 jet ecfD",
  attribute = lambda event, sample: event.delphesJet_ecfD,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfDbeta2",
  binning=[50,0,20], texX="AK8 jet ecfDbeta2",
  attribute = lambda event, sample: event.delphesJet_ecfDbeta2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfM1",
  binning=[50,0,0.35], texX="AK8 jet ecfM1",
  attribute = lambda event, sample: event.delphesJet_ecfM1,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfM2",
  binning=[50,0,0.2], texX="AK8 jet ecfM2",
  attribute = lambda event, sample: event.delphesJet_ecfM2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfM3",
  binning=[50,0,0.2], texX="AK8 jet ecfM3",
  attribute = lambda event, sample: event.delphesJet_ecfM3,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfM1beta2",
  binning=[50,0,0.35], texX="AK8 jet ecfM1beta2",
  attribute = lambda event, sample: event.delphesJet_ecfM1beta2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfM2beta2",
  binning=[50,0,0.2], texX="AK8 jet ecfM2beta2",
  attribute = lambda event, sample: event.delphesJet_ecfM2beta2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfM3beta2",
  binning=[50,0,0.2], texX="AK8 jet ecfM3beta2",
  attribute = lambda event, sample: event.delphesJet_ecfM3beta2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfN1",
  binning=[50,0,0.5], texX="AK8 jet ecfN1",
  attribute = lambda event, sample: event.delphesJet_ecfN1,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfN2",
  binning=[50,0,0.5], texX="AK8 jet ecfN2",
  attribute = lambda event, sample: event.delphesJet_ecfN2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfN3",
  binning=[50,0,5], texX="AK8 jet ecfN3",
  attribute = lambda event, sample: event.delphesJet_ecfN3,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfN1beta2",
  binning=[50,0,0.5], texX="AK8 jet ecfN1beta2",
  attribute = lambda event, sample: event.delphesJet_ecfN1beta2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfN2beta2",
  binning=[50,0,0.5], texX="AK8 jet ecfN2beta2",
  attribute = lambda event, sample: event.delphesJet_ecfN2beta2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfN3beta2",
  binning=[50,0,5], texX="AK8 jet ecfN3beta2",
  attribute = lambda event, sample: event.delphesJet_ecfN3beta2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfU1",
  binning=[50,0,0.5], texX="AK8 jet ecfU1",
  attribute = lambda event, sample: event.delphesJet_ecfU1,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfU2",
  binning=[50,0,0.04], texX="AK8 jet ecfU2",
  attribute = lambda event, sample: event.delphesJet_ecfU2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfU3",
  binning=[50,0,0.004], texX="AK8 jet ecfU3",
  attribute = lambda event, sample: event.delphesJet_ecfU3,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfU1beta2",
  binning=[50,0,0.5], texX="AK8 jet ecfU1beta2",
  attribute = lambda event, sample: event.delphesJet_ecfU1beta2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfU2beta2",
  binning=[50,0,0.04], texX="AK8 jet ecfU2beta2",
  attribute = lambda event, sample: event.delphesJet_ecfU2beta2,
  texY="Number of Events"))

plots.append(Plot( name = "delphesJet_ecfU3beta2",
  binning=[50,0,0.004], texX="AK8 jet ecfU3beta2",
  attribute = lambda event, sample: event.delphesJet_ecfU3beta2,
  texY="Number of Events"))

plots_reweight = []
if args.reweight:
    Plot.setDefaults(stack = Stack(*(stack[:1])), weight = None)#, addOverFlowBin="upper")
    for i_eft_config, eft_config in enumerate(eft_configs):
        if i_eft_config==0: continue

        name = "_".join(["%s_%3.2f"%(key, value) for key, value in eft_config['param'].iteritems()]).replace('.','p').replace('-','m')
    
        plots_reweight.append(Plot( name = "joint_weight_"+name,
          binning=[52,-.15, 5.05], texX="w(%s)/w(SM)"%name,
          attribute = lambda event, sample, i_eft_config=i_eft_config: operator.itemgetter(i_eft_config)(event.eft_weights),
          texY="Number of Events"))

        plots_reweight.append(Plot( name = "joint_reweight_"+name,
          binning=[52,-.1, 5.1], texX="w(%s)/w(SM)"%name,
          attribute = lambda event, sample, i_eft_config=i_eft_config: operator.itemgetter(i_eft_config)(event.eft_weights)*reweight[i_eft_config](event, sample),
          texY="Number of Events"))

# Text on the plots
def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if hasData else "Delphes Simulation"), 
      #(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

# draw function for plots
def drawPlots(plots, ratio=True, legend_columns=2):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory)
    plot_directory_ = os.path.join(plot_directory_, "log") if log else os.path.join(plot_directory_, "lin")
    for plot in plots:
        if  type(plot)==Plot2D:
            plotting.draw2D( plot,
                       plot_directory = plot_directory_,
                       logX = False, logY = False, logZ = log,
                       drawObjects = drawObjects(),
                       copyIndexPHP = True,
#                       oldColors = True,
                       ) 
        else:
            if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
            subtr = 0 #if args.show_derivatives else len(eft_configs)
            plotting.draw(plot,
              plot_directory = plot_directory_,
              ratio =  {'yRange':(0.8,1.2), 'histos':[(i,0) for i in range(1,len(plot.histos))]} if ratio else None,
              logX = False, logY = log, sorting = False,
              yRange = (.5, "auto") if log else (0, "auto"),
              scaling = {},
              legend =  ( (0.17,0.9-0.05*(sum(map(len, plot.histos))-subtr)/2,1.,0.9), legend_columns),
              drawObjects = drawObjects( ),
              copyIndexPHP = True,
              extensions = ["png"],
            )

plotting.fill(plots+plots_reweight, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

#color EFT
offset = 0 
for plot in plots:
    for i_eft, eft in enumerate(eft_configs):
        try:
            plot.histos[i_eft+offset][0].legendText = eft['tex']
            plot.histos[i_eft+offset][0].style      = styles.lineStyle(eft['color'],width=2)
            plot.histos[i_eft+offset][0].SetName(eft['name'])
        except IndexError:
            pass

drawPlots(plots)

rw_plots = []
for i in range(len(plots_reweight)/2):
    p1, p2 = plots_reweight[2*i:2*i+2]

    p1.histos[0][0].legendText = "nominal  (Mean: %4.3f)"%p1.histos[0][0].GetMean()
    p2.histos[0][0].legendText = "reweight (Mean: %4.3f)"%p2.histos[0][0].GetMean()
    p1.histos[0][0].style = styles.lineStyle( ROOT.kBlack, dashed=False)
    p2.histos[0][0].style = styles.lineStyle( ROOT.kBlack, dashed=True)
    rw_plots.append( 
        Plot.fromHisto( p1.name, p1.histos+p2.histos, texX = p1.texX, texY=p1.texY)
        )
    rw_plots[-1].stack=None

drawPlots(rw_plots, ratio=None, legend_columns=1)

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

syncer.sync()

