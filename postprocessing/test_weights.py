# Standard imports and batch mode
import ROOT
ROOT.gROOT.SetBatch(True)
from math                             import sqrt, cos, sin, pi, acos, cosh, sinh
import numpy as np
import imp
import os

#RootTools
from RootTools.core.standard          import *

#Analysis
from Analysis.Tools.WeightInfo        import WeightInfo
from Analysis.Tools.HyperPoly2        import HyperPoly

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--samples',            action='store',      nargs='*',  type=str, default=['TT01j2lCARef'], help="List of samples to be post-processed" )
args = argParser.parse_args()

# Logger
import HadronicSMEFT.Tools.logger as _logger
import RootTools.core.logger as _logger_rt
logger    = _logger.get_logger(   "INFO", logFile = None)
logger_rt = _logger_rt.get_logger("INFO", logFile = None)

sample_file = "$CMSSW_BASE/python/HadronicSMEFT/Samples/genCA_v1.py"
all_samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
samples = [ getattr( all_samples, s ) for s in args.samples ]
if len(samples) == 1:
    sample = samples[0] 
else:
    logger.info("Combining %i samples with a total of %i files", len(samples), len(sum( [s.files for s in samples], [])) )
    sample = FWLiteSample.combine( samples[0].name+'_comb', samples )
    sample.reweight_pkl = samples[0].reweight_pkl 
logger.debug( 'Loaded sample %s with %i files.', sample.name, len(sample.files) )

sample.reduceFiles( to=1 )

variables = []

# Load reweight pickle file if supposed to keep weights. 
weightInfo = WeightInfo( sample.reweight_pkl )
weightInfo.set_order(2) 

# coefficients for the weight parametrization
hyperPoly         = HyperPoly(2)

def interpret_weight(weight_id):
    str_s = weight_id.rstrip('_nlo').split('_')
    res={}
    for i in range(len(str_s)/2):
        res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
    return res

# Suddenly only lower case weight.id ... who on earth does such things?
weightInfo_data_lower = {k.lower():val for k, val in weightInfo.data.iteritems()}
weightInfo_data_lower.update(weightInfo.data)

# FWLite reader 
products = {
    'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
    'gen':{'type':'GenEventInfoProduct', 'label':'generator'},
}

# relocate original
#sample.copy_files( os.path.join(tmp_output_directory, "input") )

fwliteReader = sample.fwliteReader( products = products )
fwliteReader.start()

while fwliteReader.run():

    if fwliteReader.position % 100==0: logger.info("At event %i/%i", fwliteReader.position, fwliteReader.nEvents)

    nweight = weightInfo.nid
    lhe_weights = fwliteReader.products['lhe'].weights()
    weights      = []
    param_points = []
    for weight in lhe_weights:
        # Store nominal weight (First position!)
        weight_id = weight.id.rstrip('_nlo')
        if weight_id in ['rwgt_1','dummy']: 
            rw_nominal = weight.wgt
        print "weight", weight_id, ( weight_id.lower() in weightInfo_data_lower.keys()), weight.wgt
        if not weight_id.lower() in weightInfo_data_lower.keys(): 
            continue
        pos = weightInfo_data_lower[weight_id]
        #print "pos", weight.wgt, event.weight_base[pos]
        #weight_base[pos] = weight.wgt
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

    break

C    = hyperPoly.Theta[ [0,1,17] ][:, [0, 1, 17]] 
Cinv = hyperPoly.ThetaInv[ [0,1,17] ][:, [0, 1, 17]]

w   = np.array(weights)[[0,1,17]]
res = np.dot( Cinv, w) 

