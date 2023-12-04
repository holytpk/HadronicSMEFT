#!/usr/bin/env python
import ROOT
ROOT.gROOT.SetBatch(True)
import itertools
import os

from HyperPoly2 import HyperPoly

#FWLite and CMSSW tools
from DataFormats.FWLite import Events, Handle
from PhysicsTools.PythonAnalysis import *

handle = Handle('LHEEventProduct')
events = Events( ["/eos/vbc/experiments/cms/store/user/schoef/TT01j1lCARef/TT01j1lCARef/231130_094505/0000/GEN_LO_0j_102X_186.root"] )

order = 2
hyperPoly = HyperPoly(order)

isFirst = True
for i_event in range( 5 ):

    print "Event", i_event
    events.to(i_event)
    events.getByLabel("externalLHEProducer", handle)
    lhe = handle.product()

    lhe_weights = lhe.weights()
    param_points = []

    weights_at_base_points = []

    for weight in lhe_weights:

        # Store nominal weight (First position!)
        weight_id = weight.id.rstrip('_nlo')
        if weight_id in ['rwgt_1','dummy']: 
            continue 

        weights_at_base_points.append( weight.wgt )

        #interpreted_weight = interpret_weight(weight_id)

        if isFirst:
            param_point = map( lambda s: float(s.replace('m','-').replace('p','.')), weight_id.split('_')[1::2])
            param_points.append( tuple(param_point) )
            param_names = weight_id.split('_')[0::2]

    if isFirst:
        hyperPoly.initialize( param_points )
        isFirst = False

    coefficients = hyperPoly.get_parametrization( weights_at_base_points  )   

    for o in xrange(order+1):
        for i_comb, comb in enumerate( itertools.combinations_with_replacement( range(len(param_names)), o ) ):
            comb_name = ",".join( [ param_names[c] for c in comb] )
            print "%4i %15s   %+6.4f"%(i_comb, comb_name, coefficients[i_comb]/coefficients[0] )
