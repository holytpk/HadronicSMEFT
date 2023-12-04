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

# define polynomial order
order = 2
hyperPoly = HyperPoly(order)

isFirst = True
for i_event in range( 5 ):

    # read information from the event
    print "Event", i_event
    events.to(i_event)
    events.getByLabel("externalLHEProducer", handle)
    lhe = handle.product()

    lhe_weights = lhe.weights()

    base_points            = []
    weights_at_base_points = []

    # loop over the weights 
    for weight in lhe_weights:

        if weight.id in ['rwgt_1','dummy']: 
            continue 

        weights_at_base_points.append( weight.wgt )

        # at the first event, we compute the coordinates of the base points
        if isFirst:
            base_point = map( lambda s: float(s.replace('m','-').replace('p','.')), weight.id.split('_')[1::2])
            base_points.append( tuple(base_point) )
            param_names = weight.id.split('_')[0::2]

    # at the first event, we "initialize" hyperPly class -> This computes the matrix Theta and its inverse ThetaInv (class members)
    if isFirst:
        hyperPoly.initialize( base_points )
        isFirst = False

    # for EVERY event, we compute the weight-coefficients from the weights at the base points that Madgraph gives us
    coefficients = hyperPoly.get_parametrization( weights_at_base_points  )   

    combinations = []
    counter = 0
    for o in xrange(order+1):
        for _comb in itertools.combinations_with_replacement( range(len(param_names)), o ):

            comb =  [ param_names[c] for c in _comb]
            combinations.append( tuple(comb) )
            comb_name = ",".join( comb  )
            print "%4i %15s   %+6.4f relative: %+6.4f"%(counter, comb_name, coefficients[counter], coefficients[counter]/coefficients[0] )
            counter += 1

    index_const = combinations.index( tuple() )
    index_lin   = combinations.index( ('ctGRe',) )
    index_quad  = combinations.index( ('ctGRe','ctGRe') )

    # Implementation of Eq. 1 for an example parameter point ctGRe=5:
    print "Weight for ctGRe = 5:", coefficients[index_const]+ 5*coefficients[index_lin] + 5**2*coefficients[index_quad]
    
