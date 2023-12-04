#!/usr/bin/env python
''' Copy a directory and apply extra selection
'''
#
# Standard imports and batch mode
#
import ROOT
import glob
import os

# RootTools
from RootTools.core.standard             import *

from tttt.Tools.cutInterpreter           import cutInterpreter
from tttt.Tools.user                     import *
# Hello
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
#argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?')
argParser.add_argument('--overwrite',      action='store_true', help='Overwrite?')
argParser.add_argument('--sources',        action='store', nargs = "*", default=['/scratch-cbe/users/robert.schoefbeck/HadronicSMEFT/postprocessed/gen/v6/WZto2L_HT300/', '/scratch-cbe/users/robert.schoefbeck/HadronicSMEFT/postprocessed/gen/v6/DY_HT300/'])
argParser.add_argument('--target',         action='store', default='/scratch-cbe/users/$USER/HadronicSMEFT/postprocessed/gen/v6/WZandDY_selected/')
argParser.add_argument('--selection',      action='store', default='genJet_VV_p>0&&delphesJet_pt>500&&delphesJet_SDmass>70&&delphesJet_SDmass<110&&(dR_genJet_maxq1q2<0.6||(!(dR_genJet_maxq1q2>-1)))')
argParser.add_argument('--cores',          action='store',         type=int, default=-1,                  help="How many jobs to parallelize?" )
args = argParser.parse_args()

# Logger
import tttt.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

sample = Sample.fromDirectory( "sample", directory = args.sources )
target = os.path.expandvars(args.target)

try:
    os.makedirs(target)
except:
    pass

sample.copy_files(target=target, selection=args.selection,  overwrite=args.overwrite)

