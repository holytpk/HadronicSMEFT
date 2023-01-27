''' ObjectSelections
'''

max_lepton_abseta = 2.5
max_jet_abseta = 2.5

def isGoodRecoLepton( l, maxIso = None):
    return l['pt'] > 10 and abs( l['eta'] ) < max_lepton_abseta and abs( int(l['pdgId']) ) in [11,13] and (maxIso is None or l['isolationVar']<maxIso) #eta < 2.5

def isGoodRecoMuon( l, maxIso = None):
    return abs( l['pdgId'] ) == 13 and abs( l['eta'] ) < max_lepton_abseta and l['pt'] > 10 and (maxIso is None or l['isolationVar']<maxIso)#eta < 2.5

def isGoodRecoElectron( l, maxIso = None):
    return abs( l['pdgId'] ) == 11 and abs( l['eta'] ) < max_lepton_abseta and l['pt'] > 10 and (maxIso is None or l['isolationVar']<maxIso)#eta < 2.5

def isGoodRecoJet( j, pt_var = 'pt', max_jet_abseta=max_jet_abseta, minNCharged=1, minNNeutrals=1):
    return  abs( j['eta'] ) < max_jet_abseta and j[pt_var] > 30 and (j['nCharged']>=minNCharged and j['nNeutrals']>=minNNeutrals  or abs(j['eta']) > 2.4 )

def isGoodRecoPhoton( g ):
    return  abs( g['eta'] ) < 2.1 and g['pt'] > 15
