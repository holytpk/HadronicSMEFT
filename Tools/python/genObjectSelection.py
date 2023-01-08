''' ObjectSelections
'''

max_lepton_abseta = 2.5
max_jet_abseta = 2.5

def genFraction( fwlite_jet, pdgId ):
    return sum( [ fwlite_jet.getGenConstituent(i).pt() for i in range(fwlite_jet.numberOfSourceCandidatePtrs()) if abs(fwlite_jet.getGenConstituent(i).pdgId()) == pdgId ], 0 )/fwlite_jet.pt()

def genJetId( fwlite_jet ):
    return genFraction(fwlite_jet, 13)<0.80 and genFraction(fwlite_jet, 11)<0.80

def isGoodGenJet( j, max_jet_abseta=max_jet_abseta):
    ''' jet object selection
    '''
    return j['pt'] > 30 and abs( j['eta'] ) < max_jet_abseta

def isGoodGenPhoton( j ):
    ''' photon object selection
    '''
    return j['pt'] > 15 and abs( j['eta'] ) < 2.1

def isGoodGenLepton( l ):
    ''' lepton object selection
    '''
    return l['pt'] > 10 and abs( l['eta'] ) < max_lepton_abseta and abs( int(l['pdgId']) ) in [11,13] #eta < 2.5
