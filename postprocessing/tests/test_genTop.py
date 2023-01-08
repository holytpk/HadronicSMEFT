from RootTools.core.standard import *
from Analysis.Tools.WeightInfo import WeightInfo
from Analysis.Tools.HyperPoly  import HyperPoly
import ROOT
from math import cos, cosh, acos

from TMB.Tools.GenSearch              import GenSearch
from TMB.Tools.genObjectSelection     import isGoodGenJet, isGoodGenLepton, isGoodGenPhoton, genJetId
from TMB.Tools.helpers                import deltaPhi, deltaR, deltaR2, cosThetaStar, closestOSDLMassToMZ, checkRootFile

import Analysis.Tools.logger as _logger
import RootTools.core.logger as _logger_rt
logger    = _logger.get_logger(   'INFO', logFile = None)
logger_rt = _logger_rt.get_logger('INFO', logFile = None)

weightInfo = WeightInfo("/users/robert.schoefbeck/tt01j-1l-NPtHad_reweight_card.pkl")
weightInfo.set_order(2)
# get list of values of ref point in correct order
ref_point_coordinates = [weightInfo.ref_point_coordinates[var] for var in weightInfo.variables]

hyperPoly  = HyperPoly( weightInfo.order )

logger.info( "Coefficients: %i (%s), order: %i number of weights: %i", len(weightInfo.variables), ",".join(weightInfo.variables), weightInfo.order,  weightInfo.nid)

max_n = -1

import fastjet
ak8 = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.8, fastjet.E_scheme)
#sd  = fastjet.SoftDrop(,,0.8);
sd = ROOT.SoftDropWrapper(0.0, 0.1, 0.8, 200)
ns = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6 )

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

for _, args in ecfs:
    ecf.addECF( *args )

def make_pseudoJet( obj ):
    return fastjet.PseudoJet( obj.px(), obj.py(), obj.pz(), obj.energy() )

def interpret_weight(weight_id):
    str_s = weight_id.split('_')
    res={}
    for i in range(len(str_s)/2):
        res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
    return res

def varnames( vec_vars ):
    return [v.split('/')[0] for v in vec_vars.split(',')]

top_vars       =  "pt/F,eta/F,phi/F,pdgId/I,mass/F"
top_varnames   =  varnames( top_vars )

jet_read_vars       =  "pt/F,eta/F,phi/F,isMuon/I,isElectron/I,isPhoton/I"
jet_read_varnames   =  varnames( jet_read_vars )

# from here on GEN specific:
#GEN = FWLiteSample.fromFiles("GEN", ["/users/robert.schoefbeck/CMS/CMSSW_10_6_27/src/Samples/cfg/GEN_LO_0j_102X.root"])
GEN = FWLiteSample.fromDirectory("tschRefPointNoWidthRW", "/groups/hephy/cms/robert.schoefbeck/ParticleNet/GEN/t-sch-RefPoint-noWidthRW/")
GEN.reduceFiles(to=1)
GEN.reweight_pkl = "/groups/hephy/cms/robert.schoefbeck/ParticleNet/GEN/t-sch-RefPoint-noWidthRW/t-sch-RefPoint-noWidthRW_reweight_card.pkl"
#GEN = FWLiteSample.fromFiles("GEN", ["root://cms-xrd-global.cern.ch//store/user/schoef/PNet/PNet/220904_093006/0000/GEN_LO_0j_102X_27.root"])

products = {
    'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
    'gen':{'type':'GenEventInfoProduct', 'label':'generator'},
    'gp':{'type':'vector<reco::GenParticle>', 'label':("genParticles")},
    'ak8GenJets':{'type':'vector<reco::GenJet>', 'label':("ak8GenJetsNoNu")},
}

fwliteReader = GEN.fwliteReader( products = products )

fwliteReader.start()
counter=0

p_C_GEN = []
rw_GEN_debug  = [] # not needed to store base point weights

categories = [
    {'name':'e',   'func':lambda p:abs(p.pdgId())==11}, #electrons 
    {'name':'mu',  'func':lambda p:abs(p.pdgId())==13}, #muons 
    {'name':'ph',  'func':lambda p:p.pdgId()==22}, #photons 
    {'name':'chh', 'func':lambda p:abs(p.pdgId())>100 and p.charge()!=0 }, #photons 
    {'name':'neh', 'func':lambda p:abs(p.pdgId())>100 and p.charge()==0 }, #photons 
]

while fwliteReader.run( ):
    lhe_weights = fwliteReader.products['lhe'].weights()
    weights      = []
    param_points = []
    for weight in lhe_weights:
        # Store nominal weight (First position!) 
        if weight.id in ['rwgt_1','dummy']: rw_nominal = weight.wgt
        if not weight.id in weightInfo.id: continue
        pos = weightInfo.data[weight.id]
        weights.append( weight.wgt )
        interpreted_weight = interpret_weight(weight.id) 
        # weight data for interpolation
        if not hyperPoly.initialized:
            param_points.append( tuple(interpreted_weight[var] for var in weightInfo.variables) )
            logger.debug( "Weight %s -> base point %r.", weight.id, param_points[-1] ) 

    rw_GEN_debug.append( weights ) 
    # Initialize with Reference Point

    if not hyperPoly.initialized: 
        #print "ref point", ref_point_coordinates
        #for i_p, p in enumerate(param_points):
        #    print "weight", i_p, weights[i_p]/rw_nominal, p
        hyperPoly.initialize( param_points, ref_point_coordinates )
    coeff = hyperPoly.get_parametrization( weights )
    np = hyperPoly.ndof
    chi2_ndof = hyperPoly.chi2_ndof(coeff, weights)
    logger.debug( "chi2_ndof %f coeff %r", chi2_ndof, coeff )
    if chi2_ndof>10**-6: logger.warning( "chi2_ndof is large: %f", chi2_ndof )

    p_C_GEN.append( [ coeff[n] for n in xrange(hyperPoly.ndof) ] )

    counter+=1
    if counter>=max_n and max_n>0:
        break

    # All gen particles
    gp        = fwliteReader.products['gp']
    #gpPacked  = fwliteReader.products['gpPacked']

    # for searching
    search  = GenSearch( gp )

    # all gen-tops
    gen_tops = filter( lambda p:abs(p.pdgId())==6 and search.isLast(p), gp)

    # jets
    ak8GenJets = fwliteReader.products['ak8GenJets']
    genJets = filter( genJetId, ak8GenJets )

    # hadronic gen_tops
    print 
    W, hadronic_gen_top = None, None
    match_found = False
    for i_gen_top, gen_top in enumerate(gen_tops):
        t_daughters = [ gen_top.daughter(i) for i in range(gen_top.numberOfDaughters())]
        print 'top',i_gen_top,"pdgId", gen_top.pdgId(),"d_pdgid", [d.pdgId() for d in t_daughters]
        # we must always have a b
        b = next( (search.descend(b_) for b_ in t_daughters if abs(b_.pdgId()) == 5), None)
        W = next( (search.descend(W_) for W_ in t_daughters if abs(W_.pdgId()) == 24), None)
        # W resonance
        if W is not None and abs(W.daughter(0).pdgId()) in [1,2,3,4]: 
            quark1 = W.daughter(0)
            quark2 = W.daughter(1)
        # 3-body
        elif len(t_daughters)==3:
            quark1 = next( (q for q in t_daughters if abs(q.pdgId()) in [1,2,3,4]), None)
            quark2 = next( (q for q in t_daughters if abs(q.pdgId()) in [1,2,3,4] and not q.pdgId()==quark1.pdgId()), None)
        # leptonic
        else: 
            quark1, quark2 = None, None

        if quark1 is not None and quark2 is not None and b is not None: 
            hadronic_gen_top = gen_top
            hadronic_gen_top_dict = {var: getattr(hadronic_gen_top, var)() for var in top_varnames}

            # quark2 is the anti-quark
            if quark1.pdgId()<0:
                quark1, quark2 = quark2, quark1

            # 4-vectors
            vec_quark1, vec_quark2, vec_b, vec_W, vec_top = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()

            vec_quark1.SetPtEtaPhiM(quark1.pt(),quark1.eta(),quark1.phi(),quark1.mass())
            vec_quark2.SetPtEtaPhiM(quark2.pt(),quark2.eta(),quark2.phi(),quark2.mass())
            vec_b.SetPtEtaPhiM(b.pt(),b.eta(),b.phi(),b.mass())
            vec_t = vec_quark1 + vec_quark2 + vec_b
            vec_W = vec_quark1 + vec_quark2
            match_found = True 
            break

    if hadronic_gen_top is not None:

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
            phi = sign_flip*acos(n_scatter.Dot(n_decay))
        except ValueError:
            pass
            phi = -100

        boost_W = vec_W.BoostVector()

        vec_quark1.Boost(-boost_W)

        theta = float('nan')
        try:
            theta = (vec_W).Angle(vec_quark1.Vect())
        except ValueError:
            theta = -100

        #print "theta,phi",theta,phi

        # match gen jet
        matched_genJet = next( (j for j in genJets if deltaR( {'eta':j.eta(),'phi':j.phi()}, {'eta':hadronic_gen_top.eta(), 'phi':hadronic_gen_top.phi()}) < 0.6), None)
        if matched_genJet is not None:
            matched_genJet_dict = {var: getattr(matched_genJet, var)() for var in jet_read_varnames}
            #print "Jet", matched_genJet_dict
        else:
            match_found = False 

    if match_found:
        gen_particles = list(matched_genJet.getJetConstituentsQuick())
        for cat in categories:
            cat['cands'] = filter( cat['func'], gen_particles )

        assert sum( map( len, [ cat['cands'] for cat in categories] ) ) == len( gen_particles ), "Missing a gen particle in categorization!!"

        clustSeq      = fastjet.ClusterSequence( map( make_pseudoJet, matched_genJet.getJetConstituentsQuick() ), ak8 )
        sortedJets    = fastjet.sorted_by_pt(clustSeq.inclusive_jets())

        genCandsVec = ROOT.vector("TLorentzVector")()
        #genCandsPseudoVec = ROOT.vector("fastjet::PseudoJet")()
        for p in matched_genJet.getJetConstituentsQuick() :
            genCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
            #genCandsPseudoVec.push_back(fastjet.PseudoJet( p.px(), p.py(), p.pz(), p.energy() ))

        ecf.setParticles( genCandsVec )
        result = ecf.result()
        for i_ecf, (name, _) in enumerate( ecfs ):
            print name, result[i_ecf]
 
        gensdjets = sd.result( genCandsVec )
        if gensdjets.size()>=1:
            sdJet = gensdjets[0]
            sd_mass = sdJet.m() # softdrop mass
        
            print "Softdrop mass", sd_mass

            for i_subjet, subjet in enumerate(sdJet.pieces()):
                print "subjet", i_subjet, subjet.pt(), subjet.eta(), subjet.eta() - matched_genJet.eta(), subjet.phi(), subjet.phi() - matched_genJet.phi()

            maxTau = 5
            ns_tau = ns.getTau( maxTau, genCandsVec )
            for i in range(len(ns_tau)): print 'tau'+str(i), round(ns_tau[i], 3)
        else:
            assert False, ""

        print

#        if len(sortedJets)>=1
#            sd_jet = sd(sortedJets[0]);
#        genjetAK8sdmass[ngenjetAK8] = sd_jet.m();
#        Nsubjettiness nsub1_beta1(1,OnePass_WTA_KT_Axes(), UnnormalizedMeasure(1.));
#        genjetAK8tau1[ngenjetAK8] = nsub1_beta1(sd_jet);
#        Nsubjettiness nsub2_beta1(2,OnePass_WTA_KT_Axes(), UnnormalizedMeasure(1.));
#        genjetAK8tau2[ngenjetAK8] = nsub2_beta1(sd_jet);
#        Nsubjettiness nsub3_beta1(3,OnePass_WTA_KT_Axes(), UnnormalizedMeasure(1.));
#        genjetAK8tau3[ngenjetAK8] = nsub3_beta1(sd_jet);


