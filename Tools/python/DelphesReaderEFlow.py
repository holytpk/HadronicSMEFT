''' Class for reading Delphes files.
    Based on Suchitas class, now moved to RootTools reader
'''

# Standard imports
import ROOT

# Delphes Reader RootTools
from RootTools.core.DelphesReaderBase import DelphesReaderBase

class DelphesReader( DelphesReaderBase ): # version RootTools reader

    def EFlowTrack( self ):
        return self.read_collection( 'EFlowTrack', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('PID', 'pdgId'),
                ('Charge', 'charge'),
                ('D0', 'd0'), #('ErrorD0', 'd0Err'),  
                ('DZ', 'dz'), #('ErrorDZ', 'dzErr'),  
            ])

    def EFlowPhoton( self ):
        coll = self.read_collection( 'EFlowPhoton', 
            [   ('ET', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('Ehad', 'Ehad'), ('Eem', 'Eem'),  
            ])
        for c in coll:
            c['charge']=0
            c['pdgId']=22
        return coll

    def EFlowNeutralHadron( self ):
        coll = self.read_collection( 'EFlowNeutralHadron', 
            [   ('ET', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('Ehad', 'Ehad'), ('Eem', 'Eem'),  
            ])
        for c in coll:
            c['charge']=0
            c['pdgId']=130
        return coll

    def muons( self ):
        res = self.read_collection( 'Muon',
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('Charge', 'charge'), ('IsolationVar', 'isolationVar'), ('IsolationVarRhoCorr', 'isolationVarRhoCorr'),
                ('SumPtCharged', 'sumPtCharged'),  ('SumPtNeutral', 'sumPtNeutral'), ('SumPtChargedPU', 'sumPtChargedPU'),  ('SumPt', 'sumPt')
            ])
        for r in res:
            r['pdgId'] = -13*r['charge']
            r['ehadOverEem'] = float('nan')
        return res

    def electrons( self ):
        res = self.read_collection( 'Electron',
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('Charge', 'charge'), ('IsolationVar', 'isolationVar'), ('IsolationVarRhoCorr', 'isolationVarRhoCorr'),
                ('SumPtCharged', 'sumPtCharged'),  ('SumPtNeutral', 'sumPtNeutral'), ('SumPtChargedPU', 'sumPtChargedPU'),  ('SumPt', 'sumPt'),
                ('EhadOverEem','ehadOverEem')
            ])
        for r in res:
            r['pdgId'] = -11*r['charge']
        return res

    def jets( self ):
        return self.read_collection( 'Jet',
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('BTag', 'bTag'), ( 'BTagPhys', 'bTagPhys'), ('Flavor', 'flavor'),
                ('NCharged', 'nCharged'), ('NNeutrals', 'nNeutrals'),
            ])

    def genJets( self ):
        return self.read_collection( 'GenJet',
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
            ])

    #def photons( self ):
    #    return self.read_collection( 'Photon',
    #        [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
    #            ('IsolationVar', 'isolationVar'), ('IsolationVarRhoCorr', 'isolationVarRhoCorr'),
    #            ('SumPtCharged', 'sumPtCharged'),  ('SumPtNeutral', 'sumPtNeutral'), ('SumPtChargedPU', 'sumPtChargedPU'),  ('SumPt', 'sumPt'),
    #            ('EhadOverEem','ehadOverEem')
    #        ])

    def met( self ):
        return self.read_collection( 'MissingET', [('MET', 'pt'), ('Phi', 'phi')] )

    def genMet( self ):
        return self.read_collection( 'GenMissingET', [('MET', 'pt'), ('Phi', 'phi')] )


#OBJ: TObjArray  TObjArray   An array of objects : 0
# OBJ: TLeafElement  Event_  Event_ : 0 at: 0x2e7dc00
# OBJ: TLeafElement  Event.fUniqueID fUniqueID[Event_] : 0 at: 0x2e79350
# OBJ: TLeafElement  Event.fBits fBits[Event_] : 0 at: 0x2e7f1b0
# OBJ: TLeafElement  Event.Number    Number[Event_] : 0 at: 0x2e80760
# OBJ: TLeafElement  Event.ReadTime  ReadTime[Event_] : 0 at: 0x2e81d10
# OBJ: TLeafElement  Event.ProcTime  ProcTime[Event_] : 0 at: 0x2e83320
# OBJ: TLeafElement  Event.ProcessID ProcessID[Event_] : 0 at: 0x2e84930
# OBJ: TLeafElement  Event.MPI   MPI[Event_] : 0 at: 0x2e85f10
# OBJ: TLeafElement  Event.Weight    Weight[Event_] : 0 at: 0x2e874c0
# OBJ: TLeafElement  Event.CrossSection  CrossSection[Event_] : 0 at: 0x2e88ad0
# OBJ: TLeafElement  Event.CrossSectionError CrossSectionError[Event_] : 0 at: 0x2e8a140
# OBJ: TLeafElement  Event.Scale Scale[Event_] : 0 at: 0x2e8b750
# OBJ: TLeafElement  Event.AlphaQED  AlphaQED[Event_] : 0 at: 0x2e8cd30
# OBJ: TLeafElement  Event.AlphaQCD  AlphaQCD[Event_] : 0 at: 0x2e8e340
# OBJ: TLeafElement  Event.ID1   ID1[Event_] : 0 at: 0x2e8f940
# OBJ: TLeafElement  Event.ID2   ID2[Event_] : 0 at: 0x2e90f10
# OBJ: TLeafElement  Event.X1    X1[Event_] : 0 at: 0x2e924e0
# OBJ: TLeafElement  Event.X2    X2[Event_] : 0 at: 0x2e93ab0
# OBJ: TLeafElement  Event.ScalePDF  ScalePDF[Event_] : 0 at: 0x2e950b0
# OBJ: TLeafElement  Event.PDF1  PDF1[Event_] : 0 at: 0x2e966b0
# OBJ: TLeafElement  Event.PDF2  PDF2[Event_] : 0 at: 0x2e97c80
# OBJ: TLeafI    Event_size  Event_size : 0 at: 0x2e9a070
# OBJ: TLeafElement  Weight_ Weight_ : 0 at: 0x2e9dcf0
# OBJ: TLeafElement  Weight.fUniqueID    fUniqueID[Weight_] : 0 at: 0x2e9dc10
# OBJ: TLeafElement  Weight.fBits    fBits[Weight_] : 0 at: 0x2ea0590
# OBJ: TLeafElement  Weight.Weight   Weight[Weight_] : 0 at: 0x2ea2e50
# OBJ: TLeafI    Weight_size Weight_size : 0 at: 0x2ea6530
# OBJ: TLeafElement  Particle_   Particle_ : 0 at: 0x2ea8220
# OBJ: TLeafElement  Particle.fUniqueID  fUniqueID[Particle_] : 0 at: 0x2ea8140
# OBJ: TLeafElement  Particle.fBits  fBits[Particle_] : 0 at: 0x2eaab30
# OBJ: TLeafElement  Particle.PID    PID[Particle_] : 0 at: 0x2eae580
# OBJ: TLeafElement  Particle.Status Status[Particle_] : 0 at: 0x2eb0eb0
# OBJ: TLeafElement  Particle.IsPU   IsPU[Particle_] : 0 at: 0x2eb3800
# OBJ: TLeafElement  Particle.M1 M1[Particle_] : 0 at: 0x2eb6100
# OBJ: TLeafElement  Particle.M2 M2[Particle_] : 0 at: 0x2eb89e0
# OBJ: TLeafElement  Particle.D1 D1[Particle_] : 0 at: 0x2ebb2c0
# OBJ: TLeafElement  Particle.D2 D2[Particle_] : 0 at: 0x2ebdba0
# OBJ: TLeafElement  Particle.Charge Charge[Particle_] : 0 at: 0x2ec04d0
# OBJ: TLeafElement  Particle.Mass   Mass[Particle_] : 0 at: 0x2ec2e20
# OBJ: TLeafElement  Particle.E  E[Particle_] : 0 at: 0x2ec5720
# OBJ: TLeafElement  Particle.Px Px[Particle_] : 0 at: 0x2ec8000
# OBJ: TLeafElement  Particle.Py Py[Particle_] : 0 at: 0x2eca8e0
# OBJ: TLeafElement  Particle.Pz Pz[Particle_] : 0 at: 0x2ecd1c0
# OBJ: TLeafElement  Particle.P  P[Particle_] : 0 at: 0x2ecfaa0
# OBJ: TLeafElement  Particle.PT PT[Particle_] : 0 at: 0x2ed2380
# OBJ: TLeafElement  Particle.Eta    Eta[Particle_] : 0 at: 0x2ed4c60
# OBJ: TLeafElement  Particle.Phi    Phi[Particle_] : 0 at: 0x2ed7540
# OBJ: TLeafElement  Particle.Rapidity   Rapidity[Particle_] : 0 at: 0x2ed9e80
# OBJ: TLeafElement  Particle.T  T[Particle_] : 0 at: 0x2edc7c0
# OBJ: TLeafElement  Particle.X  X[Particle_] : 0 at: 0x2edf0a0
# OBJ: TLeafElement  Particle.Y  Y[Particle_] : 0 at: 0x2ee1980
# OBJ: TLeafElement  Particle.Z  Z[Particle_] : 0 at: 0x2ee4260
# OBJ: TLeafI    Particle_size   Particle_size : 0 at: 0x2ee7960
# OBJ: TLeafElement  Track_  Track_ : 0 at: 0x2ee9710
# OBJ: TLeafElement  Track.fUniqueID fUniqueID[Track_] : 0 at: 0x2ee9640
# OBJ: TLeafElement  Track.fBits fBits[Track_] : 0 at: 0x2eeace0
# OBJ: TLeafElement  Track.PID   PID[Track_] : 0 at: 0x2eec2e0
# OBJ: TLeafElement  Track.Charge    Charge[Track_] : 0 at: 0x2eed8b0
# OBJ: TLeafElement  Track.P P[Track_] : 0 at: 0x2eeee80
# OBJ: TLeafElement  Track.PT    PT[Track_] : 0 at: 0x2ef0450
# OBJ: TLeafElement  Track.Eta   Eta[Track_] : 0 at: 0x2ef1a20
# OBJ: TLeafElement  Track.Phi   Phi[Track_] : 0 at: 0x2ef2ff0
# OBJ: TLeafElement  Track.CtgTheta  CtgTheta[Track_] : 0 at: 0x2ef45f0
# OBJ: TLeafElement  Track.C C[Track_] : 0 at: 0x2ef5bf0
# OBJ: TLeafElement  Track.Mass  Mass[Track_] : 0 at: 0x2ef71c0
# OBJ: TLeafElement  Track.EtaOuter  EtaOuter[Track_] : 0 at: 0x2ef87c0
# OBJ: TLeafElement  Track.PhiOuter  PhiOuter[Track_] : 0 at: 0x2ef9df0
# OBJ: TLeafElement  Track.T T[Track_] : 0 at: 0x2efb3f0
# OBJ: TLeafElement  Track.X X[Track_] : 0 at: 0x2efc9c0
# OBJ: TLeafElement  Track.Y Y[Track_] : 0 at: 0x2efdf90
# OBJ: TLeafElement  Track.Z Z[Track_] : 0 at: 0x2eff560
# OBJ: TLeafElement  Track.TOuter    TOuter[Track_] : 0 at: 0x2f00b30
# OBJ: TLeafElement  Track.XOuter    XOuter[Track_] : 0 at: 0x2f02100
# OBJ: TLeafElement  Track.YOuter    YOuter[Track_] : 0 at: 0x2f036d0
# OBJ: TLeafElement  Track.ZOuter    ZOuter[Track_] : 0 at: 0x2f04ca0
# OBJ: TLeafElement  Track.Xd    Xd[Track_] : 0 at: 0x2f06270
# OBJ: TLeafElement  Track.Yd    Yd[Track_] : 0 at: 0x2f07840
# OBJ: TLeafElement  Track.Zd    Zd[Track_] : 0 at: 0x2f08e10
# OBJ: TLeafElement  Track.L L[Track_] : 0 at: 0x2f0a3e0
# OBJ: TLeafElement  Track.D0    D0[Track_] : 0 at: 0x2f0b9b0
# OBJ: TLeafElement  Track.DZ    DZ[Track_] : 0 at: 0x2f0cf80
# OBJ: TLeafElement  Track.Nclusters Nclusters[Track_] : 0 at: 0x2f0e5a0
# OBJ: TLeafElement  Track.dNdx  dNdx[Track_] : 0 at: 0x2f0fbc0
# OBJ: TLeafElement  Track.ErrorP    ErrorP[Track_] : 0 at: 0x2f11190
# OBJ: TLeafElement  Track.ErrorPT   ErrorPT[Track_] : 0 at: 0x2f12780
# OBJ: TLeafElement  Track.ErrorPhi  ErrorPhi[Track_] : 0 at: 0x2f13da0
# OBJ: TLeafElement  Track.ErrorCtgTheta ErrorCtgTheta[Track_] : 0 at: 0x2f15400
# OBJ: TLeafElement  Track.ErrorT    ErrorT[Track_] : 0 at: 0x2f16a30
# OBJ: TLeafElement  Track.ErrorD0   ErrorD0[Track_] : 0 at: 0x2f18020
# OBJ: TLeafElement  Track.ErrorDZ   ErrorDZ[Track_] : 0 at: 0x2f19630
# OBJ: TLeafElement  Track.ErrorC    ErrorC[Track_] : 0 at: 0x2f1ac20
# OBJ: TLeafElement  Track.ErrorD0Phi    ErrorD0Phi[Track_] : 0 at: 0x2f1c250
# OBJ: TLeafElement  Track.ErrorD0C  ErrorD0C[Track_] : 0 at: 0x2f1d8b0
# OBJ: TLeafElement  Track.ErrorD0DZ ErrorD0DZ[Track_] : 0 at: 0x2f1ef00
# OBJ: TLeafElement  Track.ErrorD0CtgTheta   ErrorD0CtgTheta[Track_] : 0 at: 0x2f20580
# OBJ: TLeafElement  Track.ErrorPhiC ErrorPhiC[Track_] : 0 at: 0x2f21c00
# OBJ: TLeafElement  Track.ErrorPhiDZ    ErrorPhiDZ[Track_] : 0 at: 0x2f23280
# OBJ: TLeafElement  Track.ErrorPhiCtgTheta  ErrorPhiCtgTheta[Track_] : 0 at: 0x2f24910
# OBJ: TLeafElement  Track.ErrorCDZ  ErrorCDZ[Track_] : 0 at: 0x2f25f70
# OBJ: TLeafElement  Track.ErrorCCtgTheta    ErrorCCtgTheta[Track_] : 0 at: 0x2f275d0
# OBJ: TLeafElement  Track.ErrorDZCtgTheta   ErrorDZCtgTheta[Track_] : 0 at: 0x2f28c60
# OBJ: TLeafElement  Track.Particle  Particle[Track_] : 0 at: 0x2f2a2c0
# OBJ: TLeafElement  Track.VertexIndex   VertexIndex[Track_] : 0 at: 0x2f2b990
# OBJ: TLeafI    Track_size  Track_size : 0 at: 0x2f2dde0
# OBJ: TLeafElement  Tower_  Tower_ : 0 at: 0x2f2fa70
# OBJ: TLeafElement  Tower.fUniqueID fUniqueID[Tower_] : 0 at: 0x2f2f9a0
# OBJ: TLeafElement  Tower.fBits fBits[Tower_] : 0 at: 0x2f310b0
# OBJ: TLeafElement  Tower.ET    ET[Tower_] : 0 at: 0x2f32730
# OBJ: TLeafElement  Tower.Eta   Eta[Tower_] : 0 at: 0x2f33d70
# OBJ: TLeafElement  Tower.Phi   Phi[Tower_] : 0 at: 0x2f353b0
# OBJ: TLeafElement  Tower.E E[Tower_] : 0 at: 0x2f369f0
# OBJ: TLeafElement  Tower.T T[Tower_] : 0 at: 0x2f38030
# OBJ: TLeafElement  Tower.NTimeHits NTimeHits[Tower_] : 0 at: 0x2f396c0
# OBJ: TLeafElement  Tower.Eem   Eem[Tower_] : 0 at: 0x2f3ad50
# OBJ: TLeafElement  Tower.Ehad  Ehad[Tower_] : 0 at: 0x2f3c390
# OBJ: TLeafElement  Tower.Etrk  Etrk[Tower_] : 0 at: 0x2f3d9d0
# OBJ: TLeafElement  Tower.Edges Edges[Tower_] : 0 at: 0x2f3f010
# OBJ: TLeafElement  Tower.Particles Particles[Tower_] : 0 at: 0x2f40860
# OBJ: TLeafI    Tower_size  Tower_size : 0 at: 0x2f44090
# OBJ: TLeafElement  EFlowTrack_ EFlowTrack_ : 0 at: 0x2f45e40
# OBJ: TLeafElement  EFlowTrack.fUniqueID    fUniqueID[EFlowTrack_] : 0 at: 0x2f45d60
# OBJ: TLeafElement  EFlowTrack.fBits    fBits[EFlowTrack_] : 0 at: 0x2f47470
# OBJ: TLeafElement  EFlowTrack.PID  PID[EFlowTrack_] : 0 at: 0x2f48b00
# OBJ: TLeafElement  EFlowTrack.Charge   Charge[EFlowTrack_] : 0 at: 0x2f4a160
# OBJ: TLeafElement  EFlowTrack.P    P[EFlowTrack_] : 0 at: 0x2f4b790
# OBJ: TLeafElement  EFlowTrack.PT   PT[EFlowTrack_] : 0 at: 0x2f4cd80
# OBJ: TLeafElement  EFlowTrack.Eta  Eta[EFlowTrack_] : 0 at: 0x2f4e3a0
# OBJ: TLeafElement  EFlowTrack.Phi  Phi[EFlowTrack_] : 0 at: 0x2f4f9d0
# OBJ: TLeafElement  EFlowTrack.CtgTheta CtgTheta[EFlowTrack_] : 0 at: 0x2f51030
# OBJ: TLeafElement  EFlowTrack.C    C[EFlowTrack_] : 0 at: 0x2f52660
# OBJ: TLeafElement  EFlowTrack.Mass Mass[EFlowTrack_] : 0 at: 0x2f53c80
# OBJ: TLeafElement  EFlowTrack.EtaOuter EtaOuter[EFlowTrack_] : 0 at: 0x2f55300
# OBJ: TLeafElement  EFlowTrack.PhiOuter PhiOuter[EFlowTrack_] : 0 at: 0x2f56990
# OBJ: TLeafElement  EFlowTrack.T    T[EFlowTrack_] : 0 at: 0x2f57fc0
# OBJ: TLeafElement  EFlowTrack.X    X[EFlowTrack_] : 0 at: 0x2f59590
# OBJ: TLeafElement  EFlowTrack.Y    Y[EFlowTrack_] : 0 at: 0x2f5ab60
# OBJ: TLeafElement  EFlowTrack.Z    Z[EFlowTrack_] : 0 at: 0x2f5c130
# OBJ: TLeafElement  EFlowTrack.TOuter   TOuter[EFlowTrack_] : 0 at: 0x2f5d760
# OBJ: TLeafElement  EFlowTrack.XOuter   XOuter[EFlowTrack_] : 0 at: 0x2f5edf0
# OBJ: TLeafElement  EFlowTrack.YOuter   YOuter[EFlowTrack_] : 0 at: 0x2f60480
# OBJ: TLeafElement  EFlowTrack.ZOuter   ZOuter[EFlowTrack_] : 0 at: 0x2f61b10
# OBJ: TLeafElement  EFlowTrack.Xd   Xd[EFlowTrack_] : 0 at: 0x2f63160
# OBJ: TLeafElement  EFlowTrack.Yd   Yd[EFlowTrack_] : 0 at: 0x2f64770
# OBJ: TLeafElement  EFlowTrack.Zd   Zd[EFlowTrack_] : 0 at: 0x2f65d80
# OBJ: TLeafElement  EFlowTrack.L    L[EFlowTrack_] : 0 at: 0x2f67370
# OBJ: TLeafElement  EFlowTrack.D0   D0[EFlowTrack_] : 0 at: 0x2f68960
# OBJ: TLeafElement  EFlowTrack.DZ   DZ[EFlowTrack_] : 0 at: 0x2f69f70
# OBJ: TLeafElement  EFlowTrack.Nclusters    Nclusters[EFlowTrack_] : 0 at: 0x2f6b5c0
# OBJ: TLeafElement  EFlowTrack.dNdx dNdx[EFlowTrack_] : 0 at: 0x2f6cc40
# OBJ: TLeafElement  EFlowTrack.ErrorP   ErrorP[EFlowTrack_] : 0 at: 0x2f6e2c0
# OBJ: TLeafElement  EFlowTrack.ErrorPT  ErrorPT[EFlowTrack_] : 0 at: 0x2f6f950
# OBJ: TLeafElement  EFlowTrack.ErrorPhi ErrorPhi[EFlowTrack_] : 0 at: 0x2f70fe0
# OBJ: TLeafElement  EFlowTrack.ErrorCtgTheta    ErrorCtgTheta[EFlowTrack_] : 0 at: 0x2f72670
# OBJ: TLeafElement  EFlowTrack.ErrorT   ErrorT[EFlowTrack_] : 0 at: 0x2f73d00
# OBJ: TLeafElement  EFlowTrack.ErrorD0  ErrorD0[EFlowTrack_] : 0 at: 0x2f75390
# OBJ: TLeafElement  EFlowTrack.ErrorDZ  ErrorDZ[EFlowTrack_] : 0 at: 0x2f76a20
# OBJ: TLeafElement  EFlowTrack.ErrorC   ErrorC[EFlowTrack_] : 0 at: 0x2f780b0
# OBJ: TLeafElement  EFlowTrack.ErrorD0Phi   ErrorD0Phi[EFlowTrack_] : 0 at: 0x2f79740
# OBJ: TLeafElement  EFlowTrack.ErrorD0C ErrorD0C[EFlowTrack_] : 0 at: 0x2f7add0
# OBJ: TLeafElement  EFlowTrack.ErrorD0DZ    ErrorD0DZ[EFlowTrack_] : 0 at: 0x2f7c460
# OBJ: TLeafElement  EFlowTrack.ErrorD0CtgTheta  ErrorD0CtgTheta[EFlowTrack_] : 0 at: 0x2f7daf0
# OBJ: TLeafElement  EFlowTrack.ErrorPhiC    ErrorPhiC[EFlowTrack_] : 0 at: 0x2f7f180
# OBJ: TLeafElement  EFlowTrack.ErrorPhiDZ   ErrorPhiDZ[EFlowTrack_] : 0 at: 0x2f80810
# OBJ: TLeafElement  EFlowTrack.ErrorPhiCtgTheta ErrorPhiCtgTheta[EFlowTrack_] : 0 at: 0x2f81ea0
# OBJ: TLeafElement  EFlowTrack.ErrorCDZ ErrorCDZ[EFlowTrack_] : 0 at: 0x2f83530
# OBJ: TLeafElement  EFlowTrack.ErrorCCtgTheta   ErrorCCtgTheta[EFlowTrack_] : 0 at: 0x2f84bc0
# OBJ: TLeafElement  EFlowTrack.ErrorDZCtgTheta  ErrorDZCtgTheta[EFlowTrack_] : 0 at: 0x2f86250
# OBJ: TLeafElement  EFlowTrack.Particle Particle[EFlowTrack_] : 0 at: 0x2f878e0
# OBJ: TLeafElement  EFlowTrack.VertexIndex  VertexIndex[EFlowTrack_] : 0 at: 0x2f88fe0
# OBJ: TLeafI    EFlowTrack_size EFlowTrack_size : 0 at: 0x2f8b460
# OBJ: TLeafElement  EFlowPhoton_    EFlowPhoton_ : 0 at: 0x2f8d170
# OBJ: TLeafElement  EFlowPhoton.fUniqueID   fUniqueID[EFlowPhoton_] : 0 at: 0x2f8d090
# OBJ: TLeafElement  EFlowPhoton.fBits   fBits[EFlowPhoton_] : 0 at: 0x2f8e7e0
# OBJ: TLeafElement  EFlowPhoton.ET  ET[EFlowPhoton_] : 0 at: 0x2f8feb0
# OBJ: TLeafElement  EFlowPhoton.Eta Eta[EFlowPhoton_] : 0 at: 0x2f91540
# OBJ: TLeafElement  EFlowPhoton.Phi Phi[EFlowPhoton_] : 0 at: 0x2f92bf0
# OBJ: TLeafElement  EFlowPhoton.E   E[EFlowPhoton_] : 0 at: 0x2f94270
# OBJ: TLeafElement  EFlowPhoton.T   T[EFlowPhoton_] : 0 at: 0x2f958c0
# OBJ: TLeafElement  EFlowPhoton.NTimeHits   NTimeHits[EFlowPhoton_] : 0 at: 0x2f96f50
# OBJ: TLeafElement  EFlowPhoton.Eem Eem[EFlowPhoton_] : 0 at: 0x2f98610
# OBJ: TLeafElement  EFlowPhoton.Ehad    Ehad[EFlowPhoton_] : 0 at: 0x2f99cd0
# OBJ: TLeafElement  EFlowPhoton.Etrk    Etrk[EFlowPhoton_] : 0 at: 0x2f9b3a0
# OBJ: TLeafElement  EFlowPhoton.Edges   Edges[EFlowPhoton_] : 0 at: 0x2f9ca70
# OBJ: TLeafElement  EFlowPhoton.Particles   Particles[EFlowPhoton_] : 0 at: 0x2f9e250
# OBJ: TLeafI    EFlowPhoton_size    EFlowPhoton_size : 0 at: 0x2fa1690
# OBJ: TLeafElement  EFlowNeutralHadron_ EFlowNeutralHadron_ : 0 at: 0x2519af0
# OBJ: TLeafElement  EFlowNeutralHadron.fUniqueID    fUniqueID[EFlowNeutralHadron_] : 0 at: 0x2519a10
# OBJ: TLeafElement  EFlowNeutralHadron.fBits    fBits[EFlowNeutralHadron_] : 0 at: 0x251b180
# OBJ: TLeafElement  EFlowNeutralHadron.ET   ET[EFlowNeutralHadron_] : 0 at: 0x1d26b00
# OBJ: TLeafElement  EFlowNeutralHadron.Eta  Eta[EFlowNeutralHadron_] : 0 at: 0x1d28190
# OBJ: TLeafElement  EFlowNeutralHadron.Phi  Phi[EFlowNeutralHadron_] : 0 at: 0x2faee40
# OBJ: TLeafElement  EFlowNeutralHadron.E    E[EFlowNeutralHadron_] : 0 at: 0x2fb0450
# OBJ: TLeafElement  EFlowNeutralHadron.T    T[EFlowNeutralHadron_] : 0 at: 0x2fb1ac0
# OBJ: TLeafElement  EFlowNeutralHadron.NTimeHits    NTimeHits[EFlowNeutralHadron_] : 0 at: 0x2fb3150
# OBJ: TLeafElement  EFlowNeutralHadron.Eem  Eem[EFlowNeutralHadron_] : 0 at: 0x2fb47e0
# OBJ: TLeafElement  EFlowNeutralHadron.Ehad Ehad[EFlowNeutralHadron_] : 0 at: 0x2fb5e70
# OBJ: TLeafElement  EFlowNeutralHadron.Etrk Etrk[EFlowNeutralHadron_] : 0 at: 0x2fb7500
# OBJ: TLeafElement  EFlowNeutralHadron.Edges    Edges[EFlowNeutralHadron_] : 0 at: 0x2fb8b90
# OBJ: TLeafElement  EFlowNeutralHadron.Particles    Particles[EFlowNeutralHadron_] : 0 at: 0x2fba2c0
# OBJ: TLeafI    EFlowNeutralHadron_size EFlowNeutralHadron_size : 0 at: 0x2fbc890
# OBJ: TLeafElement  GenJet_ GenJet_ : 0 at: 0x2fbe6c0
# OBJ: TLeafElement  GenJet.fUniqueID    fUniqueID[GenJet_] : 0 at: 0x2fbe5e0
# OBJ: TLeafElement  GenJet.fBits    fBits[GenJet_] : 0 at: 0x2fbfc90
# OBJ: TLeafElement  GenJet.PT   PT[GenJet_] : 0 at: 0x2fc1260
# OBJ: TLeafElement  GenJet.Eta  Eta[GenJet_] : 0 at: 0x2fc2830
# OBJ: TLeafElement  GenJet.Phi  Phi[GenJet_] : 0 at: 0x2fc3e00
# OBJ: TLeafElement  GenJet.T    T[GenJet_] : 0 at: 0x2fc53d0
# OBJ: TLeafElement  GenJet.Mass Mass[GenJet_] : 0 at: 0x2fc69a0
# OBJ: TLeafElement  GenJet.DeltaEta DeltaEta[GenJet_] : 0 at: 0x2fc7fc0
# OBJ: TLeafElement  GenJet.DeltaPhi DeltaPhi[GenJet_] : 0 at: 0x2fc9630
# OBJ: TLeafElement  GenJet.Flavor   Flavor[GenJet_] : 0 at: 0x2fcac70
# OBJ: TLeafElement  GenJet.FlavorAlgo   FlavorAlgo[GenJet_] : 0 at: 0x2fcc2c0
# OBJ: TLeafElement  GenJet.FlavorPhys   FlavorPhys[GenJet_] : 0 at: 0x2fcd950
# OBJ: TLeafElement  GenJet.BTag BTag[GenJet_] : 0 at: 0x2fcef80
# OBJ: TLeafElement  GenJet.BTagAlgo BTagAlgo[GenJet_] : 0 at: 0x2fd05a0
# OBJ: TLeafElement  GenJet.BTagPhys BTagPhys[GenJet_] : 0 at: 0x2fd1c10
# OBJ: TLeafElement  GenJet.TauTag   TauTag[GenJet_] : 0 at: 0x2fd3250
# OBJ: TLeafElement  GenJet.TauWeight    TauWeight[GenJet_] : 0 at: 0x2fd48a0
# OBJ: TLeafElement  GenJet.Charge   Charge[GenJet_] : 0 at: 0x2fd5ef0
# OBJ: TLeafElement  GenJet.EhadOverEem  EhadOverEem[GenJet_] : 0 at: 0x2fd7540
# OBJ: TLeafElement  GenJet.NCharged NCharged[GenJet_] : 0 at: 0x2fd8bc0
# OBJ: TLeafElement  GenJet.NNeutrals    NNeutrals[GenJet_] : 0 at: 0x2fda240
# OBJ: TLeafElement  GenJet.NeutralEnergyFraction    NeutralEnergyFraction[GenJet_] : 0 at: 0x2fdb8d0
# OBJ: TLeafElement  GenJet.ChargedEnergyFraction    ChargedEnergyFraction[GenJet_] : 0 at: 0x2fdcf60
# OBJ: TLeafElement  GenJet.Beta Beta[GenJet_] : 0 at: 0x2fde590
# OBJ: TLeafElement  GenJet.BetaStar BetaStar[GenJet_] : 0 at: 0x2fdfbb0
# OBJ: TLeafElement  GenJet.MeanSqDeltaR MeanSqDeltaR[GenJet_] : 0 at: 0x2fe1230
# OBJ: TLeafElement  GenJet.PTD  PTD[GenJet_] : 0 at: 0x2fe2860
# OBJ: TLeafElement  GenJet.FracPt   FracPt[GenJet_] : 0 at: 0x2fe3e80
# OBJ: TLeafElement  GenJet.Tau  Tau[GenJet_] : 0 at: 0x2fe5470
# OBJ: TLeafElement  GenJet.SoftDroppedJet   SoftDroppedJet[GenJet_] : 0 at: 0x2fe6aa0
# OBJ: TLeafElement  GenJet.SoftDroppedSubJet1   SoftDroppedSubJet1[GenJet_] : 0 at: 0x2fe8130
# OBJ: TLeafElement  GenJet.SoftDroppedSubJet2   SoftDroppedSubJet2[GenJet_] : 0 at: 0x2fe97c0
# OBJ: TLeafElement  GenJet.TrimmedP4    TrimmedP4[GenJet_] : 0 at: 0x2feae50
# OBJ: TLeafElement  GenJet.PrunedP4 PrunedP4[GenJet_] : 0 at: 0x2fec5c0
# OBJ: TLeafElement  GenJet.SoftDroppedP4    SoftDroppedP4[GenJet_] : 0 at: 0x2fedd20
# OBJ: TLeafElement  GenJet.NSubJetsTrimmed  NSubJetsTrimmed[GenJet_] : 0 at: 0x2fef490
# OBJ: TLeafElement  GenJet.NSubJetsPruned   NSubJetsPruned[GenJet_] : 0 at: 0x2ff0b20
# OBJ: TLeafElement  GenJet.NSubJetsSoftDropped  NSubJetsSoftDropped[GenJet_] : 0 at: 0x2ff21b0
# OBJ: TLeafElement  GenJet.ExclYmerge23 ExclYmerge23[GenJet_] : 0 at: 0x2ff3840
# OBJ: TLeafElement  GenJet.ExclYmerge34 ExclYmerge34[GenJet_] : 0 at: 0x2ff4ed0
# OBJ: TLeafElement  GenJet.ExclYmerge45 ExclYmerge45[GenJet_] : 0 at: 0x2ff6560
# OBJ: TLeafElement  GenJet.ExclYmerge56 ExclYmerge56[GenJet_] : 0 at: 0x2ff7bf0
# OBJ: TLeafElement  GenJet.Constituents Constituents[GenJet_] : 0 at: 0x2ff9280
# OBJ: TLeafElement  GenJet.Particles    Particles[GenJet_] : 0 at: 0x2ffa950
# OBJ: TLeafElement  GenJet.Area Area[GenJet_] : 0 at: 0x2ffbfc0
# OBJ: TLeafI    GenJet_size GenJet_size : 0 at: 0x2ffe3b0
# OBJ: TLeafElement  GenMissingET_   GenMissingET_ : 0 at: 0x3000060
# OBJ: TLeafElement  GenMissingET.fUniqueID  fUniqueID[GenMissingET_] : 0 at: 0x2ffff80
# OBJ: TLeafElement  GenMissingET.fBits  fBits[GenMissingET_] : 0 at: 0x3001690
# OBJ: TLeafElement  GenMissingET.MET    MET[GenMissingET_] : 0 at: 0x3002d20
# OBJ: TLeafElement  GenMissingET.Eta    Eta[GenMissingET_] : 0 at: 0x30043b0
# OBJ: TLeafElement  GenMissingET.Phi    Phi[GenMissingET_] : 0 at: 0x3005a40
# OBJ: TLeafI    GenMissingET_size   GenMissingET_size : 0 at: 0x3007ec0
# OBJ: TLeafElement  Jet_    Jet_ : 0 at: 0x3009c70
# OBJ: TLeafElement  Jet.fUniqueID   fUniqueID[Jet_] : 0 at: 0x3009bd0
# OBJ: TLeafElement  Jet.fBits   fBits[Jet_] : 0 at: 0x300b240
# OBJ: TLeafElement  Jet.PT  PT[Jet_] : 0 at: 0x300c810
# OBJ: TLeafElement  Jet.Eta Eta[Jet_] : 0 at: 0x300dde0
# OBJ: TLeafElement  Jet.Phi Phi[Jet_] : 0 at: 0x300f3b0
# OBJ: TLeafElement  Jet.T   T[Jet_] : 0 at: 0x3010980
# OBJ: TLeafElement  Jet.Mass    Mass[Jet_] : 0 at: 0x3011f50
# OBJ: TLeafElement  Jet.DeltaEta    DeltaEta[Jet_] : 0 at: 0x3013520
# OBJ: TLeafElement  Jet.DeltaPhi    DeltaPhi[Jet_] : 0 at: 0x3014af0
# OBJ: TLeafElement  Jet.Flavor  Flavor[Jet_] : 0 at: 0x30160c0
# OBJ: TLeafElement  Jet.FlavorAlgo  FlavorAlgo[Jet_] : 0 at: 0x30176c0
# OBJ: TLeafElement  Jet.FlavorPhys  FlavorPhys[Jet_] : 0 at: 0x3018cf0
# OBJ: TLeafElement  Jet.BTag    BTag[Jet_] : 0 at: 0x301a2f0
# OBJ: TLeafElement  Jet.BTagAlgo    BTagAlgo[Jet_] : 0 at: 0x301b8c0
# OBJ: TLeafElement  Jet.BTagPhys    BTagPhys[Jet_] : 0 at: 0x301ce90
# OBJ: TLeafElement  Jet.TauTag  TauTag[Jet_] : 0 at: 0x301e460
# OBJ: TLeafElement  Jet.TauWeight   TauWeight[Jet_] : 0 at: 0x301fa50
# OBJ: TLeafElement  Jet.Charge  Charge[Jet_] : 0 at: 0x3021040
# OBJ: TLeafElement  Jet.EhadOverEem EhadOverEem[Jet_] : 0 at: 0x3022660
# OBJ: TLeafElement  Jet.NCharged    NCharged[Jet_] : 0 at: 0x3023c80
# OBJ: TLeafElement  Jet.NNeutrals   NNeutrals[Jet_] : 0 at: 0x3025270
# OBJ: TLeafElement  Jet.NeutralEnergyFraction   NeutralEnergyFraction[Jet_] : 0 at: 0x30268c0
# OBJ: TLeafElement  Jet.ChargedEnergyFraction   ChargedEnergyFraction[Jet_] : 0 at: 0x3027f50
# OBJ: TLeafElement  Jet.Beta    Beta[Jet_] : 0 at: 0x3029580
# OBJ: TLeafElement  Jet.BetaStar    BetaStar[Jet_] : 0 at: 0x302ab50
# OBJ: TLeafElement  Jet.MeanSqDeltaR    MeanSqDeltaR[Jet_] : 0 at: 0x302c180
# OBJ: TLeafElement  Jet.PTD PTD[Jet_] : 0 at: 0x302d7b0
# OBJ: TLeafElement  Jet.FracPt  FracPt[Jet_] : 0 at: 0x302ed80
# OBJ: TLeafElement  Jet.Tau Tau[Jet_] : 0 at: 0x3030350
# OBJ: TLeafElement  Jet.SoftDroppedJet  SoftDroppedJet[Jet_] : 0 at: 0x3031980
# OBJ: TLeafElement  Jet.SoftDroppedSubJet1  SoftDroppedSubJet1[Jet_] : 0 at: 0x3033010
# OBJ: TLeafElement  Jet.SoftDroppedSubJet2  SoftDroppedSubJet2[Jet_] : 0 at: 0x30346a0
# OBJ: TLeafElement  Jet.TrimmedP4   TrimmedP4[Jet_] : 0 at: 0x3035d20
# OBJ: TLeafElement  Jet.PrunedP4    PrunedP4[Jet_] : 0 at: 0x30373d0
# OBJ: TLeafElement  Jet.SoftDroppedP4   SoftDroppedP4[Jet_] : 0 at: 0x3038aa0
# OBJ: TLeafElement  Jet.NSubJetsTrimmed NSubJetsTrimmed[Jet_] : 0 at: 0x303a1d0
# OBJ: TLeafElement  Jet.NSubJetsPruned  NSubJetsPruned[Jet_] : 0 at: 0x303b860
# OBJ: TLeafElement  Jet.NSubJetsSoftDropped NSubJetsSoftDropped[Jet_] : 0 at: 0x303cef0
# OBJ: TLeafElement  Jet.ExclYmerge23    ExclYmerge23[Jet_] : 0 at: 0x303e580
# OBJ: TLeafElement  Jet.ExclYmerge34    ExclYmerge34[Jet_] : 0 at: 0x303fc10
# OBJ: TLeafElement  Jet.ExclYmerge45    ExclYmerge45[Jet_] : 0 at: 0x30412a0
# OBJ: TLeafElement  Jet.ExclYmerge56    ExclYmerge56[Jet_] : 0 at: 0x3042930
# OBJ: TLeafElement  Jet.Constituents    Constituents[Jet_] : 0 at: 0x3043fc0
# OBJ: TLeafElement  Jet.Particles   Particles[Jet_] : 0 at: 0x3045610
# OBJ: TLeafElement  Jet.Area    Area[Jet_] : 0 at: 0x3046c30
# OBJ: TLeafI    Jet_size    Jet_size : 0 at: 0x3049020
# OBJ: TLeafElement  Electron_   Electron_ : 0 at: 0x304ace0
# OBJ: TLeafElement  Electron.fUniqueID  fUniqueID[Electron_] : 0 at: 0x304ac00
# OBJ: TLeafElement  Electron.fBits  fBits[Electron_] : 0 at: 0x304c2e0
# OBJ: TLeafElement  Electron.PT PT[Electron_] : 0 at: 0x304d8e0
# OBJ: TLeafElement  Electron.Eta    Eta[Electron_] : 0 at: 0x304eeb0
# OBJ: TLeafElement  Electron.Phi    Phi[Electron_] : 0 at: 0x3050480
# OBJ: TLeafElement  Electron.T  T[Electron_] : 0 at: 0x3051a50
# OBJ: TLeafElement  Electron.Charge Charge[Electron_] : 0 at: 0x3053070
# OBJ: TLeafElement  Electron.EhadOverEem    EhadOverEem[Electron_] : 0 at: 0x30546f0
# OBJ: TLeafElement  Electron.Particle   Particle[Electron_] : 0 at: 0x3055d80
# OBJ: TLeafElement  Electron.IsolationVar   IsolationVar[Electron_] : 0 at: 0x3057410
# OBJ: TLeafElement  Electron.IsolationVarRhoCorr    IsolationVarRhoCorr[Electron_] : 0 at: 0x3058aa0
# OBJ: TLeafElement  Electron.SumPtCharged   SumPtCharged[Electron_] : 0 at: 0x305a130
# OBJ: TLeafElement  Electron.SumPtNeutral   SumPtNeutral[Electron_] : 0 at: 0x305b7c0
# OBJ: TLeafElement  Electron.SumPtChargedPU SumPtChargedPU[Electron_] : 0 at: 0x305ce50
# OBJ: TLeafElement  Electron.SumPt  SumPt[Electron_] : 0 at: 0x305e4b0
# OBJ: TLeafElement  Electron.D0 D0[Electron_] : 0 at: 0x305fab0
# OBJ: TLeafElement  Electron.DZ DZ[Electron_] : 0 at: 0x3061080
# OBJ: TLeafElement  Electron.ErrorD0    ErrorD0[Electron_] : 0 at: 0x30626b0
# OBJ: TLeafElement  Electron.ErrorDZ    ErrorDZ[Electron_] : 0 at: 0x3063d40
# OBJ: TLeafI    Electron_size   Electron_size : 0 at: 0x3066190
# OBJ: TLeafElement  Photon_ Photon_ : 0 at: 0x3067e60
# OBJ: TLeafElement  Photon.fUniqueID    fUniqueID[Photon_] : 0 at: 0x3067d80
# OBJ: TLeafElement  Photon.fBits    fBits[Photon_] : 0 at: 0x3069430
# OBJ: TLeafElement  Photon.PT   PT[Photon_] : 0 at: 0x306aa00
# OBJ: TLeafElement  Photon.Eta  Eta[Photon_] : 0 at: 0x306bfd0
# OBJ: TLeafElement  Photon.Phi  Phi[Photon_] : 0 at: 0x306d5a0
# OBJ: TLeafElement  Photon.E    E[Photon_] : 0 at: 0x306eb70
# OBJ: TLeafElement  Photon.T    T[Photon_] : 0 at: 0x3070140
# OBJ: TLeafElement  Photon.EhadOverEem  EhadOverEem[Photon_] : 0 at: 0x3071770
# OBJ: TLeafElement  Photon.Particles    Particles[Photon_] : 0 at: 0x3072e00
# OBJ: TLeafElement  Photon.IsolationVar IsolationVar[Photon_] : 0 at: 0x3074490
# OBJ: TLeafElement  Photon.IsolationVarRhoCorr  IsolationVarRhoCorr[Photon_] : 0 at: 0x3075b20
# OBJ: TLeafElement  Photon.SumPtCharged SumPtCharged[Photon_] : 0 at: 0x30771b0
# OBJ: TLeafElement  Photon.SumPtNeutral SumPtNeutral[Photon_] : 0 at: 0x3078840
# OBJ: TLeafElement  Photon.SumPtChargedPU   SumPtChargedPU[Photon_] : 0 at: 0x3079ed0
# OBJ: TLeafElement  Photon.SumPt    SumPt[Photon_] : 0 at: 0x307b500
# OBJ: TLeafElement  Photon.Status   Status[Photon_] : 0 at: 0x307caf0
# OBJ: TLeafI    Photon_size Photon_size : 0 at: 0x307ef00
# OBJ: TLeafElement  Muon_   Muon_ : 0 at: 0x3080b60
# OBJ: TLeafElement  Muon.fUniqueID  fUniqueID[Muon_] : 0 at: 0x3080ab0
# OBJ: TLeafElement  Muon.fBits  fBits[Muon_] : 0 at: 0x3082130
# OBJ: TLeafElement  Muon.PT PT[Muon_] : 0 at: 0x3083700
# OBJ: TLeafElement  Muon.Eta    Eta[Muon_] : 0 at: 0x3084cd0
# OBJ: TLeafElement  Muon.Phi    Phi[Muon_] : 0 at: 0x30862a0
# OBJ: TLeafElement  Muon.T  T[Muon_] : 0 at: 0x3087870
# OBJ: TLeafElement  Muon.Charge Charge[Muon_] : 0 at: 0x3088e40
# OBJ: TLeafElement  Muon.Particle   Particle[Muon_] : 0 at: 0x308a430
# OBJ: TLeafElement  Muon.IsolationVar   IsolationVar[Muon_] : 0 at: 0x308ba80
# OBJ: TLeafElement  Muon.IsolationVarRhoCorr    IsolationVarRhoCorr[Muon_] : 0 at: 0x308d110
# OBJ: TLeafElement  Muon.SumPtCharged   SumPtCharged[Muon_] : 0 at: 0x308e7a0
# OBJ: TLeafElement  Muon.SumPtNeutral   SumPtNeutral[Muon_] : 0 at: 0x308fe30
# OBJ: TLeafElement  Muon.SumPtChargedPU SumPtChargedPU[Muon_] : 0 at: 0x30914c0
# OBJ: TLeafElement  Muon.SumPt  SumPt[Muon_] : 0 at: 0x3092af0
# OBJ: TLeafElement  Muon.D0 D0[Muon_] : 0 at: 0x30940c0
# OBJ: TLeafElement  Muon.DZ DZ[Muon_] : 0 at: 0x3095690
# OBJ: TLeafElement  Muon.ErrorD0    ErrorD0[Muon_] : 0 at: 0x3096c60
# OBJ: TLeafElement  Muon.ErrorDZ    ErrorDZ[Muon_] : 0 at: 0x3098230
# OBJ: TLeafI    Muon_size   Muon_size : 0 at: 0x309a620
# OBJ: TLeafElement  MissingET_  MissingET_ : 0 at: 0x309c2d0
# OBJ: TLeafElement  MissingET.fUniqueID fUniqueID[MissingET_] : 0 at: 0x309c1f0
# OBJ: TLeafElement  MissingET.fBits fBits[MissingET_] : 0 at: 0x309d8f0
# OBJ: TLeafElement  MissingET.MET   MET[MissingET_] : 0 at: 0x309ef30
# OBJ: TLeafElement  MissingET.Eta   Eta[MissingET_] : 0 at: 0x30a0540
# OBJ: TLeafElement  MissingET.Phi   Phi[MissingET_] : 0 at: 0x30a1b50
# OBJ: TLeafI    MissingET_size  MissingET_size : 0 at: 0x30a3f60
# OBJ: TLeafElement  ScalarHT_   ScalarHT_ : 0 at: 0x30a5c50
# OBJ: TLeafElement  ScalarHT.fUniqueID  fUniqueID[ScalarHT_] : 0 at: 0x30a5b70
# OBJ: TLeafElement  ScalarHT.fBits  fBits[ScalarHT_] : 0 at: 0x30a7250
# OBJ: TLeafElement  ScalarHT.HT HT[ScalarHT_] : 0 at: 0x30a8850
# OBJ: TLeafI    ScalarHT_size   ScalarHT_size : 0 at: 0x30aac40

