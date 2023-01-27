# https://root-forum.cern.ch/t/cannot-perform-both-dot-product-and-scalar-multiplication-on-tvector2-in-pyroot/28207
import ROOT
def fixtv2mul():
   ROOT.TVector2(1, 1) * ROOT.TVector2(2, 2)
   mvv = ROOT.TVector2.__mul__
   del ROOT.TVector2.__mul__
   ROOT.TVector2(1, 1) * 5
   mvs = ROOT.TVector2.__mul__
   def fixmul(self, other, mvv=mvv, mvs=mvs):
      if isinstance(other, self.__class__):
         return mvv(self, other)
      return mvs(self, other)
   ROOT.TVector2.__mul__ = fixmul

fixtv2mul()

def fixtv3mul():
   ROOT.TVector3(1, 1, 1) * ROOT.TVector3(2, 2, 2)
   mvv = ROOT.TVector3.__mul__
   del ROOT.TVector3.__mul__
   ROOT.TVector3(1, 1, 1) * 5
   mvs = ROOT.TVector3.__mul__
   def fixmul(self, other, mvv=mvv, mvs=mvs):
      if isinstance(other, self.__class__):
         return mvv(self, other)
      return mvs(self, other)
   ROOT.TVector3.__mul__ = fixmul

fixtv3mul()
