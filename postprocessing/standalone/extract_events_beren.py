import ROOT
import os
import glob
from array import array
from DataFormats.FWLite import Events, Handle
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.Config as edm

# save the arrays as uproot awkward arrays
import uproot3
import awkward as ak
import numpy as np 

ROOT.gROOT.SetBatch(True)

totalEvents = 0

# define the location of files 
root_files_directory = "/eos/purdue/store/user/lingqian/nanogen_eft/TT01j2l/"
#root_files = glob.glob(os.path.join(root_files_directory, "*.root"))
root_files = glob.glob(os.path.join(root_files_directory, "GEN_LO_01j_102X_1.root"))

if not os.path.exists("MiniTreeOutput"): os.makedirs("MiniTreeOutput")

h_lep_pt = ROOT.TH1F("h_lep_pt", "Lepton pT; pT (GeV);Events", 200, 0, 1000)
h_top_pt = ROOT.TH1F("h_top_pt", "Top Quark pT; pT (GeV);Events", 600, 0, 3000)
h_antitop_pt = ROOT.TH1F("h_antitop_pt", "Anti-Top Quark p_{T}; p_{T} [GeV];Events", 600, 0, 3000)
h_lep_eta = ROOT.TH1F("h_lep_eta", "eta; #eta;Events", 20, -5, 5)
h_lep_phi = ROOT.TH1F("h_lep_phi", "Azimuthal Angle; #phi;Events", 20, -ROOT.TMath.Pi(), ROOT.TMath.Pi())
h_ttbar_mass = ROOT.TH1F("h_ttbar_mass", "Invariant Mass; M (GeV);Events", 700, 0, 7000)
h_parton_multiplicity = ROOT.TH1F("h_parton_multiplicity", "Jet Multiplicity; N_{jets};Events", 20, 0, 100)
h_neu_pt = ROOT.TH1F("hMET", "MET;MET (GeV);Events", 40, 0, 200)
h_bquark_pt = ROOT.TH1F("hbquarkPt", "b-quark pT;pT (GeV);Events", 200, 0, 1000)
h_bquark_eta = ROOT.TH1F("hbquarkEta", "b-quark #eta;#eta;Events", 20, -5, 5)
h_angle_top_antitop = ROOT.TH1F("h_angle", "Angle between top and antitop;Angle (radians);Events", 50, 0, ROOT.TMath.Pi())

h_decayChannel = ROOT.TH1F("h_decayChannel", "Top Decay Channels; Channel; Events", 2, 0, 2)
h_decayChannel.GetXaxis().SetBinLabel(1, "t -> W+b")
h_decayChannel.GetXaxis().SetBinLabel(2, "Other")


h_topMultiplicity = ROOT.TH1F("h_topMultiplicity", "Top Multiplicity; N_{top};Events", 5, 0, 5)
# h_antitopMultiplicity = ROOT.TH1F("h_antitopMultiplicity", "Anti-Top Multiplicity; N_{antitop};Events", 5, 0, 5)

h_missingParticles = ROOT.TH1F("h_missingParticles", "Missing Particles; Particle Type; Events", 4, 0, 4)
h_missingParticles.GetXaxis().SetBinLabel(1, "No Top")
h_missingParticles.GetXaxis().SetBinLabel(2, "No Anti-Top")
h_missingParticles.GetXaxis().SetBinLabel(3, "No Top & No Anti-Top")

bin_edges = [-16.5, -14.5, -12.5, -10.5, 10.5, 12.5, 14.5, 16.5]
h_leptonFlavor = ROOT.TH1F("h_leptonFlavor", "Lepton Flavor; PDG ID;Events", len(bin_edges)-1, array('d', bin_edges))

h_leptonFlavor.GetXaxis().SetBinLabel(1, "tau-")
h_leptonFlavor.GetXaxis().SetBinLabel(2, "muon-")
h_leptonFlavor.GetXaxis().SetBinLabel(3, "electron-")
h_leptonFlavor.GetXaxis().SetBinLabel(5, "electron+")
h_leptonFlavor.GetXaxis().SetBinLabel(6, "muon+")
h_leptonFlavor.GetXaxis().SetBinLabel(7, "tau+")


h_nonTopMotherJets = ROOT.TH1F("h_nonTopMotherJets", "Jets without Top as Mother; Count;Events", 10, 0, 50)
h_jetMultiplicity = ROOT.TH1F("h_jetMultiplicity", "Number of Jets per Event", 10, 0, 50)

h_topMother = ROOT.TH1F("h_topMother", "Mother of Top Quarks; Mother; Events", 3, 0, 3)
h_topMother.GetXaxis().SetBinLabel(1, "qq")
h_topMother.GetXaxis().SetBinLabel(2, "gg")
h_topMother.GetXaxis().SetBinLabel(3, "Other")

h_motherPdgId = ROOT.TH1F("h_motherPdgId", "PDG ID of Top's Mother;PDG ID;Counts", 23, -6, 22)

def analyze(filename):
    
    print(filename)
    
    events = Events(filename)
    global totalEvents 
    totalEvents += events.size()
    print("Number of events in file:", events.size())
    
    handle = Handle('vector<reco::GenParticle>')
    genJetsHandle = Handle('vector<reco::GenJet>')
    
    relevant_pdgIds = {12,14,16,24,1,2,3,4,5,6,21,11,13,15}
    
    pdgIds = []  

    lep_pt = []
    lep_eta = []
    lep_phi = []    
    lep_mass = [] 
    
    b_pt = []
    b_eta = []
    b_phi = []
    b_mass = [] 
    
    neu_pt = []
    neu_eta = []
    neu_phi = []
    neu_mass = []  
    
    top_pt = [] 
    top_eta = []
    top_phi = []
    top_mass = [] 
    
    atop_pt = []
    atop_eta = []
    atop_phi = []
    atop_mass = [] 
    
    ttbar_mass = []
    ttbar_costheta = []
   
    for event in events:
        # GenParticles
        event.getByLabel("genParticles", handle)
        particles = handle.product()

        # particles = [p for p in particles if abs(p.pdgId()) in relevant_pdgIds]

        # GenJets
        event.getByLabel("ak4GenJets", genJetsHandle)
        jets = genJetsHandle.product()
        
        # print("Number of particles in event:", len(particles))
        
        tops = []
        bquarks = []
        leptons = []
        neutrinos = []
        partons = []
        non_top_mother_jet_count_j = []
        
        top = None
        antitop = None
        
        top_count = 0
        antitop_count = 0
        not_top = 0
        jet_count = 0
        non_top_mother_jet_count = 0
                
        for jet in jets:
            jet_count +=1
            
            mother = jet.mother()
            if mother and mother.pdgId() not in [6, -6]:
                non_top_mother_jet_count_j.append(jet)
        
        h_jetMultiplicity.Fill(len(jets))  
        h_nonTopMotherJets.Fill(len(non_top_mother_jet_count_j))
                
        # Iterates over each particle in the list. checks if the particle is a top of antitop. 
        # finds the mothers of the tops. makes sure that the mothers are not tops.for each particle that is top, it looks for W boson and bquark daughters
        #if both W and b quark daugthers are found, the particle is appended to the tops list. 
        #we are ensuring that we're only counting hard scatter tops and not tops from decay chain of another top. 
        
        #each top has a mother which isn't a top.
        # Each top decays into W+b.
        # iff W decays into a lepton, that lepton's pt should be greater than 30, eta is greater than 2.4. If these conditions aren't met, the loop continues to the next particle without processing the current one. 
        
        for particle in particles:
            pdgId = particle.pdgId()
            
            #if (abs(pdgId)==6): print(pdgId)  
            
            pdgIds.append(particle.pdgId())
            
            # Tops
            if abs(pdgId) == 6: 
                # Check the mothers of the top and antitop
                mother1 = particle.mother(0)
                mother2 = particle.numberOfMothers() > 1 and particle.mother(1) or None
                
                #print('mother1', mother1)
                #print('mother2', mother2)
                
                # Check if the mothers are top
                #if any([mom and abs(mom.pdgId()) == 6 for mom in [mother1, mother2]]):
                if mother1:
                    h_motherPdgId.Fill(mother1.pdgId())
                if mother2:
                    h_motherPdgId.Fill(mother2.pdgId())

                # Checking for W and b quark daughters
                w_quark_daughter = None
                b_quark_daughter = None
                for i in range(particle.numberOfDaughters()):
                    daughter = particle.daughter(i)
                    if abs(daughter.pdgId()) == 24:
                        w_quark_daughter = daughter
                        
                    elif abs(daughter.pdgId()) == 5:
                        b_quark_daughter = daughter

                if not w_quark_daughter or not b_quark_daughter:
                    continue

                has_high_pt_lepton = False
                for j in range(w_quark_daughter.numberOfDaughters()):
                    lepton_candidate = w_quark_daughter.daughter(j)

                    #print(lepton_candidate) 

                    # if abs(lepton_candidate.pdgId()) in [11, 13] and lepton_candidate.pt() > 30 and abs(lepton_candidate.eta()) < 2.4:
                    if abs(lepton_candidate.pdgId()) in [11, 13]:
                        has_high_pt_lepton = True
                        lepton = lepton_candidate 
                        h_lep_pt.Fill(lepton.pt())
                        h_lep_eta.Fill(lepton.eta())
                        h_lep_phi.Fill(lepton.phi())
                        h_leptonFlavor.Fill(lepton.pdgId())

                        lep_pt.append(lepton.pt())
                        lep_eta.append(lepton.eta())
                        lep_phi.append(lepton.phi())
                        #lep_mass.append(lepton.M()) 

                if not has_high_pt_lepton:
                    continue

                tops.append(particle)
                h_decayChannel.Fill(0)  # t -> W+b

                if particle.pdgId() == 6:
             
                    top_count += 1
                    #print("hey, there's a top and there are now " + str(top_count) + " top quarks") 
                    top_pt.append(particle.pt()) 
                    top_eta.append(particle.eta())
                    top_phi.append(particle.phi()) 

                    h_top_pt.Fill(particle.pt())
                    top = ROOT.TLorentzVector()
                    top.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.energy())
                elif particle.pdgId() == -6:
                    antitop_count += 1
                    h_antitop_pt.Fill(particle.pt())
                    atop_pt.append(particle.pt()) 
                    atop_eta.append(particle.eta())
                    atop_phi.append(particle.phi())
                    
                    antitop = ROOT.TLorentzVector()
                    antitop.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.energy())

                # for j in range(w_quark_daughter.numberOfDaughters()):
                #     lepton_candidate = w_quark_daughter.daughter(j)
                #     if abs(lepton_candidate.pdgId()) in [11, 13]:
                #         lepton = lepton_candidate
                #         h_lep_pt.Fill(lepton.pt())
                #         h_lep_eta.Fill(lepton.eta())
                #         h_lep_phi.Fill(lepton.phi())
                #         h_leptonFlavor.Fill(lepton.pdgId())

                if mother1 and mother2 and set([abs(mother1.pdgId()), abs(mother2.pdgId())]) == {21}:
                    h_topMother.Fill(1)  # gg
                elif mother1 and mother2 and set([abs(mother1.pdgId()), abs(mother2.pdgId())]).issubset({1,2,3,4,5}):
                    h_topMother.Fill(0)  # qq
                else:
                    h_topMother.Fill(2)  # Other
                
            # Partons
            if abs(pdgId) in [1, 2, 3, 4, 5, 6, 21]:
                partons.append(particle) 
            
            if abs(pdgId) == 5:
                bquarks.append(particle)
                
                b_vector = ROOT.TLorentzVector()
                b_vector.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.energy())
                
                if pdgId == 5:
                    h_bquark_pt.Fill(b_vector.Pt())
                    h_bquark_eta.Fill(b_vector.Eta())
                    
                    b_pt.append(b_vector.Pt())
                    b_eta.append(b_vector.Eta())
                    b_phi.append(b_vector.Phi())
                    #b_mass.append(b_vector.M())
        
            if abs(pdgId) in [12, 14, 16]:
                neutrino = ROOT.TLorentzVector()
                neutrino.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.energy())
                h_neu_pt.Fill(neutrino.Pt())
                
                neu_pt.append(neutrino.Pt())
                neu_eta.append(neutrino.Eta())
                neu_phi.append(neutrino.Phi())
                #neu_mass.append(neutrino.M())
                
        h_parton_multiplicity.Fill(len(partons))
        h_topMultiplicity.Fill(len(tops))
        # h_antitopMultiplicity.Fill(antitop_count)
        
        if top_count == 0:
            h_missingParticles.Fill(0)  # Filling "no top" bin

        if antitop_count == 0:
            h_missingParticles.Fill(1)  # Filling "no antitop" bin

        if top_count == 0 and antitop_count == 0:
            h_missingParticles.Fill(2)  # Filling "no top and antitop" bin
  
        if top and antitop:
            print("hey, here's a pair of tt", event)
            ttbar = top + antitop
            
            h_ttbar_mass.Fill(ttbar.M())
            h_angle_top_antitop.Fill(top.Angle(antitop.Vect()))
            
            ttbar_mass.append(ttbar.M())
            ttbar_costheta.append(top.Angle(antitop.Vect()))

    postfix = filename.replace(root_files_directory+"GEN_LO_01j_102X_", "")
    
    print(root_files_directory+"/GEN_LO_01j_102X_")
    #print("postfix", postfix) 
    newfile = uproot.recreate("MiniTreeOutput/minitree_"+postfix)     
    
    print("lep_pt", np.array(lep_pt))
    # print "lep_eta", lep_eta
    # print "lep_phi", lep_phi
    # print "b_pt", b_pt
    # print "b_eta", b_eta
    # print "b_phi", b_phi
    # print "neu_pt", neu_pt
    # print "neu_eta", neu_eta
    # print "neu_phi", neu_phi 
    # print "top_pt", top_pt
    # print "top_eta", top_eta
    # print "top_phi", top_phi
    # print "ttbar_mass", ttbar_mass
            
    newfile["tree"] = uproot3.newtree({"branch": "int32"})
    newfile["tree"].extend({"branch": np.array([1, 2, 3, 4, 5])})

total_files = len(root_files)
processed_files = 0

# newfile["pre-selection"] = {
    
#     "gen_lep_pt": lep_pt,
#     "gen_lep_eta": lep_eta,
#     "gen_lep_phi": lep_phi,
#     #"gen_lep_mass": lep_mass,
#     "gen_b_pt": b_pt,
#     "gen_b_eta": b_eta,
#     "gen_b_phi": b_phi,
#     "gen_b_mass": b_mass,
#     "gen_top_pt": top_pt,
#     "gen_top_eta": top_eta,
#     "gen_top_phi": top_phi,
#     "gen_top_mass": top_mass,    
#     "gen_neu_pt": neu_pt,
#     "gen_neu_eta": neu_eta,
#     "gen_neu_phi": neu_phi,
#     "gen_neu_mass": neu_mass, 
#     "gen_ttbar_mass": ttbar_mass
#     #"gen_ttbar_costheta": ttbar_costheta, 
    
#     }

for root_file in root_files:
    print('Analyzing: ',processed_files + 1, 'out of',total_files)
    analyze(root_file)
    processed_files += 1

canvas = ROOT.TCanvas("canvas", "Analysis Plots", 4000, 4000)
canvas.Divide(4, 5)

canvas.cd(1)
h_lep_pt.Draw()
canvas.cd(2)
h_top_pt.Draw()
canvas.cd(3)
h_antitop_pt.Draw()
canvas.cd(4)
h_lep_eta.Draw()
canvas.cd(5)
h_lep_phi.Draw()
canvas.cd(6)
h_ttbar_mass.Draw()
canvas.cd(7)
h_decayChannel.Draw()
canvas.cd(8)
h_neu_pt.Draw()
canvas.cd(9)
h_leptonFlavor.Draw()
canvas.cd(10)
h_bquark_pt.Draw()
canvas.cd(11)
h_bquark_eta.Draw()
canvas.cd(12)
h_angle_top_antitop.Draw()
canvas.cd(13)
h_parton_multiplicity.Draw()
canvas.cd(14)
h_nonTopMotherJets.Draw()
canvas.cd(15)
h_topMultiplicity.Draw()
canvas.cd(16)
h_jetMultiplicity.Draw()
canvas.cd(17)
h_topMother.Draw()
canvas.cd(18)
h_missingParticles.Draw()
canvas.cd(19)
h_motherPdgId.Draw()


canvas.SaveAs("allPlots.png")

c_lep_pt = ROOT.TCanvas("c_lep_pt", "Lepton pT Distribution", 800, 600)
h_lep_pt.Draw()
ROOT.gPad.SetLogy(1)
c_lep_pt.SaveAs("lep_ptDistribution.png")

c_top_pt = ROOT.TCanvas("c_top_pt", "Top Quark pT Distribution", 800, 600)
h_top_pt.Draw()
ROOT.gPad.SetLogy(1)
c_top_pt.SaveAs("top_ptDistribution.png")

c_antitop_pt = ROOT.TCanvas("c_antitop_pt", "Anti-Top Quark pT Distribution", 800, 600)
h_antitop_pt.Draw()
ROOT.gPad.SetLogy(1)
c_antitop_pt.SaveAs("antitop_ptDistribution.png")

c_eta = ROOT.TCanvas("c_eta", "Lepton Eta Distribution", 800, 600)
h_lep_eta.Draw()
c_eta.SaveAs("etaDistribution.png")

c_phi = ROOT.TCanvas("c_phi", "Lepton Azimuthal Angle Distribution", 800, 600)
h_lep_phi.Draw()
c_phi.SaveAs("phiDistribution.png")

c_ttbar_mass = ROOT.TCanvas("c_ttbar_mass", "Invariant Mass Distribution", 800, 600)
h_ttbar_mass.Draw()
ROOT.gPad.SetLogy(1)
c_ttbar_mass.SaveAs("ttbar_massDistribution.png")

c_decay = ROOT.TCanvas("c_decay", "Decay Channel Canvas", 800, 600)
h_decayChannel.Draw()
c_decay.SaveAs("topDecayChannel.png")

c_neu_pt = ROOT.TCanvas("cMET", "MET Distribution", 800, 600)
h_neu_pt.Draw()
ROOT.gPad.SetLogy(1)
c_neu_pt.SaveAs("METDistribution.png")

c_leptonFlavor = ROOT.TCanvas("c_leptonFlavor", "Lepton Flavor Distribution", 800, 600)
h_leptonFlavor.Draw()
c_leptonFlavor.SaveAs("leptonFlavorDistribution.png")

c_bquark_pt = ROOT.TCanvas("cbquarkPt", "b-quark pT Distribution", 800, 600)
h_bquark_pt.Draw()
ROOT.gPad.SetLogy(1)
c_bquark_pt.SaveAs("bquarkPtDistribution.png")

c_bquark_eta = ROOT.TCanvas("cbquarkEta", "b-quark Eta Distribution", 800, 600)
h_bquark_eta.Draw()
c_bquark_eta.SaveAs("bquarkEtaDistribution.png")

c_angle = ROOT.TCanvas("cangle", "Angle between top and antitop", 800, 600)
h_angle_top_antitop.Draw()
c_angle.SaveAs("angleTopAntitop.png")

c_parton_multiplicity = ROOT.TCanvas("c_parton_multiplicity", "Parton Multiplicity Distribution", 800, 600)
h_parton_multiplicity.SetFillColor(ROOT.kBlue - 10)
h_parton_multiplicity.SetLineColor(ROOT.kBlue)
h_parton_multiplicity.Draw()
ROOT.gPad.SetLogy(1)
c_parton_multiplicity.SaveAs("parton_multiplicityDistribution.png")

c_nonTopMotherJets = ROOT.TCanvas("c_nonTopMotherJets", "Jets without Top as Mother", 800, 600)
h_nonTopMotherJets.SetFillColor(ROOT.kBlue - 10)
h_nonTopMotherJets.SetLineColor(ROOT.kBlue)
h_nonTopMotherJets.Draw()
ROOT.gPad.SetLogy(1)
c_nonTopMotherJets.SaveAs("nonTopMotherJets.png")

# c_antitopMultiplicity = ROOT.TCanvas("c_antitopMultiplicity", "Anti-Top Multiplicity Distribution", 800, 600)
# h_antitopMultiplicity.SetFillColor(ROOT.kBlue - 10)
# h_antitopMultiplicity.SetLineColor(ROOT.kBlue)
# h_antitopMultiplicity.Draw()
# c_antitopMultiplicity.SaveAs("antitopMultiplicityDistribution.png")

c_topMultiplicity = ROOT.TCanvas("c_topMultiplicity", "Top Multiplicity Distribution", 800, 600)
h_topMultiplicity.SetFillColor(ROOT.kBlue - 10)
h_topMultiplicity.SetLineColor(ROOT.kBlue)
h_topMultiplicity.Draw()
ROOT.gPad.SetLogy(1)
c_topMultiplicity.SaveAs("topMultiplicityDistribution.png")

c_jetMultiplicity = ROOT.TCanvas("c_jetMultiplicity", "Number of Jets per Event", 800, 600)
h_jetMultiplicity.SetFillColor(ROOT.kBlue - 10)
h_jetMultiplicity.SetLineColor(ROOT.kBlue)
h_jetMultiplicity.Draw()
ROOT.gPad.SetLogy(1)
c_jetMultiplicity.SaveAs("jetMultiplicity.png")

c_topMother = ROOT.TCanvas("c_topMother", "Mothers of the top quark", 800, 600)
h_topMother.Draw()
c_topMother.SaveAs("topMother.png")

c_missingParticles = ROOT.TCanvas("c_missingParticles", "Missing Particles", 800, 600)
h_missingParticles.Draw()
ROOT.gPad.SetLogy(1)
c_missingParticles.SaveAs("missingpart.png")

c_motherPdgId = ROOT.TCanvas("c_motherPdgId", "PDG ID of Top's Mother", 800,600)
h_motherPdgId.Draw()
c_motherPdgId.SaveAs("motherPDG.png")

print("Total number of events:", totalEvents)
