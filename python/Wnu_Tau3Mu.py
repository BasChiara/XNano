import FWCore.ParameterSet.Config as cms
from PhysicsTools.XNano.common_cff import *


########## inputs preparation ################

# Tau -> 3mu
muonTripletForTau3Mu = cms.EDProducer('TriMuonBuilder',
    src = cms.InputTag('triMuonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('triMuonTrgSelector', 'SelectedTransientMuons'),
    beamSpot   = cms.InputTag("offlineBeamSpot"),
    lep1Selection = cms.string('isMediumMuon && ((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
    lep2Selection = cms.string('isMediumMuon && ((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
    lep3Selection = cms.string('isMediumMuon && ((abs(eta) <= 1.2 && pt > 3.5) || (abs(eta) > 1.2 && abs(eta) < 2.4 && pt > 2.0))'),
    preVtxSelection = cms.string('mass() < 3 && abs(charge()) == 1'), # selection for tau candidates pre-fit
    postVtxSelection =  cms.string('userInt("vtx_isValid")'),
    #preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
    #                             '&& mass() > 1 && charge() == 0 && userFloat("lep_deltaR") > 0.02'),
    #postVtxSelection =  cms.string('userFloat("sv_prob") > 0.01'
    #                             '&& userFloat("fitted_mass") > 2.50 && userFloat("fitted_mass") < 4.00'),
)



################################### Tables #####################################

Tau3MuTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag('muonTripletForTau3Mu', 'SelectedTriMuons'),
    cut = cms.string(""),
    name = cms.string("TauTo3Mu"),
    doc = cms.string("Tau Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        mu1_idx = uint('l1_idx'),
        mu2_idx = uint('l2_idx'),
        mu3_idx = uint('l3_idx'),
        #Tau_charge  = uint('Tau_charge'),

        vtx_prob = ufloat("vtx_prob"),
        vtx_isValid = uint("vtx_isValid"),

        fitted_wovc_pt   = ufloat('fitted_wovc_pt'),
        fitted_wovc_eta  = ufloat('fitted_wovc_eta'),
        fitted_wovc_phi  = ufloat('fitted_wovc_phi'),
        fitted_wovc_mass = ufloat("fitted_wovc_mass"),

        fitted_vc_pt   = ufloat('fitted_vc_pt'),
        fitted_vc_eta  = ufloat('fitted_vc_eta'),
        fitted_vc_phi  = ufloat('fitted_vc_phi'),
        fitted_vc_mass = ufloat("fitted_vc_mass"),
        
        mu1_pt = ufloat("mu1_pt"),
        mu1_eta = ufloat("mu1_eta"),
        mu1_phi = ufloat("mu1_phi"),
        mu1_dr = ufloat("mu1_dr"),
        mu1_trackQuality = uint("mu1_trackQuality"),
        mu2_pt = ufloat("mu2_pt"),
        mu2_eta = ufloat("mu2_eta"),
        mu2_phi = ufloat("mu2_phi"),
        mu2_dr = ufloat("mu2_dr"),
        mu2_trackQuality = uint("mu2_trackQuality"),
        mu3_pt = ufloat("mu3_pt"),
        mu3_eta = ufloat("mu3_eta"),
        mu3_phi = ufloat("mu3_phi"),
        mu3_dr = ufloat("mu3_dr"),
        mu3_trackQuality = uint("mu3_trackQuality"),
        mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = uint("mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"),
        mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = uint("mu1_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"),
        #mu1_fired_DoubleMu4_3_LowMass = uint("mu1_fired_DoubleMu4_3_LowMass"),


        mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = uint("mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"),
        mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = uint("mu2_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"),
        #mu2_fired_DoubleMu4_3_LowMass = uint("mu2_fired_DoubleMu4_3_LowMass"),

        mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1 = uint("mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1"),
        mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15 = uint("mu3_fired_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15"),
        #mu3_fired_DoubleMu4_3_LowMass = uint("mu3_fired_DoubleMu4_3_LowMass"),
    )
)


CountMuonTriplets = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag('muonTripletForTau3Mu', 'SelectedTriMuons'),
)    

########################### Sequencies  ############################

Tau3MuSequence = cms.Sequence(
    (muonTripletForTau3Mu * CountMuonTriplets)
)

Tau3MuTableSequence = cms.Sequence( Tau3MuTable )
