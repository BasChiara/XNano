import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

Path2022=["HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_v6","HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_v6","HLT_DoubleMu4_3_LowMass"]
Path = Path2022

trgTable = cms.EDProducer( "TrgBitTableProducer",
                          hltresults = cms.InputTag("TriggerResults::HLT"),
                          paths      = cms.vstring(Path),

                          # interesting paths for 2022
                          #paths      = cms.vstring(
                          #                   "HLT_Dimuon25_Jpsi",    
                          #                   "HLT_DoubleMu4_JpsiTrk_Bc",
                          #                   "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
                          #                   "HLT_DoubleMu4_LowMass_Displaced",
                          #                   "HLT_DoubleMu4_MuMuTrk_Displaced",
                          #                   "HLT_DoubleMu4_3_LowMass"
                          #                    ),
)

trgTables = cms.Sequence(trgTable)



