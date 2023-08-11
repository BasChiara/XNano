import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.triggerObjects_cff import *

triggerObjectTau3MuTable = cms.EDProducer("TriggerObjectTableProducerTau3Mu",
    name= cms.string("TrigObj"),
    src = cms.InputTag("unpackedPatTrigger"),
    selections = cms.VPSet(
        # selezione sui muoni
        cms.PSet( 
            name = cms.string("Muon"),
            id = cms.int32(13),
            sel = cms.string("type(83) && pt > 5 && coll('hltIterL3MuonCandidates')"), 
            qualityBits = cms.string("filter('hlt')"), qualityBitsDoc = cms.string("1 = Muon filters"),
        ),
    ),
)

triggerObjectTau3MuTables = cms.Sequence( unpackedPatTrigger + triggerObjectTau3MuTable )
