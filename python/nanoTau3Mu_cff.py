from __future__ import print_function
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.NanoAOD.met_cff import *
from PhysicsTools.XNano.trgbits_cff import * # modified

## for gen and trigger muon
from PhysicsTools.XNano.genparticlesT3m_cff import * # define new
from PhysicsTools.XNano.particlelevelT3m_cff import * # define new
from PhysicsTools.XNano.triggerObjectsTau3Mu_cff import * #define new
from PhysicsTools.XNano.muonsTau3mu_cff import * # define new 

## filtered input collections
#from PhysicsTools.XNano.tracksX_cff import *

## W collections
from PhysicsTools.XNano.Wnu_Tau3Mu import * #define new


nanoSequenceOnlyFullSim = cms.Sequence(triggerObjectTau3MuTables + l1bits)

nanoSequence = cms.Sequence(nanoMetadata + 
                            cms.Sequence(vertexTask) + cms.Sequence(metTablesTask) +           
                            cms.Sequence(globalTablesTask) + cms.Sequence(vertexTablesTask) +
                            triggerObjectTau3MuTables + l1bits)

nanoSequenceMC = cms.Sequence(particleLevelT3mSequence + genParticleT3mSequence + cms.Sequence(metMCTask) + 
                              cms.Sequence(globalTablesMCTask) + cms.Sequence(genWeightsTableTask) + genParticleT3mTables + lheInfoTable)

def nanoAOD_customizeMuonTriggerTau3Mu(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + muonT3mSequence + muonT3mTables)
    return process

#def nanoAOD_customizeTrackFilteredX(process):
#    process.nanoSequence = cms.Sequence( process.nanoSequence + tracksXSequence + tracksXTables)
#    return process

def nanoAOD_customizeTriggerBitsTau3Mu(process):
    process.nanoSequence = cms.Sequence( process.nanoSequence + trgTables)
    return process

def nanoAOD_customizeWnuTau3Mu(process):
    #process.nanoB0ToK0XSequence = cms.Sequence( Tau3MuSequence + Tau3MuTableSequence )
    process.nanoWnuTau3MuSequence = cms.Sequence( Tau3MuSequence + Tau3MuTableSequence + TauPlusMetSequence + TauPlusMetTableSequence)
    return process


from FWCore.ParameterSet.MassReplace import massSearchReplaceAnyInputTag
def nanoAOD_customizeMC(process):
    for name, path in process.paths.iteritems():
        # replace all the non-match embedded inputs with the matched ones
        massSearchReplaceAnyInputTag(path, 'muonTrgSelector:SelectedMuons', 'selectedMuonsMCMatchEmbedded')
        #massSearchReplaceAnyInputTag(path, 'tracksX:SelectedTracks', 'tracksXMCMatchEmbedded')

        # modify the path to include mc-specific info
        path.insert(0, nanoSequenceMC)
        path.replace(process.muonT3mSequence, process.muonT3mMC)
        #path.replace(process.tracksXSequence, process.tracksXMC)
