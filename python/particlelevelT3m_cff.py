import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.particlelevel_cff import *

particleLevelT3mSequence = cms.Sequence(mergedGenParticles + genParticles2HepMC + particleLevel)
