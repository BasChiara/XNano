import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.genparticles_cff import *

# Start with merged particles (pruned + packed),
# where pruned contain K* states, but not final states, 
# and packed contain final states (K pi).
# then you save also final states (granddaughters)
finalGenParticlesT3m = finalGenParticles.clone(
  src = cms.InputTag("mergedGenParticles"),
  select = cms.vstring(
	"drop *",
        "keep++ (abs(pdgId) == 15)",  #keep all tau(=15)
   )
)

genParticleT3mTable = genParticleTable.clone(
  src = cms.InputTag("finalGenParticlesT3m"),
  variables = cms.PSet(
      genParticleTable.variables,
      vx = Var("vx()", float, doc="x coordinate of the production vertex position, in cm", precision=10),
      vy = Var("vy()", float, doc="y coordinate of the production vertex position, in cm", precision=10),
      vz = Var("vz()", float, doc="z coordinate of the production vertex position, in cm", precision=10),
  )
)

genParticleT3mSequence = cms.Sequence(finalGenParticlesT3m)
genParticleT3mTables = cms.Sequence(genParticleT3mTable)

