from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

options = VarParsing('python')

options.register('isMC', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

options.setDefault('maxEvents', 10)
options.setDefault('tag', '10215')
options.parseArguments()

globaltag = '102X_dataRun2_v11' if not options.isMC else '94X_mcRun2_asymptotic_v3'
if options._beenSet['globalTag']:
    globaltag = options.globalTag

extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['xNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['xFullEvt', extension[options.isMC], options.tag])+'.root')
if not options.inputFiles:
    options.inputFiles = ['/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root'] if not options.isMC else \
                         ['/store/mc/RunIISummer16MiniAODv3/BdToX3872Ks_X3872ToJpsiRho_BMuonFilter_DGamma0_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/70000/34259CC8-8346-E911-ACC1-801844DEF380.root']
annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)

#from Configuration.StandardSequences.Eras import eras
#process = cms.Process('BParkNANO',eras.Run2_2018)
process = cms.Process('BParkNANO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('PhysicsTools.XNano.nanoX_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
)

process.nanoMetadata.strings.tag = annotation
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = cms.untracked.vstring(
      'drop *',
      "keep nanoaodFlatTable_*Table_*_*",     # event data
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')

from PhysicsTools.XNano.nanoX_cff import *
process = nanoAOD_customizeMuonTriggerX(process)
process = nanoAOD_customizeTrackFilteredX(process)
process = nanoAOD_customizeB0ToK0X(process)
process = nanoAOD_customizeTriggerBitsX(process)

# Path and EndPath definitions
process.nanoAOD_B0ToK0X_step = cms.Path(process.nanoSequence + process.nanoB0ToK0XSequence )

# customisation of the process.
if options.isMC:
    from PhysicsTools.XNano.nanoX_cff import nanoAOD_customizeMC
    nanoAOD_customizeMC(process)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(
                                process.nanoAOD_B0ToK0X_step,
                                process.endjob_step,
                                process.NANOAODoutput_step
                               )

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
                                   'nanoAOD_B0ToK0X_step'
                                   )
)


process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
