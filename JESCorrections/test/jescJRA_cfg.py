import fnmatch

###
### command-line arguments
###
import FWCore.ParameterSet.VarParsing as vpo
opts = vpo.VarParsing('analysis')

opts.register('skipEvents', 0,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'number of events to be skipped')

opts.register('dumpPython', None,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.string,
              'path to python file with content of cms.Process')

opts.register('numThreads', 1,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'number of threads')

opts.register('numStreams', 0,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'number of streams')

opts.register('lumis', None,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.string,
              'path to .json with list of luminosity sections')

opts.register('wantSummary', False,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.bool,
              'show cmsRun summary at job completion')

#opts.register('globalTag', None,
#              vpo.VarParsing.multiplicity.singleton,
#              vpo.VarParsing.varType.string,
#              'argument of process.GlobalTag.globaltag')

opts.register('output', 'out.root',
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.string,
              'path to output ROOT file')

opts.register('bpixMode', None,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.string,
              'if stated as noBPix removes -1.2<phi<-0.8 , else if BPix keeps only this phi region')

opts.register('verbosity', 0,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'level of output verbosity')

opts.parseArguments()

###
### HLT configuration
###

from JMETriggerAnalysis.Common.configs.HLT_dev_CMSSW_13_3_0_GRun_configDump import *

# remove cms.OutputModule objects from HLT config-dump
for _modname in process.outputModules_():
  _mod = getattr(process, _modname)
  if type(_mod) == cms.OutputModule:
    process.__delattr__(_modname)
    if opts.verbosity > 0:
      print('> removed cms.OutputModule:', _modname)

# remove cms.EndPath objects from HLT config-dump
for _modname in process.endpaths_():
  _mod = getattr(process, _modname)
  if type(_mod) == cms.EndPath:
    process.__delattr__(_modname)
    if opts.verbosity > 0:
      print('> removed cms.EndPath:', _modname)

# remove selected cms.Path objects from HLT config-dump
print('-'*108)
print('{:<99} | {:<4} |'.format('cms.Path', 'keep'))
print('-'*108)

# list of patterns to determine which paths to keep
keepPaths = [
  'MC_*Jets*',
  'MC_*AK8Calo*',
]

for _modname in sorted(process.paths_()):
  _keepPath = False
  for _tmpPatt in keepPaths:
    _keepPath = fnmatch.fnmatch(_modname, _tmpPatt)
    if _keepPath: break
  if _keepPath:
    print('{:<99} | {:<4} |'.format(_modname, '+'))
    continue
  _mod = getattr(process, _modname)
  if type(_mod) == cms.Path:
    process.__delattr__(_modname)
    print('{:<99} | {:<4} |'.format(_modname, ''))
print('-'*108)

# remove FastTimerService
if hasattr(process, 'FastTimerService'):
  del process.FastTimerService

# remove MessageLogger
#if hasattr(process, 'MessageLogger'):
#  del process.MessageLogger

#from FWCore.MessageService.MessageLogger_cfi import *
#MessageLogger = cms.Service("MessageLogger")


###
### customisations
###

## customised JME collections
#from JMETriggerAnalysis.Common.customise_hlt import *
#process = addPaths_MC_JMEPFCluster(process)
#process = addPaths_MC_JMEPFCHS(process)
#process = addPaths_MC_JMEPFPuppi(process)

## ES modules for PF-Hadron Calibrations
import os

from CondCore.CondDB.CondDB_cfi import CondDB as _CondDB
## ES modules for PF-Hadron Calibrations
process.pfhcESSource = cms.ESSource('PoolDBESSource',
   _CondDB.clone(connect = 'sqlite_file:/afs/cern.ch/user/t/tchatzis/public/PFHC2024/PFCalibration.db'),
   toGet = cms.VPSet(
     cms.PSet(
       record = cms.string('PFCalibrationRcd'),
       tag = cms.string('PFCalibration_HLT_133X_mcRun3_2024_realistic_v9'),
       label = cms.untracked.string('HLT'),
     ),
  ),
)

process.pfhcESPrefer = cms.ESPrefer('PoolDBESSource', 'pfhcESSource')
## Used to test applying the offline PFHC : 
#process.hltParticleFlow.calibrationsLabel = '' # standard label for Offline-PFHC in GT

# Test option to skip forward PFHC application (after eta = 2.5)
#process.hltParticleFlow.skipForwardCalibrations = cms.bool(True)


## Modification to select jets 
process.hltAK4PFJetsBPix = cms.EDFilter( "PFJetSelector",
    src = cms.InputTag("hltAK4PFJets"),
    filter = cms.bool(False),
    cut = cms.string("phi<-0.80 && phi>-1.20")
)

process.hltAK4PFJetsNoBPix = cms.EDFilter( "PFJetSelector",
    src = cms.InputTag("hltAK4PFJets"),
    filter = cms.bool(False),
    cut = cms.string("!(phi<-0.80 && phi>-1.20)")
)


process.HLTAK4PFJetsReconstructionSequence += process.hltAK4PFJetsBPix
process.HLTAK4PFJetsReconstructionSequence += process.hltAK4PFJetsNoBPix

process.hltAK8PFJetsBPix = cms.EDFilter( "PFJetSelector",
    src = cms.InputTag("hltAK8PFJets"),
    filter = cms.bool(False),
    cut = cms.string("phi<-0.80 && phi>-1.20")
)

process.hltAK8PFJetsNoBPix = cms.EDFilter( "PFJetSelector",
    src = cms.InputTag("hltAK8PFJets"),
    filter = cms.bool(False),
    cut = cms.string("!(phi<-0.80 && phi>-1.20)")
)


process.HLTAK8PFJetsReconstructionSequence += process.hltAK8PFJetsBPix
process.HLTAK8PFJetsReconstructionSequence += process.hltAK8PFJetsNoBPix

ak4hltCollectionName_ = 'hltAK4PFJets'
ak8hltCollectionName_ = 'hltAK8PFJets'
if opts.bpixMode == 'noBPix':
   ak4hltCollectionName_ = 'hltAK4PFJetsNoBPix'
   ak8hltCollectionName_ = 'hltAK8PFJetsNoBPix'

if opts.bpixMode == 'BPix':
   ak4hltCollectionName_ = 'hltAK4PFJetsBPix'
   ak8hltCollectionName_ = 'hltAK8PFJetsBPix'

###
### Jet Response Analyzer (JRA) NTuple
###
from JMETriggerAnalysis.JESCorrections.jescJRA_utils import addJRAPath

addJRAPath(process, genJets = 'ak4GenJetsNoNu', maxDeltaR = 0.2, moduleNamePrefix = 'ak4caloHLT'     , recoJets = 'hltAK4CaloJets'      , rho = 'hltFixedGridRhoFastjetAllCalo')
addJRAPath(process, genJets = 'ak4GenJetsNoNu', maxDeltaR = 0.2, moduleNamePrefix = 'ak4pfHLT'       , recoJets = ak4hltCollectionName_ , rho = 'hltFixedGridRhoFastjetAll')
#addJRAPath(process, genJets = 'ak4GenJetsNoNu', maxDeltaR = 0.2, moduleNamePrefix = 'ak4pfchsHLT'    , recoJets = 'hltAK4PFCHSJets'    , rho = 'hltFixedGridRhoFastjetAll')
#addJRAPath(process, genJets = 'ak4GenJetsNoNu', maxDeltaR = 0.2, moduleNamePrefix = 'ak4pfpuppiHLT'  , recoJets = 'hltAK4PFPuppiJets'  , rho = 'hltFixedGridRhoFastjetAll')

addJRAPath(process, genJets = 'ak8GenJetsNoNu', maxDeltaR = 0.4, moduleNamePrefix = 'ak8caloHLT'     , recoJets = 'hltAK8CaloJets'      , rho = 'hltFixedGridRhoFastjetAllCalo')
addJRAPath(process, genJets = 'ak8GenJetsNoNu', maxDeltaR = 0.4, moduleNamePrefix = 'ak8pfHLT'       , recoJets = ak8hltCollectionName_ , rho = 'hltFixedGridRhoFastjetAll')
#addJRAPath(process, genJets = 'ak8GenJetsNoNu', maxDeltaR = 0.4, moduleNamePrefix = 'ak8pfchsHLT'    , recoJets = 'hltAK8PFCHSJets'    , rho = 'hltFixedGridRhoFastjetAll')
#addJRAPath(process, genJets = 'ak8GenJetsNoNu', maxDeltaR = 0.4, moduleNamePrefix = 'ak8pfpuppiHLT'  , recoJets = 'hltAK8PFPuppiJets'  , rho = 'hltFixedGridRhoFastjetAll')

###
### standard options
###

# max number of events to be processed
process.maxEvents.input = opts.maxEvents

# number of events to be skipped
process.source.skipEvents = cms.untracked.uint32(opts.skipEvents)

# multi-threading settings
process.options.numberOfThreads = max(opts.numThreads, 1)
process.options.numberOfStreams = max(opts.numStreams, 0)

# show cmsRun summary at job completion
process.options.wantSummary = cms.untracked.bool(opts.wantSummary)

## update process.GlobalTag.globaltag
#if opts.globalTag is not None:
#   from Configuration.AlCa.GlobalTag import GlobalTag
#   process.GlobalTag = GlobalTag(process.GlobalTag, opts.globalTag, '')

# select luminosity sections from .json file
if opts.lumis is not None:
   import FWCore.PythonUtilities.LumiList as LumiList
   process.source.lumisToProcess = LumiList.LumiList(filename = opts.lumis).getVLuminosityBlockRange()

# create TFileService to be accessed by JRA-NTuple plugin
process.TFileService = cms.Service('TFileService', fileName = cms.string(opts.output))

# input EDM files [primary]
if opts.inputFiles:
   process.source.fileNames = opts.inputFiles
else:
   process.source.fileNames = [
     '/store/mc/Run3Winter24Digi/QCD_PT-15to7000_TuneCP5_13p6TeV_pythia8/GEN-SIM-RAW/FlatPU0to120_133X_mcRun3_2024_realistic_v9-v3/40000/0d532668-aed7-403d-943e-86f88480b3b3.root',
   ]

# input EDM files [secondary]
if not hasattr(process.source, 'secondaryFileNames'):
   process.source.secondaryFileNames = cms.untracked.vstring()
if opts.secondaryInputFiles:
   process.source.secondaryFileNames = opts.secondaryInputFiles

# dump content of cms.Process to python file
if opts.dumpPython is not None:
   open(opts.dumpPython, 'w').write(process.dumpPython())

# printouts
if opts.verbosity > 0:
   print('--- jescJRA_cfg.py ---')
   print('')
   print('option: output =', opts.output)
   print('option: reco =', opts.reco)
   print('option: dumpPython =', opts.dumpPython)
   print('')
   print('process.GlobalTag =', process.GlobalTag.dumpPython())
   print('process.source =', process.source.dumpPython())
   print('process.maxEvents =', process.maxEvents.dumpPython())
   print('process.options =', process.options.dumpPython())
   print('----------------------')
