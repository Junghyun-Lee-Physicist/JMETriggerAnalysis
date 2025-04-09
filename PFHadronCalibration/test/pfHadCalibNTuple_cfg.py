
###
### command-line arguments
###
import FWCore.ParameterSet.VarParsing as vpo
opts = vpo.VarParsing('analysis')

opts.register('skipEvents', 0,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'number of events to be skipped')

opts.register('numThreads', 1,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'number of threads')

opts.register('numStreams', 0,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'number of streams')

opts.register('dumpPython', None,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.string,
              'path to python file with content of cms.Process')

opts.register('output', 'pfHC_Online_NTuple.root',
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.string,
              'path to output ROOT file')

opts.register('verbosity', 0,
              vpo.VarParsing.multiplicity.singleton,
              vpo.VarParsing.varType.int,
              'verbosity level')

opts.parseArguments()


###
### HLT configuration
###
from JMETriggerAnalysis.Common.configs.HLT_dev_CMSSW_14_2_0_GRun_configDump import *
##from JMETriggerAnalysis.Common.configs.HLT_dev_CMSSW_14_0_0_GRun_configDump import *


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
for _modname in sorted(process.paths_()):
  _keepPath = _modname.startswith('MC_') and ('Jets' in _modname or 'MET' in _modname)
  _keepPath |= _modname.startswith('MC_ReducedIterativeTracking')
  if _keepPath:
    print('{:<99} | {:<4} |'.format(_modname, '+'))
    continue
  _mod = getattr(process, _modname)
  if type(_mod) == cms.Path:
    process.__delattr__(_modname)
    print('{:<99} | {:<4} |'.format(_modname, ''))
print('-'*108)

# remove FastTimerService
del process.FastTimerService


###
### PFHC-specific paths
###
process.particleFlowSimParticle = cms.EDProducer('PFSimParticleProducer',
  Fitter = cms.string('KFFittingSmoother'),
  MCTruthMatchingInfo = cms.untracked.bool(False),
  ParticleFilter = cms.PSet(
    EMin = cms.double(0),
    chargedPtMin = cms.double(0),
    etaMax = cms.double(5.3),
    invisibleParticles = cms.vint32(),
    protonEMin = cms.double(5000.0),
    rMax = cms.double(129.0),
    zMax = cms.double(317.0),
  ),
  Propagator = cms.string('PropagatorWithMaterial'),
  RecTracks = cms.InputTag('trackerDrivenElectronSeeds'),
  TTRHBuilder = cms.string('WithTrackAngle'),
  ecalRecHitsEB = cms.InputTag('caloRecHits','EcalRecHitsEB'),
  ecalRecHitsEE = cms.InputTag('caloRecHits','EcalRecHitsEE'),
  fastSimProducer = cms.untracked.InputTag('fastSimProducer','EcalHitsEB'),
  process_Particles = cms.untracked.bool(True),
  process_RecTracks = cms.untracked.bool(False),
  sim = cms.InputTag('g4SimHits'),
  verbose = cms.untracked.bool(False)
)

process.pfHadCalibNTuple = cms.EDAnalyzer('PFHadCalibNTuple',
  TTreeName  = cms.string('s'),
  TTreeTitle = cms.string('PFCalibration'),

  printDetailLog = cms.bool(False),

  genParticles = cms.InputTag('genParticles'),
  pfSimParticles = cms.InputTag('particleFlowSimParticle'),
  recoPFCandidates = cms.InputTag('hltParticleFlow'),

  genParticleStatus = cms.int32(1), # status code of selected GEN particles
  genParticlePdgId = cms.int32(-211), # pdgID of selected GEN particles (-211 = charged hadron)
  genParticleIsoMinDeltaR = cms.double(1.0), # min deltaR between "isolated" GEN particles and other GEN particles

  minPt = cms.double(1.00), # min pt
  minTrackP = cms.double(1.00), # min track momentum
  minTrackPt = cms.double(1.00), # min track transverse momentum
  ##minCaloEnergy = cms.double(0.5),
  minCaloEnergy = cms.double(1.00), # min ecal+hcal energy
  maxECalEnergy = cms.double(1e12), # max ecal energy

  # min nb of pixel and pixel+strip hits (per track-eta range)
  # --> which corresponds to hitPattern.numberOfValidPixelHits()
  ####minPixelHits = cms.vuint32(2, 2, 2, 2, 2),
  minPixelHits = cms.vuint32(4, 4, 4, 4, 4),
  # In run3, HLT has only 3,4 pixel hits tracks
  # https://mattermost.web.cern.ch/cms-jme-pog/pl/8oteynsbk7f39kn3zdrniidyja
  # https://mattermost.web.cern.ch/cms-jme-pog/pl/i6hogtich7bricpp5g9i954bir
  # "worst" case of 3 pixel hit:
  #     --> either 3 pixel + at least one strip hit. Usually this one strip hit is enough to improve the reconstruction from fakes Pixel tracks with 3 hits
  # 4 pixel hits case:
  #     --> these are already pretty pure.

  #minTrackerHits = cms.vuint32(14, 17, 20, 17, 10), # Offline Criterion
  minTrackerHits = cms.vuint32(2, 2, 2, 2, 2), # HLT Criterion
  # -->which corresponds to hitPattern.numberOfValidTrackerHits()
  maxEtaForMinTrkHitsCuts = cms.vdouble(1.4, 1.6, 2.0, 2.4, 2.6),

  usePFBlockElements = cms.bool(True),

  rootOutputFile = cms.string(opts.output),
)

process.pfSimParticleSeq = cms.Sequence(process.particleFlowSimParticle)
process.pfSimParticlePath = cms.Path(process.pfSimParticleSeq)
process.schedule.append(process.pfSimParticlePath)

process.pfHadCalibNTupleSeq = cms.Sequence(process.pfHadCalibNTuple)
process.pfHadCalibNTupleEndPath = cms.EndPath(process.pfHadCalibNTupleSeq)
process.schedule.append(process.pfHadCalibNTupleEndPath)


###
### options
###

# number of events
process.maxEvents.input = opts.maxEvents

# number of events to be skipped
process.source.skipEvents = cms.untracked.uint32(opts.skipEvents)

# number of threads/streams
process.options.numberOfThreads = opts.numThreads
process.options.numberOfStreams = opts.numStreams

process.options.wantSummary = cms.untracked.bool(False)

# input files
if opts.inputFiles:
  process.source.fileNames = opts.inputFiles
else:
  process.source.fileNames = [
      '/store/mc/Run3Winter25Digi/Pi_Par-E-0p2to200_PGun/GEN-SIM-RAW/NoPU_142X_mcRun3_2025_realistic_v7-v1/2820000/78bda28b-2bc3-47ab-a099-12d0bd806c91.root'
  ]


###
### output file
###
####process.TFileService = cms.Service('TFileService', fileName = cms.string(opts.output))

# dump content of cms.Process to python file
if opts.dumpPython is not None:
  process.prune()
  open(opts.dumpPython, 'w').write(process.dumpPython())
