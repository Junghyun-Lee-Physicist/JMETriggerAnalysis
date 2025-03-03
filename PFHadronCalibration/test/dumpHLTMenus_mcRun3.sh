hltGetConfiguration /dev/CMSSW_13_0_0/GRun/V16 \
   --globaltag 126X_mcRun3_2023_forPU65_v1 \
   --mc \
   --unprescale \
   --output full \
   --max-events 100 \
   --eras Run3 \
   --input /store/mc/Run3Winter23Digi/SinglePionGun_E200to500/GEN-SIM-RAW/EpsilonPU_126X_mcRun3_2023_forPU65_v1-v1/40000/53422f51-f90b-4ec6-ba57-35faab88d66a.root \
   --path MC_AK*,MC_AK*,MC_PFM*,HLTriggerF*,Status* \
   --customise HLTrigger/Configuration/customizeHLTFor2023.customizeHCALFor2023 \
   > hltMC.py

echo """
process.hltOutputFull.outputCommands = [
    'keep *',
    'drop *_hltSiPixelDigisLegacy_*_*',
    'drop *_hltSiPixelClustersLegacy_*_*',
    'drop *_hltSiPixelRecHitsFromLegacy_*_*',
    'drop *_hltEcalDigisLegacy_*_*',
    'drop *_hltEcalUncalibRecHitLegacy_*_*',
    'drop *_hltHbherecoLegacy_*_*',
    ]
""" >> hltMC.py

edmConfigDump hltMC.py > HLT_dev_CMSSW_13_0_0_GRun_configDump.py
