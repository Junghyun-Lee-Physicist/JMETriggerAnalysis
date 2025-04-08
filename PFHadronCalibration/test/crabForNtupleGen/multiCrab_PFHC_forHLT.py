import CRABClient
from CRABClient.UserUtilities import config, getUsername

jobID="25Feb2025"

config = config()

config.section_("General")
config.General.transferLogs = False
config.General.workArea = f'CrabNTupleProduction_forOnlinePFHC_{jobID}'

config.section_("JobType")
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'pfHadCalibNTuple_cfg.py' # PFHC configure file
config.JobType.outputFiles = ['pfHC_Online_NTuple.root']
#config.JobType.numCores = 4
config.JobType.allowUndistributedCMSSW = True # [ True ] for the recent version of CMSSW work
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 180
config.Data.totalUnits = -1
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 500
#config.Data.totalUnits = 6429
config.Data.publication = False

# To use the production state samples, activate below line
#config.Data.allowNonValidInputDataset = True 

config.section_("Site")
config.Site.storageSite = 'T3_KR_KNU'
#config.Site.whitelist = ['T2_CH_CERN']
#config.Site.blacklist = ['T1_US_FNAL','T2_UK_London_Brunel']

from CRABAPI.RawCommand import crabCommand
from multiprocessing import Process
import copy
import sys

def submitJob(config):
    try:
        crabCommand('submit', config = config)
    except Exception as e:
        print(f"Error submitting job {config.General.requestName}: {e}")

# Loop over datasets
for pionGun_E in [ '0p2to10', '0p2to200' ]:
    job_config = copy.deepcopy(config) # Create a fresh copy of config

    job_config.General.requestName = f"SinglePion_E_{pionGun_E}_{jobID}"
    job_config.Data.outLFNDirBase = f"/store/user/{getUsername()}/PFHC_HLT_Ntuple/{jobID}/{job_config.General.requestName}"
    #config.Data.useParent = True
    #config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'

    # Dataset selection
    if pionGun_E == '0p2to10':
        job_config.Data.inputDataset = '/Pi_Par-E-0p2to10_PGun/Run3Winter25Digi-NoPU_142X_mcRun3_2025_realistic_v7-v2/GEN-SIM-RAW'
    elif pionGun_E == '0p2to200':
        job_config.Data.inputDataset = '/Pi_Par-E-0p2to200_PGun/Run3Winter25Digi-NoPU_142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RAW'
    elif pionGun_E == '200to500':
        job_config.Data.inputDataset = '/Pi_Par-E-200to500_PGun/Run3Winter25Digi-NoPU_142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RAW'
    elif pionGun_E == '500to5000':
        job_config.Data.inputDataset = '/Pi_Par-E-500to5000_PGun/Run3Winter25Digi-NoPU_142X_mcRun3_2025_realistic_v7-v2/GEN-SIM-RAW'
    else :
        print("ERROR : There is wierd pion gun energy.. please take a look [ pionGun_E array]")
        sys.exit(2)

    p = Process(target=submitJob, args=(job_config,))
    p.start()
    p.join()

