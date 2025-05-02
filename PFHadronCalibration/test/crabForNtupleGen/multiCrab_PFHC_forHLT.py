#!/usr/bin/env python3
import os
import argparse
from CRABClient.UserUtilities import config as getConfig, getUsername
from CRABAPI.RawCommand import crabCommand

def makeConfig(energy, jobID, workArea):
    """
    Create and return a fresh CRAB config for the given pion-gun energy.
    """
    cfg = getConfig()

    # General
    cfg.section_("General")
    cfg.General.workArea    = workArea
    cfg.General.requestName = f"SinglePion_E_{energy}_{jobID}"
    cfg.General.transferLogs = False

    # JobType
    cfg.section_("JobType")
    cfg.JobType.pluginName              = 'Analysis'
    cfg.JobType.psetName                = 'pfHadCalibNTuple_cfg.py'
    cfg.JobType.outputFiles             = ['pfHC_Online_NTuple.root']
    cfg.JobType.maxMemoryMB             = 2500
    cfg.JobType.allowUndistributedCMSSW = True

    # Data
    cfg.section_("Data")
    cfg.Data.splitting    = 'FileBased'
    cfg.Data.unitsPerJob  = 2
    cfg.Data.totalUnits   = -1
    cfg.Data.publication  = False
    cfg.Data.ignoreLocality = True 

    # Site
    cfg.section_("Site")
    cfg.Site.storageSite = 'T2_CH_CERN'

    return cfg

def submit(cfg):
    """Submit a fresh CRAB job."""
    try:
        print(f"[Submit   ] {cfg.General.requestName}")
        crabCommand('submit', config=cfg)
    except Exception as e:
        print(f"[Error    ] submit {cfg.General.requestName}: {e}")

def resubmit(taskDir):
    """Resubmit failed jobs in the given CRAB project directory."""
    if not os.path.isdir(taskDir):
        print(f"[Warning  ] taskDir not found: {taskDir}")
        return
    try:
        print(f"[Resubmit ] {taskDir}")
        crabCommand('resubmit', dir=taskDir)
    except Exception as e:
        print(f"[Error    ] resubmit {taskDir}: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Submit or resubmit CRAB jobs for PFHC ntuples"
    )
    parser.add_argument(
        "--resubmit", action="store_true",
        help="If set, resubmit failed jobs instead of fresh submissions"
    )
    args = parser.parse_args()

    jobID    = "21Apr2025"
    workArea = f"NTupleProduction_forOnlinePFHC_{jobID}"
    user     = getUsername()

    # Map from pion-gun energy to its input dataset
    dataset_map = {
        '0p2to10':   '/Pi_Par-E-0p2to10_PGun/Run3Winter25Digi-NoPU_142X_mcRun3_2025_realistic_v7-v2/GEN-SIM-RAW',
        '0p2to200':  '/Pi_Par-E-0p2to200_PGun/Run3Winter25Digi-NoPU_142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RAW',
        '200to500':  '/Pi_Par-E-200to500_PGun/Run3Winter25Digi-NoPU_142X_mcRun3_2025_realistic_v7-v1/GEN-SIM-RAW',
        '500to5000': '/Pi_Par-E-500to5000_PGun/Run3Winter25Digi-NoPU_142X_mcRun3_2025_realistic_v7-v2/GEN-SIM-RAW',
    }

    for energy, dataset in dataset_map.items():
        cfg = makeConfig(energy, jobID, workArea)

        # Must specify input dataset and output LFN base
        cfg.Data.inputDataset   = dataset
        cfg.Data.outLFNDirBase  = (
            f"/store/group/phys_jetmet/{user}/"
            f"JMETriggerAnalysis/PFHC/ntuple/{jobID}/{cfg.General.requestName}"
        )

        taskDir = os.path.join(workArea, f"crab_{cfg.General.requestName}")
        if args.resubmit:
            resubmit(taskDir)
        else:
            submit(cfg)

if __name__ == "__main__":
    main()
