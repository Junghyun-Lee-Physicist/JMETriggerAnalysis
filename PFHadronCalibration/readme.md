This directory contains tools and instructions
to derive the PF-Hadron calibrations (PFHCs)
used in the High-Level Trigger (HLT).

#### Introduction

 * General description of the calibration methods,
   incl. references to the relevant documentation
   (publications, analysis notes, presentations, webpages)


#### Production of NTuples for PFHCs

The first step to derive PFHCs is to create an NTuple with the relevant information for the calibration procedure.

## CMSSW Environment Setup
For Winter25 PFHC ntuples:
```bash
ssh <user id>@lxplus8.cern.ch

cmsrel CMSSW_14_2_1
cd CMSSW_14_2_1/src
cmsenv
git cms-init
## --- You can find bellow useful additions to standard CMSSW for relevant studies ---
# Needed: Merge updates from tracking for 2025 - CA automation for patatrack params + mkFit for track building
git cms-merge-topic elusian:1501_newCAtuning 

git clone git@github.com:theochatzis/JMETriggerAnalysis.git

# Build
scram b -j 12
```

## Ntuple Production Using CRAB
> [!NOTE]
> - Before submit CRAB job, You need to check whether the current CMSSW version, the samples global tag, and the HLT menu all match!

```bash
cd ${CMSSW_BASE}/src/JMETriggerAnalysis/PFHadronCalibration/test

# Set Crab Environment Variables
source crabForNtupleGen/prepareCrabEnv.sh

# Check the HLT menu version
#   --> Compare the output to the CMSSW version and the Global Tag (GT) used for the samples
./check_global_tag.sh

# Job Submit
python3 crabForNtupleGen/multiCrab_PFHC_forHLT.py
```

PFHCs are derived using MC events simulating
the production of a single-pion without pileup (NoPU),
with the generated pion $p_{T}$ usually restricted to values below 500 GeV.
Update : The Winter25 sample extends up to 5000 GeV, while PFHC goes up to 1000 GeV.

If you want to find files from the dataset for a local test run, use a script like the one below:
```bash
${CMSSW_BASE}/src/JMETriggerAnalysis/PFHadronCalibration/test/get_test_file_from_das.sh
```

!! Add information on what to do if a sample is not available at a Tier-2 (transfer request, Rucio rule)

The first step in the calibration procedure is to create a flat ROOT NTuple
containing information on selected HLT PF-Candidates,
and the `PFSimParticle` candidates (truth-level information) associated to them.

!! Add more details on how HLT PF-Candidates are selected, and matched to `PFSimParticle` candidates.

To test the first step of the workflow, a small test NTuple can be produced locally as follows:
```
cd ${CMSSW_BASE}/src/JMETriggerAnalysis/PFHadronCalibration/test/
cmsRun pfHadCalibNTuple_cfg.py maxEvents=1000
```

which contains information

!! Add information on how to run the NTuple step with crab, HTCondor and/or other batch systems

#### Derivation of PFHCs

 * Given a set of input PFHC NTuples, how to analyse them to derive the PFHC functions
 * How to produce final `.db` file containing the PFHCs

#### Validation of PFHCs

 * Given a set of PFHCs (i.e. a `.db` file containing the relevant records),
   how to validate this output, and verify that the corrections work as expected

#### Delivery of PFHCs

 * How to upload the relevant records to the appropriate database
 * How to integrate new PFHCs in a Global-Tag (GT)

#### Available sets of PFHCs

List of the available PFHCs for HLT,
with the necessary metadata to be able to reproduce them if needed.
This should include:
 * `git` tag corresponding to the version of `JMETriggerAnalysis` used for the calibrations
 * name of the CMSSW release used
 * names of the datasets in DAS used for the calibrations
 * link to the `cmsRun` configuration file used to create the PFHC-NTuples,
 * links to the PFHC-NTuples used for the calibrations (if still available)
 * link to the relevant analysis scripts to derive and validate the calibrations
 * any further instructions specific to a certain PFHC version (if necessary)

