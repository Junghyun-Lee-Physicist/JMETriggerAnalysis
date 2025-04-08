import os
import sys

############# USER DEFINED ##################################################

jobID = "07Apr2025_forV9"

# Directory with input files
username = os.environ['USER']
first_letter = username[0]

cmssw_base = os.environ.get("CMSSW_BASE")
if not cmssw_base:
    print("Error: CMSSW_BASE environment variable is not set.")
    sys.exit(1)
# Set path for executable and related files.
# They are located in ${CMSSW_BASE}/src/JMETriggerAnalysis/JESCorrections/test/
jesc_work_dir = os.path.join(cmssw_base, "src", "JMETriggerAnalysis", "JESCorrections", "test")

storagePath = f"root://eosuser.cern.ch//eos/user/{first_letter}/{username}"
print(f"  [ Log ] : Setted storage path --> {storagePath}")
input_directory = f"{storagePath}/JESC/JESC_ntuple/250407" # Input ntuple path
output_directory = f"{storagePath}/JESC/JESC_output/250407/{jobID}" # Set JEC output directory path
print(f"  [ Log ] : n-tuple path --> {input_directory}")
print(f"  [ Log ] : output path --> {output_directory}")
o
# Set number of event for JEC
#nOfEvt = 10000000
nOfEvt = 5000000

# Set proxy name
# The proxy file will be generated in the current working directory
x509_proxy_path = os.path.join(jesc_work_dir, "condorForJEC", "proxy.cert")

bpix_categories = ['noBPix','BPix','FPix']
#flatPU_label = 'FlatPU0to80'
flatPU_label = 'FlatPU0to120'
noPU_label = 'NoPU'

doPuppiCHS = False

log_dir_name = f"htc_out_{jobID}"
log_dir_path = os.path.join(jesc_work_dir, "condorForJEC", log_dir_name)
if os.path.isdir(log_dir_path):
    print(f"Error: The directory '{log_dir_path}' already exists.")
    sys.exit(1)
else:
    os.makedirs(log_dir_path)
    print(f"Directory '{log_dir_path}' created successfully.")

log_path = os.path.join(os.getcwd(), log_dir_path)

#############################################################################

# executable file path
executable_path = os.path.join(jesc_work_dir, "fitJESCs")

# Loop over each input file
submitFileName = f"sub_jecs_{jobID}.htc"
submitFilePath = os.path.join(jesc_work_dir, "condorForJEC", submitFileName)
with open(submitFilePath, "w") as job_script:
    # Common parts in all jobs
    job_script.write("------------------------------------------------------------------------------------------------------------------\n")
    job_script.write(f"executable            = {executable_path}\n")
    job_script.write("getenv                = True\n")
    job_script.write("should_transfer_files = YES\n")
    job_script.write("when_to_transfer_output = ON_EXIT_OR_EVICT\n")
    job_script.write("output_destination      =  %s\n" %(output_directory) )
    job_script.write("MY.XRDCP_CREATE_DIR = True\n")
    job_script.write("MY.WantOS = \"el8\"\n")
    job_script.write(f"x509userproxy          = {x509_proxy_path}\n")
    job_script.write("request_cpus            = 4\n")
    job_script.write("request_memory          = 8 GB\n")
    job_script.write("------------------------------------------------------------------------------------------------------------------\n")

    for bpix_cat in bpix_categories:
        job_script.write("\n"*4)
        job_script.write(f"--- {jobID}_JESCs_{bpix_cat} ---\n")
        job_script.write(f"transfer_output_files = {jobID}_JESCs_{bpix_cat}\n")
        jet_categories = ['ak4pfHLT','ak8pfHLT']
        if doPuppiCHS:
            jet_categories += ['ak4pfPuppiHLT','ak4pfCHSHLT']
        if bpix_cat == 'noBPix':
            jet_categories += ['ak4caloHLT','ak8caloHLT']
        for jet_cat in jet_categories:
            job_script.write("\n"*2)
            job_script.write(f"arguments = -i_nopu {input_directory}/JRA_{noPU_label}{bpix_cat}.root -i_flatpu {input_directory}/JRA_{flatPU_label}{bpix_cat}.root -o {jobID}_JESCs_{bpix_cat} -b -j {jet_cat} -n {nOfEvt}\n")
            job_script.write(f"output = {log_path}/{bpix_cat}_{jet_cat}.out\n")
            job_script.write(f"error = {log_path}/{bpix_cat}_{jet_cat}.err\n")
            job_script.write(f"log = {log_path}/{bpix_cat}_{jet_cat}.log\n")
            #job_script.write("+JobFlavour = \"testmatch\"\n") # 3 days
            job_script.write("+JobFlavour = \"tomorrow\"\n") # 1 day
            #job_script.write("+JobFlavour = \"workday\"\n") # 8 hours
            #job_script.write("+JobFlavour = \"microcentury\"\n") # 1 hours
            job_script.write("queue\n")
