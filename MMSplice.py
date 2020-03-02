####################################################################################################################
####################################################################################################################
# Main code for calling MMSplice jobs.
# Author: Haiying Kong
# Last Modified: 15 November 2019
####################################################################################################################
####################################################################################################################
import os
import fnmatch
from bsub import bsub

####################################################################################################################
####################################################################################################################
# Get list of vcf files in data folder.
data_dir = '/icgc/dkfzlsdf/analysis/B240/kong/Projects/PANCSTRAT/Data/WGS/MMSplice/K20K/'
vcfs = os.listdir(data_dir)
vcfs = fnmatch.filter(vcfs, '*.vcf')

# Define folder to save results and clear up the folder.
lock_dir = '/icgc/dkfzlsdf/analysis/B240/kong/Projects/PANCSTRAT/Lock/WGS/Kipoi/MMSplice/Raw/K20K/'
file_names = os.listdir(lock_dir)
file_names = fnmatch.filter(file_names, '*.txt')
for file_name in file_names:
    os.remove(lock_dir + file_name)

# Define folder to save error and log files and clear up the folder.
err_dir = '/icgc/dkfzlsdf/analysis/B240/kong/Projects/PANCSTRAT/Err_Out/Kipoi/MMSplice/K20K/'
file_names = os.listdir(err_dir)
file_names = fnmatch.filter(file_names, '*.*')
for file_name in file_names:
    os.remove(err_dir + file_name)

# Submit jobs for each vcf file.
for vcf in vcfs:
    tumor = vcf.replace('.vcf', '')
    job_name = err_dir + tumor
    job = bsub(job_name, W='20:00', M='10G', verbose=True)
    args = tumor + ' ' + data_dir + ' ' + lock_dir
    job("module load anaconda3/2019.07; source activate kipoi-MMSplice; python /icgc/dkfzlsdf/analysis/B240/kong/Projects/PANCSTRAT/Code/WGS/Kipoi/MMSplice/Run_MMSplice/MMSplice_bsub_K20K.py" + ' ' + args)


####################################################################################################################
####################################################################################################################
