#!/bin/bash
#SBATCH --job-name=dsso_xl_gpu
#SBATCH --account=b1028              # allocation
#SBATCH --partition=b1028              # or gpu-a100
#SBATCH --gres=gpu:1                 # 1 V100 or A100
#SBATCH --time=01:00:00              # plenty for one lane
#SBATCH --cpus-per-gpu=16            # feeds I/O & pandas fallback
#SBATCH --mem=32G
#SBATCH --output=dsso_%j.out
#SBATCH --error=dsso_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yi-zhi.wang@northwestern.edu

module purge
module load python-miniconda3/4.12.0
module load cuda/12.2.2-gcc-10.4.0   # matches RAPIDS 24.06 build

# Activate RAPIDS env created beforehand
source activate rapids-24.06

# -------- paths ----------------------------------------------------------
PEP_LIST="pep_list.csv"
SPECTRA_DIR="/home/ywd617/XL_anal_GPU"         # adjust if needed
OUTPUT_CSV="DSSO_links_FDR_${SLURM_JOB_ID}.csv"
SCRIPT="dsso_link_fdr_gpu.py"              # GPU-aware script

echo "Node  : $(hostname)"
echo "GPU   : $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "Start : $(date)"

python ${SCRIPT} ${PEP_LIST} ${SPECTRA_DIR} ${OUTPUT_CSV} \
       --threads ${SLURM_CPUS_PER_GPU}

status=$?
echo "End   : $(date)  (exit $status)"
exit $status
