#!/bin/bash
#SBATCH --job-name=dsso_xl_analysis
#SBATCH --account=b1028
#SBATCH --partition=b1028
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=dsso_%j.out
#SBATCH --error=dsso_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yi-zhi.wang@northwestern.edu

# Load Python
module purge
module load python/3.12.10

# Setup environment
source $HOME/dsso_env/bin/activate || {
    python -m venv $HOME/dsso_env
    source $HOME/dsso_env/bin/activate
    pip install pandas numpy
}

# Run analysis
cd $SLURM_SUBMIT_DIR
echo "Starting at $(date) with ${SLURM_CPUS_PER_TASK} CPUs"

python DSSO_MT_XL_hybrid.py \
    pep_list.csv \
    /home/ywd617/XL_anal_S4_MI \
    S4_MI_DSSO_links_${SLURM_JOB_ID}.csv \
    --threads ${SLURM_CPUS_PER_TASK}

echo "Finished at $(date)"
