#!/bin/bash
#SBATCH --job-name=S3_MI
#SBATCH --account=b1028              # Replace with your allocation
#SBATCH --partition=b1028
#SBATCH --time=02:30:00              # 30 minutes should be sufficient
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16           # 16 threads for parallelization
#SBATCH --mem=32G                    # 16GB RAM for processing large files
#SBATCH --output=dsso_%j.out
#SBATCH --error=dsso_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yi-zhi.wang@northwestern.edu  # Replace with your email

# Load required modules
module purge
module load python/3.12.10  # or latest available

# Create virtual environment if not exists
if [ ! -d "$HOME/dsso_env" ]; then
    python -m venv $HOME/dsso_env
    source $HOME/dsso_env/bin/activate
    pip install --upgrade pip
    pip install pandas numpy
else
    source $HOME/dsso_env/bin/activate
fi

# Set working directory
cd $SLURM_SUBMIT_DIR

# Define input/output paths
PEP_LIST="pep_list.csv"
SPECTRA_DIR="/home/ywd617/XL_anal_S3_MI"  # Update this path
OUTPUT_CSV="S3_MI_DSSO_links_FDR_${SLURM_JOB_ID}.csv"
SCRIPT="DSSO_MT_XL.py"  # Or use the optimized version

# Log job information
echo "Job started on $(hostname) at $(date)"
echo "Using ${SLURM_CPUS_PER_TASK} CPUs"
echo "Working directory: $(pwd)"

# Run the DSSO analysis
python ${SCRIPT} ${PEP_LIST} ${SPECTRA_DIR} ${OUTPUT_CSV} --threads ${SLURM_CPUS_PER_TASK}

# Check exit status
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully at $(date)"
    echo "Output saved to: ${OUTPUT_CSV}"
else
    echo "Analysis failed with exit code $? at $(date)"
fi

# Deactivate virtual environment
deactivate
