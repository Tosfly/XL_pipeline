#!/bin/bash
#SBATCH --job-name=S4_BP
#SBATCH --account=b1028
#SBATCH --partition=b1028
#SBATCH --time=02:30:00              # Reduced based on actual runtime
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8            # 8 cores is usually optimal for this workload
#SBATCH --mem=4G                     # Reduced based on actual usage
#SBATCH --output=dsso_%j.out
#SBATCH --error=dsso_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=yi-zhi.wang@northwestern.edu

# Load required modules
module purge
module load python/3.12.10

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

# Define paths
PEP_LIST="pep_list.csv"
SPECTRA_DIR="/home/ywd617/XL_anal_S4_BP"
OUTPUT_CSV="S4_BP_DSSO_links_FDR_${SLURM_JOB_ID}.csv"
SCRIPT="DSSO_MT_XL_optimized.py"  # Use the optimized version

# Log job information
echo "Job started on $(hostname) at $(date)"
echo "Using ${SLURM_CPUS_PER_TASK} CPUs"
echo "Working directory: $(pwd)"

# Work directly from network storage
SPECTRA_LOCAL=${SPECTRA_DIR}

# Run timing analysis
echo "Starting analysis..."
time python ${SCRIPT} ${PEP_LIST} ${SPECTRA_LOCAL} ${OUTPUT_CSV} --threads ${SLURM_CPUS_PER_TASK}

# Check exit status
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully at $(date)"
    echo "Output saved to: ${OUTPUT_CSV}"
else
    echo "Analysis failed with exit code $? at $(date)"
fi

# Optional: Generate performance report
echo -e "\n=== Performance Summary ==="
echo "Memory usage:"
sstat -j $SLURM_JOB_ID --format=JobID,MaxRSS,AveCPU

# Deactivate virtual environment
deactivate

# Alternative submission with profiling enabled:
# Add this to your Python script call for detailed profiling
# python -m cProfile -o profile_${SLURM_JOB_ID}.prof ${SCRIPT} ${PEP_LIST} ${SPECTRA_LOCAL} ${OUTPUT_CSV} --threads ${SLURM_CPUS_PER_TASK}
# Then analyze with: python -m pstats profile_${SLURM_JOB_ID}.prof
