#!/bin/bash
#SBATCH --job-name=dsso_xl_analysis
#SBATCH --account=b1028
#SBATCH --partition=b1028
#SBATCH --time=01:00:00              # Reduced based on actual runtime
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8            # 8 cores is usually optimal for this workload
#SBATCH --mem=4G                     # Reduced based on actual usage
#SBATCH --tmp=20G                    # Request local scratch space
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
SPECTRA_DIR="/home/ywd617/XL_anal_S3_MI"
OUTPUT_CSV="S3_MI_DSSO_links_FDR_${SLURM_JOB_ID}.csv"
SCRIPT="DSSO_MT_XL_optimized.py"  # Use the optimized version

# Log job information
echo "Job started on $(hostname) at $(date)"
echo "Using ${SLURM_CPUS_PER_TASK} CPUs"
echo "Working directory: $(pwd)"
echo "Temporary directory: $TMPDIR"

# Option 1: Use local scratch for better I/O performance
# Uncomment these lines if I/O is still a bottleneck
# echo "Copying spectra files to local scratch..."
# cp -r ${SPECTRA_DIR} $TMPDIR/
# SPECTRA_LOCAL=$TMPDIR/$(basename ${SPECTRA_DIR})
# echo "Files copied to $SPECTRA_LOCAL"

# Option 2: Work directly from network storage (default)
SPECTRA_LOCAL=${SPECTRA_DIR}

# Run timing analysis
echo "Starting analysis..."
time python ${SCRIPT} ${PEP_LIST} ${SPECTRA_LOCAL} ${OUTPUT_CSV} --threads ${SLURM_CPUS_PER_TASK}

# Check exit status
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully at $(date)"
    echo "Output saved to: ${OUTPUT_CSV}"
    
    # If using local scratch, copy results back
    # if [ "$SPECTRA_LOCAL" != "$SPECTRA_DIR" ]; then
    #     cp ${OUTPUT_CSV} $SLURM_SUBMIT_DIR/
    # fi
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
