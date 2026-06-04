#!/bin/bash
# ============================================================
# STEP 2: SLURM Batch Job Script — DeepLabCut Training & Analysis
# Submit this to the scheduler with: sbatch 02_run_dlc_job.sh
# ============================================================

# --- SLURM Job Settings ---
#SBATCH --job-name=DLC_analysis          # Name shown in queue (squeue)
#SBATCH --partition=gpu                  # Use 'gpu_devel' for short tests (<30min), 'gpu' for full runs
#SBATCH --gpus=1                         # Request 1 GPU (required for DLC)
#SBATCH --cpus-per-gpu=4                 # CPU cores per GPU (for data loading)
#SBATCH --mem=32G                        # Total RAM for the job
#SBATCH --time=12:00:00                  # Max runtime: HH:MM:SS (increase for long training)
#SBATCH --output=dlc_job_%j.out          # Standard output log (%j = job ID)
#SBATCH --error=dlc_job_%j.err           # Error log (check this first if something fails)
#SBATCH --mail-type=END,FAIL             # Email you when job ends or fails
#SBATCH --mail-user=YourNetID@yale.edu   # Your Yale email

# --- Load the SAME modules used during environment setup ---
# These MUST match exactly — mismatches cause GPU/CUDA errors.
module purge
module load GCCcore/12.2.0 cuDNN/8.8.0.121-CUDA-12.0.0 miniconda

# --- Activate your DLC conda environment ---
conda activate DEEPLABCUT

# --- Navigate to your working directory ---
cd /nfs/roberts/project/pi_jac52/fc555/DeepLabCut/

# --- Print GPU info to log (useful for debugging) ---
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
nvidia-smi

# --- Run your DLC Python script ---
python 03_dlc_analysis.py

echo "Job finished at: $(date)"
