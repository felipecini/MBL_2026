#!/bin/bash
# ============================================================
# STEP 1: DeepLabCut Environment Setup — Yale HPC (Bouchet)
# Run this ONCE interactively from the login node.
# This script installs the conda environment needed for DLC.
# ============================================================

# --- 1. Request an interactive session on a compute node ---
# (Run this manually first, then proceed with steps below)
salloc --cpus-per-task=1 --time=1:00:00 --partition=devel --mem-per-cpu=32G

# --- 2. Load required HPC modules ---
# These provide the GPU drivers and Python environment manager.
# Make sure these versions match what you use when RUNNING DLC later!
module purge
module load GCCcore/12.2.0 cuDNN/8.8.0.121-CUDA-12.0.0 miniconda

# --- 3. (Optional) Remove old DLC environment if it exists ---
# Only needed if you want a clean reinstall.
# conda env remove -n DEEPLABCUT

# --- 4. Navigate to YOUR working directory ---
# Replace YourNetID with your actual Yale NetID
cd /nfs/roberts/project/pi_jac52/fc555/

# --- 5. Create the conda environment from the DLC YAML file ---
# Make sure DEEPLABCUT.yaml is present in this directory.
# You can get it from: github.com/DeepLabCut/DeepLabCut/conda-environments
conda env create -f DEEPLABCUT.yaml

echo "Done! Environment created. You can now submit DLC jobs."
