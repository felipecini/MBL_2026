#!/bin/bash

salloc --cpus-per-task=1 --time=1:00:00 --partition=devel --mem-per-cpu=32G


module purge
module load GCCcore/12.2.0 cuDNN/8.8.0.121-CUDA-12.0.0 miniconda

cd /nfs/roberts/project/pi_jac52/fc555/

conda env create -f DEEPLABCUT.yaml

echo "Done! Environment created. You can now submit DLC jobs."
