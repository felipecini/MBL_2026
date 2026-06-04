"""
STEP 3: DeepLabCut Analysis Script
This is the Python script that SLURM will run on the GPU node.

What this script does:
  1. Creates a training dataset (required after moving project to HPC)
  2. Trains the neural network on your labeled frames
  3. Evaluates how well the model trained
  4. Analyzes your behavior videos (detects body parts)
  5. Creates labeled output videos so you can visually verify results

Usage:
  This is called automatically by the SLURM script (02_run_dlc_job.sh).
  You can also run it interactively after requesting GPU resources:
    salloc --partition=gpu_devel --gpus=1 --cpus-per-gpu=2 --mem=16G --time=30:00
    conda activate DEEPLABCUT
    python 03_dlc_analysis.py
"""

import deeplabcut
import os
import time

# ============================================================
# CONFIGURE THESE PATHS FOR YOUR PROJECT
# ============================================================

# Path to your DLC config.yaml file
# This file defines your project: body parts, video paths, network settings
CONFIG_FILE = '/nfs/roberts/scratch/pi_jac52/fc555/DLC analysis/SleepTrack-FC-2026-04-20/config.yaml'

# Path(s) to the video(s) you want to analyze
# You can provide a single video file, or a list of videos:
#   VIDEO_FILES = ['/path/to/video1.avi', '/path/to/video2.avi']
# Or point to a whole folder:
#   VIDEO_FILES = '/nfs/roberts/scratch/pi_mjh24/YourNetID/videos/'
VIDEO_FILES = [
    '/nfs/roberts/scratch/pi_jac52/fc555/DLC analysis/SleepTrack-FC-2026-04-20/videos/FCephys03A_noSD_041726_2026-04-17-104643-0000_2026-04-20-142042-0000.mp4',
]

# Video file extension (used when scanning a folder)
VIDEO_TYPE = '.mp4'

# ============================================================
# STEP 1: Create training dataset
# ============================================================
# IMPORTANT: Always run this after moving your DLC project to HPC.
# It re-generates the training data using the correct Linux file paths.
# You can comment this out on subsequent runs if you haven't re-labeled.

print("\n--- Creating training dataset ---")
deeplabcut.create_training_dataset(CONFIG_FILE)

# ============================================================
# STEP 2: Train the network
# ============================================================
# This is the longest step — can take hours depending on:
#   - Number of labeled frames
#   - Number of training iterations (set in config.yaml)
#   - GPU speed
#
# Key parameters:
#   shuffle=1         -> which data shuffle to train (default: 1)
#   displayiters=100  -> print loss every N iterations
#   saveiters=1000    -> save model checkpoint every N iterations
#   maxiters=None     -> train for the full number set in config.yaml

print("\n--- Training the network ---")
deeplabcut.train_network(
    CONFIG_FILE,
    shuffle=1,
    displayiters=100,
    saveiters=1000,
    maxiters=None,   # Set a number like 50000 to cap training
    allow_growth=True  # Prevents GPU memory errors on shared nodes
)

# ============================================================
# STEP 3: Evaluate the trained network
# ============================================================
# Computes train/test error in pixels. Check the output:
#   - Training error < 5px is excellent
#   - Test error tells you how well it generalizes to new frames
# Results are saved to the evaluation-results folder in your project.

print("\n--- Evaluating the network ---")
deeplabcut.evaluate_network(CONFIG_FILE, plotting=True)

# ============================================================
# STEP 4: Analyze videos
# ============================================================
# Runs pose estimation on each video.
# Output: .h5 and .csv files with (x, y, likelihood) per body part per frame.
#
# Key parameters:
#   save_as_csv=True     -> also save results as a readable CSV
#   destfolder=None      -> save results next to the video (or set a custom path)

print("\n--- Analyzing videos ---")
deeplabcut.analyze_videos(
    CONFIG_FILE,
    VIDEO_FILES,
    videotype=VIDEO_TYPE,
    save_as_csv=True,
    destfolder=None  # Set to a custom path like '/nfs/roberts/scratch/.../results/' if needed
)

# Brief pause to ensure output files are fully written before creating videos
time.sleep(3)

# ============================================================
# STEP 5: Create labeled videos
# ============================================================
# Draws the detected body part markers on top of the original video.
# Use these to visually verify your model is tracking correctly.
#
# Key parameters:
#   draw_skeleton=True   -> connects body parts with lines (defined in config.yaml)
#   filtered=False       -> set True to use filtered predictions (smoother)

print("\n--- Creating labeled videos ---")
deeplabcut.create_labeled_video(
    CONFIG_FILE,
    VIDEO_FILES,
    videotype=VIDEO_TYPE,
    draw_skeleton=True,
    filtered=False
)

print("\n=== All done! ===")
print("Check your project folder for:")
print("  - .h5 / .csv files: raw tracking data")
print("  - _labeled.mp4 files: videos with markers drawn on")
print("  - evaluation-results/: training accuracy metrics")
