#!/bin/sh
#
# Simple "Hello World" submit script for Slurm.
#
# Replace ACCOUNT with your account name before submitting.
#
#SBATCH --job-name=ost-rec_protenix    # The job name
#SBATCH --time=0-23:59            # The time the job will take to run in D-HH:MM
#SBATCH --gres=gpu:1
#SBATCH --output logs/%x-%A.out

python3 /lus/lfs1aip2/scratch/s5h/omeir.s5h/projects/openbind/template_based_accuracy_estimation/utilX3_Py_get_receptor_ost_metrics.py -r=results_reformatted/ -f=../../fragalysis_data/A71EV2A_2026-02-04/aligned_files/ -o=ost-receptor_compare

