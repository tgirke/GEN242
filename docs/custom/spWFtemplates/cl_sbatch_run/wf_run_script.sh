#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=00:15:00 # 15 minutes
#SBATCH --mail-user=useremail@address.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="workflow_finished"
#SBATCH --partition=gen242 # Choose queue/partition from: intel, batch, highmem, gpu, short
#SBATCH --account=gen242 # The 'account' specification may not be required for every user

Rscript wf_run_script.R
