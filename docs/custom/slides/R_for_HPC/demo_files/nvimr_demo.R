##########################
## Nvim-R Demo (R code) ##
##########################
## Author: Thomas Girke
## Last update: 07-Apr-2022

## Optional: run in interactive session on node with:
# srun --x11 --partition=short --mem=2gb --cpus-per-task 4 --ntasks 1 --time 1:00:00 --pty bash -l

## Start Nvim-connected R session with \rf and then send code by pressing space bar

df <- iris # Assign built-in iris dataset to df
dim(df) # Print dimensions, here number of rows and columns
df[1:4,] # Print first 4 rows
R.version.string
sessionInfo()
