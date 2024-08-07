#!/bin/bash
#SBATCH -J snakemakePathseq
#SBATCH --output=snakemakePathseq.out

# Conda Initialization (only needed once per session)
# Assuming you've run `conda init bash` before
# If not, uncomment the next line and adjust the shell if needed
#conda init bash
source ~/.bashrc
conda activate ~/.conda/envs/pathseq/

# Check if activation was successful
if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to activate conda environment"
    exit 1
fi

# Explicitly add Conda's bin directory to PATH
export PATH="$CONDA_PREFIX/bin:$PATH"

# Print environment variables and PATH after activation
# echo "Environment variables:"
# env | grep -i conda
# echo "PATH:"
# echo $PATH

snakemake -p -v --conda-frontend conda

echo "RAN SNAKEMAKE"

# conda deactivate