#!/bin/bash
#SBATCH --job-name=job
#SBATCH --output=logfile.log
#SBATCH --error=error_file.err
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G

module purge
module load miniconda

source /home/apps/MLDL/DL-CondaPy3/etc/profile.d/conda.sh
conda activate spv_env

cd $SLURM_SUBMIT_DIR

which python
python --version

python SPV_June_18.py
