#!/bin/bash
#PBS -N job
#PBS -o logfile.log
#PBS -e error_file.err
#PBS -l walltime=36:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb

cd $PBS_O_WORKDIR

source /home/apps/anaconda3/bin/activate
conda activate myenv

which python
python --version

python SPV_Aug_4.py
