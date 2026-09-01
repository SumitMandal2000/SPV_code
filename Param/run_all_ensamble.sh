#!/bin/bash

l=32
N1=1
N=100
pin_all=(0 2 4 6 8 10)
p0=3.70
v0=0.30
# Set the base directory for your files
WORK_DIR="$(pwd)"

# Path to the original files
SPV_PY="$WORK_DIR/SPV_June_18.py"
JOB_SH="$WORK_DIR/job.sh"

# Loop over each p0 value
for pin in "${pin_all[@]}"; do
    echo "Processing for pin=${pin}"

    # Inner loop for i from 1 to N
    for i in $(seq $N1 $N); do
        #echo "Processing iteration i=$i for p0=$p0"

        # 1. Modify the SPV_31_10.py file, change p0 in line 148, and save with a new name
        cp "$SPV_PY" "$WORK_DIR/SPV_p0_${p0}_v0_${v0}_pin_${pin}_${i}.py"
        sed -i "534s/^l=.*/l=$l/" "$WORK_DIR/SPV_p0_${p0}_v0_${v0}_pin_${pin}_${i}.py"
        sed -i "537s/^p0=.*/p0=$p0/" "$WORK_DIR/SPV_p0_${p0}_v0_${v0}_pin_${pin}_${i}.py"
        sed -i "538s/^v0=.*/v0=$v0/" "$WORK_DIR/SPV_p0_${p0}_v0_${v0}_pin_${pin}_${i}.py"
        sed -i "540s/^pinning_percent=.*/pinning_percent=$pin/" "$WORK_DIR/SPV_p0_${p0}_v0_${v0}_pin_${pin}_${i}.py"
        sed -i "545s/^set_i=.*/set_i=$i/" "$WORK_DIR/SPV_p0_${p0}_v0_${v0}_pin_${pin}_${i}.py"

        # 3. Modify the job.sh file, update the Python script name, and save with a new name
        cp "$JOB_SH" "$WORK_DIR/job_p0_${p0}_v0_${v0}_pin_${pin}_${i}.sh"
        sed -i "2s|^#SBATCH --job-name=.*|#SBATCH --job-name=${p0}_${v0}_${pin}_${i}|" "$WORK_DIR/job_p0_${p0}_v0_${v0}_pin_${pin}_${i}.sh"
        sed -i "21s|^python .*|python $WORK_DIR/SPV_p0_${p0}_v0_${v0}_pin_${pin}_${i}.py|" "$WORK_DIR/job_p0_${p0}_v0_${v0}_pin_${pin}_${i}.sh"

        # 4. Submit the job using qsub
        echo "Submitting job_p0_${p0}_v0_${v0}_pin_${pin}_${i}.sh to qsub"
        qsub "$WORK_DIR/job_p0_${p0}_v0_${v0}_pin_${pin}_${i}.sh"

        echo "Job for p0=$p0, v0=$v0 and i=$i submitted."
    done
done

echo "All jobs completed."
