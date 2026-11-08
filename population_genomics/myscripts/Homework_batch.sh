#!/bin/bash
#SBATCH --job-name=Admixture_HW
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Anna.E.Smith@uvm.edu

###change to where your data is
cd /gpfs1/home/a/s/asmit168/eco_genomics_2026/population_genomics/mydata/Homework_2

###this is where you should output your data 
OUTDIR=/gpfs1/home/a/s/asmit168/eco_genomics_2026/population_genomics/myresults

echo "This is my out directory"
echo $OUTDIR

module load admixture

INPUT="Lg_LD_0.5.bed"
K_MIN=1
K_MAX=10
THREADS=1   # we can increase this when running as a batch job

# start with an empty results file (>> below adds to the end of a file,
# so without this, re-running the loop would add duplicate rows)
rm -f cv_errors.tsv

for K in $(seq ${K_MIN} ${K_MAX}); do
    admixture --cv -j${THREADS} "${INPUT}" ${K} | tee "admixture_K${K}.out"
    CV=$(grep "CV error" "admixture_K${K}.out" | awk '{print $NF}')
    echo -e "${K}\t${CV}" >> cv_errors.tsv
done

