#!/bin/bash
#SBATCH --job-name=pca
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=./logs/%x.%j.out
#SBATCH --error=./logs/%x.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=asmit168@uvm.edu

## change to where data is
/users/a/s/asmit168/eco_genomics_2026/population_genomics/mydata

## this is where output your data
OUTDIR=/users/a/s/asmit168/eco_genomics_2026/population_genomics/myresults/Population_Structure

echo $VCF
### Call VCF
VCF=$DATA_DIR/Lottia_gigantea_pl_flt3.renamed.vcf.gz

## Load plink
module load plink

## LD thin first
plink --vcf $VCF --double-id --allow-extra-chr \
    --set-missing-var-ids @:# \
    --indep-pairwise 50 10 0.2 \
    --out $OUTDIR/Lg_LD_0.2

## Now lets run PCA
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
    --extract $OUTDIR/Lg_LD_0.2.prune.in \
    --make-bed --pca --out $OUTDIR/Lg_LD_0.2

echo "my script is done"