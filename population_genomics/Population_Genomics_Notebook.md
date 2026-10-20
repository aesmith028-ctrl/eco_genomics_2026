# Population Genomics Notebook

**Course:** Intro to Ecological Genomics - Fall 2026

**Name:** Anna Smith

------------------------------------------------------------------------

## 10/13/2026 - Working with VCF files

-   Describe the main steps that take raw sequencing reads to a filtered VCF

-   Read the header, fixed columns and genotype fields of a VCF

-   Query a VCF with bcftools (header, samples, records and regions)

-   Calculate per-individual missingness and depth with vcftools

-   Write, submit and monitor a batch job on the VACC

**Working Directory**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/population_genomics`

**Input Files**

`none`

**Output Files**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/population_genomics/Population_Genomics_Notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

-   `load ecogen-rlibs`

**Scripts:**

`none; all code run in terminal`

**Code:**

``` bash
#### What is a VCF? ####
# A VCF (Variant Call Format) file stores the positions in the genome where variants were found among our samples, along with every individual’s genotype at those positions. It has two parts:

# a header: lines starting with ## that record how the file was made and define every field used in it
# a body: one line per variant
# The full specification is at https://samtools.github.io/hts-specs/VCFv4.2.pdf.

## Pipeline ##

# raw sequences (.fastq)
# Quality Control (QC, fastp)
# Alignment/mapping (bwa)
# Mapping Quality (picard)
# variant calling (samtools/bcftools)
# SNP filtering (bcftools)

#### From reads to alignments ####
# 1. Quality control and trimming. Raw reads arrive from the sequencer as FASTQ files. We check read quality with FastQC (MultiQC combines many FastQC reports into one), then remove adapters and trim low-quality bases with a trimmer such as fastp (fastp also does the same as FastQC, but it is sometimes good to do both!).
# 2. Mapping. The cleaned reads are mapped to the reference genome (or assembled de novo if there is no reference). We used BWA (the Burrows–Wheeler Aligner), which indexes the reference using the Burrows–Wheeler transform so it can work out where millions of short reads belong quickly and with little memory.
# 3. Cleaning up the alignments. We then remove reads that would bias the variant calls, for example PCR duplicates (identical copies of the same DNA fragment made during library amplification) and reads that map poorly.

# mapping = approximate, alignment = exact #

#### CODING DONE IN BASH ####
 # copy over 4 files from /gpfs1/cl/biol3990/PopulationGenomics/VCF to /gpfs1/home/a/s/asmit168/eco_genomics_2026/population_genomics/mydata

# run zcat Lottia_gigantea_biallelic.vcf.gz | head -n 20 
  # There is one ##contig line for every sequence in the reference genome. Our reference, the owl limpet (Lottia gigantea) assembly, is scaffold-level: roughly 360 Mb split into thousands of scaffolds with no assembled chromosomes (4469 to be precise), so the header contains thousands of these lines.
  # A contig is a continuous stretch of DNA sequence assembled from overlapping reads, with no gaps. Contigs are ordered and joined into scaffolds (which can contain gaps of unknown sequence, written as Ns), and in the best assemblies scaffolds are placed onto chromosomes. How far an assembly gets depends on its quality. The VCF header calls every reference sequence a “contig”, whichever of these it actually is.

#### Reading the header with bcftools ####
  # A better way to look at a VCF is with our first tool, bcftools. Load it on the VACC first, and always record which version you are using:
module load gcc
module load bcftools

bcftools --version   # bcftools 1.19, Using htslib 1.19.1 (ALWAYS RECORD VERSIONS)
bcftools -h          # -h or --help usually gets you the manual

bcftools head Lottia_gigantea_biallelic.vcf.gz

# Again we get a big string of contig lines. Let’s pipe (|) the output into grep to remove them:
bcftools head Lottia_gigantea_biallelic.vcf.gz | grep -v "##contig="
  # grep (“global regular expression print”) searches text for a pattern and prints the lines that match. The -v option inverts this, printing only the lines that don’t match.

# Better. Now you can see how the file was made and what each field means. Each ##INFO and ##FORMAT line defines one of the tags used in the body of the file:
bcftools view -h Lottia_gigantea_biallelic.vcf.gz | grep "^##INFO"
bcftools view -h Lottia_gigantea_biallelic.vcf.gz | grep "^##FORMAT"

# bcftools also writes every command it runs into the header (in lines containing Command=), so you can often reconstruct how a file was made and filtered:
bcftools view -h Lottia_gigantea_biallelic.vcf.gz | grep "Command="

#### Reading a variant record ####
  # Let’s look at the first few SNPs to see what our data look like. The -H option drops the header so we only see variant records:
  
# View the first SNP
bcftools view -H Lottia_gigantea_biallelic.vcf.gz | head -n 1

# the first 5 records; less -S stops long lines wrapping (press q to exit)
bcftools view -H Lottia_gigantea_biallelic.vcf.gz | head -n 5 | less -S
```

**Filtering The VCF**

| Filter | Threshold | Why |
|------------------|------------------|------------------------------------|
| Variant type | Biallelic SNPs only | Indels and sites with more than two alleles are harder to call reliably, and many are sequencing errors |
| Site Quality (QUAL) | ≥ 30 | A Phred score of 30 means a 1 in 1,000 chance that the call is wrong |
| Genotype Depth | ≥ 3 reads | Genotypes supported by very few reads are unreliable, so they are set to missing |
| Missingness | Remove sites missing in \> 30% of individuals | Sites with lots of missing data add noise and can bias downstream analyses |
| Minor Allele Frequency (MAF) | Commonly 1% or 5% | Very rare alleles are hard to distinguish from sequencing errors; 5% is the stricter cut-off |

------------------------------------------------------------------------

## 10/15/2026 - Working with VCF files

-   Describe the main steps that take raw sequencing reads to a filtered VCF

-   Read the header, fixed columns and genotype fields of a VCF Query a VCF with bcftools (header, samples, records and regions)

-   Calculate per-individual missingness and depth with vcftools

-   Write, submit and monitor a batch job on the VACC

**Working Directory**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/population_genomics`

**Input Files**

`none`

**Output Files**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/population_genomics/myscripts/Batch1.sh`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

-   `load ecogen-rlibs`

**Scripts:**

`none; all code run in terminal`

**Code:**

``` bash
# Where did we get the data
  # SRA - sequence read archive

#### fixing wrong data from tuesday ####
# delete previous data from tuesday
# cp (4 data) ~/eco_genomics_2026/population_genomics/mydata

#### working in terminal with VCF data ####
module purge
module load gcc
module load bcftools
module load vcftools

bcftools --version   # bcftools 1.19 Using htslib 1.19.1
bcftools -h          # -h or --help usually gets you the manual

bcftools head Lottia_gigantea_biallelic.vcf.gz

## Query tool ##

# bcftools query -f prints only the fields you ask for, which is how we will get data out of a VCF and into other programs. Fields are written as %NAME, \t is a tab, \n is a new line, and anything in square brackets is repeated for every sample:

# scaffold, position, alleles and every individual's genotype for the first 5 SNPs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
    Lottia_gigantea_p1_flt3.renamed.vcf.gz | head -n 5

#### Running a Batch Job ####
# type nano in terminal to get text editor
  # gets you more memory so does stuff fast
  # first lines tell computer specifics of job and where it is stored, partition will always be general (big computer), how many nodes are needed (always 1 node), cpu per task may change depending on data size (1 should be fine), CHANGE MEMORY DEPENDING ON DATA (small task = less than 5 gigs of data), time should change depending on analysis, output and error stay same, mail type emails when job starts and ends, TEMPLATE BELOW to use every time, 
```

**BATCH TEMPLATE**

#!/bin/bash

#SBATCH --job-name=vcf_stats

#SBATCH --partition=general #SBATCH --nodes=1

#SBATCH --cpus-per-task=1 #SBATCH --mem=40G

#SBATCH --time=30:00:00 #SBATCH --output=./logs/%x.%j.out

#SBATCH --error=./logs/%x.%j.err

#SBATCH --mail-type=ALL

#SBATCH --[mail-user=NETID\@uvm.edu](mailto:mail-user=NETID@uvm.edu){.email}

------------------------------------------------------------------------

## 10/20/2026 - Population Structure

-   explain what linkage disequilibrium (LD) is and why we prune it before looking at population structure

-   LD-prune a VCF with plink

-   run a principal component analysis (PCA) with plink and plot it in R

-   run ADMIXTURE across a range of K and use cross-validation to choose which values to look at

-   interpret PCA and ADMIXTURE plots, and recognise their common pitfalls

**Working Directory**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/population_genomics`

**Input Files**

`none`

**Output Files**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/Population_Genomics_Notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

-   `load ecogen-rlibs`

**Scripts:**

`none; all code run in terminal`

**Code:**

``` bash
#### Setting Up ####
# First, make a folder for today’s work and move into it.

# Then store the location of the VCF in a variable, so we don’t have to type the full path every time:
mkdir -p ~/eco_genomics_2026/population_genomics/myresults/Population_Structure 

# where the VCF lives (make sure this is the right path!)
DATA_DIR=~/eco_genomics_2026/population_genomics/mydata
VCF=$DATA_DIR/Lottia_gigantea_pl_flt3.renamed.vcf.gz

# check the variable holds what you expect, and that the file is really there
echo $VCF
ls -lh $VCF

cd $DATA_DIR
pwd ##check you are in the right place
```
