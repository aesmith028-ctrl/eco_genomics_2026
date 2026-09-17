# Transcriptomics Notebook

**Course:** Intro to Ecological Genomics - Fall 2026

**Name:** Anna Smith

------------------------------------------------------------------------

## 9/15/2026 - Setting up lab notebook and learning markdown

-   setting up transcriptomics notebook

-   Learn how to take notes in markdown

-   push notes to github

**Working Directory**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics`

**Input Files**

`none`

**Output Files**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/Transcriptomics_Notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`none`

**Code:**

``` r
print("Hello world")
```

**Table:**

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |

![](markdown_cheat_sheet.png){width="588"}

**Notes/Observations**:

oh cool graph, yay! And some interpretation of what you are seeing, maybe what the next step could be?

**Next Steps?**

-   Lets dive into the project!

------------------------------------------------------------------------

## 9/17/2026 - Introduce the study system

-   Understand the experimental design and data and questions that can be addressed with these data.

-   Understand what a .fastq file is.

-   Understand the work flow or “pipeline” for processing and analyzing RNAseq data.

-   Conceptualize a “counts” file.

**Working Directory**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics`

**Input Files**

`none`

**Output Files**

`/gpfs1/home/a/s/asmit168/eco_genomics_2026/transcriptomics/Transcriptomics_Notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts:**

`none`

**Code:**

``` r
# Useful commands for later in shell access of VACC

# to get to class directory
cd /gpfs1/cl/biol3990

# to get list of what is in directory (ls short, ll long)
ls
ll

# take the output of zcat (reads a gzipped file; AA...qz file) and give me (|) this number of lines (-n 8) - output is the first few lines of the quality of the reads
zcat AA_F0_Rep3_2_clean.fq.gz | head -n 8

#DO NOT run the above file without adding | (pipe) aka filtering it thru a secondary command.... or else computer may explode

# tell me how many lines are in this file
zcat AA_F0_Rep3_2_clean.fq.gz | wc -l
```
