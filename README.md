# Combining-RNA-seq-and-ChIP-seq-analysis-using-Galaxy
This is an example workflow that combines ChIP-seq and RNA-seq analysis primarily using Galaxy. We are using a [dataset](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA729780&o=acc_s%3Aa) from [Boersma et al's study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9306810/) investigating the role of ODO1 in the production of FVBPs in Petunia flowers. FVBPs (floral volatile benzoid and phenylpropanoid compounds) are specialized metabolites that produce scent bouquets in Petunia flowers. ODO1 have been found to be a master regulator of metabolic FVBP synthesis and emission. However, the way in which ODO1 regulate the synthesis and emission of FVBPs in the Petunia was unknown. 

Thus, the  authors used *ChIP-seq* to identify the direct genomic targets of the ODO1 transcription factor in petunia flowers. Two transgenic petunia flowers were used :  
- pODO1:GFP‐ODO1
  -  uses the natural ODO1 promoter (pODO1)
  - crucial for observing when and where the protein is normally present—such as its diurnal oscillation (with levels peaking in the evening), reflecting the actual biological scenario in petunia flowers.

- 35S:GFP‐ODO1
  - 35S is a constitutive promoter from the cauliflower mosaic virus
  - This ensures that GFP‐ODO1 is expressed at high levels in a more uniform manner, regardless of the endogenous regulatory cues.
 
Using both the native promoter construct (pODO1:GFP‐ODO1) and the overexpression construct (35S:GFP‐ODO1) allowed the researchers to compare the binding profiles. Only those binding peaks detected in both setups were considered high-confidence, thereby minimizing the risk of false positives.

The authors used *RNA-seq* to determine whether the genes bound by ODO1 (as shown by ChIP‐seq) actually exhibited changes in expression when ODO1 was suppressed. This validation helps distinguish between direct targets (whose expression changes are directly linked to ODO1 binding) and indirect effects. RNA-seq was performed on WT and odo1i (ODO1 RNAi knockdown line) transgenic petunia flowers


## Table of contents

- ChIP-seq analysis
  - [Step 1: Import data](#step-1-import-data)
  - [Step 2: Quality control using FastQC](#step-2-quality-control-using-fastqc)
 - [Step 3: Trim using Trimmomatic ](#step-3-trim-using-trimmomatic)
 - [Step 4: Map reads to Petunia genome using Bowtie2](#step-4-map-reads-to-petunia-genome-using-bowtie2)
 - [Step 5: Merge Input files using MergeSamFiles](#step-5-Merge-Input-files-using-MergeSamFiles)
 - [Step 6: Find common peaks between the 2 petunia lines using Bedtools Intersect intervals](#step-6-Find-common-peaks-between-the-2-petunia-lines-using-Bedtools-Intersect-intervals)
 - [Step 7: bedtools Intersect intervals to find common peaks between the 2 petunia lines](#step-7-bedtools-Intersect-intervals-to-find-common-peaks-between-the-2-petunia-lines)
 - [Step 8: Gene Ontology](#step-8-gene-ontology)

- RNA-seq analysis
  - [Step 1: Import data](#step-1-import-data)
  - [Step 2: Quality control using FastQC](#step-2-quality-control-using-fastqc)
  - [Step 3: Trim using Trimmomatic ](#step-3-trim-using-trimmomatic)
  - [Step 4: Align reads to Petunia genome using HISAT2](#step-4-Align-reads-to-Petunia-genome-using-HISAT2)
  - [Step 5: Count the number of aligned reads that overlapp Petunia gff annotation file using Htseq-counts](#step-5-Count-the-number-of-aligned-reads-that-overlapp-Petunia-gff-annotation-file-using-Htseq-counts) 
  - [Step 6: Find differentially expressed genes using DESeq2](#step-6-Find-differentially-expressed-genes-using-DESeq2) 
  - [Step 7: Gene Ontology](#step-7-gene-ontology)

- Combining ChIP-seq and RNA-seq analyses
  - [Step 1: Create 3 files of filtered gff based off of the 3 gene lists](#step-1-Create-3-files-of-filtered-gff-based-off-of-the-3-gene-lists) 
  - [Step 2: Motif analysis using memeChIP](#step-9-motif-analysis-using-memeChIP)
  -  


### ChIP-seq analysis
#### Step 1: Import data
From this [website](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA729780&o=acc_s%3Aa), click on the boxes next to the SRR numbers for the Input and IP data for Ctcf homozygous mutation (9 and 15 on the list) and then press the the galaxy button shown in the picture below. Do the same for the Input and IP data for Wildtype (10 and 16 on the list)


#### Step 2: Quality control using ```FastQC```

#### Step 3: Trim using ```Trimmomatic``` 

#### Step 4: Map reads to Petunia genome using ```Bowtie2```

#### Step 5: Merge Input files using ```MergeSamFiles``` 

#### Step 6: Find peaks using ```MACS2 callpeak``` 

#### Step 7: Find common peaks between the 2 petunia lines using ```Bedtools Intersect intervals``` 

#### Step 8: Gene Ontology



### RNA-seq analysis
#### Step 1: Import data

#### Step 2: Quality control using ```FastQC```

#### Step 3: Trim using ```Trimmomatic``` 

#### Step 4: Align reads to Petunia genome using ```HISAT2```

#### Step 5: Count the number of aligned reads that overlapp Petunia gff annotation file using ```Htseq-counts```

#### Step 6: Find differentially expressed genes using ```DESeq2``` 

#### Step 7: Gene Ontology




### Combining ChIP-seq and RNA-seq analyses
#### Step 1: Create 3 files of filtered gff based off of the 3 gene lists

#### Step 2: Motif analysis using ```memeChIP```







