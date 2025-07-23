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
  -  [Step 7: bedtools Intersect intervals to find common peaks between the 2 petunia lines](#step-7-bedtools-Intersect-intervals-to-find-common-peaks-between-the-2-petunia-lines)
  -  [Step 8: Gene Ontology](#step-8-gene-ontology)

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
We will first download the IP dataset, from this [website](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA729780&o=acc_s%3Aa), click on the boxes next to the SRR number for the 	
ChIP-seq data (SRR14528049 and SRR14528050) and then press the the galaxy button shown in the picture below. This will automatically bring you to Galaxy (Note that you will have to create a Galaxy account to get enough memory to do this analysis). Name this SRA collection ***SRA(ODO1)***. 

Then, we will download the background control(input) dataset. From this [website](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA650505&o=acc_s%3Aa), click on the boxes for SRR12442821,  SRR12442822 and SRR12442825,  then press the Galaxy button. Name this SRA collection ***SRA(Input)***. 

Then go to ```tools``` → ```Get Data``` → ```Download and Extract Reads in FASTQ format from NCBI SRA```

Use the following settings:
* ```select input type```: list of SRA accession, one per line
* Under ```sra accession list```, input your SRA collection (i.e. ***SRA(ODO1)***/ ***SRA(Input)***)
* ```select output format```: gzip compressed fastqc
Then press ```Run Tool```. Run it twice, once for each SRA collection in your history


After it has finished running, for the ODO1 dataset, you should see ***a list with 2 fastqsanger.gz pairs*** under ***Single-end data (fastq-dump)*** and ***a list with 0 datasets*** under ***Paired-end data (fastq-dump)***. For, the input dataset, you will see ***a list with 3 fastqsanger.gz pairs*** under ***Paired-end data (fastq-dump)*** and ***a list with 0 datasets*** under ***Single-end data (fastq-dump)***

![SRA](chipseq_img/1-SRA.png)

#### Step 2: Quality control using ```FastQC```

Run fastQC on ***Single-end data (ODO1)*** and ***Paired-end data (Input)***. From the fastQC report, we can see that the ODO1 dataset have a high percentage of polyA sequence and the Input dataset have a high percentage of PolyG sequence. 

![ODO1_fastQC_before_trimming](chipseq_img/2-ODO1_fastQC_before_trimming.png)
![Input_fastQC_before_trimming](chipseq_img/3-Input_fastQC_before_trimming.png)

#### Step 3: Trim using ```Trimmomatic``` 

Run ```Trimmomatic``` twice on each dataset collection

For the ODO1 dataset:
##### For the 1st run:
- input: ```Single-end or paired-end reads?```: single-end
  - ***SRA(ODO1)***
- Use the following settings:
     - ```Perform initial ILLUMINACLIP step?``` : Yes
          - ```Select standard adapter sequences or provide custom?```: Standard
     - ```Adapter sequence to use```: TruSeq3 (single-ended, for MiSeq and HiSeq)

##### For the 2nd run:
- input: ```Single-end or paired-end reads?```: single-end
  - ***SRA(ODO1)***
- Use the following settings:
     - ```Perform initial ILLUMINACLIP step?``` : Yes
          - ```Select standard adapter sequences or provide custom?```: custom
     - ```Adapter sequence```:
  ```
          > polyA
            AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
  ```

 
For the Input dataset:
##### For the 1st run:
- input: ```Single-end or paired-end reads?```: paired-end (as a collection)
  - ***SRA(Input)***
- Use the following settings:
   ```Perform initial ILLUMINACLIP step?``` : Yes
          - ```Select standard adapter sequences or provide custom?```: Standard
     - ```Adapter sequence to use```: TruSeq3 (paired-ended, for MiSeq and HiSeq)

##### For the 2nd run:
- input: ```Single-end or paired-end reads?```: paired-end (as a collection)
  - ***SRA(Input)***
- Use the following settings:
     - ```Perform initial ILLUMINACLIP step?``` : Yes
          - ```Select standard adapter sequences or provide custom?```: custom
     - ```Adapter sequence```:
  ```
          > polyA
            GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
  ```


For the output of trimmomatic for the Input dataset, we will get paired and unpaired, we will use the paired output. 

![trimmomatic_outputs](chipseq_img/4-trimmomatic_outputs.png)

Now we will only be using ***trimmomatic2 (ODO1)*** and ***trimmomatic2 (Input): paired***

Then run fastQC on ***trimmomatic2 (ODO1)*** and ***trimmomatic2 (Input): paired*** to see whether trimmomatic has succesfully trimmed out the adapter sequence and polyG sequence. Now in the fastQC report, we can see that we have trimmed out the adapter sequences. 

![Input_fastQC_after_trimming](chipseq_img/5-Input_fastQC_after_trimming.png)
![ODO1_fastQC_after_trimming](chipseq_img/6-ODO1_fastQC_after_trimming.png)


#### Step 4: Map reads to Petunia genome using ```Bowtie2```
Download the petunia genome fasta file  and upload it to Galaxy

Run ```Bowtie2``` twice: Once with ***trimmomatic2 (ODO1)*** as the input and once with ***trimmomatic2 (Input): paired*** as the input 

Use the following settings:
- ```Will you select a reference genome from your history or use a built-in index?```: Use a genome from history and build indes
     -```Select reference genome```: your petunia genome
- ```Select analysis mode```
     -```Do you want to use presets?```: No, just use defaults
- ```Save the bowtie2 mapping statistics to the history``` : Yes
  
you will get 2 Bowtie2 outputs for each dataset: alignments and mapping stats.I named the outputs ***Bowtie2(Input):alignments***, ***Bowtie2(Input):mapping stats***, ***Bowtie2(ODO1):alignments***, ***Bowtie2(ODO1):mapping stats***. 
From the mapping stats output you can see how mappy percentage of reads were mapped. Generally a percentage above 80% is considered good. 


#### Step 5: Merge Input files using ```MergeSamFiles``` 
We have 3 Input samples in our dataset colleciton so we want to combine them into 1 file for cleaner, consolidated background control in MACS2 peak calling

Use the following settings:
- ```Select SAM/BAM dataset or dataset collection``: Bowtie2(Input):aligments

I named the output: ***MergeSamFiles on input***


#### Step 6: Find peaks using ```MACS2 callpeak``` 
```MACS2 callpeak``` is used to identify enriched regions of DNA — called "peaks" — from ChIP-seq data. It basically looks for places in the genome where there are many sequencing reads aligned to it which indicates where CTCF binds. 


Use these settings:
-```ChIP-Seq Treatment File``` : ***Bowtie2(ODO1):alignments***
- ```Do you have a Control File?```: yes
   - ```ChIP-Seq Control File``` : ***MergeSamFiles on input***


Use the following settings:
- ```Format of Input Files``` : Single-end BAM
- ```Effective genome size``` :1200000000


#### Step 7: Find common peaks between the 2 petunia lines using ```Bedtools Intersect intervals``` 
First we have to extract the 2 datasets from ***MACS2 callpeak (narrow peaks)***

Run ```Extract dataset``:
- ```Input List``` : ***MACS2 callpeak (narrow peaks)***
     - ```How should a dataset be selected?```: Select by index
     -``` Element index```: 0 and 1 //Run once with element index as 0 and run once with element index as 1.

For the output with element index as 0 (SRR14528049), name it ***MACS2 callpeak (35S:GFP-ODO1)***. For the output with element index as 1 (SRR14528050), name it ***MACS2 callpeak (pODO1:GFP-ODO1)***

Then we want to find the peaks that are common between the 2 lines: pODO1:GFP-ODO1 and 35S:GFP-ODO1. 
Use the following settings:
- ```File A to intersect with B ```: ***MACS2 callpeak (35S:GFP-ODO1)***
-  ```Combined or separate output files ```: One output file per 'input B' file
  -  ```File B to intersect with A ```: ***MACS2 callpeak (pODO1:GFP-ODO1)***
-  ```Calculation based on strandedness? ```: Overlaps on either strand

I named the output ***bedtools Intersect intervals on pODO1:GFP-ODO1 and 35S:GFP-ODO1***

Then we will run ```Bedtools Intersect intervals```  again to find the intersection between the peaks found in ***bedtools Intersect intervals on pODO1:GFP-ODO1 and 35S:GFP-ODO1*** and the Petunia GFF annotation. This is  done to identify which genes are likely direct targets of ODO1.


#### Step 8: Gene Ontology



### RNA-seq analysis
#### Step 1: Import data
We will first download the IP dataset, from this [website](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA729780&o=acc_s%3Aa), click on the boxes next to the SRR numbers for the 	
RNA-seq data (SRR14528051, SRR145280512, SRR14528053, SRR14528054, SRR14528055, SRR14528056) and then press the the galaxy button shown in the picture below.

Then go to ```tools``` → ```Get Data``` → ```Download and Extract Reads in FASTQ format from NCBI SRA```

Use the following settings:
* ```select input type```: list of SRA accession, one per line
* ```sra accession list```: ***SRA***
* ```select output format```: gzip compressed fastqc
Then press ```Run Tool```. Run it twice, once for each SRA collection in your history

You will see ***a list with 6 fastqsanger.gz pairs*** under ***Paired-end data (fastq-dump)*** and ***a list with 0 datasets*** under ***Single-end data (fastq-dump)***. 
![SRA](rnaseq_img/1-SRA.png)

#### Step 2: Quality control using ```FastQC```
Run fastQC on ***Paired-end data (fastq-dump)***. From the fastQC report, we can see the adapter content is very low which is good, however, we can still run trimmomatic to remove the small percentage of Ilumina universal adapter. 

![fastQC_before_trimming](rnaseq_img/2-fastQC_before_trimming.png)

#### Step 3: Trim using ```Trimmomatic``` 
Use the following settings: 
- input: ```Single-end or paired-end reads?```: paired-end (as a collection)
  - ***SRA(Input)***
- Use the following settings:
   ```Perform initial ILLUMINACLIP step?``` : Yes
          - ```Select standard adapter sequences or provide custom?```: Standard
     - ```Adapter sequence to use```: TruSeq3 (paired-ended, for MiSeq and HiSeq)
 
we are only going to use the apired output( ***Trimmomatic:paired***)

Then run fastQC on  ***Trimmomatic:paired*** to see whether trimmomatic has succesfully trimmed out the adapter sequence. Now in the fastQC report, we can see that we have trimmed out the adapter sequences. 
![fastQC_after_trimming](rnaseq_img/2-fastQC_after_trimming.png)

#### Step 4: Align reads to Petunia genome using ```HISAT2```
HISAT2 employs a graph FM index and hierarchical indexing strategy which means that it employs two types of indexes: (1) one global FM index representing the whole genome, and (2) many separate local FM indexes for small regions collectively covering the genome. 

The way it works is that it builds a genome index and a splicing graph index (from a GTF/GFF). This lets HISAT2 search efficiently through known exon-exon junctions and novel junctions. When a read fails to align contiguously, HISAT2 attempts to align across known or inferred introns by referencing its splicing graph. It can match reads even if they’re split across multiple exons, using either annotation or read evidence. This makes HISAT2 extremely efficient for genomes with alternative splicing

Use the following settings:
- ```Source for the reference genome```: cDNA.fasta file
- ```Is this a single or paired library```: Paired-end Dataset Collection
  - ```Paired Collection ```: ***Trimmomatic:paired***
 
There are 2 outputs ***HISAT2: Mapping summary*** and ***HISAT2: aligned reads (BAM)***
From ***HISAT2: Mapping summary*** you can see what percentage of the reads were mapped. 


#### Step 5: Count the number of aligned reads that overlapp Petunia gff annotation file using ```Htseq-counts```
HTSeq-count takes an alignment file (typically in SAM or BAM format) containing sequencing reads mapped to a reference genome and an annotation file (usually in GTF or GFF format) that describes genomic features (such as genes or exons). Its main function is to count the number of reads that map to each genomic feature. HTSeq-count outputs 2-column tab-delimited text file where the first column is the Gene ID and the 2nd column is the counts. There are 5 rows in the bottoms that do not correspond to geneIDs 
 - __no_feature	->  reads that didn’t overlap any feature (e.g. intergenic)
 - __ambiguous	-> reads that overlapped multiple genes and can’t be assigned
 - __alignment_not_unique ->	reads that map to multiple locations in the genome
 - __not_aligned	-> reads that did not map at all
 - __too_low_aQual	-> reads with low mapping quality (if you enabled filtering)

Use the following settings:
- ```Aligned SAM/BAM File``` : ***HISAT2: aligned reads (BAM)***
- ```GFF/GTF File```: the petunia genome gff file
- ```Mode``` : Union
- ```Minimum alignment quality```: 30
- ```ID Attribute```: ID

I named the output ***htseq-count***


#### Step 6: Find differentially expressed genes using ```DESeq2``` 
Before we run ```DESeq2```  we are going to extract all of the different sata from 

#### Step 7: Gene Ontology




### Combining ChIP-seq and RNA-seq analyses
#### Step 1: Create 3 files of filtered gff based off of the 3 gene lists

#### Step 2: Motif analysis using ```memeChIP```







