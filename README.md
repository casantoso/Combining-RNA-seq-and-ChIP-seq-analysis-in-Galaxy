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

Run fastQC on ***SRA(ODO1)*** and ***SRA(Input)***. From the fastQC report, we can see that the ODO1 dataset have a high percentage of polyA sequence and the Input dataset have a high percentage of PolyG sequence. 

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
     - ```Adapter sequence to use```: TrueSeq3 (single-ended, for MiSeq and HiSeq)

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
     - ```Adapter sequence to use```: TrueSeq3 (single-ended, for MiSeq and HiSeq)

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
Download the petunia genome

Run ```Bowtie2``` twice: Once with ***trimmomatic2 (ODO1)*** as the input and once with ***trimmomatic2 (Input): paired*** as the input 

Use the following settings:
- ```Will you select a reference genome from your history or use a built-in index?```: Use a built-in genome index
     -```Select reference genome```: Mouse (mus musculus) : mm10
- ```Select analysis mode```
     -```Do you want to use presets?```: Very sensitive end-to-end
  

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







