README for MixClone 1.1.1
=========================


INTRODUCTION
============

MixClone is a comprehensive software package to study the subclonal 
structures of tumor genomes, including subclonal cellular prevalences 
estimation, allelic configuration estimation, absolute copy number 
estimation and a few visualization tools. It takes next-generation 
sequencing data of paired normal-tumor samples as input and integrates 
sequence information from both somatic copy number alterations and allele
frequencies within a unified probabilistic framework. If you have any 
questions, please email yil8@uci.edu



INSTALL
=======

Prerequisites
-------------
**Mandatory** 

* Python (2.7). [Python 2.7.3](http://www.python.org/download/releases/2.7.3/) is recommended.

* [Numpy](http://www.numpy.org/)(>=1.6.1). You can download the source of Numpy from [here](http://sourceforge.net/projects/numpy/files/).

* [Scipy](http://www.scipy.org/)(>=0.10). You can download the source of Scipy from [here](http://sourceforge.net/projects/scipy/files/).

* [Pysam](https://code.google.com/p/pysam/)(>=0.7). To install Pysam, you also need to install [Cython](http://cython.org/) first. 


**Optional** 
* [matplotlib](http://matplotlib.org/)(>=1.2.0) is required for a few visualization tools.


Although not mandatory, Linux system is recommended. Also, [samtools](http://samtools.sourceforge.net/) is not required by MixClone, but can be useful for creating bam, bam index and fasta index files which are required by the pysam module of MixClone. 

Install from source
-------------------
Download the compressed source file MixClone-*.tar.gz and do as follows:

```
$ tar -xzvf MixClone-*.tar.gz
$ cd MixClone-*
$ python setup.py install
```

If you prefer to install MixClone other than the default directory, you can also use this command:
```
$ python setup.py install --prefix /home/yili/
```

There is also a `bin/` folders under MixClone-*. The `bin/` folder contains useful utilities, such as the R code to run [BICseq](http://compbio.med.harvard.edu/Supplements/PNAS11.html) and the python script to convert BICseq results to BED file. You can copy these two folders somewhere easily accessible.



USAGE
=====


Overview
--------
MixClone is composed of three modules: 
* `preprocess`. Preprocess the reads aliments of paired normal-tumor samples in BAM format, the tumor genome segmentation file in BED format, and produce the *.MixClone.input.pkl file as output, which will be used for running the model.
 
* `run_model`. Take the *.MixClone.input.pkl as input, estimate the subclonal cellular prevalence, the absolute copy number and the allelic configuration of each segment, and produce the *.MixClone.output.pkl file as output, which will be used for postprocessing. If the user runs the model without specifying the number of subclonal populations, MixClone will run the model five times with subclonal number ranges from 1 to 5, and recommend the most likely model.

* `postprocess`. Take the *.MixClone.output.pkl file as input, and extract various output files, including the segments file with estimated parameters, the allele counts file, the subclonal analysis summary and a few plots.

The general workflow of MixClone is this
![alt tag](https://github.com/uci-cbcl/MixClone/blob/gh-pages/images/workflow.png)



Tumor genome segmentation
-------------------------
MixClone requires a segmentation file of the tumor genome in BED format before running the package. We used [BICseq](http://compbio.med.harvard.edu/Supplements/PNAS11.html) in the original paper. To run a BICseq analysis, you can copy the commands in `bin/BICseq.R` (Li, Y., Xie, X. 2014, Bioinformatics) and paste them in a R interative shell. Or you can also run the R script from the command line:
```
$ R CMD BATCH bin/BICseq.R
```
Note that,`normal.bam` and `tumor.bam` must be in the same directory where you run the command. The R script will output a segments file `segments.BICseq`. Then you can use the other script `bin/BICseq2bed.py` (Li, Y., Xie, X. 2014, Bioinformatics) to convert the segments file into BED format:
```
$ BICseq2bed.py segments.BICseq segments.bed --seg_length 1000000
```

**--seg_length** Only convert segments with length longer than the threshold.



Preprocess
----------
This part of README is based on [JoinSNVMix](https://code.google.com/p/joint-snv-mix/wiki/runningOld). To preprocess the paired cancer sequencing data, execute:
```
$ MixClone.py preprocess REFERENCE_GENOME.fasta SEGMENTS.bed NORMAL.bam TUMOUR.bam INPUT_BASENAME --min_depth 20 --min_base_qual 10 --min_map_qual 10 --process_num 10
```

**REFERENCE_GENOME.fasta** The path to the fasta file that the paired BAM files aligned to. Currently, only the
[UCSC](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) and [ENSEMBL](http://uswest.ensembl.org/info/website/upload/bed.html) chromosome format are supported. Note that the index file should be generated for the reference genome. This can be done by running samtools as follows:

`$ samtools faidx REFERENCE_GENOME.fasta`

**SEGMENTS.bed** The BED file for the tumor genome segmentation.

**NORMAL.bam** The BAM file for the normal sample. The BAM index file should be generated for this file and named NORMAL.bam.bai. This can be done by running

`$ samtools index NORMAL.bam`

**TUMOUR.bam** The bam file for the tumour sample. As for the normal this file needs to be indexed.

**INPUT_BASENAME** The base name of the preprocessed input file to be created.

**--min_depth** Minimum depth in both normal and tumor sample required to use a site in the analysis.

**--min_base_qual** Minimum base quality required for each base.

**--min_map_qual** Minimum mapping quality required for each base.

**--process_num** Number of processes to launch for preprocessing.



Run model
---------
After the preprocessed input file is created, we can run the probabilistic model of MixClone by execute:
```
$ MixClone.py run_model INPUT_BASENAME OUTPUT_BASENAME --max_copynumber 6 --subclone_num 2 --max_iters 30 --stop_value 1e-6
```
**INPUT_BASENAME** The base name of the preprocessed input file created in the preprocess step.

**OUTPUT_BASENAME** The base name of the output file with model parameters estimated to be created.

**--max_copynumber** The maximum copy number of each segment allows to take.

**--subclone_num** The number of subclones within the tumor sample. If not provided, go through [1, 5] and select the most likely model.

**--max_iters** Maximum number of iterations for training.

**--stop_value** Stop value of the EM algorithm for training. If the change of log-likelihood is lower than this value, stop training.



Postprocess
-----------
After the output file with model parameters estimated, we can extract various result files from the output file by execute:
```
$ MixClone.py postprocess OUTPUT_BASENAME
```

**BASENAME** The base name of the output file created in the run_model step.



Output files
------------
**\*.MixClone.segments** The segments file. It contains the genomic and subclonal information of each segment. The definition of each column in a *.MixClone.segments file is listed here:

| Column           | Definition                                                              | 
| :--------------- | :---------------------------------------------------------------------- | 
| seg_name         | Name of the segment                                                     |      
| chrom            | Chromosome of the segment                                               |  
| start            | Start position of the segment                                           |
| end              | End position of the segment                                             |
| normal_reads_num | Count of reads mapped to the segment in the normal sample               |
| tumor_reads_num  | Count of reads mapped to the segment in the tumor sample               |
| LOH_frac         | Fraction of LOH sites in the segment                                    |
| LOH_status       | FALSE -> no LOH; TRUE -> significant LOH; UNCERTAIN -> medium level LOH |
| log2_ratio       | Log2 ratio between tumor_reads_num and normal_reads_num                 |
| copy_number      | Estimated absolute copy number of the segment                           |
| allele_type      | Estimated allelic configuration of the segment                          |
| subclone_prev    | Estimated subclonal cellular prevalence of the segment                  |
| subclone_cluster | Estimated subclonal cluster label of the segment                        |

**\*.MixClone.counts** The allele counts file. It contains the allelic counts information of sites, which are heterozygous SNP sites in the normal genome. The definition of each column in a *.MixClone.counts file is listed here:

| Column    | Definition                                         | 
| :-------- | :------------------------------------------------- | 
| seg_name  | Name of the segment                                |      
| normal_A  | Count of bases match A allele in the normal sample |
| normal_B  | Count of bases match B allele in the normal sample |
| tumor_A   | Count of bases match A allele in the tumor sample  |
| tumor_B   | Count of bases match B allele in the tumor sample  |
| chrom     | Chromosome of the site                             |
| pos       | Genomic position of the site                       |

**\*.MixClone.summary** The summary file about the subclonal analysis, including log-likelihood, subclonal labels and the corresponding cellular prevalences.

**\*.MixClone.heatmap** The folder of BAF heat maps plotted for each segment (Li, Y., Xie, X. 2014, Bioinformatics).

**\*.MixClone.segplot.png** The plot of estimated subclonal cellular prevalences and absolute copy numbers of all the segments of non-diploid allelic configuration.

**\*.MixClone.model_selection** The summary about model selection if `--subclone_num` is not specified. It includes the log-likelihood and related change of each model with different number of subclones. Although MixClone selects the most likely model based on a heuristic described in the original paper, users can select the best model use their own judgement based on the log-likelihood related information.

