# atacseq

## Introduction
This is a shell pipeline for automatically analyzing all paired-end ATAC-seq data in a folder.

It includes QC, trimming, contamination screening, mapping, filtering reads, peak calling, generating BigWig files and summary steps.

Downstream analyses, such as merging and quantifying chromatin accessibility for multiple samples, different peak analysis, peak annotation and pathway enrichment analaysis, are case specific, so are not included in this pipeline.

## Configuration
The paths of some genome files and tools need to be modified accordingly:
* Modify the paths of genome index files, ENCODE blacklist bed files and chromosome size files.
* Make sure **fastqc**, **multiqc**, **trim_galore**, **fastq_screen**, **bowtie2**, **samtools**, **bedtools**, **macs2**, **bedGraphToBigWig** are installed and can be ran by just typing their names, otherwise modify the paths of these tools in this pipeline script.
* Install **Picard** tools and set an environment variable/alias named "PICARD" for the path of picard.jar script, otherwise manually change "**$PICARD**" in this pipeline script to the full path of picard.jar script.
* Add summarize_ATACseq.pl to environment PATH or put it in /usr/local/bin/, otherwise modify the path of this perl file in this pipeline script.

## Usage
```
run_ATACseq.sh -g <hg19 or mm10> -i <path_of_inputs> -o <path_of_outputs> -t <threads>
```

## Arguments
* *-g*: set hg19 or mm10, default is hg19.
* *-i*: directory of input fastq or fastq.gz files.
* *-o*: directory for output files.
* *-t*: set average number of threads for **each** sample in parallel analysis, must be integer, default is 1.

## Outputs
* A **qc** folder including QC results.
* A **process** folder including a summary table and subfolders of analysis results for each sample.

## Contact
[Ting Wang](http://wt2015-github.github.io/) ([email](wang9ting@gmail.com)).
