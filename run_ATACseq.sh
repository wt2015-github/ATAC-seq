#!/bin/bash
#Ting Wang, 201703

manual="
Shell script for automatically analyzing all paired-end ATAC-seq data in a folder.
QC with fastqc at the beginning and use multiqc to summarize qc figures.
Usage:
    run_ATACseq.sh -g <hg19 or mm10> -i <path_of_inputs> -o <path_of_outputs> -t <threads>
Options:
    -g    set hg19 or mm10, default is hg19
    -i    directory of input fastq or fastq.gz files
    -o    directory for output files
    -t    average number of threads for each sample, must be integer, default is 1
"

if [[ $# -le 0 ]]; then
    echo "${manual}" >&2
    exit 1
fi

while getopts :g:i:o:t: ARGS
do
case $ARGS in
    g)
        genome=$OPTARG
        ;;
    i)
        pathin=$OPTARG
        ;;
    o)
        pathout=$OPTARG
        ;;
    t)
        threads=$OPTARG
        ;;
    :)
        echo "no value for option: $OPTARG"
        echo "${manual}" >&2
        exit 1
        ;;
    *)
        echo "unknow option: $OPTARG"
        echo "${manual}" >&2
        exit 1
        ;;
esac
done

echo ${genome:="hg19"} >/dev/null
echo ${threads:=1} >/dev/null
echo ${pathin:="nopathin"} >/dev/null
echo ${pathout:="nopathout"} >/dev/null

if [[ ${pathin} == "nopathin" ]]; then
    echo "=== please set directory of inputs!! === ${manual}" >&2
    exit 1
elif [[ ! -d ${pathin} ]]; then
    echo "=== input directory does not exist!! === ${manual}" >&2
    exit 1
fi

if [[ ${pathout} == "nopathout" ]]; then
    echo "=== please set directory of outputs!! === ${manual}" >&2
    exit 1
elif [[ ! -d ${pathout} ]]; then
    echo "=== output directory does not exist!! === ${manual}" >&2
    exit 1
fi

(( ${threads} )) 2>/dev/null
if [[ $? != 0 || ! ${threads} -ge 1 ]]; then
    echo "=== thread number should be a positive integer!! === ${manual}" >&2
    exit 1
fi
echo ""
echo "=== set ${threads} threads for each sample in parallel analysis ==="

if [[ ${genome} == "hg19" ]]; then
    echo "=== set reference genome: hg19 ==="
    genome_Bowtie2_index="/home1/data/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
    blacklist_bed="/home1/data/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/ENCODE_blacklist_hg19_sorted.bed"
    chromsizes="/home1/data/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/hg19.chrom.sizes"
    gsize="hs"
elif [[ ${genome} == "mm10" ]]; then
    echo "=== set reference genome: mm10 ==="
    genome_Bowtie2_index="/home1/data/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
    blacklist_bed="/home1/data/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/ENCODE_blacklist_mm10_sorted.bed"
    chromsizes="/home1/data/iGenomes/Mus_musculus/UCSC/mm10/Annotation/Genes/mm10.chrom.sizes"
    gsize="mm"
else
    echo "=== supported reference genomes are hg19 or mm10 === ${manual}" >&2
    exit 1
fi

fastq_files=(`find ${pathin} -maxdepth 1 -name "*.fastq*" -type f | sort`)
if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "=== there is no fastq or fastq.gz file in input directory!! === ${manual}" >&2
    exit 1
elif (( ${#fastq_files[@]}%2 == 1 )); then
    echo "=== the number of fastq files is not even, not a paired-end data set!! === ${manual}" >&2
    exit 1
else
    echo "=== there are ${#fastq_files[@]} raw fastq files in input directory ==="
fi

echo "=== quality control with FastQC ==="
mkdir ${pathout}/fastqc
for file in ${fastq_files[@]}; do fastqc $file -t ${threads} -o ${pathout}/fastqc/ & done
wait
fastqc_files=(`ls ${pathout}/fastqc/*fastqc.zip`)
if [[ ${#fastq_files[@]} -ne ${#fastqc_files[@]} ]]; then
    echo "=== some samples not pass fastqc!!===" >&2
    exit 1
fi
multiqc -o ${pathout}/qc/ ${pathout}/fastqc/
mv ${pathout}/fastqc ${pathout}/qc/

echo "=== parallel analyzing paired-end ATAC-seq data, see log.all.txt in each subfolder ==="
mkdir ${pathout}/process
for ((i=0; i<${#fastq_files[@]}; i+=2)); do
    subfolder=$(echo ${fastq_files[$i]} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1)
    mkdir ${pathout}/process/${subfolder}
    echo "=== analyzing for sample ${subfolder} ==="
    {
        echo "=== trimming data and QC for sample ${subfolder} ==="
        trim_galore ${fastq_files[$i]} ${fastq_files[$i+1]} --suppress_warn --paired --gzip --o ${pathout}/process/${subfolder}/
        echo "=== screening contamination for sample ${subfolder} ==="
        fastq_screen ${pathout}/process/${subfolder}/*.fq.gz --subset 1000000 --threads ${threads} --quiet
        echo "=== mapping to reference genome for sample ${subfolder} ==="
        bowarr=(${pathout}/process/${subfolder}/*.fq.gz)
        bowtie2 -x ${genome_Bowtie2_index} -X 650 -I 20 -p ${threads} -t --very-sensitive -1 ${bowarr[0]} -2 ${bowarr[1]} -S ${pathout}/process/${subfolder}/bowtie2align.sam 2> ${pathout}/process/${subfolder}/bowtie2align.log
        echo "=== remove low MAPQ, chrM, sort and remove duplicates for sample ${subfolder} ==="
        perl -lane 'print $_ if $F[2] ne "chrM"' ${pathout}/process/${subfolder}/bowtie2align.sam | samtools view -bSq 20 - > ${pathout}/process/${subfolder}/bowtie2align.nochrM.bam
        rm ${pathout}/process/${subfolder}/bowtie2align.sam
        samtools sort -o ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.bam ${pathout}/process/${subfolder}/bowtie2align.nochrM.bam
        samtools flagstat ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.bam > ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.flagstat
        java -Xmx10g -XX:ParallelGCThreads=${threads} -jar $PICARD MarkDuplicates I=${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.bam O=${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.bam M=${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.metrics REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/home1/wangt5/tmp/
        samtools flagstat ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.bam > ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.flagstat
        echo "=== remove ENCODE blacklist regions for sample ${subfolder} ==="
        bedtools intersect -v -abam ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.bam -b ${blacklist_bed} > ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.bam
        samtools flagstat ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.bam > ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.flagstat
        samtools index ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.bam
        echo "=== call peaks with macs2 BAMPE mode for sample ${subfolder} ==="
        samtools sort -n -o ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.bam ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.bam
        macs2 callpeak -g ${gsize} -f BAMPE -t ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.bam -n ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2 -B --trackline --SPMR 2> ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2.log
        echo "=== generate bigwig file for sample ${subfolder} ==="
        sed -i '1d' ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_peaks.narrowPeak
        head -n 1 ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_treat_pileup.bdg > ${pathout}/process/${subfolder}/tmp.header
        sed -i '1d' ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_treat_pileup.bdg
        sort -k1,1 -k2,2n ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_treat_pileup.bdg > ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_treat_pileup.bdg.sorted
        cat ${pathout}/process/${subfolder}/tmp.header ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_treat_pileup.bdg.sorted > ${pathout}/process/${subfolder}/tmp.txt
        mv ${pathout}/process/${subfolder}/tmp.txt ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_treat_pileup.bdg.sorted
        bedGraphToBigWig ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_treat_pileup.bdg.sorted ${chromsizes} ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_treat_pileup.bw
        rm ${pathout}/process/${subfolder}/tmp.header
        rm ${pathout}/process/${subfolder}/bowtie2align.nochrM.sorted.rmdup.rmblacklist.sortname.macs2_treat_pileup.bdg.sorted
    } > ${pathout}/process/${subfolder}/log.all.txt 2>&1 &
    done
wait

# Generate summary table
echo "=== generate summary table ==="
#summarize_ATACseq.pl is added in /usr/local/bin/
summarize_ATACseq.pl ${pathout}/process/ ${pathout}/process/summary.txt

echo "=== Finished ==="

