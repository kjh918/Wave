

#! /bin/bash

# Author  : kangsm
# Version/ModifyDate History
# v1.0.0001 / 20250905
# v1.1.0001 / 20250930

### COMMAND OPTIONS --------------------------------------------------------------------------------
    usage() { 
        echo "Usage: cbNIPT_run.NGS_Processing.sh [ OPTIONS ]...
        --SeqID            = SeqID. REQUIRED. 
        --SampleID         = SampleID. if NULL, SeqID will be used as SampleID
        --Threads          = ncpu cores. default = 10
        --PicoplexGold     = WGA kit is Picoplex-gold ot not. default = "No"
        --ReadGroupID      = RGID for BAM header. use Flowcell-ID.
        --RunFastqc        = run FastQC or not. 'Yes' or 'No'
        --RunFastqTrimming = run fastp fastq-trimming. 'Yes' or 'No'
        --RunUnalignedBAM  = run uBAM pipeline or not. 'Yes' or 'No'
        --RunLocalReAlign  = run local re-alignment or not. 'Yes' or 'No'
        -h | --help        = print this usage
        "
    }
###
### PARSE CMD OPTIONS ------------------------------------------------------------------------------
    ARGS=$(getopt -a -o h: --long SeqID:,SampleID:,Threads:,PicoplexGold:,ReadGroupID:,RunFastqc:,RunFastqTrimming:,RunUnalignedBAM:,RunLocalReAlign:,help -- "$@" )
    VALID_ARGS=$?
    if [ "$VALID_ARGS" != "0" ]; then 
        usage >&2 
        exit 2
    fi

    eval set -- "$ARGS"
    while :
    do
    case "$1" in
        --SeqID )
            case "$2" in
                -* | --* | "") echo "no value at $1. REQUIRED." >&2 ; exit 2 ;;
                *)  SeqID=$2 ; shift 2 ;;
            esac ;;
        --SampleID )
            case "$2" in
                -* | --* | "") echo "no value at $1. default = SeqID " ; shift ;;
                *)  SampleID=$2 ; shift 2 ;;
            esac ;;
        --Threads )
            case "$2" in
                -* | --* | "") echo "no value at $1. default = 10 " ; shift ;;
                *)  Threads=$2 ; shift 2 ;;
            esac  ;;
        --PicoplexGold )
            case "$2" in
                -* | --* | "") echo "no value at $1. default = No " ; shift ;;
                *)  PicoplexGold=$2 ; shift 2 ;;
            esac ;;
        --ReadGroupID )
            case "$2" in
                -* | --* | "") echo "no value at $1. REQUIRED." >&2 ; exit 2 ;;
                *)  ReadGroupID=$2 ; shift 2 ;;
            esac ;;
        --RunFastqc )
            case "$2" in
                -* | --* | "") echo "no value at $1. default = No " ; shift ;;
                *)  RunFastqc=$2 ; shift 2 ;;
            esac ;;
        --RunFastqTrimming )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = No" ; shift ;;
                *)  RunFastqTrimming=$2 ; shift 2 ;;
            esac ;;
        --RunUnalignedBAM )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = No" ; shift ;;
                *)  RunUnalignedBAM=$2 ; shift 2 ;;
            esac ;;
        --RunLocalReAlign )
            case "$2" in
                -* | --* | "") echo "no value at $1, default = No" ; shift ;;
                *)  RunLocalReAlign=$2 ; shift 2 ;;
            esac ;;
        -h | --help )
        usage >&2 ; exit 2 ;;   
        --) shift ; break ;;
        *)  usage >&2 ; exit 2 ;;
    esac
    done
###
### debug params ###
    # SeqID="cbNIPT_PubDB_01_02"
    # SampleID="cbnipt.pubdb.01_trisomy7.sample02.wgs"
    # Threads=10
    # PicoplexGold="No"
    # ReadGroupID="SRR27376808"
    # RunFastqc="No"
    # RunFastqTrimming="No"
    # RunUnalignedBAM="No"
    # RunLocalReAlign="No"
###
### DEFAULT PARAMS ---------------------------------------------------------------------------------
    if [ -z "$Threads" ]; then Threads=10; fi
    if [ -z "$SampleID" ]; then SampleID=${SeqID}; fi
    if [ -z "$PicoplexGold" ]; then PicoplexGold="No"; fi
    if [ -z "$RunFastqc" ]; then RunFastqc="No"; fi
    if [ -z "$RunFastqTrimming" ]; then RunFastqTrimming="No"; fi
    if [ -z "$RunUnalignedBAM" ]; then RunUnalignedBAM="Yes"; fi
    if [ -z "$RunLocalReAlign" ]; then RunLocalReAlign="Yes"; fi
    #---------------------------------------------------------------------------
    BaseDir=/data/cbNIPT
    BatchID=processed_dataset
    ReadGroupPlatform=ILLUMINA
    ReadGroupLibrary=WGS
    ReadGroupCenter=GENCURIX
    RawFastqDir=/data/cbNIPT/fastq_raw
    TrimFastqDir=/data/cbNIPT/fastq_trimmed
    BamToBedDir=/data/cbNIPT/bamToBeds
###
### SAMPLE-DEPEND PARAMS ---------------------------------------------------------------------------
    # BAM DIR
    BamDir=${BaseDir}/${BatchID}/${SeqID}/bam
    if [ ${RunUnalignedBAM} = "No" ]; then
        BamDir=${BamDir}_bwa_only
    fi
    if [ ! -d ${BamDir} ]; then mkdir -p ${BamDir}; fi
    # QC RESULT DIR
    qcResDir=${BaseDir}/${BatchID}/${SeqID}/qcfiles
    if [ ! -d ${qcResDir} ]; then 
        mkdir -p ${qcResDir}
    fi
    # TEMP DIR
    TmpDir=${BamDir}/tmp
    if [ ! -d ${TmpDir} ]; then mkdir -p ${TmpDir}; fi
###
### RESOURCES --------------------------------------------------------------------------------------
    BwaIndex=/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa
    ReferenceFasta=/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa
    KnownSnp=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz
    KnownIndel1=/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz
    KnownIndel2=/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
###
### RUN FASTQC -------------------------------------------------------------------------------------
    if [ ${RunFastqc} = "Yes" ]; then
        singularity exec -B /storage,/data /storage/images/fastqc-0.12.1.sif fastqc --extract --threads 8 --outdir ${qcResDir} \
            ${RawFastqDir}/${SeqID}_R1.fastq.gz ${RawFastqDir}/${SeqID}_R2.fastq.gz
    fi
###
### RUN FASTQ-TRIMMING -----------------------------------------------------------------------------
    if [ ${RunFastqTrimming} = "Yes" ]; then
        if [ ${PicoplexGold} = "Yes" ]; then
            singularity exec -B /storage,/data /storage/images/fastp-0.23.4.sif fastp --thread ${Threads} \
                --in1 ${RawFastqDir}/${SeqID}_R1.fastq.gz --in2 ${RawFastqDir}/${SeqID}_R2.fastq.gz \
                --out1 ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz --out2 ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz \
                --json ${qcResDir}/${SeqID}.fastp.json --html ${qcResDir}/${SeqID}.fastp.html \
                --trim_poly_g --detect_adapter_for_pe --adapter_sequence "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" --adapter_sequence_r2 "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
                --length_required 100 --average_qual 10 --qualified_quality_phred 15 \
                --trim_front1 14 --trim_front2 14 
        else
            singularity exec -B /storage,/data /storage/images/fastp-0.23.4.sif fastp --thread ${Threads} \
                --in1 ${RawFastqDir}/${SeqID}_R1.fastq.gz --in2 ${RawFastqDir}/${SeqID}_R2.fastq.gz \
                --out1 ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz --out2 ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz \
                --json ${qcResDir}/${SeqID}.fastp.json --html ${qcResDir}/${SeqID}.fastp.html \
                --trim_poly_g --detect_adapter_for_pe --adapter_sequence "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" --adapter_sequence_r2 "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
                --length_required 100 --average_qual 10 --qualified_quality_phred 15 
        fi
    fi
###
### ALIGN TO REFERECE GENOME -----------------------------------------------------------------------
    if [ ${RunUnalignedBAM} = "Yes" ]; then
        touch ${BamDir}/uBAM_and_BWA_Pipeline
        # unaligned-bam
        java -XX:ParallelGCThreads=14 -Xmx16384m -jar /storage/apps/bin/picard.jar FastqToSam \
            --FASTQ ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz --FASTQ2 ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz \
            --SAMPLE_NAME ${SeqID} --OUTPUT ${BamDir}/${SeqID}.fastqtosam.bam \
            --READ_GROUP_NAME ${ReadGroupID} --PLATFORM ${ReadGroupPlatform} --LIBRARY_NAME ${ReadGroupLibrary} --SEQUENCING_CENTER ${ReadGroupCenter} \
            --TMP_DIR ${TmpDir}
        # aligned-bam
        singularity exec -B /storage,/data /storage/images/bwa-0.7.17.sif bwa mem -M -t ${Threads} -Y -L 50,50 \
            -R "@RG\tID:${ReadGroupID}\tPL:${ReadGroupPlatform}\tLB:${ReadGroupLibrary}\tSM:${SeqID}\tCN:${ReadGroupCenter}" ${BwaIndex} \
            ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz > ${BamDir}/${SeqID}.bwa.mem.sam 
        # merge bam alignment
        java -XX:ParallelGCThreads=14 -Xmx16384m -jar /storage/apps/bin/picard.jar MergeBamAlignment \
            --UNMAPPED_BAM ${BamDir}/${SeqID}.fastqtosam.bam --ALIGNED_BAM ${BamDir}/${SeqID}.bwa.mem.sam \
            --REFERENCE_SEQUENCE ${ReferenceFasta} --OUTPUT ${BamDir}/${SeqID}.primary.bam \
            --CREATE_INDEX true --MAX_INSERTIONS_OR_DELETIONS -1 --CLIP_ADAPTERS false --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --ATTRIBUTES_TO_RETAIN XS --EXPECTED_ORIENTATIONS FR --EXPECTED_ORIENTATIONS RF    
    else
        touch ${BamDir}/BWA_only_Pipeline
        # run bwa
        singularity exec -B /storage,/data /storage/images/bwa-0.7.17.sif bwa mem -M -t ${Threads} -Y -L 50,50 \
            -R "@RG\tID:${ReadGroupID}\tPL:${ReadGroupPlatform}\tLB:${ReadGroupLibrary}\tSM:${SeqID}\tCN:${ReadGroupCenter}" ${BwaIndex} \
            ${TrimFastqDir}/${SeqID}.trimmed_R1.fastq.gz ${TrimFastqDir}/${SeqID}.trimmed_R2.fastq.gz | samtools view -bS -o ${BamDir}/${SeqID}.primary.bam -
    fi
###
### BAM SORT ADN INDEX -----------------------------------------------------------------------------
    samtools sort -@ ${Threads} -o ${BamDir}/${SeqID}.sorted.bam ${BamDir}/${SeqID}.primary.bam
    samtools index -b -@ ${Threads} ${BamDir}/${SeqID}.sorted.bam
###
### MARK DUPLICATES --------------------------------------------------------------------------------
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk MarkDuplicates --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --INPUT ${BamDir}/${SeqID}.sorted.bam --OUTPUT ${BamDir}/${SeqID}.sorted.dedup.bam \
        --METRICS_FILE ${qcResDir}/${SeqID}.mark.duplicates.metrics.txt --CREATE_INDEX true --REMOVE_SEQUENCING_DUPLICATES true
###
### LOCAL RE-ALIGNMENT -----------------------------------------------------------------------------
    if [ ${RunLocalReAlign} = "Yes" ]; then
        # create targets
        singularity exec -B /storage,/data /storage/images/gatk-3.8-1.sif java -Xmx16384m -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator \
            -R ${ReferenceFasta} -I ${BamDir}/${SeqID}.sorted.dedup.bam -o ${qcResDir}/${SeqID}.realign.target.intervals \
            -known ${KnownIndel1} -known ${KnownIndel2} -nt ${Threads}
        # run re-align   
        singularity exec -B /storage,/data /storage/images/gatk-3.8-1.sif java -Xmx16384m -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
            -R ${ReferenceFasta} -targetIntervals ${qcResDir}/${SeqID}.realign.target.intervals -known ${KnownIndel1} -known ${KnownIndel2} \
            -I ${BamDir}/${SeqID}.sorted.dedup.bam -o ${BamDir}/${SeqID}.sorted.dedup.realign.bam
    fi
###
### BASE RECALIBRATION -----------------------------------------------------------------------------
    if [ ${RunLocalReAlign} = "Yes" ]; then
        RecalInputBAM=${BamDir}/${SeqID}.sorted.dedup.realign.bam
        touch ${BamDir}/local_realign_run
    else
        RecalInputBAM=${BamDir}/${SeqID}.sorted.dedup.bam
        touch ${BamDir}/local_realign_NOT_run
    fi
    # create recal table
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk BaseRecalibrator --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --input ${RecalInputBAM} --reference ${ReferenceFasta} --output ${qcResDir}/${SeqID}.recal.table.txt \
        --known-sites ${KnownSnp} --known-sites ${KnownIndel1} --known-sites ${KnownIndel2} 
    # run bqsr
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk ApplyBQSR --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --input ${RecalInputBAM} --bqsr-recal-file ${qcResDir}/${SeqID}.recal.table.txt --output ${BamDir}/${SeqID}.recal.bam 
###
### Filter BAM
    # filter : mapq > 20
    #        : proper-paired only
    #        : remove unmapped and secondary alignments
    #        : softclip length < 20
    #        : mismatch < 12
    
    samtools view -b -h -q 20 -f 0x2 -F 0x100 -F 0x4 -e "sclen < 20" -e "[NM] < 12" --threads ${Threads} ${BamDir}/${SeqID}.recal.bam > ${BamDir}/${SeqID}.recal.filtered.bam
    samtools index --threads ${Threads} ${BamDir}/${SeqID}.recal.filtered.bam
###
### CREATE analysisReady BAM LINK
    ln -Tsf  ${BamDir}/${SeqID}.recal.filtered.bam ${BamDir}/${SeqID}.analysisReady.bam
    ln -Tsf  ${BamDir}/${SeqID}.recal.filtered.bam.bai ${BamDir}/${SeqID}.analysisReady.bam.bai
###
### REMOVE TEMP FILES ------------------------------------------------------------------------------
    mv ${BamDir}/${SeqID}.recal.bai ${BamDir}/${SeqID}.recal.bam.bai
    rm ${BamDir}/${SeqID}.bwa.mem.sam
    rm ${BamDir}/${SeqID}.fastqtosam.bam
    rm ${BamDir}/${SeqID}.primary.bam
    rm ${BamDir}/${SeqID}.sorted.bam
    rm ${BamDir}/${SeqID}.sorted.bam.bai
    rm ${BamDir}/${SeqID}.sorted.bai
    rm ${BamDir}/${SeqID}.sorted.dedup.bam
    rm ${BamDir}/${SeqID}.sorted.dedup.bai
    if [ -f ${BamDir}/${SeqID}.sorted.dedup.realign.bam ]; then 
        rm ${BamDir}/${SeqID}.sorted.dedup.realign.bam
        rm ${BamDir}/${SeqID}.sorted.dedup.realign.bai
    fi
    rm -r ${BamDir}/tmp
###
### BAM QC -----------------------------------------------------------------------------------------
    # sequencing artifacts
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif gatk CollectSequencingArtifactMetrics --java-options "-XX:ParallelGCThreads=14 -Xmx16384m" \
        --INPUT ${BamDir}/${SeqID}.analysisReady.bam \
        --OUTPUT ${qcResDir}/${SeqID}.artifacts.txt \
        --FILE_EXTENSION .txt --REFERENCE_SEQUENCE ${ReferenceFasta}
    # alignment summary
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif java -XX:ParallelGCThreads=14 -Xmx16384m -jar /gatk/gatk-package-4.4.0.0-local.jar CollectAlignmentSummaryMetrics \
        --INPUT ${BamDir}/${SeqID}.analysisReady.bam \
        --OUTPUT ${qcResDir}/${SeqID}.alignment.summary.metrics.txt
    # wgs mteric
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif java -XX:ParallelGCThreads=14 -Xmx16384m -jar /gatk/gatk-package-4.4.0.0-local.jar CollectWgsMetrics \
        --INPUT ${BamDir}/${SeqID}.analysisReady.bam \
        --OUTPUT ${qcResDir}/${SeqID}.collect.wgs.metrics.txt \
        --REFERENCE_SEQUENCE ${ReferenceFasta}
    # insert size metric
    export LC_ALL=en_US.UTF-8
    singularity exec -B /storage,/data /storage/images/gatk-4.4.0.0.sif java -XX:ParallelGCThreads=14 -Xmx16384m -jar /gatk/gatk-package-4.4.0.0-local.jar CollectInsertSizeMetrics \
        --INPUT ${BamDir}/${SeqID}.analysisReady.bam \
        --OUTPUT ${qcResDir}/${SeqID}.insert.size.metrics.txt \
        --Histogram_FILE ${qcResDir}/${SeqID}.insert.size.histogram.pdf \
        --REFERENCE_SEQUENCE ${ReferenceFasta}
    # mosdepth (coverage)
    singularity exec -B /storage,/data /storage/images/mosdepth-0.3.6.sif /opt/mosdepth --threads ${Threads} --no-per-base --fast-mode \
        --by 100000 --mapq 20 ${qcResDir}/${SeqID} ${BamDir}/${SeqID}.analysisReady.bam 
    # copy duplicates metric
    #cp ${BamDir}/${SeqID}.mark.duplicates.metrics.txt ${qcResDir}/
    # MULTI-QC Summarization
    singularity exec -B /storage,/data /storage/images/multiqc-1.16.sif multiqc --force \
        --filename ${SeqID}.QC.Results --outdir ${qcResDir} --config /storage/home/kangsm/runScripts/NGS_config.MultiQC_Custom.yaml --data-dir ${qcResDir}
    # resutl summary     
    Rscript /storage/home/kangsm/runScripts/cbNIPT/GCX.cbNIPT_run.NGS.Stats.Processing.R --SeqID ${SeqID} --SampleID ${SampleID} --DatabaseImport TRUE
###
### BAM TO BED BamToBedDir
    bedtools bamtobed -i ${BamDir}/${SeqID}.analysisReady.bam > ${BamDir}/${SeqID}.bam.bed
    bgzip -@ ${Threads} ${BamDir}/${SeqID}.bam.bed
    ln -Tsf ${BamDir}/${SeqID}.bam.bed.gz ${BamDir}/${SampleID}.bed.gz
    cp ${BamDir}/${SeqID}.bam.bed.gz ${BamToBedDir}/${SampleID}.bed.gz
###
    


