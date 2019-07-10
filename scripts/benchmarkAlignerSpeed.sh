#####################################
#   AUTHOR: ALEXANDER LACHMANN      #
#   DATE: 5/16/2019                 #
#   Mount Sinai School of Medicine  #
#####################################

# This script will align an SRA file and create a transcript and gene quantification file
# some aligners have built in mechanisms to create gene counts, others require additional tools to create counts from aligned reads

thread_num=16

center() {
  termwidth="$(tput cols)"
  padding="$(printf '%0.1s' ={1..500})"
  printf '%*.*s %s %*.*s\n' 0 "$(((termwidth-2-${#1})/2))" "$padding" "$1" 0 "$(((termwidth-1-${#1})/2))" "$padding"
}

#small helper function for downloads
downloadSRA(){
    center "Download SRA file"
    rm -f sradata/*
    fasterq-dump \
        --split-files \
        --outdir $2 \
        -e $thread_num -p -t "sradata" \
        $1
}

# align list or SRA files with a defined species
alignFile(){
    input=$1
    while IFS= read -r line
    do
        SRA_FILE=$line
        
        printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' =
        center $SRA_FILE
        downloadSRA $SRA_FILE "sradata"
        count=`find sradata -name "${SRA_FILE}*" | wc -l`
        
        if [ "$count" -eq "1" ];
        then
            alignSingleSRA $SRA_FILE $2
        else
            alignPairedSRA $SRA_FILE $2
        fi
        
        rm sradata/${SRA_FILE}*
        
        rm -f quant/star/${SRA_FILE}/${SRA_FILE}*.bam
        rm -f quant/star/${SRA_FILE}/${SRA_FILE}*.mate1
        rm -f quant/star/${SRA_FILE}/${SRA_FILE}*.mate2
        rm -f quant/star/${SRA_FILE}/${SRA_FILE}*J.out.tab
        rm -f quant/hisat2/${SRA_FILE}/${SRA_FILE}*.sam
        rm -f quant/bwa-mem2/${SRA_FILE}/${SRA_FILE}*.sam
        
    done < "$input"
}

# this function will align reads for single read data
# fastq-dump will create one file for the corresponding reads
alignSingleSRA(){
    SRA_FILE=$1
    SPECIES=$2
    
    center "Salmon"
    { time {
        echo "-----------" >> time/${SRA_FILE}_time.txt
        echo salmon >> time/${SRA_FILE}_time.txt
        salmon quant -i "index/salmon/salmon_${SPECIES}_96" -l A \
            -r "sradata/${SRA_FILE}.fastq" \
            -p $thread_num -q --validateMappings -o quant/salmon/$SRA_FILE
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    center "kallisto"
    { time {
        echo "-----------" >> time/${SRA_FILE}_time.txt
        echo kallisto >> time/${SRA_FILE}_time.txt
        kallisto quant -i index/kallisto/kallisto_${SPECIES}_96.idx \
            -t $thread_num -o quant/kallisto/$SRA_FILE \
            --single -l 180 -s 20 "sradata/${SRA_FILE}.fastq"
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    center "STAR"
    mkdir -p quant/star/$SRA_FILE
    { time {
        echo "-----------" >> time/${SRA_FILE}_time.txt
        echo star >> time/${SRA_FILE}_time.txt
        STAR \
            --genomeDir "index/star/${SPECIES}_96" \
            --limitBAMsortRAM 10000000000 \
            --runThreadN $thread_num \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFileNamePrefix quant/star/$SRA_FILE/$SRA_FILE \
            --readFilesIn "sradata/${SRA_FILE}.fastq" \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outSAMmode Full \
            --quantMode GeneCounts \
            --limitIObufferSize 50000000
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    center "HISAT2"
    mkdir -p quant/hisat2/$SRA_FILE
    { time {
        echo "-----------" >> time/${SRA_FILE}_time.txt
        echo hisat2 >> time/${SRA_FILE}_time.txt
        hisat2 \
            -x "index/hisat2/${SPECIES}_96/${SPECIES}" \
            -U "sradata/${SRA_FILE}.fastq" \
            -p $thread_num \
            -S "quant/hisat2/${SRA_FILE}/${SRA_FILE}.sam"
        featureCounts -T $thread_num -a "reference/${SPECIES}_96.gtf" -o "quant/hisat2/${SRA_FILE}/${SRA_FILE}.tsv" "quant/hisat2/${SRA_FILE}/${SRA_FILE}.sam"
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    mkdir -p quant/bwa/$SRA_FILE
    { time {
        bwa mem \
            -t $thread_num \
            "index/bwa/${SPECIES}_96/${SPECIES}_96" \
            "sradata/${SRA_FILE}_1.fastq" \
            "sradata/${SRA_FILE}_2.fastq" \
            > "quant/bwa/${SRA_FILE}/${SRA_FILE}.sam"
        featureCounts -T $thread_num -a "reference/${SPECIES}_96.gtf" -o "quant/bwa/${SRA_FILE}/${SRA_FILE}.tsv" "quant/bwa/${SRA_FILE}/${SRA_FILE}.sam"
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    Rscript --vanilla scripts/aggregatecounts.r $SRA_FILE $SPECIES
    aws s3 cp "quant/combined/${SRA_FILE}.tsv" "s3://mssm-genecount-combined/${SRA_FILE}.tsv"
}

# this function will align reads for paired end data
# fastq-dump will create two files for the corresponding reads
alignPairedSRA(){
    SRA_FILE=$1
    SPECIES=$2
    
    center "Salmon"
    
    { time {
        echo "-----------" >> time/${SRA_FILE}_time.txt
        echo salmon >> time/${SRA_FILE}_time.txt
        salmon quant -i "index/salmon/salmon_${SPECIES}_96" -l A \
            -1 "sradata/${SRA_FILE}_1.fastq" \
            -2 "sradata/${SRA_FILE}_2.fastq" \
            -p $thread_num -q --validateMappings -o quant/salmon/$SRA_FILE
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    center "kallisto"
    { time {
        echo "-----------" >> time/${SRA_FILE}_time.txt
        echo kallisto >> time/${SRA_FILE}_time.txt
        kallisto quant -i index/kallisto/kallisto_${SPECIES}_96.idx \
            "sradata/${SRA_FILE}_1.fastq" \
            "sradata/${SRA_FILE}_2.fastq" \
            -t $thread_num -o quant/kallisto/$SRA_FILE
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    center "STAR"
    mkdir -p quant/star/$SRA_FILE
    
    { time {
        echo "-----------" >> time/${SRA_FILE}_time.txt
        echo star >> time/${SRA_FILE}_time.txt
        STAR \
            --genomeDir "index/star/${SPECIES}_96" \
            --limitBAMsortRAM 10000000000 \
            --runThreadN $thread_num \
            --outSAMstrandField intronMotif \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFileNamePrefix quant/star/$SRA_FILE/$SRA_FILE \
            --readFilesIn "sradata/${SRA_FILE}_1.fastq" "sradata/${SRA_FILE}_2.fastq" \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outSAMmode Full \
            --quantMode GeneCounts \
            --limitIObufferSize 50000000
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    center "HISAT2"
    mkdir -p quant/hisat2/$SRA_FILE
    { time {
        echo "-----------" >> time/${SRA_FILE}_time.txt
        echo hisat2 >> time/${SRA_FILE}_time.txt
        hisat2 \
            -x "index/hisat2/${SPECIES}_96/${SPECIES}" \
            -1 "sradata/${SRA_FILE}_1.fastq" \
            -2 "sradata/${SRA_FILE}_2.fastq" \
            -p $thread_num \
            -S "quant/hisat2/${SRA_FILE}/${SRA_FILE}.sam"
        featureCounts -T $thread_num -a "reference/${SPECIES}_96.gtf" -o "quant/hisat2/${SRA_FILE}/${SRA_FILE}.tsv" "quant/hisat2/${SRA_FILE}/${SRA_FILE}.sam"
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    center "BWA"
    mkdir -p quant/bwa/$SRA_FILE
    { time {
        echo "-----------" >> time/${SRA_FILE}_time.txt
        echo bwa >> time/${SRA_FILE}_time.txt
        bwa mem \
            -t $thread_num \
            "index/bwa/${SPECIES}_96/${SPECIES}_96" \
            "sradata/${SRA_FILE}_1.fastq" \
            "sradata/${SRA_FILE}_2.fastq" \
            > "quant/bwa/${SRA_FILE}/${SRA_FILE}.sam"
        featureCounts -T $thread_num -a "reference/${SPECIES}_96.gtf" -o "quant/bwa/${SRA_FILE}/${SRA_FILE}.tsv" "quant/bwa/${SRA_FILE}/${SRA_FILE}.sam"
    } 2> temp.txt ; } 2>> time/${SRA_FILE}_time.txt
    
    Rscript --vanilla scripts/aggregatecounts.r $SRA_FILE $SPECIES
    pigz "quant/combined/${SRA_FILE}.tsv"
    aws s3 cp "quant/combined/${SRA_FILE}.tsv.gz" "s3://mssm-genecount-combined/${SRA_FILE}.tsv.gz"
}


testtime(){
    { time {
        echo "-----------" >> time/${1}_time.txt
        echo bwa >> time/${1}_time.txt
        echo "haha"
    } 2> temp.txt ; } 2>> time/${1}_time.txt
    
    { time {
        echo "-----------" >> time/${1}_time.txt
        echo bwa >> time.txt
        echo "haha"
    } 2> temp.txt ; } 2>> time/${1}_time.txt
    
    { time {
        echo "-----------" >> time/${1}_time.txt
        echo bwa >> time.txt
        echo "haha"
    } 2> temp.txt ; } 2>> time/${1}_time.txt
    
    { time {
        echo "-----------" >> time/${1}_time.txt
        echo bwa >> time/${1}_time.txt
        echo "haha"
    } 2> temp.txt ; } 2>> time/${1}_time.txt
}

# First we need to download a file
# this can be done using sra-tools, there is a couple of functions that perform this task
# fastq-dump is the standard way, it will however cache the SRA file in a folder which will create an issue if yoiu are downloading many files
# fasterq-dump generally performs much faster but causes issues with HISAT2 which can't read the generated output

# create folder in which we want to move the sequencing data
mkdir -p sradata
mkdir -p quant/combined
mkdir -p quant/salmon/
mkdir -p quant/kallisto/
mkdir -p quant/star/
mkdir -p quant/hisat2/

# when using high thread count in STAR there will be a large number of open files. This can lead to errors
# set ulimit -n 1024 and it should work with many threads
ulimit -n 1024


# before we start downloading data we should disable the caching of SRA files
# when downloading many files this will quickly end up filling up all disk space otherwise
# this command below magically works to disable caching
mkdir -p ~/.ncbi
echo '/repository/user/cache-disabled = "true"' > ~/.ncbi/user-settings.mkfg


SRA_FILE="SRR2972202"
time downloadSRA $SRA_FILE "sradata"
time alignPairedSRA $SRA_FILE human

rm -f time/${SRA_FILE}_time.txt
echo "---- file ${SRA_FILE} -------">> time/${SRA_FILE}_time.txt
for i in 1 2 3 4 5 6 8 12 16
do
    thread_num=$i
    echo "---- thread count $i -------">> time/${SRA_FILE}_time.txt
    time alignPairedSRA $SRA_FILE human
done


