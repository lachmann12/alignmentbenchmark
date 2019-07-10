#####################################
#   AUTHOR: ALEXANDER LACHMANN      #
#   DATE: 5/16/2019                 #
#   Mount Sinai School of Medicine  #
#####################################

SPECIES=$1
AWS_BUCKET=$2

downloadIndex(){
    mkdir -p index/$1
    curl https://s3.amazonaws.com/$3/$2 -o index/$1/$2
    tar -I pigz -xvf index/$1/$2
    rm index/$1/$2
}

# this script will download all files needed to align reads. Alternatively the script create_index can be run
# this would built indeces from scratch. The indeces for this project are all publically available at

# 1) Download all data required for alignment or transcript quantification
#    We are using the Ensembl genomes and annotations in this project

mkdir -p reference
mkdir -p index/salmon
mkdir -p index/kallisto
mkdir -p index/star/${SPECIES}_96
mkdir -p index/hisat2/${SPECIES}_96
mkdir -p index/bwa/${SPECIES}_96

# fetch GTF files annotationg the raw DNA sequences
# they are needed by the alignment algorithms STAR and HISAT2 to count transcript reads
if [ $SPECIES = "human" ]; then
    url="ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz"
else
    url="ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz"
fi

curl $url -o reference/${SPECIES}_96.gtf.gz &

# the script is now waiting for the downloads to finish before continuing
# the jobs -p command lists forked jobs (commands followed by &)
for job in `jobs -p`
do
echo $job
    wait $job || let "FAIL+=1"
done

# unpack some files that can't be used when zipped
pigz -d reference/${SPECIES}_96.gtf.gz

# this will download a large amout of data and will take some time, this will still be much faster than building the 
# indeces from scratch
downloadIndex "salmon" "salmon_${SPECIES}_96.tar.gz" ${AWS_BUCKET}
downloadIndex "kallisto" "kallisto_${SPECIES}_96.idx.tar.gz" ${AWS_BUCKET}
downloadIndex "star" "star_${SPECIES}_96.tar.gz" ${AWS_BUCKET}
downloadIndex "hisat2" "hisat2_${SPECIES}_96.tar.gz" ${AWS_BUCKET}
downloadIndex "bwa" "bwa_${SPECIES}_96.tar.gz" ${AWS_BUCKET}






