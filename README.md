# Alignment Benchmark

This is a benchmark tool to test the performance of genomic alignment algorithms. The tool is encapsulated in a docker container and all functionality can be opened in Jupyter notebooks. 
The Jupyter notebooks are either shell scripts or R scipts. There are currently 5 different aligners that can be tested. These algorithms are BWA, STAR, HISAT2, kallisto, and Salmon. There are currently 4 notebooks for different types of benchmarks:

- Build index files: This Jupyter notebook will build index files required to run the alignment. There is also the option of downloading prebuild index files for human and mouse reference genomes

- Align SRA files listed in a text file: This Jupyter Notebook will download and align all SRA files form a file

- Alignment Speed Benchmark: This Jupyter notebook will run alignment of a specified SRA file on the 5 alignment algorithms and writes the elaped time to a file

- Process gene expression: This is a Jupyter Notebook with an R script that aggregates gene count data produced by the 5 aligners and wirtes them into text files for further processing

## Running Jupyter Notebooks

To run the Jupyter Notebooks in the maayanlab/alignmentbenchmark docker is required. To run the docker image locally in a unix environment use the following command:

docker run -p 8889:8888 maayanlab/alignmentbenchmark

This will launch a virtual machine with a Jupyter Notebook webserver exposed at the port 8889 of the host machine. To access the notbooks the URL is: http://localhost:8889/tree.

The resources required to perform alignment varies with the alignment algorithms used. STAR needs about 20GB of memory. BWA requires more than 40GB of memory to build the index file. Because of this we recommend running the Jupyter notebook in a cloud compute environment such as Amazon AWS. A good instance to use is c5.4xlarge from AWS. It has 32GB of memory and 16 vCPUs. The current cost of such an image is $0.68. This instance does however not provide enough memory to build the BWA index. Consider downloading the prebuild indeces in this case.

## SRA data

The support file folder contains a list of human and mouse SRA ids. It also contains a signatures.tsv file with information about gene expression studies from which differential gene expression can be copmputed.

## Adding additional alignment tools

To preconfigure the docker image with all dependencies we rely on conda to install software. In the DockerAlign dockerfile additional aligners can be installed. For this add commands such as RUN conda install -y star. The Jupyter Notebook files will also have to be modified accordingly to support additional algorithms.
