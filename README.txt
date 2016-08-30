AHCG PIPELINE README 
BY HARSHMI SHAH

Download 
- VirtualBox : https://www.virtualbox.org/wiki/Downloads
- Basespace Native App VM : https://da1s119xsxmu0.cloudfront.net/sites/developer/native/nativeappsvm/BaseSpace%20Native%20App%20VM%20(phix%20only)%20v9.ova

Import Basespace Native App VM in Virtual Box

Connect via PuTTY
- Username : basespace@localhost
- Port : 2222
- Password : basespace

Clone ahcg_pipeline from GitHub
- git clone https://github.com/shashidhar22/ahcg_pipeline.git
- Open the folder and download lib files : git pull origin master (for Python3, GATK, Picard, Bowtie & Trimmomatic)

Installations:
- Samtools : sudo apt-get install samtools	
- Java
	sudo apt-get install software-properties-common python-software-properties
	sudo add-apt-repository ppa:webupd8team/java
	sudo apt-get update
	sudo apt-get install oracle-java8-installer
	java -version

Downloads:
- Ref genome and dbsnp: www.prism.gatech.edu/~sravishankar9/resources.tar.gz
- Test data:
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq

Getting index files
- Bowtie :  ./lib/bowtie2-2.2.9/bowtie2-build -f ./resources/genome/hg19.fa hg19
- Fasta index using Samtools : samtools faidx ./resources/genome/hg19.fa
- Genome dict file using picard : java -jar ./lib/picard.jar CreateSequenceDictionary R=./resources/genome/hg19.fa O=hg19.dict 

Run the script:
python3 ahcg_pipeline.py -t ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b ./lib/bowtie2-2.2.9/bowtie2 -p ./lib/picard.jar -g ./lib/GenomeAnalysisTK.jar -i /path/to/folder/test*.fastq -w /path/to/folder/hg19 -d ./resources/dbsnp/dbsnp_138.hg19.vcf -r ./resources/genome/hg19.fa -a ./lib/Trimmomatic-0.36/adapters/TruSeq2-PE.fa -o <outputdir>


### GIT commands
git add <file>
git commit -m "comment"
git push origin master
