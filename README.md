HARSHMI SHAH


Download 
- VirtualBox :
``` {sh} 
https://www.virtualbox.org/wiki/Downloads 
```
- Basespace Native App VM 
``` {sh} 
https://da1s119xsxmu0.cloudfront.net/sites/developer/native/nativeappsvm/BaseSpace%20Native%20App%20VM%20(phix%20only)%20v9.ova 
```


Import Basespace Native App VM in Virtual Box
- Open File tab -> Import Appliance -> Browse and select Basespace's .ova file
- Start


Connect via PuTTY
- Username : basespace@localhost
- Port : 2222
- Password : basespace


Requirements:
- Python3 v3.4.1
- GATK v3.4
- Picard v2.6.0 
- Bowtie v2.2.9
- Trimmomatic v0.36)


Clone AHCG Pipeline from GitHub and get all requirements
``` {sh}
git clone https://github.com/shashidhar22/ahcg_pipeline.git
git pull origin master 
```


Installations:
- Samtools 
``` {sh} 
sudo apt-get install samtools	
```
- Java
``` {sh}
	sudo apt-get install software-properties-common python-software-properties
	sudo add-apt-repository ppa:webupd8team/java
	sudo apt-get update
	sudo apt-get install oracle-java8-installer
	java -version
```


Downloads:
- Reference genome and dbsnp: 
Download: 
``` {sh}
wget  www.prism.gatech.edu/~sravishankar9/resources.tar.gz
```
Extract: 
``` {sh}
tar -xvzf resources.tar.gz
```


- Test data:

	- Download: 
``` {sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
```
	- Extract: 
```{sh}
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
```
	- Build a test set: 
``` {sh}
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```


Getting index files
- Bowtie :
``` {sh}
./lib/bowtie2-2.2.9/bowtie2-build -f ./resources/genome/hg19.fa hg19
```
- Fasta index using Samtools : 
``` {sh}
samtools faidx ./resources/genome/hg19.fa
```
- Genome dict file using picard : 
``` {sh}
java -jar ./lib/picard.jar CreateSequenceDictionary R=./resources/genome/hg19.fa O=hg19.dict 
```


Run the script:
- Help : 
``` {sh}
python3 ahcg_pipeline.py -h
```
- Command :
```{sh}
python3 ahcg_pipeline.py 
-t ./lib/Trimmomatic-0.36/trimmomatic-0.36.jar
		Path to Trimmomatic 
-b ./lib/bowtie2-2.2.9/bowtie2 
		Path to Bowtie
-p ./lib/picard.jar 
		Path to Picard
-g ./lib/GenomeAnalysisTK.jar
		Path to GATK 
-i /path/to/folder/test*.fastq
		Path to read files 
-w /path/to/folder/hg19 
		Path to Bowtie indexes
-d ./resources/dbsnp/dbsnp_138.hg19.vcf
		Path to dbSNP vcf 
-r ./resources/genome/hg19.fa 
		Path to reference genome
-a ./lib/Trimmomatic-0.36/adapters/TruSeq2-PE.fa
		Path to adapter sequence 
-o <outputdir>
		Path to output directory
```


Change the remote URL for GIT Repository
- Fork the original repository to personal GitHub
- Change the remote URL for GIT Repository by adding the URL of personal repository in .git/config file


GIT Ignore
- Use the .gitignore file to add files/directories that you want to ignore when updating the git repository 


Set up name & email address
``` {sh}
git config --global user.email johndoe@example.com
git config --global user.name "John Doe"
```


GIT commands to upload files/folder
``` {sh}
git add <filename or * for everything>
git commit -m "comment describing the upload"
git push origin master
```


Get gene annotation file 
``` {sh}
wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt
```


Locate BRCA1 gene and variants
```{sh}
grep BRCA1 hg19_refGene.txt
``` 
- Found 6 variants
- Find which one is best
	- Go to https://dnasu.org/DNASU/Home.do
	- Search BRCA and select reference HsCD00022337
	- Go to Reference Sequence Alignment
	- Based on this alignment NM_007294 was the best 
- Use the NM_007294 to get genome coordinates and create a bed file
	- grep NM_007294 hg19_refGene.txt > test.txt
	- use the script bed.py to convert it to bed file
- Get fasta sequences from bedtools
	- Install bedtools & get fasta :
``` {sh}
sudo apt-get install bedtools
bedtools getfasta -s -fo brcafa.fa -fi ./resources/genome/hg19.fa -bed brca1.bed
```


Extracting reads mapped to region of interest
- Download four bam files of  NA12878 exome from GIAB ftp website
``` {sh}
ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/
```
- Use Samtools to get the region of interest i.e. BRCA1
``` {sh}
samtools view -L <bed file> -b -o <output bam file> <input bam file>
```
- Convert bam to fastq
``` {sh}
bedtools bamtofastq -i <bam file> -fq <fastq r1> -fq2 <fastq r2>
```


Verification of variant call:
- Extract the lines corresponding to the genes in breastcancer_genes.txt from reference file (hg19_refGene.txt)
``` {sh}
awk '{print "\\<" $2 "\\>" }' breastcancer_genes.txt > genelist.txt
grep -f genelist.txt hg19_refGene.txt > bed.txt
```
- Use bed.py to create bed file using the genes in bed.txt
``` {sh}
./bed.py -i bed.txt -o outputfile.bed
```
- Compare the variants.vcf file with the outputfile.bed 
``` {sh}
bedtools intersect -header -wa -a variants.vcf -b outputfile.bed > outfile1
```
- Compare the GIAB vcf file with outputfile.bed
``` {sh}
add "chr" to get correct match
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf > giab_with_chr.vcf
bedtools intersect -header -wa -a giab_with_chr.vcf -b outputfile.bed > outfile2
```
- Compare the intersected files
``` {sh}
bedtools intersect -header -a outputfile1 -b outfile2 > outfile
```
