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


Test data:

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


VariantRecalibrator
- Download the bundle files
``` {sh}
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/hapmap_3.3.hg19.sites.vcf.idx.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.gz
ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
```
- Gunzip them all
- Bgzip them
- Run the variant recalibrator
``` {sh}
java -Xmx4g -jar ./lib/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R ./resources/genome/hg19.fa \
-input ./oct_4/variants.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./oct_4/hapmap_3.3.hg19.sites.vcf.gz \
-resource:omni,known=false,training=true,truth=false,prior=12.0 ./oct_4/1000G_omni2.5.hg19.sites.vcf.gz \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 ./oct_4/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./resources/dbsnp/dbsnp_138.hg19.vcf \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP -recalFile output.recal -tranchesFile output.tranches -rscriptFile output.plots.R
```
- Run GATK again to get vcf file
``` {sh}
java -jar ./lib/GenomeAnalysisTK.jar \ 
    -T ApplyRecalibration \ 
    -R ./path/to/reference.fa \ 
    -input ./path/to/raw_variants.vcf \ 
    -mode SNP \ 
    --ts_filter_level 99.0 \ 
    -recalFile recalibrate_SNP.recal \ 
    -tranchesFile recalibrate_SNP.tranches \ 
    -o recalibrated_snps_raw_indels.vcf
``` 

Coverage 
``` {sh}
grep 'NM_007298' bcoc_padded.bed > brca1.bed
samtools view -L brca1.bed data/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam -b > new.bam
bedtools genomecov -ibam new.bam -bga > na12878.bga.bed
bedtools intersect -loj -a brca1.bed -b na12878.bga.bed -bed > brca1.join_final.bed
awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$8,$8,$4,$10,$6)}' brca1.join_final.bed \
| sed -E -e 's/^chr//' > brca1.final.bed

bedtools intersect -a brca1.final.bed -b brca_clinical_nonbenign_xref.vcf -wo > brca_clinical_nonbenign_final.bed

awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$10,$6)}' brca1.join_final.bed > brca1.depths.bed
python cov.py brca1.depths.bed brca_depth.txt
Rscript draw_depth.R brca_depth.txt brca_depth.png
```

DCM
```{sh}
#shrink clinvar to just DCM genes
bedtools intersect -a clinvar.vcf.gz -b dcm_gene_list.bed -header > clinvar_allfrombed.vcf

#recalibrate all variants using VQSR
vcf_vqrs.sh patient2_variants.vcf
vcf_apply_recal.sh patient2_variants.vcf

#shrink variants to just DCM genes
bedtools intersect -a patient2_variants_recal.vcf -b dcm_gene_list.bed -header > patient2_dcm_final.vcf

#match variants to clinvar
bedtools intersect -b patient2_dcm_final.vcf -a clinvar_allfrombed.vcf -header > patient2_intersect_clinvar.vcf

#generate simple report on findings
python3 parse_clnsig.py -i patient2_intersect_clinvar.vcf.gz 2>&1 | tee patient2_simple_report.txt
cut -c 24- patient2_simple_report.txt
```

Master Script README
```{sh}

README: README_master_script

Usage: 
./master_script.sh <output_dir_name> <patient_bam_file> <gene_list_bed> <clinvar vcf>

Example:
./master_script.sh results ./Patient1_RG_MD_IR_BQ.bam ./dcm_gene_list.bed ./clinvar.vcf

Requirements:
1. Python3 v3.4.1
	- logging
	- vcf (pyvcf)
	- click
2. GATK v3.4
3. R v3.1.2
	- ggplot2
4. Bedtools v2.26.0
5. Samtools v0.1.19
5. Imagemagick 
	- convert

Steps:
1. Variant Calling using GATK's HaplotypeCaller
2. Recalibrating using GATK's VQSR ie VariantRecalibrator & ApplyRecalibration
3. Cross referencing variants file with clinvar file using bedtools
4. Coverage Calculation using samtools & bedtools
5. Report Generation using convert

Directory Structure:
/home/basespace/ahcg_pipeline
	- lib
		* bowtie2-2.2.9
			-bowtie2
		* GATK ie GenomeAnalysisTK.jar
	- resources
		* genome
			- hg19.fa
		*dbsnp
			- dbsnp_138.hg19.vcf
	- master_script.sh
	- clinvar.vcf
	- Patient1_RG_MD_IR_BQ.bam
	- Patient1_RG_MD_IR_BQ.bai
	- oct_4
		* hapmap_3.3.hg19.sites.vcf.gz
		* 1000G_omni2.5.hg19.sites.vcf.gz
		* 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
	- cov.py
	- dcm_gene_list.bed
	- draw_depth.R
	- parse_clnsig.py
	- output_dir_name
```

Master Script (master_script.sh)
```{sh}
#!/bin/bash

if [ "$#" -ne 4 ]
then
    echo "usage: $0 output_dir patient.bam gene_list.bed clinvar.vcf.gz"
    exit 1
fi

#### VARIANT CALLING
BASE_DIR="/home/basespace/ahcg_pipeline"
GATK=$BASE_DIR/lib/GenomeAnalysisTK.jar


mkdir $1
PATIENT=$(basename $2)
echo "$PATIENT" > $1/patientname.txt
echo "Running HaplotypeCaller"
java -jar $GATK -T HaplotypeCaller -R $BASE_DIR/resources/genome/hg19.fa -I $2 --dbsnp $BASE_DIR/resources/dbsnp/dbsnp_138.hg19.vcf -o $1/variants.vcf -nct 1 -gt_mode DISCOVERY

#### VARIANT RECALIBERATION 
VARIANTS="$1/variants.vcf"
RECAL=$(dirname $VARIANTS)/output.recal
TRANCH=$(dirname $VARIANTS)/output.tranches

echo "Running VQSR"

java -Xmx4g -jar $GATK \
-T VariantRecalibrator \
-R $BASE_DIR/resources/genome/hg19.fa \
-input $VARIANTS \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $BASE_DIR//oct_4/hapmap_3.3.hg19.sites.vcf.gz \
-resource:omni,known=false,training=true,truth=false,prior=12.0 $BASE_DIR/oct_4/1000G_omni2.5.hg19.sites.vcf.gz \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 $BASE_DIR/oct_4/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $BASE_DIR/resources/dbsnp/dbsnp_138.hg19.vcf \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP -recalFile $RECAL -tranchesFile $TRANCH -rscriptFile output.plots.R


suffix=_recal.vcf
outname="$(basename $VARIANTS | sed -e 's/\.vcf$//')$suffix"

echo "Applying VQSR"
java -jar $GATK -T ApplyRecalibration -R $BASE_DIR/resources/genome/hg19.fa -input $VARIANTS -mode SNP --ts_filter_level 99.0 -recalFile $RECAL -tranchesFile $TRANCH -o $1/$outname


#### CROSS REFERENCING VARIANT FILE WITH CLINVAR DATA
echo "Cross referencing with clinvar"
bedtools intersect -a $4 -b $3 -header > $1/clinvar_allfrombed.vcf
bedtools intersect -a $1/$outname -b $3 -header > $1/patient_dcm_final.vcf
bedtools intersect -b $1/patient_dcm_final.vcf -a $1/clinvar_allfrombed.vcf -header > $1/patient_intersect_clinvar.vcf
python3 parse_clnsig.py -i $1/patient_intersect_clinvar.vcf 2>&1 | tee $1/patient_simple_report.txt
cut -c 24- $1/patient_simple_report.txt

#### COVERAGE CALCULATION
out=$(basename $VARIANTS | sed -e 's/\.bam$//').bga.bed
final=$(basename $VARIANTS | sed -e 's/\.bam$//').join_final.bed
depths=$(basename $VARIANTS | sed -e 's/\.bam$//').depths.bed

echo "Calculating coverage"
samtools view -L $3 $2 -b > $1/new.bam
bedtools genomecov -ibam $1/new.bam -bga > $1/$out
bedtools intersect -loj -F 0.10 -a $3 -b $1/$out -bed > $1/$final
awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$6,$7,$4,$8)}' $1/$final > $1/$depths

for gene in $(cut -f4 $3 | sort -u | xargs)
	do
		echo "Making graphs for $gene"
		grep $gene $1/$depths > $1/${gene}_raw.txt
		python cov.py $1/${gene}_raw.txt $1/${gene}.txt
		./draw_depth.R $1/${gene}.txt $1/${gene}.png
	done

echo "Creating report"
convert  $1/patientname.txt $1/patient_simple_report.txt $1/*.png $1/patient_report_final.pdf
echo "Done!"
```
