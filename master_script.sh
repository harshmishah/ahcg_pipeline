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
#echo ".Depth of Coverage" >> $1/depth-chart.txt
for gene in $(cut -f4 $3 | sort -u | xargs)
	do
		echo "Making graphs for $gene"
		grep $gene $1/$depths > $1/${gene}_raw.txt
		python cov.py $1/${gene}_raw.txt $1/${gene}.txt
		#xvfb-run --server-args="-screen 0 1024x768x24" ./draw_depth.R $1/${gene}.txt
		./draw_depth.R $1/${gene}.txt $1/${gene}.png
		#echo "image::${gene}coverage.png[height=408]" >> $1/depth-chart.txt
	done

echo "Creating report"
convert  $1/patientname.txt $1/patient_simple_report.txt $1/*.png $1/patient_report_final.pdf
echo "Done!"


