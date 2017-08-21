#!/bin/bash

#------------------------------------------------------------------------------------------------------------------------
# Name		: reproduceResults_complete.sh
# Author	: Mathijs A Sanders
# Version	: 1.0
# Copyright	: None
# Description	: reproduceResults_complete.sh -i input_directory -o output_directory  [-g reference_genome] [-h]
#------------------------------------------------------------------------------------------------------------------------

display_usage() {
cat << HEREDOC
Usage: reproduceResults_complete.sh -i input_directory -m wgbs_directory -o output_directory [-g reference_genome] [-h]

Required arguments:

	-i <input_directory>		Input directory containing the JAR file and preprocessed (epi)genomic feature interval files
	-m <wgbs_directory>		Input directory containing the CD34+ coverage and methylation bigWig files (Spencer et al., https://wustl.app.box.com/s/6fvifntz6golh4vu52epkfqxeyd0vqmk, CD34_[12456].cpg.filtered.(coverage/meth).bw)
	-o <output_directory>		Output directory for storing all downloaded input files, determined tables and RMR output files
	
Optional arguments:
	
	-g <reference_genome>		Input FASTA file used for aligning whole genome sequencing data (Optional)

HEREDOC
}
#--------------------------------------------------------------------------------------------------
# Get input parameters
#--------------------------------------------------------------------------------------------------


while getopts ':i:m:o:g:h' ARG
do
	case $ARG in

		i)
			INPUT_PREFIX=$OPTARG
		;;
		
		m)
			METHYLATION_PREFIX=$OPTARG
		;;

		o)
			OUTPUT_PREFIX=$OPTARG
		;;

		g)
			REF=$OPTARG
		;;

		h)
			display_usage
			exit 0
		;;

		\?)
			echo "Error: Invalid option -$OPTARG (-h for help)"
			display_usage
			exit -1
		;;
	esac
done

if [[ -z $INPUT_PREFIX ]]
then
	echo "Please specify an input directory containing the RoaMeR JAR file and preprocessed (epi)genomic feature interval files (-h for help)"
	exit -1
fi

if [[ -z $METHYLATION_PREFIX ]]
then
	echo "Please specify the input directory containing the CD34+ coverage and methylation bigWig files (Spencer et al)"
	exit -1
fi

if [[ -z $OUTPUT_PREFIX ]]
then
	echo "Please specify an output directory (-h for help)"
	exit -1
fi

if [[ ! -f ${INPUT_PREFIX}/RoaMeR.jar ]]; then
	echo "Please specify the input directory containing the RoaMeR JAR file"
	exit -1
else
	JAR=${INPUT_PREFIX}/RoaMeR.jar
fi

if [[ -z $REF ]] 
then
	mkdir -p ${OUTPUT_PREFIX}/Reference
	for i in {1..22} X
	do
		wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz -P ${OUTPUT_PREFIX}/Reference/
		zcat ${OUTPUT_PREFIX}/Reference/chr${i}.fa.gz >> ${OUTPUT_PREFIX}/Reference/hg19.fa
	done
	export REF=${OUTPUT_PREFIX}/Reference/hg19.fa
	if [[ ! -f $REF ]] 
	then
		echo "Error: Downloading and compiling the reference sequence file has failed"
		exit -1
	else
		rm -f ${OUTPUT_PREFIX}/Reference/*.fa.gz
	fi
elif [[ ! -f $REF ]] 
then
	echo "Error: Specified FASTQ file \"$g\" is not a valid file or does not exist"
	exit -1
fi

#------------------------------------------------------------------------------------------------
# Download R-code, decompress and reproduce Figures 1a,b,d; 2d and Supplemantary Figure 3a
#------------------------------------------------------------------------------------------------

DIR=${OUTPUT_PREFIX}/plotReproduction/

if [[ ! -d $DIR ]]
then
        wget https://gitlab.wehi.edu.au/flensburg.c/codeForPaper/raw/master/Ranalysis.tar.gz -P ${OUTPUT_PREFIX}/
        tar -xzf ${OUTPUT_PREFIX}/Ranalysis.tar.gz && rm ${OUTPUT_PREFIX}/Ranalysis.tar.gz
        if [[ ! -d $DIR ]]
        then
                echo "Error: Unable to download and decompress the R code"
                exit -1
        fi
fi

CDIR=$(pwd)
cd ${OUTPUT_PREFIX}/plotReproduction/
Rscript runAll.R

if [[ $? -ne 0 ]]
then
        echo "Error: Downloading dependencies and running the R code has failed"
        exit -1
else
        cd $CDIR
fi


#------------------------------------------------------------------------------------------------
# Construct genome-wide CG table used for extracting methylation data and calculating RMR scores
#------------------------------------------------------------------------------------------------

mkdir -p ${OUTPUT_PREFIX}/db

java -Xmx29G -jar $JAR --analysis substrate --inputFasta $REF --outputDatabase ${OUTPUT_PREFIX}/db/cg_table.txt --substrate CG --bufferSize 10 --excludeContigs $(pwd)/Auxiliary/ExcludeChromosomes.txt --verbose

if [[ $? -ne 0 ]]
then
	echo "Error: Creating the substrate table has failed"
	exit -1
fi

#------------------------------------------------------------------------------------------------
# Extract methylation data from the bigWig files and construct a 5mCG table
#------------------------------------------------------------------------------------------------

mkdir -p ${OUTPUT_PREFIX}/methylation/

if [[ ! -f ${METHYLATION_PREFIX}/CD34_1.cpg.filtered.coverage.bw ]]
then
	echo "Error: Could not locate CD34_1.cpg.filtered.coverage.bw"
	exit -1
fi

if [[ ! -f ${METHYLATION_PREFIX}/CD34_1.cpg.filtered.meth.bw ]]
then
	echo "Error: Could not locate CD34_1.cpg.filtered.meth.bw"
	exit -1
fi


if [[ ! -f ${METHYLATION_PREFIX}/CD34_2.cpg.filtered.coverage.bw ]]
then
	echo "Error: Could not locate CD34_2.cpg.filtered.coverage.bw"
	exit -1
fi

if [[ ! -f ${METHYLATION_PREFIX}/CD34_2.cpg.filtered.meth.bw ]]
then
	echo "Error: Could not locate CD34_2.cpg.filtered.meth.bw"
	exit -1
fi

if [[ ! -f ${METHYLATION_PREFIX}/CD34_4.cpg.filtered.coverage.bw ]]
then
	echo "Error: Could not locate CD34_4.cpg.filtered.coverage.bw"
	exit -1
fi


if [[ ! -f ${METHYLATION_PREFIX}/CD34_4.cpg.filtered.meth.bw ]]
then
	echo "Error: Could not locate CD34_4.cpg.filtered.meth.bw"
	exit -1
fi

if [[ ! -f ${METHYLATION_PREFIX}/CD34_5.cpg.filtered.coverage.bw ]]
then
	echo "Error: Could not locate CD34_5.cpg.filtered.coverage.bw"
	exit -1
fi

if [[ ! -f ${METHYLATION_PREFIX}/CD34_5.cpg.filtered.meth.bw ]]
then
	echo "Error: Could not locate CD34_5.cpg.filtered.meth.bw"
	exit -1
fi

if [[ ! -f ${METHYLATION_PREFIX}/CD34_6.cpg.filtered.coverage.bw ]]
then
	echo "Error: Could not locate CD34_6.cpg.filtered.coverage.bw"
	exit -1
fi

if [[ ! -f ${METHYLATION_PREFIX}/CD34_6.cpg.filtered.meth.bw ]]
then
	echo "Error: Could not locate CD34_6.cpg.filtered.meth.bw"
	exit -1
fi

#Build triplet 5mCG table

java -Xmx29G -jar $JAR --analysis extractmethylation --coverageFiles ${METHYLATION_PREFIX}/CD34_1.cpg.filtered.coverage.bw,${METHYLATION_PREFIX}/CD34_2.cpg.filtered.coverage.bw,${METHYLATION_PREFIX}/CD34_4.cpg.filtered.coverage.bw,${METHYLATION_PREFIX}/CD34_5.cpg.filtered.coverage.bw,${METHYLATION_PREFIX}/CD34_6.cpg.filtered.coverage.bw --methylationFiles ${METHYLATION_PREFIX}/CD34_1.cpg.filtered.meth.bw,${METHYLATION_PREFIX}/CD34_2.cpg.filtered.meth.bw,${METHYLATION_PREFIX}/CD34_4.cpg.filtered.meth.bw,${METHYLATION_PREFIX}/CD34_5.cpg.filtered.meth.bw,${METHYLATION_PREFIX}/CD34_6.cpg.filtered.meth.bw --outputMethylation ${OUTPUT_PREFIX}/methylation/CD34_meth_triplet.txt --cgDatabase ${OUTPUT_PREFIX}/db/cg_table.txt --precedingBases 1 --succeedingBases 0 --verbose --excludeContigs ${INPUT_PREFIX}/Auxiliary/excludeContigs.txt

if [[ $? -ne 0 ]]
then
	echo "Error: Building of the triplet 5mCG database has failed"
	exit -1
fi

#Build quadruplet 5mCG table

java -Xmx29G -jar $JAR --analysis extractmethylation --coverageFiles ${METHYLATION_PREFIX}/CD34_1.cpg.filtered.coverage.bw,${METHYLATION_PREFIX}/CD34_2.cpg.filtered.coverage.bw,${METHYLATION_PREFIX}/CD34_4.cpg.filtered.coverage.bw,${METHYLATION_PREFIX}/CD34_5.cpg.filtered.coverage.bw,${METHYLATION_PREFIX}/CD34_6.cpg.filtered.coverage.bw --methylationFiles ${METHYLATION_PREFIX}/CD34_1.cpg.filtered.meth.bw,${METHYLATION_PREFIX}/CD34_2.cpg.filtered.meth.bw,${METHYLATION_PREFIX}/CD34_4.cpg.filtered.meth.bw,${METHYLATION_PREFIX}/CD34_5.cpg.filtered.meth.bw,${METHYLATION_PREFIX}/CD34_6.cpg.filtered.meth.bw --outputMethylation ${OUTPUT_PREFIX}/methylation/CD34_meth_quadruplet.txt --cgDatabase ${OUTPUT_PREFIX}/db/cg_table.txt --precedingBases 1 --succeedingBases 1 --verbose --excludeContigs ${INPUT_PREFIX}/Auxiliary/excludeContigs.txt

if [[ $? -ne 0 ]]
then
	echo "Error: Building of the quadruplet 5mCG database has failed"
	exit -1
fi

#------------------------------------------------------------------------------------------------
# Calculate RMR and RMR-CG scores for the (epi)genomic features
#------------------------------------------------------------------------------------------------


mkdir -p ${OUTPUT_PREFIX}/results/Genome_triplet_2B

java -Xmx29G -jar $JAR --analysis rmr --database ${OUTPUT_PREFIX}/methylation/CD34_meth_triplet.txt --bedfiles ${INPUT_PREFIX}/Intervals/Genome/Genome.bed --namesBedFiles Genome --vcfs ${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-1.vcf,${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-2.vcf --sampleNames WEHI-AML-1,WEHI-AML-2 --outputStatistics ${OUTPUT_PREFIX}/results/Genome_triplet_2B/Genome_triplet.txt --outputContext ${OUTPUT_PREFIX}/results/Genome_triplet_2B/Genome_triplet_context.txt --verbose

if [[ $? -ne 0 ]]
then
	echo "Error: Error calculating RMR scores for genome-wide triplets"
	exit -1
fi

mkdir -p ${OUTPUT_PREFIX}/results/Features_2C

java -Xmx29G -jar $JAR --analysis rmr --database ${OUTPUT_PREFIX}/methylation/CD34_meth_quadruplet.txt --bedfiles ${INPUT_PREFIX}/Intervals/Features/cpgIsland.bed,${INPUT_PREFIX}/Intervals/Features/Promoter.bed,${INPUT_PREFIX}/Intervals/Features/Promoter_flanking.bed,${INPUT_PREFIX}/Intervals/Features/Ensembl_exon.bed,${INPUT_PREFIX}/Intervals/Features/Ensembl_intron.bed --namesBedFiles CpGisland,Promoter,Promoter_flanking,Ensembl_exon,Ensembl_intron --vcfs ${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-1.vcf,${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-2.vcf --sampleNames WEHI-AML-1,WEHI-AML-2 --outputStatistics ${OUTPUT_PREFIX}/results/Features_2C/Features.txt --verbose

if [[ $? -ne 0 ]]
then
	echo "Error: Error calculating RMR scores for the (epi)genomic features"
	exit -1
fi

mkdir -p ${OUTPUT_PREFIX}/results/Genome_quadruplet_2E

java -Xmx29G -jar $JAR --analysis rmr --database ${OUTPUT_PREFIX}/methylation/CD34_meth_quadruplet.txt --bedfiles ${INPUT_PREFIX}/Intervals/Genome/Genome.bed --namesBedFiles Genome --vcfs ${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-1.vcf,${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-2.vcf --sampleNames WEHI-AML-1,WEHI-AML-2 --outputStatistics ${OUTPUT_PREFIX}/results/Genome_quadruplet_2E/Genome_quadruplet.txt --outputContext ${OUTPUT_PREFIX}/results/Genome_quadruplet_2E/Genome_quadruplet_context.txt --verbose

if [[ $? -ne 0 ]]
then
	echo "Error: Error calculating RMR scores for genome-wide quadruplets"
	exit -1
fi


mkdir -p ${OUTPUT_PREFIX}/results/Replication_Domain_context_2F

java -Xmx29G -jar $JAR --analysis rmr --database ${OUTPUT_PREFIX}/methylation/CD34_meth_quadruplet.txt --bedfiles ${INPUT_PREFIX}/Intervals/ReplicationDomains/Earliest_Equal.bed,${INPUT_PREFIX}/Intervals/ReplicationDomains/Early_Equal.bed,${INPUT_PREFIX}/Intervals/ReplicationDomains/Late_Equal.bed,${INPUT_PREFIX}/Intervals/ReplicationDomains/Latest_Equal.bed --namesBedFiles Earliest,Early,Late,Latest --vcfs ${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-1.vcf,${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-2.vcf --sampleNames WEHI-AML-1,WEHI-AML-2 --outputStatistics ${OUTPUT_PREFIX}/results/Replication_Domain_context_2F/Replication_domain.txt --outputContext ${OUTPUT_PREFIX}/results/Replication_Domain_context_2F/Replication_domain_context.txt --verbose

if [[ $? -ne 0 ]]
then
	echo "Error: Error calculating RMR scores for the replication domains"
	exit -1
fi

mkdir -p ${OUTPUT_PREFIX}/results/Strand_bias_Supp5A

java -Xmx29G -jar $JAR --analysis rmr --database ${OUTPUT_PREFIX}/methylation/CD34_meth_quadruplet.txt --bedfiles ${INPUT_PREFIX}/Intervals/Strand/Ensembl_gene_body_strand.bed --namesBedFiles Strand --vcfs ${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-1.vcf,${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-2.vcf --sampleNames WEHI-AML-1,WEHI-AML-2 --outputStatistics ${OUTPUT_PREFIX}/results/Strand_bias_Supp5A/Strand_bias.txt --verbose --stranded

if [[ $? -ne 0 ]]
then
	echo "Error: Error calculating RMR scores for the transcriptional strand bias"
	exit -1
fi

mkdir -p ${OUTPUT_PREFIX}/results/Expression_strand_bias_Supp5B

java -Xmx29G -jar $JAR --analysis rmr --database ${OUTPUT_PREFIX}/methylation/CD34_meth_quadruplet.txt --bedfiles ${INPUT_PREFIX}/Intervals/Expression_Strand/Expression_none.bed,${INPUT_PREFIX}/Intervals/Expression_Strand/Expression_lowest.bed,${INPUT_PREFIX}/Intervals/Expression_Strand/Expression_low.bed,${INPUT_PREFIX}/Intervals/Expression_Strand/Expression_high.bed,${INPUT_PREFIX}/Intervals/Expression_Strand/Expression_highest.bed --namesBedFiles None,Lowest,Low,High,Highest --vcfs ${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-1.vcf,${INPUT_PREFIX}/Intervals/VCF/WEHI-AML-2.vcf --sampleNames WEHI-AML-1,WEHI-AML-2 --outputStatistics ${OUTPUT_PREFIX}/results/Expression_strand_bias_Supp5B/Expression_strand_bias.txt --verbose --stranded

if [[ $? -ne 0 ]]
then
	echo "Error: Error calculating RMR scores for the expression bins"
	exit -1
fi

mkdir -p ${OUTPUT_PREFIX}/results/gnomAD_Supp5C

if [[ ! -f ${INPUT_PREFIX}/Intervals/gnomAD/gnomad_0.0001_0.001.vcf ]]
then
	if [[ ! -f ${INPUT_PREFIX}/Intervals/gnomAD/gnomad_0.0001_0.001.vcf.gz ]]
	then
		echo "Error: gnomAD source files do not exist"
		exit -1
	else
		gunzip -c ${INPUT_PREFIX}/Intervals/gnomAD/gnomad_0.0001_0.001.vcf.gz > ${INPUT_PREFIX}/Intervals/gnomAD/gnomad_0.0001_0.001.vcf
	fi
fi

java -Xmx29G -jar $JAR --analysis rmr --database ${OUTPUT_PREFIX}/methylation/CD34_meth_quadruplet.txt --bedfiles ${INPUT_PREFIX}/Intervals/ReplicationDomains/Earliest_Equal.bed,${INPUT_PREFIX}/Intervals/ReplicationDomains/Early_Equal.bed,${INPUT_PREFIX}/Intervals/ReplicationDomains/Late_Equal.bed,${INPUT_PREFIX}/Intervals/ReplicationDomains/Latest_Equal.bed --namesBedFiles Earliest,Early,Late,Latest --vcfs ${INPUT_PREFIX}/Intervals/gnomAD/gnomad_0.0001_0.001.vcf --sampleNames gnomAD --outputStatistics ${OUTPUT_PREFIX}/results/gnomAD_Supp5C/gnomAD_Replication_domain.txt --outputContext ${OUTPUT_PREFIX}/results/gnomAD_Supp5C/gnomAD_Replication_domain_context.txt --verbose

if [[ $? -ne 0 ]]
then
	echo "Error: Error calculating RMR scores for the replication domains using the rare gnomAD variants"
	exit -1
fi

exit 0

