while getopts i:h:o:s:d:c:v:f: flag
do
    case "${flag}" in
        i) illumina_vcf=${OPTARG};;
        h) hifi_vcf=${OPTARG};;
        o) ont_vcf=${OPTARG};;
        s) sample=${OPTARG};;
        d) output_directory=${OPTARG};;
		c) variant_caller=${OPTARG};;
        v) version=${OPTARG};;
		f) filter=${OPTARG};;
    esac
done

# script to generate benchmarking VCFs

echo "illumina_vcf: $illumina_vcf";
echo "hifi_vcf: $hifi_vcf";
echo "ont_vcf: $ont_vcf";
echo "sample: $sample";
echo "output_directory: $output_directory";
echo "variant_caller: $variant_caller";
echo "variant caller options: deepsomatic or clairs"
echo "variant caller version: $version";
echo "filter: $filter";
echo "filter options: 'filter4', 'orthogonal_technology'"


set -o pipefail
set -e
set -u

# input files should be filtered for somatic-only

if [[ $filter == "filter4" ]]
then
	#########################################
	# generate filter4 benchmark 
	#########################################

	## merge
	bcftools merge --force-samples --threads 16 ${illumina_vcf} ${hifi_vcf} ${ont_vcf} | bgzip > ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_ONT_somaticOnly_merged.vcf.gz
	bcftools index -t ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_ONT_somaticOnly_merged.vcf.gz


	## filter4 on somatic only merged

	TOOL=/private/groups/patenlab/jimin/scripts/deepsomatic/vcf_intersection_complex_v2.py

	MERGE_VCF=${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_ONT_somaticOnly_merged.vcf.gz
	OUTPUT=${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_ONT_somaticOnly_filter4.vcf.gz

	python3 ${TOOL} -v ${MERGE_VCF} -i ${MERGE_VCF}.tbi -f 'filter4' -o ${OUTPUT}
	bcftools index -t ${OUTPUT}
elif [[ $filter == "orthogonal_technology" ]]
then
	################################################
	# generate orthogonal technology benchmark 
	################################################

	#################################################### MERGE ####################################################

	# WGS + PACBIO ###################
	bcftools merge --force-samples --threads 4 ${illumina_vcf} ${hifi_vcf} | bgzip > ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_somaticOnly_merged.vcf.gz
	bcftools index -t ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_somaticOnly_merged.vcf.gz

	# PACBIO + ONT ###################
	bcftools merge --force-samples --threads 4 ${hifi_vcf} ${ont_vcf} | bgzip > ${output_directory}/${sample}_${variant_caller}_${version}_PacBio_ONT_somaticOnly_merged.vcf.gz
	bcftools index -t ${output_directory}/${sample}_${variant_caller}_${version}_PacBio_ONT_somaticOnly_merged.vcf.gz

	# WGS + ONT ###################
	bcftools merge --force-samples --threads 4 ${illumina_vcf} ${ont_vcf} | bgzip > ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_ONT_somaticOnly_merged.vcf.gz
	bcftools index -t ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_ONT_somaticOnly_merged.vcf.gz

	#################################################### INTERSECTION ####################################################

	if [[ $variant_caller == "deepsomatic" ]]
	then
		FILTER='two_tech'
	elif [[ $variant_caller == "clairs" ]]
	then
		FILTER='two_tech_clairs'
	fi

	# WGS + PACBIO ###################
	VCF=${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_somaticOnly_merged.vcf.gz

	python3 /private/groups/patenlab/jimin/scripts/deepsomatic/working/vcf_comprehensive_update.py --vcf ${VCF} --vcfindex ${VCF}.tbi \
	--filter ${FILTER} --outfile ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_somaticOnly_intersection.vcf

	bgzip ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_somaticOnly_intersection.vcf
	bcftools index -t ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_PacBio_somaticOnly_intersection.vcf.gz


	# PACBIO + ONT ###################
	VCF=${output_directory}/${sample}_${variant_caller}_${version}_PacBio_ONT_somaticOnly_merged.vcf.gz

	python3 /private/groups/patenlab/jimin/scripts/deepsomatic/working/vcf_comprehensive_update.py --vcf ${VCF} --vcfindex ${VCF}.tbi \
	--filter ${FILTER} --outfile ${output_directory}/${sample}_${variant_caller}_${version}_PacBio_ONT_somaticOnly_intersection.vcf

	bgzip ${output_directory}/${sample}_${variant_caller}_${version}_PacBio_ONT_somaticOnly_intersection.vcf
	bcftools index -t ${output_directory}/${sample}_${variant_caller}_${version}_PacBio_ONT_somaticOnly_intersection.vcf.gz


	# WGS + ONT ###################
	VCF=${output_directory}/${sample}_${variant_caller}_${version}_Illumina_ONT_somaticOnly_merged.vcf.gz

	python3 /private/groups/patenlab/jimin/scripts/deepsomatic/working/vcf_comprehensive_update.py --vcf ${VCF} --vcfindex ${VCF}.tbi \
	--filter ${FILTER} --outfile ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_ONT_somaticOnly_intersection.vcf

	bgzip ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_ONT_somaticOnly_intersection.vcf
	bcftools index -t ${output_directory}/${sample}_${variant_caller}_${version}_Illumina_ONT_somaticOnly_intersection.vcf.gz
fi



















