while getopts b:m:t:s:p:o:d:j:k:x:y: flag
do
    case "${flag}" in
        b) bed_file=${OPTARG};;
        m) merged_vcf=${OPTARG};;
        t) truth_vcf=${OPTARG};;
        s) sample=${OPTARG};;
        p) output_prefix=${OPTARG};;
        o) output_directory=${OPTARG};;
		d) date=${OPTARG};;
        j) sv_bed1=${OPTARG};;
        k) sv_bed2=${OPTARG};;
        x) other_regions=${OPTARG};;
        y) name_other_regions=${OPTARG};;

    esac
done

echo "bed_file: $bed_file";
echo "merged_vcf: $merged_vcf";
echo "truth_vcf: $truth_vcf";
echo "sample: $sample";
echo "output_prefix: $output_prefix";
echo "output_directory: $output_directory";
echo "date: $date";

# optional parameters
if [[ "$sv_bed1" != "" ]] ; then
echo "sv_bed1: $sv_bed1"; fi

if [[ "$sv_bed2" != "" ]] ; then
echo "sv_bed2: $sv_bed2"; fi

if [[ "$other_regions" != "" ]] ; then
echo "other_regions: $other_regions"; fi

if [[ "$name_other_region" != "" ]] ; then
echo "name_other_regions: $name_other_regions"; fi


set -o pipefail
set -e
set -u

mkdir -p ${output_directory}

# 1. subtract "confusing regions" and generate "high conf regions bed"
/private/groups/patenlab/jimin/scripts/deepsomatic/subtract_confusing_regions.sh \
-c ${bed_file} \
-m ${merged_vcf} \
-t ${truth_vcf} \
-s ${sample} \
-p ${output_prefix} \
-o ${output_directory} \
-d ${date}

# output from subtract_confusing_regions.sh
HIGH_CONF=${output_directory}/${sample}/${date}_${sample}_${output_prefix}_highconf.bed

# 2. subtract SV regions and SD regions

SD_BED=/private/groups/patenlab/jimin/data/BED/GIAB/GRCh38_segdups.bed

if [[ -z "$sv_bed1" ]]
then
    bedtools subtract -a ${HIGH_CONF} -b ${SD_BED} > ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSDs.bed
    CURR_BED=${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSDs.bed
elif [[ -z "$sv_bed2" ]]
then
    bedtools subtract -a ${HIGH_CONF} -b ${sv_bed1} | bedtools subtract -a stdin -b ${SD_BED} > ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs.bed
    CURR_BED=${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs.bed
else
    bedtools subtract -a ${HIGH_CONF} -b ${sv_bed1} | bedtools subtract -a stdin -b ${sv_bed2} | bedtools subtract -a stdin -b ${SD_BED} > ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs.bed
    CURR_BED=${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs.bed
fi

# 3. subtract any other regions

if [[ "$other_regions" != "" ]]
then
    if [[ "${CURR_BED}" == ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs.bed ]]
    then
        bedtools subtract -a ${CURR_BED} -b ${other_regions} > ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs_${name_other_regions}.bed

        # update CURR_BED
        CURR_BED=${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs_${name_other_regions}.bed
    elif [[ "${CURR_BED}" == ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSDs.bed ]]
    then
        bedtools subtract -a ${CURR_BED} -b ${other_regions} > ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSDs_${name_other_regions}.bed

        # update CURR_BED
        CURR_BED=${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSDs_${name_other_regions}.bed
fi


# 4. generate autosome only verions (chr1-22)
# (if statements just for correct file name)
if [[ "${CURR_BED}" == ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs.bed ]]
then
    grep -v "chrX" ${CURR_BED} | grep -v "chrY" > ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs_chr1_22.bed
elif [[ "${CURR_BED}" == ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSDs.bed ]]
then
    grep -v "chrX" ${CURR_BED} | grep -v "chrY" > ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSDs_chr1_22.bed
elif [[ "${CURR_BED}" == ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs_${name_other_regions}.bed ]]
then 
    grep -v "chrX" ${CURR_BED} | grep -v "chrY" > ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSVs_SDs_${name_other_regions}_chr1_22.bed
elif
    [[ "${CURR_BED}" == ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSDs_${name_other_regions}.bed ]]
then
    grep -v "chrX" ${CURR_BED} | grep -v "chrY" > ${output_directory}/${date}_${sample}_${output_prefix}_highconf_minusSDs_${name_other_regions}_chr1_22.bed
fi