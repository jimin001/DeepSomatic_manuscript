# test benchmark generating scripts

################
# 578
################

illumina_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/WGS/1017_578_DeepSomatic_v1.7.0_WGS_wg_somaticOnly.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/PacBio/1107_578_DeepSomatic_v1.7.0_PacBio_wg_somaticOnly.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/ONT/1017_578_DeepSomatic_v1.7.0_ONT_wg_somaticOnly.vcf.gz
sample=578
output_directory=/private/groups/patenlab/jimin/data/VCF/04042025_test_benchmark_script
variant_caller=deepsomatic
version=v1.7.0
filter="orthogonal_technology"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}

# bcftools view 578_deepsomatic_v1.7.0_Illumina_ONT_somaticOnly_intersection.vcf.gz | wc -l
# 12448

# bcftools view 578_deepsomatic_v1.7.0_Illumina_ONT_somaticOnly_merged.vcf.gz | wc -l
# 19915

# bcftools view 578_v1.7.0_Illumina_ONT_somaticOnly_intersection.vcf.gz | wc -l
# 12448

# bcftools view 1017_578_DeepSomatic_v1.7.0_Illumina_ONT_somaticOnly_merged.vcf.gz | wc -l
# 19915

################
# 1437
################

illumina_vcf=/private/groups/patenlab/jimin/data/VCF/0701_final_docker/WGS/0701_1437_DeepSomatic_WGS_wg_somaticOnly.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/data/VCF/0701_final_docker/PacBio/0701_1437_DeepSomatic_PacBio_wg_somaticOnly.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/data/VCF/0701_final_docker/ONT/0701_1437_DeepSomatic_ONT_wg_somaticOnly.vcf.gz
sample=1437
output_directory=/private/groups/patenlab/jimin/data/VCF/04042025_test_benchmark_script
variant_caller=deepsomatic
version=v1.7.0
filter="filter4"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}

# bcftools view 0704_1437_DeepSomatic_Illumina_PacBio_ONT_somaticOnly_filter4.vcf.gz | wc -l
# 92632

# bcftools view 1437_deepsomatic_v1.7.0_Illumina_PacBio_ONT_somaticOnly_filter4.vcf.gz | wc -l
# 92632

#################################################################################################################################################

bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/578_working/orthogonal_tech_beds/578_MajorityCallable_Illumina+ONT_concat_sort_merge_minusChrUn_random.bed
merged_vcf=/private/groups/patenlab/jimin/data/VCF/04042025_test_benchmark_script/578_deepsomatic_v1.7.0_Illumina_ONT_somaticOnly_merged.vcf.gz
truth_vcf=/private/groups/patenlab/jimin/data/VCF/04042025_test_benchmark_script/578_deepsomatic_v1.7.0_Illumina_ONT_somaticOnly_intersection.vcf.gz
sample=578
output_prefix=deepsomatic_v1.7.0_Illumina_ONT
output_directory=/private/groups/patenlab/jimin/data/BED/04042025_test_benchmark_script
date=04042025
sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_Hs578T.bp1000.bed
sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_Hs578Bst.bp1000.bed
other_regions=""
name_other_regions=""

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
-b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
-j ${sv_bed1} -k ${sv_bed2} -x "" -y ""

# 50066 04042025_578_deepsomatic_v1.7.0_Illumina_ONT_highconf_minusSVs_SDs.bed
# 47400 04042025_578_deepsomatic_v1.7.0_Illumina_ONT_highconf_minusSVs_SDs_chr1_22.bed

# 47400 578_Illumina+ONT_highconf_selfgenprior_minusSVs_SDs.bed


bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/578_working/orthogonal_tech_beds/578_MajorityCallable_HiFi+ONT_concat_sort_merge_minusChrUn_random.bed
merged_vcf=/private/groups/patenlab/jimin/data/VCF/04042025_test_benchmark_script/578_deepsomatic_v1.7.0_PacBio_ONT_somaticOnly_merged.vcf.gz
truth_vcf=/private/groups/patenlab/jimin/data/VCF/04042025_test_benchmark_script/578_deepsomatic_v1.7.0_PacBio_ONT_somaticOnly_intersection.vcf.gz
sample=578
output_prefix=deepsomatic_v1.7.0_PacBio_ONT
output_directory=/private/groups/patenlab/jimin/data/BED/04042025_test_benchmark_script
date=04042025
sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_Hs578T.bp1000.bed
sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_Hs578Bst.bp1000.bed
other_regions=""
name_other_regions=""

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
-b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
-j ${sv_bed1} -k ${sv_bed2} -x "" -y ""

# 48555 04042025_578_deepsomatic_v1.7.0_PacBio_ONT_highconf_minusSVs_SDs_chr1_22.bed
# 48555 578_HiFi+ONT_highconf_selfgenprior_minusSVs_SDs.bed


















