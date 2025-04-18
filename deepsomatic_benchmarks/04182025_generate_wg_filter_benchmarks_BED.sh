
# check that all files exist in this naming format
for sample in 1395 1437 1937 1954 2009
do
    bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_Illumina_HiFi_ONT_MajorityCallable_merge_minusChrUn_random_minusSD_minusSVs_chr1_22_sort_merge.bed
    ls -lh $bed_file
done

 # files with LoH regions
for sample in 1395 2009
do
    if [[ $sample -eq 1395 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        other_regions=/private/groups/patenlab/jimin/data/high_conf_regions_working/1395_working/subtract_regions_chr6chr16.bed
        name_other_regions="LoH"
    elif [[ $sample -eq 1437 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        other_regions=""
        name_other_regions=""
    elif [[ $sample -eq 1937 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        other_regions=""
        name_other_regions=""
    elif [[ $sample -eq 1954 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        other_regions=""
        name_other_regions=""
    elif [[ $sample -eq 2009 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        other_regions=/private/groups/patenlab/jimin/data/high_conf_regions_working/2009_working/subtract_LoH_regions.bed
        name_other_regions="LoH"
    fi

    bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_Illumina_HiFi_ONT_MajorityCallable_merge_minusChrUn_random_minusSD_minusSVs_chr1_22_sort_merge.bed
    merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/MERGED/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_merged.vcf.gz
    truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_filter4.vcf.gz
    output_prefix=deepsomatic_v1.8.0_filter4
    output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
    date=04142025
    
    /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
    -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
    -j ${sv_bed1} -k ${sv_bed2} -x ${other_regions} -y ${name_other_regions}
done

for sample in 1437 1937 1954
do
    if [[ $sample -eq 1395 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        other_regions=/private/groups/patenlab/jimin/data/high_conf_regions_working/1395_working/subtract_regions_chr6chr16.bed
        name_other_regions="LoH"
    elif [[ $sample -eq 1437 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        other_regions=""
        name_other_regions=""
    elif [[ $sample -eq 1937 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        other_regions=""
        name_other_regions=""
    elif [[ $sample -eq 1954 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        other_regions=""
        name_other_regions=""
    elif [[ $sample -eq 2009 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        other_regions=/private/groups/patenlab/jimin/data/high_conf_regions_working/2009_working/subtract_LoH_regions.bed
        name_other_regions="LoH"
    fi

    bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_Illumina_HiFi_ONT_MajorityCallable_merge_minusChrUn_random_minusSD_minusSVs_chr1_22_sort_merge.bed
    merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/MERGED/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_merged.vcf.gz
    truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_filter4.vcf.gz
    output_prefix=deepsomatic_v1.8.0_filter4
    output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
    date=04142025
    
    /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
    -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
    -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
done

# check that all files exist in this naming format
for sample in 578 HG008
do
    bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_Illumina_HiFi_ONT_MajorityCallable_merge_minusChrUn_random_minusSD_minusSVs_chr1_22_sort_merge.bed
    ls -lh $bed_file
done

for sample in 578 HG008
do
    if [[ $sample -eq 578 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_Hs578T.bp1000.bed
        sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_Hs578Bst.bp1000.bed
    elif [[ $sample -eq HG008 ]]
    then
        sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HG008.bp1000.bed
    fi

    bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_Illumina_HiFi_ONT_MajorityCallable_merge_minusChrUn_random_minusSD_minusSVs_chr1_22_sort_merge.bed
    merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/MERGED/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_merged.vcf.gz
    truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_filter4.vcf.gz
    output_prefix=deepsomatic_v1.8.0_filter4

    output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
    date=04142025

    other_regions=""
    name_other_regions=""

    /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
    -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
    -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
done

# check that all files exist in this naming format
for sample in UPN237
do
    bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_Illumina_HiFi_ONT_MajorityCallable_merge_minusChrUn_random_minusSD_minusSVs_chr1_22_sort_merge.bed
    ls -lh $bed_file
done

# glioblastoma old tumor version
sample=UPN237

sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_ont_60x.bp1000.bed
sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_hifi.bp1000.bed

bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_Illumina_HiFi_ONT_MajorityCallable_merge_minusChrUn_random_minusSD_minusSVs_chr1_22_sort_merge.bed
merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/MERGED/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_merged.vcf.gz
truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_filter4.vcf.gz
output_prefix=deepsomatic_v1.8.0_filter4

output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
date=04142025

other_regions=""
name_other_regions=""

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
-b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
-j ${sv_bed1} -k ${sv_bed2} -x "" -y ""


# new Illumina tumor version
sample=UPN237_new_tumor

sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_ont_60x.bp1000.bed
sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_hifi.bp1000.bed

bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/UPN237_working/UPN237_Illumina_HiFi_ONT_MajorityCallable_merge_minusChrUn_random_minusSD_minusSVs_chr1_22_sort_merge.bed
merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/MERGED/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_merged.vcf.gz
truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/filter4_benchmark/${sample}_deepsomatic_v1.8.0_Illumina_PacBio_ONT_somaticOnly_filter4.vcf.gz
output_prefix=deepsomatic_v1.8.0_filter4

output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
date=04142025

other_regions=""
name_other_regions=""

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
-b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
-j ${sv_bed1} -k ${sv_bed2} -x "" -y ""


mv 04142025_1395_deepsomatic_v1.8.0_filter4_highconf_minusSVs_SDs_LoH_chr1_22.bed 04142025_1395_deepsomatic_v1.8.0_filter4_highconf_minusSVs_SDs_chr1_22.bed
mv 04142025_2009_deepsomatic_v1.8.0_filter4_highconf_minusSVs_SDs_LoH_chr1_22.bed 04142025_2009_deepsomatic_v1.8.0_filter4_highconf_minusSVs_SDs_chr1_22.bed
















