# generate high confidence regions BED files, excluding "confusing regions" derived from orthogonal technology benchmark VCFs

for two_tech in "PacBio_ONT"
do
    for sample in 1395 1437 1937 1954 2009
    do
        if [[ $sample -eq 1395 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        elif [[ $sample -eq 1437 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        elif [[ $sample -eq 1937 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        elif [[ $sample -eq 1954 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        elif [[ $sample -eq 2009 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        fi

        bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_MajorityCallable_${two_tech}_concat_sort_merge_minusChrUn_random.bed
        merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_merged.vcf.gz
        truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_intersection.vcf.gz
        output_prefix=deepsomatic_v1.8.0_${two_tech}
        output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
        date=04142025
        
        other_regions=""
        name_other_regions=""

        /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
        -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
        -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
    done
done


for two_tech in "Illumina_ONT"
do
    for sample in 1395 1437 1937 1954 2009
    do
        if [[ $sample -eq 1395 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        elif [[ $sample -eq 1437 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        elif [[ $sample -eq 1937 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        elif [[ $sample -eq 1954 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        elif [[ $sample -eq 2009 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        fi

        bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_MajorityCallable_WGS_ONT_concat_sort_merge_minusChrUn_random.bed
        merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_merged.vcf.gz
        truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_intersection.vcf.gz
        output_prefix=deepsomatic_v1.8.0_${two_tech}
        output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
        date=04142025
        
        other_regions=""
        name_other_regions=""

        /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
        -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
        -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
    done
done

for two_tech in "Illumina_PacBio"
do
    for sample in 1395 1437 1937 1954 2009
    do
        if [[ $sample -eq 1395 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        elif [[ $sample -eq 1437 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        elif [[ $sample -eq 1937 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        elif [[ $sample -eq 1954 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}BL.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_HCC${sample}.bp1000.bed
        elif [[ $sample -eq 2009 ]]
        then
            sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_BL${sample}.bp1000.bed
            sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_H${sample}.bp1000.bed
        fi

        bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/${sample}_MajorityCallable_WGS_PacBio_concat_sort_merge_minusChrUn_random.bed
        merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_merged.vcf.gz
        truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_intersection.vcf.gz
        output_prefix=deepsomatic_v1.8.0_${two_tech}
        output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
        date=04142025
        
        other_regions=""
        name_other_regions=""

        /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
        -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
        -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
    done
done

# subtract LoH regions

# 1395 ############################

sample=1395
for two_tech in "Illumina_PacBio" "Illumina_ONT" "PacBio_ONT"
do
	cd /private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer/${two_tech}
	BED=04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs.bed
	LOH=/private/groups/patenlab/jimin/data/high_conf_regions_working/1395_working/subtract_regions_chr6chr16.bed
	bedtools subtract -a ${BED} -b ${LOH} > 04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed
done

for two_tech in "Illumina_PacBio" "Illumina_ONT" "PacBio_ONT"
do
	cd /private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer/${two_tech}
	mv 04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed 04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs.bed 
done

sample=1395
for two_tech in "Illumina_PacBio" "Illumina_ONT" "PacBio_ONT"
do
	cd /private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer/${two_tech}
	BED=04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed
	LOH=/private/groups/patenlab/jimin/data/high_conf_regions_working/1395_working/subtract_regions_chr6chr16.bed
	bedtools subtract -a ${BED} -b ${LOH} > 04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed
	wc -l 04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed
	mv 04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed 04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed
	wc -l 04142025_1395_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed
done

# 2009 ############################

sample=2009
for two_tech in "Illumina_PacBio" "Illumina_ONT" "PacBio_ONT"
do
	cd /private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer/${two_tech}
	BED=04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed
	LOH=/private/groups/patenlab/jimin/data/high_conf_regions_working/2009_working/subtract_LoH_regions.bed
	bedtools subtract -a ${BED} -b ${LOH} > 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed
	wc -l 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed
	mv 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed
	wc -l 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed
done

sample=2009
for two_tech in "Illumina_PacBio" "Illumina_ONT" "PacBio_ONT"
do
	cd /private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer/${two_tech}
	BED=04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs.bed
	LOH=/private/groups/patenlab/jimin/data/high_conf_regions_working/2009_working/subtract_LoH_regions.bed
	bedtools subtract -a ${BED} -b ${LOH} > 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed
	wc -l 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed
	mv 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_LoH.bed 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs.bed
	wc -l 04142025_${sample}_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs.bed
done


for two_tech in "Illumina_PacBio"
do
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

        bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/orthogonal_tech_beds/${sample}_MajorityCallable_Illumina+HiFi_concat_sort_merge_minusChrUn_random.bed
        merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_merged.vcf.gz
        truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_intersection.vcf.gz
        output_prefix=deepsomatic_v1.8.0_${two_tech}

        output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
        date=04142025

        other_regions=""
        name_other_regions=""

        /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
        -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
        -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
    done
done

two_tech="Illumina_PacBio"
grep -v "chrM" 04142025_HG008_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs.bed > 04142025_HG008_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_temp.bed
mv 04142025_HG008_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_temp.bed 04142025_HG008_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs.bed

grep -v "chrM" 04142025_HG008_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed > 04142025_HG008_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22_temp.bed
mv 04142025_HG008_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22_temp.bed 04142025_HG008_deepsomatic_v1.8.0_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed


for two_tech in "Illumina_ONT"
do
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

        bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/orthogonal_tech_beds/${sample}_MajorityCallable_Illumina+ONT_concat_sort_merge_minusChrUn_random.bed
        merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_merged.vcf.gz
        truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_intersection.vcf.gz
        output_prefix=deepsomatic_v1.8.0_${two_tech}

        output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
        date=04142025

        other_regions=""
        name_other_regions=""

        /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
        -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
        -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
    done
done


grep -v "chrM" 04142025_HG008_deepsomatic_v1.8.0_Illumina_ONT_highconf_minusSVs_SDs.bed > 04142025_HG008_deepsomatic_v1.8.0_Illumina_ONT_highconf_minusSVs_SDs_temp.bed
mv 04142025_HG008_deepsomatic_v1.8.0_Illumina_ONT_highconf_minusSVs_SDs_temp.bed 04142025_HG008_deepsomatic_v1.8.0_Illumina_ONT_highconf_minusSVs_SDs.bed

grep -v "chrM" 04142025_HG008_deepsomatic_v1.8.0_Illumina_ONT_highconf_minusSVs_SDs_chr1_22.bed > 04142025_HG008_deepsomatic_v1.8.0_Illumina_ONT_highconf_minusSVs_SDs_chr1_22_temp.bed
mv 04142025_HG008_deepsomatic_v1.8.0_Illumina_ONT_highconf_minusSVs_SDs_chr1_22_temp.bed 04142025_HG008_deepsomatic_v1.8.0_Illumina_ONT_highconf_minusSVs_SDs_chr1_22.bed



for two_tech in "PacBio_ONT"
do
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

        bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/orthogonal_tech_beds/${sample}_MajorityCallable_HiFi+ONT_concat_sort_merge_minusChrUn_random.bed
        merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_merged.vcf.gz
        truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_intersection.vcf.gz
        output_prefix=deepsomatic_v1.8.0_${two_tech}

        output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
        date=04142025

        other_regions=""
        name_other_regions=""

        /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
        -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
        -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
    done
done

sample=UPN237
for two_tech in "Illumina_PacBio" "Illumina_ONT" "PacBio_ONT"
do
    sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_ont_60x.bp1000.bed
    sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_hifi.bp1000.bed

    bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/UPN237_MajorityCallable_${two_tech}_concat_sort_merge_minusChrUn_random.bed
    merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_merged.vcf.gz
    truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_v1.8.0_${two_tech}_somaticOnly_intersection.vcf.gz
    output_prefix=deepsomatic_v1.8.0_${two_tech}

    output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
    date=04142025

    other_regions=""
    name_other_regions=""

    /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
    -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
    -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
done

# new Illumina tumor version
sample=UPN237_new_tumor
for two_tech in "Illumina_PacBio" "Illumina_ONT"
do
    sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_ont_60x.bp1000.bed
    sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_hifi.bp1000.bed

    bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/UPN237_working/UPN237_MajorityCallable_${two_tech}_concat_sort_merge_minusChrUn_random.bed
    merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/UPN237_new_tumor_deepsomatic_v1.8.0_${two_tech}_somaticOnly_merged.vcf.gz
    truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/UPN237_new_tumor_deepsomatic_v1.8.0_${two_tech}_somaticOnly_intersection.vcf.gz
    output_prefix=deepsomatic_v1.8.0_${two_tech}

    output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
    date=04142025

    other_regions=""
    name_other_regions=""

    /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
    -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
    -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
done


# 5/2/2025 update samples 578, HG008 and UPN237 new tumor with v17_rc0_07012024 docker for PacBio and ONT

version=v17_rc0_07012024
for two_tech in "Illumina_PacBio"
do
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

        bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/orthogonal_tech_beds/${sample}_MajorityCallable_Illumina+HiFi_concat_sort_merge_minusChrUn_random.bed
        merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_${version}_${two_tech}_somaticOnly_merged.vcf.gz
        truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_${version}_${two_tech}_somaticOnly_intersection.vcf.gz
        output_prefix=deepsomatic_${version}_${two_tech}

        output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
        date=05022025

        other_regions=""
        name_other_regions=""

        /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
        -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
        -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
    done
done

grep -v "chrM" ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs.bed > ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_temp.bed
mv ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_temp.bed ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs.bed

grep -v "chrM" ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed > ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_chr1_22_temp.bed
mv ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_chr1_22_temp.bed ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed


for two_tech in "Illumina_ONT"
do
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

        bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/orthogonal_tech_beds/${sample}_MajorityCallable_Illumina+ONT_concat_sort_merge_minusChrUn_random.bed
        merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_${version}_${two_tech}_somaticOnly_merged.vcf.gz
        truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_${version}_${two_tech}_somaticOnly_intersection.vcf.gz
        output_prefix=deepsomatic_${version}_${two_tech}

        output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
        date=05022025

        other_regions=""
        name_other_regions=""

        /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
        -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
        -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
    done
done

grep -v "chrM" ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs.bed > ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_temp.bed
mv ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_temp.bed ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs.bed

grep -v "chrM" ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed > ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_chr1_22_temp.bed
mv ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_chr1_22_temp.bed ${date}_HG008_deepsomatic_${version}_${two_tech}_highconf_minusSVs_SDs_chr1_22.bed


for two_tech in "PacBio_ONT"
do
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

        bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/${sample}_working/orthogonal_tech_beds/${sample}_MajorityCallable_HiFi+ONT_concat_sort_merge_minusChrUn_random.bed
        merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_${version}_${two_tech}_somaticOnly_merged.vcf.gz
        truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/${sample}_deepsomatic_${version}_${two_tech}_somaticOnly_intersection.vcf.gz
        output_prefix=deepsomatic_${version}_${two_tech}

        output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
        date=05022025

        other_regions=""
        name_other_regions=""

        /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
        -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
        -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
    done
done

# new Illumina tumor version
sample=UPN237_new_tumor
for two_tech in "Illumina_PacBio" "Illumina_ONT" "PacBio_ONT"
do
    sv_bed1=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_ont_60x.bp1000.bed
    sv_bed2=/private/groups/patenlab/jimin/data/BED/SVs/severus_all_UPN237_hifi.bp1000.bed

    bed_file=/private/groups/patenlab/jimin/data/high_conf_regions_working/UPN237_working/UPN237_MajorityCallable_${two_tech}_concat_sort_merge_minusChrUn_random.bed
    merged_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/UPN237_new_tumor_deepsomatic_${version}_${two_tech}_somaticOnly_merged.vcf.gz
    truth_vcf=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model/${two_tech}/UPN237_new_tumor_deepsomatic_${version}_${two_tech}_somaticOnly_intersection.vcf.gz
    output_prefix=deepsomatic_${version}_${two_tech}

    output_directory=/private/groups/patenlab/jimin/data/BED/DeepSomatic_v1.8.0_fixed_illumina_multicancer
    date=05022025

    other_regions=""
    name_other_regions=""

    /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_bed.sh \
    -b ${bed_file} -m ${merged_vcf} -t ${truth_vcf} -s ${sample} -p ${output_prefix} -o ${output_directory} -d ${date} \
    -j ${sv_bed1} -k ${sv_bed2} -x "" -y ""
done

