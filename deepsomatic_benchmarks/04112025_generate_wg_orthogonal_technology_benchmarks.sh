################################
# training cell lines
################################
1437 1937 2009
# 1954 is re-running
for sample in 1954
do
  illumina_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/04072025_multicancer_wo_spp_${sample}_Illumina_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/04072025_multicancer_wo_spp_${sample}_Illumina_wg.somatic_only.vcf.gz
  hifi_vcf=/private/groups/patenlab/jimin/data/VCF/0701_final_docker/PacBio/0701_${sample}_DeepSomatic_PacBio_wg_somaticOnly.vcf.gz
  ont_vcf=/private/groups/patenlab/jimin/data/VCF/0701_final_docker/ONT/0701_${sample}_DeepSomatic_ONT_wg_somaticOnly.vcf.gz

  output_directory=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model
  variant_caller=deepsomatic
  version=v1.8.0
  filter="orthogonal_technology"

  /private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
  -i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}
done

################
# 578
################
sample=578

illumina_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/04072025_multicancer_wo_spp_${sample}_Illumina_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/04072025_multicancer_wo_spp_${sample}_Illumina_wg.somatic_only.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/PacBio/1107_578_DeepSomatic_v1.7.0_PacBio_wg_somaticOnly.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/ONT/1017_578_DeepSomatic_v1.7.0_ONT_wg_somaticOnly.vcf.gz

output_directory=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model
variant_caller=deepsomatic
version=v1.8.0
filter="orthogonal_technology"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}

################
# HG008
################
sample=HG008

illumina_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/04072025_multicancer_wo_spp_${sample}_Illumina_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/04072025_multicancer_wo_spp_${sample}_Illumina_wg.somatic_only.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/PacBio/0907_HG008_N-P_DeepSomatic_v1.7.0_PacBio_wg_somaticOnly.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/ONT/0907_HG008_N-P_DeepSomatic_v1.7.0_ONT_wg_somaticOnly.vcf.gz

output_directory=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model
variant_caller=deepsomatic
version=v1.8.0
filter="orthogonal_technology"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}


################
# UPN237
################
sample=UPN237

illumina_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/04072025_multicancer_wo_spp_${sample}_Illumina_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/04072025_multicancer_wo_spp_${sample}_Illumina_wg.somatic_only.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/PacBio/012225_UPN237_DeepSomatic_v1.7.0_PacBio_wg_somaticOnly.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/ONT/020525_UPN237_DeepSomatic_v1.7.0_ONT_60x_wg_somaticOnly.vcf.gz

output_directory=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model
variant_caller=deepsomatic
version=v1.8.0
filter="orthogonal_technology"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}


######################
# UPN237 new tumor
######################
sample=UPN237_new_tumor

illumina_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/04072025_multicancer_wo_spp_UPN237_Illumina_new_tumor_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/04072025_multicancer_wo_spp_UPN237_Illumina_new_tumor_wg.somatic_only.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/PacBio/012225_UPN237_DeepSomatic_v1.7.0_PacBio_wg_somaticOnly.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/data/VCF/DeepSomatic_v1.7.0/ONT/020525_UPN237_DeepSomatic_v1.7.0_ONT_60x_wg_somaticOnly.vcf.gz

output_directory=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model
variant_caller=deepsomatic
version=v1.8.0
filter="orthogonal_technology"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}


################
# COLO829
################
sample=COLO829
illumina_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/04072025_multicancer_wo_spp_COLO829_Illumina_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/04072025_multicancer_wo_spp_COLO829_Illumina_wg.somatic_only.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/DeepSomatic_v1.8.0_COLO829_PacBio_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/DeepSomatic_v1.8.0_COLO829_PacBio_wg.somatic_only.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/DeepSomatic_v1.8.0_COLO829_ONT_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/DeepSomatic_v1.8.0_COLO829_ONT_wg.somatic_only.vcf.gz

output_directory=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model
variant_caller=deepsomatic
version=v1.8.0
filter="orthogonal_technology"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}


# 5/2/2025 update samples 578, HG008 and UPN237 new tumor with v17_rc0_07012024 docker for PacBio and ONT


################
# 578
################
sample=578

illumina_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/04072025_multicancer_wo_spp_${sample}_Illumina_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/04072025_multicancer_wo_spp_${sample}_Illumina_wg.somatic_only.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/DeepSomatic_v17_rc0_07012024_578_PacBio_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/DeepSomatic_v17_rc0_07012024_578_PacBio_wg.somatic_only.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/DeepSomatic_v17_rc0_07012024_578_ONT_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/DeepSomatic_v17_rc0_07012024_578_ONT_wg.somatic_only.vcf.gz

output_directory=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model
variant_caller=deepsomatic
version=v17_rc0_07012024
filter="orthogonal_technology"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}

################
# HG008
################
sample=HG008

illumina_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/04072025_multicancer_wo_spp_${sample}_Illumina_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/04072025_multicancer_wo_spp_${sample}_Illumina_wg.somatic_only.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/DeepSomatic_v17_rc0_07012024_HG008_PacBio_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/DeepSomatic_v17_rc0_07012024_HG008_PacBio_wg.somatic_only.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/DeepSomatic_v17_rc0_07012024_HG008_ONT_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/DeepSomatic_v17_rc0_07012024_HG008_ONT_wg.somatic_only.vcf.gz

output_directory=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model
variant_caller=deepsomatic
version=v17_rc0_07012024
filter="orthogonal_technology"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}

######################
# UPN237 new tumor
######################
sample=UPN237_new_tumor

illumina_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/04072025_multicancer_wo_spp_UPN237_Illumina_new_tumor_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/04072025_multicancer_wo_spp_UPN237_Illumina_new_tumor_wg.somatic_only.vcf.gz
hifi_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/DeepSomatic_v17_rc0_07012024_UPN237_PacBio_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/DeepSomatic_v17_rc0_07012024_UPN237_PacBio_wg.somatic_only.vcf.gz
ont_vcf=/private/groups/patenlab/jimin/GITHUB/run_workflows/workflows/DeepSomatic/DeepSomatic_v17_rc0_07012024_UPN237_ONT_wg/analysis/DeepSomatic_outputs/DeepSomatic.postProcess/DeepSomatic_v17_rc0_07012024_UPN237_ONT_wg.somatic_only.vcf.gz


output_directory=/private/groups/patenlab/jimin/data/VCF/two_technology_truthset/BENCHMARK_fixed_multicancer_wgs_model
variant_caller=deepsomatic
version=v17_rc0_07012024
filter="orthogonal_technology"

/private/groups/patenlab/jimin/GITHUB/DeepSomatic_manuscript/deepsomatic_benchmarks/generate_benchmarks_vcf.sh \
-i ${illumina_vcf} -h ${hifi_vcf} -o ${ont_vcf} -s ${sample} -d ${output_directory} -c ${variant_caller} -v ${version} -f ${filter}










