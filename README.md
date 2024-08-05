# DeepSomatic_manuscript
Collection of scripts used for bam, vcf and bed manipulation described in the DeepSomatic manuscript.

# Convert VCF to BED
Script used to convert variant regions in VCF to BED format.
```
Usage: python3 vcf_to_bed_v4.py -v ${VCF} -i ${VCF}.tbi -t 'deepsomatic' -o ${BED}
```

# Filter for high confidence variants
Script used to filter for somatic variants that meet a certain criteria.
```
Usage: python3 vcf_intersection_complex_v2.py -v ${VCF} -i ${VCF}.tbi -f 'filter4' -o filtered.vcf.gz
```

# Titration Scripts
Scripts used to generate titration bams for DeepSomatic analysis. 

`tumor_purity_titration.sh` and `normal_purity_titration.sh` are almost identical, only difference is the naming convention for output files.

### tumor_purity_titration.sh
This script is an adaptation of `purity_titration.sh`, and only creates tumor purity titrations and allows user to specify which tumor purity percentages to titrate input bams to. 
This script uses samtools with 60 threads.

run locally:
```
Usage: ./tumor_purity_titration.sh \
   -t <tumor_bam> \
   -n <normal_bam> \
   -c <tumor_coverage> \
   -q <normal_coverage> \
   -g <tumor_goal_total_coverage> \
   -p <platform> \
   -s <sample> \
   -o <output_directory> \
   -l <tumor_percent_list>
```
### normal_purity_titration.sh
This script is an adaptation of `purity_titration.sh`, and only creates normal purity titrations and allows user to specify which normal purity percentages to titrate input bams to. 
This script uses samtools with 30 threads.

run locally:
```
Usage: ./normal_purity_titration.sh \
   -t <tumor_bam> \
   -n <normal_bam> \
   -c <tumor_coverage> \
   -q <normal_coverage> \
   -x <normal_goal_total_coverage> \
   -p <platform> \
   -s <sample> \
   -o <output_directory> \
   -l <normal_percent_list>
```
## Split Bams
These scripts were used as a pre-processing step before performing purity titrations, in order to ensure that there would not be duplicate reads in our titration set and evaluation set. 

These scripts with split a bam file into two sub-bams, so that there are no overlapping reads in each of the two split sub-bams. 
Two scripts are identical except for naming conventions specified for "normal" or "tumor" bams. These scripts use samtools with 60 threads and require large disk space and memory due to handling files in SAM format during intermediate steps.

### split_bam_normal.sh
run locally:
```
Usage: ./split_bam_normal.sh \
   -n <normal_bam> \
   -q <normal_coverage> \
   -g <normal_goal_coverage> \ # Desired coverage for one of the two sub-bams. The other sub-bam will be the remainder of coverage from <normal_coverage>.
   -e <normal_evaluation_bam> \ # One of the two sub-bams, if already present, if not already present use "None". Coverage of <normal_evaluation_bam> must match <normal_goal_coverage>
   -p <platform> \
   -s <sample> \
   -o <output_directory>
```

### split_bam_tumor.sh
run locally:
```
Usage: ./split_bam_tumor.sh \
   -n <tumor_bam> \
   -q <tumor_coverage> \
   -g <tumor_goal_coverage> \ # Desired coverage for one of the two sub-bams. The other sub-bam will be the remainder of coverage from <tumor_coverage>.
   -e <tumor_evaluation_bam> \ # One of the two sub-bams, if already present, if not already present use "None". Coverage of <tumor_evaluation_bam> must match <tumor_goal_coverage>
   -p <platform> \
   -s <sample> \
   -o <output_directory>
```
