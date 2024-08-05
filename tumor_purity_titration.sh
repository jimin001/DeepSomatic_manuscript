while getopts t:n:c:q:g:p:s:o:l: flag
do
    case "${flag}" in
        t) tumor_bam=${OPTARG};;
        n) normal_bam=${OPTARG};;
        c) tumor_coverage=${OPTARG};;
        q) normal_coverage=${OPTARG};;
        g) tumor_goal_total_coverage=${OPTARG};;
        p) platform=${OPTARG};;
        s) sample=${OPTARG};;
        o) output_directory=${OPTARG};;
        l) tumor_percent_list=${OPTARG};;
    esac
done

echo "Tumor bam: $tumor_bam";
echo "Normal bam: $normal_bam";
echo "Coverage of input tumor bam: $tumor_coverage";
echo "Coverage of input normal bam: $normal_coverage";
echo "Goal total coverage for tumor bams: $tumor_goal_total_coverage";
echo "Sequencing platform: $platform";
echo "Cell line sample: $sample";
echo "Output directory: $output_directory";
echo "List of tumor purity percentages to titrate bams with, as decimal: $tumor_percent_list"

# Set the exit code of a pipeline to that of the rightmost command
# to exit with a non-zero status, or zero if all commands of the pipeline exit
set -o pipefail
# cause a bash script to exit immediately when a command fails
set -e
# cause the bash shell to treat unset variables as an error and exit immediately
set -u

echo "---------------------------------------------------------------------"
echo "Creating purity titration bams at : $tumor_percent_list tumor purity"
echo "-------------------------------------"
echo ""
for tumor_percent in $tumor_percent_list
do
    echo "tumor percent: $tumor_percent"

    normal_percent=$( bc <<< "scale=2; 1-$tumor_percent" )
    echo "normal percent: $normal_percent"

    tumor_goal_coverage_temp=$( bc <<< "scale=1;$tumor_goal_total_coverage * $tumor_percent" )
    normal_goal_coverage_temp=$( bc <<< "scale=1;$tumor_goal_total_coverage * $normal_percent" )

    tumor_goal_coverage=${tumor_goal_coverage_temp%.*}
    normal_goal_coverage=${normal_goal_coverage_temp%.*}

    echo "tumor goal coverage: $tumor_goal_coverage"
    echo "normal goal coverage: $normal_goal_coverage"

    tumor_fraction=$( bc <<< "scale=3;$tumor_goal_coverage/$tumor_coverage")
    echo "tumor fraction: $tumor_fraction"

    normal_fraction=$( bc <<< "scale=3;$normal_goal_coverage/$normal_coverage + 42")
    echo "normal fraction: $normal_fraction"
    echo "-------------------------------------"

    echo "Downsampling tumor bam to ${tumor_goal_coverage}x coverage..."
    time samtools view -@60 -b -s ${tumor_fraction} ${tumor_bam} > ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam
    time samtools index -@60 ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam

    echo "Downsampling normal bam to ${normal_goal_coverage}x coverage..."
    time samtools view -@60 -b -s ${normal_fraction} ${normal_bam} > ${output_directory}/${sample}_${platform}_Normal.GRCh38.${normal_goal_coverage}x.bam
    time samtools index -@60 ${output_directory}/${sample}_${platform}_Normal.GRCh38.${normal_goal_coverage}x.bam

    echo "Mergeing tumor and normal bams..."
    time samtools merge -@60 ${output_directory}/${sample}_${platform}.GRCh38.${tumor_goal_coverage}xT_${normal_goal_coverage}xN.bam ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam ${output_directory}/${sample}_${platform}_Normal.GRCh38.${normal_goal_coverage}x.bam
    time samtools index -@60 ${output_directory}/${sample}_${platform}.GRCh38.${tumor_goal_coverage}xT_${normal_goal_coverage}xN.bam

    echo "Checking coverage with samtools depth for ${sample}_${platform}.GRCh38.${tumor_goal_coverage}xT_${normal_goal_coverage}xN.bam : "
    time samtools depth ${output_directory}/${sample}_${platform}.GRCh38.${tumor_goal_coverage}xT_${normal_goal_coverage}xN.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
    echo "-------------------------------------"

done
