while getopts n:q:g:e:p:s:o: flag
do
    case "${flag}" in
        n) tumor_bam=${OPTARG};;
        q) tumor_coverage=${OPTARG};;
        g) tumor_goal_coverage=${OPTARG};;
        e) tumor_evaluation_bam=${OPTARG};;
        p) platform=${OPTARG};;
        s) sample=${OPTARG};;
        o) output_directory=${OPTARG};;
    esac
done

echo "tumor bam: $tumor_bam";
echo "Coverage of input tumor bam: $tumor_coverage";
echo "Goal coverage to subset tumor bam by: $tumor_goal_coverage";
echo "One part of split_bam if a partial bam is already present, else put None (string) for this flag: $tumor_evaluation_bam"
echo "Sequencing platform: $platform";
echo "Cell line sample: $sample";
echo "Output directory: $output_directory";

echo "---------------------------------------------------------------------"
echo "Splitting tumor bam into: evaluation.bam (${tumor_goal_coverage}x), spikein.bam (total_tumor_coverage - ${tumor_goal_coverage}x)"
echo "-------------------------------------"

set -o pipefail
set -e
set -u

############################
if [[ ${tumor_evaluation_bam} == "None" ]]
then
        echo "1. Downsampling tumor bam to ${tumor_goal_coverage}x coverage..."
        tumor_fraction=$( bc <<< "scale=3;$tumor_goal_coverage/$tumor_coverage")
        echo "tumor fraction: $tumor_fraction"

        samtools view -@60 -b -s ${tumor_fraction} ${tumor_bam} > ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam
        samtools index -@60 ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam

        echo "Checking coverage with samtools depth for ${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam : "
        samtools depth ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
        tumor_evaluation_bam=${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x.bam

else
        echo "Use existing evalution bam"
        echo "1. Use ${tumor_goal_coverage}x downsampled tumor.bam"

        echo ${tumor_evaluation_bam}

        ############################

        echo "Checking coverage with samtools depth for ${tumor_evaluation_bam} : "
        samtools depth ${tumor_evaluation_bam} | awk '{sum+=$3} END { print "Average = ",sum/NR}'
fi
############################

echo "2. Extracting readnames of full coverage bam.."
time samtools view -@60 ${tumor_bam} | cut -f1 > ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_coverage}x_READNAMES.txt

echo "Check number of lines in: ${sample}_${platform}_Tumor.GRCh38.${tumor_coverage}x_READNAMES.txt "
wc -l ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_coverage}x_READNAMES.txt

echo "Extracting readnames of ${tumor_goal_coverage}x.bam..."
time samtools view -@60 ${tumor_evaluation_bam} | cut -f1 > ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x_READNAMES.txt

echo "Check number of lines in: ${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x_READNAMES.txt"
wc -l ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x_READNAMES.txt

############################

echo "3. Obtaining readnames NOT in evaluation set..."
time grep -Fxvf ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x_READNAMES.txt ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_coverage}x_READNAMES.txt \
> ${output_directory}/${sample}_${platform}_Tumor.GRCh38.SPIKEIN_READNAMES.txt

############################
echo "4. Creating SPIKEIN.bam..."
spikein_coverage=$(( $tumor_coverage - $tumor_goal_coverage ))

samtools view -H ${tumor_bam} > ${output_directory}/temp_header.txt
time samtools view -@60 ${tumor_bam} | grep -Fwf ${output_directory}/${sample}_${platform}_Tumor.GRCh38.SPIKEIN_READNAMES.txt \
> ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${spikein_coverage}x_SPIKEIN.sam

cat ${output_directory}/temp_header.txt ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${spikein_coverage}x_SPIKEIN.sam > ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${spikein_coverage}x_SPIKEIN_withheader.sam
samtools view -b -@60 ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${spikein_coverage}x_SPIKEIN_withheader.sam > ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${spikein_coverage}x_SPIKEIN.bam

# remove intermediate sam files
rm ${output_directory}/*.sam

# index bam file
samtools index -@60 ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${spikein_coverage}x_SPIKEIN.bam

echo "Checking coverage with samtools depth for ${sample}_${platform}_Tumor.GRCh38.${spikein_coverage}x_SPIKEIN.bam"
samtools depth ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${spikein_coverage}x_SPIKEIN.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'

############################
echo "5. Compressing text files to save space..."
pigz ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_coverage}x_READNAMES.txt
pigz ${output_directory}/${sample}_${platform}_Tumor.GRCh38.${tumor_goal_coverage}x_READNAMES.txt
pigz ${output_directory}/${sample}_${platform}_Tumor.GRCh38.SPIKEIN_READNAMES.txt
