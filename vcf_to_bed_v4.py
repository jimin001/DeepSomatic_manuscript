import pysam
import argparse

"""
running command:

VCF=/Users/jimin/CGL/nanopore_somatic/VCF/1205_1937_SomaticOnly_DeepSomatic_Illumina_PacBio_ONT_merged.vcf.gz
VCF_IDX=/Users/jimin/CGL/nanopore_somatic/VCF/1205_1937_SomaticOnly_DeepSomatic_Illumina_PacBio_ONT_merged.vcf.gz.tbi
OUTFILE=/Users/jimin/CGL/nanopore_somatic/VCF/1205_1937_SomaticOnly_DeepSomatic_Illumina_PacBio_ONT_merged.bed

python3 vcf_to_bed_script.py -v $VCF -i $VCF_IDX -o $OUTFILE
"""
parser = argparse.ArgumentParser()
parser.add_argument('--vcf', '-v', type=str, required=True, help='vcf file path')
parser.add_argument('--vcfindex', '-i', type=str, help='vcf index file path')
parser.add_argument('--vcftype', '-t', type=str, help='vcf type: deepsomatic or severus')
parser.add_argument('--basepairs', '-b', type=int, help='SV only: number of basepairs to extract upstream and downstream of SV')
parser.add_argument('--outfile', '-o', type=str, help='output file path')

args = parser.parse_args()

class DeepSomaticSampleData:
    """
    Store information for GT:AP:GQ:DP:AD:VAF:REP.
    """
    def __init__(self, samples):
        self.samples = samples

        for sample in samples:
            self.sample = sample
            self.GT = self.sample[0]
            self.GQ = self.sample[1]
            self.DP = self.sample[2]
            self.AD = self.sample[3]
            self.VAF = self.sample[4]
            self.PL = self.sample[5]

class VCFHandler:
    """
    VCF Handler class.
    """
    def __init__(self, input_vcf_path, input_vcf_index_path):
        self.vcf = pysam.VariantFile(input_vcf_path, 'r', index_filename=input_vcf_index_path)
        self.header = self.vcf.header
        self.contig = None

    def get_header(self):
        """
        Get header of VCF.
        """
        return self.header

    def get_records(self):
        """
        Read in vcf line by line in specified contig region
        """
        # returns iterator object
        return self.vcf.fetch(self.contig)

    def get_record_attributes(self, records, attribute):
        """
        Get specified feature of vcf record.
        """
        # returns generator object of specified attribute
        for record in records:
            if attribute == 'chrom':
                yield record.chrom
            elif attribute == 'pos':
                yield record.pos
            elif attribute == 'id':
                yield record.id
            elif attribute == 'ref':
                yield record.ref
            elif attribute == 'alt':
                yield record.alts
            elif attribute == 'qual':
                yield record.qual
            elif attribute == 'filter':
                yield record.filter.keys()
            elif attribute == 'info':
                yield record.info.keys()
            # elif attribute == 'format':
                # yield record.format.keys()
            # elif attribute == 'sample':
                # yield [y.values() for y in record.samples.itervalues()]

    def get_record_samples(self, samples_list, type):
        """
        Return generator object of specified data type of sample.
        :param samples_list: samples attribute generator object or list
        :param type: GT:AP:GQ:DP:AD:VAF:REP
        :return:
        """
        for s in samples_list:
            samples = DeepSomaticSampleData(s)
            if type == 'GT':
                yield samples.GT
            elif type == 'GQ':
                yield samples.GQ
            elif type == 'DP':
                yield samples.DP
            elif type == 'AD':
                yield samples.AD
            elif type == 'VAF':
                yield samples.VAF
            elif type == 'PL':
                yield samples.PL
            else:
                return 'Incorrect sample type'

    def vcf_to_bed(self, records):
        for record in records:
            chrom = record.chrom
            pos = record.pos
            ref = record.ref

            pos1 = int(pos) - 1
            pos2 = int(pos) + len(ref) - 1
            reformat_string = chrom + '\t' + str(pos1) + '\t' + str(pos2) + '\n'

            yield reformat_string

    def vcf_to_bed_sv(self, records, bp):
        for record in records:
            chrom = record.chrom
            pos = record.pos
            info = record.info.values()
            sv_type = record.info['SVTYPE']
            sv_length = record.info['SVLEN']

            if sv_type == "DEL" or sv_type == "INS" or sv_type == "DUP":
                start_pos1 = int(pos) - bp
                start_pos2 = int(pos) + bp

                end_pos1 = int(pos) + int(sv_length) - bp
                end_pos2 = int(pos) + int(sv_length) + bp

                start_bed_string = chrom + '\t' + str(start_pos1) + '\t' + str(start_pos2) + '\n'
                end_bed_string = chrom + '\t' + str(end_pos1) + '\t' + str(end_pos2) + '\n'
                yield_string = start_bed_string + end_bed_string

            elif sv_type == "INV" or sv_type == "BND":
                start_pos1 = int(pos) - bp
                start_pos2 = int(pos) + bp

                yield_string = chrom + '\t' + str(start_pos1) + '\t' + str(start_pos2) + '\n'

            yield yield_string

if __name__ == '__main__':
    vcf_file = args.vcf
    vcf_index = args.vcfindex
    out_file = args.outfile
    vcf_type = args.vcftype
    basepairs = args.basepairs
    vcf_h = VCFHandler(vcf_file, vcf_index)

    records = vcf_h.get_records()

    if vcf_type == 'deepsomatic':
        bed_strings = vcf_h.vcf_to_bed(records)
    elif vcf_type == 'severus':
        bed_strings = vcf_h.vcf_to_bed_sv(records, bp=basepairs)
    else:
        print('enter correct vcf caller type')

    out = open(out_file, 'w')

    for region in bed_strings:
        out.write(region)
