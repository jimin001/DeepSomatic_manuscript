import pysam
import argparse

"""
Take in union vcf's and filter out for different criteria:
    1. Illumina must be present + one other technology
    2. At least 2/3 present for SNPs + 3/3 present for indels
    3. At least 2/3 present for all records

Input VCF must be an UNION VCF with samples merged in specific order: ILLUMINA, HIFI, ONT

Running command:
python3 vcf_intersection_complex.py -v $UNION_VCF -i $VCF_IDX -f <filter> -o <output file path/name>

"""
parser = argparse.ArgumentParser()
parser.add_argument('--vcf', '-v', type=str, required=True, help='vcf file path')
parser.add_argument('--vcfindex', '-i', type=str, help='vcf index file path')
parser.add_argument('--filter', '-f', type=str, required=True, help='which filter, options: [filter1 = Illumina must be present + one other technology, \
filter2 = 2/3 present for SNPs + 3/3 present for indels, filter3, filter4, filter4_updated, minvaf]')
parser.add_argument('--minvaf', '-m', type=float,
                    help='minimum vaf, any records with vaf less than this value will be filtered out')

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
            elif attribute == 'format':
                yield record.format.keys()
            elif attribute == 'sample':
                yield [y.values() for y in record.samples.itervalues()]

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

    def filter1(self, records):
        for record in records:
            sample_info = [y.values() for y in record.samples.itervalues()]
            illumina_sample = sample_info[0]
            hifi_sample = sample_info[1]
            ont_sample = sample_info[2]

            illumina_bool = False
            hifi_bool = False
            ont_bool = False

            if None not in illumina_sample:
                illumina_bool = True
            if None not in hifi_sample:
                hifi_bool = True
            if None not in ont_sample:
                ont_bool = True

            if illumina_bool and (hifi_bool or ont_bool):
                yield record

    def filter2(self, records):
        for record in records:
            sample_info = [y.values() for y in record.samples.itervalues()]
            ref = record.ref
            alt = record.alts[0]
            illumina_sample = sample_info[0]
            hifi_sample = sample_info[1]
            ont_sample = sample_info[2]

            illumina_bool = False
            hifi_bool = False
            ont_bool = False

            if None not in illumina_sample:
                illumina_bool = True
            if None not in hifi_sample:
                hifi_bool = True
            if None not in ont_sample:
                ont_bool = True

            is_indel = False
            if len(ref) > 1 or len(alt) > 1:
                is_indel = True

            bool_list = [illumina_bool, hifi_bool, ont_bool]
            sample_count = bool_list.count(True)

            if is_indel and sample_count == 3:
                yield record

            elif not is_indel and sample_count >= 2:
                yield record

    def filter3(self, records):
        for record in records:
            sample_info = [y.values() for y in record.samples.itervalues()]

            illumina_sample = sample_info[0]
            hifi_sample = sample_info[1]
            ont_sample = sample_info[2]

            illumina_bool = False
            hifi_bool = False
            ont_bool = False

            if None not in illumina_sample:
                illumina_bool = True
            if None not in hifi_sample:
                hifi_bool = True
            if None not in ont_sample:
                ont_bool = True

            bool_list = [illumina_bool, hifi_bool, ont_bool]
            sample_count = bool_list.count(True)

            if sample_count >= 2:
                yield record

    def filter4(self, records):
        for record in records:
            sample_info = [y.values() for y in record.samples.itervalues()]
            illumina_sample = sample_info[0]
            hifi_sample = sample_info[1]
            ont_sample = sample_info[2]
            ref = record.ref
            alt = record.alts[0]

            illumina_bool = False
            hifi_bool = False
            ont_bool = False

            if None not in illumina_sample:
                illumina_bool = True
            if None not in hifi_sample:
                hifi_bool = True
            if None not in ont_sample:
                ont_bool = True

            is_indel = False
            if len(ref) > 1 or len(alt) > 1:
                is_indel = True

            if is_indel:
                if illumina_bool and (hifi_bool or ont_bool):
                    yield record

            elif not is_indel:
                bool_list = [illumina_bool, hifi_bool, ont_bool]
                sample_count = bool_list.count(True)

                if sample_count >= 2:
                    yield record

    def filter4_updated(self, records):
        for record in records:
            sample_info = [y.values() for y in record.samples.itervalues()]
            illumina_som_sample = sample_info[0]
            illumina_germ_sample = sample_info[1]

            hifi_som_sample = sample_info[2]
            hifi_germ_sample = sample_info[3]

            ont_som_sample = sample_info[4]
            ont_germ_sample = sample_info[5]

            ref = record.ref
            alt = record.alts[0]

            illumina_som_bool = False
            hifi_som_bool = False
            ont_som_bool = False

            if None not in illumina_som_sample:
                illumina_som_bool = True
            if None not in hifi_som_sample:
                hifi_som_bool = True
            if None not in ont_som_sample:
                ont_som_bool = True

            illumina_germ_bool = False
            hifi_germ_bool = False
            ont_germ_bool = False

            if None not in illumina_germ_sample:
                illumina_germ_bool = True
            if None not in hifi_germ_sample:
                hifi_germ_bool = True
            if None not in ont_germ_sample:
                ont_germ_bool = True

            is_indel = False
            if len(ref) > 1 or len(alt) > 1:
                is_indel = True

            if is_indel:
                if illumina_som_bool and (hifi_som_bool or ont_som_bool):
                    yield record
                if illumina_germ_bool and (hifi_germ_bool or ont_germ_bool):
                    yield record

            elif not is_indel:
                som_bool_list = [illumina_som_bool, hifi_som_bool, ont_som_bool]
                som_sample_count = som_bool_list.count(True)
                if som_sample_count >= 2:
                    yield record

                germ_bool_list = [illumina_germ_bool, hifi_germ_bool, ont_germ_bool]
                germ_sample_count = germ_bool_list.count(True)
                if germ_sample_count >= 2:
                    yield record
                
    def filter_min_vaf(self, records, min_vaf):
        for record in records:
            sample_info = [y.values() for y in record.samples.itervalues()]
            sample_obj = DeepSomaticSampleData(sample_info)
            vaf = sample_obj.VAF[0]
            if vaf >= min_vaf:
                yield record

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
    filter = args.filter
    min_vaf = args.minvaf

    vcf_h = VCFHandler(vcf_file, vcf_index)
    out_file = pysam.VariantFile(out_file, 'w', header=vcf_h.header)

    records = vcf_h.get_records()

    filter_criteria1 = vcf_h.filter1(records)
    filter_criteria2 = vcf_h.filter2(records)
    filter_criteria3 = vcf_h.filter3(records)
    filter_criteria4 = vcf_h.filter4(records)
    filter_criteria4_updated = vcf_h.filter4_updated(records)
    
    filter_minvaf = vcf_h.filter_min_vaf(records, min_vaf)

    if filter == 'filter1':
        for record in filter_criteria1:
            out_file.write(record)

    elif filter == 'filter2':
        for record in filter_criteria2:
            out_file.write(record)

    elif filter == 'filter3':
        for record in filter_criteria3:
            out_file.write(record)

    elif filter == 'filter4':
        for record in filter_criteria4:
            out_file.write(record)
            
    elif filter == 'filter4_updated':
        for record in filter_criteria4_updated:
            out_file.write(record)

    elif filter == 'minvaf':
        for record in filter_minvaf:
            out_file.write(record)

