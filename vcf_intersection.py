import pysam
import argparse

"""
Take in union vcf's and filter out for different filtering criteria:
    1. filter1: Illumina must be present + one other technology
    2. filter2: At least 2/3 present for SNPs + 3/3 present for indels
    3. filter3: At least 2/3 present for all records
    4. filter4: At least 2/3 present for SNPs + Illumina and one other technology present for indels
    5. filter4_orthogonal: At least 2/3 present for SNPs + Illumina and one other technology present for indels applied to
    Strelka2 and ClairS union set

Input VCF must be an UNION VCF with samples merged in specific order: ILLUMINA, HIFI, ONT

Running command:
python3 vcf_intersection.py -v $UNION_VCF -i $UNION_VCF_IDX -f <filter> -o <output file path/name>

"""
parser = argparse.ArgumentParser()
parser.add_argument('--vcf', '-v', type=str, required=True, help='vcf file path')
parser.add_argument('--vcfindex', '-i', type=str, help='vcf index file path')
parser.add_argument('--filter', '-f', type=str, required=True, help='choose filtering criteria')
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
        """
        Filter1: Illumina must be present + one other technology
        """
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
        """
        Filter2: At least 2/3 technologies present for SNPs + 3/3 technologies present for indels
        """
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
        """
        At least 2/3 technologies present for all records
        """
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
        """
        Filter4: (At least 2/3 technologies present for SNVs) + (Illumina + at least 1 other technology present for INDELs)
        """
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

    def filter4_orthogonal(self, records):
        """
        Apply "filter4" to a union set of Strelka2 and ClairS variants.
        """
        for record in records:
            sample_info = [y.values() for y in record.samples.itervalues()]
            strelka2_snv_normal = sample_info[0]
            strelka2_snv_tumor = sample_info[1]
            strelka2_indel_normal = sample_info[2]
            strelka2_indel_tumor = sample_info[3]
            clairs_hifi_snv = sample_info[4]
            clairs_hifi_indel = sample_info[5]
            clairs_ont_snv = sample_info[6]
            clairs_ont_indel = sample_info[7]

            ref = record.ref
            alt = record.alts[0]

            illumina_snv_bool = False
            illumina_indel_bool = False
            hifi_snv_bool = False
            hifi_indel_bool = False
            ont_snv_bool = False
            ont_indel_bool = False

            # check that not all info entries are None
            if len(list(i for i in strelka2_snv_tumor if isinstance(i, int))) > 0:
                illumina_snv_bool = True
            if len(list(i for i in strelka2_indel_tumor if isinstance(i, int))) > 0:
                illumina_indel_bool = True
            if len(list(i for i in clairs_hifi_snv if isinstance(i, int))) > 0:
                hifi_snv_bool = True
            if len(list(i for i in clairs_hifi_indel if isinstance(i, int))) > 0:
                hifi_indel_bool = True
            if len(list(i for i in clairs_ont_snv if isinstance(i, int))) > 0:
                ont_snv_bool = True
            if len(list(i for i in clairs_ont_indel if isinstance(i, int))) > 0:
                ont_indel_bool = True

            is_indel = False

            if len(ref) > 1 or len(alt) > 1:
                is_indel = True

            # Illumina + 1 other technology needs to be called for INDELs
            if is_indel:
                indel_bool_list = [illumina_indel_bool, hifi_indel_bool, ont_indel_bool]
                indel_count = indel_bool_list.count(True)

                if illumina_indel_bool and (hifi_indel_bool or ont_indel_bool):
                    yield record

            # 2/3 technologies must be called for SNPs
            elif not is_indel:
                bool_list = [illumina_snv_bool, hifi_snv_bool, ont_snv_bool]
                sample_count = bool_list.count(True)

                if sample_count >= 2:
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
            
    elif filter == 'filter4_orthogonal':
        for record in filter_criteria4_updated:
            out_file.write(record)

    elif filter == 'minvaf':
        for record in filter_minvaf:
            out_file.write(record)

