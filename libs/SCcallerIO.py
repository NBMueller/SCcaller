#!/usr/bin/env python

# Standard libraries
import os 
import re
import time
import gzip
from subprocess import Popen, PIPE
# Additional libraries
from pysam import FastaFile, AlignmentFile  # 0.15.1
# SCcaller internal libraries
from libs.OutLineStruct import MULTIPLEGENOTYPE, NOTENOUGHVARIANTS

# Fixed default values
PHREDSCORE = 33

# ------------------------------------------------------------------------------
# INPUT 
# ------------------------------------------------------------------------------

def parse_fasta(fasta_file):
    """ Parse the reference fasta file into the format [(name,head,tail), ...].
    name refers to the contig, head to the index of the first non "N" base, and 
    tail to the index of the last non "N" base

    :param fasta_file_name:
    :return:
    """

    # check catalog
    catalog_file = '{}.catalog'.format(os.path.splitext(fasta_file)[0])
    if os.path.exists(catalog_file) \
            and os.path.getmtime(catalog_file) > os.path.getmtime(fasta_file):
        with open(catalog_file, 'r') as fp:
            catalog = fp.read().split(';')
        my_result = map(
            lambda x: map(
                lambda y: y if x.split(',').index(y) == 0 else int(y),
                    x.split(',')),
            catalog
        )
    # parse fasta
    else:  
        with open(fasta_file) as fp:
            fasta = fp.read().split('>')

        my_result = []
        for target in fasta:
            if target == '':
                continue

            first_line_idx = target.find('\n')
            first_line = target[:first_line_idx]

            if not ' ' in first_line:
                name = first_line
            else:
                name = first_line[:first_line.index(' ')]

            target = target[first_line_idx + 1:].replace('\n', '')
            target_len = len(target)

            try:
                head = re.search('[^N]', target).start()
            except AttributeError:
                head = target_len
                tail = 0
            else:
                for j in range(target_len -1, head, -1):
                    if target[j] != 'N':
                        break
                tail = j + 1
            my_result.append([name, head, tail])

    return my_result


def load_snp_pos(name, snp_info, vcf_file_name):
    """ Read the position information of the specified name chromosome from the 
        vcf file, and store it as a list
    """
    snp_pos = []
    curr_info = [i for i in snp_info if i[0] == name]

    if vcf_file_name.endswith('.gz'):
        file_stream = gzip.open(vcf_file_name, 'rb')
    else:
        file_stream = open(vcf_file_name, 'r')

    if curr_info:
        with file_stream as f:
            for i in curr_info:
                f.seek(i[1])
                lines = f.read(i[2]).splitlines()
                snp_pos.extend([int(line.split('\t')[1]) for line in lines])
    return snp_pos


def parse_snp_info(my_args):
    """ Parse the vcf file into [[name,head,length],[name,head,length]...] format
    """
    snp_file = my_args.snp_in

    has_shown_info = False
    base_file, file_type = os.path.splitext(snp_file)

    # read vcf catalog
    catalog_file = '{}.catalog'.format(base_file)
    if os.path.exists(catalog_file) \
            and os.path.getmtime(catalog_file) > os.path.getmtime(snp_file):
        with open(catalog_file, 'r') as fp:
            catalog = fp.read().split(';')
        result = map(
            lambda x: map(lambda y: y if x.split(',').index(y) == 0 else int(y),
                x.split(',')),
            catalog
        )
        res = [map(lambda y: y if x.split(',').index(y) == 0 else int(y), x.split(',')) \
            for x in catalog]
        import pdb; pdb.set_trace()
    # parse vcf
    else:
        result = []
        name = ''
        head = 0
        if file_type == '.gz':
            file_stream = gzip.open(snp_file, 'rb')
        else:
            file_stream = open(snp_file, 'r')

        with file_stream as f:
            while True:
                line = f.readline()
                if line == '':
                    result.append([name, head, f.tell() - head])
                    break
                if line[0] == '#' or line == '\n':
                    continue
                
                cols = line.split('\t')
                if name != cols[0]:
                    curr_idx = f.tell() - len(line)
                    if name == '':
                        head = curr_idx
                    else:
                        tail = curr_idx
                        result.append([name, head, tail - head])
                        head = tail
                    name = cols[0]

                if my_args.snp_type == 'hsnp' and not has_shown_info \
                        and len(cols) > 8:
                    tmp_str = '\t'.join(cols[9:])
                    tmp_list = re.findall('0\\|0|0/0|1\\|1|1/1|2/2|2\\|2', tmp_str)
                    if len(tmp_list) > 0:
                        print('\n>>>Please confirm the input VCF only ' \
                            'contains heterozygote loci in bulk!<<<\n')
                        has_shown_info = True

    return result


def data_generator(my_args, name, start, stop, is_bulk):
    if my_args.engine == 'samtools':
        generator = data_generator_samtools
    else:
        generator = data_generator_pysam
    return generator(my_args, name, start, stop, is_bulk)


def data_generator_pysam(my_args, name, start, stop, is_bulk):
    fasta_file = FastaFile(my_args.fasta)
    ref = fasta_file.fetch(name, start, stop)

    my_arg = {'fastafile': fasta_file, 'stepper': 'samtools',
        'adjust_capq_threshold': 50, 'contig': name, 'start': start,
        'stop': stop, 'min_mapping_quality': 0 if is_bulk else 20,
        'min_base_quality': 13,}

    if is_bulk:
        bam_file = AlignmentFile(my_args.bulk, 'rb')
    else:
        bam_file = AlignmentFile(my_args.bam, 'rb')

    read_bases_list = []
    for pileup_column in bam_file.pileup(**my_arg):
        pos = pileup_column.reference_pos

        if pos >= stop:
            break
        if pos < start:
            continue

        read_bases_list = pileup_column.get_query_sequences(mark_matches=True,
            mark_ends=True, add_indels=True)

        read_bases = ''.join(read_bases_list).upper()
        n = pileup_column.get_num_aligned()
        if n == 0:
            read_bases = '*'
            base_q = '*'
            map_q = '*'
        else:
            base_q = ''.join([chr(int(i) + PHREDSCORE) \
                for i in pileup_column.get_query_qualities()])
            map_q = ''.join([chr(int(i) + PHREDSCORE) \
                for i in pileup_column.get_mapping_qualities()])

        yield [name, pos, ref[pos - start], str(n), read_bases, base_q, map_q]

    yield None


def data_generator_mpileup(pileup_file):
    pileup_stream = open(pileup_file, 'r')
    for line in pileup_stream:
        el = line.strip().split('\t')
        base_q = ''.join([chr(int(i) + PHREDSCORE) for i in [40] * int(el[3])])
        yield [el[0], int(el[1]), el[2], el[3], el[4], base_q, base_q]
    yield None


def data_generator_samtools(my_args, name, start, stop, is_bulk):
    if is_bulk:
        cmd_str = 'samtools mpileup -C50  -f {} -s {} -r {}:{}-{}' \
            .format(my_args.fasta, my_args.bam, name, start, stop)
    else:
        cmd_str = 'samtools mpileup -C50  -f {} -q 40 -s {} -r {}:{}-{}' \
            .format(my_args.fasta, my_args.bam, name, start, stop)
    process_samtools = Popen([cmd_str], shell=True, stdout=PIPE)
    while 1:
        pileup = []
        if not read_line_from_process(process_samtools, '#', '\t', 7, pileup):
            break
        else:
            yield pileup
    yield None


def read_line_from_process(process_handle, ch, spliter, columns, data_out):
    """ Read data from the process's stdout, filter out the comments, split by spliter, and check the format
    :param process_handle: handle of the process
    :param ch: comment character
    :param spliter:
    :param columns: input. Greater than 0: column number of the data. Others, Unlimited
    :param data_out: output data (should clear buffer before using)
    :return: True, Success. Others, processes have been terminated, and stdout can't read the data.
    """
    while 1:
        buf = []
        read_flag = read_line_from_file(process_handle.stdout, ch, spliter, 
            columns, buf)
        if not read_flag:

            # samtools has been terminated
            if not process_handle.poll() is None:
                return False
            else:
                continue
        else:
            buf[4] = buf[4].upper()
            data_out.extend(buf)
            return True


def read_line_from_file(from_file, ch, spliter, columns, data_out):
    """
    Read data from file, filter out comments, and split by spliter, check format
    :param from_file: file pointer
    :param ch: comment character
    :param spliter:
    :param columns: input. Greater than 0: column number of the data. Others, Unlimited
    :param data_out: output data (should clear buffer before using)
    :return: True, Success. Others, end of file.
    """
    while 1:
        buf = from_file.readline()
        if len(buf) == 0:
            return False
        if buf[0] == ch:
            continue
        buf = buf.strip('\n').split(spliter)
        if columns > 0:
            if len(buf) != columns:
                continue
            else:
                break
        else:
            break
    buf[1] = int(buf[1])
    data_out.extend(buf)
    return True



# ------------------------------------------------------------------------------
# OUTPUT 
# ------------------------------------------------------------------------------

def write_vcf(my_args, version='2.0.0_NB'):
    """ Write the vcf file
    """
    head_str = '##fileformat=VCFv4.1\n' \
        '##fileDate={}\n' \
        '##source=SCcaller_v{}\n' \
        '##reference=file:{}\n' \
        '##INFO=<ID=NS,Number=1,Type=Integer,' \
            'Description="Number of Samples With Data">\n' \
        '##FILTER=<ID={},Description="Multiple genotypes">\n' \
        '##FILTER=<ID={},' \
            'Description="Number of variant supporting reads smaller than {}">\n' \
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' \
        '##FORMAT=<ID=AD,Number=.,Type=Integer,' \
            'Description="Allelic depths for the ref and alt alleles">\n' \
        '##FORMAT=<ID=BI,Number=1,Type=Float,Description="Amplification Bias">\n' \
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n' \
        '##FORMAT=<ID=PL,Number=G,Type=Integer,' \
            'Description="Phred scaled likelihood of: REF/REF, REF/ALT, ALT/ALT">\n' \
        .format(time.strftime('%Y:%m:%d-%H:%M:%S', time.localtime()), version,
            my_args.fasta, MULTIPLEGENOTYPE, NOTENOUGHVARIANTS, my_args.minvar)

    if my_args.bulk != '':
        head_str += '##FORMAT=<ID=SO,Number=1,Type=String,' \
            'Description="Whether it is a somatic mutation.">\n' \
            '##FORMAT=<ID=BN,Number=.,Type=String,' \
            'Description="Bulk Normal information for position.">\n'

    eta_file = '{}.eta'.format(my_args.output)
    with open(eta_file, 'r') as f:
        etas = f.read()
    head_str += etas
    os.remove(eta_file)

    if my_args.pileup:
        sc_name = os.path.basename(my_args.pileup).split('.')[0]
    else:
        sc_name = os.path.basename(my_args.bam).split('.')[0]
    head_str += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n' \
        .format(sc_name)

    body_file = '{}.body'.format(my_args.output)
    with open(body_file, 'r') as f:
        body_str = f.read()
    os.remove(body_file)

    with open(my_args.output, 'w') as fp_out:
        if my_args.format == 'vcf':
            fp_out.write(head_str)
        fp_out.write(body_str)


if __name__ == '__main__':
    print('Here be dragons...')