#!/usr/bin/env python

### Single Cell Caller (SCcaller) - Identify single nucleotide variations (SNV) from single cell sequencing data
# Copyright (C) 2016  Dong, et. al.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU Affero General Public License for more details.

### updates
# v2.0.0, 2019.04.01, allowing parallele processing, optimizing I/O, optimizing pipeline, output in vcf format, and fixing bugs
# v1.21, 2018.08.18, fixing bugs
# v1.2, 2017.05.01, allowing INDEL calling
# v1.1.3, 2017.01.09, users can change the min mapQ, default to 40
# v1.1.2, 2016.12.30, fixing bugs
# v1.1.1, 2016.12.29, updating read_mpileup to consider indels
# v1.1, 2016.07.25, fixing bugs
# v1.0, 2016.04.26, cleanup and release version of v0.0.4
# v0.0.4, 2016.04.26, fixing bugs in readling mpileup file
# v0.0.3, 2016.04.22, read_mpilup function returns mindepth fails before returning reference genotype
# v0.0.3, 2016.04.22, default mapQ change from 20 to 40
# v0.0.2, 2016.04.19, fixing bugs - jump mpileup file column not fit problem.
# v0.0.1, 2016.03, adding likelihood ratio test based on null distribution from the data resampling.

# Standard libraries
import sys
import os
import re
import time
import random
import copy
import argparse
import logging
from subprocess import Popen, PIPE
import multiprocessing as mp
from collections import Counter
# Additional libraries
import pysam  # 0.15.1
import numpy as np
# SCcaller internal libraries
from libs.BigForewordList import BigForewordList
from libs.OutLineStruct import OutLineStruct, MULTIPLEGENOTYPE, NOTENOUGHVARIANTS


if float(sys.version[:3]) != 2.7:
    print('\nWARNING: Python 3.X not tested for vcf .catalog file (for dbsnp) '
        'and bed output (instead of vcf)!')


# regular expression of indels
INSERTION = "\+[0-9]+[ACGTNacgtn]+"
DELETION = "\-[0-9]+[ACGTNacgtn]+"
INDEL = "\+[0-9]+[ACGTNacgtn]+|\-[0-9]+[ACGTNacgtn]+"

# Fixed queue indicators
WORKVAR = 1
WORKCOVERAGE = 2
WORKVCF = 3
WORKDONE = 'done'
WORKERROR = 'error'
WORKEXIT = 'exit'

# Fixed default values
PHREDSCORE = 33
DEFAULTBASEQUALITY = 30
CUTOFFNUM = 100000

LOG_FORMAT = "%(asctime)s %(levelname)-8s %(process)d [%(funcName)s]" \
    "%(message)s(%(filename)s line%(lineno)d)"
VERSION = "2.0.0_NB"


class QueueData:
    def __init__(self, worker_id, work_type, outline=None, coverage=None,
            mode=""):
        self.worker_id = worker_id  # type: int
        self.work_type = work_type

        self.outline = outline
        if coverage == None:
            self.coverage = []
        else:
            self.coverage = coverage
        self.mode = mode


    def __str__(self):
        out_str = 'QueueData:\n\tID: {}\n\tType: {}\n\tMode: {}\n\t' \
                'Outline: {}\n\tCoverage: {}\n\t' \
            .format(self.worker_id, self.work_type, self.mode, self.outline, 
                self.coverage)
        return out_str


class GoldenHetero:
    def __init__(self, name, pos, refcount, altcount):
        self.name = name  # type: str
        self.pos = pos  # type: int
        self.refcount = refcount  # type: int
        self.altcount = altcount  # type: int


    def __str__(self):
        return '{}:{}\tAD:{},{}' \
            .format(self.name, self.pos, self.refcount, self.altcount)


class PileupDataStruct:
    def __init__(self, name, pos, reference_base, variant,
                reference_allele_num, variant_allele_num, genotype_num,
                read_info_list, so=None, genotype_class=-1, variant_all=""):
        # type: (str, int, str, str, int, int, int, list[tuples], str, int) -> object
        self.name = name  # type: str
        self.pos = pos  # type: int
        self.reference_base = reference_base  # type: str
        self.variant = variant  # type: str
        self.reference_allele_num = reference_allele_num  # type: int
        self.variant_allele_num = variant_allele_num  # type: int
        self.genotype_num = genotype_num  # type: int
        self.read_info_list = read_info_list  # type: list[tuples]
        self.so = so  # type: str  # reasoning  bulk ref -- "True"   bulk var -- "False"  else -- "NA"
        self.genotype_class = genotype_class  # type: int  # 0 unknown 1:0,0 0,1 1,1      2:1,1 1,2 2,2
        self.variant_all = variant_all


    def __str__(self):
        return '{}:{}\tref:{}({}),alt:{}({})' \
            .format(self.name, self.pos, self.reference_base,
                self.reference_allele_num, self.variant, self.variant_allele_num)


def compress(l, f):
    return [j for i, j in enumerate(l) if f[i]]


def parse_args():
    parser = argparse.ArgumentParser(
        prog='SCcaller',
        usage='''\npython {}
            [-h] [-d WKDIR] [-l LAMB] [--bias NUM] [--minvar NUM]
            [--mapq NUM] [--min_depth NUM] [--RD NUM] [--null NUM]
            [--bulk BAM] [--bulk_min_depth NUM] [--bulk_min_mapq NUM]
            [--bulk_min_var NUM] [--format {{bed,vcf}}] [--head NUM]
            [--tail NUM] [-e {{pysam, samtools}}] [--cpu_num NUM] [-w NUM]
            [-n NUM] -t {{dbsnp,hsnp}} -b BAM -f FASTA -s
            SNP_IN -o OUTPUT'''.format(__file__),
        description='''SCcaller v{};
            Xiao Dong, biosinodx@gmail.com, xiao.dong@einstein.yu.edu;
            Yujue Wang, spsc83@gmail.com'''.format(VERSION)
    )

    parser.add_argument("-b", "--bam", type=str, required=True, 
        help="Bamfile of a single cell")
    parser.add_argument("-f", "--fasta", type=str, required=True, 
        help="Fasta file of reference genome")
    parser.add_argument("-s", "--snp_in", type=str, required=True,
        help="Candidate snp input file, either from dbsnp data or heterozygous "
            "snp (hsnp) data of the bulk, for known heterogous call. "
            "file type: bed (1-based) or vcf.")
    parser.add_argument("-o", "--output", type=str, required=True, 
        help='Output file name')
    parser.add_argument("-t", "--snp_type", type=str, required=True,
        choices=["dbsnp", "hsnp"],
        help="SNP type for --snp_in. When choosing dbsnp, "
            "--bulk bulk_bamfile is required.")
    parser.add_argument("-d", "--wkdir", type=str, default="./",
        help="Work dir. Default: ./")
    # parser.add_argument("-h", "--help", action="store_true")
    parser.add_argument("-l", "--lamb", type=int, default=10000,
        help="Lambda for bias estimation. Default: 10000.")
    parser.add_argument("--bias", type=float, default=0.75,
        help="Default theta (bias) for SNVs whose theta cannot be estimated. "
            "Default: 0.75.")
    parser.add_argument("--minvar", type=int, default=4,
        help="Min. # variant supporting reads. Default: 4")
    parser.add_argument("-mvf", "--minvarfrac", type=float, default=0.2,
        help='Min. fraction of variant supporting reads for a 0/1 genotype. '
            'Default: 0.2')
    parser.add_argument("--mapq", type=int, default=40,
        help="Min. mapQ. Default: 40")
    parser.add_argument("--min_depth", type=int, default=10,
        help="Min. # reads. Default: 10")
    parser.add_argument("--RD", type=int, default=20,
        help="Min. read depth of known heterogous SNP called from bulk when "
            "choosing -t dbsnp. Default: 20. Recommended: 10,15,20, depending "
            "on average read depth.")
    parser.add_argument("--null", type=float, default=0.03,
        help="Min allelic fraction considered. Default=0.03.")
    parser.add_argument("-e", "--engine", type=str, choices=["pysam", "samtools"],
        default="pysam", help="Pileup engine. Default: pysam")
    parser.add_argument("-w", "--work_num", type=int, default=-1,
        help="Deprecated (defined automatically.")
    parser.add_argument("-n", "--cpu_num", type=int, default=1,
        help="Num. processes. Default: 1")
    parser.add_argument("--head", type=int, default=1,
        help="First chromosome as sorted as in fasta file to analyze (1-based). "
            "Default: 1")
    parser.add_argument("--tail", type=int, default=-1,
        help="Last chromosome as sorted as in fasta file to analyze (1-based). "
        "Default: -1")
    parser.add_argument("--format", type=str, choices=["bed", "vcf"],
        default="vcf", help="Output file format. Default: vcf")
    parser.add_argument("--coverage", action="store_true", default=False,
        help="use \"--coverage\" to generate the coverage file at the same time")
    parser.add_argument("--bulk", type=str, default="",
        help="bamfile of bulk DNA sequencing")
    parser.add_argument("--bulk_min_var", type=int, default=1,
        help="Min. num. variant supporting reads for bulk. Default: 1")
    parser.add_argument("--bulk_min_mapq", type=int, default=20,
        help="Min. mapQ for bulk. Default: 20")
    parser.add_argument("--bulk_min_depth", type=int, default=20,
        help="Min. reads for bulk. Default: 20")
    parser.add_argument("--debug", action='store_true',
        help="Turn of multiprocessing/threading for debugging. Default: False")

    args = parser.parse_args()

    os.chdir(args.wkdir)

    # check mandatory files
    if args.snp_type == "dbsnp" and args.bulk == "":
        raise IOError("{}: When choosing dbsnp, --bulk file is required." \
            .format(__file__))
    for req_file in (args.bam, args.fasta, args.snp_in):
        if not os.path.exists(req_file):
            raise IOError('{}: error: file [{}] does not exist.' \
                .format(__file__, req_file))

    # check result file
    if os.path.exists(args.output): #TAMA
        os.remove(args.output) #TAMA
    for out_ending in (".coverage", ".reasoning"):
        if os.path.exists(get_my_filename(args.output, out_ending)):
            os.remove(get_my_filename(args.output, out_ending))

    return args


def parse_fasta(fasta_file):
    """
    Parse the fasta file into the format [[name,head,tail],[name,head,tail]...]
    :param fasta_file_name:
    :return:
    """

    # check catalog
    catalog_file = "{}.catalog".format(os.path.splitext(fasta_file)[0])
    if os.path.exists(catalog_file) \
            and os.path.getmtime(catalog_file) > os.path.getmtime(fasta_file):
        with open(catalog_file, "r") as fp:
            catalog = fp.read().split(";")
        my_result = map(
            lambda x: map(
                lambda y: y if x.split(",").index(y) == 0 else int(y),
                    x.split(",")),
            catalog
        )
    # parse fasta
    else:  
        with open(fasta_file) as fp:
            fasta = fp.read().split(">")

        my_result = []
        for target in fasta:
            if target == "":
                continue

            first_line_idx = target.find("\n")
            first_line = target[:first_line_idx]

            if not " " in first_line:
                name = first_line
            else:
                name = first_line[:first_line.index(' ')]

            target = target[first_line_idx + 1:].replace("\n", "")
            target_len = len(target)

            try:
                head = re.search('[^N]', target).start() + 1
            except AttributeError:
                head = target_len
                tail = 0
            else:
                for j in range(target_len -1, 0, -1):
                    if target[j] != "N":
                        break
                tail = j + 1 

            my_result.append([name, head, tail])

    return my_result


def parse_indel(indel_str, indel_list_out):
    """
    parse indel string, return the indel string and length of indel string
    :param indel_str:indel string (the '+' and '-' in indel string doesn't affect the result)
    :param indel_list_out:indel string (should be clear before)
    :return:
    >0 return indel string length
    =0 in case of indel and svn
    <0 in case of invalid format of indel
    """
    i = 0
    j = 0
    len_indel_str = len(indel_str)
    for i in range(len_indel_str):
        if indel_str[i] == "+" or indel_str[i] == "-":
            continue
        else:
            break
    # more than 1 '+' or '-'
    if i > 1:
        return -1
    buf = indel_str[i:]
    len_buf = len(buf)
    for j in range(len_buf):
        if buf[j].isalpha():
            break

    if len_indel_str - i - j > int(buf[:j]):
        indel_list_out.extend(indel_str[0:i + j + int(buf[:j])])
        return 0
    elif len_indel_str - i - j < int(buf[:j]):
        return -2
    indel_list_out.extend(indel_str[0:i + j + int(buf[:j])])
    return int(buf[:j]) + j + i


def remove_head_end(pileup):
    # remove '$' which is not head mapQ
    str1 = re.sub("(?<!\^)\$", "", pileup[4])
    str2 = re.sub("\^\S{1}", "", str1)

    missing_num = int(pileup[3]) - len(str2.replace("I", ""))
    if missing_num == 0:
        return str2
    # look for pattern like '^*^'
    total_found_num = 0
    while 1:
        tmp = re.findall("\^\S{1}\^", str1)
        if len(tmp) == 0:
            break
        total_found_num += 1
        # fill in the lost data
        head_index = str1.find(tmp[0])
        str1 = str1[:head_index + 1] + pileup[6][get_head_MQ(str1, head_index)] \
            + str1[head_index + 1:]

    if total_found_num == missing_num:
        return re.sub("\^\S{1}", "", str1)

    num_found = 0
    # look for pattern like '^'
    for i in range(len(str1) - 1):
        MQ_head = pileup[6][get_head_MQ(str1, i)]
        if str1[i] == "^" and str1[i + 1] != MQ_head:
            num_found += 1
            str1 = str1[:i + 1] + MQ_head + str1[i + 1:]

    total_found_num += num_found
    if total_found_num == missing_num:
        return re.sub("\^\S{1}", "", str1)

    logging.critical("Cannot handle this pileup: "
            "name={} pos={} readbase_len={} readbase={} map_q={}" \
        .format(pileup[0], pileup[1], pileup[3], pileup[4], pileup[6])
    )
    return "???"


def get_head_MQ(read_bases, head_index):
    """
    calculate the mapQ of read head (^) from readbase string without I
    :param read_bases:
    :param head_index:
    :return:
    """
    counter = 0
    for i, read_base in enumerate(read_bases):
        if i == head_index:
            break
        if read_base in [".", ",", "A", "T", "G", "C"]:
            if i == 0:
                counter += 1
            elif read_bases[i - 1] != "^":
                counter += 1
    return counter


def compress_read_bases(read_bases, my_filter):
    """
    compress read bases but leave the I (indel)
    :type read_bases: str
    :type my_filter: list
    """
    for i, base in enumerate(read_bases):
        if base == "I":
            my_filter.insert(i, True)
    return "".join(compress(read_bases, my_filter))


def get_reference_variant_allele_num(read_bases):
    v_num = 0
    r_num = 0
    for i, base in enumerate(read_bases):
        if base in [".", ","]:
            r_num += 1
        if base in ["A", "T", "C", "G"]:
            v_num += 1
        if i > 0 and base in ["I", "i"] and read_bases[i - 1] in [".", ","]:
            r_num -= 1
            v_num += 1
    return r_num, v_num


def rebuild_read_base_list(read_base_without_i, indel_name_list, indel_count_list):
    """
    add indel to the end of the readbase without I, return as list
    gather
    :param read_base_without_i:
    :param indel_name_list:
    :param indel_count_list:
    :return:
    """
    read_base_list = []
    read_base_list.extend(list(read_base_without_i))
    tmp = [[indel_name_list[list_index] \
                for i in range(indel_count_list[list_index])]
            for list_index in range(len(indel_count_list))]
    for i in tmp:
        read_base_list.extend(i)
    return read_base_list


def read_mpileup(pileup, rm_minvar, min_mapq, rm_mindepth, is_gh, worker_id):
    """
    screenfor mutations, filter out alignments that do not meet the requirements of mapQ, and convert the data format
    :param is_gh:
    :type is_gh: bool
    :type rm_mindepth: int
    :type min_mapq: int
    :type rm_minvar: int
    :type pileup: list
    :param pileup: 1 line data of mpileup
    :param rm_minvar: variant supporting reads.
    :param min_mapq: minimun mapQ
    :param rm_mindepth: minimun read depth
    :return:
    not empty PileupDataStruct: mutation data
    -1: return reference genome type and for skipping this line, but should be count in coverage
    []: not mutation, should not be count in coverage
    """

    # indel read overlaps
    if "*" in pileup[4]:
        return []  

    if "-" in pileup[4] or "+" in pileup[4]:
        indel_list = re.findall(INDEL, pileup[4])
        if indel_list:
            tmp = Counter(indel_list).most_common()
             # multiple different indel calls
            if len(tmp) > 1:
                return []
            rm_indel = tmp[0][0]
            indel_name_list = [tmp[0][0]]
            indel_count_list = [tmp[0][1]]
            
            result = []
            # indel reads followed by a SNV read, e.g. ....-2AAAA.....
            if parse_indel(rm_indel, result) == 0:
                return []

            pileup[4] = pileup[4].replace(rm_indel, "I")
    else:
        rm_indel = ""
        indel_name_list = []
        indel_count_list = []
        
    # remove head (^) and tail ($)
    pileup[4] = remove_head_end(pileup)

    # filter out alignments that do not meet the requirements of mapQ
    map_q_filter = [ord(i) - PHREDSCORE >= min_mapq for i in pileup[6]]
    read_base_with_i = compress_read_bases(pileup[4], map_q_filter)
    read_base_without_i = re.sub("I", "", read_base_with_i)

    if len(read_base_without_i) < rm_mindepth:
        return []

    if is_gh:
        return read_base_without_i

    var_names = ["A", "C", "G", "T", "I"]
    var_counts = [read_base_with_i.count(element) for element in var_names]

    # return reference genome type and for skipping this line
    if np.max(var_counts) < rm_minvar:
        return -1

    var_max = var_names[np.argmax(var_counts)]
    if var_max == 'I':
        var_max = rm_indel
        if len(rm_indel) == 0:
            logging.critical("({}:{}) has 'I' but no rm_indel" \
                .format(pileup[0], pileup[1]))

    read_base_with_d = "".join([j if map_q_filter[i] else "D" \
        for i, j in enumerate(pileup[4])])
    r_num, v_num = get_reference_variant_allele_num(read_base_with_d)
    
    read_base_list_final = rebuild_read_base_list(read_base_without_i,
        indel_name_list, indel_count_list)
    
    base_q = [ord(j) - PHREDSCORE for i,j in enumerate(pileup[5]) \
        if map_q_filter[i]]
    num_of_i = len(read_base_with_i) - len(read_base_without_i)
    base_q.extend([DEFAULTBASEQUALITY] * num_of_i)

    read_info_list = [(j, base_q[i]) for i, j in enumerate(readbase_list_with_i)]

    return PileupDataStruct(pileup[0], pileup[1], pileup[2], var_max,
        r_num, v_num, 1, read_info_list, 1)


def choose(candidate_index, current_basis, choice_out, basis_filter):
    """

    choose indexs from the candidate_index according to current_basis
    :param candidate_index: the length equal to the number of 'True's in basis_filter
    :type values: list[int]
    """
    values = compress(current_basis, basis_filter)
    tmp = Counter(values).most_common()
    value_list = [i[0] for i in tmp]
    count_list = [i[1] for i in tmp]
    sorted_value_list = copy.copy(value_list)
    sorted_value_list.sort(reverse=True)
    candidate = []
    if len(choice_out) >= 2:
        return candidate

    # Only one of the maximum
    if count_list[value_list.index(sorted_value_list[0])] == 1:
        new_choice_0 = candidate_index[values.index(sorted_value_list[0])]
        choice_out.append(new_choice_0)
        basis_filter[new_choice_0] = False

        if len(choice_out) < 2:
            # only one second largest value
            if count_list[value_list.index(sorted_value_list[1])] == 1:
                new_choice_1 = candidate_index[values.index(sorted_value_list[1])]
                choice_out.append(new_choice_1)
                basis_filter[new_choice_1] = False
            else:
                new_candidates = [candidate_index[i] for i,j in enumerate(values) \
                    if j == sorted_value_list[1]]
                candidate.extend(new_candidates)

    else:
        new_candidates = [candidate_index[i] for i,j in enumerate(values) \
            if j == sorted_value_list[0]]
        candidate.extend(new_candidates)
    for i in range(len(basis_filter)):
        if i not in candidate:
            basis_filter[i] = False
    return candidate


def choose_random(candidate, num):
    result = []
    for i in range(num):
        rdm_i = np.random.randint(len(candidate))
        rdm_candidate = candidate.pop(rdm_i)
        result.append(rdm_candidate)
    return result


def get_variant_info(var_names, basis, read_bases_filtered, geno_class):
    # type: (list[str], list[list[int]], str, str) -> list
    """
    calculate the variant's name_str, ref_num and alt_num
    if variant is more than 2, choose 2 according to basis. If basis doesn't work, choose randomly.
    :param basis: [reads, map_q, base_q]
    :param var_names:
    :param ead_bases_filtered:
    :return: [variant_str, ref_num, alt_num]
    """
    if len(var_names) == 0:
        var_info = ("", 0, 0)
    elif len(var_names) == 1:
        r_num = 0
        v_num = 0
        for read_base in read_bases_filtered:
            if read_base in [".", ","]:
                r_num += 1
            elif read_base in ["A", "T", "C", "G"] or "+" in read_base \
                    or "-" in read_base:
                v_num += 1

        var_info = (var_names[0], r_num, v_num)
    else:
        choice = []
        pos_index = range(len(var_names))
        basis_filter = [True for i in var_names]
        for current_basis in basis:
            pos_index = choose(pos_index, current_basis, choice, basis_filter)
            if len(choice) >= 2:
                break
        if len(choice) < 2:
            choice.extend(choose_random(pos_index, 2 - len(choice)))

        i1, i2 = choice[:2]
        if geno_class == 2:
            var_info = ("{},{}".format(var_names[i1], var_names[i2]),
                basis[0][i1], basis[0][i2])
        else:
            ref_num = read_bases_filtered.count(".") \
                + read_bases_filtered.count(",")
            var_info = (var_names[i1], ref_num, basis[0][i1])
    return var_info


def split_readbase(read_bases, indel_list):
    """
    ignore the ',' and '.' before indel
    ..+1A ---> ['.', '+1A']
    :param readbase_with_i:
    :param indel_list:
    :return:
    """
    tmp_list = []
    index = 0
    for read_base in read_bases:
        if read_base == "I":
            tmp_list[-1] = indel_list[index]
            index += 1
        else:
            tmp_list.append(read_base)
    return tmp_list


def read_mpileup_vcf(pileup, my_args):
    """
    screenfor mutations for vcf format output, filter out alignments that do not meet the requirements of mapQ, and convert the data format
    :type pileup: list
    :param pileup: one line of mpileup file
    :return:
    not empty PileupDataStruct: data
    -1: return reference genome type and for skipping this line, but should be count in coverage
    []: should be skip.
    """

    # indel read overlaps
    if "*" in pileup[4]:
        return []
    
    indel_name_list = []
    indel_count_list = []
    indel_list = []
    if "-" in pileup[4] or "+" in pileup[4]:
        indel_list = re.findall(INDEL, pileup[4])
        if indel_list:
            # check data: eliminate the invalid indel like "+9AA", 
            #   split the indel and the mutation together, 
            #   for example, convert +2AACT to +2AA
            counter = 0
            while 1:
                if counter == len(indel_list):
                    break
                result = []
                ret = parse_indel(indel_list[counter], result)
                if ret < 0:
                    indel_list.remove(indel_list[counter])
                elif ret == 0:
                    indel_list[counter] = "".join(result)
                    counter += 1
                else:
                    counter += 1
            indel_count_list = Counter(indel_list).most_common()
            indel_name_list = [i[0] for i in indel_count_list]
            indel_count_list = [i[1] for i in indel_count_list]
            for i in indel_name_list:
                pileup[4] = pileup[4].replace(i, "I")

    # remove head(^) and tail($)
    pileup[4] = remove_head_end(pileup)
    
    # screen according to mapQ
    map_q_filter = [ord(i) - PHREDSCORE >= my_args.mapq for i in pileup[6]]
    base_q = [ord(j) - PHREDSCORE for i, j in enumerate(pileup[5]) \
        if map_q_filter[i]] 
    map_q = [ord(j) - PHREDSCORE for i, j in enumerate(pileup[6]) \
        if map_q_filter[i]] 

    readbase_list_with_i = split_readbase(pileup[4], indel_list)
    read_bases_filtered = compress(readbase_list_with_i, map_q_filter)

    i_filter = ["+" not in i and "-" not in i for i in read_bases_filtered]
    read_base_list_without_i = compress(read_bases_filtered, i_filter)

    # less min_depth
    if len(read_bases_filtered) < my_args.min_depth:
        return []
    
    # Count occurence of all possible mutations/bases
    var_names = ["A", "C", "G", "T"]
    var_counts = [read_base_list_without_i.count(i) for i in var_names]
    # Filter out mutations/bases that were not observed
    name_filter = [i > 0 for i in var_counts]
    var_names = compress(var_names, name_filter)
    var_counts = compress(var_counts, name_filter)

    var_names.extend(indel_name_list)
    var_counts.extend(indel_count_list)

    # Only reference bases reported
    if not var_counts:
        return -1

    # if max(var_counts) < my_args.minvar :
    #     return -1  # reference genome type, skip this line

    var_MQ = []
    var_BQ = []
    for var_name in var_names:
        q_filter = [j == var_name for j in read_bases_filtered]
        var_MQ.append(sum(compress(map_q, q_filter)))
        var_BQ.append(sum(compress(base_q, q_filter)))
    
    geno_class = get_geno_class(read_bases_filtered, map_q, base_q, var_counts,
        var_MQ, var_BQ)

    variant_str, r_num, v_num = get_variant_info(var_names, 
        [var_counts, var_MQ, var_BQ], read_bases_filtered, geno_class)
    sort_by_name(var_names, var_counts, var_MQ, var_BQ)

    var_all_str = ",".join(var_names)
    read_info_list = [(j, base_q[i]) for i, j in enumerate(readbase_list_with_i)]
    
    return PileupDataStruct(pileup[0], pileup[1], pileup[2], variant_str, r_num,
        v_num, len(var_names), read_info_list, None, geno_class, var_all_str)


def sort_by_name(var_names, var_counts, var_MQ, var_BQ):
    len_names = len(var_names)
    for i in range(len_names):
        for j in range(len_names):
            index = len_names - 1 - j - i
            if index == 0:
                break
            if var_counts[index] > var_counts[index - 1] \
                    or var_MQ[index] > var_MQ[index - 1] \
                    or var_BQ[index] >= var_BQ[index - 1]:
                for ex_list in [var_names, var_counts, var_MQ, var_BQ]:
                    ex_list[index], ex_list[index - 1] = \
                        ex_list[index - 1], ex_list[index]
    return var_names


def get_geno_class(read_bases_filtered, map_q, base_q, var_counts, var_MQ, var_BQ):
    if len(var_counts) < 2:
        return 1
    second_num = sorted(var_counts)[-2]

    ref_num = read_bases_filtered.count(".") + read_bases_filtered.count(",")
    if ref_num > second_num:
        return 1
    elif ref_num < second_num:
        return 2
    else:
        read_base_filter = [i in [".", ","] for i in read_bases_filtered]
        sum_ref_map_q = sum(compress(map_q, read_base_filter))
        try:
            v_map_q = var_MQ[var_counts.index(second_num)]
        except:
            import pdb; pdb.set_trace()

        if sum_ref_map_q > v_map_q:
            return 1
        elif sum_ref_map_q < v_map_q:
            return 2
        else:
            sum_ref_base_q = sum(compress(base_q, read_base_filter))
            v_base_q = var_BQ[var_counts.index(second_num)]
            
            if sum_ref_base_q >= v_base_q:
                return 1
            return 2


def get_golden_hetero(my_args, pileup, worker_id, min_map_q=20):
    """
    handle the data after vcf screen, store the result in q2
    :param pileup:
    :param min_map_q: mapQ
    :return: None  should skip
    """

    read_base = read_mpileup(pileup, 0, min_map_q, my_args.min_depth, True,
        worker_id)
    if not read_base or read_base == -1:
        return None
    ref_count = len([1 for i in read_base if i in [pileup[2], ".", ","]])
    alt_count = len(read_base) - ref_count

    if my_args.snp_type == 'dbsnp' \
            and (alt_count < 2 or ref_count < 2 or len(read_base) < my_args.RD):
        return None
    return GoldenHetero(pileup[0], pileup[1],  ref_count, alt_count)


def window_data_one_chromosome(center_pos, q2_list, lamb, is_list_ended,
            current_pos):
    """
    window_data from one chromosome data
    :type current_pos: int
    :param current_pos: current handling position.
    :param is_list_ended:
    :rtype: BigForewordList
    :type lamb: int
    :type q2_list: BigForewordList
    :type center_pos: int
    :param center_pos: pos of data in q3
    :param data_list: data after GH in q2(GoldenHetero)
    :param lamb: half of window width
    :return:
        None can not window yet.
        others data windowed
    """
    if q2_list.is_empty(): # is empty
        if is_list_ended:
            return []
        else:
            if current_pos - center_pos > lamb:
                return []
            return None

    # now data_list is not empty, looking for left edge
    if q2_list.get_last_element().pos < center_pos - lamb:
        if is_list_ended:
            return []
        if current_pos - center_pos > lamb:
            return []
        return None
    counter = 0  # type: int
    len_data_list = q2_list.len()
    while 1:
        if counter >= len_data_list:
            if is_list_ended:
                return []
            else:
                if current_pos - center_pos > lamb:
                    return []
                return None
        if q2_list.get_current_element().pos >= center_pos - lamb:
            left_edge = q2_list.get_current_element()
            if left_edge.pos > center_pos + lamb:
                return []
            break
        else:
            q2_list.move_foreword()
        counter += 1

    # right edge
    if q2_list.get_current_element().pos > center_pos + lamb:
        return []
    counter = 0  # type: int
    while 1:
        if counter == len_data_list:
            return []
        if q2_list.get_element(-1 - counter).pos <= center_pos + lamb:
            if is_list_ended or counter > 0:
                right_edge = q2_list.get_element(-1 - counter)
                break
            else:
                if current_pos - center_pos > lamb:
                    right_edge = q2_list.get_last_element()
                    break
                else:
                    return None
        else:
            counter += 1

    # double check edges
    if int(right_edge.pos) < int(left_edge.pos):
        return []
    return q2_list.filter(lambda x: left_edge.pos <= x.pos <= right_edge.pos)


# <INFO, NB> Calculation of GQ value (phred-scaled genotype quality score)
# Replaced with np array for clarity and performance
def get_gq(log_probs):
    # type: (float, float, float) -> str
    while 1:
        probs = 10 ** log_probs
        if probs.sum() - probs.min() != 0:
            break
        else:
            log_probs += 1

    norm_probs = probs / probs.sum()
    value = 1 - norm_probs.max()

    if value == 0:
        return 99
    else:
        return int(round(-10 * np.log10(value)))


def get_pl(log_lh):
    pl_raw = -10 * log_lh
    # <DONE, NB>
    pl = pl_raw - pl_raw.min()
    # <BUG, Original> Dont normalize Phred scores to smallest = 0
    # pl = pl_raw
    return '{:.0f},{:.0f},{:.0f},{:.0f}'.format(*pl)


def differential(my_args, q2_list, rel_pileup, q5, worker_id, q2_list_is_end,
        current_pos):
    """
    window the data in q2, send result to q5.
    :type my_format: str
    :param my_format:
    :type q2_list_is_end: bool
    :param q2_list_is_end:
    :type worker_id: int
    :param worker_id:
    :param artifact:
    :param default_bias:
    :param q5:
    :type current_pos: int
    :param current_pos: current data pos from source
    :type rel_pileup: PileupDataStruct
    :param rel_pileup:
    :param q2_list: data queue to be windowed
    :type q2_list: BigForewordList
    :param lamb: half of window width
    :type lamb: int
    :return:
    True data is windowed, send result to W
    False can not window yet
    """

    tracked_data = window_data_one_chromosome(rel_pileup.pos, q2_list,
        my_args.lamb, q2_list_is_end, current_pos)
    
    if tracked_data is None:
        return False

    bias = bias_estimator(rel_pileup.pos, tracked_data, my_args)

    # <INFO, NB> Calculation of basic values for GQ and PL
    # rr = Sequencing  Noise (rr_u = log10(rr))
    # ra = Amplification Artefact/Error (ra_u = log10(ra))
    # rm = Heterozygous SNV (rm_u = log10(rm))
    # mm = Homozygous SNV (mm_u = log10(mm))
    log_lh = sc_caller(rel_pileup, bias, my_args.null)
    lh = 10 ** log_lh

    if my_args.format == "bed":
        worker_type = WORKVAR
        gq = None
        pl = None
    else:
        worker_type = WORKVCF
        # <DONE, NB>
        gq = get_gq(np.array([np.log10(lh[0] + lh[1]), log_lh[2], log_lh[3]]))
        # <BUG, Original> Ignoring the sequencing  noise
        # gq = get_gq(log_lh[1:])
        pl = get_pl(log_lh)

    outline = OutLineStruct(rel_pileup.name, rel_pileup.pos, 
        rel_pileup.reference_base, rel_pileup.variant, 
        rel_pileup.reference_allele_num, rel_pileup.variant_allele_num, lh,
        rel_pileup.so, rel_pileup.variant_all, bias, gq, pl, 
        rel_pileup.genotype_num, rel_pileup.genotype_class)
 
    q5.put(QueueData(worker_id, worker_type, outline=outline), block=False)
    return True


def get_so_source(bulk_pileup_source, my_args):
    # type: (GeneratorExit, int, int, int) -> GeneratorExit
    """
    When pos in bulk greater than data_pos, should not continue to read bulk.
    should use the current bulkpos to compare the next data_pos.
    Build the generator, retain the bulkpos data, and achieve the above functions
    :param bulk_pileup_source: bulk generator
    :return:
        Datapos is not in bulk or bulk has no data: noCoverageInControl
        Datapos in bulk does not have enough depth: lessmindepth
        Datapos is ref in bulk: refgenotype
        Datapos is var in bulk: varreads
        -1 Pysam crushed and should be recalculated
    """
    should_read = True
    pos = -1
    while 1:
        if should_read:
            bulk_pileup_list = next(bulk_pileup_source)
            if bulk_pileup_list == -1:
                logging.info("Use samtools engine instead!")
                yield -1
            if bulk_pileup_list is None:
                while 1:
                    pos = yield "noCoverageInControl2"

        if int(bulk_pileup_list[1]) < pos:
            should_read = True
            continue
        elif int(bulk_pileup_list[1]) > pos:
            should_read = False
            pos = yield "noCoverageInControl"
        else:
            should_read = True
            pos = yield read_bulk_mpileup(bulk_pileup_list, my_args)


def read_bulk_mpileup(pileup, my_args):
    if "*" in pileup[4] or "-" in pileup[4] or "+" in pileup[4]:
        return "indel"
    pileup[4] = remove_head_end(pileup)
    MQ_filter = [ord(i) - PHREDSCORE >= my_args.bulk_min_mapq for i in pileup[6]]
    read_base = compress_read_bases(pileup[4], MQ_filter)
    if len(read_base) < my_args.bulk_min_depth:
        return "lessmindepth"
    maxvar = max([read_base.count('A'), read_base.count('C'),
        read_base.count('G'), read_base.count('T')])
    if maxvar < my_args.bulk_min_var:
        return "refgenotype"
    return "varreads"


def bias_estimator(pos, tracked_data, my_args):
    """
    calculate bias
    :type lamb: int
    :rtype: float
    :type tracked_data: list[GoldenHetero]
    :type default: float
    :type pos: int
    :param pos: center of window
    :param tracked_data: windowed data
    :param lamb: half of window width
    :param default: default value
    :return:
    """

    be_kwy = []
    be_kw = []
    for i in tracked_data:
        try:
            be_tmp1 = float(i.refcount + i.altcount)  # type: float
            be_tmp2 = max(i.refcount, i.altcount) / be_tmp1
        except AttributeError:
            be_tmp1 = float(i.reference_allele_num + i.variant_allele_num)
            be_tmp2 = max(i.reference_allele_num, i.variant_allele_num) / be_tmp1

        if be_tmp1 <= 0:
            continue

        be_tmp = bias_kernel(int(pos), int(i.pos), my_args.lamb)  # K in the formula
        be_kwy.append(be_tmp * be_tmp1 * be_tmp2)
        be_kw.append(be_tmp * be_tmp1)
    # Nadaraya-Watson kernel-weighted average
    if len(be_kwy) > 0 and sum(be_kw) > 0:
        return sum(be_kwy) / sum(be_kw)
    # return args.bias when no neighboring heterozygous base
    return my_args.bias


def bias_kernel(bk_x0, bk_xi, lamb):
    # Equation 2
    if -lamb < bk_x0 - bk_xi < lamb:
        return 0.75 * (1 - (float(bk_x0 - bk_xi) / lamb) ** 2)
    else:
        return 0.0


def sc_caller(sc_candidate, sc_bias, min_frac):
    """ Calculate RR, RA, RM, MM from data in q3 and bias.
    """
    if "," in sc_candidate.variant:
        ref, mut = sc_candidate.variant.split(",")
    else:
        ref = sc_candidate.reference_base
        mut = sc_candidate.variant

    # <DONE, NB>
    f_h0 = 0.125 * sc_bias
    # <BUG, Original> Ignoring the bias
    # f_h0 = 0.125

    f_h1 = sc_bias
    f_h2 = 1
    if sc_candidate.reference_allele_num > sc_candidate.variant_allele_num:
        f_h1 = 1 - f_h1
        # <TODO, NB> Confirm if Eq. 6 from the manuscript is wrong
        # f_h2 = 0
        # <BUG, Original> Ignoring the 1 - F_K(\theta) if ref count > alt count
        f_h2 = 1
  
    log_lh = np.zeros(4)
    for idx, f in enumerate([min_frac, f_h0, f_h1, f_h2]):
        log_lh[idx] = np.log10(
            [P_b_GG(i, ref, mut, f) for i in sc_candidate.read_info_list]).sum()

    return log_lh


def P_b_GG(bp, ref, mut, f):
    """ Equation (7), calculate lg Probability value
    """
    e = 10 ** (-bp[1] / 10.0)

    if bp[0] in [',', '.', ref]:
        a = f * e / 3 + (1 - f) * (1 - e)
    elif bp[0] == mut:
        a = (1 - f) * e / 3 + f * (1 - e)
    else:
        a = e / 3

    return a


def get_my_filename(output, suffix):
    file_name_base = os.path.splitext(os.path.basename(output))[0]
    file_name = '{}{}'.format(file_name_base, suffix)
    return os.path.join(os.path.dirname(output), file_name)


def calculate_eta(list_var_buf, list_var_tag):
    allele_nums = []
    bias_list = []
    for i, var_tag in enumerate(list_var_tag):
        if not var_tag:
            break
        for q5_item in list_var_buf[i]:
            allele_nums.append(q5_item.outline.total_num)
            bias_list.append(q5_item.outline.bias) 
    
    allele_nums = allele_nums[:CUTOFFNUM]
    bias_list = bias_list[:CUTOFFNUM]

    LLR = np.zeros(len(allele_nums))
    for i, allele_num in enumerate(allele_nums):
        f_artifact = 0.125 * bias_list[i] / 0.5
        # bin (depth(number of trials),prob_success)
        alt = np.random.binomial(allele_num, f_artifact)
        ref = allele_num - alt
        L_filter = (1 - f_artifact) ** ref * f_artifact ** alt
        ## random select major allele
        major = np.random.randint(2)
        if major == 0:
            f_true = 0.5 * 0.5 / bias_list[i]
        if major == 1:
            f_true = 0.5 * bias_list[i] / 0.5
        L_true = (1 - f_true) ** ref * f_true ** alt
        ## if L_filter/true is 0, assign a very small value
        if L_filter == 0:
            L_filter = 10 ** -100
        if L_true == 0:
            L_true = 10 ** -100

        ratio = L_filter / L_true
        if ratio != 0:
            LLR[i] = -np.log10(ratio)
        else:
            LLR[i] = 10000

    co_001 = np.percentile(LLR, 99)
    result = 10 ** -co_001
    return result


def write_result(q5, args, name):
    """
    Receive the data in q5, organize it, and write the result file.
    :return:
    """
    logging.info("Waiting for data coming")
    
     # Initialize data for each worker
    list_var_buf = [[] for i in range(args.work_num)]
    list_var_tag = [False for i in range(args.work_num)]
    list_coverage_buf = [[] for i in range(args.work_num)]
    
    reasoning_file = get_my_filename(args.output,
        "_{:0>2d}to{:0>2d}.reasoning".format(args.head, args.tail))
    cov_file = get_my_filename(args.output,
        "_{:0>2d}to{:0>2d}.coverage".format(args.head, args.tail))
    body_file = '{}.body'.format(args.output)
    eta_file = '{}.eta'.format(args.output)

    eta = -1
    current_handling_worker = 1

    with open(body_file, "a") as fp_out, \
            open(reasoning_file, "a") as fp_reasoning, \
            open(eta_file, "a") as fp_eta:
        while 1:
            msg_q5 = q5.get(block=True)
            w_id = msg_q5.worker_id

            # from main
            if msg_q5.mode == WORKEXIT:
                # try to write coverage file
                merged_coverage_list = merge_coverage(list_coverage_buf)
                if merged_coverage_list:
                    cov_list = map(lambda z: "\t".join(z),
                        map(lambda x: map(lambda y: str(y), x), 
                            merged_coverage_list))
                    import pdb; pdb.set_trace()
                    with open(cov_file, "a") as fp_cov:
                        fp_cov.write("\n".join(cov_list))
                if eta == -1:
                    if get_current_cutoff_num(list_var_buf, list_var_tag) > 0:
                        eta = calculate_eta(list_var_buf, list_var_tag)
                        for i, var_buf in enumerate(list_var_buf):
                            if var_buf:
                                write_do([j.outline for j in var_buf],
                                    msg_q5.work_type, i, fp_out, eta,
                                    fp_reasoning, args)
                            else:
                                logging.info("worker{} write len = 0".format(i))
                    else:
                        logging.info("Nothing need to be written.")
                fp_eta.write("##contig=<ID={},eta={}>\n".format(name, eta))
                break

            elif msg_q5.work_type == WORKVAR or msg_q5.work_type == WORKVCF:

                # record tag and data
                if msg_q5.mode == WORKDONE:
                    list_var_tag[w_id] = True
                    logging.info("chr{} worker{} done".format(name, w_id))
                    # Try to calculate eta, write data
                    if get_current_cutoff_num(list_var_buf, list_var_tag) \
                            >= CUTOFFNUM and eta == -1:
                        eta = calculate_eta(list_var_buf, list_var_tag)
                    if eta != -1:
                        for i, var_buf in enumerate(list_var_buf):
                            if i < current_handling_worker - 1:
                                continue
                            if list_var_tag[i]:
                                if var_buf:
                                    write_do([j.outline for j in var_buf],
                                        msg_q5.work_type, i + 1, fp_out, eta,
                                        fp_reasoning, args)
                                else:
                                    logging.info("Chr{} worker{} write len=0" \
                                        .format(name, i + 1))
                                current_handling_worker += 1
                            else:
                                break
                elif msg_q5.mode == WORKERROR:
                    logging.info("worker{} pysam crashed. Data cleaned" \
                        .format(w_id))
                    list_var_tag[w_id] = False
                    list_var_buf[w_id] = []
                    list_coverage_buf[w_id] = []
                else:
                    list_var_buf[w_id].append(msg_q5)

            else:
                list_coverage_buf[w_id].append(msg_q5.coverage)

    if args.bulk == "" or msg_q5.work_type == WORKVCF:
        os.remove(reasoning_file)

    logging.info("Quit!")


def get_current_cutoff_num(list_var_buf, list_var_tag):
    cutoff_counter = 0
    for i, var_tag in enumerate(list_var_tag):
        if not var_tag:
            break
        cutoff_counter += len(list_var_buf[i])
    return cutoff_counter


def merge_coverage(list_coverage_buf):
    # list_coverage_buf: [[coverage1, coverage2...]]
    result_list = []
    for coverage_buf in list_coverage_buf:
        if not result_list:
            result_list.extend(coverage_buf)
            continue
        if not coverage_buf:
            continue

        if coverage_buf[0][1] == result_list[-1][2] + 1:
            result_list[-1][2] = coverage_buf[0][2]
            result_list.extend(coverage_buf[1:])
        else:
            result_list.extend(coverage_buf)
    return result_list


def write_do(data_list, work_type, worker_id, fp_out, eta, fp_reasoning, my_args):
    if work_type == WORKVAR:
        res_filter = [i.he != 0 and i.ae / i.he < eta \
                and my_args.min_var_frac * i.var_num > i.ref_num \
                and 2 * i.sn < i.he \
            for i in data_list]
        data_list = [j for i, j in enumerate(data_list) if res_filter[i]]

        if len(data_list) > 0:
            fp_out.write("\n".join([i.get_varcall_str() for i in data_list]))
            fp_out.write("\n")
        if my_args.bulk != "" and data_list:
            fp_reasoning.write("\n".join(
                [i.get_reason_str() for i in data_list]))
            fp_reasoning.write("\n")
    elif work_type == WORKVCF:
        fp_out.write("\n".join(
            [i.get_vcf_str(eta, my_args.minvar, my_args.minvarfrac) \
                for i in data_list]))
        fp_out.write("\n")

    logging.debug("worker{} write len={}".format(worker_id, len(data_list)))


def data_generator(my_args, name, start, stop, is_bulk):
    if my_args.engine == "samtools":
        generator = data_generator_samtools
    else:
        generator = data_generator_pysam
    return generator(my_args, name, start, stop, is_bulk)


def data_generator_pysam(my_args, name, start, stop, is_bulk):
    fasta_file = pysam.FastaFile(my_args.fasta)
    str_ref = fasta_file.fetch(name, start - 1, stop + 1)

    my_arg = {"fastafile": fasta_file, "stepper": "samtools",
        "adjust_capq_threshold": 50, "contig": name, "start": start,
        "stop": stop, "min_mapping_quality": 0 if is_bulk else 40}

    if is_bulk:
        bam_file = pysam.AlignmentFile(my_args.bulk, "rb")
    else:
        bam_file = pysam.AlignmentFile(my_args.bam, "rb")
    read_bases_list = []
    for pileup_column in bam_file.pileup(**my_arg):
        pos = pileup_column.reference_pos + 1

        if pos > stop:
            break
        if pos < start:
            continue

        try:
            read_bases_list = pileup_column \
                .get_query_sequences(mark_matches=True, mark_ends=True,
                    add_indels=True)
        except Exception as e:
            logging.debug("Pysam crashed! Unexpected Error: {}".format(e))
            yield -1

        read_bases = ''.join(read_bases_list)
        if len(read_bases) == 0:
            read_bases = "*"

        base_qualities = "".join([chr(int(i) + 33) \
            for i in pileup_column.get_query_qualities()])
        if len(base_qualities) == 0:
            base_qualities = "*"

        mapq = "".join([chr(int(i) + 33) \
            for i in pileup_column.get_mapping_qualities()])
        if len(mapq) == 0:
            mapq = "*"

        result = [name, pos, str_ref[pos - start], 
            str(pileup_column.get_num_aligned()), read_bases.upper(),
            base_qualities, mapq]

        yield result
    yield None


def data_generator_samtools(my_args, name, start, stop, is_bulk):
    if is_bulk:
        cmd_str = "samtools mpileup -C50  -f {} -s {} -r {}:{}-{}" \
            .format(my_args.fasta, my_args.bam, name, start, stop)
    else:
        cmd_str = "samtools mpileup -C50  -f {} -q 40 -s {} -r {}:{}-{}" \
            .format(my_args.fasta, my_args.bam, name, start, stop)
    process_samtools = Popen([cmd_str], shell=True, stdout=PIPE)
    while 1:
        pileup = []
        if not read_line_from_process(process_samtools, "#", "\t", 7, pileup):
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
        buf = buf.strip("\n").split(spliter)
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


def control(my_args, list_vcf, q5, name, head, stop, worker_id):
    """
    :param stop: Unexpanded region
    :type my_args: object
    :type list_vcf: list
    :type name: str
    :type q5: object (main queue)
    :type worker_id: int
    :type stop: int
    :type head: int
    """
    logging.info("worker {} begin! name={} head={} tail={} len={} engine={}" \
        .format(worker_id, name, head, stop, stop - head, my_args.engine))

    rel_pileups = BigForewordList([])
    q2_list = BigForewordList([])
    q2_list_is_end = False
    my_vcf_list = BigForewordList(list_vcf)
    
    copy_vcf_list = copy.copy(list_vcf)
    # Expanding the edge
    if worker_id != 0:
        my_start = head - my_args.lamb
    else:
        my_start = head

    if worker_id != my_args.work_num -1:
        my_stop = stop + my_args.lamb
    else:
        my_stop = stop
    total_work = my_stop - my_start

    # data source
    pileup_source = data_generator(my_args, name, my_start, my_stop, False)
    # bulk source
    if my_args.bulk != '':
        bulk_pileup_source = data_generator(my_args, name, head, stop, True)
        so_source = get_so_source(bulk_pileup_source, my_args)
        next(so_source)

    if my_args.debug:
        print('w{:<3d}\t0/{:<6}'.format(worker_id, total_work))

    coverage_list = []  # list of [name str, head int, stop int]
    counter = 0
    main_index = 0
    while 1:
        pileup = next(pileup_source)

        if pileup == -1:  # pysam crashed
            logging.info("Use samtools engine instead!")
            error_data = QueueData(worker_id, WORKVAR, mode=WORKERROR)
            q5.put(error_data, block=False)
            my_args.engine = "samtools"
            control(my_args, copy_vcf_list, q5, name, head, stop, worker_id)
            return
        elif pileup is None:
            if my_args.coverage and len(coverage_list) > 0:
                done_data = QueueData(worker_id, WORKCOVERAGE, mode=WORKDONE)
                q5.put(done_data, block=False)
                del coverage_list[0]
            break

        # append the relevant pileups into rel_pileups
        if head <= pileup[1] < stop:
            # Check if pileup is relevant/passes filters
            if my_args.format == "bed":
                rel_pileup = read_mpileup(copy.copy(pileup),
                    my_args.minvar, my_args.mapq, my_args.min_depth, False,
                    worker_id) 
            else:
                rel_pileup = read_mpileup_vcf(copy.copy(pileup), my_args)

            if rel_pileup and rel_pileup != -1:
                if my_args.bulk != "":
                    #  calculate so
                    result = so_source.send(rel_pileup.pos)
                    if result == -1:
                        logging.info("Use samtools engine instead!")
                        error_data = QueueData(worker_id, WORKVAR, mode=WORKERROR)
                        q5.put(error_data, block=False)
                        my_args.engine = "samtools"
                        control(my_args, copy_vcf_list, q5, name, head, 
                            stop, worker_id)
                        return
                    else:
                        rel_pileup.so = result
                else:
                    rel_pileup.so = ""
                rel_pileups.append(rel_pileup)

            # calculate coverage
            if my_args.coverage:
                if rel_pileup and (my_args.bulk == "" \
                        or (rel_pileup == -1 or (rel_pileup != -1 \
                            and rel_pileup.so in ["", "varreads", "refgenotype"]))):

                    if not coverage_list:  # is empty
                        coverage_list.append([name, pileup[1], pileup[1]])
                    elif pileup[1] == coverage_list[-1][2] + 1:
                        coverage_list[-1][2] = pileup[1]
                    else:
                        coverage_list.append([name, pileup[1], pileup[1]])
                        new_coverage = coverage_list.pop(0)
                        coverage_data = QueueData(worker_id, WORKCOVERAGE, 
                            coverage=new_coverage)
                        q5.put(coverage_data, block=False)

        # handle the head data in rel_pileups
        if not rel_pileups.is_empty():
            diff_success = differential(my_args, q2_list, 
                rel_pileups.get_current_element(), q5, worker_id, q2_list_is_end,
                pileup[1])
            if diff_success:
                rel_pileups.move_foreword()

        if my_args.debug:
            counter += 1
            if counter == 10000:
                counter = 0
                main_index += 1
                print("w{:<3d} {:>6}/{:<6}"\
                    .format(worker_id, main_index * 10000, total_work))
        # vcf screen
        if not q2_list_is_end:
            while not my_vcf_list.is_empty():
                if my_vcf_list.get_current_element() < pileup[1]:
                    my_vcf_list.move_foreword()
                else:
                    break
            if pileup[1] == my_vcf_list.get_current_element():
                ret = get_golden_hetero(my_args, pileup, worker_id)
                if ret is not None:
                    q2_list.append(ret)
            if my_vcf_list.is_empty() \
                    or pileup[1] >= my_vcf_list.get_last_element():
                q2_list_is_end = True

    while not rel_pileups.is_empty():
        diff_success = differential(my_args, q2_list,
            rel_pileups.get_current_element(), q5, worker_id, True, my_stop)
        if diff_success:
            rel_pileups.move_foreword()

    worker_type = WORKVAR if my_args.format == "bed" else WORKVCF
    done_data = QueueData(worker_id, worker_type, mode=WORKDONE)
    q5.put(done_data, block=False)

    if my_args.debug:
        print('w{:<3d} I am done!'.format(worker_id))
    logging.info("worker {} Quit!".format(worker_id))


def load_vcf(name, vcf_info, vcf_file_name):
    """ Read the position information of the specified name chromosome from the 
        vcf file, and store it as a list
    """
    vcf_list = []
    curr_info = [i for i in vcf_info if i[0] == name]
    if curr_info:
        with open(vcf_file_name) as f:
            for i in curr_info:
                f.seek(i[1])
                lines = f.read(i[2]).splitlines()
                vcf_list.extend([int(line.split('\t')[1]) for line in lines])
    return vcf_list


def parse_vcf(vcf_file, snp_type):
    """ Parse the vcf file into [[name,head,length],[name,head,length]...] format
    """
    has_shown_info = False
    base_file, file_type = os.path.splitext(vcf_file)

    # read vcf catalog
    catalog_file = "{}.catalog".format(base_file)
    if os.path.exists(catalog_file) \
            and os.path.getmtime(catalog_file) > os.path.getmtime(vcf_file):
        with open(catalog_file, "r") as fp:
            catalog = fp.read().split(";")
        result = map(
            lambda x: map(lambda y: y if x.split(",").index(y) == 0 else int(y),
                x.split(",")),
            catalog
        )
        res = [map(lambda y: y if x.split(",").index(y) == 0 else int(y), x.split(",")) \
            for x in catalog]
        import pdb; pdb.set_trace()
    # parse vcf
    else:
        if file_type == '.gz':
            raise IOError('VCF file {} needs to be unzipped'.format(vcf_file))
        result = []
        name = ""
        head = 0
        with open(vcf_file, 'r') as f:
            while True:
                line = f.readline()
                if line == "":
                    result.append([name, head, f.tell() - head])
                    break
                if line[0] == "#" or line == "\n":
                    continue
                
                cols = line.split("\t")
                if name != cols[0]:
                    curr_idx = f.tell() - len(line)
                    if name == "":
                        head = curr_idx
                    else:
                        tail = curr_idx
                        result.append([name, head, tail - head])
                        head = tail
                    name = cols[0]

                if snp_type == "hsnp" and not has_shown_info and len(cols) > 8:
                    tmp_str = "\t".join(cols[9:])
                    tmp_list = re.findall("0\\|0|0/0|1\\|1|1/1|2/2|2\\|2", tmp_str)
                    if len(tmp_list) > 0:
                        logging.info(">>>Please confirm the input VCF only " \
                            "contains heterozygote loci in bulk!<<<")
                        has_shown_info = True

    return result


def main(my_args):
    # Init logging
    log_file = get_my_filename(my_args.output, "_sccaller_{:0>2d}to{:0>2d}.log" \
        .format(my_args.head, my_args.tail))
    logging.basicConfig(filename=log_file, level=logging.DEBUG,
        format=LOG_FORMAT, filemode="w")

    logging.info("Welcome to SCcaller v{}".format(VERSION))

    if my_args.debug:
        my_args.cpu_num = 1
        my_args.work_num = 1

    if my_args.debug:
        print("\nRunning SCcaller v{} in debugging mode".format(VERSION))
        print('parsing vcf...')
    logging.info("parsing SNP vcf...")
    vcf_info = parse_vcf(my_args.snp_in, my_args.snp_type)

    if my_args.debug:
        print('parsing fasta...')
    logging.info("parsing reference fasta...")
    fasta_info = parse_fasta(my_args.fasta)

    if my_args.tail == -1 or my_args.tail > len(fasta_info):
        my_args.tail = len(fasta_info)

    logging.info("args: {}".format(my_args))

    queue = mp.Manager().Queue()
    # Start the operation process
    for j, fasta_j in enumerate(fasta_info):
        if j + 1 < my_args.head:
            continue
        if j + 1 > my_args.tail:
            break

        name, start, stop = fasta_j
        
        if my_args.debug:
            print('loading vcf...')
        logging.info("loading vcf...")
        list_vcf = load_vcf(name, vcf_info, my_args.snp_in)
        
        if my_args.debug:
            print('SCcaller v2.0.0 is handling chromosome {}...'.format(name))
        logging.debug("name={}; list_vcf len={}".format(name, len(list_vcf)))

        # Calculate the amount of tasks for each process
        no_bases = stop - start
        step = int(np.ceil(no_bases / float(my_args.work_num)))
        if step < my_args.lamb:
            logging.info("work_num={} is too large for chr={}. Using 1 instead." \
                .format(my_args.work_num, name))
            my_args.work_num = 1
            step = no_bases

        # Start the write file process
        proc_write_result = mp \
            .Process(target=write_result, args=(queue, my_args, name))
        proc_write_result.daemon = True
        proc_write_result.start()

        if my_args.work_num > 1:
            process_pool = mp.Pool(processes=my_args.cpu_num)

        for woker_id in range(my_args.work_num):
            head = start + woker_id * step
            tail = head + step
            if woker_id != 0:
                if (head - my_args.lamb) < 0:
                    err_str = "lamb={} is too large. Fasta: name={} ({} - {})." \
                        .format(my_args.lamb, name, start, stop)
                    logging.critical(err_str)
                    raise RuntimeError(err_str)

            vcf_worker = [i for i in list_vcf \
                if head - my_args.lamb <= i <= tail + my_args.lamb]

            if my_args.work_num > 1:
                process_pool.apply_async(control, 
                    (my_args, vcf_worker, queue, name, head, tail, woker_id))
            else:
                control(my_args, vcf_worker, queue, name, head, tail, woker_id)
                
        if my_args.work_num > 1:
            process_pool.close()
            process_pool.join()

        # Exit the W process
        worker_type = WORKVAR if my_args.format == "bed" else WORKVCF
        queue.put(QueueData(0, worker_type, mode=WORKEXIT), block=True)
        proc_write_result.join()

    write_vcf(my_args)

    logging.info("W quit. All done.")


def write_vcf(my_args):
    """ Write the vcf file
    """
    head_str = "##fileformat=VCFv4.1\n" \
        "##fileDate={date}\n" \
        "##source=SCcaller_v{v}\n" \
        "##reference=file:{ref}\n" \
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=" \
        "\"Number of Samples With Data\">\n" \
        "##FILTER=<ID={mg},Description=\"Multiple genotype\">\n" \
        "##FILTER=<ID={nev}{mv},Description=\"Number of variant supporting " \
        "reads <{mv}\">\n" \
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" \
        "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"" \
        "Allelic depths for the ref and alt alleles in the order listed\">\n" \
        "##FORMAT=<ID=BI,Number=1,Type=Float," \
        "Description=\"Amplification Bias\">\n" \
        "##FORMAT=<ID=GQ,Number=1,Type=Integer," \
        "Description=\"Genotype Quality\">\n" \
        "##FORMAT=<ID=FPL,Number=4,Type=Integer,Description=\"" \
        "sequencing noise, amplification artifact, heterozygous SNV and " \
        "homozygous SNV respectively\">\n" \
            .format(date=time.strftime("%Y:%m:%d-%H:%M:%S", time.localtime()),
                ref=my_args.fasta, mg=MULTIPLEGENOTYPE, v=VERSION,
                nev=NOTENOUGHVARIANTS, mv=my_args.minvar)

    if my_args.bulk != "":
        head_str += "##FORMAT=<ID=SO,Number=1,Type=String," \
            "Description=\"Whether it is a somatic mutation.\">\n"

    eta_file = '{}.eta'.format(my_args.output)
    with open(eta_file, 'r') as f:
        etas = f.read()
    head_str += etas
    os.remove(eta_file)

    sc_name = os.path.basename(my_args.bam).split('.')[0]
    head_str += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n" \
        .format(sc_name)

    body_file = '{}.body'.format(my_args.output)
    with open(body_file, "r") as f:
        body_str = f.read()
    os.remove(body_file)

    with open(my_args.output, "w") as fp_out:
        if my_args.format == "vcf":
            fp_out.write(head_str)
        fp_out.write(body_str)


if __name__ == "__main__":
    start = time.time()
    args = parse_args()
    main(args)
    print("Duration = {}s".format(str(time.time() - start)))