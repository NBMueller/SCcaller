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


import sys
import os
import re
from datetime import datetime
import random
import copy
import argparse
import logging
import multiprocessing as mp
from collections import Counter
# Additional libraries
import numpy as np
# SCcaller internal libraries
import libs.SCcallerIO as io
from libs.BigForewordList import BigForewordList
from libs.OutLineStruct import OutLineStruct


if float(sys.version[:3]) != 2.7:
    print('\nWARNING: Python 3.X not tested for vcf .catalog file (for dbsnp) '
        'and bed output (instead of vcf)!\n')


# regular expression of indels
INSERTION = "\+[0-9]+[ACGTNacgtn]+"
DELETION = "\-[0-9]+[ACGTNacgtn]+"
INDEL = "\+[0-9]+[ACGTNacgtn]+|\-[0-9]+[ACGTNacgtn]+"

# Fixed queue indicators
WORKVAR = 1
WORKCOVERAGE = 2
WORKVCF = 3
WORKDONE = 'done'
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
    def __init__(self, name, pos, ref_num, var_num):
        self.name = name  # type: str
        self.pos = pos  # type: int
        self.ref_num = ref_num  # type: int
        self.var_num = var_num  # type: int


    def __str__(self):
        return '{}:{}\tAD:{},{}' \
            .format(self.name, self.pos, self.ref_num, self.var_num)


class PileupDataStruct:
    def __init__(self, name, pos, ref_base, var_base, ref_num, var_num, gt_num,
                read_infos, gt_class=-1, var_all="", so=None):
        # type: (str, int, str, str, int, int, int, list[tuples], str, int) -> object
        self.name = name  # type: str
        self.pos = pos  # type: int
        self.ref_base = ref_base  # type: str
        self.var_base = var_base  # type: str
        self.ref_num = ref_num  # type: int
        self.var_num = var_num  # type: int
        self.gt_num = gt_num  # type: int
        self.read_infos = read_infos  # type: list[tuples]
        self.gt_class = gt_class  # type: int  # 0 unknown 1:0,0 0,1 1,1      2:1,1 1,2 2,2
        self.var_all = var_all
        self.so = so  # type: str  # reasoning  bulk ref -- "True"   bulk var -- "False"  else -- "NA"


    def __str__(self):
        return '{}:{}\tref:{}({}),alt:{}({})' \
            .format(self.name, self.pos, self.ref_base, self.ref_num,
                self.var_base, self.var_num)


def parse_args():
    parser = argparse.ArgumentParser(
        prog='SCcaller',
        description='''SCcaller v{}_NB;
            Xiao Dong, biosinodx@gmail.com, xiao.dong@einstein.yu.edu;
            Yujue Wang, spsc83@gmail.com

            (Modifications & refactoring: nico.borgsmueller@bsse.ethz.ch)''' \
                .format(VERSION)
    )
    # Input related arguments
    parser.add_argument("-b", "--bam", type=str, required=True, 
        help="BAM file of a single cell sequencing experiment.")
    parser.add_argument("-f", "--fasta", type=str, required=True, 
        help="FASTA file of the reference genome.")
    parser.add_argument("-s", "--snp_in", type=str, required=True,
        help="Candidate snp input file, either from dbsnp data or heterozygous "
            "snp (hsnp) data of the bulk, for known heterozygous calls. "
            "File type: bed (1-based) or vcf.")
    parser.add_argument("-t", "--snp_type", type=str, required=True,
        choices=["dbsnp", "hsnp"],
        help="SNP type for --snp_in. When choosing dbsnp, "
            "--bulk bulk_bamfile is required.")
    parser.add_argument("-d", "--wkdir", type=str, default="./",
        help="Work dir. Default: ./")
    # Cutoff/threshold related arguments
    parser.add_argument("--minvar", type=int, default=2,
        help="Min. # variant supporting reads. Default: 2")
    parser.add_argument("-mvf", "--minvcalarfrac", type=float, default=0.2,
        help='Min. fraction of variant supporting reads for a 0/1 genotype. '
            'Default: 0.2')
    parser.add_argument('-a', '--lrt_alpha', default=0.05, type=float,
        help='Significance level for calculating eta with LRT via simulations')
    parser.add_argument("--mapq", type=int, default=40,
        help="Min. mapQ. Default: 40")
    parser.add_argument("--min_depth", type=int, default=10,
        help="Min. # reads. Default: 10")
    # Paramter related arguments
    parser.add_argument("-l", "--lamb", type=int, default=10000,
        help="Lambda for bias estimation. Default: 10000.")
    parser.add_argument("--bias", type=float, default=0.75,
        help="Default theta (bias) for SNVs whose theta cannot be estimated. "
            "Default: 0.75.")
    parser.add_argument("--null", type=float, default=0.03,
        help="Min allelic fraction considered. Default=0.03.")
    
    parser.add_argument("--head", type=int, default=1,
        help="First chromosome as sorted as in fasta file to analyze (1-based). "
            "Default: 1")
    parser.add_argument("--tail", type=int, default=-1,
        help="Last chromosome as sorted as in fasta file to analyze (1-based). "
        "Default: -1")
    # Output related arguments
    parser.add_argument("-o", "--output", type=str, default='', 
        help='Output file name. Default: <BAM_WITHOUT_EXTENSION>.calls.vcf|bed')
    parser.add_argument("--format", type=str, choices=["vcf", "bed"],
        default="vcf", help="Output file format. Default: vcf")
    parser.add_argument("--coverage", action="store_true", default=False,
        help="use \"--coverage\" to generate the coverage file at the same time")
    # Bulk related arguments
    parser.add_argument("--bulk", type=str, default="",
        help="bamfile of bulk normal DNA sequencing")
    parser.add_argument("--bulk_min_var", type=int, default=1,
        help="Min. num. variant supporting reads for bulk. Default: 1")
    parser.add_argument("--bulk_min_mapq", type=int, default=20,
        help="Min. mapQ for bulk. Default: 20")
    parser.add_argument("--bulk_min_depth", type=int, default=20,
        help="Min. reads for bulk. Default: 20")
    parser.add_argument("--RD", type=int, default=15,
        help="Min. read depth of known heterozygous SNPs called from bulk when "
            "choosing -t dbsnp. Default: 20. Recommended: 10,15,20, depending "
            "on average read depth.")
    # Miscellaneous arguments
    parser.add_argument("-e", "--engine", type=str, choices=["pysam", "samtools"],
        default="pysam", help="Pileup engine. Default: pysam")
    parser.add_argument("-w", "--work_num", type=int, default=100,
        help='Number of threads to use. If you run into memory issues, try '
            'increasing the number. Default: 100')
    parser.add_argument("-n", "--cpu_num", type=int, default=1,
        help="Num. processes. Default: 1")
    parser.add_argument("--debug", action='store_true',
        help="Turn of multiprocessing/threading for debugging. Default: False")
    parser.add_argument("--seed", type=int, default=-1,
        help="Seed to use for random number generation. Default: random")

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
    if args.output == '':
        args.output = '{}.calls.{}' \
            .format(os.path.splitext(args.bam)[0], args.format)

    if os.path.exists(args.output):
        print('\n>>>WARNING: Output file {} exists and will be overwritten<<<\n' \
            .format(args.output))
    for out_ending in (".coverage", ".reasoning"):
        if os.path.exists(get_my_filename(args.output, out_ending)):
            os.remove(get_my_filename(args.output, out_ending))

    if args.format == 'bed':
        args.worker_type = WORKVAR 
    else:
        args.worker_type = WORKVCF

    return args


def compress(l, f):
    return [j for i, j in enumerate(l) if f[i]]


def parse_indel(indel_str, indel_list_out):
    """ parse indel string, return the indel string and length of indel string
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
    """ calculate the mapQ of read head (^) from readbase string without I
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
    """ add indel to the end of the readbase without I, return as list
    """
    read_base_list = []
    read_base_list.extend(list(read_base_without_i))
    tmp = [[indel_name_list[list_index] \
                for i in range(indel_count_list[list_index])]
            for list_index in range(len(indel_count_list))]
    for i in tmp:
        read_base_list.extend(i)
    return read_base_list


def read_mpileup(pileup, rm_minvar, min_mapq, rm_mindepth, is_gh):
    """ screenfor mutations, filter out alignments that do not meet the 
        requirements of mapQ, and convert the data format
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
    MQ_filter = [ord(i) - PHREDSCORE >= min_mapq for i in pileup[6]]
    for i, base in enumerate(pileup[4]):
        if base == "I":
            MQ_filter.insert(i, True)
    read_base_with_i = "".join(compress(pileup[4], MQ_filter))
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

    read_base_with_d = "".join([j if MQ_filter[i] else "D" \
        for i, j in enumerate(pileup[4])])
    ref_num, var_num = get_reference_variant_allele_num(read_base_with_d)
    
    read_base_list_final = rebuild_read_base_list(read_base_without_i,
        indel_name_list, indel_count_list)
    
    base_q = [ord(j) - PHREDSCORE for i,j in enumerate(pileup[5]) \
        if MQ_filter[i]]
    num_of_i = len(read_base_with_i) - len(read_base_without_i)
    base_q.extend([DEFAULTBASEQUALITY] * num_of_i)

    read_infos = [(j, base_q[i]) for i, j in enumerate(readbase_list_with_i)]

    return PileupDataStruct(pileup[0], pileup[1], pileup[2], var_max, ref_num,
        var_num, 1, read_infos)


def choose(candidate_index, current_basis, choice_out, basis_filter):
    """ choose indexs from the candidate_index according to current_basis
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
    """ calculate the variant's name_str, ref_num and var_num
        if variant is more than 2, choose 2 according to basis.
        If basis doesn't work, choose randomly.
    :param basis: [reads, map_q, base_q]
    :param var_names:
    :param ead_bases_filtered:
    :return: [variant_str, ref_num, var_num]
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
    MQ_filter = [ord(i) - PHREDSCORE >= my_args.mapq for i in pileup[6]]
    base_q = [ord(j) - PHREDSCORE for i, j in enumerate(pileup[5]) \
        if MQ_filter[i]] 
    map_q = [ord(j) - PHREDSCORE for i, j in enumerate(pileup[6]) \
        if MQ_filter[i]] 

    readbase_list_with_i = split_readbase(pileup[4], indel_list)
    read_bases_filtered = compress(readbase_list_with_i, MQ_filter)

    i_filter = ["+" not in i and "-" not in i for i in read_bases_filtered]
    read_base_list_without_i = compress(read_bases_filtered, i_filter)

    # less min_depth
    if len(read_bases_filtered) < my_args.min_depth:
        return []
    
    # Count occurence of all possible mutations/bases
    var_names = ["A", "C", "G", "T"]
    var_counts = [read_base_list_without_i.count(i) for i in var_names]
    # Filter out mutations/bases that were not observed
    name_filter = [i > my_args.minvar for i in var_counts]
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

    var_all_str = ",".join(sorted(var_names))
    read_infos = [(j, base_q[i]) for i, j in enumerate(readbase_list_with_i)]
    return PileupDataStruct(pileup[0], pileup[1], pileup[2], variant_str, r_num,
        v_num, len(var_names), read_infos,
        gt_class=geno_class, var_all=var_all_str)


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

    read_base = read_mpileup(pileup, 0, min_map_q, my_args.min_depth, True)
    if not read_base or read_base == -1:
        return None
    ref_count = len([1 for i in read_base if i in [pileup[2], ".", ","]])
    alt_count = len(read_base) - ref_count

    if my_args.snp_type == 'dbsnp' \
            and (alt_count < 2 or ref_count < 2 or len(read_base) < my_args.RD):
        return None
    return GoldenHetero(pileup[0], pileup[1],  ref_count, alt_count)


def window_data_one_chromosome(center_pos, gh_list, lamb, gh_list_end,
            current_pos):
    """
    window_data from one chromosome data
    :type current_pos: int
    :param current_pos: current handling position.
    :param gh_list_end:
    :rtype: BigForewordList
    :type lamb: int
    :type gh_list: BigForewordList
    :type center_pos: int
    :param center_pos: pos of data in q3
    :param data_list: data after GH in q2(GoldenHetero)
    :param lamb: half of window width
    :return:
        None can not window yet.
        others data windowed
    """
    if gh_list.is_empty(): # is empty
        if gh_list_end:
            return []
        else:
            if current_pos - center_pos > lamb:
                return []
            return None

    # now data_list is not empty, looking for left edge
    if gh_list.get_last_element().pos < center_pos - lamb:
        if gh_list_end:
            return []
        if current_pos - center_pos > lamb:
            return []
        return None
    counter = 0  # type: int
    len_data_list = gh_list.len()
    while 1:
        if counter >= len_data_list:
            if gh_list_end:
                return []
            else:
                if current_pos - center_pos > lamb:
                    return []
                return None
        if gh_list.get_current_element().pos >= center_pos - lamb:
            left_edge = gh_list.get_current_element()
            if left_edge.pos > center_pos + lamb:
                return []
            break
        else:
            gh_list.move_foreword()
        counter += 1

    # right edge
    if gh_list.get_current_element().pos > center_pos + lamb:
        return []
    counter = 0  # type: int
    while 1:
        if counter == len_data_list:
            return []
        if gh_list.get_element(-1 - counter).pos <= center_pos + lamb:
            if gh_list_end or counter > 0:
                right_edge = gh_list.get_element(-1 - counter)
                break
            else:
                if current_pos - center_pos > lamb:
                    right_edge = gh_list.get_last_element()
                    break
                else:
                    return None
        else:
            counter += 1

    # double check edges
    if int(right_edge.pos) < int(left_edge.pos):
        return []
    return gh_list.filter(lambda x: left_edge.pos <= x.pos <= right_edge.pos)


def differential(my_args, gh_list, rel_pileup, queue, worker_id, gh_list_end,
        current_pos):
    """
    window the data in gh_list, send result to queue.
    :type my_format: str
    :param my_format:
    :type gh_list_end: bool
    :param gh_list_end:
    :type worker_id: int
    :param worker_id:
    :param artifact:
    :param default_bias:
    :param queue:
    :type current_pos: int
    :param current_pos: current data pos from source
    :type rel_pileup: PileupDataStruct
    :param rel_pileup:
    :param gh_list: data queue to be windowed
    :type gh_list: BigForewordList
    :param lamb: half of window width
    :type lamb: int
    :return:
    True data is windowed, send result to W
    False can not window yet
    """

    tracked_data = window_data_one_chromosome(rel_pileup.pos, gh_list,
        my_args.lamb, gh_list_end, current_pos)
    
    if tracked_data is None:
        return False

    if tracked_data:
        bias = bias_estimator(rel_pileup.pos, tracked_data, my_args.lamb)
    else:
        bias = my_args.bias

    # <INFO, NB> Calculation of basic values for GQ and PL
    # rr = Sequencing  Noise (rr_u = log10(rr))
    # ra = Amplification Artefact/Error (ra_u = log10(ra))
    # rm = Heterozygous SNV (rm_u = log10(rm))
    # mm = Homozygous SNV (mm_u = log10(mm))
    # <DONE, NB>
    # First two values, rr and ra, are added to get a single likelihood for REF/REF
    log_lh = sc_caller(rel_pileup, bias, my_args.null)

    outline = OutLineStruct(rel_pileup.name, rel_pileup.pos, 
        rel_pileup.ref_base, rel_pileup.var_base, rel_pileup.ref_num,
        rel_pileup.var_num, log_lh, rel_pileup.so, rel_pileup.var_all, bias, 
        rel_pileup.gt_num, rel_pileup.gt_class)
 
    queue.put(QueueData(worker_id, my_args.worker_type, outline=outline),
        block=False)
    return True


def so_generator(bulk_pileup_source, my_args):
    """
    :return:
        pos is not in bulk or bulk has no data: noCoverageInControl
        pos in bulk does not have enough depth: shallow
        pos is indel in bulk: indel
        pos is ref in bulk: ref
        pos is var in bulk: germline
    """
    should_read = True
    pos = -1
    while 1:
        if should_read:
            bulk_pileup_list = next(bulk_pileup_source)
            if bulk_pileup_list is None:
                while 1:
                    pos = yield "noCoverage"

        if int(bulk_pileup_list[1]) < pos:
            should_read = True
            continue
        elif int(bulk_pileup_list[1]) > pos:
            should_read = False
            pos = yield "noCoverage"
        else:
            should_read = True
            pos = yield read_bulk_mpileup(bulk_pileup_list, my_args)


def read_bulk_mpileup(pileup, my_args):
    if "*" in pileup[4] or "-" in pileup[4] or "+" in pileup[4]:
        return "indel"
    pileup[4] = remove_head_end(pileup)
    MQ_filter = [ord(i) - PHREDSCORE >= my_args.bulk_min_mapq for i in pileup[6]]
    for i, base in enumerate(pileup[4]):
        if base == "I":
            MQ_filter.insert(i, True)
    read_base = "".join(compress(pileup[4], MQ_filter))
    # Coverage in bulk/dbsnp too shallow
    if len(read_base) < my_args.bulk_min_depth:
        return "shallow({})".format(len(read_base))

    var_counts = sorted([(read_base.count(i), i) for i in  ["A", "C", "G", "T"]])
    # Only reference supporintg reads in bulk normal (smaller than threshold)
    if var_counts[-1][0] < my_args.bulk_min_var:
        return 'ref'
    # Variant supporting reads in bulk normal -> germline mutation
    if var_counts[-2][0] > 0:
        return 'var({}{},{}{})' \
            .format(var_counts[-1][0], var_counts[-1][1], 
                var_counts[-2][0], var_counts[-2][1])
    else:
        return 'var({}{})'.format(*var_counts[-1])


def bias_estimator(pos, tracked_data, lamb):
    K = np.zeros((len(tracked_data), 2))

    # Equation 2
    for i, data in enumerate(tracked_data):
        be_tmp1 = float(data.ref_num + data.var_num)
        be_tmp2 = max(data.ref_num, data.var_num) / be_tmp1

        t  = pos - data.pos
        if -lamb <= t <= lamb:
            # Equation 3
            be_tmp = 0.75 * (1 - (float(t) / lamb) ** 2)
            # <DONE, NB>
            K[i, 0] = be_tmp * be_tmp2
            K[i, 1] = be_tmp
            # <BUG, Original> Adding an additional factor be_tmp1 to 
            #   nominator and denominator
            # K[i, 0] = be_tmp * be_tmp1 * be_tmp2
            # K[i, 1] = be_tmp * be_tmp1

    K_sum = K.sum(axis=0)
    # Nadaraya-Watson kernel-weighted average
    if K_sum[1] == 0:
        return 1
    else:
        return K_sum[0] / K_sum[1]


def sc_caller(sc_candidate, sc_bias, min_frac):
    """ Calculate RR, RA, RM, MM from data in q3 and bias.
    """
    if "," in sc_candidate.var_base:
        ref, mut = sc_candidate.var_base.split(",")
    else:
        ref = sc_candidate.ref_base
        mut = sc_candidate.var_base

    # <DONE, NB>
    f_h0 = 0.125 * sc_bias
    # <BUG, Original> Ignoring the bias
    # f_h0 = 0.125

    f_h1 = sc_bias
    f_h2 = 1
    if sc_candidate.ref_num > sc_candidate.var_num:
        f_h1 = 1 - f_h1
        # <TODO, NB> Confirm if Eq. 6 from the manuscript is wrong
        # f_h2 = 0
        # <BUG, Original> Ignoring the 1 - F_K(\theta) if ref count > alt count
        f_h2 = 1
  
    f_v = np.array([min_frac, f_h0, f_h1, f_h2])
    lh_m = np.zeros((len(sc_candidate.read_infos), 4))
    for idx, bp in enumerate(sc_candidate.read_infos):
        # Equation (7), calculate lg Probability value
        e = 10 ** (-bp[1] / 10.0)
        if bp[0] in [',', '.', ref]:
            lh_m[idx] = f_v * e / 3 + (1 - f_v) * (1 - e)
        elif bp[0] == mut:
            lh_m[idx] = (1 - f_v) * e / 3 + f_v * (1 - e)
        else:
            lh_m[idx] = e / 3

    log_lh = np.log10(lh_m).sum(axis=0)
    # Add noise and artifact likelihoods in log space for REF/REF log-likelihood
    if log_lh[0] > log_lh[1]:
        log_lh_wt = log_lh[0] + np.log10(1 + 10 ** (log_lh[1] - log_lh[0]))
    else:
        log_lh_wt = log_lh[1] + np.log10(1 + 10 ** (log_lh[0] - log_lh[1]))

    return np.array([log_lh_wt, log_lh[2], log_lh[3]])


def get_my_filename(output, suffix):
    file_name_base = os.path.splitext(os.path.basename(output))[0]
    file_name = '{}{}'.format(file_name_base, suffix)
    return os.path.join(os.path.dirname(output), file_name)


# def calculate_eta(list_var_buf, list_var_tag, alpha=0.05):
#     allele_nums = []
#     bias_list = []
#     for i, var_tag in enumerate(list_var_tag):
#         if not var_tag:
#             break
#         for queue_data in list_var_buf[i]:
#             allele_nums.append(queue_data.outline.total_num)
#             bias_list.append(queue_data.outline.bias) 
    
#     allele_nums = allele_nums[:CUTOFFNUM]
#     bias_list = bias_list[:CUTOFFNUM]

#     LLR = np.zeros(len(allele_nums))
#     for i, allele_num in enumerate(allele_nums):
#         f_artifact = 0.125 * bias_list[i] / 0.5
#         # bin (depth(number of trials),prob_success)
#         alt = np.random.binomial(allele_num, f_artifact)
#         ref = allele_num - alt
#         L_H0 = (1 - f_artifact) ** ref * f_artifact ** alt
#         # random select major allele
#         major = np.random.randint(2)
#         if major == 0:
#             f_true = 0.5 * 0.5 / bias_list[i]
#         else:
#             f_true = 0.5 * bias_list[i] / 0.5
#         L_H1 = (1 - f_true) ** ref * f_true ** alt
#         # if L_H0/L_H1 is 0, assign a very small value
#         if L_H0 == 0:
#             L_H0 = 10 ** -100
#         if L_H1 == 0:
#             L_H1 = 10 ** -100

#         ratio = L_H0 / L_H1
#         if ratio != 0:
#             LLR[i] = -np.log10(ratio)
#         else:
#             LLR[i] = 10000

#     co = np.percentile(LLR, int(100 - alpha * 100))
#     result = 10 ** -co

#     return result


def calculate_eta(list_var_buf, list_var_tag, alpha=0.05):
    idx = 0
    LLR = np.zeros(CUTOFFNUM)

    for i, var_tag in enumerate(list_var_tag):
        if not var_tag or idx >= CUTOFFNUM:
            break

        for queue_data in list_var_buf[i]:
            read_count = queue_data.outline.total_num
            bias = queue_data.outline.bias

            f_artifact = 0.125 * bias
            # Sample read counts under H_0 from a binomial
            alt = np.random.binomial(read_count, f_artifact)
            ref = read_count - alt
            L_H0 = ref * np.log10(1 - f_artifact) + alt * np.log10(f_artifact)
        
            # Random select major allele
            if np.random.random() < 0.5:
                f_mut = 1 - bias
            else:
                f_mut = bias
            L_H1 = ref * np.log10(1 - f_mut) + alt * np.log10(f_mut)

            try:
                LLR[idx] = L_H0 - L_H1
            except IndexError:
                break
            else:
                idx += 1

    return 10 ** np.percentile(LLR[:idx], alpha * 100)


def write_result(queue, args, name):
    """
    Receive the data in queue, organize it, and write the result file.
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
    current_handling_worker = 0

    with open(body_file, "a") as fp_out, \
            open(reasoning_file, "a") as fp_reasoning, \
            open(eta_file, "a") as fp_eta:
        while 1:
            queue_data = queue.get(block=True)
            w_id = queue_data.worker_id

            # from main
            if queue_data.mode == WORKEXIT:
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
                        eta = calculate_eta(list_var_buf, list_var_tag,
                            args.lrt_alpha)
                        for i, var_buf in enumerate(list_var_buf):
                            if var_buf:
                                write_do([j.outline for j in var_buf],
                                    queue_data.work_type, i, fp_out, eta,
                                    fp_reasoning, args)
                            else:
                                logging.info("worker{} write len = 0".format(i))
                    else:
                        logging.info("Nothing need to be written.")
                fp_eta.write("##contig=<ID={},eta={}>\n".format(name, eta))
                break

            elif queue_data.work_type == WORKVAR \
                    or queue_data.work_type == WORKVCF:

                # record tag and data
                if queue_data.mode == WORKDONE:
                    list_var_tag[w_id] = True
                    logging.info("chr{} worker{} done".format(name, w_id))
                    # Try to calculate eta, write data
                    if get_current_cutoff_num(list_var_buf, list_var_tag) \
                            >= CUTOFFNUM and eta == -1:
                        eta = calculate_eta(list_var_buf, list_var_tag,
                            args.lrt_alpha)
                    if eta != -1:
                        for i, var_buf in enumerate(list_var_buf):
                            if i < current_handling_worker:
                                continue
                            elif list_var_tag[i]:
                                if var_buf:
                                    write_do([j.outline for j in var_buf],
                                        queue_data.work_type, i, fp_out, eta,
                                        fp_reasoning, args)
                                else:
                                    logging.info("Chr{} worker{} write len=0" \
                                        .format(name, i))
                                current_handling_worker += 1
                            else:
                                break
                else:
                    list_var_buf[w_id].append(queue_data)

            else:
                list_coverage_buf[w_id].append(queue_data.coverage)

    if args.bulk == "" or queue_data.work_type == WORKVCF:
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
            fp_out.write('\n'.join([i.get_varcall_str() for i in data_list]))
            fp_out.write('\n')
        if my_args.bulk != "" and data_list:
            fp_reasoning.write('\n'.join(
                [i.get_reason_str() for i in data_list]))
            fp_reasoning.write('\n')
    elif work_type == WORKVCF:
        fp_out.write('\n'.join(
            [i.get_vcf_str(eta, my_args.minvar, my_args.minvarfrac) \
                for i in data_list]))
        fp_out.write('\n')

    logging.debug('worker{} write len={}'.format(worker_id, len(data_list)))


def control(my_args, snp_pos_subset, queue, name, head, stop, worker_id):
    """
    :param stop: Unexpanded region
    :type my_args: object
    :type snp_pos_subset: list
    :type name: str
    :type queue: object (main queue)
    :type worker_id: int
    :type stop: int
    :type head: int
    """
    logging.info('worker {} start! name={} head={} tail={} len={} engine={}' \
        .format(worker_id, name, head, stop, stop - head, my_args.engine))

    snp_list = BigForewordList(snp_pos_subset)
    rel_pileups = BigForewordList([])
    gh_list = BigForewordList([])
    gh_list_end = False
    
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

    # Pileup of single cell [tumor] sample
    pileup_source = io.data_generator(my_args, name, my_start, my_stop, False)

    if my_args.bulk != '':
        # Pileup of bulk normal sample
        bulk_pileup_source = io.data_generator(my_args, name, head, stop, True)
        so_source = so_generator(bulk_pileup_source, my_args)
        next(so_source)

    if my_args.debug:
        print('w{:<3d}\t0/{:<6}'.format(worker_id, total_work))

    coverage_list = []  # list of [name str, head int, stop int]
    counter = 0
    main_index = 0
    while 1:
        pileup = next(pileup_source)
        if pileup is None:
            if my_args.coverage and len(coverage_list) > 0:
                done_data = QueueData(worker_id, WORKCOVERAGE, mode=WORKDONE)
                queue.put(done_data, block=False)
                del coverage_list[0]
            break

        # append the relevant pileups into rel_pileups
        if head <= pileup[1] < stop:
            # Check if pileup is relevant/passes filters
            if my_args.worker_type == WORKVAR:
                rel_pileup = read_mpileup(copy.copy(pileup),
                    my_args.minvar, my_args.mapq, my_args.min_depth, False) 
            else:
                rel_pileup = read_mpileup_vcf(copy.copy(pileup), my_args)

            if rel_pileup and rel_pileup != -1:
                if my_args.bulk != '':
                    # calculate so
                    result = so_source.send(rel_pileup.pos)
                    rel_pileup.so = result
                else:
                    rel_pileup.so = ''
                rel_pileups.append(rel_pileup)

            # calculate coverage
            if my_args.coverage:
                if rel_pileup and (my_args.bulk == '' \
                        or (rel_pileup == -1 or (rel_pileup != -1 \
                            and (rel_pileup.so == '' or rel_pileup.so == 'ref' \
                                or rel_pileup.so.startswith('var'))))):

                    if not coverage_list:  # is empty
                        coverage_list.append([name, pileup[1], pileup[1]])
                    elif pileup[1] == coverage_list[-1][2] + 1:
                        coverage_list[-1][2] = pileup[1]
                    else:
                        coverage_list.append([name, pileup[1], pileup[1]])
                        new_coverage = coverage_list.pop(0)
                        coverage_data = QueueData(worker_id, WORKCOVERAGE, 
                            coverage=new_coverage)
                        queue.put(coverage_data, block=False)

        # handle the head data in rel_pileups
        if not rel_pileups.is_empty():
            diff_success = differential(my_args, gh_list, 
                rel_pileups.get_current_element(), queue, worker_id,
                    gh_list_end, pileup[1])
            if diff_success:
                rel_pileups.move_foreword()

        if my_args.debug:
            counter += 1
            if counter == 10000:
                counter = 0
                main_index += 1
                print('w{:<3d} {:>6}/{:<6}'\
                    .format(worker_id, main_index * 10000, total_work))
        # Screen SNP list
        if not gh_list_end:
            while not snp_list.is_empty():
                if snp_list.get_current_element() < pileup[1]:
                   snp_list.move_foreword()
                else:
                    break
            # current SC tumor pileup has corresponding entry in SNP list
            if pileup[1] == snp_list.get_current_element():
                golden_hetero = get_golden_hetero(my_args, pileup, worker_id)
                if golden_hetero is not None:
                    gh_list.append(golden_hetero)
            if snp_list.is_empty() \
                    or pileup[1] >= snp_list.get_last_element():
                gh_list_end = True

    while not rel_pileups.is_empty():
        diff_success = differential(my_args, gh_list,
            rel_pileups.get_current_element(), queue, worker_id, True, my_stop)
        if diff_success:
            rel_pileups.move_foreword()

    done_data = QueueData(worker_id, my_args.worker_type, mode=WORKDONE)
    queue.put(done_data, block=False)

    if my_args.debug:
        print('w{:<3d} I am done!'.format(worker_id))
    logging.info('worker {} stop!'.format(worker_id))


def main(my_args):
    # Init logging
    log_file = get_my_filename(my_args.output, '_sccaller_{:0>2d}to{:0>2d}.log' \
        .format(my_args.head, my_args.tail))
    logging.basicConfig(filename=log_file, level=logging.DEBUG,
        format=LOG_FORMAT, filemode='w')

    logging.info('Welcome to SCcaller v{}'.format(VERSION))

    if my_args.debug:
        my_args.cpu_num = 1
        my_args.work_num = 1
        if my_args.seed == -1:
            my_args.seed = 1
        print('\nRunning SCcaller v{} in debugging mode'.format(VERSION))
        print('parsing vcf...')

    if my_args.seed == -1:
        my_args.seed = np.random.randint(0, 2 ** 32 - 1)
    np.random.seed(my_args.seed)

    logging.info('parsing SNP vcf...')
    snp_info = io.parse_snp_info(my_args)

    if my_args.debug:
        print('parsing fasta...')
    logging.info('parsing reference fasta...')
    fasta_info = io.parse_fasta(my_args.fasta)

    if my_args.tail == -1 or my_args.tail > len(fasta_info):
        my_args.tail = len(fasta_info)

    # Write arguments to log file for reproducability
    logging.info('args: {}'.format(my_args))

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
        logging.info('loading vcf...')
        snp_pos = np.array(io.load_snp_pos(name, snp_info, my_args.snp_in))

        if my_args.debug:
            print('SCcaller v2.0.0 is handling chromosome {}...'.format(name))
        logging.debug('name={}; #snps={}'.format(name, len(snp_pos)))

        # Calculate the amount of tasks for each process
        no_bases = stop - start
        step = int(np.ceil(no_bases / float(my_args.work_num)))
        if step < my_args.lamb:
            logging.info('work_num={} is too large for chr={}. Using 1 instead.' \
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
                    err_str = 'lamb={} is too large. Fasta: name={} ({} - {}).' \
                        .format(my_args.lamb, name, start, stop)
                    logging.critical(err_str)
                    raise RuntimeError(err_str)

            snp_pos_subset = [i for i in snp_pos \
                if head - my_args.lamb <= i <= tail + my_args.lamb]
            
            if my_args.work_num > 1:
                process_pool.apply_async(control, 
                    (my_args, snp_pos_subset, queue, name, head, tail, woker_id))
            else:
                control(my_args, snp_pos_subset, queue, name, head, tail, woker_id)
                
        if my_args.work_num > 1:
            process_pool.close()
            process_pool.join()

        # Exit the W process
        queue.put(QueueData(0, my_args.worker_type, mode=WORKEXIT), block=True)
        proc_write_result.join()

    io.write_vcf(my_args, VERSION)
    logging.info('W quit. All done.')


if __name__ == '__main__':
    start = datetime.now()
    args = parse_args()
    print('Start running SCcaller v{} at {:%Y%m%d_%H:%M:%S}' \
        .format(VERSION, start))
    main(args)
    print('Stop running SCcaller (duration={} s)' \
        .format(str(datetime.now() - start)))