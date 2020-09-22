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
import argparse
import os
from subprocess import Popen, PIPE
import multiprocessing as mp
import re
from collections import Counter
from itertools import compress
import logging
import time
from functools import wraps
import pysam  # 0.15.1
import copy
import random
import numpy as np

from BigForewordList import BigForewordList


if float(sys.version[:3]) != 2.7:
    print("CRITICAL: Python version must be 2.7!\n")
    sys.exit(1)


# regular expression of indels
INSERTION = "\+[0-9]+[ACGTNacgtn]+"
DELETION = "\-[0-9]+[ACGTNacgtn]+"
INDEL = "\+[0-9]+[ACGTNacgtn]+|\-[0-9]+[ACGTNacgtn]+"
PHREDSCORE = 33  #
INVALIDNUM = -1

VARCALLFORMAT = "bed"
VCFFORMAT = "vcf"
ENGINESAMTOOLS = "samtools"
ENGINEPYSAM = "pysam"
MULTIPLEGENOTYPE = "multiple-genotype"
NOTENOUGHVARIANTS = "No.variants<"

WORKVAR = 1
WORKCOVERAGE = 2
WORKVCF = 3
STOPSIGN = -1

DEFAULTBASEQUALITY = 30
CUTOFFNUM = 100000


V_ROWS = 26
LOG_FORMAT = "%(asctime)s %(levelname)-8s %(process)d [%(funcName)s]" \
    "%(message)s(%(filename)s line%(lineno)d)"
VERSION = "2.0.0_NB"



class DataInQ7:
    def __init__(self, name, time_float):
        # type: (str, float)
        self.name = name 
        self.time_float = time_float


def send_name_time(should_send, q7):
    def send(f):
        @wraps(f)
        def decorated(*args, **kwargs):
            if should_send:
                start_ = time.time()
            ret_ = f(*args, **kwargs)
            if should_send:
                q7.put(DataInQ7(f.__name__, time.time() - start_), block=False)
            return ret_

        return decorated

    return send


# @send_name_time(should_send=should_analyze, q7=q7)
class VcfInfo:
    def __init__(self, gt, ad, bi, gq, pl, qual, genotype_num, genotype_class):
        """

        :type qual: str
        :type pl: str
        :type gq: str
        :type bi: float
        :type genotype_num: int
        :type ad: str
        :type so: str
        :type gt: str
        """
        self.gt = gt

        self.ad = ad
        self.bi = bi
        self.gq = gq
        self.pl = pl
        self.qual = qual
        self.genotype_num = genotype_num
        self.genotype_class = genotype_class


class OutLineStruct:
    def __init__(self, name, pos, ref, var, ref_num, var_num, bias,
                sn, ae, he, ho, vcf_info, so, var_all):
        # type: (str, int, str, str, int, int, float, float, float, float, float, VcfInfo) -> OutLineStruct
        self.name = name
        self.pos = pos
        self.ref = ref
        self.var = var
        self.ref_num = ref_num  # 6
        self.var_num = var_num  # 7
        self.bias = bias
        self.sn = sn  # 10
        self.ae = ae  # 11
        self.he = he  # 12
        self.ho = ho  # 13
        self.so = so
        self.vcf_info = vcf_info
        self.var_all = var_all


    def __str__(self):
        out_str = '{}\t{}\t{}\t{}\t{}\t{}\t{}' \
            .format(self.name, self.pos, self.ref, self.var, self.ref_num,
                self.var_num, self.so)
        return out_str


class DataInQ5:
    def __init__(self, outline_struct, worker_id, work_type, coverage, done=False):
        self.outlineStruct = outline_struct  # type: OutLineStruct
        self.worker_id = worker_id  # type: int
        self.work_type = work_type
        self.coverage = coverage  # type: list  # [name str, start int, stop int]
        self.done = done


    def __str__(self):
        out_str = 'DataInQ5 Object:\n\tID: {}\n\tType: {}\n\tOutlineStruct: {}\n' \
                '\tCoverage: {}\n\tDone: {}' \
            .format(self.worker_id, self.work_type, self.outlineStruct, 
                self.coverage, self.done)
        return out_str


    def is_done(self):
        return self.done


    def set_done(self):
        self.done = True


class DataInQ2:
    def __init__(self, name, coordinate_1, refcount, altcount):
        # type: (str, int, int, int) -> DataInQ2
        self.name = name  # type: str
        self.coordinate_1 = coordinate_1  # type: int
        self.refcount = refcount  # type: int
        self.altcount = altcount  # type: int


    def __str__(self):
        out_str = '{}:{}\tAD:{},{}' \
            .format(self.name, self.coordinate_1, self.refcount, self.altcount)
        return out_str



class ReadInfo:
    def __init__(self, read_base, base_quality):
        # type: (str, str) -> None
        self.read_base = read_base  # type: str
        self.base_quality = base_quality  # type: int


    def __str__(self):
        return '{} (BQ: {})'.format(self.read_base, self.base_quality)


class PileupDataStruct:
    def __init__(self, name, coordinate_1, reference_base, variant,
                reference_allele_num, variant_allele_num, genotype_num,
                read_info_list, so=None, genotype_class=-1, variant_all=""):
        # type: (str, int, str, str, int, int, int, list[ReadInfo], str, int) -> object
        self.name = name  # type: str
        self.coordinate_1 = coordinate_1  # type: int
        self.reference_base = reference_base  # type: str
        self.variant = variant  # type: str
        self.reference_allele_num = reference_allele_num  # type: int
        self.variant_allele_num = variant_allele_num  # type: int
        self.genotype_num = genotype_num  # type: int
        self.read_info_list = read_info_list  # type: list[ReadInfo]
        self.so = so  # type: str  # reasoning  bulk ref -- "True"   bulk var -- "False"  else -- "NA"
        self.genotype_class = genotype_class  # type: int  # 0 unknown 1:0,0 0,1 1,1      2:1,1 1,2 2,2
        self.variant_all = variant_all


    def __str__(self):
        out_str = '{}:{}\tref:{}({}),alt:{}({})' \
            .format(self.name, self.coordinate_1, self.reference_base,
                self.reference_allele_num, self.variant, self.variant_allele_num)
        return out_str



def my_print(msg, name='main'):
    pass
    print('{}'.format(msg))


# @send_name_time(should_send=should_analyze, q7=q7)
def parse_args():
    """
    handle the parameters.
    :return: argparse.Namespace
    """

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
    parser.add_argument("-e", "--engine", type=str,
        choices=[ENGINEPYSAM, ENGINESAMTOOLS], default=ENGINEPYSAM,
        help="Pileup engine. Default: pysam")
    parser.add_argument("-w", "--work_num", type=int, default=100,
        help="# splits per chromosome for multi-process computing. Default: 100.")
    parser.add_argument("-n", "--cpu_num", type=int, default=1,
        help="Num. processes. Default: 1")
    parser.add_argument("--head", type=int, default=1,
        help="First chromosome as sorted as in fasta file to analyze (1-based). "
            "Default: 1")
    parser.add_argument("--tail", type=int, default=-1,
        help="Last chromosome as sorted as in fasta file to analyze (1-based). "
        "Default: -1")
    parser.add_argument("--format", type=str, choices=[VARCALLFORMAT, VCFFORMAT],
        default=VCFFORMAT, help="Output file format. Default: vcf")
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
        raise IOError("{}: When choosing dbsnp, --bulk bulk_bamfile is required." \
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
                for j in xrange(target_len -1, 0, -1):
                    if target[j] != "N":
                        break
                tail = j + 1 

            my_result.append([name, head, tail])

    return my_result


def parse_indel(indel_str, indel_list_out):
    # type: (str, list) -> int
    """
    parse indel string, return the indel string and length of indel string
    :param indel_str:indel string (the '+' and '-' in indel string doesn't affect the result)
    :param indel_list_out:indel string (should be clear before)
    :return:
    >0 return indel string length
    =0 in case of indel and svn
    <0 in case of invalid format of indel
    eg:
        rr = []
        ret = parse_indel("+2AAC",rr)   #it will print:
        print ret                       #0
        print rr                        #['+', '2', 'A', 'A']

        rr = []
        ret = parse_indel("+2AA",rr)    #it will print:
        print ret                       #4
        print rr                        #['+', '2', 'A', 'A']
    """
    i = 0  # type: int
    j = 0  # type: int
    len_indel_str = len(indel_str)
    for i in xrange(len_indel_str):
        if indel_str[i] == "+" or indel_str[i] == "-":
            continue
        else:
            break
    # more than 1 '+' or '-'
    if i > 1:
        return -1
    buf = indel_str[i:]
    len_buf = len(buf)
    for j in xrange(len_buf):
        if buf[j].isalpha():
            break

    if len_indel_str - i - j > int(buf[:j]):
        indel_list_out.extend(indel_str[0:i + j + int(buf[:j])])
        return 0
    elif len_indel_str - i - j < int(buf[:j]):
        return -2
    indel_list_out.extend(indel_str[0:i + j + int(buf[:j])])
    return int(buf[:j]) + j + i


# @send_name_time(should_send=should_analyze, q7=q7)
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
    for i in xrange(len(str1) - 1):
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


def compress_read_base(read_base, my_filter):
    """
    compress read base but leave the I (indel)
    :type read_base: str
    :type my_filter: list
    """
    counter = 0
    for i in read_base:
        if i == "I":
            my_filter.insert(counter, True)
        counter += 1
    return "".join(list(compress(read_base, my_filter)))


def get_reference_variant_allele_num(read_base):
    v_num = 0
    r_num = 0
    for i, base in enumerate(read_base):
        if base in [".", ","]:
            r_num += 1
        if base in ["A", "T", "C", "G"]:
            v_num += 1
        if i > 0 and base in ["I", "i"] and read_base[i - 1] in [".", ","]:
            r_num -= 1
            v_num += 1
    return [r_num, v_num]


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
                for i in xrange(indel_count_list[list_index])]
            for list_index in xrange(len(indel_count_list))]
    for i in tmp:
        read_base_list.extend(i)
    return read_base_list


def read_mpileup(mpileup_list, rm_minvar, rm_minmap_q, rm_mindepth, is_gh, worker_id):
    """
    screenfor mutations, filter out alignments that do not meet the requirements of mapQ, and convert the data format
    :param is_gh:
    :type is_gh: bool
    :type rm_mindepth: int
    :type rm_minmap_q: int
    :type rm_minvar: int
    :type mpileup_list: list
    :param mpileup_list: 1 line data of mpileup
    :param rm_minvar: variant supporting reads.
    :param rm_minmap_q: minimun mapQ
    :param rm_mindepth: minimun read depth
    :return:
    not empty PileupDataStruct: mutation data
    -1: return reference genome type and for skipping this line, but should be count in coverage
    []: not mutation, should not be count in coverage
    """
    if "*" in mpileup_list[4]:
        return []  # indel read overlaps

    indel_name_list = []
    indel_count_list = []
    rm_indel = ""
    if ("-" in mpileup_list[4]) or ("+" in mpileup_list[4]):
        indel_list = re.findall(INDEL, mpileup_list[4])
        if indel_list:
            tmp = Counter(indel_list).most_common()
            if len(tmp) > 1:
                return []  # multiple different indel calls
            indel_name_list = map(lambda x: x[0], tmp)
            indel_count_list = map(lambda x: x[1], tmp)
            rm_indel = tmp[0][0]
            result = []
            if parse_indel(rm_indel, result) == 0:  # indel reads followed by a SNV read, e.g. ....-2AAAA.....
                return []
            mpileup_list[4] = mpileup_list[4].replace(rm_indel, "I")
    # remove head (^) and tail ($)
    mpileup_list[4] = remove_head_end(mpileup_list)

    # filter out alignments that do not meet the requirements of mapQ
    data_filter = [ord(i) - PHREDSCORE >= rm_minmap_q for i in mpileup_list[6]]
    if not is_gh:
        base_quality = map(lambda x: ord(x) - PHREDSCORE,
                           list(compress(mpileup_list[5], data_filter)))  # type: list

    read_base_with_i = compress_read_base(mpileup_list[4], data_filter)  # type: str
    read_base_with_d = "".join(
        [mpileup_list[4][i] if data_filter[i] else "D" \
            for i in xrange(len(mpileup_list[4]))])  # type: str
    read_base_without_i = re.sub("I", "", read_base_with_i)  # type: str
    if len(read_base_without_i) < rm_mindepth:
        return []

    if is_gh:
        return read_base_without_i

    maxvarp = ["A", "C", "G", "T", "I"]
    maxvar = [read_base_with_i.count(element) for element in maxvarp]
    i_index_max = maxvar.index(max(maxvar))
    if 4 == i_index_max:
        if len(rm_indel) == 0:
            logging.critical("has I but no rm_indel name = {} pos = {}" \
                .format(mpileup_list[0], mpileup_list[1]))
        maxvarp[4] = rm_indel

    if maxvar[i_index_max] < rm_minvar:
        return -1  # return reference genome type and for skipping this line
    r_num, v_num = get_reference_variant_allele_num(read_base_with_d)
    num_of_i = len(read_base_with_i) - len(read_base_without_i)
    read_base_list_final = rebuild_read_base_list(read_base_without_i, indel_name_list, indel_count_list)
    base_quality.extend([DEFAULTBASEQUALITY for i in xrange(num_of_i)])  # type: list

    rm_out = PileupDataStruct(mpileup_list[0],  # name str
                              mpileup_list[1],  # pos  int
                              mpileup_list[2],  # reference_base  str
                              maxvarp[i_index_max],  # variant  str
                              r_num,  # reference_allele_num  int
                              v_num,  # variant_allele_num  int
                              1,  # genotype_num  int
                              map(lambda x: ReadInfo(read_base_list_final[x], base_quality[x]),
                                  xrange(len(read_base_list_final))),
                              1)  # genotype class
    return rm_out


def choose(candidate_index_list, current_basis, choice_out, basis_filter):
    """

    choose indexs from the candidate_index_list according to current_basis
    :param candidate_index_list: the length equal to the number of 'True's in basis_filter
    :type values: list[int]
    """
    values = list(compress(current_basis, basis_filter))  # compressed basis
    tmp = Counter(values).most_common()
    value_list = map(lambda x: x[0], tmp)  # value
    count_list = map(lambda x: x[1], tmp)  # count of value
    sorted_value_list = copy.copy(value_list)
    sorted_value_list.sort(reverse=True)  # descending sorted values
    candidate = []
    if len(choice_out) >= 2:
        return candidate
    if count_list[value_list.index(sorted_value_list[0])] == 1:  # Only one of the maximum
        choice_out.append(candidate_index_list[values.index(sorted_value_list[0])])
        basis_filter[candidate_index_list[values.index(sorted_value_list[0])]] = False

        if len(choice_out) < 2:
            if count_list[value_list.index(sorted_value_list[1])] == 1:  # only one second largest value
                choice_out.append(candidate_index_list[values.index(sorted_value_list[1])])
                basis_filter[candidate_index_list[values.index(sorted_value_list[1])]] = False
            else:
                candidate.extend(map(lambda y: candidate_index_list[y],
                                     filter(lambda x: values[x] == sorted_value_list[1], xrange(len(values)))))

    else:
        candidate.extend(map(lambda y: candidate_index_list[y],
                             filter(lambda x: values[x] == sorted_value_list[0], xrange(len(values)))))
    for i in xrange(len(basis_filter)):
        if i not in candidate:
            basis_filter[i] = False
    return candidate


def choose_random(candidate, num):
    result = []
    for i in xrange(num):
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
        if read_base != "I":
            tmp_list.append(read_base)
        else:
            tmp_list[-1] = indel_list[index]
            index += 1
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
            indel_name_list = map(lambda x: x[0], indel_count_list)
            indel_count_list = map(lambda x: x[1], indel_count_list)
            for i in indel_name_list:
                pileup[4] = pileup[4].replace(i, "I")

    # remove head(^) and tail($)
    pileup[4] = remove_head_end(pileup)
    
    # screen according to mapQ
    map_q_filter = [ord(i) - PHREDSCORE >= my_args.mapq for i in pileup[6]]
    base_q = map(lambda x: ord(x) - PHREDSCORE, list(compress(pileup[5], map_q_filter)))
    map_q = map(lambda x: ord(x) - PHREDSCORE, list(compress(pileup[6], map_q_filter)))

    readbase_list_with_i = split_readbase(pileup[4], indel_list)
    read_bases_filtered = list(compress(readbase_list_with_i, map_q_filter))

    i_filter = ["+" not in i and "-" not in i for i in read_bases_filtered]
    read_base_list_without_i = list(compress(read_bases_filtered, i_filter))

    # less min_depth
    if len(read_bases_filtered) < my_args.min_depth:
        return []
    
    # Count occurence of all possible mutations/bases
    var_names = ["A", "C", "G", "T"]
    var_counts = [read_base_list_without_i.count(i) for i in var_names]
    # Filter out mutations/bases that were not observed
    name_filter = [i > 0 for i in var_counts]
    var_names = list(compress(var_names, name_filter))
    var_counts = list(compress(var_counts, name_filter))

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
        var_MQ.append(sum(list(compress(map_q, q_filter))))
        var_BQ.append(sum(list(compress(base_q, q_filter))))
    
    geno_class = get_geno_class(read_bases_filtered, map_q, base_q, var_counts,
        var_MQ, var_BQ)

    variant_str, r_num, v_num = get_variant_info(var_names, 
        [var_counts, var_MQ, var_BQ], read_bases_filtered, geno_class)
    sort_by_name(var_names, var_counts, var_MQ, var_BQ)

    var_all_str = ",".join(var_names)
    read_info_list = [ReadInfo(j, base_q[i]) \
        for i, j in enumerate(readbase_list_with_i)]
    
    return PileupDataStruct(pileup[0], pileup[1], pileup[2], variant_str, r_num,
        v_num, len(var_names), read_info_list, None, geno_class, var_all_str)


def sort_by_name(var_names, v_count_list, v_map_q_list, v_base_q_list):
    len_names = len(var_names)
    for i in xrange(len_names):
        for j in xrange(len_names):
            index = len_names - 1 - j - i
            if index == 0:
                break
            if v_count_list[index] > v_count_list[index - 1] \
                    or v_map_q_list[index] > v_map_q_list[index - 1] \
                    or v_base_q_list[index] >= v_base_q_list[index - 1]:
                for ex_list in [var_names, v_count_list, v_map_q_list, v_base_q_list]:
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
        sum_ref_map_q = sum(list(compress(map_q, read_base_filter)))
        v_map_q = var_MQ[var_counts.index(second_num)]
        if sum_ref_map_q > v_map_q:
            return 1
        elif sum_ref_map_q < v_map_q:
            return 2
        else:
            sum_ref_base_q = sum(list(compress(base_q, read_base_filter)))
            v_base_q = var_BQ[var_counts.index(second_num)]
            if sum_ref_base_q >= v_base_q:
                return 1
            return 2


# @send_name_time(should_send=should_analyze, q7=q7)
def golden_hetero(pileup, min_map_q, min_depth, snp_type, read_depth, worker_id):
    """
    handle the data after vcf screen, store the result in q2
    :param pileup:
    :param min_depth: Known heterogous SNP call required read depths. RD
    :param snp_type: Known snp type input for known heterogous SNP call
    :param read_depth: minimum reads
    :param min_map_q: mapQ
    :return: None  should skip
    """
    read_base = read_mpileup(pileup, 0, min_map_q, min_depth, True, worker_id)  # type: str
    if not read_base or read_base == -1:
        return None
    ref_count = len([1 for i in read_base if i in [pileup[2], ".", ","]])
    alt_count = len(read_base) - ref_count

    if snp_type == 'dbsnp' \
            and (alt_count < 2 or ref_count < 2 or len(read_base) < read_depth):
        return None
    return DataInQ2(pileup[0], pileup[1],  ref_count, alt_count)


def window_data_one_chromosome(center_pos, q2_list, lamb, is_list_ended,
            current_pos, worker_id):
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
    :param data_list: data after GH in q2(DataInQ2)
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
    if q2_list.get_last_element().coordinate_1 < center_pos - lamb:
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
        if q2_list.get_current_element().coordinate_1 >= center_pos - lamb:
            left_edge = q2_list.get_current_element()  # type: DataInQ2
            if left_edge.coordinate_1 > center_pos + lamb:
                return []
            break
        else:
            q2_list.move_foreword()
        counter += 1

    # right edge
    if q2_list.get_current_element().coordinate_1 > center_pos + lamb:
        return []
    counter = 0  # type: int
    while 1:
        if counter == len_data_list:
            return []
        if q2_list.get_element(-1 - counter).coordinate_1 <= center_pos + lamb:
            if is_list_ended or counter > 0:
                right_edge = q2_list.get_element(-1 - counter)
                break
            else:
                if current_pos - center_pos > lamb:
                    right_edge = q2_list.get_last_element()  # type: DataInQ2
                    break
                else:
                    return None
        else:
            counter += 1

    # double check edges
    if int(right_edge.coordinate_1) < int(left_edge.coordinate_1):
        return []
    return q2_list.filter(lambda x: left_edge.coordinate_1 <= x.coordinate_1 <= right_edge.coordinate_1)


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
        return "99"
    else:
        return str(round(-10 * np.log10(value)))


def get_pl(rr_u, ra_u, rm_u, mm_u):
    pl_raw = -10 * np.array([rr_u, ra_u, rm_u, mm_u])
    pl = pl_raw - pl_raw.min()
    result = '{:.0f},{:.0f},{:.0f},{:.0f}'.format(*pl)
    return result


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
    :type data_in_q3_curr: PileupDataStruct
    :param data_in_q3_curr:
    :param q2_list: data queue to be windowed
    :type q2_list: BigForewordList
    :param lamb: half of window width
    :type lamb: int
    :return:
    True data is windowed, send result to W
    False can not window yet
    """

    # TODO <NB> Check length of q2_list before and after call!
    # if rel_pileup.coordinate_1 == 68044 and current_pos == 83302:
    #     import pdb; pdb.set_trace()
    tracked_data = window_data_one_chromosome(rel_pileup.coordinate_1,
        q2_list, my_args.lamb, q2_list_is_end, current_pos, worker_id)  # type: List[DataInQ2]
    
    # print([q2_list.len(), q2_list.pos, len(q2_list.my_list)])
    # print(rel_pileup.coordinate_1, q2_list.is_empty(), q2_list_is_end,
                                              # current_pos, worker_id)


    if tracked_data is None:
        return False

    bias = bias_estimator(rel_pileup.coordinate_1, tracked_data, my_args.lamb,
        my_args.bias)

    # <INFO, NB> Calculation of basic values for GQ and PL
    # rr = Sequencing  Noise (rr_u = log10(rr))
    # ra = Amplification Artefact/Error (ra_u = log10(ra))
    # rm = Heterozygous SNV (rm_u = log10(rm))
    # mm = Homozygous SNV (mm_u = log10(mm))
    [rr_u, ra_u, rm_u, mm_u] = sc_caller(rel_pileup, bias, my_args.null)  # 10 11 12 13

    rr = 10 ** rr_u  # type: float
    ra = 10 ** ra_u  # type: float
    rm = 10 ** rm_u  # type: float
    mm = 10 ** mm_u  # type: float

    if my_args.format == VARCALLFORMAT:
        vcf_info = None
        worker_type = WORKVAR
    else:
        worker_type = WORKVCF
        gq = get_gq(np.array([np.log10(rr + ra), rm_u, mm_u]))  # type: str

        vcf_info = VcfInfo(
            "",
            "{},{}".format(rel_pileup.reference_allele_num,
                rel_pileup.variant_allele_num),  # AD
            np.round(bias, 3),  # BI
            gq,  # GQ
            get_pl(rr_u, ra_u, rm_u, mm_u),  # type:str  # PL
            gq,  # QUAL
            rel_pileup.genotype_num,
            rel_pileup.genotype_class
        )

    outline = OutLineStruct(
        rel_pileup.name,  # name (str)
        rel_pileup.coordinate_1,  # pos (int)
        rel_pileup.reference_base,  # ref (str)
        rel_pileup.variant,  # var (str)
        rel_pileup.reference_allele_num,  # ref_num (int)
        rel_pileup.variant_allele_num, # var_num (int)
        bias,  # bias (float)
        rr,  # sn (float)
        ra,  # ae (float)
        rm,  # he (float)
        mm,  # ho (float)
        vcf_info, # vcf_info
        rel_pileup.so, # sp
        rel_pileup.variant_all #var_all
    )

    q5.put(DataInQ5(outline, worker_id, worker_type, []), block=False)
    return True


def read_bulk_mpileup(pileup, rm_minvar, rm_minmapQ, rm_mindepth):
    if "*" in pileup[4]:
        return "indel"
    if ("-" in pileup[4]) or ("+" in pileup[4]):
        if re.findall(INDEL, pileup[4]):
            return "indel"

    # remove head(^) and tail($)
    pileup[4] = remove_head_end(pileup)
    # mapQ
    map_q_filter = [ord(i) - PHREDSCORE >= rm_minmapQ for i in pileup[6]]
    read_base = compress_read_base(pileup[4], map_q_filter)
    maxvar = [read_base.count('A'), read_base.count('C'), read_base.count('G'),
        read_base.count('T')]

    if len(read_base) < rm_mindepth:
        return "lessmindepth"
    elif max(maxvar) < rm_minvar:
        return "refgenotype"  # return reference genome type and for skipping this line
    else:
        return "varreads"


def get_so_source(bulk_pileup_source, bulk_minvar, bulk_min_mapq, bulk_mindepth):
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
            bulk_pileup_list = bulk_pileup_source.next()
            if bulk_pileup_list == -1:  # pysam crashed
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
            bulk_pileup_list[4] = bulk_pileup_list[4].upper()
            should_read = True
            pos = yield read_bulk_mpileup(bulk_pileup_list, bulk_minvar,
                bulk_min_mapq, bulk_mindepth)


# @send_name_time(should_send=should_analyze, q7=q7)
def bias_estimator(pos, tracked_data, lamb, default):
    """
    calculate bias
    :type lamb: int
    :rtype: float
    :type tracked_data: list[DataInQ2]
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

        be_tmp = bias_kernel(int(pos), int(i.coordinate_1), lamb)  # K in the formula
        be_kwy.append(be_tmp * be_tmp1 * be_tmp2)
        be_kw.append(be_tmp * be_tmp1)
    # Nadaraya-Watson kernel-weighted average
    if len(be_kwy) > 0 and sum(be_kw) > 0:
        return sum(be_kwy) / sum(be_kw)
    # return args.bias when no neighboring heterozygous base
    return default


def bias_kernel(bk_x0, bk_xi, lamb):
    # Eq. 2
    if -lamb < bk_x0 - bk_xi < lamb:
        return 0.75 * (1 - (float(bk_x0 - bk_xi) / lamb) ** 2)
    else:
        return 0.0


def sc_caller(sc_candidate, sc_bias, sc_artifact):
    """
    Calculate RR, RA, RM, MM from data in q3 and bias.
    :type sc_bias: float
    :type sc_candidate: PileupDataStruct
    :param sc_candidate: data in q3
    :param sc_bias: Estimated bias in this position
    :param sc_artifact:
    :return:
    [RR, RA, RM, MM]
    """
    try:
        if "," in sc_candidate.variant:
            ref, mut = sc_candidate.variant.split(",")
        else:
            ref = sc_candidate.reference_base
            mut = sc_candidate.variant
    except:
        import pdb; pdb.set_trace()

    f_h0 = 0.125 * sc_bias
    f_h1 = sc_bias
    f_h2 = 1
    if sc_candidate.reference_allele_num > sc_candidate.variant_allele_num:
        f_h1 = 1 - f_h1
        f_h2 = 0

    rr_u = sum(map(lambda x: P_b_GG(x, ref, mut, sc_artifact), sc_candidate.read_info_list))
    ra_u = sum(map(lambda x: P_b_GG(x, ref, mut, f_h0), sc_candidate.read_info_list))
    rm_u = sum(map(lambda x: P_b_GG(x, ref, mut, f_h1), sc_candidate.read_info_list))
    mm_u = sum(map(lambda x: P_b_GG(x, ref, mut, f_h2), sc_candidate.read_info_list))

    return [rr_u, ra_u, rm_u, mm_u]


def P_b_GG(bp, ref, mut, f):
    """
    Formula (7), calculate lg Probability value
    :type bp: ReadInfo (read_base and base_quality)
    :type ref: str
    :type mut: str
    :type f: float
    :param bp: data of 1 read
    :param ref: reference base
    :param mut: mismatch
    :param f: {0, 0.125, bias (or 1-bias depending on mut>ref or ref>mut), 1} for {ref/ref, ref/artifacts, ref/mut; mut/mut}
    :return: lg Probability value
    """
    e = 10 ** (-bp.base_quality / 10.0)

    if bp.read_base in [',', '.', ref]:
        a = f * e / 3 + (1 - f) * (1 - e)
    elif bp.read_base == mut:
        a = (1 - f) * e / 3 + f * (1 - e)
    else:
        a = e / 3

    return np.log10(a)


def get_my_filename(output, suffix, prefix=None):
    base = os.path.basename(output)
    if not prefix:
        prefix =  os.path.dirname(output)
    prefix = os.path.dirname(output)
    return '{}{}{}'.format(prefix, os.path.splitext(base)[0], suffix)


def calculate_eta(list_var_buf, list_var_tag):
    """
    :param list_var_buf: list[list[OutLineStruct]]
    :return:
    """
    allele_nums = []
    bias_list = []
    for i, var_tag in enumerate(list_var_tag):
        if not var_tag:
            break
        for q5_item in list_var_buf[i]:
            allele_nums.append(q5_item.outlineStruct.var_num \
                + q5_item.outlineStruct.ref_num)
            bias_list.append(q5_item.outlineStruct.bias) 
    
    allele_nums = allele_nums[:CUTOFFNUM]
    bias_list = bias_list[:CUTOFFNUM]

    LLR = np.zeros(len(allele_nums))
    for i, allele_num in enumerate(allele_nums):
        f_artifact = 0.125 * bias_list[i] / 0.5
        alt = np.random.binomial(allele_num, f_artifact) #bin (depth(number of trials),prob_success)
        ref = allele_num - alt
        L_filter = (1 - f_artifact) ** ref * f_artifact ** alt
        ## random select major allele
        major = np.random.binomial(1, .5)
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
            msg_q5 = q5.get(block=True)  # type: DataInQ5
            w_id = msg_q5.worker_id

            # from main
            if w_id == 0:
                # try to write coverage file
                merged_coverage_list = merge_coverage(list_coverage_buf)
                if merged_coverage_list:
                    with open(cov_file, "a") as fp_cov:
                        fp_cov.write("\n".join(
                            map(lambda z: "\t".join(z),
                                map(lambda x: map(lambda y: str(y), x),
                                    merged_coverage_list)
                            )
                        ))
                if eta == -1:
                    if get_current_cutoff_num(list_var_buf, list_var_tag) > 0:
                        eta = calculate_eta(list_var_buf, list_var_tag)
                        for i, var_buf in enumerate(list_var_buf):
                            if var_buf:
                                write_do([j.outlineStruct for j in var_buf],
                                    msg_q5.work_type, i + 1, fp_out, eta,
                                    args.bulk, fp_reasoning, args.minvar)
                            else:
                                logging.info("worker{} write len = 0" \
                                    .format(i + 1))
                    else:
                        logging.info("Nothing need to be written.")
                fp_eta.write("##contig=<ID={},eta={}>\n".format(name, eta))
                break

            if msg_q5.work_type == WORKVAR or msg_q5.work_type == WORKVCF:

                # record tag and data
                if msg_q5.is_done():
                    list_var_tag[msg_q5.worker_id - 1] = True
                    logging.info("chr{} worker{} done".format(name, w_id))
                    # Try to calculate eta, write data
                    if get_current_cutoff_num(list_var_buf, list_var_tag) >= CUTOFFNUM \
                            and eta == -1:
                        eta = calculate_eta(list_var_buf, list_var_tag)
                    if eta != -1:
                        for i, var_buf in enumerate(list_var_buf):
                            if i < current_handling_worker - 1:
                                continue
                            if list_var_tag[i]:
                                if var_buf[i]:
                                    write_do([j.outlineStruct for j in var_buf],
                                        msg_q5.work_type, i + 1, fp_out, eta,
                                        args.bulk, fp_reasoning, args.minvar)
                                else:
                                    logging.info("Chr{} worker{} write len=0" \
                                        .format(name, i + 1))
                                current_handling_worker += 1
                            else:
                                break
                elif msg_q5.outlineStruct == "pysam crashed":
                    logging.info("worker{} pysam crashed. Data cleaned" \
                        .format(w_id))
                    list_var_tag[w_id - 1] = False
                    list_var_buf[w_id - 1] = []
                    list_coverage_buf[w_id - 1] = []
                else:
                    list_var_buf[w_id - 1].append(msg_q5)  # record data

            else:
                list_coverage_buf[w_id - 1].append(msg_q5.coverage)

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
    result_list = []  # [coverage1, coverage2...]  coverage: [name str ,start int ,stop int]
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


def write_do(data_list, work_type, worker_id, fp_out, eta, bulk, fp_reasoning,
            min_var):
    """

    :type data_list: list[OutLineStruct]
    """
    if work_type == WORKVAR:
        res_filter = [i.he != 0 and i.ae / i.he < eta \
                and 7 * i.var_num > i.ref_num and 2 * i.sn < i.he \
            for i in data_list]
        data_list = [j for i, j in enumerate(data_list) if res_filter[i]]

        logging.debug("worker{} write len={}".format(worker_id, len(data_list)))
        if len(data_list) > 0:
            fp_out.write("\n".join([pack_varcall_str(i) for i in data_list]))
            fp_out.write("\n")
        if bulk != "" and data_list:
            fp_reasoning.write("\n".join([
                "{}\t{}\t{}\t{}\t{}\t{}\t{}" \
                    .format(i.name, i.pos, i.ref, i.var, i.ref_num, i.var_num, i.so) \
                for i in data_list])
            )
            fp_reasoning.write("\n")
    elif work_type == WORKVCF:
        logging.info("worker{} write len = {}" .format(worker_id, len(data_list)))
        fp_out.write("\n".join([pack_vcf_str(i, eta, min_var) for i in data_list]))
        fp_out.write("\n")


def pack_vcf_str(outline, eta, min_var):
    """
    Assemble data of vcf according to the outline
    :type outline: OutLineStruct
    """
    outline.vcf_info.gt = get_gt(outline, eta)
    if outline.so != "":
        format_str = "GT:SO:AD:BI:GQ:FPL"
        cell_str = "{}:{}:{}:{}:{}:{}" \
            .format(outline.vcf_info.gt, transform_vcf_so(outline.so),
                outline.vcf_info.ad, outline.vcf_info.bi,
                outline.vcf_info.gq, outline.vcf_info.pl)
    else:
        format_str = "GT:AD:BI:GQ:FPL"
        cell_str = "{}:{}:{}:{}:{}" \
            .format(outline.vcf_info.gt, outline.vcf_info.ad,
                outline.vcf_info.bi, outline.vcf_info.gq,
                outline.vcf_info.pl)

    if outline.vcf_info.genotype_class == 1:
        filter_str = get_filter(outline.var, outline.vcf_info.genotype_num,
            min_var, outline.var_num)
    else:
        filter_str = get_filter(outline.var, outline.vcf_info.genotype_num,
            min_var, outline.ref_num)
    
    result = "\t".join([outline.name, str(outline.pos), ".", outline.ref,
        outline.var_all, outline.vcf_info.qual, filter_str, "NS=1", 
        format_str, cell_str])
    return result


def transform_vcf_so(so_in):
            if so_in == "refgenotype":
                return "True"
            elif so_in == "varreads":
                return "False"
            else:
                return "NA"

def get_gt(outline, eta):
    AE = outline.ae / eta
    SN = 2 * outline.sn
    HE = outline.he
    HO = outline.ho
    likelihood_list = [AE, SN, HE, HO]

    max_likelihood = max(likelihood_list)
    if max_likelihood in [AE, SN]:
        result_str = "0/0"
    elif max_likelihood == HE:
        # <TODO, NB> Why: 7 * var > ref ???
        result_str = "0/1"
        # if 7 * outline.var_num > outline.ref_num:
        #     result_str = "0/1"
        # else:
        #     result_str = "0/0"
    else:
        result_str = "1/1"

    if outline.vcf_info.genotype_class == 2:
        if result_str == "0/0":
            result_str = "1/1"
        elif result_str == "0/1":
            result_str = "1/2"
        else:
            result_str = "1/1"

    return result_str


def get_filter(var_str, genotype_num, min_var, var_num):
    comment = ""
    if len(var_str.split(",")) < genotype_num:
        comment += MULTIPLEGENOTYPE
    if var_num < min_var:
        if len(comment) == 0:
            comment = "{}{}".format(NOTENOUGHVARIANTS, min_var)
        else:
            comment += ",{}{}".format(NOTENOUGHVARIANTS, min_var)
    if len(comment) == 0:
        comment = "."
    return comment


def pack_varcall_str(outline):
    # type: (OutLineStruct) -> str
    return "{0}\t{1}\t{1}\t{2}\t{3}\t{4}\t{5}\tPASS\t{6}\t{7}\t{8}\t{9}\t{10}" \
        .format(outline.name, outline.pos, outline.ref, outline.var,
            outline.ref_num, outline.var_num, outline.bias, outline.sn,
            outline.ae, outline.he, outline.ho)


def data_generator(my_args, name, start, stop, is_bulk):
    if my_args.engine == ENGINESAMTOOLS:
        generator = data_generator_samtools
    else:
        generator = data_generator_pysam
    return generator(my_args, name, start, stop, is_bulk)


# type: (str, str, str, int, int, str) -> GeneratorExit
def data_generator_pysam(my_args, name, start, stop, is_bulk):
    """

    :param stop:
    :param name:
    :param start:
    :type fasta_file: str
    :type bam_file: str
    """

    fasta_file = pysam.FastaFile(my_args.fasta)
    str_ref = fasta_file.fetch(name, start - 1, stop + 1)

    my_arg = {"fastafile": fasta_file, "stepper": "samtools",
        "adjust_capq_threshold": 50, "contig": name, "start": start,
        "stop": stop, "min_mapping_quality": 0 if is_bulk else 40}

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
        cmd_str = "samtools mpileup -C50  -f {0} -s {1} -r {2}:{3}-{4}" \
            .format(my_args.fasta, my_args.bam, name, start, stop)
    else:
        cmd_str = "samtools mpileup -C50  -f {0} -q 40 -s {1} -r {2}:{3}-{4}" \
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
    # type: (Popen, str, str, int, list) -> int
    """
    Read data from the process's stdout, filter out the comments, split by spliter, and check the format
    :param process_handle: handle of the process
    :param ch: comment character
    :param spliter:
    :param columns: input. Greater than 0: column number of the data. Others, Unlimited
    :param data_out: output data (should clear buffer before using)
    :return: True, Success. Others, processes have been terminated, and stdout can't read the data.
    """
    while 1:
        buf = []
        if not read_line_from_file(process_handle.stdout, ch, spliter, columns, buf):
            if not process_handle.poll() is None:  # samtools has been terminated
                return False
            else:
                continue
        else:
            buf[4] = buf[4].upper()
            data_out.extend(buf)
            return True


def read_line_from_file(from_file, ch, spliter, columns, data_out):
    # type: (file, str, str, int, list) -> int
    """
    Read data from file, filter out comments, and split by spliter, check format
    :param from_file: file pointer
    :param ch: comment character
    :param spliter:
    :param columns: input. Greater than 0: column number of the data. Others, Unlimited
    :param data_out: output data (should clear buffer before using)
    :return: True, Success. Others, end of file.
    eg:
    data = []
    with open("test") as file:
    ret = read_data_from_file(file,"#","\t",8,data)
    """
    while 1:
        buf = from_file.readline()
        if len(buf) == 0:
            return False
        if buf[0] == ch:
            continue
        buf = buf.strip("\n")
        buf2 = buf.split(spliter)
        if columns > 0:
            if len(buf2) != columns:
                continue
            else:
                break
        else:
            break
    buf2[1] = int(buf2[1])
    data_out.extend(buf2)
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
    if worker_id != 1:
        my_start = head - my_args.lamb
    else:
        my_start = head

    if worker_id != my_args.work_num:
        my_stop = stop + my_args.lamb
    else:
        my_stop = stop
    total_work = my_stop - my_start

    # data source
    pileup_source = data_generator(my_args, name, my_start, my_stop, False)
    # bulk source
    if my_args.bulk != '':
        bulk_pileup_source = data_generator(my_args, name, head, stop, True)
        so_source = get_so_source(bulk_pileup_source, my_args.bulk_min_var,
            my_args.bulk_min_mapq, my_args.bulk_min_depth)
        x = so_source.next()

    row = (worker_id - 1) % (V_ROWS - 2) + 2
    col = (worker_id - 1) / (V_ROWS - 2) * 21 + 1
    my_print("\x1b[{};{}Hw{:<3d}\t0/{:<6}".format(row, col, worker_id, total_work))

    coverage_list = []  # list of [name str, head int, stop int]
    counter = 0
    main_index = 0
    while 1:
        pileup = pileup_source.next()

        if pileup == -1:  # pysam crashed
            logging.info("Use samtools engine instead!")
            q5.put(DataInQ5("pysam crashed", worker_id, WORKVAR, []), block=False)
            my_args.engine = ENGINESAMTOOLS
            control(my_args, copy_vcf_list, q5, name, head, stop, worker_id)
            return
        elif pileup is None:
            if my_args.coverage and len(coverage_list) > 0:
                q5_item = DataInQ5(None, worker_id, None, None, True)
                q5.put(q5_item, block=False)
                del coverage_list[0]
            break

        # append the relevant pileups into rel_pileups
        if head <= pileup[1] < stop:
            # Check if pileup is relevant/passes filters
            if my_args.format == VARCALLFORMAT:
                rel_pileup = read_mpileup(copy.copy(pileup),
                    my_args.minvar, my_args.mapq, my_args.min_depth, False,
                    worker_id) 
            else:
                rel_pileup = read_mpileup_vcf(copy.copy(pileup), my_args)

            if rel_pileup and rel_pileup != -1:
                if my_args.bulk != "":
                    #  calculate so
                    result = so_source.send(rel_pileup.coordinate_1)
                    if result == -1:
                        logging.info("Use samtools engine instead!")
                        q5.put(DataInQ5("pysam crashed", worker_id, WORKVAR, []),
                            block=False)
                        my_args.engine = ENGINESAMTOOLS
                        control(my_args, copy_vcf_list, q5, name, head, stop, worker_id)
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
                        q5.put(DataInQ5("", worker_id, WORKCOVERAGE, new_coverage),
                            block=False)

        # handle the head data in rel_pileups
        if not rel_pileups.is_empty():
            diff_success = differential(my_args, q2_list, 
                rel_pileups.get_current_element(), q5, worker_id, q2_list_is_end,
                pileup[1])
            if diff_success:
                rel_pileups.move_foreword()

        counter += 1
        if counter == 10000:
            counter = 0
            main_index += 1
            my_print("\x1b[{};{}Hw{:<3d} {:>6}/{:<6}"
                .format(row, col, worker_id, main_index * 10000, total_work),
                "control")
        # vcf screen
        if not q2_list_is_end:
            while not my_vcf_list.is_empty():
                if my_vcf_list.get_current_element() < pileup[1]:
                    my_vcf_list.move_foreword()
                else:
                    break
            if pileup[1] == my_vcf_list.get_current_element():
                ret = golden_hetero(pileup, 20, my_args.min_depth, my_args.snp_type,
                                    my_args.RD, worker_id)  # type: DataInQ2
                if ret is not None:
                    q2_list.append(ret)
            if my_vcf_list.is_empty() or pileup[1] >= my_vcf_list.get_last_element():
                q2_list_is_end = True

    while not rel_pileups.is_empty():
        diff_success = differential(my_args, q2_list,
            rel_pileups.get_current_element(), q5, worker_id, True, my_stop)
        if diff_success:
            rel_pileups.move_foreword()

    worker_type = WORKVAR if my_args.format == VARCALLFORMAT else WORKVCF
    q5_item = DataInQ5(None, worker_id, worker_type, None, done=True)
    q5.put(q5_item, block=False)

    my_print("\x1b[{};{}H\x1b[1;32mw{:<3d} I'm done!\x1b[0m" \
        .format(row, col, worker_id))
    logging.info("worker {} Quit!".format(worker_id))


def load_vcf(name, vcf_info, vcf_file_name):
    """
    Read the position information of the specified name chromosome from the vcf file, and store it as a list
    :param name: chromosome name
    :param vcf_info:
    :param vcf_file_name: vcf file name
    :return: POS information of vcf (list[int])
    """
    vcf_list = []
    curr_info = filter(lambda x: x[0] == name, vcf_info)
    if curr_info:
        with open(vcf_file_name) as fp:
            for i in curr_info:
                fp.seek(i[1])
                vcf_list.extend(
                    [int(i.split('\t')[1]) for i in fp.read(i[2]).splitlines()])
    return vcf_list


def parse_vcf(vcf_file, need_check):
    """
    Parse the vcf file into [[name,head,length],[name,head,length]...] format
    Here head, length is the file pointer
    :param vcf_file:
    :return: the vcf infomation [[name,head,length],[name,head,length]...]
    """
    has_shown_info = False
    base_file, file_type = os.path.splitext(vcf_file)

    catalog_file = "{}.catalog".format(base_file)
    if os.path.exists(catalog_file) \
            and os.path.getmtime(catalog_file) > os.path.getmtime(vcf_file):
        # read vcf catalog
        with open(catalog_file, "r") as fp:
            catalog = fp.read().split(";")
        result = map(
            lambda x: map(lambda y: y if x.split(",").index(y) == 0 else int(y),
                x.split(",")),
            catalog
        )
    else:
        if file_type == '.gz':
            raise IOError('VCF file {} needs to be unzipped'.format(vcf_file))
        # parse vcf
        result = []
        name = ""
        head = 0
        with open(vcf_file, "rb") as fp:
            while True:
                line = fp.readline()
                if line == "":
                    result.append([name, head, fp.tell() - head])
                    break
                if line[0] == "#" or line == "\n":
                    continue
                
                cols = line.split("\t")
                if name != cols[0]:
                    curr_idx = fp.tell() - len(line)
                    if name == "":
                        head = curr_idx
                    else:
                        tail = curr_idx
                        result.append([name, head, tail - head])
                        head = tail
                    name = cols[0]

                if need_check and not has_shown_info and len(cols) > 8:
                    tmp_str = "\t".join(cols[9:])
                    tmp_list = re.findall("0\\|0|0/0|1\\|1|1/1|2/2|2\\|2", tmp_str)
                    if len(tmp_list) > 0:
                        logging.info(">>>Please confirm the input VCF only " \
                            "contains heterozygote loci in bulk!!<<<")
                        has_shown_info = True

    return result


def main(my_args):
    # Init logging
    log_file = get_my_filename(my_args.output, "_{:0>2d}to{:0>2d}.log" \
        .format(my_args.head, my_args.tail), "sc_")
    logging.basicConfig(filename=log_file, level=logging.DEBUG,
        format=LOG_FORMAT, filemode="w")

    logging.info("Welcome to SCcaller v{}".format(VERSION))

    if my_args.debug:
        my_args.cpu_num = 1
        my_args.work_num = 1

    my_print("parsing vcf...")
    logging.info("parsing SNP vcf...")
    vcf_info = parse_vcf(my_args.snp_in, my_args.snp_type == "hsnp")

    my_print("parsing fasta...")
    logging.info("parsing reference fasta...")
    fasta_info = parse_fasta(my_args.fasta)

    if my_args.tail == -1 or my_args.tail > len(fasta_info):
        my_args.tail = len(fasta_info)

    logging.info("args: {}".format(my_args))

    my_print("\x1b[?25l")

    q5 = mp.Manager().Queue()
    # Start the operation process
    for j, fasta_j in enumerate(fasta_info):
        if j + 1 < my_args.head:
            continue
        if j + 1 > my_args.tail:
            break

        name, start, stop = fasta_j
        
        my_print("loading vcf...")
        logging.info("loading vcf...")
        list_vcf = load_vcf(name, vcf_info, my_args.snp_in)
        
        my_print("\x1b[2J\x1b[0;0HSCcaller v2.0.0 is handling chromosome {}...\x1b[0J" \
            .format(name))
        logging.debug("name={}; list_vcf len={}".format(name, len(list_vcf)))

        # Calculate the amount of tasks for each process
        no_bases = stop - start
        step = int(np.ceil(no_bases / float(my_args.work_num)))
        if step < my_args.lamb:
            logging.info("work_num={} is too large for {}. Use 1 instead. lamb={}" \
                .format(my_args.work_num, name, my_args.lamb))
            my_args.work_num = 1
            step = no_bases

        # Start the write file process
        proc_write_result = mp \
            .Process(target=write_result, args=(q5, my_args, name))
        proc_write_result.daemon = True
        proc_write_result.start()

        if not my_args.debug:
            process_pool = mp.Pool(processes=my_args.cpu_num)

        for woker_id in xrange(my_args.work_num):
            head = start + woker_id * step
            tail = head + step
            if woker_id != 0:
                if (head - my_args.lamb) < 0:
                    err_str = "lamb={} is too large. Fasta: name={} from {} to {}." \
                        .format(my_args.lamb, name, start, stop)
                    logging.critical(err_str)
                    raise RuntimeError(err_str)

            vcf_worker = [i for i in list_vcf \
                if head - my_args.lamb <= i <= tail + my_args.lamb]

            if my_args.debug:
                control(my_args, vcf_worker, q5, name, head, tail, woker_id + 1)
            else:
                process_pool.apply_async(control, 
                    (my_args, vcf_worker, q5, name, head, tail, woker_id + 1))

        if not my_args.debug:
            process_pool.close()
            process_pool.join()

        # Exit the W process
        worker_type = WORKVAR if my_args.format == VARCALLFORMAT else WORKVCF
        q5.put(DataInQ5(None, 0, worker_type, None), block=True)
        proc_write_result.join()

    write_vcf(my_args)

    my_print("\x1b[{};0H\x1b[?25h".format(my_args.work_num + 2))
    logging.info("W quit. All done.")


def write_vcf(args):
    """
    write the vcf file head
    :param my_format:
    :param fasta_file:
    :return:
    """

    head_str = "##fileformat=VCFv4.1\n" \
        "##fileDate={date}\n" \
        "##source=SCcaller_v{v}\n" \
        "##reference=file:{ref}" \
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
        "homozygous SNV respectively\">\n".format(
            date=time.strftime("%Y%m%d", time.localtime()), ref=args.fasta,
            mg=MULTIPLEGENOTYPE, v=VERSION, nev=NOTENOUGHVARIANTS, mv=args.minvar)

    if args.bulk != "":
        head_str += "##FORMAT=<ID=SO,Number=1,Type=String," \
            "Description=\"Whether it is a somatic mutation.\">\n"

    with open('{}.eta'.format(args.output), 'r') as f:
        etas = f.read()
    head_str += etas

    sc_name = os.path.basename(args.bam).split('.')[0]
    head_str += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n" \
        .format(sc_name)

    with open('{}.body'.format(args.output), "r") as f:
        body_str = f.read()

    with open(args.output, "w") as fp_out:
        if args.format == VCFFORMAT:
            fp_out.write(head_str)
        fp_out.write(body_str)


if __name__ == "__main__":
    start = time.time()
    args = parse_args()
    main(args)
    print("Duration = {}s".format(str(time.time() - start)))
