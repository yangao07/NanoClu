import collections
import os
import re
import sys

import pysam as ps
import utils


def find_all_occ(string, character):
    index = [x.start() for x in re.finditer(character, string)]
    return index


def get_ref_op_length(cigartuples):
    # get op length for MDNP=XB
    # input: cigartuples is a list of tuples
    # return: ref_len
    op_len = 0
    for c in cigartuples:
        if c[0] == 0 or c[0] == 2 or c[0] == 3 or c[0] == 7 or c[0] == 8:
            op_len += c[1]
    return op_len


def get_read_op_length(cigartuples):
    # get op length for MDNP=XB
    # input: cigartuples is a list of tuples
    # return: read_len
    op_len = 0
    for c in cigartuples:
        if c[0] == 0 or c[0] == 1 or c[0] == 4 or c[0] == 7 or c[0] == 8:
            op_len += op[1]
    return op_len


def get_op_length(cigartuples):
    # get op length for MIDNSHP=XB
    # input: cigartuples is a list of tuples
    # return: Counter({op:op_len})
    op_len = collections.Counter()
    op_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9}
    for op in op_dict:
        op_all = [i for i in cigartuples if i[0] == op_dict[op]]
        op_len[op] = sum(i[1] for i in op_all)
    return op_len


def clus_by_splice_site(in_bam):
    # 1. identify first/last splice site for each bam record
    filter_bam_fp = ps.AlignmentFile(in_bam, 'rb')
    single_exon_set = []
    splice_site_set = []

    utils.format_time(__name__, 'Clustering alignments\n')

    clu_cnt = 0
    bam_iter = filter_bam_fp.fetch()
    for record in bam_iter:
        query_name = record.query_name
        query_seq = record.query_alignment_sequence
        ref_start = record.reference_start + 1  # 1-bae exonic base
        ref_end = record.reference_end + 1  # 1-base exonic base
        cigar = record.cigartuples

        N_idx = [i for i, item in enumerate(cigar) if item[0] == 3]  # 'N':3
        if len(N_idx) > 0:  # multi-exon read
            first_site = ref_start + get_ref_op_length(cigar[:N_idx[0]]) - 1  # 1-base exonic base
            last_site = ref_end - get_ref_op_length(cigar[N_idx[-1] + 1:])  # 1-base exonic base
            splice_site_set.append(
                (query_name, record.reference_name, record.is_reverse, first_site, last_site, query_seq))
        else:  # single-exon read
            single_exon_set.append(
                (query_name, record.reference_name, record.is_reverse, ref_start, ref_end, query_seq))

    utils.format_time(__name__, 'Read clutsers: %d\n' % clu_cnt)

    print 'multi-exons'
    for i in splice_site_set:
        print i
    print 'single-exon'
    for i in single_exon_set:
        print i

    filter_bam_fp.close()
    return
