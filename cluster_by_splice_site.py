import collections as cl
import math
import os
import re
import sys
import threading as td

import gffutils as gu
import numpy as np
import pysam as ps

import utils


def find_all_occ(string, character):
    index = []
    for x in re.finditer(character, string):
        index.append(x.start())
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
            op_len += c[1]
    return op_len


def get_op_length(cigartuples):
    # get op length for MIDNSHP=XB
    # input: cigartuples is a list of tuples
    # return: Counter({op:op_len})
    op_len = cl.Counter()
    op_dict = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9}
    for op in op_dict:
        op_all = [i for i in cigartuples if i[0] == op_dict[op]]
        op_len[op] = sum(i[1] for i in op_all)
    return op_len


# return: 1-base, exonic base position
def cigar_to_splice_site(ref_start, cigar, min_intron_len):  # ref_start: 1-base
    N_idx = [(i, item[1]) for i, item in enumerate(cigar) if item[0] == 3 and item[1] >= min_intron_len]  # 'N':3
    if len(N_idx) > 0:
        five_site = []
        three_site = []
        cur_pos = ref_start
        last_N = -1
        for i, item in N_idx:
            ref_len = get_ref_op_length(cigar[last_N + 1:i])
            last_N = i

            five_site.append(cur_pos + ref_len - 1)
            three_site.append(cur_pos + ref_len + item)
            cur_pos += (ref_len + item)
        return five_site, three_site
    else:
        return [], []


# maxmize: square(sum(site cnt)) / site bin cnt
def splice_site_bin_DP(site, bin_size, bin_dis):
    if not site:
        return [], []

    solu_avg = cl.defaultdict(tuple)
    solu_bin = cl.defaultdict(list)

    for i, site_i in enumerate(site):
        # solu_avg: object, bin number, total site number
        solu_avg[(i, i)] = (1.0, 1, 1)
        solu_avg[(i, i - 1)] = (0, 0, 0)
        # solu_bin: (start, end):[bin tuples(bin1_start, bin2_end), (bin2) ... ]
        solu_bin[(i, i)] = (site_i, site_i)  # start,end  1-base

    for i, site_i in enumerate(site[1:], 1):
        # calculate solu[(j,i)]
        # every bin in (j, i-1) is pre-calculated
        for j in range(i - 1, -1, -1):
            # 0. exclude site i
            ex_solu = solu_avg[(j, i - 1)]
            ex_solu_bin = list(solu_bin[(j, i - 1)])

            # 1. include site i
            # 1.1 add site i to last site bin
            # last_bin = solu_last_bin[(j, i - 1)]
            # update last site bin, update all other bin
            new_last_bin_start = i
            new_other_bin_end = j - 1
            for k in range(0, i):
                if site_i - site[k] < bin_size:
                    new_last_bin_start = k
                    break
            for k in range(new_last_bin_start - 1, j - 1, -1):
                if site[new_last_bin_start] - site[k] > bin_dis:
                    new_other_bin_end = k
                    break
            new_last_avg = (math.pow(i - new_last_bin_start + 1.0, 2), 1, i - new_last_bin_start + 1.0)
            other_avg = solu_avg[(j, new_other_bin_end)]  # update_all_other_bin
            in_add_solu = (math.pow(other_avg[2] + new_last_avg[2], 2) / (other_avg[1] + new_last_avg[1]),
                           other_avg[1] + new_last_avg[1], other_avg[2] + new_last_avg[2])
            in_add_solu_bin = list(solu_bin[(j, new_other_bin_end)])
            in_add_solu_bin.append((site[new_last_bin_start], site_i))

            # 1.2 add site i to a new bin
            new_other_bin_end = j - 1
            for k in range(i - 1, j - 1, -1):
                if site_i - site[k] > bin_dis:
                    new_other_bin_end = k
                    break
            other_avg = solu_avg[(j, new_other_bin_end)]
            in_new_solu = (math.pow(other_avg[2] + 1, 2) / (other_avg[1] + 1),
                           other_avg[1] + 1, other_avg[2] + 1)
            in_new_solu_bin = list(solu_bin[(j, new_other_bin_end)])
            in_new_solu_bin.append((site_i, site_i))

            if ex_solu[0] >= in_add_solu[0] and ex_solu[0] >= in_new_solu[0]:
                solu_avg[(j, i)] = ex_solu
                solu_bin[(j, i)] = ex_solu_bin
            elif in_add_solu[0] >= ex_solu[0] and in_add_solu[0] >= in_new_solu[0]:
                solu_avg[(j, i)] = in_add_solu
                solu_bin[(j, i)] = in_add_solu_bin
            else:  # in_new_solu[0]
                solu_avg[(j, i)] = in_new_solu
                solu_bin[(j, i)] = in_new_solu_bin
    # print solu_avg[(0, len(site) - 1)], solu_bin[(0, len(site) - 1)]
    return solu_avg, solu_bin


def get_max_bin(site, bin_size, min_cnt):
    max_bin = (-1, -1)
    max_idx = (-1, -1)
    max_bin_std = float('inf')
    max_bin_cnt = 0
    for start_i, start_site in enumerate(site):
        bin_cnt = 1
        _end_i = start_i
        for end_i, end_site in enumerate(site[start_i + 1:], start_i + 1):
            if end_site - start_site < bin_size:
                bin_cnt += 1
                _end_i = end_i
            else:
                break
        bin_std = np.std(site[start_i:_end_i + 1])
        if bin_cnt > max_bin_cnt or bin_cnt == max_bin_cnt and bin_std < max_bin_std:
            max_bin_cnt = bin_cnt
            max_bin = (start_site, site[_end_i])
            max_idx = (start_i, _end_i)
            max_bin_std = bin_std

    if max_bin_cnt >= min_cnt:
        return max_bin, max_idx, max_bin_cnt
    else:
        return [], [], 0


def check_bin_dis(site_bin, max_bin, bin_dis):
    if not site_bin:
        return True
    for start, end in site_bin:
        if (end + bin_dis >= max_bin[0] > end) \
                or (start - bin_dis <= max_bin[1] < start):
            return False
    return True


# greedily find bin with most sites/second most sites/ ..., until reach min_cnt threshold
def splice_site_bin_greedy(site, bin_size, bin_dis, min_cnt):
    if not site:
        return []
    site_bin = []
    while site:
        max_bin, max_idx, max_bin_size = get_max_bin(site, bin_size, min_cnt)
        if not max_bin:
            break
        if check_bin_dis(site_bin, max_bin, bin_dis):
            site_bin.append(max_bin)

        site = site[0:max_idx[0]] + site[max_idx[1] + 1:len(site)]
    return site_bin


# TODO: TSS/TES
# XXX: ignore +/- strand
def gen_splice_site_from_gtf(gtf_db, rname, min_site, max_site):
    five_site = []
    three_site = []
    for exon in gtf_db.region((rname, min_site, max_site), featuretype='exon'):
        five_site.append(exon.stop)
        three_site.append(exon.start)
    return list(set(five_site)), list(set(three_site))


def gen_se_site_from_gtf(gtf_db, rname, min_site, max_site):
    start_site = []
    end_site = []
    for transcript in gtf_db.region((rname, min_site, max_site), featuretype='transcript'):
        start_site.append(transcript.start)
        end_site.append(transcript.stop)
    return list(set(start_site)), list(set(end_site))


# bin_size: with gtf: [bin_size--,anno_site,--bin_size]; without gtf: [bin_size]
# bin_dis:  [site bin]--bin_dis--[site bin]
# min_cnt:  size of [site bin] >= min_cnt
def infer_splice_site(gtf_db, rname, five_ss, three_ss, bin_size, bin_dis, min_cnt):
    five_site_bin = []  # (start, end) : count
    three_site_bin = []
    novel_five_ss = []
    novel_three_ss = []
    five_ss.sort()
    three_ss.sort()
    if gtf_db:
        # supervised clustering
        # 1. cluster site into nearest annotated site or novel site_bin
        # 2. can NOT cluster, add new novel site_bin
        # generate gtf site
        max_site = max(five_ss[-1], three_ss[-1])
        min_site = min(five_ss[0], three_ss[0])
        anno_five_site, anno_three_site = gen_splice_site_from_gtf(gtf_db, rname, min_site, max_site)
        for five_s in five_ss:
            anno_five_s = min(anno_five_site, key=lambda x: abs(x - five_s))
            if abs(anno_five_s - five_s) <= bin_size:
                five_site_bin.append((anno_five_s - bin_size, anno_five_s + bin_size))
            elif abs(anno_five_s - five_s) > bin_dis:
                novel_five_ss.append(five_s)
        for three_s in three_ss:
            anno_three_s = min(anno_three_site, key=lambda x: abs(x - three_s))
            if abs(anno_three_s - three_s) <= bin_size:
                three_site_bin.append((anno_three_s - bin_size, anno_three_s + bin_size))
            elif abs(anno_three_s - three_s) > bin_dis:
                novel_three_ss.append(three_s)
        five_site_bin.extend(splice_site_bin_greedy(novel_five_ss, bin_size, bin_dis, min_cnt))
        three_site_bin.extend(splice_site_bin_greedy(novel_three_ss, bin_size, bin_dis, min_cnt))
    else:
        # unsupervised clustering
        five_site_bin = splice_site_bin_greedy(five_ss, bin_size, bin_dis, min_cnt)
        three_site_bin = splice_site_bin_greedy(three_ss, bin_size, bin_dis, min_cnt)

    return five_site_bin, three_site_bin


def infer_start_end_site(gtf_db, rname, start_ss, end_ss, bin_size, bin_dis, min_cnt):
    start_site_bin = []
    end_site_bin = []
    novel_start_ss = []
    novel_end_ss = []
    start_ss.sort()
    end_ss.sort()

    if gtf_db:
        # supervised clustering
        # 1. cluster site into nearest annotated site or novel site_bin
        # 2. can NOT cluster, add new novel site_bin
        # generate gtf site
        max_site = max(start_ss[-1], end_ss[-1])
        min_site = min(start_ss[0], end_ss[0])
        anno_start_site, anno_end_site = gen_se_site_from_gtf(gtf_db, rname, min_site, max_site)
        for start_s in start_ss:
            anno_start_s = min(anno_start_site, key=lambda x: abs(x - start_s))
            if abs(anno_start_s - start_s) <= bin_size:
                start_site_bin.append((anno_start_s - bin_size, anno_start_s + bin_size))
            elif abs(anno_start_s - start_s) > bin_dis:
                novel_start_ss.append(start_s)
        for end_s in end_ss:
            anno_end_s = min(anno_end_site, key=lambda x: abs(x - end_s))
            if abs(anno_end_s - end_s) <= bin_size:
                end_site_bin.append((anno_end_s - bin_size, anno_end_s + bin_size))
            elif abs(anno_end_s - end_s) > bin_dis:
                novel_end_ss.append(end_s)
        start_site_bin.extend(splice_site_bin_greedy(novel_start_ss, bin_size, bin_dis, min_cnt))
        end_site_bin.extend(splice_site_bin_greedy(novel_end_ss, bin_size, bin_dis, min_cnt))
    else:
        # unsupervised clustering
        start_site_bin = splice_site_bin_greedy(start_ss, bin_size, bin_dis, min_cnt)
        end_site_bin = splice_site_bin_greedy(end_ss, bin_size, bin_dis, min_cnt)

    return start_site_bin, end_site_bin


def infer_site_bin(gtf_db, rname, start_ss, end_ss, five_ss, three_ss,
                   bin_size, bin_dis, min_bin_cnt):
    (start_site_bin, end_site_bin) = ([], [])
    (five_site_bin, three_site_bin) = ([], [])

    if len(start_ss) > 0:
        start_site_bin, end_site_bin = infer_start_end_site(gtf_db, rname, start_ss, end_ss, bin_size, bin_dis,
                                                            min_bin_cnt)

    if len(five_ss) > 0:
        five_site_bin, three_site_bin = infer_splice_site(gtf_db, rname, five_ss, three_ss, bin_size, bin_dis,
                                                          min_bin_cnt)
    return start_site_bin, end_site_bin, five_site_bin, three_site_bin


# site_seq: ((start, end), (five_site ...), (three_site ...)) : [seq]
def assign_splice_site_bin(site_seq, five_site_bin, three_site_bin):
    bin_seq_set = {}  # (first_bin, last_bin) : seq
    # assign site to bin
    for site, seq in site_seq.items():
        # site[0] : (start, end)
        # site[1] : (five_site1, five_site2 ...)
        # site[2] : (three_site1, three_site2 ...)

        site_bin = []
        for i in range(len(site[1])):
            five_bin = [five_bin for five_bin in five_site_bin if five_bin[1] >= site[1][i] >= five_bin[0]]
            three_bin = [three_bin for three_bin in three_site_bin if three_bin[1] >= site[2][i] >= three_bin[0]]
            if len(five_bin) == 0 or len(three_bin) == 0:
                continue
            # print five_bin[0], three_bin[0]
            site_bin.append(five_bin[0])
            site_bin.append(three_bin[0])

        site_bin_tuple = tuple(site_bin)
        if len(site_bin_tuple) == 0:
            continue

        if site_bin_tuple in bin_seq_set:
            bin_seq_set[site_bin_tuple].extend(seq)
        else:
            bin_seq_set[site_bin_tuple] = seq
    return bin_seq_set


# site_seq: ((start, end), (five_site ...), (three_site ...)) : [seq]
def assign_first_last_site_bin(site_seq, five_site_bin, three_site_bin):
    bin_seq_set = cl.defaultdict(list)  # (first_bin, last_bin) : seq
    # assign site to bin
    for site, seq in site_seq.items():
        # site[0] : (start, end)
        # site[1] : (five_site1, five_site2 ...)
        # site[2] : (three_site1, three_site2 ...)
        first_bin = [five_bin for five_bin in five_site_bin if five_bin[1] >= site[1][0] >= five_bin[0]]
        if len(first_bin) != 1:
            continue
        last_bin = [three_bin for three_bin in three_site_bin if three_bin[1] >= site[2][-1] >= three_bin[0]]
        if len(last_bin) != 1:
            continue

        first_bin = first_bin[0]
        last_bin = last_bin[0]

        if (first_bin, last_bin) in bin_seq_set:
            bin_seq_set[(first_bin, last_bin)].extend(seq)
        else:
            bin_seq_set[(first_bin, last_bin)] = seq
    return bin_seq_set


def assign_se_site_bin(site_seq, start_site_bin, end_site_bin):
    bin_seq_set = cl.defaultdict(list)  # (first_bin, last_bin) : seq
    # assign site to bin
    for site, seq in site_seq.items():
        # site[0] : (start, end)
        # site[1] : (five_site1, five_site2 ...)
        # site[2] : (three_site1, three_site2 ...)
        start_bin = [start_bin for start_bin in start_site_bin if start_bin[1] >= site[0][0] >= start_bin[0]]
        if len(start_bin) != 1:
            continue
        end_bin = [end_bin for end_bin in end_site_bin if end_bin[1] >= site[0][1] >= end_bin[0]]
        if len(end_bin) != 1:
            continue

        start_bin = start_bin[0]
        end_bin = end_bin[0]

        if (start_bin, end_bin) in bin_seq_set:
            bin_seq_set[(start_bin, end_bin)].extend(seq)
        else:
            bin_seq_set[(start_bin, end_bin)] = seq
    return bin_seq_set


def assign_se_sp_site_bin(site_seq, start_site_bin, end_site_bin, five_site_bin, three_site_bin):
    bin_seq_set = cl.defaultdict(list)  # (first_bin, last_bin) : seq
    # assign site to bin
    for site, seq in site_seq.items():
        # site[0] : (start, end)
        # site[1] : (five_site1, five_site2 ...)
        # site[2] : (three_site1, three_site2 ...)
        start_bin = [start_bin for start_bin in start_site_bin if start_bin[1] >= site[0][0] >= start_bin[0]]
        if len(start_bin) != 1:
            continue
        end_bin = [end_bin for end_bin in end_site_bin if end_bin[1] >= site[0][1] >= end_bin[0]]
        if len(end_bin) != 1:
            continue

        start_bin = start_bin[0]
        end_bin = end_bin[0]

        splice_site_bin = []
        for i in range(len(site[1])):
            five_bin = [five_bin for five_bin in five_site_bin if five_bin[1] >= site[1][i] >= five_bin[0]]
            three_bin = [three_bin for three_bin in three_site_bin if three_bin[1] >= site[2][i] >= three_bin[0]]
            if len(five_bin) == 0 or len(three_bin) == 0:
                continue
            # print five_bin[0], three_bin[0]
            splice_site_bin.append(five_bin[0])
            splice_site_bin.append(three_bin[0])
        if len(splice_site_bin) == 0:
            continue

        site_bin = [start_bin, end_bin]
        site_bin.extend(splice_site_bin)
        site_bin = tuple(site_bin)
        if len(site_bin) == 0:
            continue

        if site_bin in bin_seq_set:
            bin_seq_set[site_bin].extend(seq)
        else:
            bin_seq_set[site_bin] = seq
    return bin_seq_set


def write_clu_seq(bin_seq_set, min_clu_size, rname, base_dir):
    clu_n = 0
    clu_seq_n = 0
    for site_bin, seq in bin_seq_set.items():
        if len(seq) >= min_clu_size:
            fa_fn = base_dir + '/' + rname + '_' + \
                    str(np.mean(site_bin[0])) + '_' + str(np.mean(site_bin[1]))
            fa_fp = open(fa_fn, 'w')

            for i, s in enumerate(seq):
                fa_fp.write('>' + os.path.basename(fa_fn) + '_' + str(i) + '\n')
                fa_fp.write(s + '\n')
            fa_fp.close()
            clu_n += 1
            clu_seq_n += len(seq)
    return clu_n, clu_seq_n


def gen_clu_seq(clu_mode, min_clu_size, rname, base_dir, site_seq, start_site_bin, end_site_bin,
                five_site_bin, three_site_bin):
    # 1. assign each record's first/last splice site with five_site_bin/three_site_bin
    if clu_mode == 1:  # first and end
        # (frist_bin, last_bin) : [seq]
        bin_seq_set = assign_first_last_site_bin(site_seq, five_site_bin, three_site_bin)
        clu_n, clu_seq_n = write_clu_seq(bin_seq_set, min_clu_size, rname, base_dir)
    elif clu_mode == 2 or clu_mode == 3:
        # (start_bin, end_bin) : [seq]
        bin_seq_set = assign_se_site_bin(site_seq, start_site_bin, end_site_bin)
        clu_n, clu_seq_n = write_clu_seq(bin_seq_set, min_clu_size, rname, base_dir)
    elif clu_mode == 4:
        # ([],[]) : [seq]
        bin_seq_set = assign_splice_site_bin(site_seq, five_site_bin, three_site_bin)
        clu_n, clu_seq_n = write_clu_seq(bin_seq_set, min_clu_size, rname, base_dir)
    elif clu_mode == 5:
        # ((start_bin, end_bin), ([], [])) : [seq]
        bin_seq_set = assign_se_sp_site_bin(site_seq, start_site_bin, end_site_bin, five_site_bin,
                                            three_site_bin)
        clu_n, clu_seq_n = write_clu_seq(bin_seq_set, min_clu_size, rname, base_dir)
    else:
        sys.stderr.write('Unknown clu-mode: %d\n' % clu_mode)
        sys.exit(1)

    return clu_n, clu_seq_n


def infer_cluster(shared_site_seq, min_clu_size):
    cluster_site_seq = {}
    last_site_seq = []
    seq = []
    clu_n = 0
    clu_seq_n = 0
    for i, site_seq in enumerate(shared_site_seq):
        if len(last_site_seq) != 0 and len(set(site_seq[:-1]).intersection(set(last_site_seq[:-1]))) < 1:
            if len(seq) >= min_clu_size:
                cluster_site_seq[str(last_site_seq[0])] = seq
                clu_n += 1
            seq = []

        seq.append(site_seq[-1])
        clu_seq_n += 1
        last_site_seq = site_seq
    if len(seq) >= min_clu_size:
        cluster_site_seq[str(last_site_seq[0])] = seq
        clu_n += 1

    return cluster_site_seq, clu_n, clu_seq_n


def write_clu_seq_with_shared_site(rname, base_dir, shared_site_seq):
    for name, seq in shared_site_seq.items():
        fa_fn = base_dir + '/' + rname + '_' + name
        fa_fp = open(fa_fn, 'w')

        for i, s in enumerate(seq):
            fa_fp.write('>' + os.path.basename(fa_fn) + '_' + str(i) + '\n')
            fa_fp.write(s + '\n')
        fa_fp.close()
    return


# for each bam bundle:
# 1. infer splice site
# 2. assign each record's (first/last) splice site with a inferred splice-site bin(20 bp)
# 3. cluster by first 5' and last 3' splice site TODO: use different criteria to cluster read
# 4. TODO: keep or discard clipping sequence
def clu_by_splice_site_core(gtf_db, clu_mode, min_intron_len, bin_size, bin_dis, min_bin_cnt, min_clu_size, bam_bundle,
                            rname, base_dir):
    if not bam_bundle:
        return
    bam_bundle_read_cnt = 0
    # 0. generate splice site from bam record
    if clu_mode == 0:  # default mode
        five_site_set = cl.Counter()
        three_site_set = cl.Counter()
        site_seq = []
        shared_site_seq = []
        for record in bam_bundle:
            bam_bundle_read_cnt += 1
            query_seq = record.query_sequence[record.query_alignment_start:record.query_alignment_end]
            ref_start = record.reference_start + 1  # 1-bae exonic base
            # ref_end = record.reference_end + 1  # 1-base exonic base
            cigar = record.cigartuples

            (five_site, three_site) = cigar_to_splice_site(ref_start, cigar, min_intron_len)

            site = []
            for five, three in zip(five_site, three_site):
                five_site_set[five] += 1
                three_site_set[three] += 1
                site.append(five)
                site.append(three)
            site.append(query_seq)
            site_seq.append(tuple(site))

        # TODO: annotation splice site min cnt
        if min_bin_cnt < 1:
            min_bin_cnt = bam_bundle_read_cnt * min_bin_cnt
        if min_clu_size < 1:
            min_clu_size = bam_bundle_read_cnt * min_clu_size

        shared_five_site = [five for five in five_site_set if five_site_set[five] >= min_bin_cnt]
        shared_three_site = [three for three in three_site_set if three_site_set[three] >= min_bin_cnt]

        for ss in site_seq:
            shared_site = []
            for i, site in enumerate(ss[:-1]):
                if i % 2 == 0 and site in shared_five_site:
                    shared_site.append(site)
                elif i % 2 != 0 and site in shared_three_site:
                    shared_site.append(site)
            if len(shared_site) > 0:
                shared_site.append(ss[-1])
                shared_site_seq.append(tuple(shared_site))

        # [shared_site, shared_site, ..., seq]
        shared_site_seq.sort()
        # cluster into clusters
        cluster_site_seq, clu_n, clu_seq_n = infer_cluster(shared_site_seq, min_clu_size)
        # write cluster seq
        write_clu_seq_with_shared_site(rname, base_dir, cluster_site_seq)
    else:
        start_site_set = []
        end_site_set = []
        five_site_set = []
        three_site_set = []
        site_seq = {}
        for record in bam_bundle:
            bam_bundle_read_cnt += 1
            query_seq = record.query_sequence[record.query_alignment_start:record.query_alignment_end]
            ref_start = record.reference_start + 1  # 1-bae exonic base
            ref_end = record.reference_end + 1  # 1-base exonic base
            cigar = record.cigartuples

            # record => bin_seq = {((start, end), (five_site), (three_site)) : query_seq}
            (start, end) = (ref_start, ref_end)
            five_site = []
            three_site = []
            # XXX: use strand or XS of most clustered read to represent strand of clustered sequences
            # start/end site
            if clu_mode == 2 or clu_mode == 3 or clu_mode == 5:
                start_site_set.append(ref_start)
                end_site_set.append(ref_end)

            # internal splice site, five_site, three_site: 1-base exonic position
            if clu_mode == 0 or clu_mode == 1 or clu_mode == 4 or clu_mode == 5:
                (five_site, three_site) = cigar_to_splice_site(ref_start, cigar, min_intron_len)
                five_site_set.extend(five_site)
                three_site_set.extend(three_site)

            if ((start, end), tuple(five_site), tuple(three_site)) in site_seq:
                site_seq[((start, end), tuple(five_site), tuple(three_site))].extend(query_seq)
            else:
                site_seq[((start, end), tuple(five_site), tuple(three_site))] = [query_seq]

        if min_bin_cnt < 1:
            min_bin_cnt = bam_bundle_read_cnt * min_bin_cnt
        if min_clu_size < 1:
            min_clu_size = bam_bundle_read_cnt * min_clu_size
        # 1. generate TSS/TES/splice site bin with/without annotation
        start_site_bin, end_site_bin, five_site_bin, three_site_bin = \
            infer_site_bin(gtf_db, rname, start_site_set, end_site_set,
                           five_site_set, three_site_set, bin_size,
                           bin_dis, min_bin_cnt)

        # 2. generate clusters
        clu_n, clu_seq_n = gen_clu_seq(clu_mode, min_clu_size, rname, base_dir, site_seq,
                                       start_site_bin, end_site_bin,
                                       five_site_bin, three_site_bin)

    return clu_n, clu_seq_n


def clu_thread(id, bam_bundle_list_n, bam_bundle_list, bam_bundle_rname, bam_bundle_clu_out, gtf_db, clu_mode,
               min_intron_len, bin_size, bin_dis, min_bin_cnt, min_clu_size, base_dir):
    global CLU_THREAD_I
    global CLU_THREAD_LK

    while 1:
        CLU_THREAD_LK.acquire()
        list_i = CLU_THREAD_I
        CLU_THREAD_I += 1
        CLU_THREAD_LK.release()

        if bam_bundle_list_n <= list_i:
            break

        bam_bundle = bam_bundle_list[list_i]
        rname = bam_bundle_rname[list_i]

        n, seq_n = clu_by_splice_site_core(gtf_db, clu_mode, min_intron_len, bin_size, bin_dis, min_bin_cnt, min_clu_size,
                                           bam_bundle, rname, base_dir)
        bam_bundle_clu_out[list_i] = (n, seq_n)
    return


def clu_by_splice_site(logfp, gtf, clu_mode, min_intron_len, bin_size, bin_dis, min_bin_cnt, min_clu_size, in_bam,
                       clu_folder,
                       thread_n):
    global CLU_THREAD_I
    global CLU_THREAD_LK

    filter_bam_fp = ps.AlignmentFile(in_bam, 'rb')
    gtf_db = None
    if gtf is not None:
        gtf_db_fn = os.path.dirname(in_bam) + '/' + gtf + '.db'
        if not os.path.isfile(gtf_db_fn):
            try:
                gtf_db = gu.create_db(gtf, gtf_db_fn)
            except:
                logfp.write('Error in parsing %s\n' % gtf)
                sys.exit(IOError)
        else:
            try:
                gtf_db = gu.FeatureDB(gtf_db_fn)
            except:
                logfp.write('Error in parsing %s\n' % gtf_db_fn)
                sys.exit(IOError)

    utils.format_time(logfp, __name__, 'Clustering alignments ... \n')

    bam_iter = filter_bam_fp.fetch()

    bundle_thd = 0
    last_record = None
    bam_bundle_chunk = 1000
    bam_bundle = []

    bam_bundle_list = []
    bam_bundle_rname = []
    bam_bundle_list_n = 0
    threads = []
    CLU_THREAD_LK = td.Lock()
    clu_n = 0
    clu_seq_n = 0

    if not os.path.isdir(clu_folder):
        utils.exec_cmd(logfp, 'gen_clu_seq', 'mkdir %s 2> /dev/null' % clu_folder)
    else:
        utils.exec_cmd(logfp, 'gen_clu_seq', 'rm %s/* -rf 2> /dev/null' % clu_folder)

    for record in bam_iter:
        # separate bam bundle if cur_start - last_end > thd
        if last_record and (record.reference_id != last_record.reference_id
                            or record.reference_start - last_record.reference_end > bundle_thd):

            bam_bundle_list.append(bam_bundle)
            bam_bundle_rname.append(last_record.reference_name)
            bam_bundle_list_n += 1

            if bam_bundle_list_n == bam_bundle_chunk:
                if not os.path.isdir(clu_folder + '/' + str(clu_n / 1000)):
                    utils.exec_cmd(logfp, 'gen_clu_seq', 'mkdir %s/%d 2> /dev/null' % (clu_folder, clu_n / 1000))
                base_dir = clu_folder + '/' + str(clu_n / 1000)

                bam_bundle_clu_out = [(0, 0)] * bam_bundle_list_n
                CLU_THREAD_I = 0

                for i in range(thread_n):
                    t = td.Thread(target=clu_thread, args=(i, bam_bundle_list_n, bam_bundle_list, bam_bundle_rname,
                                                           bam_bundle_clu_out, gtf_db, clu_mode, min_intron_len, bin_size,
                                                           bin_dis, min_bin_cnt, min_clu_size, base_dir,))
                    threads.append(t)
                    t.start()
                for t in threads:
                    t.join()

                clu_n += sum(first for first, second in bam_bundle_clu_out)
                clu_seq_n += sum(second for first, second in bam_bundle_clu_out)

                bam_bundle_list = []
                bam_bundle_rname = []
                bam_bundle_list_n = 0

            bam_bundle = []

        # add to bundle
        bam_bundle.append(record)
        last_record = record

    bam_bundle_list.append(bam_bundle)
    bam_bundle_rname.append(last_record.reference_name)
    bam_bundle_list_n += 1
    if not os.path.isdir(clu_folder + '/' + str(clu_n / 1000)):
        utils.exec_cmd(logfp, 'gen_clu_seq', 'mkdir %s/%d 2> /dev/null' % (clu_folder, clu_n / 1000))
    base_dir = clu_folder + '/' + str(clu_n / 1000)

    bam_bundle_clu_out = [(0, 0)] * bam_bundle_list_n
    CLU_THREAD_I = 0

    for i in range(thread_n):
        t = td.Thread(target=clu_thread, args=(i, bam_bundle_list_n, bam_bundle_list, bam_bundle_rname,
                                               bam_bundle_clu_out, gtf_db, clu_mode, min_intron_len, bin_size, bin_dis,
                                               min_bin_cnt, min_clu_size, base_dir,))
        threads.append(t)
        t.start()
    for t in threads:
        t.join()

    clu_n += sum(first for first, second in bam_bundle_clu_out)
    clu_seq_n += sum(second for first, second in bam_bundle_clu_out)

    utils.format_time(logfp, __name__, 'Read clutsers: %d %d\n' % (clu_n, clu_seq_n))
    filter_bam_fp.close()
    return
