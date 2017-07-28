import os
import sys

import gffutils as gu
import pysam as ps
import utils


def filter_align_rate(record, align_rate):
    read_len = record.infer_query_length(always=True) + 0.0
    query_len = record.query_alignment_length
    return query_len / read_len < align_rate


def filter_iden_rate(record, iden_rate):
    query_len = record.query_alignment_length + 0.0
    del_all = [i for i in record.cigartuples if i[0] == 2]  # 'D'
    del_len = sum(i[1] for i in del_all)
    match_len = query_len - record.get_tag('NM') + del_len
    score = match_len
    return match_len / query_len < iden_rate, score


def filter_rRNA(record, rRNA_db):
    rname = record.reference_name
    ref_start = record.reference_start + 1  # 1-base
    ref_end = record.reference_end + 1  # 1-base
    for i in rRNA_db.features_of_type('transcript'):
        if i.chrom is rname and not (ref_start > i.stop or ref_end < i.start):
            return True
    return False


def filter_bam(samtools, in_sam, filter_sort_bam, thread, align_rate, iden_rate, best_rate, rRNA):
    in_fp = ps.AlignmentFile(in_sam, 'r')
    filtered_bam = in_sam + '.filter'
    filter_fp = ps.AlignmentFile(filtered_bam, 'wb', template=in_fp)
    if rRNA is not None:
        rRNA_db_fn = os.path.dirname(in_sam) + '/' + rRNA + '.db'
        if not os.path.isfile(rRNA_db_fn):
            try:
                rRNA_db = gu.create_db(rRNA, rRNA_db_fn)
            except:
                sys.stderr.write('Error in parsing %s\n' % (rRNA))
                sys.exit(IOError)
        else:
            try:
                rRNA_db = gffutils.FeatureDB(rRNA_db_fn)
            except:
                sys.stderr.write('Error in parsing %s\n' % (rRNA_db_fn))
                sys.exit(IOError)

    last_qname = ''
    b_score = 0
    s_score = 0
    cnt = 0
    utils.format_time(__name__, 'Filtering alignments\n')
    for record in in_fp.fetch():
        if record.is_unmapped:
            continue
        # 1. aligned bases
        if filter_align_rate(record, align_rate):
            continue
        # 2. identical rate
        filter_iden, score = filter_iden_rate(record, iden_rate)
        if filter_iden:
            continue
        # 3. rRNA
        if rRNA and filter_rRNA(record, rRNA_db):
            continue
        # 4. best
        if record.query_name is last_qname:
            if score > b_score:
                b_record = record
                s_score = b_score
                b_score = score
            elif score > s_score:
                s_score = score
        else:
            if last_qname is not '' and s_score < best_rate * b_score:
                filter_fp.write(b_record)
                cnt += 1
            b_record = record
            b_score = score
            s_score = 0
            last_qname = record.query_name

    if last_qname is not '' and s_score < best_rate * b_score:
        filter_fp.write(b_record)
        cnt += 1

    in_fp.close()
    filter_fp.close()

    utils.format_time(__name__, 'Filtered alignments: ' + str(cnt) + '\n')

    utils.exec_cmd(__name__, '%s sort -@ %d %s > %s' % (samtools, thread, filtered_bam, filter_sort_bam))
    utils.exec_cmd(__name__, '%s index %s' % (samtools, filter_sort_bam))
    utils.exec_cmd(__name__, 'rm %s' % filtered_bam)
