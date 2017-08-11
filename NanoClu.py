#!/usr/bin/python
import argparse
import os
import sys

import cluster_by_splice_site as cbss
import filter_bam as fb
import filter_fastq as ff
import run_poa as rp
import update_gtf as ug
import utils


def parse_argv(args):
    # 1. filter fastq
    argv = {}
    argv['in_fq'] = args.file
    argv['output'] = args.output
    argv['quality'] = args.quality_cutoff
    argv['length'] = args.len
    argv['filter_fq'] = args.output + '/' + os.path.basename(args.file) + '_filter.fq'
    # 2. align with GMAP
    if args.gmap_bin is None:
        argv['gmap_bin'] = 'gmap'
    else:
        argv['gmap_bin'] = args.gmap_bin
    if args.gmap_genome_dir is None:
        argv['gmap_D'] = ''
    else:
        argv['gmap_D'] = '-D ' + args.gmap_genome_dir
    argv['gmap_db'] = args.gmap_db
    argv['min_intron'] = args.min_intron
    argv['gmap_thread'] = args.gmap_thread
    argv['gmap_sam'] = argv['filter_fq'] + '.gmap.sam'
    argv['gmap_log'] = argv['filter_fq'] + '.gmap.log'

    # 3. filter bam
    argv['align_rate'] = args.align_rate
    argv['match_rate'] = args.match_rate
    argv['best_rate'] = args.best_rate
    argv['rRNA'] = args.rRNA
    argv['gmap_filter_sort_bam'] = argv['filter_fq'] + '.gmap.filter.sort.bam'
    # 4 generate consensus with poa
    if args.poa_bin is None:
        argv['poa_bin'] = 'poa'
    else:
        argv['poa_bin'] = args.poa_bin
    argv['samtools'] = args.samtools
    argv['gtf'] = args.gtf
    argv['clu_mode'] = args.clu_mode
    argv['splice_site_bin_size'] = args.bin_size
    argv['splice_site_bin_dis'] = args.bin_dis
    argv['splice_site_bin_min_cnt'] = float(args.min_bin_cnt)
    argv['min_clu_size'] = float(args.min_clu_size)
    argv['clu_folder'] = argv['output'] + '/clu_seq' + '_M' + str(argv['clu_mode']) + '_b' + str(
        argv['splice_site_bin_min_cnt']) + '_c' + str(argv['min_clu_size']) + '/'
    # align consensus seq with GAMP
    argv['cons'] = argv['output'] + '/' + 'M' + str(argv['clu_mode']) + '_b' + str(
        argv['splice_site_bin_min_cnt']) + '_c' + str(argv['min_clu_size']) + '_cons.fa'
    argv['cons_sam'] = argv['cons'] + '.gmap.sam'
    argv['cons_log'] = argv['cons'] + '.gmap.log'
    # update gtf
    argv['out_gtf'] = args.update_gtf
    # log
    argv['log'] = argv['output'] + '/' + 'M' + str(argv['clu_mode']) + '_b' + str(
        argv['splice_site_bin_min_cnt']) + '_c' + str(argv['min_clu_size']) + '.log'
    argv['logfp'] = open(argv['log'], 'w')
    return argv


def filter_read(argv):
    logfp = argv['logfp']
    ff.filter_fastq(logfp, argv['in_fq'], argv['filter_fq'], argv['length'], argv['quality'])


def align_filtered_reads_with_gmap(argv):
    min_intron_len = 25
    logfp = argv['logfp']
    utils.exec_cmd(logfp, 'align_with_gamp', '%s %s -d %s -t %d --min-intronlength %d -f samse -O %s > %s 2> %s' %
                   (argv['gmap_bin'], argv['gmap_D'], argv['gmap_db'], argv['gmap_thread'],
                    min_intron_len, argv['filter_fq'], argv['gmap_sam'], argv['gmap_log']))


def filter_bam(argv):
    samtools = argv['samtools']
    gmap_sam = argv['gmap_sam']
    align_rate = argv['align_rate']
    match_rate = argv['match_rate']
    best_rate = argv['best_rate']
    rRNA = argv['rRNA']
    gmap_filter_sort_bam = argv['gmap_filter_sort_bam']
    gmap_thread = argv['gmap_thread']
    intron_N = 1  # filter out non-spliced reads
    min_intron_len = 25
    logfp = argv['logfp']

    fb.filter_bam(logfp, samtools, gmap_sam, gmap_filter_sort_bam, gmap_thread, align_rate, match_rate,
                  best_rate, min_intron_len, intron_N, rRNA)


def gen_cluster(argv):
    gtf = argv['gtf']
    min_intron = argv['min_intron']
    gmap_filter_sort_bam = argv['gmap_filter_sort_bam']
    clu_mode = argv['clu_mode']
    bin_site_dis = argv['splice_site_bin_dis']
    bin_site_size = argv['splice_site_bin_size']
    min_bin_cnt = argv['splice_site_bin_min_cnt']
    min_clu_size = argv['min_clu_size']
    clu_folder = argv['clu_folder']
    thread_n = argv['gmap_thread']
    logfp = argv['logfp']

    cbss.clu_by_splice_site(logfp, gtf, clu_mode, min_intron, bin_site_size, bin_site_dis, min_bin_cnt, min_clu_size,
                            gmap_filter_sort_bam, clu_folder, thread_n)


def gen_cons_with_poa(argv):
    poa = argv['poa_bin']
    thread_n = argv['gmap_thread']
    clu_folder = argv['clu_folder']
    cons = argv['cons']
    logfp = argv['logfp']

    rp.run_poa(logfp, poa, thread_n, clu_folder, cons)


def align_cons_with_gmap(argv):
    min_intron_len = 25
    cons = argv['cons']
    cons_sam = argv['cons_sam']
    cons_log = argv['cons_log']
    logfp = argv['logfp']

    utils.exec_cmd(logfp, 'align_with_gamp', '%s %s -d %s -t %d --min-intronlength %d -f samse -O %s > %s 2> %s' %
                   (argv['gmap_bin'], argv['gmap_D'], argv['gmap_db'], argv['gmap_thread'],
                    min_intron_len, cons, cons_sam, cons_log))


def update_gtf(argv):
    old_gtf = argv['gtf']
    cons_sam = argv['cons_sam']
    rRNA = argv['rRNA']
    out_gtf = argv['out_gtf']
    logfp = argv['logfp']

    ug.update_gtf(logfp, cons_sam, rRNA, old_gtf, out_gtf)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    raw_filter_group = parser.add_argument_group('Filter raw reads')
    raw_filter_group.add_argument('-f', '--file', type=str, required=True, help='input file, fasta/fastq format')
    raw_filter_group.add_argument('-o', '--output', type=str, required=True, help='output folder')
    raw_filter_group.add_argument('-l', '--len', type=int, default=200,
                                  help='minimum read length [200]')
    raw_filter_group.add_argument('-q', '--quality-cutoff', type=int, default=9,
                                  help='minimum average phred-score, only if input file is fastq format [9]')

    bam_filter_group = parser.add_argument_group('Filter bam alignemt')
    bam_filter_group.add_argument('-a', '--align-rate', default=0.67, type=float,
                                  help='minimum fraction of aligned bases [0.80]')
    argument = bam_filter_group.add_argument('-m', '--match-rate', default=0.75, type=float,
                                             help='minimum fraction of identically matched bases [0.85]')
    bam_filter_group.add_argument('-b', '--best-rate', default=0.98, type=float,
                                  help='maximum secondary/best alignment score to retain record [0.80]')
    add_argument = bam_filter_group.add_argument('-r', '--rRNA', metavar='rRNA.gtf', type=str,
                                                 help='GTF annotation of rRNA to filter out')

    gmap_group = parser.add_argument_group('Align with GMAP')
    gmap_group.add_argument('-g', '--gmap-bin', type=str, default='gmap',
                            help='path of GMAP [gmap]')
    gmap_group.add_argument('-D', '--gmap-genome-dir', type=str, help='GMAP genome directory')
    gmap_group.add_argument('-d', '--gmap-db', type=str, required=True, help='GMAP reference genome database')
    gmap_group.add_argument('-i', '--min-intron', type=int, default=25, help='minimum intron length [25]')
    gmap_group.add_argument('-t', '--gmap-thread', type=int, default=1,
                            help='number of threads to run GMAP [1]')

    poa_group = parser.add_argument_group('Generating consensus with POA')
    poa_group.add_argument('-p', '--poa-bin', type=str, default='poa', help='path of POA [poa]')
    poa_group.add_argument('--samtools', type=str, default='samtools', help='path of samtools [samtools]')
    poa_group.add_argument('-G', '--gtf', type=str, help='use annotation file to guide generating consensus sequence')
    poa_group.add_argument('-M', '--clu-mode', type=int, default=0, help='specify strategy to cluster reads [0]. ' +
                                                                         '(0): first and last splice site; ' +
                                                                         '(1): first or last splice site; ' +
                                                                         '(2): start and end site; ' +
                                                                         '(3): start or end site; ' +
                                                                         '(4): all internal splice site; ' +
                                                                         '(5): start, end site and all internal splice site')
    poa_group.add_argument('--bin-size', type=int, default=20,
                           help='inferred splice site bin size [20]')
    poa_group.add_argument('--bin-dis', type=int, default=30,
                           help='minimum distance of two adjacent splice site bins [30]')
    poa_group.add_argument('--min-bin-cnt', type=str, default='2',
                           help='minimum number or fraction of splice sites to retain the site bin [2]')
    poa_group.add_argument('--min-clu-size', type=str, default='0.02',
                           help='minimum number or fraction of cluster reads to retain the read cluster [0.02]')

    gtf_group = parser.add_argument_group('Update gtf with consensus alignment')
    gtf_group.add_argument('-u', '--update-gtf', type=str, required=True, help='updated gtf file')

    args = parser.parse_args()

    if args.file is None or args.gmap_db is None:
        parser.print_help()
        sys.exit(1)

    # 0. parse argument
    argv = parse_argv(args)

    # 1. filter fastq based on length and phred-score
    filter_read(argv)

    # 2. pesudo-align with nanoclu
    # pseudo_align_with_nanoclu(argv)
    # 2. align filtered reads with GMAP
    align_filtered_reads_with_gmap(argv)

    # 3. filter GMAP alignment
    filter_bam(argv)  # filter with 1.best/secondary, 2.aligned ratio, 3. identical bases ratio

    # 4. sort bam record with splice site
    gen_cluster(argv)  # TODO: merge splice site into bins

    # 5. generate consensus with poa
    gen_cons_with_poa(argv)

    # 6. align consensus sequence with GMAP
    align_cons_with_gmap(argv)

    # 7. generate updated GTF with cons alignment
    # update_gtf(argv)
