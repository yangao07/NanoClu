#!/usr/bin/python
import argparse
import os
import sys

import cluster_by_splice_site as cbss
import filter_bam as fb
import filter_fastq as ff
import run_poa as rp


def parse_argv(args):
    # 1. filter fastq
    argv = {}
    argv['in_fq'] = args.file
    argv['quality'] = args.quality_cutoff
    argv['length'] = args.len
    argv['filter_fq'] = args.file + '_filter.fq'
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
    argv['gmap_thread'] = args.gmap_thread
    argv['gmap_sam'] = argv['filter_fq'] + '.gmap.sam'
    argv['gmap_log'] = argv['filter_fq'] + '.gmap.log'

    # 3. filter bam
    argv['gmap_bam'] = argv['filter_fq'] + '.gmap.bam'
    argv['gmap_sort_bam'] = argv['filter_fq'] + '.gmap.sort.bam'
    argv['gmap_filter_bam'] = argv['filter_fq'] + '.gmap.filter.bam'
    # 4 generate consensus with poa
    if args.poa_bin is None:
        argv['poa_bin'] = 'poa'
    else:
        argv['poa_bin'] = args.poa_bin
    argv['samtools'] = args.samtools
    argv['poa_cpus'] = args.poa_cpus
    return argv


def filter_read(argv):
    ff.filter_fastq(argv['in_fq'], argv['filter_fq'], argv['length'], argv['quality'])


def align_with_gmap(argv):
    min_intron_len = 25
    sys.stdout.write('%s %s -d %s -t %d --min-intronlength %d -f samse -O %s > %s 2> %s\n'
                     % (argv['gmap_bin'], argv['gmap_D'], argv['gmap_db'], argv['gmap_thread'],
                        min_intron_len, argv['filter_fq'], argv['gmap_sam'], argv['gmap_log']))
    os.system('%s %s -d %s -t %d --min-intronlength %d -f samse -O %s > %s 2> %s\n'
              % (argv['gmap_bin'], argv['gmap_D'], argv['gmap_db'], argv['gmap_thread'],
                 min_intron_len, argv['filter_fq'], argv['gmap_sam'], argv['gmap_log']))


def filter_bam(argv):
    samtools = argv['samtools']
    gmap_sam = argv['gmap_sam']
    gmap_bam = argv['gmap_bam']
    gmap_sort_bam = argv['gmap_sort_bam']
    gmap_filter_bam = argv['gmap_filter_bam']
    gmap_thread = argv['gmap_thread']

    sys.stdout.write('%s view -bS %s > %s\n' % (samtools, gmap_sam, gmap_bam))
    os.system('%s view -bS %s > %s' % (samtools, gmap_sam, gmap_bam))
    sys.stdout.write('%s sort -@ %d %s > %s\n' % (samtools, gmap_thread, gmap_bam, gmap_sort_bam))
    os.system('%s sort -@ %d %s > %s' % (samtools, gmap_thread, gmap_bam, gmap_sort_bam))
    sys.stdout.write('%s index %s\n' % (samtools, gmap_sort_bam))
    os.system('%s index %s' % (samtools, gmap_sort_bam))

    fb.filter_bam(samtools, gmap_sort_bam, gmap_filter_bam)


def gen_cluster(argv):
    gmap_filter_bam = argv['gmap_filter_bam']
    cbss.clus_by_splice_site(gmap_filter_bam)


def gen_cons_with_poa(argv):
    poa = argv['poa_bin']
    poa_cpus = argv['poa_cpus']
    rp.run_poa(poa, poa_cpus)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    raw_filter_group = parser.add_argument_group('Filter raw reads')
    raw_filter_group.add_argument('-f', '--file', type=str, help='input file, fasta/fastq format')
    raw_filter_group.add_argument('-l', '--len', type=int, default=200,
                                  help='minimum read length [200]')
    raw_filter_group.add_argument('-q', '--quality_cutoff', type=int, default=9,
                                  help='minimum average phred-score, only if input file is fastq format [9]')

    gmap_group = parser.add_argument_group('Align with GMAP')
    gmap_group.add_argument('-b', '--gmap_bin', type=str, default='gmap', help='path of GMAP [gmap]')
    gmap_group.add_argument('-D', '--gmap_genome_dir', type=str, help='GMAP genome directory')
    gmap_group.add_argument('-d', '--gmap_db', type=str, help='GMAP reference genome database')
    gmap_group.add_argument('-t', '--gmap_thread', type=int, default=1,
                            help='number of threads to run GMAP [1]')

    poa_group = parser.add_argument_group('Generating consensus with POA')
    poa_group.add_argument('-p', '--poa_bin', type=str, default='poa', help='path of POA [poa]')
    poa_group.add_argument('-s', '--samtools', type=str, default='samtools', help='path of samtools [samtools]')
    poa_group.add_argument('-c', '--poa_cpus', type=int, default=1,
                           help='number of processors used to run POA at one time')

    args = parser.parse_args()

    if args.file is None or args.gmap_db is None:
        parser.print_help()
        sys.exit(1)

    # 0. parse argument
    argv = parse_argv(args)

    # 1. filter fastq based on length and phred-score
    filter_read(argv)

    # 2. align with GMAP
    align_with_gmap(argv)

    # 3. filter GMAP alignment
    filter_bam(argv) #TODO: filter with 1.best/secondary, 2.aligned ratio, 3. identical bases ratio

    # 4. sort bam record with splice site
    gen_cluster(argv) #TODO: merge splice site into bins

    # 5. generate consensus with poa
    gen_cons_with_poa(argv)
