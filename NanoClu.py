import argparse
import os
import sys

import cluster_by_splice_site as cbss
import filter_fastq as ff
import run_poa as rp

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    raw_filter_group = parser.add_argument_group('Filter raw reads')
    raw_filter_group.add_argument('-f', '--file', type=str, help='input file, fasta/fastq format')
    raw_filter_group.add_argument('-l', '--len', type=int, default=200, help='minimum read length [200]')
    raw_filter_group.add_argument('-q', '--quality_cutoff', type=int, default=9, help='minimum average phred-score, only if input file is fastq format [9]')

    gmap_group = parser.add_argument_group('Align with GMAP')
    gmap_group.add_argument('-b', '--gmap_bin_dir', type=str, help='bin folder of GMAP')
    gmap_group.add_argument('-D', '--gmap_genome_dir', type=str, help='GMAP genome directory')
    gmap_group.add_argument('-d', '--gmap_db', type=str, help='GMAP reference genome database')
    gmap_group.add_argument('-t', '--gmap_thread', default=1, type=int, help='number of threads to run GMAP [1]')

    poa_group = parser.add_argument_group('Generating consensus with POA')
    poa_group.add_argument('-p', '--poa_bin_dir', type=str, help='bin folder of POA')
    poa_group.add_argument('-c', '--poa_cpus', type=int, default=1, help='number of processors used to run POA at one time')

    args = parser.parse_args()

    if args.file is None or args.gmap_db is None:
        parser.print_help()
        sys.exit(1)


    # 1. filter fastq based on length and phred-score
    in_fq = args.file
    quality = args.quality_cutoff
    length = args.len
    filter_fq = os.path.dirname(in_fq)+'/filter.fq'
    ff.filter_fastq(in_fq, filter_fq, length, quality)

    # 2. align with GMAP
    if args.gmap_bin_dir is None:
        gmap_bin = 'gmap'
    else:
        gmap_bin = args.gmap_bin_dir+'/gmap'

    if args.gmap_genome_dir is None:
        gmap_D = ''
    else:
        gmap_D = '-D ' + args.gmap_genome_dir
    gmap_db = args.gmap_db
    gmap_thread = args.gmap_thread
    gmap_sam = filter_fq+'.gamp.sam'
    gmap_log = filter_fq+'.gmap.log'
    min_intron_len = 25
    sys.stdout.write('%s %s -d %s -t %d --min-intronlength %d -f samse -O %s > %s 2> %s\n'
              %(gmap_bin, gmap_D, gmap_db, gmap_thread, min_intron_len, filter_fq, gmap_sam, gmap_log))
    os.system('%s %s -d %s -t %d --min-intronlength %d -f samse -O %s > %s 2> %s'
              %(gmap_bin, gmap_D, gmap_db, gmap_thread, min_intron_len, filter_fq, gmap_sam, gmap_log))

    # 3. samtools sort GMAP alignment
    samtools = 'samtools'
    gmap_bam = filter_fq+'.gamp.bam'
    gmap_sort_bam = filter_fq+'.gamp.sort.bam'
    sys.stdout.write('%s view -bS %s > %s\n' % (samtools, gmap_sam, gmap_bam))
    os.system('%s view -bS %s > %s' % (samtools, gmap_sam, gmap_bam))
    sys.stdout.write('%s sort -@ %d %s > %s\n' % (samtools, gmap_thread, gmap_bam, gmap_sort_bam))
    os.system('%s sort -@ %d %s > %s' % (samtools, gmap_thread, gmap_bam, gmap_sort_bam))

    # 4. sort bam record with splice site
    full_isoform_list = ''
    first_last_list = ''
    first_or_last_list = ''
    cbss.clus_by_splice_site(gmap_sort_bam, '')

    # 5. run POA for each cluster
    if args.poa_bin_dir is Nono:
        poa='poa'
    else:
        poa=args.poa_bin_dir+'/poa'
    poa_cpus = args.poa_cpus
    rp.run_poa(full_isoform_list, poa_cpus, '')
