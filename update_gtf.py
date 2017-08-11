import utils


def update_gtf(fp, cons_sam, rRNA, old_gtf, out_gtf):
    cons_filter_bam = cons_sam + '.filter.bam'
    gtools = 'gtools'
    samtools = 'samtools'
    out_gtf_stat = out_gtf + '.stat'
    utils.exec_cmd(fp, 'filter_cons_bam', '%s filter %s %s | %s sort > %s' % (gtools, cons_sam, rRNA,
                                                                          samtools, cons_filter_bam))
    utils.exec_cmd(fp, 'update_gtf', '%s update-gtf %s %s > %s 2> %s' % (gtools, cons_filter_bam,
                                                                     old_gtf, out_gtf, out_gtf_stat))

    return
