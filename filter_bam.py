import os
import sys

import pysam


def filter_bam(samtools, in_bam, out_bam):
    in_fp = pysam.AlignmentFile(in_bam, 'rb')
    out_fp = pysam.AlignmentFile(out_bam, 'wb', template=in_fp)
    for record in in_fp.fetch():
        out_fp.write(record)

    sys.stdout.write('%s index %s' % (samtools, out_bam))
    in_fp.close()
    out_fp.close()

    os.system('%s index %s' % (samtools, out_bam))

