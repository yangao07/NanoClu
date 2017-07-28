import argparse
import os
import sys

import numpy as np
import pysam as ps
import utils


def filter_fastq(infile, outfile, length, quality):
    try:
        fin = ps.FastxFile(infile)
    except:
        sys.stderr.write('Error in parsing read file %s\n' % (infile))
        sys.exit(IOError)
    try:
        fout = open(outfile, 'w')
    except:
        sys.stderr.write('Error in opening output read file %s\n' % (outfile))
        sys.exit(IOError)

    utils.format_time(__name__, 'Filtering reads\n')

    cnt = 0
    for entry in fin:
        if len(entry.sequence) >= length:
            if entry.quality is not None:
                nq = []
                for q in entry.quality:
                    nq.append(ord(q) - 33)
                if np.average(nq) >= quality:
                    fout.write(str(entry) + '\n')
                    cnt += 1
            else:
                fout.write(str(entry) + '\n')
                cnt += 1
    utils.format_time(__name__, 'Filtered reads: %s\n' % cnt)

    fin.close()
    fout.close()
