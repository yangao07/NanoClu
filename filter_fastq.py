#import sys, os, argparse
import numpy as np
import pysam as ps


def filter_fastq(infile, outfile, length=200, quality=9):
    with ps.FastxFile(infile) as fin, open(outfile, mode='w') as fout:
        for entry in fin:
            if len(entry.sequence) >= length:
                if not entry.quality == None:
                    nq = []
                    for q in entry.quality:
                        nq.append(ord(q) - 33)
                    # print np.average(nq)
                    if np.average(nq) >= quality:
                        fout.write(str(entry)+'\n')
                else:
                    fout.write(str(entry)+'\n')
    fin.close()
    fout.close()
