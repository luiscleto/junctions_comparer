from __future__ import print_function
import sys
import csv
import os
from itertools import groupby


def e_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def represents_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def split_csv(filename, outputdir):
    chromosomes_found = set()
    for key, rows in groupby(csv.reader(open(filename), delimiter='\t'), lambda r: r[0]):
        nk = str(key).replace("chr", "")
        if len(nk) == 1 and represents_int(nk):
            nk = "0" + nk
        chromosomes_found.add(nk)
        with open(os.path.join(outputdir, os.path.splitext(os.path.basename(filename))[0] + "_%s.bed" % nk), "w+") as output:
            for row in rows:
                output.write("\t".join(row) + "\n")
    return chromosomes_found
