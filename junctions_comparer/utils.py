from __future__ import print_function
import sys
import csv
import os
from itertools import groupby


def e_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def split_csv(filename, outputdir):
    for key, rows in groupby(csv.reader(open(filename), delimiter='\t'),
                             lambda row: row[0]):
        with open(os.path.join(outputdir, os.path.splitext(os.path.basename(filename))[0] + "_%s.bed" % key), "w+") as output:
            for row in rows:
                output.write("\t".join(row) + "\n")
