from __future__ import print_function
import sys
import csv
import os
from itertools import groupby


def e_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def file_len(file_name):
    i = 0
    with open(file_name) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def represents_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def split_csv(filename, output_dir):
    chromosomes_found = set()
    for key, rows in groupby(csv.reader(open(filename), delimiter='\t'), lambda r: r[0]):
        nk = str(key).replace("chr", "")
        if len(nk) == 1 and represents_int(nk):
            nk = "0" + nk
        chromosomes_found.add(nk)
        with open(os.path.join(output_dir, os.path.splitext(os.path.basename(filename))[0] + "_%s.bed" % nk), "w+") as output:
            for row in rows:
                output.write("\t".join(row) + "\n")
    return chromosomes_found


# Print iterations progress
def print_progress (iteration, total, prefix = '', suffix = '', decimals = 2, bar_length = 40):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : number of decimals in percent complete (Int)
        barLength   - Optional  : character length of bar (Int)
    """
    filled_length    = int(round(bar_length * iteration / float(total)))
    percents        = round(100.00 * (iteration / float(total)), decimals)
    bar             = '#' * filled_length + '-' * (bar_length - filled_length)
    sys.stdout.write('\r%s [%s] %s%s %s' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("")
