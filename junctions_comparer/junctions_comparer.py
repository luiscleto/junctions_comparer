#!/usr/bin/python

from utils import *


__temp_dir__ = "sj_tmp"

if len(sys.argv) < 3:
    e_print('Insufficient parameters provided!')
    print('Usage:')
    print('\tpython junctions_comparer.py <sample1.bed> <sample2.bed> [<sample3.bed> ...]')


def process_samples(file_list):
    for f in file_list:
        split_csv(f,__temp_dir__)

if not os.path.exists(__temp_dir__):
    os.makedirs(__temp_dir__)
process_samples(sys.argv[1:])
