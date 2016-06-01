#!/usr/bin/python

from utils import *

if len(sys.argv) < 3:
    e_print('Insufficient parameters provided!')
    print('Usage:')
    print('\tpython junctions_comparer.py <sample1.bed> <sample2.bed> [<sample3.bed> ...]')

