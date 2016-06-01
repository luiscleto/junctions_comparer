#!/usr/bin/python

from utils import *


__temp_dir__ = "sj_tmp"
__out_dir__ = "sj_out"
__out_file__ = "results.csv"
__out_delimiter__ = ","

if len(sys.argv) < 3:
    e_print('[ERROR] Insufficient parameters provided!')
    print('Usage:')
    print('\tpython junctions_comparer.py <sample1.bed> <sample2.bed> [<sample3.bed> ...]')


def write_chromosome_junctions_to_file(chromosome, junctions):
    writer = csv.writer(open(os.path.join(__out_dir__, __out_file__), 'a'), delimiter=__out_delimiter__, lineterminator='\n')
    for key, value in junctions.items():
        row = [key]
        row.extend(value)
        writer.writerow(row)


def read_junctions(samples, chromosomes):
    for c in chromosomes:
        chromosome_junctions = {}
        sample_index = 0
        for s in map(lambda sam: os.path.splitext(os.path.basename(sam))[0], samples):
            file_name = os.path.join(__temp_dir__, s + "_" + c + ".bed")
            with open(file_name, 'rb') as csv_file:
                csv_reader = csv.reader(csv_file, delimiter='\t')
                for row in csv_reader:
                    junc_id = row[0] + "_" + row[1] + "_" + row[2]
                    if junc_id not in chromosome_junctions.keys():
                        chromosome_junctions[junc_id] = [0] * len(samples)
                    chromosome_junctions[junc_id][sample_index] += 1
            sample_index += 1
        write_chromosome_junctions_to_file(c, chromosome_junctions)


def clean_files(samples, chromosomes):
    for c in chromosomes:
        for s in map(lambda sam: os.path.splitext(os.path.basename(sam))[0], samples):
            file_name = os.path.join(__temp_dir__, s + "_" + c + ".bed")
            os.remove(file_name)


def process_samples(file_list):
    # Initialize csv file
    with open(os.path.join(__out_dir__, __out_file__), "w+") as o:
        samples = map(lambda sam: os.path.splitext(os.path.basename(sam))[0], file_list)
        o.write("id"+__out_delimiter__+__out_delimiter__.join(samples)+"\n")
    # Create smaller files with chromosome-specific junctions
    print("[INFO] Splitting BED files by chromosome...")
    chromosomes_found = set()
    for f in file_list:
        print("\t[INFO] Splitting sample: " + os.path.splitext(os.path.basename(f))[0])
        chromosomes_found.add(split_csv(f, __temp_dir__))
        print("\t[DONE]")
    print("[DONE]")
    if len(chromosomes_found) == 0:
        e_print('\t[ERROR] No chromosomes found in junction files')
        sys.exit(-1)
    chromosomes = sorted(list(chromosomes_found))
    # Per chromosome, process junctions in file and add them to .csv file
    print("[INFO] Collecting junction information by chromosome...")
    read_junctions(file_list, chromosomes)
    print("[DONE]")
    print("[INFO] Cleaning temporary files...")
    clean_files(samples, chromosomes)
    print("[DONE]")
    print("[INFO] Finished without (apparent) errors")


if not os.path.exists(__temp_dir__):
    os.makedirs(__temp_dir__)
if not os.path.exists(__out_dir__):
    os.makedirs(__out_dir__)
process_samples(sorted(sys.argv[1:]))
