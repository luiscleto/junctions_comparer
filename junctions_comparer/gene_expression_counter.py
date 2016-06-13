#!/usr/bin/python
import errno

from utils import *
from gtf_features import *
import argparse


## strand -> chromosome -> gene_id -> GeneDescription
__gene_set__ = {'+': {}, '-': {}}
## strand -> chromosome -> [(pos_start, pos_end, gene_id)] (list of gene id and their positions) (deleted after computing disjoint intervals)
__gene_intervals__ = {'+': {}, '-': {}}
## strand -> chromosome -> ([pos_start], [pos_end], [[gene_ids]]) (last element is a list of all possible genes for that interval (due to overlaps))
## it is a (sorted) tuple of lists to facilitate binary search
__gene_disjoint_intervals__ = {'+': {}, '-': {}}
## sample -> unknown_read_count
__unknown_reads_per_sample__ = {}
## sample -> total_read_count
__total_reads_per_sample__ = {}


parser = argparse.ArgumentParser(description="Gene read counter", usage='''gene_expression_counter.py [-h] [-nb | -q] [-t TEMP_DIR] [-o OUT_DIR]
                                  [-rl RESULTS_DELIMITER] [-u UNKNOWN_ID]
                                  [-a {stranded,unstranded}]
                                  gtf_file reads_file reads_file [reads_file ...]
''')
group = parser.add_mutually_exclusive_group()
group.add_argument("-nb", "--no-bars", help="disables loading bars (use when redirecting output to file)", action="store_true")
group.add_argument("-q", "--quiet", action="store_true")
parser.add_argument("-t", "--temp-dir", help="directory where to keep temporary files (default: gec_tmp)")
parser.add_argument("-o", "--out-dir", help="directory where to keep output files (default: gec_out)")
parser.add_argument("-rl", "--results-delimiter", help="delimiter for output TSV file (default: \\t)")
parser.add_argument("-u", "--unknown-id", help="name for unknown genes (default: UNKNOWN)")
parser.add_argument("-a", "--analysis", choices=['stranded','unstranded'], help="type of analysis (default: unstranded)")
parser.add_argument("gtf_file", help="reference annotation file in GTF format")
parser.add_argument('reads_file1', nargs=1, metavar='reads_file', help="read annotation file in custom format (chromosome<tab>read_start<tab>read_end<tab>strand)")
parser.add_argument('reads_file2', nargs='+', metavar='reads_file', help=argparse.SUPPRESS)

parser.set_defaults(temp_dir="gec_tmp", out_dir="gec_out", results_delimiter="\t", unknown_id="UNKNOWN", analysis="unstranded")
args = parser.parse_args()

UNKNOWN_GENE_ID = args.unknown_id
__temp_dir__ = args.temp_dir
__out_dir__ = args.out_dir
__out_delimiter__ = args.results_delimiter
__out_extension__ = ".gene_expression.tsv"
__stranded_analysis__ = args.analysis == "stranded"


reads_file = args.reads_file1
reads_file.extend(args.reads_file2)


def generate_disjoint_gene_intervals(strand, chromosome):
    global __gene_disjoint_intervals__
    __gene_disjoint_intervals__[strand][chromosome] = list()
    current_start=-1
    current_ids=list()
    for s,e,g_id in __gene_intervals__[strand][chromosome]:
        if current_start == -1:
            current_start = s
            current_end = e
            current_ids.append(g_id)
        elif s <= current_end:
            current_ids.append(g_id)
            if e > current_end:
                current_end = e
        else:
            __gene_disjoint_intervals__[strand][chromosome].append((current_start, current_end, current_ids))
            current_start = s
            current_end = e
            current_ids = [g_id]
    if current_start != -1:
        __gene_disjoint_intervals__[strand][chromosome].append((current_start, current_end, current_ids))
    __gene_disjoint_intervals__[strand][chromosome] = zip(*__gene_disjoint_intervals__[strand][chromosome])


def read_reference_genome(filename, sample_list):
    global __gene_set__
    global __gene_intervals__
    print('[INFO] Reading GTF file...')
    num_rows = file_len(filename)
    with open(filename, 'rb') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        row_count = 0
        for row in csv_reader:
            if row_count % 1000 == 0:
                print_progress(row_count, num_rows, '\t')
            row_count += 1
            if len(row) > 0 and str(row[0]).startswith("#"):
                continue
            elif len(row) == 0:
                continue
            if str(row[GTFIndices.feature]) != 'gene':
                continue
            gen = GeneDescription(row, sample_list)
            if gen.chromosome not in __gene_set__[gen.strand]:
                __gene_set__[gen.strand][gen.chromosome] = {}
                __gene_intervals__[gen.strand][gen.chromosome] = list()
            if gen.id not in __gene_set__[gen.strand][gen.chromosome]:
                __gene_set__[gen.strand][gen.chromosome][gen.id] = gen
                __gene_intervals__[gen.strand][gen.chromosome].append((gen.start_position, gen.end_position, gen.id))
            elif __gene_set__[gen.strand][gen.chromosome][gen.id].start_position != gen.start_position \
                    or __gene_set__[gen.strand][gen.chromosome][gen.id].end_position != gen.end_position:
                e_print("\t[WARNING] Duplicate gene ID found in GTF file with different positions: " + gen.id)
            else:
                e_print("\t[WARNING] Duplicate gene ID found in GTF file: " + gen.id)
        print_progress(num_rows, num_rows, '\t')
        print("\t[INFO] Sorting gene lists...")
        for strand in ['+', '-']:
            for chromosome in __gene_intervals__[strand].keys():
                __gene_intervals__[strand][chromosome] = sorted(__gene_intervals__[strand][chromosome], key=lambda x: (x[0], x[1]))
                generate_disjoint_gene_intervals(strand, chromosome)
                del __gene_intervals__[strand][chromosome]
            __gene_intervals__[strand].clear()
        print("\t[DONE]")
    print('[DONE]')


def find_included_genes(start, end, chromosome, strand):
    if chromosome not in __gene_disjoint_intervals__[strand]:
        return [UNKNOWN_GENE_ID]

    try:
        i = find_le(__gene_disjoint_intervals__[strand][chromosome][0], start)
    except ValueError:
        return [UNKNOWN_GENE_ID]

    try:
        last_i = find_le(__gene_disjoint_intervals__[strand][chromosome][0], end)
    except ValueError:
        last_i = len(__gene_disjoint_intervals__[strand][chromosome][0]) - 1

    gene_set = []
    while i <= last_i:
        gene_set.extend(filter(lambda g_id: max(__gene_set__[strand][chromosome][g_id].start_position, start) <= min(__gene_set__[strand][chromosome][g_id].end_position, end),
                          __gene_disjoint_intervals__[strand][chromosome][2][i]))
        i += 1

    if not gene_set:
        return [UNKNOWN_GENE_ID]

    return gene_set


def write_read_counts_to_file(chromosome_gene_read_count, s, c):
    writer = csv.writer(open(os.path.join(__out_dir__, s + __out_extension__), 'a'), delimiter=__out_delimiter__,lineterminator='\n')
    for key, value in chromosome_gene_read_count.items():
        row = [key, value, c]
        writer.writerow(row)


def process_reads(samples, chromosomes):
    for c in chromosomes:
        print('\t[INFO] Collecting info for chromosome ' + c + "...")
        sample_index = 0
        for s in map(lambda sam: os.path.splitext(os.path.basename(sam))[0], samples):
            print('\t\tProcessing sample: ' + s + "...")
            file_name = os.path.join(__temp_dir__, s + "_" + c + ".bed")
            if not os.path.exists(file_name):
                e_print("\t\t\t[WARNING] Sample " + s + " does not contain any reads for chromosome " + c + "!")
            else:
                chromosome_gene_read_count = {}
                num_rows = file_len(file_name)
                __total_reads_per_sample__[s] += num_rows
                with open(file_name, 'rb') as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter='\t')
                    row_no = 0
                    for row in csv_reader:
                        if row_no % 10000 == 0:
                            print_progress(row_no, num_rows, '\t\t\t')
                        true_start = int(row[1])
                        true_end = int(row[2])
                        strand = str(row[3])
                        if __stranded_analysis__:
                            gene_set = find_included_genes(true_start, true_end, c, strand)
                        else:
                            gene_set = filter(lambda x: x!= UNKNOWN_GENE_ID,find_included_genes(true_start, true_end, c, "+"))
                            gene_set.extend(filter(lambda x: x!= UNKNOWN_GENE_ID,find_included_genes(true_start, true_end, c, "-")))
                            if not gene_set:
                                gene_set = [UNKNOWN_GENE_ID]
                        for g_id in gene_set:
                            if g_id == UNKNOWN_GENE_ID:
                                __unknown_reads_per_sample__[s] += 1
                            elif g_id not in chromosome_gene_read_count:
                                chromosome_gene_read_count[g_id] = 1
                            else:
                                chromosome_gene_read_count[g_id] += 1
                        row_no += 1
                    print_progress(num_rows, num_rows, '\t\t\t')
                    write_read_counts_to_file(chromosome_gene_read_count, s, c)
            sample_index += 1
        print('\t[DONE]')


def clean_files(samples, chromosomes):
    for c in chromosomes:
        for s in samples:
            file_name = os.path.join(__temp_dir__, s + "_" + c + ".bed")
            if os.path.exists(file_name):
                os.remove(file_name)
    try:
        os.rmdir(__temp_dir__)
    except OSError as ex:
        if ex.errno == errno.ENOTEMPTY:
            e_print("\t[WARNING] Could not delete temporary directory! Directory not empty.")


def write_summary(samples):
    for s in samples:
        writer = csv.writer(open(os.path.join(__out_dir__, s + __out_extension__), 'a'), delimiter=__out_delimiter__,
                            lineterminator='\n')
        row1 = [UNKNOWN_GENE_ID, __unknown_reads_per_sample__[s], "_"]
        row2 = ["TOTAL", __total_reads_per_sample__[s], "_"]
        writer.writerow(row1)
        writer.writerow(row2)


def process_samples(file_list):
    # Initialize csv file
    samples = map(lambda sam: os.path.splitext(os.path.basename(sam))[0], file_list)
    global __total_reads_per_sample__
    __total_reads_per_sample__ = map_to_dict(lambda x:x, lambda _:0, samples)
    global __unknown_reads_per_sample__
    __unknown_reads_per_sample__ = map_to_dict(lambda x:x, lambda _:0, samples)

    for s in samples:
        with open(os.path.join(__out_dir__, s + __out_extension__), "w+") as o:
            o.write("gene_id" + __out_delimiter__ + "read_count" + __out_delimiter__ + "chr\n")

    # Create smaller files with chromosome-specific junctions
    print("[INFO] Splitting read files by chromosome...")
    chromosomes_found = set()
    for f in file_list:
        print("\t[INFO] Splitting sample: " + os.path.splitext(os.path.basename(f))[0])
        chromosomes_found |= split_csv(f, __temp_dir__)
        print("\t[DONE]")
    if len(chromosomes_found) == 0:
        e_print('\t[ERROR] No chromosomes found in read files')
        sys.exit(-1)
    print("[DONE]")
    chromosomes = sorted(list(chromosomes_found))
    # Per chromosome, process junctions in file and add them to .csv file
    print("[INFO] Collecting read information by chromosome...")
    process_reads(file_list, chromosomes)
    print("[DONE]")
    # Save summary for each sample
    print("[INFO] Summarizing read counts...")
    write_summary(samples)
    print("[DONE]")
    # Delete temporary files
    print("[INFO] Cleaning temporary files...")
    clean_files(samples, chromosomes)
    print("[DONE]")
    print("[INFO] Finished without (apparent) errors")


if not os.path.exists(__temp_dir__):
    os.makedirs(__temp_dir__)
if not os.path.exists(__out_dir__):
    os.makedirs(__out_dir__)
read_reference_genome(args.gtf_file, map(lambda sam: os.path.splitext(os.path.basename(sam))[0],reads_file))
process_samples(sorted(reads_file))
