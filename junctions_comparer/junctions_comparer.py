#!/usr/bin/python
import errno

from bed_features import BEDIndices
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
## strand -> chromosome -> ([pos_start], [pos_end]) (tuple of lists to allow easy binary search. Same index in each list corresponds to same exon.
__exon_list___ = {'+': {}, '-': {}}
## strand -> chromosome -> ([pos_start], [pos_end], [(pos_start, pos_end)]) (last element is a list of all possible exons for that interval (due to overlaps))
__exon_disjoint_intervals__ = {'+': {}, '-': {}}
__total_reads_per_sample__ = {}


parser = argparse.ArgumentParser(description="Junction files analyser", usage='''junctions_comparer.py [-h] [-nb | -q] [-n MIN_READS] [-t TEMP_DIR]
                             [-o OUT_DIR] [-r RESULTS_FILE]
                             [-rl RESULTS_DELIMITER] [-u UNKNOWN_ID]
                             [-gl GENE_DELIMITER]
                             gtf_file bed_file bed_file [bed_file ...]
''')
group = parser.add_mutually_exclusive_group()
group.add_argument("-nb", "--no-bars", help="disables loading bars (use when redirecting output to file)", action="store_true")
group.add_argument("-q", "--quiet", action="store_true")
parser.add_argument("-n", "--min-reads", type=int, help="minimum number of reads supporting junction to pass filter (default: 3)")
parser.add_argument("-t", "--temp-dir", help="directory where to keep temporary files (default: sj_tmp)")
parser.add_argument("-o", "--out-dir", help="directory where to keep output files (default: sj_out)")
parser.add_argument("-r", "--results-file", help="name for final result file (default: results.csv)")
parser.add_argument("-rl", "--results-delimiter", help="delimiter for output CSV file (default: ,)")
parser.add_argument("-u", "--unknown-id", help="name for unknown genes (default: UNKNOWN)")
parser.add_argument("-gl", "--gene-delimiter", help="delimiter for gene lists in output (default: |)NOTE: MUST NOT BE THE SAME AS -rl")
parser.add_argument("gtf_file", help="reference annotation file in GTF format")
parser.add_argument('bed_file1', nargs=1, metavar='bed_file', help="junction annotation file in BED format")
parser.add_argument('bed_file2', nargs='+', metavar='bed_file', help=argparse.SUPPRESS)

parser.set_defaults(min_reads=3, temp_dir="sj_tmp", out_dir="sj_out", results_file="results.csv", results_delimiter=",", unknown_id="UNKNOWN", gene_delimiter="|")
args = parser.parse_args()

UNKNOWN_GENE_ID = args.unknown_id
__temp_dir__ = args.temp_dir
__out_dir__ = args.out_dir
__out_file__ = args.results_file
__out_delimiter__ = args.results_delimiter
__gene_list_delimiter__ = args.gene_delimiter
__min_reads__ = 3

__out_file_filtered__ = __out_file__.rsplit('.', 1)[0] + '.filtered.' + __out_file__.rsplit('.', 1)[1]

bed_files = args.bed_file1
bed_files.extend(args.bed_file2)

class __SpliceTypes:
    def __init__(self):
        pass
    canonical = "C"
    exon_skip = "ES"
    intron_inclusion = "IN"
    alt_3_prime = "3'"
    alt_5_prime = "5'"
    unknown = "UK"


SpliceTypes = __SpliceTypes()


def find_gene(position_s, chromosome, strand):
    if chromosome not in __gene_disjoint_intervals__[strand]:
        return [UNKNOWN_GENE_ID]

    try:
        i = find_le(__gene_disjoint_intervals__[strand][chromosome][0], position_s)
    except ValueError:
        return [UNKNOWN_GENE_ID]

    gene_set = filter(lambda g_id: __gene_set__[strand][chromosome][g_id].start_position <= position_s <= __gene_set__[strand][chromosome][g_id].end_position,
                      __gene_disjoint_intervals__[strand][chromosome][2][i])

    if not gene_set:
        return [UNKNOWN_GENE_ID]

    return gene_set


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


def generate_disjoint_exon_intervals(strand, chromosome):
    global __exon_disjoint_intervals__
    __exon_disjoint_intervals__[strand][chromosome] = list()
    current_start = -1
    current_exons = list()
    for s, e in __exon_list___[strand][chromosome]:
        if current_start == -1:
            current_start = s
            current_end = e
            current_exons.append((s,e))
        elif s == current_exons[-1][0] and e == current_exons[-1][1]:
            continue # ignore duplicates
        elif s <= current_end:
            current_exons.append((s,e))
            if e > current_end:
                current_end = e
        else:
            __exon_disjoint_intervals__[strand][chromosome].append((current_start, current_end, current_exons))
            current_start = s
            current_end = e
            current_exons = [(s,e)]
    if current_start != -1:
        __exon_disjoint_intervals__[strand][chromosome].append((current_start, current_end, current_exons))
    __exon_disjoint_intervals__[strand][chromosome] = zip(*__exon_disjoint_intervals__[strand][chromosome])


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
            elif str(row[GTFIndices.feature] == 'exon'):
                chrom = str(row[GTFIndices.seq_name]).replace("chr", "")
                if len(chrom) == 1 and represents_int(chrom):
                    chrom = "0" + chrom
                elif chrom == "M":
                    chrom = "MT"

                if chrom not in __exon_list___[row[GTFIndices.strand]]:
                    __exon_list___[row[GTFIndices.strand]][chrom] = list()
                __exon_list___[row[GTFIndices.strand]][chrom].append((int(row[GTFIndices.start]), int(row[GTFIndices.end])))
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
        print("\t[INFO] Sorting gene and exon lists...")
        for strand in ['+', '-']:
            for chromosome in __exon_list___[strand].keys():
                __exon_list___[strand][chromosome] = sorted(__exon_list___[strand][chromosome], key=lambda x: x[0])
                generate_disjoint_exon_intervals(strand, chromosome)
                __exon_list___[strand][chromosome] = zip(*__exon_list___[strand][chromosome])
            for chromosome in __gene_intervals__[strand].keys():
                __gene_intervals__[strand][chromosome] = sorted(__gene_intervals__[strand][chromosome], key=lambda x: (x[0], x[1]))
                generate_disjoint_gene_intervals(strand, chromosome)
                del __gene_intervals__[strand][chromosome]
            __gene_intervals__[strand].clear()
        print("\t[DONE]")
    print('[DONE]')


def write_chromosome_junctions_to_file(junctions):
    writer = csv.writer(open(os.path.join(__out_dir__, __out_file__), 'a'), delimiter=__out_delimiter__, lineterminator='\n')
    writer2 = csv.writer(open(os.path.join(__out_dir__, __out_file_filtered__), 'a'), delimiter=__out_delimiter__, lineterminator='\n')
    for key, value in junctions.items():
        row = [key]
        row.extend(value)
        writer.writerow(row)
        if any(v >= __min_reads__ for v in value[:-3]):
            writer2.writerow(row)


def find_next_exon(chrom, strand, pos_end):
    try:
        i = find_gt(__exon_list___[strand][chrom][0], pos_end)
        return __exon_list___[strand][chrom][0][i], __exon_list___[strand][chrom][1][i]
    except ValueError:
        return -1, -1


def find_next_exons(chrom, strand, pos_end):
    # try:
    s,e = find_next_exon(chrom, strand, pos_end)
    if s == -1 and e == -1:
        return [(s,e)]

    i = find_le(__exon_disjoint_intervals__[strand][chrom][0], e)
    return __exon_disjoint_intervals__[strand][chrom][2][i]
    # except ValueError: THIS SHOULD NOT BE POSSIBLE
    #     return [(-1, -1)]


def find_splice_for_end_junction(chrom, strand, junc_end, possible_first_exons):
    exon_start = junc_end+1
    three_prime_found = False
    for s,e in possible_first_exons: # find canonical first
        possible_nexts = filter(lambda coords: coords[0] < exon_start <= coords[1],find_next_exons(chrom,strand,e))
        if (not possible_nexts) or (len(possible_nexts) == 1 and possible_nexts[0][0] == -1):
            continue
        exact_nexts = filter(lambda coords: exon_start == coords[1], possible_nexts)
        if exact_nexts:
            return SpliceTypes.canonical
        three_prime_nexts = filter(lambda coords: exon_start > coords[0], possible_nexts)
        if three_prime_nexts:
            three_prime_found = True

    if three_prime_found:
        return SpliceTypes.alt_3_prime

    # Try to find exon skip event or handle cases where first event is intron_inclusion
    prefix = (SpliceTypes.exon_skip + __gene_list_delimiter__) if possible_first_exons else ""

    possible_nexts = filter(lambda coords: coords[0] <= exon_start < coords[1], find_next_exons(chrom, strand, junc_end))
    if (not possible_nexts) or (len(possible_nexts) == 1 and possible_nexts[0][0] == -1):
        return SpliceTypes.intron_inclusion
    exact_nexts = filter(lambda coords: exon_start == coords[0], possible_nexts)
    if exact_nexts:
        return prefix + SpliceTypes.canonical
    three_prime_nexts = filter(lambda coords: exon_start < coords[1], possible_nexts)
    if three_prime_nexts:
        return prefix + SpliceTypes.alt_3_prime

    return SpliceTypes.intron_inclusion


def determine_splice_type(chrom, strand, junc_start, junc_end):
    if chrom not in __exon_disjoint_intervals__[strand]:
        e_print("\t\t\t[WARNING] Chromosome " + chrom + " found in BED file but not in GTF!")
        # to avoid "spamming" warning
        __exon_list___[strand][chrom] = ([], [])
        __exon_disjoint_intervals__[strand][chrom] = ([],[],[])
        return SpliceTypes.intron_inclusion + __gene_list_delimiter__ + SpliceTypes.intron_inclusion

    try:
        i = find_le(__exon_disjoint_intervals__[strand][chrom][0], junc_start)
        possible_exon_set = filter(lambda coords: coords[0] < junc_start <= coords[1], __exon_disjoint_intervals__[strand][chrom][2][i])
        precise_exons = filter(lambda coords: junc_start == coords[1], possible_exon_set)
        if precise_exons:
            return SpliceTypes.canonical + __gene_list_delimiter__ + find_splice_for_end_junction(chrom, strand, junc_end, precise_exons)
        five_prime_exons = filter(lambda coords: junc_start > coords[0], possible_exon_set)
        if five_prime_exons:
            return SpliceTypes.alt_5_prime + __gene_list_delimiter__ + find_splice_for_end_junction(chrom, strand, junc_end, five_prime_exons)
    except ValueError:
        pass

    return SpliceTypes.intron_inclusion + __gene_list_delimiter__ + find_splice_for_end_junction(chrom, strand, junc_end, [])


def read_junctions(samples, chromosomes):
    for c in chromosomes:
        print('\t[INFO] Collecting info for chromosome ' + c + "...")
        chromosome_junctions = {}
        sample_index = 0
        for s in map(lambda sam: os.path.splitext(os.path.basename(sam))[0], samples):
            print('\t\tProcessing sample: ' + s + "...")
            file_name = os.path.join(__temp_dir__, s + "_" + c + ".bed")
            if not os.path.exists(file_name):
                e_print("\t\t\t[WARNING] Sample " + s + " does not contain any junction reads for chromosome " + c + "!")
            else:
                num_rows = file_len(file_name)
                __total_reads_per_sample__[s] += num_rows
                with open(file_name, 'rb') as csv_file:
                    csv_reader = csv.reader(csv_file, delimiter='\t')
                    row_no = 0
                    for row in csv_reader:
                        if row_no % 10000 == 0:
                            print_progress(row_no, num_rows, '\t\t\t')
                        block_sizes = map(int, str(row[BEDIndices.block_sizes]).split(","))
                        true_start = int(row[BEDIndices.start]) + block_sizes[0]
                        true_end = int(row[BEDIndices.end]) - block_sizes[1]
                        junc_id = row[BEDIndices.chromosome] + "_" + row[BEDIndices.strand] + "_" + str(true_start) + "_" + str(true_end)
                        if junc_id not in chromosome_junctions:
                            chromosome_junctions[junc_id] = [0] * len(samples)
                            gen1_ids = find_gene(true_start, row[BEDIndices.chromosome], row[BEDIndices.strand])
                            gen2_ids = find_gene(true_end, row[BEDIndices.chromosome], row[BEDIndices.strand])
                            chromosome_junctions[junc_id].extend([__gene_list_delimiter__.join(gen1_ids),
                                                                  __gene_list_delimiter__.join(gen2_ids),
                                                                  determine_splice_type(row[BEDIndices.chromosome], row[BEDIndices.strand], true_start, true_end)])
                        gen1_ids = chromosome_junctions[junc_id][len(samples)].split(__gene_list_delimiter__)
                        gen2_ids = chromosome_junctions[junc_id][len(samples)+1].split(__gene_list_delimiter__)
                        for gen1_id in gen1_ids:
                            if gen1_id != UNKNOWN_GENE_ID:
                                __gene_set__[row[BEDIndices.strand]][row[BEDIndices.chromosome]][gen1_id].gene_reads_by_sample[s] += 1
                            else:
                                __unknown_reads_per_sample__[s] += 1
                        for gen2_id in filter(lambda x: x not in gen1_ids, gen2_ids):
                            if gen2_id != UNKNOWN_GENE_ID:
                                __gene_set__[row[BEDIndices.strand]][row[BEDIndices.chromosome]][gen2_id].gene_reads_by_sample[s] += 1
                            else:
                                __unknown_reads_per_sample__[s] += 1
                        chromosome_junctions[junc_id][sample_index] += 1
                        row_no += 1
                    print_progress(num_rows, num_rows, '\t\t\t')
            sample_index += 1
        write_chromosome_junctions_to_file(chromosome_junctions)
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
    file_handlers = map_to_dict(lambda x:x, lambda y:open(os.path.join(__out_dir__, y + ".read_counts.txt"),"w+"), samples)
    for k, f in file_handlers.items():
        f.write("##Summary for sample " + k + "\n")
        f.write("##This file indicates the number of reads in this sample (of junctions) that mapped to a certain gene as well as the TOTAL_READS value\n")
        f.write("##Format: <gene_id>=<number_of_reads>\n")
    for strand, chromosome_to_id in __gene_set__.items():
        for id_to_gene in chromosome_to_id.values():
            for gene_id, gen in id_to_gene.items():
                for sample, count in gen.gene_reads_by_sample.items():
                    if count > 0:
                        file_handlers[sample].write(gene_id+"="+str(count) + "\n")
    for sample, count in __unknown_reads_per_sample__.items():
        if count > 0:
            file_handlers[sample].write("UNKNOWN=" + str(count) + "\n")
    for k, f in file_handlers.items():
        f.write("TOTAL_READS="+str(__total_reads_per_sample__[k])+"\n")
        f.close()


def process_samples(file_list):
    # Initialize csv file
    samples = map(lambda sam: os.path.splitext(os.path.basename(sam))[0], file_list)
    global __total_reads_per_sample__
    __total_reads_per_sample__ = map_to_dict(lambda x:x, lambda _:0, samples)
    global __unknown_reads_per_sample__
    __unknown_reads_per_sample__ = map_to_dict(lambda x:x, lambda _:0, samples)
    with open(os.path.join(__out_dir__, __out_file__), "w+") as o:
        o.write("id"+__out_delimiter__+__out_delimiter__.join(samples)
                + __out_delimiter__+"gene_ids (start)"
                + __out_delimiter__ + "gene_ids (end)"
                + __out_delimiter__ + "type\n")
    with open(os.path.join(__out_dir__, __out_file__), "w+") as o:
        o.write("id" + __out_delimiter__ + __out_delimiter__.join(samples)
                + __out_delimiter__ + "gene_ids (start)"
                + __out_delimiter__ + "gene_ids (end)"
                + __out_delimiter__ + "type\n")
    with open(os.path.join(__out_dir__, __out_file_filtered__), "w+") as o:
        o.write("id" + __out_delimiter__ + __out_delimiter__.join(samples)
                + __out_delimiter__ + "gene_ids (start)"
                + __out_delimiter__ + "gene_ids (end)"
                + __out_delimiter__ + "type\n")
    # Create smaller files with chromosome-specific junctions
    print("[INFO] Splitting BED files by chromosome...")
    chromosomes_found = set()
    for f in file_list:
        print("\t[INFO] Splitting sample: " + os.path.splitext(os.path.basename(f))[0])
        chromosomes_found |= split_csv(f, __temp_dir__)
        print("\t[DONE]")
    if len(chromosomes_found) == 0:
        e_print('\t[ERROR] No chromosomes found in junction files')
        sys.exit(-1)
    print("[DONE]")
    chromosomes = sorted(list(chromosomes_found))
    #chromosomes = ['01']
    # Per chromosome, process junctions in file and add them to .csv file
    print("[INFO] Collecting junction information by chromosome...")
    read_junctions(file_list, chromosomes)
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
read_reference_genome(args.gtf_file, map(lambda sam: os.path.splitext(os.path.basename(sam))[0],bed_files))
process_samples(sorted(bed_files))
