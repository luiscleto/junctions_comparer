#!/usr/bin/python
import errno

from bed_features import BEDIndices
from utils import *
from gtf_features import *

UNKNOWN_GENE_ID = "UNKNOWN"

__temp_dir__ = "sj_tmp"
__out_dir__ = "sj_out"
__out_file__ = "results.csv"
__out_delimiter__ = ","
__gene_list_delimiter__ = "|"
## strand -> chromosome -> gene_id -> GeneDescription
__gene_set__ = {'+': {}, '-': {}}
## strand -> chromosome -> ([pos_start], [pos_end], [gene_id]) (tuple of lists (sorted by pos_start and secondly by pos_end) to allow easy binary search. Same index in each list corresponds to same gene)
__gene_intervals__ = {'+': {}, '-': {}}
## sample -> unknown_read_count
__unknown_reads_per_sample__ = {}
## strand -> chromosome -> sort_type -> ([pos_start], [pos_end]) (tuple of lists to allow easy binary search. Same index in each list corresponds to same exon.
__exon_list___ = {'+': {}, '-': {}}
__sorted_by_start_key___ = 'sort_start'
__sorted_by_end_key___ = 'sort_end'


if len(sys.argv) < 4:
    e_print('[ERROR] Insufficient parameters provided!')
    print('Usage:')
    print('\tpython junctions_comparer.py <ref_genome.gtf> <sample1.bed> <sample2.bed> [<sample3.bed> ...]')


class __SpliceTypes:
    def __init__(self):
        pass
    canonical = 0
    non_canonical = 1

SpliceTypes = __SpliceTypes()


def find_gene(position_s, chromosome, strand):
    if chromosome not in __gene_intervals__[strand]:
        return [UNKNOWN_GENE_ID]
    gene_set = []
    try:
        i = find_le(__gene_intervals__[strand][chromosome][0], position_s)
    except ValueError:
        return [UNKNOWN_GENE_ID]
    while __gene_intervals__[strand][chromosome][0][i] <= position_s and i >= 0:
        s = __gene_intervals__[strand][chromosome][0][i]
        e = __gene_intervals__[strand][chromosome][1][i]
        gene_id = __gene_intervals__[strand][chromosome][2][i]
        i -= 1
        if s <= position_s <= e:
            gene_set.append(gene_id)
        elif len(gene_set) > 0:
            break
    if len(gene_set) == 0:
        return [UNKNOWN_GENE_ID]
    return gene_set


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
                if chrom not in __exon_list___[row[GTFIndices.strand]]:
                    __exon_list___[row[GTFIndices.strand]][chrom] = {__sorted_by_start_key___: list(), __sorted_by_end_key___: list()}
                __exon_list___[row[GTFIndices.strand]][chrom][__sorted_by_start_key___].append((int(row[GTFIndices.start]), int(row[GTFIndices.end])))
                __exon_list___[row[GTFIndices.strand]][chrom][__sorted_by_end_key___].append((int(row[GTFIndices.start]), int(row[GTFIndices.end])))
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
                __exon_list___[strand][chromosome][__sorted_by_start_key___] = zip(*sorted(__exon_list___[strand][chromosome][__sorted_by_start_key___], key=lambda x: x[0]))
                __exon_list___[strand][chromosome][__sorted_by_end_key___] = zip(*sorted(__exon_list___[strand][chromosome][__sorted_by_end_key___], key=lambda x: x[1]))
            for chromosome in __gene_intervals__[strand].keys():
                __gene_intervals__[strand][chromosome] = zip(*sorted(__gene_intervals__[strand][chromosome], key=lambda x: (x[0], x[1])))
        print("\t[DONE]")
    print('[DONE]')


def write_chromosome_junctions_to_file(junctions):
    writer = csv.writer(open(os.path.join(__out_dir__, __out_file__), 'a'), delimiter=__out_delimiter__, lineterminator='\n')
    for key, value in junctions.items():
        row = [key]
        row.extend(value)
        writer.writerow(row)

def find_next_exon(chrom, strand, pos_end):
    try:
        i = find_gt(__exon_list___[strand][chrom][__sorted_by_start_key___][0], pos_end)
        return __exon_list___[strand][chrom][__sorted_by_start_key___][0][i], __exon_list___[strand][chrom][__sorted_by_start_key___][1][i]
    except ValueError:
        return -1, -1

def determine_splice_type(chrom, strand, junc_start, junc_end):
    if chrom not in __exon_list___[strand]:
        e_print("\t\t\t[WARNING] Chromosome " + chrom + " found in BED file but not in GTF!")
        # to avoid "spamming" warning
        __exon_list___[strand][chrom][__sorted_by_start_key___] = ([],[])
        __exon_list___[strand][chrom][__sorted_by_end_key___] = ([],[])
        return SpliceTypes.non_canonical

    try:
        i = index(__exon_list___[strand][chrom][__sorted_by_end_key___][1], junc_start)
    except ValueError:
        return SpliceTypes.non_canonical
    while __exon_list___[strand][chrom][__sorted_by_end_key___][1][i] == junc_start:
        s = __exon_list___[strand][chrom][__sorted_by_end_key___][0][i]
        e = __exon_list___[strand][chrom][__sorted_by_end_key___][1][i]
        i += 1
        if e < junc_start:
            continue
        elif e == junc_start:
            next_ex_start, next_ex_end = find_next_exon(chrom, strand, e)
            if next_ex_start == -1 and next_ex_end == -1:
                continue # TODO handle this special case (no known exons after junction start exon)
            if junc_end+1 == next_ex_start:
                return SpliceTypes.canonical
            # TODO detect types of non canonical junctions
        else:
            break
    return SpliceTypes.non_canonical


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
    total_reads_per_sample = map_to_dict(lambda x:x, lambda _:0, samples)
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
                        total_reads_per_sample[sample] += count
    for sample, count in __unknown_reads_per_sample__.items():
        if count > 0:
            file_handlers[sample].write("UNKNOWN=" + str(count) + "\n")
            total_reads_per_sample[sample] += count
    for k, f in file_handlers.items():
        f.write("TOTAL_READS="+str(total_reads_per_sample[k])+"\n")
        f.close()


def process_samples(file_list):
    # Initialize csv file
    samples = map(lambda sam: os.path.splitext(os.path.basename(sam))[0], file_list)
    global __unknown_reads_per_sample__
    __unknown_reads_per_sample__ = map_to_dict(lambda x:x, lambda _:0, samples)
    with open(os.path.join(__out_dir__, __out_file__), "w+") as o:
        o.write("id"+__out_delimiter__+__out_delimiter__.join(samples)
                + __out_delimiter__+"gene_ids (start)"
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
read_reference_genome(sys.argv[1], map(lambda sam: os.path.splitext(os.path.basename(sam))[0],sys.argv[2:]))
process_samples(sorted(sys.argv[2:]))
