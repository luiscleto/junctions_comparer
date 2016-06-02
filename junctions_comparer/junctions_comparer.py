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
## chromosome -> gene_id -> GeneDescription
__gene_set__ = {}
## chromosome -> [(pos_start, pos_end, gene_id)]
__gene_intervals__ = {}
## unknown_read_count
__unknown_reads_per_sample__ = {}

if len(sys.argv) < 4:
    e_print('[ERROR] Insufficient parameters provided!')
    print('Usage:')
    print('\tpython junctions_comparer.py <ref_genome.gtf> <sample1.bed> <sample2.bed> [<sample3.bed> ...]')


def find_gene(positionS, chromosome):
    if chromosome not in __gene_intervals__:
        return UNKNOWN_GENE_ID
    for s,e,gene_name in __gene_intervals__[chromosome]:
        if s <= positionS < e:
            return gene_name
    return UNKNOWN_GENE_ID


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
            elif str(row[GTFIndices.feature]) != 'gene':
                continue
            gen = GeneDescription(row, sample_list)
            if gen.chromosome not in __gene_set__:
                __gene_set__[gen.chromosome] = {}
                __gene_intervals__[gen.chromosome] = list()
            if gen.id not in __gene_set__[gen.chromosome]:
                __gene_set__[gen.chromosome][gen.id] = gen
                __gene_intervals__[gen.chromosome].append((gen.start_position, gen.end_position, gen.id))
            elif __gene_set__[gen.chromosome][gen.id].start_position != gen.start_position \
                    or __gene_set__[gen.chromosome][gen.id].end_position != gen.end_position:
                e_print("\t[WARNING] Duplicate gene ID found in GTF file with different positions: " + gen.id)
            else:
                e_print("\t[WARNING] Duplicate gene ID found in GTF file: " + gen.id)
        print_progress(num_rows, num_rows, '\t')
    print('[DONE]')


def write_chromosome_junctions_to_file(junctions):
    writer = csv.writer(open(os.path.join(__out_dir__, __out_file__), 'a'), delimiter=__out_delimiter__, lineterminator='\n')
    for key, value in junctions.items():
        row = [key]
        row.extend(value)
        writer.writerow(row)


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
                        junc_id = row[BEDIndices.chromosome] + "_" + row[BEDIndices.start] + "_" + row[BEDIndices.end]
                        if junc_id not in chromosome_junctions:
                            chromosome_junctions[junc_id] = [0] * len(samples)
                            gen1_id = find_gene(row[BEDIndices.start], row[BEDIndices.chromosome])
                            gen2_id = find_gene(row[BEDIndices.end], row[BEDIndices.chromosome])
                            if gen1_id != ("%s" % UNKNOWN_GENE_ID):
                                gen1_name = __gene_set__[row[BEDIndices.chromosome]][gen1_id].name
                            else:
                                gen1_name = gen1_id
                            if gen2_id != ("%s" % UNKNOWN_GENE_ID):
                                gen2_name = __gene_set__[row[BEDIndices.chromosome]][gen2_id].name
                            else:
                                gen2_name = gen2_id
                            chromosome_junctions[junc_id].extend([gen1_name, gen1_id, gen2_name, gen2_id])
                        gen1_id = chromosome_junctions[junc_id][len(samples)+1]
                        gen2_id = chromosome_junctions[junc_id][-1]
                        if gen1_id != UNKNOWN_GENE_ID:
                            __gene_set__[row[BEDIndices.chromosome]][gen1_id].gene_reads_by_sample[s] += 1
                        else:
                            __unknown_reads_per_sample__[s] += 1
                        if gen2_id != UNKNOWN_GENE_ID and gen2_id != gen1_id:
                            __gene_set__[row[BEDIndices.chromosome]][gen2_id].gene_reads_by_sample[s] += 1
                        chromosome_junctions[junc_id][sample_index] += 1
                        row_no += 1
                    print_progress(num_rows, num_rows, '\t\t\t')
            sample_index += 1
        write_chromosome_junctions_to_file(chromosome_junctions)
        print('\t[DONE]')


def clean_files(samples, chromosomes):
    for c in chromosomes:
        for s in map(lambda sam: os.path.splitext(os.path.basename(sam))[0], samples):
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
    for id_to_gene in __gene_set__.values():
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
                + __out_delimiter__ + "gene_name (start)"
                + __out_delimiter__+"gene_id (start)"
                + __out_delimiter__ + "gene_name (end)"
                + __out_delimiter__ + "gene_id (end)\n")
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
