# Junctions Comparer
A simple Python tool to compare junctions across several samples

#Usage 

```
usage: junctions_comparer.py [-h] [-nb | -q] [-n MIN_READS] [-t TEMP_DIR]
                             [-o OUT_DIR] [-r RESULTS_FILE]
                             [-rl RESULTS_DELIMITER] [-u UNKNOWN_ID]
                             [-gl GENE_DELIMITER] [-a {stranded,unstranded}]
                             gtf_file bed_file bed_file [bed_file ...]
```

# Arguments
```
positional arguments:
  gtf_file              reference annotation file in GTF format
  bed_file              junction annotation file in BED format
```

# Options
```
optional arguments:
  -h, --help            show this help message and exit
  -nb, --no-bars        disables loading bars (use when redirecting output to
                        file)
  -q, --quiet
  -n MIN_READS, --min-reads MIN_READS
                        minimum number of reads supporting junction to pass
                        filter (default: 3)
  -t TEMP_DIR, --temp-dir TEMP_DIR
                        directory where to keep temporary files (default:
                        sj_tmp)
  -o OUT_DIR, --out-dir OUT_DIR
                        directory where to keep output files (default: sj_out)
  -r RESULTS_FILE, --results-file RESULTS_FILE
                        name for final result file (default: results.csv)
  -rl RESULTS_DELIMITER, --results-delimiter RESULTS_DELIMITER
                        delimiter for output CSV file (default: ,)
  -u UNKNOWN_ID, --unknown-id UNKNOWN_ID
                        name for unknown genes (default: UNKNOWN)
  -gl GENE_DELIMITER, --gene-delimiter GENE_DELIMITER
                        delimiter for gene lists in output (default: |)NOTE:
                        MUST NOT BE THE SAME AS -rl
  -a {stranded,unstranded}, --analysis {stranded,unstranded}
                        type of analysis (default: unstranded)
  ```
