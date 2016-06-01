from utils import map_to_dict, represents_int


class GTFIndices:
    def __init__(self):
        pass
    seq_name = 0
    source = 1
    feature = 2
    start = 3
    end = 4
    score = 5
    strand = 6
    frame = 7
    attribute = 8


class GeneDescription:
    def __init__(self, row, sample_list):
        attr_map = GeneAttributes(str(row[GTFIndices.attribute])).attr_dict
        self.name = "unnamed" if "gene_name" not in attr_map else attr_map["gene_name"]
        self.id = "unknown" if "gene_id" not in attr_map else attr_map["gene_id"]
        self.id = self.id.split(".")[0]
        self.start_position = row[GTFIndices.start]
        self.end_position = row[GTFIndices.end]
        self.gene_reads_by_sample = map_to_dict(lambda x:x, lambda _:0, sample_list)
        self.strand = row[GTFIndices.strand]
        self.chromosome = str(row[GTFIndices.seq_name]).replace("chr", "")
        if len(self.chromosome) == 1 and represents_int(self.chromosome):
            self.chromosome = "0" + self.chromosome

class GeneAttributes:
    def __init__(self, attribute_str):
        pair_list = map(lambda x:x.lstrip().split(" "), filter(lambda y: len(y) > 1, attribute_str.split(";")))
        self.attr_dict = map_to_dict(lambda x: x[0],
                                     lambda x: x[1].strip("\""),
                                     pair_list)
