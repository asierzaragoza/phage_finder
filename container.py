import pickle
from Bio.SeqIO.FastaIO import SimpleFastaParser


class Container():
    def __init__(self):
        self.saved_args = None
        self.input_dict = None
        self.contig_length_dict = None
        self.script_step = 0
        self.blastn_results = []
        self.seed_contigs = set([])
        self.seed_contigs_in_sets = set([])
        self.overlap_df = None
        self.sets_calculated = None
        self.large_sets_finished = False
        self.duo_sets_finished = False

    def save_pickle(self):
        with open(self.saved_args.outdir + 'container.pickle', 'wb') as out_handle:
            pickle.dump(self, out_handle)

    def build_length_dict(self):
        contig_list = self.build_list_from_input_dict('contig_file')
        self.contig_length_dict = {k:{} for k in [x[0] for x in contig_list]}
        for contig_file_tuple in contig_list:
            length_dict = {}
            with open(contig_file_tuple[1]) as in_handle:
                for id, seq in SimpleFastaParser(in_handle):
                    length_dict[id.split(' ')[0]] = len(seq)
            self.contig_length_dict[contig_file_tuple[0]] = length_dict

    def build_list_from_input_dict(self, target_key):
        final_list = []
        for key in self.input_dict.keys():
            if target_key in self.input_dict[key]:
                final_list.append((key, self.input_dict[key][target_key]))
        return final_list

    def write_overlap_tsv(self):
        self.overlap_df.to_csv(self.saved_args.outdir + 'overlap_df.tsv', index=False, sep='\t')

