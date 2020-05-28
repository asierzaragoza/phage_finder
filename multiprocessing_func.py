from blast_parser import blast_parser
import subprocess, os
from Bio.SeqIO.FastaIO import SimpleFastaParser


def worker_blastn(args, item_tuple):
    result = blast_parser.run_blastn(item_tuple[0][1], item_tuple[1][1], blast_outfile=args.outdir + item_tuple[0][0] +
                                        '--' + item_tuple[1][0] + '__blastn', threads=str(args.proc), min_identity=str(95.0))
    if result is not None:
        return (item_tuple, result)
    else:
        return item_tuple


def worker_hmmsearch(args, input_tuple):
    #input_tuple: 0=output name, 1=hmm file, 2=seqdb
    p = subprocess.call([args.hmmer_path + 'hmmsearch', '--domtblout', args.outdir + input_tuple[0] + '.domtblout',
                             input_tuple[1], input_tuple[2]], stdout = open(os.devnull, 'w'))
    if p == 0: #call returns the exit, 0 is usually 'everything OK'
        return input_tuple[0]
    else:
        return None


def worker_extract_contigs(contigs_to_extract, fasta_file):
    result_fastas = []
    with open(fasta_file) as in_handle:
        for id, seq in SimpleFastaParser(in_handle):
            if id.split(' ')[0] in contigs_to_extract:
                result_fastas.append('>' + id.split(' ')[0] + '\n' + seq + '\n')
    return result_fastas
