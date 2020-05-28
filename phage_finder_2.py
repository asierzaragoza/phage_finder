#!/usr/bin/python3
import argparse, os, pickle
import blast_parser.blast_parser as blast_parser
from container import Container
import phage_finder_helper

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, help='tab-separated file')
parser.add_argument('-m', '--hmm', default = './phage_clusters_filtered.reduced.hmm', help='hmm file to search')
parser.add_argument('-j', '--jobs', default=1, type=int, help='number of simultaneous processes (default=1)')
parser.add_argument('-n', '--proc', default=1, type=int, help='number of processes per run (default=1)')
parser.add_argument('-o', '--outdir', required=True, help='output directory')
parser.add_argument('--blast_path', default='', dest='blast_path', help='location of blast binaries (default = $PATH')
parser.add_argument('--hmmer_path', default='', dest='hmmer_path', help='location of HMMER binaries (default = $PATH')
parser.add_argument('--continue', dest='_continue', action='store_true', default=False, help='continue previous run?')

args = parser.parse_args()
args.outdir = args.outdir.rstrip('/') + '/'



#Check blast binaries location
if args.blast_path != '':
    args.blast_path = args.blast_path.rstrip('/') + '/'
else:
    print('No --blast_path option selected, using whatever is on $PATH')
if not blast_parser.check_blast_path(args.blast_path):
    exit()

#Check hmmer binaries location
if args.hmmer_path != '':
    args.hmmer_path = args.hmmer_path.rstrip('/') + '/'
else:
    print('No --hmmer_path option selected, using whatever is on $PATH')

#Create outdir if it does not exist
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

#Create a new container object or load one if there is one available
container_object = None
if args._continue == True and os.path.isfile(args.outdir + 'container.pickle'):
    with open(args.outdir + 'container.pickle', 'rb') as in_handle:
        container_object = pickle.load(in_handle)
        args = container_object.saved_args
elif args._continue == True:
    print('\'container.pickle\' cannot be found in {}, but --continue was selected'.format(args.outdir))
    exit()
else:
    container_object = Container()
    container_object.saved_args = args

#Main Loop
script_step = 0
function_list = [phage_finder_helper.validate_input, phage_finder_helper.run_hmmsearch,
                 phage_finder_helper.filter_domtblout, phage_finder_helper.aggregate_contig_data,
                 phage_finder_helper.run_blastn_combinations, phage_finder_helper.find_overlapped_contigs,
                 phage_finder_helper.cluster_overlapped_contigs, phage_finder_helper.extract_single_contigs]


for i in range(container_object.script_step, len(function_list)):
    print('Running step {} / {}: {}'.format(i, len(function_list), function_list[i].__name__))
    return_code = function_list[i](args, container_object)
    if return_code != 'OK':
        exit()
    container_object.save_pickle()
print('finished!')