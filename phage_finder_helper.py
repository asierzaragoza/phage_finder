import pandas as pd
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
import multiprocessing
import itertools
import functools
import math
import multiprocessing_func
import os
from blast_parser import blast_parser

#todo: lacks test_file_path
#todo: lacks save_pickle
#todo: lacks get_length_from_contig
#todo: add in test_file_path a check to see if the file is larger than 0 bytes


def write_blast_dict(Blast_Hit, blast_file):
    p_aln_size = Blast_Hit.match_length * 2 / sum(Blast_Hit.family.parent_lengths)
    merged_contig_size = sum(Blast_Hit.family.parent_lengths) - Blast_Hit.match_length

    blast_dict = {'seq_left':Blast_Hit.parents[0],
            'seq_right':Blast_Hit.parents[1], 'length_left':Blast_Hit.family.parent_lengths[0],
            'length_right': Blast_Hit.family.parent_lengths[1], 'aln_length':Blast_Hit.match_length,
            'aln_start_left':Blast_Hit.seq1pos[0], 'aln_end_left':Blast_Hit.seq1pos[1],
            'aln_start_right':Blast_Hit.seq2pos[0], 'aln_end_right':Blast_Hit.seq2pos[1],
            'p_aln_size':p_aln_size, 'merged_seq_size':merged_contig_size, 'similarity':Blast_Hit.similarity,
            'gaps':Blast_Hit.gaps, 'identity':Blast_Hit.identity, 'file_tuple':blast_file}
    return blast_dict


def test_file_path(file_path):
    #check if the path provided actually points to a file
    if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
        return True
    else:
        return False


def get_length_from_contig(id, length_dict):
    try:
        return length_dict[id]
    except KeyError:
        return None


#Main Functions
def validate_input(args, container):
    validate_input.name = 'validate_input'
    def test_separator(protein_file, separator):
        # 1) Check the separator divides all fasta headers in 2 parts
        # 2) Find out the order of the contig and protein parts in the header
        total_proteins = set([])
        total_contigs = set([])

        with open(protein_file) as in_handle:
            for id, seq in SimpleFastaParser(in_handle):
                split_id = id.split(separator)
                assert (len(split_id) == 2), 'header {} in file {} does not work with separator {}'.format(id,protein_file,separator)
                total_proteins.add(split_id[0])
                total_contigs.add(split_id[1])
        assert (len(total_proteins) != len(total_contigs)), 'Protein file {} does not look properly formatted (same number of protein IDs and contig IDs)'.format(protein_file)
        if len(total_proteins) > len(total_contigs):
            return False
        else:
            return True

    def test_contig_protein_coherency(contig_file, protein_file, separator, reverse_protein):
        #check that the contig headers for the protein file and the contig file are the same
        contig_header_set = set()
        protein_header_set = set()
        with open(contig_file) as in_handle:
            for id, seq in SimpleFastaParser(in_handle):
                contig_header_set.add(id.split(' ')[0])
        with open(protein_file) as in_handle:
            for id, seq in SimpleFastaParser(in_handle):
                splitHeader = id.rstrip('\n').split(separator)
                if reverse_protein == False:
                    protein_header_set.add(splitHeader[-1])
                else:
                    protein_header_set.add(splitHeader[0])
        if contig_header_set >= protein_header_set:
            return True
        else:
            print('\tWARNING: files {} & {} do not have the same contig headers. They will not be included.'.format(contig_file, protein_file))
            return False

    input_dict = {}
    with open(args.input) as in_handle:
        for i, line in enumerate(in_handle):
            splitLine = line.rstrip('\n').split('\t')
            if len(splitLine) != 3:
                continue
            if not test_file_path(splitLine[0]) and test_file_path(splitLine[1]):
                continue
            candidate_dict = {'contig_file': splitLine[0], 'protein_file': splitLine[1], 'separator': splitLine[2],
                              'reverse_protein': test_separator(splitLine[1], splitLine[2])}
            if not test_contig_protein_coherency(splitLine[0], splitLine[1], splitLine[2], candidate_dict['reverse_protein']):
                continue
            input_dict['dataset_' + str(i)] = candidate_dict

        print('\tinput file successfully validated, {} datasets added'.format(len(input_dict.keys())))
        container.input_dict = input_dict
        print('\tbuilding contig_length dict...')
        container.build_length_dict()
        container.script_step += 1
        return 'OK'


def run_hmmsearch(args, container):
    #Find previous runs
    previous_runs = container.build_list_from_input_dict('hmmsearch_file')
    validated_prev_runs = set()
    if len(previous_runs) > 0:
        print('\tFound {} previous searches. Validating hmmscan files...'.format(len(previous_runs)))
        for prev in previous_runs:
            if test_file_path(prev[1]):
                validated_prev_runs.add(prev[0])
        print('\tValidation finished, skipping {} hmmsearch runs'.format(len(validated_prev_runs)))

    #Add to the input list ONLY the runs not in previous_runs
    input_list = [(x, args.hmm, container.input_dict[x]['protein_file']) for x in list(container.input_dict.keys()) if x not in validated_prev_runs]

    #Prepare multiprocessing pool, run hmmsearch. imap means the results will be sent as soon as they are finished
    print('\tRunning hmmsearch with {} files...'.format(len(input_list)))
    partial_worker = functools.partial(multiprocessing_func.worker_hmmsearch, args)
    p_pool = multiprocessing.Pool(processes=args.jobs)

    for result in p_pool.imap(partial_worker, input_list):
        if result != None:
            print('\t', result, 'finished')
            container.input_dict[result]['hmmsearch_file'] = args.outdir + result + '.domtblout'
            container.save_pickle()
    p_pool.close()
    p_pool.join()

    #Check that all runs from the input_list have finished. We have to create input_list again becase we don't add previous runs
    input_list = container.build_list_from_input_dict('protein_file')
    result_list = container.build_list_from_input_dict('hmmsearch_file')

    if len(result_list) == len(input_list):
        container.script_step += 1
        return 'OK'
    else:
        print('\tNot all tasks were completed. Remaining datasets: {}'.format(set([x[0] for x in input_list]) -
                                                                              set([x[0] for x in result_list])))
        return None


def filter_domtblout(args, container):
    #Format hmmsearch result files & filter hits
    hmmsearch_list = container.build_list_from_input_dict('hmmsearch_file')
    print('\tFormatting {} hmmsearch files...'.format(len(hmmsearch_list)))
    #Format hmmsearch results
    for hmmsearch_file_tuple in hmmsearch_list:
        header = ['target_name', 'accession', 'tlen', 'query_name', 'accession', 'qlen', 'full_evalue', 'full_score',
                  'full_bias', '#', 'of', 'c-evalue', 'i-evalue', 'dom_score', 'dom_bias', 'hmm_from', 'hmm_to',
                  'ali_from', 'ali_to', 'env_from', 'env_to', 'acc', 'target_desc', 'query_aln_len', 'target_aln_len',
                  'tlen_qlen_ratio']
        out_handle = open(args.outdir + str(hmmsearch_file_tuple[0]) + '.domtblout_tformed', 'w')
        out_handle.write('\t'.join(header) + '\n')
        with open(hmmsearch_file_tuple[1], 'r') as in_handle:
            for line in in_handle:
                if line[0] == '#':
                    continue
                splitLine = [x for x in line.rstrip('\n').split(' ') if x != '']
                splitLine = splitLine[:22] + [' '.join(splitLine[22:])]
                query_alnLen = round((int(splitLine[18]) - int(splitLine[17])) / int(splitLine[5]), 2)
                splitLine.append(str(query_alnLen))
                target_alnLen = round((int(splitLine[16]) - int(splitLine[15])) / int(splitLine[2]), 2)
                splitLine.append(str(target_alnLen))
                tlen_qlen_ratio = round(int(splitLine[2]) / int(splitLine[5]), 2)
                splitLine.append(str(tlen_qlen_ratio))
                out_handle.write('\t'.join(splitLine) + '\n')
        out_handle.close()

        #Filter hmmsearch results by evalue, aln_length, difference in size between query and target
        domtblout_df = pd.read_csv(args.outdir + str(hmmsearch_file_tuple[0]) + '.domtblout_tformed', sep='\t')
        domtblout_df = domtblout_df[domtblout_df['i-evalue'] <= 0.00001]
        domtblout_df = domtblout_df[(domtblout_df['query_aln_len'] >= 0.7) & (domtblout_df['target_aln_len'] >= 0.7)]
        domtblout_df = domtblout_df[(domtblout_df['tlen_qlen_ratio'] >= 0.5) & (domtblout_df['tlen_qlen_ratio'] <= 1.5)]
        domtblout_df.drop(labels='accession.1', axis=1, inplace=True)

        #Separate protein from contig (both are mixed in target_name)
        separator = container.input_dict[str(hmmsearch_file_tuple[0])]['separator']
        if container.input_dict[str(hmmsearch_file_tuple[0])]['reverse_protein'] is False:
            domtblout_df['contig'] = domtblout_df['target_name'].apply(lambda x: x.split(separator)[-1])
            domtblout_df['protein'] = domtblout_df['target_name'].apply(lambda x: x.split(separator)[0])
        else:
            domtblout_df['contig'] = domtblout_df['target_name'].apply(lambda x: x.split(separator)[0])
            domtblout_df['protein'] = domtblout_df['target_name'].apply(lambda x: x.split(separator)[-1])

        #Keep only 1 hit per protein
        new_domtblout = pd.DataFrame(columns=domtblout_df.columns.values)
        for id, group in domtblout_df.groupby('target_name'):
            group.sort_values(by='full_score', ascending=False, inplace=True)
            new_domtblout = new_domtblout.append(group.iloc[0, :])

        new_domtblout.to_csv(args.outdir + str(hmmsearch_file_tuple[0]) + '.domtblout_filtered', sep='\t', index=False)
        container.input_dict[hmmsearch_file_tuple[0]]['domtblout_dfs'] = args.outdir + str(hmmsearch_file_tuple[0]) + '.domtblout_filtered'

    #Check all files have been processed
    result_list = container.build_list_from_input_dict('domtblout_dfs')
    if len(hmmsearch_list) == len(result_list):
        container.script_step += 1
        return 'OK'
    else:
        print('\tNot all tasks were completed. Remaining datasets: {}'.format(set([x[0] for x in hmmsearch_list]) -
                                                                              set([x[0] for x in result_list])))
        return None


def aggregate_contig_data(args, container):
    '''Select the contigs to use in later steps.'''
    domtblout_list = container.build_list_from_input_dict('domtblout_dfs')
    for domtblout_file in domtblout_list:
        #Aggregate all hit data by contig
        contig_list = []
        domtblout_df = pd.read_csv(domtblout_file[1], sep='\t')
        for id, group_df in domtblout_df.groupby('contig'):
            contig_list.append({'contig_name': id, 'n_of_prots': len(group_df.index),
                                'length': get_length_from_contig(id, container.contig_length_dict[domtblout_file[0]])})
        new_df = pd.DataFrame(contig_list)

        #new_df has the dataframe, with columns 'contig_name', 'n_of_prots' & 'length'. We can add some filtering here
        #Remember that when working with pandas the truth-values used by python are considered untrustworthy, so you have to use
        #the numpy equivalents whenever possible, hence the use of & and np.maximum here
        new_df = new_df[(new_df['length'] >= 7500) & (new_df['length'] <= 500000) & (new_df['n_of_prots'] >= np.maximum(2,new_df['length']/20000))]

        #Add candidate contigs to seed_contigs
        container.seed_contigs = container.seed_contigs|set(list(zip(itertools.repeat(domtblout_file[0]), new_df['contig_name'].tolist())))
        new_df.to_csv(args.outdir + domtblout_file[0] + '.domtblout_contig_melted', sep='\t', index=False)
        container.input_dict[domtblout_file[0]]['contig_melted_file'] = args.outdir + domtblout_file[0] + '.domtblout_contig_melted'

    result_list = container.build_list_from_input_dict('contig_melted_file')
    if len(domtblout_list) == len(result_list):
        container.script_step += 1
        return 'OK'
    else:
        print('\tNot all tasks were completed. Remaining datasets: {}'.format(set([x[0] for x in domtblout_list]) -
                                                                              set([x[0] for x in result_list])))
        return None


def extract_seed_contigs(args, container):
    container.seed_contigs = list(container.seed_contigs)
    container.seed_contigs.sort(key=lambda x:x[0])
    container.seed_contigs = tuple(container.seed_contigs)
    out_handle = open(args.outdir + 'seed_contigs.fasta', 'w')
    for dataset, seed_contigs in itertools.groupby(container.seed_contigs, lambda x: x[0]):
        contigs_to_extract = set([x[1] for x in seed_contigs])
        print('\t',dataset, len(contigs_to_extract))
        with open(container.input_dict[dataset]['contig_file']) as in_handle:
            for id, seq in SimpleFastaParser(in_handle):
                if id in contigs_to_extract:
                    out_handle.write('>' + id + '\n' + seq + '\n')
    out_handle.close()
    container.seed_contigs = set(container.seed_contigs)
    container.script_step += 1
    return 'OK'


def run_blastn_combinations(args, container):
    '''Run blastn between all combinations of datasets.'''
    # Find previous runs
    print('starting run_blast_combinations')
    previous_runs = container.blastn_results
    validated_prev_runs = []
    if len(previous_runs) > 0:
        print('\tFound {} previous runs. Validating blastn files...'.format(len(previous_runs)))
        for prev in previous_runs:
            if test_file_path(prev[1]):
                # rstrip removes instances of the characters given, not of the complete string
                # e.g. if I do 'dog.fna__blastn'.rstrip('__blastn), the result will be 'dog.f' bc 'na' is in '__blastn'

                validated_prev_runs.append(set(prev[0]))
                # if prev not in container.blastn_results: container.blastn_results.append(tuple([run_tuple, prev]))
        print('\tValidation finished, found {} finished blastn runs'.format(len(validated_prev_runs)))

    # Get all combinations needeed...
    contig_file_list = container.build_list_from_input_dict('contig_file')
    combination_no = math.factorial(len(contig_file_list) + 1) / (
                math.factorial(2) * math.factorial(len(contig_file_list) - 1))
    combination_set = set(itertools.combinations_with_replacement(contig_file_list, 2))
    print(combination_set)

    # Then remove the combinations run previously. These are combinations (the order does not matter), so build sets
    # with the datasets and remove those that match with the ones from previous_runs
    items_to_remove = set([])
    for item in combination_set:
        #print(item)
        item_set = set([x[0] for x in item])
        #print(item_set)
        #print(validated_prev_runs)
        if item_set in validated_prev_runs:
            items_to_remove.add(item)
    print('\tcombinations removed after accounting for existing blastn files: {}'.format(len(items_to_remove)))
    combination_set -= items_to_remove
    print('\tcombinations left: {}'.format(len(combination_set)))

    partial_worker = functools.partial(multiprocessing_func.worker_blastn, args)

    p_pool = multiprocessing.Pool(processes=args.jobs)
    for result in p_pool.imap(partial_worker, combination_set):
        print('\t', result)
        if len(result) == 2:
            container.blastn_results.append(result)
            container.save_pickle()
    p_pool.close()
    p_pool.join()
    if len(container.blastn_results) == combination_no:
        container.script_step += 1
        return 'OK'
    else:
        print('\tERROR: Some blastn runs could not be completed. Please run the script again.')
        return None


def find_overlapped_contigs(args, container):
    '''Find overlap pairs between the blastn results'''

    def check_side(blast):
        '''We are only interested in overlapping blasts between sequences (so one contig is the continuation of the other).
           '''
        assert (blast.family.parent_lengths != (-1, -1))
        seq1_mean = int((blast.seq1pos[0] + blast.seq1pos[1]) / 2)
        seq2_mean = int((blast.seq2pos[0] + blast.seq2pos[1]) / 2)
        seq1_parent_len = blast.family.parent_lengths[0]
        seq2_parent_len = blast.family.parent_lengths[1]
        if (seq1_parent_len - seq1_mean) > seq1_mean and (seq2_parent_len - seq2_mean) < seq2_mean:
            # seq1 is Start, seq2 is end
            if 1 in blast.seq1pos and seq2_parent_len in blast.seq2pos:
                return True
        elif (seq1_parent_len - seq1_mean) < seq1_mean and (seq2_parent_len - seq2_mean) > seq2_mean:
            # seq1 is End, seq2 is Start
            if 1 in blast.seq2pos and seq1_parent_len in blast.seq1pos:
                return True
        else:
            return False

    blast_dict_list = []
    for blastn_result in container.blastn_results:
        blast_families = blast_parser.parallel_parse_blast_file(blastn_result[-1], min_identity=97.0, min_aln_len=250, nproc=args.jobs)

        for i, family in enumerate(blast_families):
            family.remove_own_hits()
            family.equalize()
            if family.is_duplicate():  # remove duplicates: if the family is duplicated we just move to the next family
                continue
            family.merge_blast_list(50)
            family.blast_list.sort(key=lambda x: x.match_length, reverse=True)
            for blast in family.blast_list:
                if blast.gaps > 10:  # Remove matches with > 10 gaps
                    continue

                if check_side(blast):  # check if it is a side blast, if it is store it in the dict
                    blast_dict_list.append(write_blast_dict(blast, blastn_result[-1]))
                    break

    container.overlap_df = pd.DataFrame(blast_dict_list)
    #container.write_overlap_tsv()
    container.script_step += 1
    return 'OK'


def cluster_overlapped_contigs(args, container):
    '''Group overlapping pairs into groups, then extract the groups into proper '''
    result_list = []
    #Calculate overlap groups
    if container.sets_calculated is not None:
        print('\tOverlapping sets found in \'run.pickle\', skipping step')
        result_list = container.sets_calculated
    else:
        print('\tCalculating overlapping sets...')
        overlap_contig_pairs = container.overlap_df.apply(lambda x: (x['seq_left'], x['seq_right']), axis=1).tolist()
        overlap_dataset = set([frozenset(x) for x in overlap_contig_pairs])

        while (overlap_dataset):
            nset = set(overlap_dataset.pop())
            check = len(overlap_dataset)
            while check:
                check = False
                for iset in overlap_dataset.copy():
                    if nset.intersection(iset):
                        check = True
                        overlap_dataset.remove(iset)
                        nset.update(iset)
            result_list.append(tuple(nset))
        container.sets_calculated = result_list
        container.save_pickle()

    #get a set of the seed_contigs to compare with the sets:
    seed_contig_set = set([x[1] for x in container.seed_contigs])
    print('number of seed_contigs: {}'.format(len(seed_contig_set)))
    #Process large sets
    total_large_sets = 0
    stored_large_sets = 0
    if container.large_sets_finished == True:
        print('\t\'large_sets_finished\' flag is TRUE, skipping')
    else:
        print('\tProcessing large sets')
        if not os.path.isdir(args.outdir + 'large_sets'): os.mkdir(args.outdir + 'large_sets')
        print('candidate large_sets: {}'.format(len([x for x in result_list if len(x) > 2])))
        for i, large_set in enumerate([x for x in result_list if len(x) > 2]):
            total_large_sets += 1
            if len(seed_contig_set.intersection(large_set)) > 0:
                stored_large_sets += 1
                for contig in seed_contig_set.intersection(large_set):
                    container.seed_contigs_in_sets.add(contig)
                #Extract the datasets that contain the contigs present in this set
                large_set_df = container.overlap_df[container.overlap_df['seq_left'].isin(large_set) &
                                                    container.overlap_df['seq_right'].isin(large_set)]
                blastn_file_set = set(large_set_df['file_tuple'].tolist())
                blastn_file_set = set(list(itertools.chain.from_iterable([[x.split('/')[-1].split('--')[0],
                                                                           x.split('/')[-1].split('--')[-1].split('__')[0]] for x in blastn_file_set])))

                fasta_files_to_extract = [container.input_dict[x]['contig_file'] for x in blastn_file_set]


                p_pool = multiprocessing.Pool(processes=args.jobs)
                iterable = zip(itertools.repeat(large_set), fasta_files_to_extract)
                extracted_contigs = p_pool.starmap(multiprocessing_func.worker_extract_contigs, iterable)
                extracted_contigs = list(set(itertools.chain.from_iterable(extracted_contigs)))

                if len(extracted_contigs) != len(large_set):
                    print('Did not find all contigs.')
                    print('extracted contigs:', len(extracted_contigs))
                    print('large set:', len(large_set))
                    exit()
                else:
                    with open(args.outdir + 'large_sets/'  + 'large_set_' + str(i) + '.fna', 'w') as out_handle:
                        out_handle.write(''.join(extracted_contigs))
                    large_set_df.to_csv(args.outdir + 'large_sets/' + 'large_set_' + str(i) + '.tsv', sep = '\t', index=False)
                p_pool.terminate()
        container.large_sets_finished = True
        print('Out of {} sets, {} sets were saved'.format(total_large_sets, stored_large_sets))
        container.save_pickle()

    if container.duo_sets_finished == True:
        print('\t\'duo_sets_finished\' flag is TRUE, skipping')
    else:
        print('\tProcessing duo sets')
        if not os.path.isdir(args.outdir + 'duo_sets'): os.mkdir(args.outdir + 'duo_sets')
        print('candidate dup_sets: {}'.format(len([x for x in result_list if len(x) == 2])))
        #Process duo sets
        total_duo_sets = 0
        stored_duo_sets = 0
        for i, duo_set in enumerate([x for x in result_list if len(x) == 2]):
            total_duo_sets += 1
            if len(seed_contig_set.intersection(duo_set)) > 0:
                stored_duo_sets += 1
                for contig in seed_contig_set.intersection(duo_set):
                    container.seed_contigs_in_sets.add(contig)
                duo_set_df = container.overlap_df[container.overlap_df['seq_left'].isin(duo_set) &
                                                    container.overlap_df['seq_right'].isin(duo_set)]

                blastn_file_set = set(duo_set_df['file_tuple'].tolist())
                blastn_file_set = set(list(itertools.chain.from_iterable([[x.split('/')[-1].split('--')[0],
                                                                           x.split('/')[-1].split('--')[-1].split('__')[
                                                                               0]] for x in blastn_file_set])))

                fasta_files_to_extract = [container.input_dict[x]['contig_file'] for x in blastn_file_set]

                p_pool = multiprocessing.Pool(processes=args.jobs)
                iterable = zip(itertools.repeat(duo_set), fasta_files_to_extract)
                extracted_contigs = p_pool.starmap(multiprocessing_func.worker_extract_contigs, iterable)
                extracted_contigs = list(set(itertools.chain.from_iterable(extracted_contigs)))

                if len(extracted_contigs) != len(duo_set):
                    print('Did not find all contigs.')
                    print('extracted contigs:', len(extracted_contigs))
                    print('duo set:', len(duo_set))
                else:
                    with open(args.outdir + 'duo_sets/' + 'duo_set_' + str(i) + '.fna', 'w') as out_handle:
                        out_handle.write(''.join(extracted_contigs))
                    duo_set_df.to_csv(args.outdir + 'duo_sets/' + 'duo_set_' + str(i) + '.tsv', sep = '\t', index=False)
                p_pool.terminate()
        container.duo_sets_finished = True
        print('Out of {} sets, {} sets were saved'.format(total_duo_sets, stored_duo_sets))
        container.save_pickle()
    container.write_overlap_tsv()
    container.script_step += 1
    return 'OK'


def extract_single_contigs(args, container):
    out_handle = open(args.outdir + 'seed_contigs.fna', 'w')
    single_contig_list = [x for x in container.seed_contigs if x[1] not in container.seed_contigs_in_sets]
    single_contig_list.sort(key=lambda x: x[1])
    out_handle = open(args.outdir + 'single_seed_contigs.fna', 'w')
    for dataset, contigs_to_extract in itertools.groupby(single_contig_list, lambda x:x[0]):
        contigs_to_extract = set([x[1] for x in contigs_to_extract])
        with open(container.input_dict[dataset]['contig_file']) as in_handle:
            for id, seq in SimpleFastaParser(in_handle):
                if id in contigs_to_extract:
                    out_handle.write('>' + id + '\n' + seq +'\n')
    out_handle.close()
    container.script_step += 1
    return 'OK'
