#Initialize variables
import subprocess, platform, os, resource, itertools, multiprocessing


class Blast_Hit():
    def __init__(self):
        self.parents = (None, None)
        self.parent_lengths = (0, 0)
        self.seq1pos = (-1, -1)
        self.seq2pos = (-1, -1)

        self.mismatches = None
        self.gaps = None
        self.identity = None
        self.similarity = None
        self.match_length = None
        self.evalue = None
        self.inverted = None
        self.family = None

    def add_family(self, blast_family):
        assert set(self.parents) == set(blast_family.parents)
        self.family = blast_family

    def reverse_blast(self):
        if self.inverted == True:
            if self.seq1pos[0] > self.seq1pos[1]:
                self.seq1pos = (self.seq1pos[1], self.seq1pos[0])
            elif self.seq2pos[0] > self.seq2pos[1]:
                self.seq2pos = (self.seq2pos[1], self.seq2pos[0])
            self.inverted = False
        else:
            self.seq1pos = (self.seq1pos[1], self.seq1pos[0])
            self.inverted = True
        return self

    def build_from_blast(self, line):
        '''qseqid sseqid length qlen slen pident positive mismatch gaps qstart qend sstart send evalue bitscore'''
        blast_line = line.rstrip('\n').split('\t')
        self.parents = (blast_line[0], blast_line[1])
        self.parent_lengths = (int(blast_line[3]), int(blast_line[4]))
        self.match_length = int(blast_line[2])
        self.seq1pos = (int(blast_line[9]), int(blast_line[10]))
        self.seq2pos = (int(blast_line[11]), int(blast_line[12]))
        self.identity = float(blast_line[5])
        self.similarity = float(blast_line[6])
        self.mismatches = int(blast_line[7])
        self.gaps = int(blast_line[8])
        self.bitscore = float(blast_line[14])

        #Process evalue
        if 'e' in blast_line[13]:
            evalue_split = blast_line[13].split('e')
            if evalue_split[1][0] == '+':
                self.evalue = float(evalue_split[0]) * (10 ^ int(evalue_split[1][1:]))
            elif evalue_split[1][0] == '-':
                self.evalue = float(evalue_split[0]) * (10 ^ int(evalue_split[1][1:]) * -1)
            else:
                print('I do not understand this E-Value: {}'.format(blast_line[13]))
                raise ValueError
        else:
            self.evalue = float(blast_line[13])

        #Process inverted
        if self.seq1pos[0] > self.seq1pos[1]:
            self.inverted = True
        else:
            self.inverted = False

        return self

class Blast_Family():
    def __init__(self, type, filter_dict, parent_list, parent_lengths):
        self.parents = parent_list
        self.parent_lengths = parent_lengths
        assert type in ['blastn', 'tblastx']
        self.blast_list = []
        self.filter_dict = filter_dict

    def __iter__(self):
        return self.blast_list

    def __contains__(self, item):
        return item in self.blast_list

    def add_blast(self, blast_hit):
        assert set(self.parents) == set(blast_hit.parents)
        self.blast_list.append(blast_hit)

    def sort(self, by='seq1pos'):
        assert by in ['seq1pos', 'match_len']
        if by == 'seq1pos':
            self.blast_list.sort(key= lambda Blast_Hit : Blast_Hit.seq1pos)
        elif sortBy == 'matchLen':
            self.blast_list.sort(key = lambda Blast_Hit: Blast_Hit.match_length)

    def equalize(self):
        for Blast_Hit in self.blast_list:
            if Blast_Hit.parents[0] != self.parents[0]:
                newSeq2 = Blast_Hit.seq1pos
                newSeq1 = Blast_Hit.seq2pos
                Blast_Hit.parents = self.parents
                Blast_Hit.seq1pos = newSeq1
                Blast_Hit.seq2pos = newSeq2

    def is_duplicate(self, identity_threshold = 99.5, size_threshold = 0.01):
        #print('{} / {}: checking duplicates... '.format(self.parents[0], self.parents[1]), end = '')
        '''Check if the Blast_Family is between 2 identical sequences. For that, we get all blasts over a specified
        similarity_threshold, then sum the length of all blasts over that size and check if it is close or similar to
        the length of any of the two parent sequences.
        #todo: Might have problems with multiple blasts to the same zone (e.g. ribosomal proteins). Should check later'''
        assert -1 not in self.parent_lengths

        total_blast_len = sum([x.match_length for x in self.blast_list if x.identity >= identity_threshold])

        '''Old code does something weird here (looks like either I didn't bother to check if the OR operator was exclusive
        in python, or this was some strange hack to evade some kind of error... If the function starts doing funny stuff,
         start looking here'''
        seq_1_equal = (1-size_threshold) <= total_blast_len / self.parent_lengths[0] <= (1+size_threshold)
        seq_2_equal = (1-size_threshold) <= total_blast_len / self.parent_lengths[1] <= (1+size_threshold)

        if seq_1_equal or seq_2_equal:
            #print('equal = TRUE;\t{}\t{}\t{}'.format(total_blast_len, self.parent_lengths[0], self.parent_lengths[1]))
            return True
        else:
            #print('equal = FALSE;\t{}\t{}\t{}'.format(total_blast_len, self.parent_lengths[0], self.parent_lengths[1]))
            return False

    def merge_blast_list(self, threshold):
        def merge_blasts(blast_to_merge_list):
            pident = round(sum([x.identity for x in blast_to_merge_list]) / len(blast_to_merge_list), 2)
            ppos = round(sum([x.similarity for x in blast_to_merge_list]) / len(blast_to_merge_list), 2)
            bitscore = round(sum([x.bitscore for x in blast_to_merge_list] )/ len(blast_to_merge_list), 1)
            evalue = sum([x.evalue for x in blast_to_merge_list]) / len(blast_to_merge_list)
            gaps = sum([x.gaps for x in blast_to_merge_list])
            mismatches = sum([x.mismatches for x in blast_to_merge_list])
            qstart = blast_to_merge_list[0].seq1pos[0]
            qend = blast_to_merge_list[-1].seq1pos[1]
            sstart = blast_to_merge_list[0].seq2pos[0]
            send = blast_to_merge_list[-1].seq2pos[1]
            match_length = qend - qstart
            new_blast_line = '\t'.join(list(map(str, [blast_to_merge_list[0].family.parents[0], blast_to_merge_list[0].family.parents[1], match_length,
                                            blast_to_merge_list[0].family.parent_lengths[0], blast_to_merge_list[0].family.parent_lengths[1],
                                            pident, ppos, mismatches, gaps, qstart, qend, sstart, send, evalue, bitscore])))

            return Blast_Hit().build_from_blast(new_blast_line)


        def find_merge_candidates(blast_list, threshold):
            merge_candidate_list = []
            in_merged_list = []
            for i in range(0, len(blast_list) -1):
                fst_blast = blast_list[i]
                snd_blast = blast_list[i + 1]
                pos_1_dtce = snd_blast.seq1pos[0] - (fst_blast.seq1pos[1] + 0.01)
                pos_2_dtce = snd_blast.seq2pos[0] - (fst_blast.seq2pos[1] + 0.01)

                #first check if the blast hits overlap. Id they do, add them to the candidate merge list
                if pos_1_dtce < 0 and pos_2_dtce < 0:
                    merge_candidate_list.append([fst_blast, snd_blast])
                    in_merged_list.append(fst_blast)
                    in_merged_list.append(snd_blast)
                    continue
                #if they don't, check that the dtce between blasts is between the specified parameters
                elif pos_1_dtce <= threshold and pos_2_dtce <= threshold:
                    merge_candidate_list.append([fst_blast, snd_blast])
                    in_merged_list.append(fst_blast)
                    in_merged_list.append(snd_blast)

            #remove merged candidates from the blast_list:
            final_merged_list = list(set(blast_list) - set(in_merged_list))
            # Merge concatenated pairs (merge pairs that have a blast hit in common)
            i = 0
            while i < len(merge_candidate_list) - 1:
                if merge_candidate_list[i][-1] == merge_candidate_list[i + 1][0]:
                    new_merge_candidate_list = merge_candidate_list[i][:-1] + merge_candidate_list[i + 1]
                    merge_candidate_list[i] = new_merge_candidate_list
                    merge_candidate_list.pop(i + 1)
                    i = 0
                    continue
                else:
                    i += 1

            #Now merge the candidate blasts together

            for candidate_list in merge_candidate_list:
                final_merged_list.append(merge_blasts(candidate_list))

            return final_merged_list

        self.equalize()
        std_blast_list = [x for x in self.blast_list if x.inverted == False]
        std_blast_list.sort(key=lambda Blast_Hit: Blast_Hit.seq1pos)
        std_blast_list = find_merge_candidates(std_blast_list, threshold)
        for blast_hit in std_blast_list:
            blast_hit.add_family(self)
            #print(blast_hit.family.parent_lengths)

        rev_blast_list = [x.reverse_blast() for x in self.blast_list if x.inverted == True]
        rev_blast_list = find_merge_candidates(rev_blast_list, threshold)
        for blast_hit in rev_blast_list:
            blast_hit.add_family(self)
            #print(blast_hit.family.parent_lengths)

        self.blast_list = std_blast_list + rev_blast_list

    def remove_own_hits(self):
        '''Iterate over the blast list and remove all BlastHits that are equal. We define equal as having the same positions
        for seq1 and seq2.'''

        clean_blast_list = []
        frozenset_list = []
        own_hits = 0
        for Blast_Hit in self.blast_list:
            for blast_hit_frozenset in frozenset_list:
                if frozenset(Blast_Hit.seq1pos + Blast_Hit.seq2pos) == blast_hit_frozenset:
                    own_hits += 1
                    break
            else:
                clean_blast_list.append(Blast_Hit)
                frozenset_list.append(frozenset(Blast_Hit.seq1pos + Blast_Hit.seq2pos))

        #print('{} duplicate hits removed, family now contains {} hits'.format(own_hits, len(clean_blast_list)))
        self.blast_list = clean_blast_list

def run_tblastx(fasta_file_search, fasta_file_db=None, blast_matrix=None, blast_outfile = 'blast_seqs.temp.blast', blast_path = '', threads=4):
    print('running makeblastdb... ')
    if not fasta_file_db:
        makeblastdb_proc = subprocess.Popen(
            [blast_path + 'makeblastdb', '-in', fasta_file_search, '-out', 'dbTemp', '-dbtype', 'nucl'])
        makeblastdb_proc.wait()
    else:
        makeblastdb_proc = subprocess.Popen(
            [blast_path + 'makeblastdb', '-in', fasta_file_db, '-out', 'dbTemp', '-dbtype', 'nucl'])
        makeblastdb_proc.wait()
    print('done!')

    # Run blast
    print('running tblastx...')
    blast_process = None
    executable_name = None

    if platform.system() == 'Windows': executable_name = 'tblastx.exe'
    else: executable_name = 'tblastx'

    blast_process = subprocess.Popen([blast_path + executable_name, '-matrix', blast_matrix, '-query', fasta_file, '-db',
                                    'dbTemp', '-out', blast_outfile, '-num_threads', threads, '-outfmt', '6 qseqid sseqid length qlen '
                                    'slen pident ppos mismatch gaps qstart qend sstart send evalue bitscore'])

    # There is no single way to limit memory usage directly in both unix and windows, so I'll need to write a watchdog
    # function to constantly check memory usage and kill the program if it gets out of hand -> with multiprocessing?
    #print(system_utilities.find_process_by_pid(blast_process.pid))
    blast_process.wait()
    print('done!')

    os.remove(os.path.join(os.getcwd() + os.sep,'dbTemp.nhr'))
    os.remove(os.path.join(os.getcwd() + os.sep,'dbTemp.nin'))
    os.remove(os.path.join(os.getcwd() + os.sep,'dbTemp.nsq'))
    return (blast_outfile)

def run_blastn(fasta_file_search, fasta_file_db=None, blast_outfile = 'blast_seqs.temp.blast', blast_path = '', threads = '4', min_identity = '0.0'):
    #if fasta_file_db = None, do a self blast
    if not fasta_file_db:
        fasta_file_db = fasta_file_search
    devnull = open(os.devnull, 'wb')
    #run makeblastdb
    makeblastdb_proc = subprocess.Popen([blast_path + 'makeblastdb', '-in', fasta_file_db, '-out', 'dbTemp_' +
                        str(os.getpid()), '-dbtype', 'nucl'], stdout=devnull, stderr = subprocess.PIPE)
    error = makeblastdb_proc.stderr.read()
    makeblastdb_proc.communicate()
    if makeblastdb_proc.returncode != 0:
        print('makeblastdb failed. The error message is shown below:\n{}'.format(error.decode('UTF-8')))
        return None
    devnull.close()

    # Run blastn
    blast_process = None
    executable_name = None
    if platform.system() == 'Windows':
        executable_name = 'blastn.exe'
    else:
        executable_name = 'blastn'

    blast_process = subprocess.Popen([blast_path + executable_name, '-query', fasta_file_search, '-db', 'dbTemp_' + str(os.getpid()),
                                    '-out', blast_outfile, '-num_threads', threads, '-perc_identity', min_identity, '-outfmt', '6 qseqid sseqid length qlen '
                                    'slen pident ppos mismatch gaps qstart qend sstart send evalue bitscore'], stderr = subprocess.PIPE)
    error = blast_process.stderr.read()
    #print(system_utilities.find_process_by_pid(blast_process.pid))
    blast_process.communicate()
    if blast_process.returncode != 0:
        print('blastn failed. The error message is shown below:\n{}'.format(error.decode('UTF-8')))
        for file in [x for x in os.listdir(os.getcwd()) if x.split('.')[0] == 'dbTemp_' + str(os.getpid())]:
            os.remove(os.path.join(os.getcwd() + os.sep, file))
        return None
    else:
        print('blastn finished successfully!')
        for file in [x for x in os.listdir(os.getcwd()) if x.split('.')[0] == 'dbTemp_' + str(os.getpid())]:
            os.remove(os.path.join(os.getcwd() + os.sep, file))
        return (blast_outfile)

def check_blast_path(blast_path=''):
    print('Checking blast_path... ', end='')
    p=None
    #Check that the executable is found in the path provided
    if blast_path != '':
        p = subprocess.Popen([os.path.join(blast_path + os.sep, 'blastn'), '-version'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    else: #if no blast_path is given, check that it is in the path
        p = subprocess.Popen(['blastn', '-version'], stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    output = p.stdout.read()
    error = p.stderr.read()
    p.communicate()
    if p.returncode != 0 and blast_path != '':
        print('blast_path test failed (exit code =/=0). Please check that the blast_path provided ({}) is correct.\n '
              'The output of stderror is shown below:\n{}'.format(blast_path,error.decode('UTF-8')))
        return False
    elif p.returncode != 0:
        print('Unable to find blast in the PATH.')
        return False

    #check blastn version
    split_output = output.decode('UTF-8').split('\n')[0].split(' ')
    version_numbers = list(map(int, split_output[1][:-1].split('.')))
    if version_numbers[0] >= 2 and version_numbers[1] >= 6:
        print('BLAST+ version is {}'.format(split_output[1]))
        return True
    else:
        print('BLAST+ version is outdated ({}). Please update your version of BLAST+.'.format(split_output[1]))
        return False

#Basic filtering: self hits, min length, min identity + group blasthits into families
def parse_blast_file(blast_file, type = 'blastn', min_identity = 0, min_ppos = 0,  min_aln_len = 0):
    cause_dict = {'Self hits': 0, 'Low identity': 0, 'Low similarity':0, 'Small match': 0}
    filter_dict = {'min_identity':min_identity, 'min_similarity':min_ppos, 'min_aln_len': min_aln_len}
    total_blast_hits = 0
    accepted_blast_hits = 0
    blast_family_list = []
    blast_hit_list = []

    with open(blast_file, 'r') as in_handle:
        for line in in_handle:
            if len(line.split('\t')) != 15:
                continue
            else:
                total_blast_hits += 1

                new_blast_hit = Blast_Hit().build_from_blast(line)
                # Remove self-hits
                if new_blast_hit.parents[0] == new_blast_hit.parents[1]:
                    cause_dict['Self hits'] += 1
                    continue
                # Remove low identity hits
                elif new_blast_hit.identity < min_identity:
                    cause_dict['Low identity'] += 1
                    continue
                elif new_blast_hit.similarity < min_ppos:
                    cause_dict['Low similarity'] += 1
                    continue
                # Remove small hits
                elif new_blast_hit.match_length < min_aln_len:
                    cause_dict['Small match'] += 1
                    continue
                else:
                    accepted_blast_hits += 1
                    blast_hit_list.append(tuple([set(new_blast_hit.parents), new_blast_hit]))
    blast_hit_list.sort(key=lambda x: x[0])

    for key, group in itertools.groupby(blast_hit_list, lambda x: x[0]):
        new_blast_family = Blast_Family(type, filter_dict, key, (0, 0))
        for blast_hit in [x[1] for x in list(group)]:
            new_blast_family.add_blast(blast_hit)
            blast_hit.add_family(new_blast_family)
        new_blast_family.parents = new_blast_family.blast_list[0].parents
        new_blast_family.parent_lengths = new_blast_family.blast_list[0].parent_lengths
        new_blast_family.equalize()
        blast_family_list.append(new_blast_family)

    print('{} blast hits accepted out of {}, grouped into {} families'.format(accepted_blast_hits,
                                                    total_blast_hits, len(blast_family_list)))
    print('Self hits: {}\tLow identity (<{}): {}\tLow similarity (<{}): {}\tSmall matches(<{}): {}'.format(
        cause_dict['Self hits'], min_identity, cause_dict['Low identity'], min_ppos, cause_dict['Low similarity'],
        min_aln_len, cause_dict['Small match']))
    return blast_family_list

def parse_blast_file_worker(line_chunk, filter_dict):
    return_list = []
    for line in line_chunk:
        if len(line.split('\t')) != 15:
            continue
        else:
            new_blast_hit = Blast_Hit().build_from_blast(line)
            if new_blast_hit.parents[0] == new_blast_hit.parents[1]:
                continue
            elif new_blast_hit.identity < filter_dict['min_identity']:
                continue
            elif new_blast_hit.similarity < filter_dict['min_similarity']:
                continue
            elif new_blast_hit.match_length < filter_dict['min_aln_len']:
                continue
            else:
                return_list.append(tuple([set(new_blast_hit.parents), new_blast_hit]))
    return return_list

def parallel_parse_blast_file(blast_file, blast_type = 'blastn', min_identity = 0, min_ppos = 0,  min_aln_len = 0, nproc = 1, chunk_size = 250000000):

    blast_family_list = []


    def get_line_chunk(file_handle, chunk_size):
        line_chunk = file_handle.readlines(chunk_size)
        if not line_chunk:
            line_chunk = file_handle.readlines()
        yield line_chunk

    filter_dict = {'min_identity': min_identity, 'min_similarity': min_ppos, 'min_aln_len': min_aln_len}


    p_pool = multiprocessing.Pool(nproc)


    in_handle = open(blast_file)
    result_list = p_pool.starmap(parse_blast_file_worker, zip(get_line_chunk(in_handle, chunk_size), itertools.repeat(filter_dict)))
    p_pool.close()
    p_pool.join()

    result_list = list(itertools.chain(result_list))[0]
    print('result_list_len', len(result_list))
    result_list.sort(key= lambda x: x[0])

    for key, group in itertools.groupby(result_list, lambda x: x[0]):
        new_blast_family = Blast_Family(blast_type, filter_dict, key, (0, 0))
        for blast_hit in [x[1] for x in list(group)]:
            new_blast_family.add_blast(blast_hit)
            blast_hit.add_family(new_blast_family)
        new_blast_family.parents = new_blast_family.blast_list[0].parents
        new_blast_family.parent_lengths = new_blast_family.blast_list[0].parent_lengths
        new_blast_family.equalize()
        blast_family_list.append(new_blast_family)

    print('blast_family_len', len(blast_family_list))
    return blast_family_list

#The following code only executes when its run as a script
if __name__ == '__main__':
    pass






