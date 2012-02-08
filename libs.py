import re
import subprocess
from os import listdir, path, makedirs
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def load_multifasta(seqfile):
    """Load multiple records from Fasta file."""
    input_handle = open(seqfile, 'rU')
    multifasta = SeqIO.parse(input_handle, 'fasta')
    fasta_list = list(multifasta)
    for record in fasta_list:
        assert record.id
    return fasta_list

def from_dir(ori_dir, pattern):
    """Load filenames in a directory using a regex."""
    contents = listdir(ori_dir)
    filenames = []
    for item in contents:
        match = re.match(pattern, item)
        if match:
            filenames.append(item)
    return filenames

def run_primer3(file):
    """Make external call to Primer3."""
    comps = ["primer3_core", file]
    cline = " ".join(comps)
    try:
        child = subprocess.Popen(str(cline), stdout=subprocess.PIPE,
                                 shell=True)
        output, error = child.communicate()
    except: raise
    else: return output

def make_boulderIO(seq_id, template, task, left, right, pick_num,
                   interval, lead, acc, spacing,
                   t_start, t_len, x_start, x_len):
    return [
        'SEQUENCE_ID='+str(seq_id),
        'SEQUENCE_TEMPLATE='+template,
        'SEQUENCE_TARGET='+str(t_start)+','+str(t_len),
        'SEQUENCE_EXCLUDED_REGION='+str(x_start)+','+str(x_len),
        'PRIMER_TASK='+task,
        'PRIMER_NUM_RETURN='+str(pick_num),
        'PRIMER_SEQUENCING_INTERVAL='+str(interval),
        'PRIMER_SEQUENCING_LEAD='+str(lead),
        'PRIMER_SEQUENCING_ACCURACY='+str(acc),
        'PRIMER_SEQUENCING_SPACING='+str(spacing),
        'PRIMER_PICK_LEFT_PRIMER='+str(left),
        'PRIMER_PICK_RIGHT_PRIMER='+str(right),
        'PRIMER_MAX_NS_ACCEPTED=0',
        'PRIMER_EXPLAIN_FLAG=1',
        '='
    ]

def parse_boulder_single(raw_output, ord_res):
    """Parse output of Primer3 for single primer."""
    output_list = raw_output.split("\n")
    if int(output_list[ord_res][-1]) > 0:
        # retrieve primer location
        sequence = output_list[20][output_list[20].find('=')+1:]
        position = output_list[21][output_list[21].find('=')+1:].split(",")
        start = int(position[0])
        length = int(position[1])
    else:
        # should parse for error but this will do for now
        sequence, start, length = False, False, False
    return sequence, start, length

def parse_boulder_mass(raw_output):
    """Parse output of Primer3 for multiple primers."""
    output_list = raw_output.split("\n")
    fwd_seqs = [x[x.find('=')+1:]
                for x in output_list if x.find('SEQUENCE=') != -1
                                      and x.find('LEFT') != -1]
    rev_seqs = [x[x.find('=')+1:]
                for x in output_list if x.find('SEQUENCE=') != -1
                                      and x.find('RIGHT') != -1]
    fwd_locs = [x[x.find('=')+1:].split(",")
                for x in output_list if x.find('PRIMER_LEFT') != -1
                                      and x.find('END') == -1
                                      and x.find('NUM') == -1
                                      and x.find('PENALTY') == -1
                                      and x.find('TM') == -1
                                      and x.find('GC') == -1
                                      and x.find('SEQUENCE') == -1
                                      and x.find('EXPLAIN') == -1
                                      and x.find('SELF') == -1 ]
    rev_locs = [x[x.find('=')+1:].split(",")
                for x in output_list if x.find('PRIMER_RIGHT') != -1
                                      and x.find('END') == -1
                                      and x.find('NUM') == -1
                                      and x.find('PENALTY') == -1
                                      and x.find('TM') == -1
                                      and x.find('GC') == -1
                                      and x.find('SEQUENCE') == -1
                                      and x.find('EXPLAIN') == -1
                                      and x.find('SELF') == -1 ]
    try:
        assert len(fwd_seqs) == len(fwd_locs)
    except AssertionError:
        print "\nERROR in parsed FWD counts\n"
        print fwd_seqs
        print fwd_locs
    try:
        assert len(rev_seqs) == len(rev_locs)
    except AssertionError:
        print "\nERROR in parsed REV counts\n"
        print rev_seqs
        print rev_locs
    return fwd_seqs, fwd_locs, rev_seqs, rev_locs

def edge_primer_maker(contig, c_count, g_name):
    """Design outward sequencing primers."""
    name_base = g_name+"_C"+str(c_count)+"_OUT_"
    task = 'pick_sequencing_primers'
    left_record_txt = make_boulderIO(name_base+'LT', str(contig.seq),
                                     task, 0, 1, 1, 100, 10, 50, 200,
                                     0, 150, 250, len(contig)-250)
    right_record_txt = make_boulderIO(name_base+'RT', str(contig.seq),
                                      task, 1, 0, 1, 100, 10, 50, 200,
                                      len(contig)-150, 150,
                                      0, len(contig)-250)
    # make tempfile
    filename = "tempfile.txt"
    # make left
    open(filename, 'w').write("\n".join(left_record_txt))
    left_result_raw = run_primer3(filename)
    seq, start, length = parse_boulder_single(left_result_raw, 16)
    left_primer = {'name': name_base+'LT',
                   'seq': seq,
                   'start': start,
                   'len': length,
                   'direction': 'rev'}
    # make right
    open(filename, 'w').write("\n".join(right_record_txt))
    right_result_raw = run_primer3(filename)
    seq, start, length = parse_boulder_single(right_result_raw, 15)
    right_primer = {'name': name_base+'RT',
                    'seq': seq,
                    'start': start,
                    'len': length,
                    'direction': 'fwd'}
    return [left_primer, right_primer]

def inner_primers_maker(contig, c_count, g_name): # aim for overlap by 250
    """Design inner sequencing primers."""
    name_base = g_name+"_C"+str(c_count)+"_IN_"
    task = 'pick_sequencing_primers'
    primers = []
    # make tempfile
    filename = "tempfile.txt"
    # run mass primer design
    spacing = 1100
    if len(contig.seq) < spacing:
        spacing = len(contig.seq)/2
        num_return = 5
    else:
        num_return = len(contig.seq)/(spacing-250)
    record_txt = make_boulderIO(name_base+'ALL', str(contig.seq),
                                task, 1, 1, num_return,
                                250, 20, 20, spacing,
                                100, len(contig)-200, 0, 0)
    open(filename, 'w').write("\n".join(record_txt))
    mass_results_raw = run_primer3(filename)
    fwd_seqs, fwd_locs, rev_seqs, rev_locs = \
                                        parse_boulder_mass(mass_results_raw)
    # unpack primer data (should refactor)
    counter = 0
    while counter < len(fwd_seqs):
        fwd_primer = {'name': name_base+'FWD_'+str(counter+1),
                      'seq': fwd_seqs[counter],
                      'start': int(fwd_locs[counter][0]),
                      'len': int(fwd_locs[counter][1]),
                      'direction': 'fwd'}
        primers.append(fwd_primer)
        counter +=1
    counter = 0
    while counter < len(rev_seqs):
        rev_primer = {'name': name_base+'REV_'+str(counter+1),
                      'seq': rev_seqs[counter],
                      'start': int(rev_locs[counter][0]),
                      'len': int(rev_locs[counter][1]),
                      'direction': 'rev'}
        primers.append(rev_primer)
        counter +=1
    return primers

def ensure_dir(dir_path):
    """Check that the directory exists; if not, create it."""
    abs_path = path.abspath(dir_path)
    if not path.exists(abs_path):
        try: makedirs(abs_path)
        except Exception as message:
            status = 1
            # TODO: make graceful fail or request input if interactive mode
        else:
            message = 'created path'
            status = 0
    else:
        message = 'path exists'
        status = 0
    report = {'message': message, 'status': status}
    return abs_path, report

def write_genbank(filename, seqrecords):
    """Write GenBank file."""
    output_handle = open(filename, 'w')
    counter = SeqIO.write(seqrecords, output_handle, 'genbank')
    output_handle.close()
    return counter

def annot_primer(primer):
    """Generate a primer annotation."""
    if primer['direction'] == 'fwd':
        start_pos = primer['start']
        end_pos = primer['start']+primer['len']
        strand_pos = 1
    else:
        start_pos = primer['start']-primer['len']+1
        end_pos = primer['start']+1
        strand_pos = -1
    feat_loc = FeatureLocation(start_pos, end_pos)
    feature = SeqFeature(location=feat_loc,
                         strand = strand_pos,
                         id=primer['name'],
                         type='primer_bind')
    return feature
