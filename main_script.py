import re
from sys import argv
from libs import load_multifasta, from_dir, edge_primer_maker, \
    inner_primers_maker, ensure_dir, write_genbank, annot_primer
from Bio.Alphabet import generic_dna


print "\n",\
"##################################################\n",\
"### NovoReSeq v. 0.1                           ###\n",\
"### Copyright 2012 Geraldine A. Van der Auwera ###\n",\
"##################################################\n"\

if len(argv) > 1 and argv[1] == '-h':
    print "Basic usage: \n",\
    "$ python main_script.py data_dir\n"
    exit()

data_dir = "data/"+argv[1]+"/"
gbk_dir = "data/"+argv[1]+"/genbank/"

ensure_dir(gbk_dir)

status_file = open(data_dir+"dataset_info.txt", 'w')

genomes = from_dir(data_dir, re.compile(r'.*\.fas.*'))

primers = []

for genome in genomes:
    contigs = load_multifasta(data_dir+genome)

    g_name = genome[:genome.find(".fas")]
    genome_dir = gbk_dir+"/"+g_name+"/"
    ensure_dir(genome_dir)

    statline = "\n"+g_name+"\n"
    status_file.write(statline)
    print "\t", statline

    c_count = 0
    for contig in contigs:

        c_count +=1
        statline = "\t".join(["C"+str(c_count), contig.id])
        status_file.write(statline)
        print "\t", statline,

        if len(contig.seq) > 300:
            # outward primers
            edge_primers = edge_primer_maker(contig, c_count, g_name)
            for primer in edge_primers:
                primers.append(primer)
                # add primer to sequence file
                contig.features.append(annot_primer(primer))
            # inner primers
            inner_primers = inner_primers_maker(contig, c_count, g_name)
            for primer in inner_primers:
                primers.append(primer)
                # add primer to sequence file
                contig.features.append(annot_primer(primer))
            status_file.write("\tOK\n")
            print 'OK'
        else:
            status_file.write("\tskipped (too short)\n")
            print "skipped (too short)"

        # make a genbank file of the contig
        gbk_file = "".join([genome_dir, g_name, "_C", str(c_count), ".gbk"])
        contig.name = "".join([g_name, "_C", str(c_count)])
        contig.id = "".join([g_name, "_C", str(c_count)])
        contig.seq.alphabet = generic_dna
        write_genbank(gbk_file, contig)

    print "\n"

# create output file
output_file = open(data_dir+'primers.txt', 'w')
for primer in primers:
    #print primer
    output_file.write("\t".join([primer['name'], primer['seq'], "\n"]))
output_file.close()
