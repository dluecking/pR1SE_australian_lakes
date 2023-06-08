import argparse
import subprocess
from Bio import SeqIO, Entrez
import re
import shutil
import os
import time

################################################################################
# constants
EMAIL = "dom.luecking@gmail.com"
API_KEY = "ed38a68e4cac02507e4bc585e8913bab5a08"
MAX_ITERATIONS = 5

# variables
known_accessions = []
iteration_counter = 0

# Define command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="Input FASTA file")
parser.add_argument("output_file", help="Output FASTA file")
parser.add_argument("--db", help="Path to database")
parser.add_argument("--evalue", type=float, default=0.01, help="E-value cutoff for Diamond")
parser.add_argument("--bitscore", type=float, default=50, help="Bitscore cutoff for Diamond")
parser.add_argument("--threads", default=16, help="Number of threads to use for Diamond")
args = parser.parse_args()


################################################################################
def get_record(accession):
    Entrez.email = EMAIL
    Entrez.api = API_KEY
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text", api_key=Entrez.api)
    record = SeqIO.read(handle, "fasta")
    return record

def get_new_hit_accessions(blast_result, known_accessions):
    new_hit_accessions = set()
    with open(blast_out_file) as blast_result:
        for line in blast_result:
            fields = line.strip().split('\t')
            # field 10 = evalue
            # field 11 = bitscore
            if float(fields[10]) <= args.evalue and float(fields[11]) >= args.bitscore and fields[1] not in known_accessions:
                new_hit_accessions.add(fields[1])
                
    l = len(new_hit_accessions)
    print(f'Number of new accessions: {l}')
    return(new_hit_accessions)

################################################################################
# get ORF
ORF = re.findall('ORF\d*', args.input_file)[0]
# get DB
DB = os.path.basename(args.db)

print(args)

# copy input file to output file
tmp_file = args.output_file + ".tmp"
shutil.copy(args.input_file, tmp_file)


while iteration_counter < MAX_ITERATIONS:
    # blastp
    blast_out_file = f'../02_blastp_results/{ORF}_it{iteration_counter}_{DB}_blast_results.tsv'
    subprocess.run(['diamond', 'blastp', '-d', args.db, '-q', tmp_file, '-o', blast_out_file, '--faster', '--threads', args.threads])

    # add accessions present in input file to known_accessions
    for record in SeqIO.parse(tmp_file, 'fasta'):
        known_accessions.append(record.id)
    
    # get the new hits based on the blastp result, excluding the new 
    new_hit_accessions = get_new_hit_accessions(blast_out_file, known_accessions)

    # get the NCBI record (including sequence) for each new hit
    new_records = []
    for acc in new_hit_accessions:
        new_records.append(get_record(acc))
        time.sleep(1)
        
    
    if new_records:
        # append new records to tmpfile
        with open(tmp_file, "a") as tmpfile:
            SeqIO.write(new_records, tmpfile, "fasta")
    else:
        print("No new accessions found!")
        break

    iteration_counter += 1
    print("##########################################################")
    print(f"Following interation: {iteration_counter}")

print('No new hits found or MAX_ITERATION was reached. Tterative blastp complete.')
shutil.copy(tmp_file, args.output_file)
