import argparse
from Bio import SeqIO

# Create the argument parser
parser = argparse.ArgumentParser(description="Compare accessions in two FASTA files.")
parser.add_argument("file1", help="Path to the first FASTA file")
parser.add_argument("file2", help="Path to the second FASTA file")

# Parse the command line arguments
args = parser.parse_args()

def get_accessions(file):
    accessions = set()
    with open(file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            accession = record.id
            accessions.add(accession)
    return accessions

def compare_accessions(file1, file2):
    accessions_file1 = get_accessions(file1)
    accessions_file2 = get_accessions(file2)
    unique_accessions = accessions_file1 - accessions_file2
    for accession in unique_accessions:
        print(accession)

# Call the compare_accessions function with the provided file paths
compare_accessions(args.file1, args.file2)
