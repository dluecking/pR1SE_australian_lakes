from Bio import SeqIO
import pandas as pd
import os
import re

# Read contig_df
contig_df = pd.read_csv("novel_hosts_df.tsv")

# Read all input fasta files and keep relevant hits
relevant_seqs = []

fasta_files = [os.path.join("../tal/4_Dom/", file) for file in os.listdir("../tal/4_Dom/") if file.endswith(".fasta")]

for fasta_file in fasta_files:
    for record in SeqIO.parse(fasta_file, "fasta"):
        short_name = re.sub("\_length.*", "", record.name)
        if short_name in contig_df["contig"].values and record.name not in [record.name for record in relevant_seqs]:
            relevant_seqs.append(record)
            print(short_name)

# Save to file
SeqIO.write(relevant_seqs, "novel_host_sequences.fasta", "fasta")
