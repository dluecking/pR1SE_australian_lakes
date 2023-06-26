import Bio.Entrez as Entrez

def get_ipg_nuccore(prot_acc):
    # download the record from IPG
    ipg_record = Entrez.efetch(id=prot_acc, db='ipg').read().decode()
    
    # Get the nucleotide ID from the first line of the resulting table
    nt_acc = ipg_record.split("\n")[1].split("\t")[2]
    
    # Get the nucleotide record
    nuc_record = Entrez.efetch(id=nt_acc, db='nuccore', rettype='fasta').read()
    
    # Parse into a dict - key is ID, value is sequence
    fasta = {nuc_record.split("\n")[0]: "".join(nuc_record.split("\n")[1:])}
    return (fasta)

# Step 1: Read accs
file = open("accs_to_download.txt", "r")
unique_accessions = file.readlines()

# Step 2: download genome
for accession in unique_accessions:
    dictionary = get_ipg_nuccore(accession)
    
    genome_acc = list(dictionary.keys())[0]
    genome_acc = genome_acc.split(" ")[0]
    genome_acc = genome_acc.replace(">", "")
    print(genome_acc)
    # Step 3: Write key-value pairs to a file
    output_filename = f'../known_hosts/{genome_acc}.fasta'
    with open(output_filename, 'w') as fasta_file:
        for key, value in dictionary.items():
            fasta_file.write(f">{key}\n{value}\n")
