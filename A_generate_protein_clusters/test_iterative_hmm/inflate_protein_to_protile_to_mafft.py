import argparse
import os
import subprocess
import shutil
import re


# Step 1: Parse input and output file paths using argparse
parser = argparse.ArgumentParser(description='Protein sequence analysis script')
parser.add_argument('--protein_file', required=True, help='Path to protein fasta file (e.g. ORF42.faa)')
parser.add_argument('--output_file', required=True, help='Path of the final alignment file (e.g. path/to/file/ORF42.msa)')
parser.add_argument('--database_path', required=True, help='Path to the database')
parser.add_argument('--keep_intermediate', action='store_true', help='Keep intermediate files')
args = parser.parse_args()

# Extract the ORF name from the protein file name
orf_name = re.findall("ORF\d*", args.protein_file)[0]

# Create a temporary folder for intermediate files
tmp_folder = f'tmp_{orf_name}'
os.makedirs(tmp_folder, exist_ok=True)

try:
    # Step 1: Align protein sequences using MAFFT
    alignment_file = f'{tmp_folder}/{orf_name}_alignment1.msa'
    mafft_command = f'mafft --localpair --reorder {args.protein_file} > {alignment_file}'
    print('# STEP 1 #################################')
    print('Running MAFFT command:', mafft_command)
    subprocess.run(mafft_command, shell=True)

    # Step 2: Create HMM profile from the alignment
    profile_file = f'{tmp_folder}/{orf_name}_profile.hmm'
    hmmbuild_command = f'hmmbuild --amino {profile_file} {alignment_file}'
    print('# STEP 2 #################################')
    print('Running hmmbuild command:', hmmbuild_command)
    subprocess.run(hmmbuild_command, shell=True)

    # Step 3: Search HMM profile against the database
    hmmsearch_result_file = f'{tmp_folder}/{orf_name}_hmmsearch_result.domtableout.txt'
    hmmsearch_command = f'hmmsearch --domtblout {hmmsearch_result_file} --cpu 16 {profile_file} {args.database_path}'
    print('# STEP 3 #################################')
    print('Running hmmsearch command:', hmmsearch_command)
    subprocess.run(hmmsearch_command, shell=True)

    # Step 4: Run Rscript to filter and download new accessions
    filter_and_download_script = 'filter_and_download_new_accessions.R'
    output_file = f'{tmp_folder}/{orf_name}_output.fasta'
    # the script needs four arguments:
    # Rscript filter_and_download_new_accessions.R <ORF> <DOMTBLOUT_FILE> <ORIGINAL_PROTEIN_IN_FILE> <PROTEIN_OUT_FILE>   
    rscript_command = f'Rscript {filter_and_download_script} {orf_name} {hmmsearch_result_file} {args.protein_file} {output_file}'
    print('# STEP 4 #################################')
    print('Running Rscript command:', rscript_command)
    subprocess.run(rscript_command, shell=True)

    # Step 5: Align the resulting protein file again using MAFFT
    final_alignment_file = f'{orf_name}_final_alignment.msa'
    mafft_command = f'mafft --localpair --reorder {output_file} > {final_alignment_file}'
    print('# STEP 5 #################################')
    print('Running MAFFT command for final alignment:', mafft_command)
    subprocess.run(mafft_command, shell=True)

    # Move the final alignment file to the specified output path
    print('# DONE ###################################')
    print(f'We are done with {orf_name}')
    final_output_file = os.path.join(args.output_path, final_alignment_file)
    shutil.move(final_alignment_file, args.output_path)

    print('Script execution complete!')

finally:
    # Delete the temporary folder unless the --keep-intermediate flag is provided
    if not args.keep_intermediate:
        shutil.rmtree(tmp_folder)

