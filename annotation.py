import os
import glob
import logging
from Bio import SeqIO
import pandas as pd
from datetime import datetime
import logging
from concurrent.futures import ThreadPoolExecutor
import blast_handler

# Create a dictionary to store BLAST results
blast_results = {}

# Create a dictionary to store unique records
unique_records = {}

def collect_unique_sequences(directory):
    # Walk through the directory and its subdirectories
    for dirpath, dirnames, filenames in os.walk(directory):
        # For each file in the current directory
        for filename in filenames:
            # If the file is a .fasta file
            if filename.endswith('.fasta'):
                # Construct the full file path
                fasta_file_path = os.path.join(dirpath, filename)
                # Read sequences from fasta file
                sequences = list(SeqIO.parse(open(fasta_file_path),'fasta'))
                for record in sequences:
                    # Add record to the dictionary with sequence as the key
                    unique_records[str(record.seq)] = record

def run_blast_on_unique_sequences():
    try:
        # Create a pool of threads
        with ThreadPoolExecutor() as executor:
            # For each unique record, run BLAST
            futures = []
            for seq, rec in unique_records.items():
                if seq not in blast_results:
                    # If not, run BLAST and store the result
                    future = executor.submit(run_blast, rec)
                    futures.append(future)

            # Wait for all BLAST searches to finish
            for future in futures:
                result_list = future.result()
                # Use sequence as the key in blast_results
                blast_results[seq].extend(result_list)

    except Exception as e:
        logging.error(f"Error: {e}")

def run_blast(record):
    logging.info(f"{datetime.now()} - Running BLAST for {record.id}, this may take a while...")
    # Run BLAST and parse the result
    blast_record = blast_handler.blast_parse(record)

    result_list = []
    # Check there are at least 5 alignments
    if len(blast_record.alignments) >= 5:
        # Get the top 5 hits
        top_hits = blast_record.alignments[:5]

        for hit in top_hits:
            top_hsp = hit.hsps[0]

            # Check e-value
            if top_hsp.expect < 0.01:
                # Extract required data
                sequence_id = record.id
                protein_hit = hit.hit_def
                bit_score = top_hsp.bits
                e_value = top_hsp.expect
                query_cover = (top_hsp.align_length / len(record.seq)) * 100
                perc_identity = (top_hsp.identities / top_hsp.align_length) * 100
                accession = hit.accession

                result_list.append([sequence_id, protein_hit, bit_score, e_value, query_cover, perc_identity, accession])
    
    logging.info(f"{datetime.now()} - BLAST for {record.id} finished running.")
    return result_list

def create_excel_files(directory):
    try:
        fasta_files = glob.glob(os.path.join(directory, "*.fasta"))

        for fasta_file in fasta_files:
            output_file_name = os.path.splitext(fasta_file)[0] + "_annotation.xlsx"
            logging.debug(f"{datetime.now()} - Processing {fasta_file}...")

            # Read sequences from fasta file
            sequences = list(SeqIO.parse(open(fasta_file),'fasta'))

            # Create a list to store the results for this file
            result_list = []

            # For each sequence in the file, get the BLAST result from the dictionary
            for record in sequences:
                result_list.extend(blast_results[str(record.seq)])

            # Create a DataFrame from the results
            df = pd.DataFrame(result_list, columns=["Fasta header", "Functional annotation", "Bit-score", "E-value", "Query Cover", "Percent Identity", "Accession"])

            # Write the DataFrame to an Excel file
            df.to_excel(output_file_name, index=False)

    except Exception as e:
        logging.error(f"Error: {e}")

def blast_search(p):

    collect_unique_sequences(p.temps)
    run_blast_on_unique_sequences()

    for plastic in p.plastic_list:
        create_excel_files(os.path.join(p.temps, plastic))
