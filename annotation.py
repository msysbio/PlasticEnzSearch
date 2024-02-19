import os
import glob
import logging
import Bio.Blast.NCBIWWW
import Bio.Blast.NCBIXML
from Bio import SeqIO
import pandas as pd
from datetime import datetime

error2="""
 __ __  __  __  __  
|_ |__)|__)/  \|__) 
|__| \ | \ \__/| \  
                                     
"""

def run_blast_on_sequences(fasta_file_path, output_excel_file_path):
    # Create a logger
    logging.basicConfig(filename='annotation.log', level=logging.ERROR)
    try:
        # Check if the file is empty
        if os.path.getsize(fasta_file_path) == 0:
            print(f"{datetime.now()} - The file {fasta_file_path} is empty. Skipping...")
            return

        # Read sequences from fasta file
        sequences = SeqIO.parse(open(fasta_file_path),'fasta')

        result_list = []

        for record in sequences:
            print(f"{datetime.now()} - Running BLAST for {record.id}, this may take a while...")
            # Run BLAST and parse the result
            result_handle = Bio.Blast.NCBIWWW.qblast("blastp", "nr", record.seq)
            blast_record = Bio.Blast.NCBIXML.read(result_handle)
            
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
                        print(f"{datetime.now()} - BLAST for {record.id} finished running with success! :-)")

        # Create a DataFrame from the results
        df = pd.DataFrame(result_list, columns=["Fasta header", "Functional annotation", "Bit-score", "E-value", "Query Cover", "Percent Identity", "Accession"])

        # Write the DataFrame to an Excel file
        df.to_excel(output_excel_file_path, index=False)

    except Exception as e:
        logging.error(f"Error in run_blast_on_sequences function: {e}")
        print(error2)
        print("The annotation step was interrupted :-( , check the annotation.log file for more information")


def blast_search_in_directory(directory):
    try:
        fasta_files = glob.glob(os.path.join(directory, "*_hmm_output.fasta"))

        for fasta_file in fasta_files:
            output_file_name = os.path.splitext(fasta_file)[0] + "_annotation.xlsx"
            print(f"{datetime.now()} - Processing {fasta_file}...")
            run_blast_on_sequences(fasta_file, output_file_name)
    
    except Exception as e:
        logging.error(f"Error in blast_search_in_directory function: {e}")
        print(error2)
        print("The annotation step was interrupted :-( , check the annotation.log file for more information")

blast_search_in_directory("/Users/u0145079/Library/CloudStorage/OneDrive-KULeuven/Desktop/PlasticDaphnia/Plastic_tool/Plastic_tool/temps/")