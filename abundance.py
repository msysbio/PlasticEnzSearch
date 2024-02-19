#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 22 15:09:30 2023

@author: u0145079
"""
from Bio import SeqIO
from Bio.Seq import Seq
import os
from utilities import check_dependencies,spinning_cursor_task
import threading
import subprocess
import time
import os

def check_files_and_folders(args):
    # Check if output directory exists
    if not os.path.isdir(args.output):
        print(f"Output directory {args.output} does not exist.")
        return False
    
    # Check if temps directory exists
    global temps_folder
    temps_folder = os.path.join(args.output, 'temps')
    if not os.path.isdir(temps_folder):
        print(f"Temps directory {temps_folder} does not exist.")
        return False
    
    # Check if contigs file exists
    if not os.path.isfile(args.contigs):
        print(f"Contigs file {args.contigs} does not exist.")
        return False
    
    # Check if bam/sam files exist
    if ',' in args.mappings:
        mapping_files = args.mappings.split(',')
    else:
        mapping_files = [args.mappings]
    
    for mapping_file in mapping_files:
        if not os.path.isfile(mapping_file.strip()):
            print(f"Mapping file {mapping_file.strip()} does not exist.")
            return False

    # Check if aa and nt files exist
    contigs_file = os.path.basename(args.contigs)
    contigs_base = os.path.basename(args.contigs).split(".")[0]
    aa_file = os.path.join(temps_folder, f"{contigs_base}.faa")
    nt_file = os.path.join(temps_folder, f"{contigs_base}.ffn")
    
    if not os.path.isfile(aa_file):
        print(f"AA file {aa_file} does not exist.")
        return False
    
    if not os.path.isfile(nt_file):
        print(f"NT file {nt_file} does not exist.")
        return False
    
    # If all checks passed
    return True

def search_for_hits(hits_file_path, genes_file_path, output_file_path):
    # Read all sequence IDs from hits.fa into a list
    hits_ids = [rec.id for rec in SeqIO.parse(hits_file_path, "fasta")]

    # Open output file
    with open(output_file_path, "w") as output_file:
        # Open genes.ffn and iterate through each sequence
        with open(genes_file_path, "r") as genes_file:
            for record in SeqIO.parse(genes_file, "fasta"):
                # If the sequence ID (without the extra details) is in the list of hits
                # write it to the output file with original header from genes.ffn
                if '_'.join(record.id.split('_')[0:3]) in hits_ids:
                    SeqIO.write(record, output_file, "fasta")

def correct_ffn_file(input_file, output_file):
    with open(input_file, 'r') as original, open(output_file, 'w') as corrected:
        for record in SeqIO.parse(original, 'fasta'):
            header = record.description.split(' ')
            contig_id, start, end, strand = header[0], int(header[2]), int(header[4]), int(header[6])
            if strand == -1:
                record.seq = record.seq.reverse_complement()
            record.id = record.description  # Important: Ensure record.id is the same as record.description
            SeqIO.write(record, corrected, "fasta")  # Use SeqIO.write to write the record in FASTA format

def ffn_to_saf(input_ffn_file, output_saf_file):
    with open(input_ffn_file, "r") as ffn, open(output_saf_file, "w") as saf:
        # Write the SAF header
        saf.write("\t".join(["GeneID", "Chr", "Start", "End", "Strand"]) + "\n")
        
        # Parse each FASTA entry in the input file
        for record in SeqIO.parse(ffn, "fasta"):
            # Extract information from the record description
            fields = record.description.split("#")
            gene_id = fields[0].strip()
            start = fields[1].strip()
            end = fields[2].strip()
            strand = fields[3].strip()
            if strand == '1':
                strand = '+'
            elif strand == '-1':
                strand = '-'
            
            # Remove last suffix from Chr
            chr = "_".join(gene_id.split("_")[:-1])
            
            # Write the SAF entry
            saf.write("\t".join([gene_id, chr, start, end, strand]) + "\n")

def annotation(args):
    if not check_files_and_folders(args):
        return
    
    check_dependencies("featureCounts")
    
    # Get the base name of the contigs file
    contigs_base = os.path.basename(args.contigs).split(".")[0]

    # Determine the plastic types to be considered
    motif_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "motifs")
    if isinstance(args.plastic, str) and args.plastic.lower() != "all":
        plastic_types = args.plastic.split(',')
    elif isinstance(args.plastic, str) and args.plastic.lower() == "all":
        plastic_types = [os.path.splitext(file)[0] for file in os.listdir(motif_dir) if file.endswith('.hmm')]

    # Process each plastic type
    for plastic_type in plastic_types:
        temp_folder_path = os.path.join(args.output, 'temps', plastic_type)
        fasta_files = [file for file in os.listdir(temp_folder_path) if file.endswith('.fasta')]
        
        if fasta_files:
            # Construct file paths
            print(fasta_files)
            hits_file_path = os.path.join(temp_folder_path, fasta_files[0])
            genes_file_path = os.path.join(temps_folder, f"{contigs_base}.ffn")
            print(genes_file_path)
            output_file_base = os.path.basename(fasta_files[0]).split(".")[0]
            output_file_path = os.path.join(temp_folder_path, f"{output_file_base}.fnn")
            
            # Process hits and correct the .ffn file
            search_for_hits(hits_file_path, genes_file_path, output_file_path)
            corrected_file_path = os.path.join(temp_folder_path, f"{output_file_base}-corrected.fnn")
            correct_ffn_file(output_file_path, corrected_file_path)
            
            # Convert corrected .ffn file to .saf
            saf_file_path = os.path.join(temp_folder_path, f"{output_file_base}.saf")
            ffn_to_saf(corrected_file_path, saf_file_path)
        else:
            print(f"No .fasta files found in {temp_folder_path}. Skipping this folder.")
            continue

        
        #use FeatureCount to quantify
        
        #SAM/BAM file may be a file list so you will need to loop through them
        #one by one
        
        #There should one .saf file per plastic folder.      
        
        # Split the mappings argument into a list of individual file paths
        mapping_files = args.mappings.split(',')
        for mapping_file in mapping_files:
            mapping_file = mapping_file.strip()  # Remove any leading/trailing whitespace
            fc_input = saf_file_path
            fc_output = os.path.join(temp_folder_path, f"{plastic_type}_{os.path.basename(mapping_file)}_counts.out")
            fc_command = f"featureCounts -a {fc_input} -F SAF -o {fc_output} {mapping_file}"
            log_file = os.path.join(temp_folder_path, f"{os.path.basename(mapping_file)}_featureCounts.log")
        
            # Specify the file to capture program output
            program_output_file = os.path.join(temp_folder_path, f"{os.path.basename(mapping_file)}_featureCounts.out")
        
            task_done = threading.Event()
        
            with open(log_file, 'w') as f, open(program_output_file, 'w') as p_out:
                process = subprocess.Popen(fc_command, shell=True, stdout=p_out, stderr=f)
        
            t = threading.Thread(target=spinning_cursor_task, args=(task_done,'featureCounts'))
            t.start()
        
            while process.poll() is None:
                time.sleep(0.1)
            task_done.set()
            t.join()

            #Calculate statistsc on the abundance data: average, standard div, sum.

            # Folder where fc_output files are stored
            fc_folder_path = temp_folder_path  # Replace with your folder path
            tsv_output_file = os.path.join(temp_folder_path, 'mapping_summary.tsv')  # Replace with your desired output file path
            
            with open(tsv_output_file, 'w') as tsv_file:
                tsv_file.write('\t'.join(["sample name", "reads mapped", "total reads", "proportion"]) + '\n')
            
                # Iterate through all files in the directory
                for filename in os.listdir(fc_folder_path):
                    # Check if the file is a featureCounts output file
                    if filename.endswith("_counts.out.summary"):
                        print(f"Processing file: {filename}")  # Debugging print
                        # Open the file
                        with open(os.path.join(fc_folder_path, filename), 'r') as fc_output:
                            # Initialize counts
                            reads_mapped = 0
                            total_reads = 0
            
                            # Iterate through each line in the file
                            for line in fc_output:
                                # Ignore lines that are not part of the summary
                                if not line.startswith('Assigned') and not line.startswith('Unassigned'):
                                    continue
                                
                                # Extract count from the line
                                count = int(line.split('\t')[1].strip())
                                
                                # Add the count to the appropriate variable
                                if line.startswith('Assigned'):
                                    reads_mapped = count
                                total_reads += count
            
                            # Calculate proportion
                            proportion = reads_mapped / total_reads
            
                            # Get sample name from filename
                            sample_name = filename.split('_counts.out.summary')[0]
            
                            # Write row to TSV file
                            tsv_file.write('\t'.join([sample_name, str(reads_mapped), str(total_reads), str(proportion)]) + '\n')

        
        
        
        
        
        
        
        
        
        
        
        
    
