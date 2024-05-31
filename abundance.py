from Bio import SeqIO
import os
import utilities
import subprocess
import time
import multiprocessing
from functools import partial
import logging

def check_translate_result(p):

    # Check if aa and nt files exist
    contigs_file = os.path.basename(p.contigs)
    contigs_base = contigs_file.split(".")[0]
    aa_file = os.path.join(p.temps, f"{contigs_base}.faa")
    nt_file = os.path.join(p.temps, f"{contigs_base}.ffn")
    
    if not os.path.isfile(aa_file):
        logging.warning(f"AA file {aa_file} does not exist.")
        return False
    
    if not os.path.isfile(nt_file):
        logging.warning(f"NT file {nt_file} does not exist.")
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

def annotation(p):
    if not check_translate_result(p):
        return

    # Determine the plastic types to be considered
    motif_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "motifs")
    if isinstance(p.plastic, str) and p.plastic.lower() != "all":
        plastic_types = p.plastic.split(',')
    elif isinstance(p.plastic, str) and p.plastic.lower() == "all":
        plastic_types = [os.path.splitext(file)[0] for file in os.listdir(motif_dir) if file.endswith('.hmm')]
    
    #create tsv files for each mapping file
    for mapping_file in p.mappings.split(','):
        sample_name = mapping_file.split('.')[0].split('/')[-1]
        tsv_output_file = os.path.join(p.temps, f'{sample_name}.tsv')
        with open(tsv_output_file, 'w') as tsv_file:
            tsv_file.write('\t'.join(["plastic name", "reads mapped", "total reads", "proportion", "rpkm"]) + '\n')


    # Create a new function that has `p` already filled in
    featurecounts_p = partial(featurecounts, p=p)

    # Start the tasks
    pool = multiprocessing.Pool(processes=utilities.get_logical_cores() if p.multiprocessing else 1)
    pool.map(featurecounts_p, plastic_types)

    # Wait for the tasks to finish
    pool.close()
    pool.join()

def featurecounts(plastic_type, p):
    try:
        # Get the base name of the contigs file
        contigs_base = os.path.basename(p.contigs).split(".")[0]

        temp_folder_path = os.path.join(p.temps, plastic_type)
        fasta_files = [file for file in os.listdir(temp_folder_path) if file.endswith('.fasta')]
        
        if fasta_files:
            # Construct file paths
            hits_file_path = os.path.join(temp_folder_path, fasta_files[0])
            genes_file_path = os.path.join(p.temps, f"{contigs_base}.ffn")
            output_file_base = os.path.basename(fasta_files[0]).split(".")[0]
            output_file_path = os.path.join(temp_folder_path, f"{output_file_base}.fnn")
            
            # Process hits and correct the .ffn file
            search_for_hits(hits_file_path, genes_file_path, output_file_path)
            corrected_file_path = os.path.join(temp_folder_path, f"{output_file_base}-corrected.fnn")
            correct_ffn_file(output_file_path, corrected_file_path)
            
            # Convert corrected .ffn file to .saf
            saf_file_path = os.path.join(temp_folder_path, f"{output_file_base}.saf")
            ffn_to_saf(corrected_file_path, saf_file_path)


            logging.info(f"Starting {plastic_type} featurecounts.")

            #use FeatureCount to quantify
            
            #SAM/BAM file may be a file list so you will need to loop through them
            #one by one
            
            #There should one .saf file per plastic folder.      
            
            # Split the mappings argument into a list of individual file paths
            mapping_files = p.mappings.split(',')
            for mapping_file in mapping_files:

                sample_name = mapping_file.split('.')[0].split('/')[-1]
                mapping_file = mapping_file.strip()  # Remove any leading/trailing whitespace
                fc_input = saf_file_path
                fc_output = os.path.join(temp_folder_path, f"{plastic_type}_{os.path.basename(mapping_file)}_counts.out")
                fc_command = f"featureCounts -a {fc_input} -F SAF -o {fc_output} {mapping_file}"
                log_file = os.path.join(temp_folder_path, f"{os.path.basename(mapping_file)}_featureCounts.log")
            
                # Specify the file to capture program output
                program_output_file = os.path.join(temp_folder_path, f"{os.path.basename(mapping_file)}_featureCounts.out")
            
                with open(log_file, 'w') as f, open(program_output_file, 'w') as p_out:
                    process = subprocess.Popen(fc_command, shell=True, stdout=p_out, stderr=f)

                while process.poll() is None:
                    time.sleep(0.1)

                #Calculate statistsics on the abundance data: average, standard div, sum.

                # Initialize counts
                reads_mapped = 0
                total_reads = 0
                total_rpk = 0

                # Folder where fc_output files are stored
                fc_folder_path = temp_folder_path            
                # Iterate through all files in the directory
                for filename in os.listdir(fc_folder_path):
                    # Check if the file is a featureCounts output summary file
                    if filename.endswith("_counts.out.summary") and sample_name in filename:
                        # Open the file
                        with open(os.path.join(fc_folder_path, filename), 'r') as fc_output:

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


                    # Check if the file is a featureCounts output file
                    if filename.endswith("_counts.out") and sample_name in filename:
                        # Open the file
                        with open(os.path.join(fc_folder_path, filename), 'r') as fc_output:
                            # Skip the header lines
                            for _ in range(2):
                                next(fc_output)

                            # Iterate through each line in the file
                            for line in fc_output:

                                # Extract count and length from the line
                                fields = line.split('\t')
                                reads = int(fields[6].strip())
                                gene_length = int(fields[5].strip())

                                # Calculate RPK and add it to the list
                                rpk = reads / (gene_length / 1e3) if gene_length != 0 else 0
                                total_rpk += rpk

                            
                # Calculate proportion
                proportion = reads_mapped / total_reads if total_reads != 0 else 0
                # Format the proportion in scientific notation
                proportion = "{:.2e}".format(proportion)

                # Calculate rpkm
                rpkm = total_rpk / (total_reads / 1e6) if total_reads != 0 else 0

                # Open the TSV output file in append mode
                tsv_output_file = os.path.join(p.temps, f'{sample_name}.tsv')
                with open(tsv_output_file, 'a+') as tsv_file:
                    # Write row to TSV file
                    tsv_file.write('\t'.join([plastic_type, str(reads_mapped), str(total_reads), str(proportion), str(rpkm)]) + '\n')
        else:
            logging.warning(f"No .fasta files found in {temp_folder_path}. Skipping this folder.")
            pass
    except Exception as e:
        logging.error(f"Error running featurecounts for {plastic_type}: {e}")
