import os
import glob
from blast_handler import blast_local
import concurrent.futures


import os
import glob
from blast_handler import blast_local
import concurrent.futures

def run_blast(dir):
    """
    Runs a local BLAST search for each fasta file in the specified directory.

    The results of each BLAST search are saved as an XML file in the same directory.

    Args:
        dir (str): The directory containing the fasta files. This should be the output directory of pathmanager but can also be used more generally.
    """
    fasta_files = glob.glob(os.path.join(dir, "*.fasta"))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for fasta_file in fasta_files:
            filename = os.path.basename(fasta_file)
            output_file = os.path.join(dir, f"{filename}.xml")

            future = executor.submit(blast_local, fasta_file, output_file)
            futures.append(future)

        # Wait for all the blast jobs to complete
        concurrent.futures.wait(futures)



