import os
import subprocess
from Bio import SearchIO, SeqIO
import time
import threading
import utilities
import multiprocessing
from functools import partial

#use pprodigal if installed, parallelizes prodigal
def get_prodigal_command(p):
    contigs_base = os.path.basename(p.contigs).split(".")[0]
    aa_file = os.path.join(p.temps, f"{contigs_base}.faa")
    nt_file = os.path.join(p.temps, f"{contigs_base}.ffn")

    cores = utilities.run_in_parallel("pprodigal")
    if cores == 1:
        return f"prodigal -i {p.contigs} -a {aa_file} -p meta -d {nt_file}"
    else:
        return f"pprodigal -i {p.contigs} -a {aa_file} -p meta -d {nt_file}" 


def run_prodigal(p):
    if isinstance(p.plastic, str) and p.plastic != "all":
        plastic_names = p.plastic.split(',')
    elif isinstance(p.plastic, str) and p.plastic == "all":
        plastic_names = p.all_plastics
              
    if p.plastic.lower() == "all":
        # Create a sub-directory for each plastic type
        for plastic_name in plastic_names:
            temp_dir = os.path.join(os.path.abspath(p.output), "temps", plastic_name.lower())
            if os.path.exists(temp_dir):
                raise ValueError("ERROR: Temps folder already exists, please remove it and run the program again.")
            else:
                os.makedirs(temp_dir, exist_ok=True)
    else:
        temp_dir = os.path.join(os.path.abspath(p.output), "temps")
        if os.path.exists(temp_dir):
            raise ValueError("ERROR: Temps folder already exists, please remove it and run the program again.")
        else:
            os.makedirs(temp_dir, exist_ok=True)
            
    

    prodigal_command = get_prodigal_command(p)
    
    prodigal_log_file = os.path.join(p.output, "temps", "prodigal.log")
    
    task_done = threading.Event()

    with open(prodigal_log_file, 'w') as f:
        process = subprocess.Popen(prodigal_command, shell=True, stdout=f, stderr=f)

    t = threading.Thread(target=utilities.spinning_cursor_task, args=(task_done,'prodigal '))
    t.start()
    
    while process.poll() is None:
        time.sleep(0.1)
    task_done.set()
    t.join()
    
    print("\nprodigal finished running. Prodigal logs saved to {}".format(prodigal_log_file))


def run_hmmer(p):
    if isinstance(p.plastic, str) and p.plastic != "all":
        plastic_names = p.plastic.split(',')
    elif isinstance(p.plastic, str) and p.plastic == "all":
        plastic_names = p.all_plastics

    # Create a new function that has `p` already filled in
    run_hmmer_thread_p = partial(run_hmmer_thread, p=p)

    task_done = threading.Event()
    t = threading.Thread(target=utilities.spinning_cursor_task, args=(task_done,'hmmer'))
    t.start()

    # Start the tasks
    pool = multiprocessing.Pool(processes=4)
    pool.map(run_hmmer_thread_p, plastic_names)

    # Wait for the tasks to finish
    pool.close()
    pool.join()
    task_done.set()
    t.join()
    

def run_hmmer_thread(plastic_name, p):
    try:
        contigs_base = os.path.basename(p.contigs).split(".")[0]
        aa_file = os.path.join(p.output, "temps", f"{contigs_base}.faa")
        temp_dir = os.path.join(p.temps, plastic_name)

        incT = None
        with open(p.bitscores) as f:
            for line in f:
                if plastic_name in line.lower():
                    incT = float(line.split(":")[1].strip())
                    break
        if incT is None:
            raise Exception(f"\n ERROR: No bitscore found for specified plastic type: {plastic_name} in  {p.bitscores}\n ")
    
        hmm_input = os.path.join(p.motif, f"{plastic_name}.hmm")
        hmm_output = os.path.join(temp_dir, f"{contigs_base}_{plastic_name}_HMMER.out")
        hmmer_command = f"hmmsearch --tformat fasta -T {incT} --tblout {hmm_output} {hmm_input} {aa_file}"
        log_file = os.path.join(temp_dir, f"{plastic_name}_hmmsearch.log")
        
        # specify the file to capture program output
        program_output_file = os.path.join(temp_dir, f"{contigs_base}_{plastic_name}_hmmsearch.out")
        
        task_done = threading.Event()
        
        with open(log_file, 'w') as f, open(program_output_file, 'w') as p_out:
            process = subprocess.Popen(hmmer_command, shell=True, stdout=p_out, stderr=f)
        
        while process.poll() is None:
            time.sleep(0.1)

        print("\nhmmsearch finished running. Results saved to {}".format(hmm_output))
        print("hmmsearch logs saved to {}".format(log_file))
        print("hmmsearch program output saved to {}".format(program_output_file))
        
        try:
            qresults = SearchIO.read(hmm_output, "hmmer3-tab")
            hits = [result.id for result in qresults]
        except Exception as e:
            print(f"Error reading HMMER output: {e}")
            hits = []

        records = list(SeqIO.parse(aa_file, "fasta"))
        output = os.path.join(temp_dir, f"{contigs_base}_{plastic_name}_hmm_output.fasta")

        if os.path.exists(output):
            os.remove(output)

        for record in records:
            if record.id in hits and len(record) > 10:
                with open(output, "a") as f:
                    f.write(">" + str(record.id) + "\n" + str(record.seq) + "\n")
    except Exception as e:
        print(f"Error running HMMER for {plastic_name}: {e}")

    

