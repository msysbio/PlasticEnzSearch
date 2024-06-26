import os
import shutil

def hmm_fetch(output_dir):
    # Prompt the user for plastic motif type
    motif_type = input("Enter the plastic motif type: ")

    # Define source and destination for copying the file
    src = os.path.join(os.getcwd(), "motifs", f"{motif_type}.hmm")
    dest = output_dir

    # Copy the hmm file to the output directory
    shutil.copy2(src, dest)

    # Open and read the bitscores.txt file
    with open(os.path.join(os.getcwd(), "motifs", "bitscores.txt")) as f:
        lines = f.readlines()

    # Find the line with the correct motif and print the bitscore
    for line in lines:
        motif_name, bitscore = line.strip().split(':')
        if motif_name.lower().strip() == f"{motif_type}.hmm".lower():
            print(f"Bitscore for {motif_type}: {bitscore}")
            break
    else:
        print("Motif not found in bitscores.txt.")

