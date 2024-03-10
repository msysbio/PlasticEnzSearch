import os
import sys
import utilities

def check_directory_exists(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return directory

def check_file_exists(file):
    if not os.path.isfile(file):
        raise FileNotFoundError(f"File {file} does not exist.")
    return file

def check_files_and_folders(output, contigs, mappings):
    # Check if output directory exists
    output = check_directory_exists(output)
    
    # Check if temps directory exists
    temps_folder = check_directory_exists(os.path.join(output, 'temps'))
    
    # Check if contigs file exists
    contigs = check_file_exists(contigs)
    
    # Check if bam/sam files exist
    mapping_files = mappings.split(',') if ',' in mappings else [mappings]
    for mapping_file in mapping_files:
        check_file_exists(mapping_file.strip())

    # If all checks passed
    return temps_folder

def check_arg(arg, arg_name):
    if arg is None:
        error = f"\nError: No {arg_name} provided\n"
        error += f"Use the flag --{arg_name} to provide it.\n"
        print(error)
        sys.exit(1)

def fetch_motifs():
    motif_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "motifs")
    check_directory_exists(motif_dir)

    bitscores_file = os.path.join(motif_dir, 'bitscores.txt')
    check_file_exists(bitscores_file)
    
    return motif_dir, bitscores_file


class PathManager:
    def __init__(self, args):
        self._args = args
        self._temps = check_files_and_folders(args.output, args.contigs, args.mappings)
        self._motif, self._bitscores = fetch_motifs()
        self._all_plastics = [os.path.splitext(file)[0] for file in os.listdir(self.motif) if file.endswith('.hmm')]
        
        utilities.check_dependencies()
        
    @property
    def output(self):
        return self._args.output
    
    @property
    def contigs(self):
        return self._args.contigs
    
    @property
    def plastic(self):
        return self._args.plastic
    
    @property
    def mappings(self):
        return self._args.mappings
    
    @property
    def temps(self):
        return self._temps
    
    @property
    def motif(self):
        return self._motif
    
    @property
    def bitscores(self):
        return self._bitscores
    
    @property
    def all_plastics(self):
        return self._all_plastics
    