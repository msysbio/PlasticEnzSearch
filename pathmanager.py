import os
import sys
import utilities
import logging

def check_directory_exists(directory):
    """Check if a directory exists, create it if it doesn't.

    Args:
        directory (str): Path to the directory.

    Returns:
        str: Path to the directory.
    """
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return directory

def check_file_exists(file):
    """Check if a file exists, raise an error if it doesn't.

    Args:
        file (str): Path to the file.

    Returns:
        str: Path to the file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not os.path.isfile(file):
        raise FileNotFoundError(f"File {file} does not exist.")
    return file

def check_files_and_folders(output, contigs, mappings):
    """Check if the necessary files and folders exist.

    Args:
        output (str): Path to the output directory.
        contigs (str): Path to the contigs file.
        mappings (str): Comma-separated string of paths to the mapping files.

    Returns:
        str: Path to the temporary folder.
    """

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
    """Check if an argument is provided, exit the program if it isn't.

    Args:
        arg (str): The argument to check.
        arg_name (str): The name of the argument.
    """
    if arg is None:
        error = f"\nError: No {arg_name} provided\n"
        error += f"Use the flag --{arg_name} to provide it.\n"
        logging.error(error)
        sys.exit(1)

def fetch_motifs():
    """Fetch the hmm motifs.

    Returns:
        tuple: A tuple containing the path to the motif directory and the bitscores file.
    """
    motif_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "motifs")
    check_directory_exists(motif_dir)

    bitscores_file = os.path.join(motif_dir, 'bitscores.txt')
    check_file_exists(bitscores_file)
    
    return motif_dir, bitscores_file

def set_cores(cores=2, max_cores=True):
    """Set the number of cores to use.

    Args:
        cores (int, optional): The number of cores to use. Defaults to 2 if max_cores is False.
        max_cores (bool, optional): If True, use the maximum number of cores. Defaults to True.
    """
    if max_cores:
        utilities.set_max_cores()
    else:
        utilities.set_cores(cores)


class PathManager:
    """Class to manage the paths used in the program.

    Args:
        args (argparse.Namespace): The parsed command-line arguments.
    """
    def __init__(self, args):
        self._args = args
        self._temps = check_files_and_folders(args.output, args.contigs, args.mappings)
        self._motif, self._bitscores = fetch_motifs()
        self._all_plastics = [os.path.splitext(file)[0] for file in os.listdir(self.motif) if file.endswith('.hmm')]
        
        utilities.check_dependencies()
        set_cores(max_cores=self.multiprocessing)
        
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
    def plastic_list(self):
        plastic_arg = self._args.plastic
        if isinstance(plastic_arg, str) and plastic_arg != "all":
            plastic_names = plastic_arg.split(',')
        elif isinstance(plastic_arg, str) and plastic_arg == "all":
            plastic_names = self._all_plastics
        return plastic_names
    
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
    
    @property
    def multiprocessing(self):
        return self._args.p if hasattr(self._args, 'p') else True
    
    @property
    def force_overwrite(self):
        return self._args.f if hasattr(self._args, 'f') else True
    