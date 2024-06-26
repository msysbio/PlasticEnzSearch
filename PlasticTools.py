import argparse
import translate_search
from abundance import annotation
from annotation import blast_search
from blast_annotation import run_blast
from quantify_hmm import quantify_hmm
from pathmanager import PathManager
import traceback
from database_operations import database_fetch
from hmm_operations import hmm_fetch
from output import create_html, remove_temps
import logging
import threading


def main(args, debug=False, blast=True):
    """Main function to run PlasticTools.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        debug (bool, optional): If True, enables debug logging. Defaults to False.
        blast (bool, optional): If True, enables blast search. Defaults to True.
    """


    # Create an instance of PathManager with the provided arguments
    p = PathManager(args)

    # Set up logging
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)
    
    # extract ORFs from the contigs and translate them into protein code using PRODIGAL
    logging.info('Extracting ORFs using PRODIGAL...')
    translate_search.run_prodigal(p)
    
    # Search extracted ORFs with HMM motifs using HMMER
    logging.info('Searching for HMM hits...')
    translate_search.run_hmmer(p)

    if blast:
        # Create a thread to run blast search in the background
        logging.info('Blasting reads...')
        blast_thread = threading.Thread(target=run_blast, args=(p.output,))
        blast_thread.start()

    # Continue with annotation
    logging.info('Calculating abundances...')
    annotation(p)

    if blast:
        # Wait for blast to finish
        blast_thread.join()

    # Create html
    logging.info('Creating html report...')
    create_html(p)

    # Remove temporary files
    logging.info('Removing temporary files...')
    remove_temps(p, debug) # debug=True to keep temporary files and only moving fasta, tsv and html files to the output directory
    



if __name__ == "__main__":
    """Entry point of the tool for the command-line interface."""

    # Print the logo
    logo = """
       ___ _           _   _        __          
      / _ | | __ _ ___| |_(_) ___  /___ __  ____
     / /_)| |/ _` / __| __| |/ __|/_\| '_ \|_  /
    / ___/| | (_| \__ | |_| | (__//__| | | |/ / 
    \/    |_|\__,_|___/\__|_|\___\__/|_| |_/___|
    """
    print(logo)

    # Create the parser
    parser = argparse.ArgumentParser(prog='PlasticTools', add_help=False)

    # Define required and optional argument groups
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    # Add required arguments
    required.add_argument('--output', required=True, help='Provide the output directory where all temporary files and outputs will be saved')
    required.add_argument('--plastic', required=True, help='Provide type of plastic searched (PLA,PET,nylon...)')
    required.add_argument('--contigs', required=True, help='Provide contigs file path')
    required.add_argument('--mappings', help='Provide path to BAM/SAM files separated by the comma')
 
    # Add optional arguments
    optional.add_argument('--b', action='store_false', help='Disable blast')
    optional.add_argument('--p', action='store_true', help='Enable multiprocessing')
    optional.add_argument('--f', action='store_true', help='Force overwrite existing files')
    optional.add_argument('--d', action='store_true', help='Enable debug logging')
 
    optional.add_argument('-v', '--version', action='version', version='%(prog)s 1.0', help="Show the version number and exit")
    optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
              help='Show this help message and exit')
 
    subparsers = parser.add_subparsers(dest='command')
    database_fetch_parser = subparsers.add_parser('database-fetch', help="Fetches data from the database")
    hmm_fetch_parser = subparsers.add_parser('hmm-fetch', help="Fetches HMM for plastics")
 
    args = parser.parse_args()
    try:
        if args.command == 'database-fetch':
            plastic_type = input("Enter the plastic type: ")
            database_fetch(plastic_type, args.output)
        elif args.command == 'hmm-fetch':
            hmm_fetch(args.output)
        else:
            main(args, debug=args.d, blast=not args.b)
       
    except Exception as e:
        with open('program_crash_data.txt', 'w') as f:
            traceback.print_exc(file=f)
            f.write(str(e))
        print(str(e))



