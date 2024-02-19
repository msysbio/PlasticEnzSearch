import argparse
import sys
import traceback
import os
from database_operations import database_fetch
from hmm_operations import hmm_fetch
from translate_search import translate_search
from quantify_hmm import quantify_hmm

logo = """
   ___ _           _   _        __          
  / _ | | __ _ ___| |_(_) ___  /___ __  ____
 / /_)| |/ _` / __| __| |/ __|/_\| '_ \|_  /
/ ___/| | (_| \__ | |_| | (__//__| | | |/ / 
\/    |_|\__,_|___/\__|_|\___\__/|_| |_/___|
                                            
"""
print(logo)

# Create the parser
parser = argparse.ArgumentParser(prog='PlasticEnz Tool', add_help=False)

required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument('--output', required=True, help='Provide the output directory where all temporary files and outputs will be saved')
required.add_argument('--plastic', required=True, help='Provide type of plastic searched (PLA,PET,nylon...)')
required.add_argument('--contigs', required=True, help='Provide contigs file path')
required.add_argument('--mappings', help='Provide path to BAM/SAM files separated by the comma')

# Add the help argument back in
optional.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                      help='Show this help message and exit')

# Create subparsers for the "database-fetch" and "hmm-fetch" commands
subparsers = parser.add_subparsers(dest='command')
database_fetch_parser = subparsers.add_parser('database-fetch', help="Fetches data from the database")
hmm_fetch_parser = subparsers.add_parser('hmm-fetch', help="Fetches HMM for plastics")

# Parse the command-line arguments
args = parser.parse_args()

try:
    # Check which command was specified and call the respective function
    if args.command == 'database-fetch':
        plastic_type = input("Enter the plastic type: ")
        database_fetch(plastic_type, args.output)
    elif args.command == 'hmm-fetch':
        hmm_fetch(args.output)
    else:  # if no command is given, run the main function
        print('Running translate_search...')
        translate_search(args)
        print('Running quantify_hmm...')
        quantify_hmm(args)
        #print('Running BLAST searches...')
        #blast_search_in_directory(os.path.join(args.output, 'temps'))
except Exception as e:
    with open('program_crash_data.txt', 'w') as f:
        traceback.print_exc(file=f)
        f.write(str(e))
    print(str(e))


