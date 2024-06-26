import shutil
import os
import sys
import time
import multiprocessing
import logging

MAX_CORES = True
CORES = 2

def check_dependencies():
    # List of commands to check
    commands = ["hmmsearch", "featureCounts", "prodigal", "samtools"]
    
    for command in commands:
        add_to_path(command)

# Find the full path to the command
def get_path(command):
    return shutil.which(command)


def add_to_path(command):
    command_path = get_path(command)

    if command_path:
        # Get the directory containing the command
        command_dir = os.path.dirname(command_path)

        # Split the current PATH into a list of directories
        current_path = os.environ['PATH'].split(os.pathsep)

        # Only add the command's directory to the PATH if it's not already there
        if command_dir not in current_path:
            os.environ['PATH'] = command_dir + os.pathsep + os.environ['PATH']
    else:
        raise ValueError(f"ERROR: {command} command not found")


# Returns the number of logical cores on the cpu
def get_logical_cores():
    if MAX_CORES:
        return multiprocessing.cpu_count()
    else:
        return min(multiprocessing.cpu_count(), CORES)
    
def set_cores(cores):
    MAX_CORES = False
    CORES = cores

def set_max_cores():
    MAX_CORES = True
    

# If the command used to run in parallel is in the PATH, return the number of logical cores, otherwise return 1
def run_in_parallel(command):
    try:
        command_path = get_path(command)
        if command_path is not None:
            logical_cores = get_logical_cores()
            if isinstance(logical_cores, int):
                return logical_cores
        return 1
    except Exception:
        return 1



def spinning_cursor():
    while True:
        for cursor in '|/-\\':
            yield cursor

def spinning_cursor_task(task_done,program):
    spinner = spinning_cursor()
    logging.info(f"Starting {program} in the background. This may take a while...")
    while not task_done.is_set():
        sys.stdout.write(next(spinner))  # write the next character
        sys.stdout.flush()                # flush stdout buffer (actual character display)
        sys.stdout.write('\b')            # erase the last written char
        time.sleep(0.1)
    logging.info(f"{program} finished running.")