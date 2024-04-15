import os
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.offline as pyo
import shutil

# Function to create a pie chart and return its HTML representation
def create_plot(df):
    fig = px.bar(df, x='mapping', y='log reads mapped', color='plastic name', pattern_shape='plastic name', title='Plastic Reads Mapped', 
                 hover_data={'reads mapped': True, 'log reads mapped': False}, barmode='group')
    return fig

def create_html(p):
    # Create an empty list to store all the dataframes
    combined_dfs = []

    for filename in os.listdir(p.temps):
        if filename.endswith(".tsv"):
            name = filename.split(".tsv")[0]

            # Read the TSV file
            df = pd.read_csv(os.path.join(p.temps, filename), sep='\t')
            # Scale the values logarithmically
            df['log reads mapped'] = df['reads mapped'].apply(lambda x: np.log(x) if x != 0 else 0)
            df['mapping'] = name

            # Append the current dataframe to the list of dataframes
            combined_dfs.append(df)

    # Concatenate all the dataframes in the list into a single dataframe
    combined_df = pd.concat(combined_dfs)

    # Write the HTML string to an HTML file
    html_file = os.path.join(p.output, 'abundances.html')
    pyo.plot(create_plot(combined_df), filename=html_file, auto_open=False)


def move_files(file_extensions, source, destination):
    for root, dirs, files in os.walk(source):
        for file in files:
            if any(file.endswith(ext) for ext in file_extensions):
                src = os.path.join(root, file)
                dst = os.path.join(destination, file)
                shutil.copy(src, dst)  # Overwrites the destination file if it exists
                #os.remove(src)  # Deletes the source file

def remove_temps(p, debug=False):

    # Move tsv and fasta files from p.temps to p.output
    move_files(['.tsv', '.fasta'], p.temps, p.output)

    # Delete p.temps folder
    if not debug:
        shutil.rmtree(p.temps)

