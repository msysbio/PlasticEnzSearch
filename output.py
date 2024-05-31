import os
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.offline as pyo
from plotly.offline import plot
from plotly.graph_objs import Layout
from plotly.subplots import make_subplots
import shutil

# Function to create a bar chart and return its HTML representation
def create_plot(df, y, hover_y = None):
    if hover_y is None:
        return px.bar(df, x='sample', y=y, color='plastic name', pattern_shape='plastic name', title=y, barmode='group')
    else:
        return px.bar(df, x='sample', y=hover_y, color='plastic name', pattern_shape='plastic name', title='Plastic Reads Mapped', 
                 hover_data={y: True, hover_y: False}, barmode='group')
    


def create_fig(df):
        
    # Initialize a dictionary to keep track of trace names
    trace_names = {}

    # Create a list of figures
    figs = [
        create_plot(df, 'reads mapped', 'log reads mapped'),
        create_plot(df, 'proportion'),
        create_plot(df, 'rpkm')
    ]
    n = len(figs)  # Number of figures
    fig = make_subplots(rows=n, cols=1)

    # Add traces to the subplot
    for i, fig_i in enumerate(figs):
        for trace in fig_i.data:
            # If the trace name is already in the legend, hide it
            if trace.name in trace_names:
                trace.showlegend = False
            else:
                trace_names[trace.name] = True
            fig.add_trace(trace, row=i+1, col=1)

    # Set the labels for the x- and y-axis for each subplot
    ylabels = ['Log Reads Mapped', 'Proportion', 'RPKM']
    xlabels = ['Sample', 'Sample', 'Sample']
    for i in range(n):
        fig.update_yaxes(title_text=ylabels[i], row=i+1, col=1)
        fig.update_xaxes(title_text=xlabels[i], row=i+1, col=1) 

    # Set a single legend for the entire subplot
    fig.update_layout(showlegend=True)

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
            df['sample'] = name

            # Append the current dataframe to the list of dataframes
            combined_dfs.append(df)

    # Concatenate all the dataframes in the list into a single dataframe
    combined_df = pd.concat(combined_dfs)

    fig = create_fig(combined_df)

    # Generate the Plotly figure's div element
    plot_div = plot(fig, include_plotlyjs=False, output_type='div')

    # Create the HTML content with the title, text, and Plotly figure
    html_content = f"""
    <html>
    <head>
    <title>PlasticTools</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>
    <h1>Enzyme Abundance</h1>
    <p>These figures show the abundance of potential plastic degrading enzymes in the samples given to the <i>PlasticTools</i> program.</p>
    {plot_div}
    </body>
    </html>
    """

    # Save the HTML content to a file
    html_file = os.path.join(p.output, 'abundances.html')
    with open(html_file, 'w') as f:
        f.write(html_content)
        



def move_files(file_extensions, source, destination):
    for root, dirs, files in os.walk(source):
        for file in files:
            if any(file.endswith(ext) for ext in file_extensions):
                src = os.path.join(root, file)
                dst = os.path.join(destination, file)
                shutil.copy(src, dst)  # Overwrites the destination file if it exists
                #os.remove(src)  # Deletes the source file(should be removed anyway by remove_temps function)

def remove_temps(p, debug=False):

    # Move tsv and fasta files from p.temps to p.output
    move_files(['.tsv', '.fasta'], p.temps, p.output)

    # Delete p.temps folder
    if not debug:
        shutil.rmtree(p.temps)

