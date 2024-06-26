import sqlite3
import os

def database_fetch(plastic_type, output_dir):
    # Connect to the SQLite database
    database_path = os.path.join(os.getcwd(), "PlasticEnzymes.db")  
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    # SQL query
    query = f"""
    SELECT Main_copy.Genus, Main_copy.Species, Main_copy.Strain, Main_copy.Enzyme, Sequences.SEQUENCE
    FROM Main_copy
    JOIN Sequences ON Main_copy.Sequences_id = Sequences.Seq_Pk
    JOIN Plastic ON Main_copy.Plastic_id = Plastic.Plastic_id
    WHERE Plastic.PLASTIC = ?;
    """

    # Execute the query
    cursor.execute(query, (plastic_type,))

    # Fetch all rows from the executed SQL query
    rows = cursor.fetchall()

    # Close the connection
    conn.close()

    # Save the fetched data in FASTA format
    output_file_path = os.path.join(output_dir, f"{plastic_type}_sequences.fasta")

    with open(output_file_path, 'w') as output_file:
        for row in rows:
            # Prepare the first part of the header (row[0] and row[1]), replace None with an empty string
            first_part = ' '.join([str(item).replace(' ', ' ') if item is not None else '' for item in row[:2]])
            
            # Prepare the second part of the header (row[2] and row[3]), skip None values
            second_part = '|'.join([str(item).replace(' ', ' ') for item in row[2:4] if item is not None])
    
            # Combine the two parts of the header
            header = f"{first_part}|{second_part}"
    
            # Prepare the sequence
            if row[4] is not None:
                sequence = row[4].replace(' ', '').replace('\t', '').replace('\n', '')
            else:
                sequence = ''
                #print(row)
    
            output_file.write(f">{header}\n{sequence}\n")


def write_all_records_to_file(output_dir):
    # Connect to the SQLite database
    database_path = os.path.join(os.getcwd(), "PlasticEnzymes.db")
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    # SQL query to fetch all records
    query = "SELECT * FROM Main_copy"

    # Execute the query
    cursor.execute(query)

    # Fetch all rows from the executed SQL query
    rows = cursor.fetchall()

    # Close the connection
    conn.close()

    # Write all records to a text file
    output_file_path = os.path.join(output_dir, "all_records.txt")
    with open(output_file_path, 'w') as output_file:
        for row in rows:
            output_file.write(str(row) + '\n')

def print_example_sequence():
    # Connect to the SQLite database
    database_path = os.path.join(os.getcwd(), "PlasticEnzymes.db")  
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    # SQL query to fetch the first sequence
    query = "SELECT SEQUENCE FROM Sequences LIMIT 1"

    # Execute the query
    cursor.execute(query)

    # Fetch the first row from the executed SQL query
    row = cursor.fetchone()

    # Print the sequence
    if row is not None:
        print("Example Sequence: ", row[0])
    else:
        print("No sequence found in the database.")

    # Close the connection
    conn.close()

def print_all_tables():
    # Connect to the SQLite database
    database_path = os.path.join(os.getcwd(), "PlasticEnzymes.db")  
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    # SQL query to fetch all table names
    query = "SELECT name FROM sqlite_master WHERE type='table';"

    # Execute the query
    cursor.execute(query)

    # Fetch all rows from the executed SQL query
    rows = cursor.fetchall()

    # Print the table names
    print("Tables in the database:")
    for row in rows:
        print(row[0])

    # Close the connection
    conn.close()

def print_unique_plastic_values():
    # Connect to the SQLite database
    database_path = os.path.join(os.getcwd(), "PlasticEnzymes.db")  
    conn = sqlite3.connect(database_path)
    cursor = conn.cursor()

    # SQL query to fetch all unique PLASTIC values
    query = "SELECT DISTINCT PLASTIC FROM Plastic"

    # Execute the query
    cursor.execute(query)

    # Fetch all rows from the executed SQL query
    rows = cursor.fetchall()

    # Print the unique PLASTIC values
    print("Unique PLASTIC values:")
    for row in rows:
        print(row[0])

    # Close the connection
    conn.close()

def main():
    # Define the plastic types and the output directory
    plastic_types = ['PBAT', 'Nylon', 'n-alkanes', 'PET', 'PE', 'PHB', 'PLA', 'PCL', 'PHA', 'PHO', 'PBSA', 'MHET', 'NR']
    output_dir = '/home/jasper/Thesis/database_output'

    # Fetch the data for each plastic type
    for plastic_type in plastic_types:
        database_fetch(plastic_type, output_dir)

    '''
        # Print the output
        output_file_path = os.path.join(output_dir, f"{plastic_type}_sequences.fasta")
        with open(output_file_path, 'r') as output_file:
            print(output_file.read())

    # Write all records from the database to a file
    write_all_records_to_file(output_dir)

    # Print an example sequence
    print_example_sequence()

    # Print all table names
    print_all_tables()
    
    # Print unique PLASTIC values
    print_unique_plastic_values()
    '''

if __name__ == "__main__":
    main()


