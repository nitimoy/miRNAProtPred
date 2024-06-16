from flask import Flask, render_template, request, send_file, jsonify
import re
import pandas as pd
import io
from werkzeug.utils import secure_filename
import psycopg2
import time
import random
import string

app = Flask(__name__)

# Database credentials
db_params = {
    'host': 'localhost',
    'database': 'db_name',
    'user': 'db_user',
    'password': 'db_password'
}

# Establish a connection
conn = psycopg2.connect(**db_params)
    
# Genetic code dictionary
genetic_code = {
    "A": "GCU", "R": "CGU", "N": "AAU", "D": "GAU",
    "C": "UGU", "Q": "CAA", "E": "GAA", "G": "GGU",
    "H": "CAU", "I": "AUU", "L": "UUA", "K": "AAA",
    "M": "AUG", "F": "UUU", "P": "CCU", "S": "UCU",
    "T": "ACU", "W": "UGG", "Y": "UAU", "V": "GUU",
    "*": "UAA"
}

#boyer_moore_search algo
def preprocess_bad_character_table(pattern):
    table = {}
    for i in range(len(pattern)):
        table[pattern[i]] = len(pattern) - i - 1
    return table

def boyer_moore_search(text, pattern):
    m = len(pattern)
    n = len(text)
    if m > n:
        return []

    table = preprocess_bad_character_table(pattern)
    occurrences = []

    i = m - 1
    while i < n:
        j = m - 1
        while j >= 0 and text[i] == pattern[j]:
            i -= 1
            j -= 1

        if j == -1:
            occurrences.append(i + 1)

        bad_char_skip = table.get(text[i], m)
        i += max(m - j, bad_char_skip)

    return occurrences

#All the transcription
def transcription_translation(sequence):
    if all(nucleotide in "ATGC" for nucleotide in sequence):
        input_format = "DNA"
    elif all(nucleotide in "AUGC" for nucleotide in sequence):
        input_format = "RNA"
    else:
        input_format = "Protein"

    if input_format == "Protein":
        # Transcribe protein to RNA
        rna_sequence = protein_to_rna(sequence)
        # Translate RNA to DNA
        dna_sequence = rna_to_dna(rna_sequence)
    elif input_format == "RNA":
        # Translate RNA to DNA
        dna_sequence = rna_to_dna(sequence)
    elif input_format == "DNA":
        # Directly use DNA sequence
        dna_sequence = sequence
    else:
        raise ValueError("Invalid input format")

    return dna_sequence

def protein_to_rna(protein_sequence):
    rna_sequence = ""
    for amino_acid in protein_sequence:
        if amino_acid in genetic_code:
            codon = genetic_code[amino_acid]
            rna_sequence += codon
    return rna_sequence

def rna_to_dna(rna_sequence):
    dna_sequence = ""
    for nucleotide in rna_sequence:
        if nucleotide == "U":
            dna_sequence += "T"
        else:
            dna_sequence += nucleotide
    return dna_sequence

def rows(conn):
    cursor = conn.cursor()
    sql_query = "SELECT * FROM mirdb"
    cursor.execute(sql_query)
    rows = cursor.fetchall()
    cursor.close()
    return rows

def search_miRNA_id_in_db(miRNA_id, connection):
    cursor = connection.cursor()
    sql_query = "SELECT * FROM mirdb WHERE mir_id = %s"
    cursor.execute(sql_query, (miRNA_id,))
    result = cursor.fetchall()
    cursor.close()
    return result

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')

@app.route('/about_us', methods=['GET', 'POST'])
def about_us():
    return render_template('about_us.html')

@app.route('/contact_us', methods=['GET', 'POST'])
def contact_us():
    return render_template('contact_us.html')

@app.route('/utr_prime', methods=['GET', 'POST'])
def utr_prime():
    if request.method == 'POST':
        sequence = request.form['sequence']
        job_id = request.form['job_id']  # Get the user-provided Job Id
        sequence = re.sub(r'\s+', '', sequence).upper()
        dna_sequence = transcription_translation(sequence)
        dna_sequence = dna_sequence.lower()
        
        # Extract the subset of the DNA sequence
        subset_dna_sequence = dna_sequence

        # Convert fetched rows to a DataFrame
        column_labels = ["r_id", "Description", "HumanmiRNAID","Accession", "Sequence", "seed_1","seed_2", "seed_3"]  
        df = pd.DataFrame(rows(conn), columns=column_labels)

        # Identify the correct # Identify the correct column labels
        column_labels = ["seed_1","seed_2", "seed_3"]

        # Extract data from the identified column labels
        all_target_sequences = []
        for column_label in column_labels:
            if column_label in df.columns:
                target_sequences_from_column = df[column_label].dropna().tolist()
                all_target_sequences.extend(target_sequences_from_column)

        # Find occurrences of each target sequence in the given DNA sequence using Boyer-Moore
        target_sequence_matches = []
        for target_sequence in all_target_sequences:
            occurrences = boyer_moore_search(subset_dna_sequence, target_sequence)
            if occurrences:
                target_sequence_matches.append((target_sequence, occurrences))

        # Convert fetched rows to a DataFrame for generating new Excel File
        column_labels = ["r_id", "Description", "Human_miRNA_ID","Accession_ID", "Sequence", "seed_1","seed_2", "seed_3"]  
        df = pd.DataFrame(rows(conn), columns=column_labels)

        # Create a new DataFrame to store matching rows with two additional columns
        matching_rows_df = pd.DataFrame(columns=df.columns.tolist() + ["Seed", "Position"])

        # Set to keep track of already added rows
        added_rows = set()

        # Add matching rows to the new DataFrame with "Seed" and "Position" columns
        for target_sequence, occurrences in target_sequence_matches:
            matching_rows = df[df.isin([target_sequence]).any(axis=1)]
            matching_rows = matching_rows[~matching_rows.index.isin(added_rows)]  # Exclude already added rows
            
            if not matching_rows.empty:
                for occurrence in occurrences:
                    # Create a new row for each occurrence
                    occurrence_row = matching_rows.copy()
                    occurrence_row["Seed"] = target_sequence
                    occurrence_row["Position"] = occurrence + 1
                    
                    matching_rows_df = pd.concat([matching_rows_df, occurrence_row])
                    
                added_rows.add(str(target_sequence_matches))
                added_rows.update(matching_rows.index)
        
        #For displaying in the webpage in ascending order
        display_rows = matching_rows_df.sort_values(by="Position", ascending=True)
        
        if display_rows.empty:
            no_results_message = "No matching rows found based on your input."
            return render_template('utr_prime.html', no_results_message=no_results_message)

        
        # Remove the first column before saving to Excel
        matching_rows_df = matching_rows_df.iloc[:, 1:]  # This line removes the first column
          # Exclude "seed1," "seed2," and "seed3" columns
        matching_rows_df = matching_rows_df.drop(["seed_1", "seed_2", "seed_3"], axis=1) 

        # Sort the DataFrame based on the "Position" column in ascending order
        matching_rows_df = matching_rows_df.sort_values(by="Position", ascending=True)
        
        # Use the Job Id as the filename if provided, or generate a unique filename
        if job_id:
            filename = f"downloads/{job_id}.xlsx"
        else:
            # Generate a unique filename based on timestamp and random string
            timestamp = int(time.time())
            random_string = ''.join(random.choices(string.ascii_letters, k=5))
            filename = f"downloads/matching_rows_{timestamp}_{random_string}.xlsx"

        
        matching_rows_df.to_excel(filename, index=False)

        return render_template('utr_prime.html', display_rows=display_rows.to_dict(orient='records'), excel_path=filename, job_id=job_id)
    return render_template('utr_prime.html', display_rows=None)

@app.route('/sequence_finder', methods=['GET', 'POST'])
def sequence_finder():
    if request.method == 'POST':
        sequence = request.form['sequence']
        job_id = request.form['job_id']  # Get the user-provided Job Id
        sequence = re.sub(r'\s+', '', sequence).upper()
        dna_sequence = transcription_translation(sequence)
        dna_sequence = dna_sequence.lower()
        
        # Extract the subset of the DNA sequence
        subset_dna_sequence = dna_sequence

        # Convert fetched rows to a DataFrame
        column_labels = ["r_id", "Description", "Human miRNA ID","Accession ID", "Sequence", "seed_1","seed_2", "seed_3"] 
        df = pd.DataFrame(rows(conn), columns=column_labels)

        # Identify the correct # Identify the correct column labels
        column_labels = ["seed_1","seed_2", "seed_3"]

        # Extract data from the identified column labels
        all_target_sequences = []
        for column_label in column_labels:
            if column_label in df.columns:
                target_sequences_from_column = df[column_label].dropna().tolist()
                all_target_sequences.extend(target_sequences_from_column)

        # Find occurrences of each target sequence in the given DNA sequence using Boyer-Moore
        target_sequence_matches = []
        for target_sequence in all_target_sequences:
            occurrences = boyer_moore_search(subset_dna_sequence, target_sequence)
            if occurrences:
                target_sequence_matches.append((target_sequence, occurrences))

        # Convert fetched rows to a DataFrame for generating new Excel File
        column_labels = ["r_id", "Description", "Human_miRNA_ID","Accession_ID", "Sequence", "seed_1","seed_2", "seed_3"]  
        df = pd.DataFrame(rows(conn), columns=column_labels)

        # Create a new DataFrame to store matching rows with two additional columns
        matching_rows_df = pd.DataFrame(columns=df.columns.tolist() + ["Seed", "Position"])

        # Set to keep track of already added rows
        added_rows = set()

        # Add matching rows to the new DataFrame with "Seed" and "Position" columns
        for target_sequence, occurrences in target_sequence_matches:
            matching_rows = df[df.isin([target_sequence]).any(axis=1)]
            matching_rows = matching_rows[~matching_rows.index.isin(added_rows)]  # Exclude already added rows
            
            if not matching_rows.empty:
                for occurrence in occurrences:
                    # Create a new row for each occurrence
                    occurrence_row = matching_rows.copy()
                    occurrence_row["Seed"] = target_sequence
                    occurrence_row["Position"] = occurrence + 1
                    
                    matching_rows_df = pd.concat([matching_rows_df, occurrence_row])
                    
                added_rows.add(str(target_sequence_matches))
                added_rows.update(matching_rows.index)

        display_rows = matching_rows_df

        if display_rows.empty:
            no_results_message = "No matching rows found based on your input."
            return render_template('sequence_finder.html', no_results_message=no_results_message)
        elif (len(matching_rows_df.head(int(len(matching_rows_df))))) < 10:
            display_rows = matching_rows_df
        # Select first 20% of rows
        else: 
            display_rows = matching_rows_df.head(int(len(matching_rows_df) * 0.2))
        
        # Remove the first column before saving to Excel
        matching_rows_df = matching_rows_df.iloc[:, 1:]  # This line removes the first column

        # Sort the DataFrame based on the "Occurrences" column in ascending order
        matching_rows_df = matching_rows_df.sort_values(by="Position", ascending=True)

          # Exclude "seed1," "seed2," and "seed3"
        matching_rows_df = matching_rows_df.drop(["seed_1", "seed_2", "seed_3"], axis=1) 
        
        # Use the Job Id as the filename if provided, or generate a unique filename
        if job_id:
            filename = f"downloads/{job_id}.xlsx"
        else:
            # Generate a unique filename based on timestamp and random string
            timestamp = int(time.time())
            random_string = ''.join(random.choices(string.ascii_letters, k=5))
            filename = f"downloads/matching_rows_{timestamp}_{random_string}.xlsx"

        
        matching_rows_df.to_excel(filename, index=False)

        return render_template('sequence_finder.html', display_rows=display_rows.to_dict(orient='records'), excel_path=filename, job_id=job_id)
    return render_template('sequence_finder.html', display_rows=None)


@app.route('/validator', methods=['GET', 'POST'])
def validator():
    if request.method == 'POST':
        sequence = request.form['sequence']
        miRNA_ids_input = request.form['miRNA_ids']
        sequence = re.sub(r'\s+', '', sequence).upper()
        dna_sequence = transcription_translation(sequence)
        dna_sequence = dna_sequence.lower()
        # Split miRNA IDs by commas and ensure "hsa-" is added if not already present
        miRNA_ids = []
        for id in miRNA_ids_input.split(','):
            id = id.strip()
            if not id.startswith('hsa-'):
                id = 'hsa-' + id
            miRNA_ids.append(id)

        results = []
        for miRNA_id in miRNA_ids:
            search_result = search_miRNA_id_in_db(miRNA_id, conn)
            
            if search_result:
                miRNA_info = {
                    'id': miRNA_id,
                    'found': True,
                    'seed_match': False
                }
                seeds = search_result[0][-3:]
                found_match = any(seed in dna_sequence for seed in seeds)
                if found_match:
                    miRNA_info['seed_match'] = True
                
                results.append(miRNA_info)
            else:
                results.append({'id': miRNA_id, 'found': False})
        
        return render_template('validator.html', results=results)
    
    return render_template('validator.html', results=[])

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
