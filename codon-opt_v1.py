
import pandas as pd
import openpyxl
import numpy as np
import os
from collections import Counter
import math
import random
import re

def reverse_complement(sequence):
    # Define complement mapping
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    
    # Reverse the sequence and map each base to its complement
    reversed_complemented_sequence = "".join(complement[base] for base in reversed(sequence))
    
    return reversed_complemented_sequence
def calculate_gc_content(dna_sequence):
    """
    Calculates the GC content of a DNA sequence.

    Parameters:
    dna_sequence (str): The DNA sequence to analyze.

    Returns:
    float: The GC content as a percentage of the sequence.
    """
    # Count occurrences of 'G' and 'C' in the DNA sequence
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    
    # Calculate GC content as a percentage
    gc_content = (gc_count / len(dna_sequence)) * 100
    return gc_content


# Define valid characters
valid_bases = ["G", "A", "T", "C", "*"]
valid_AA = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
# Define stop codons
stop_codons = ["TAA", "TAG", "TGA"]
# Define a codon-to-amino-acid dictionary for translation
codon_table = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}

# Define codon usage profiles for optimization
codon_usage_profiles = {
    1: "Human",
    2: "EColi",
    3: "Yeast",
    4: "Mouse",
    5: "Insect"
}


# Define file paths for each codon usage profile (replace with your actual file paths)
# Get the directory where the script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
print("Directory")
print(script_dir)
codon_usage_files = {
    1: os.path.join(script_dir, "human.xlsx"),
    2: os.path.join(script_dir, "ecoli.xlsx"),
    3: os.path.join(script_dir, "yeast.xlsx"),
    4: os.path.join(script_dir, "mouse.xlsx"),
    5: os.path.join(script_dir, "insect.xlsx")
}


#Start the program

while True:
    # Introduction
    print("Basti's Codon Optimization Tool.")
    print("Please provide only the DNA or Amino Acid Sequence in a text file.")

    # Define the path
    path = input("Enter path to Sequence.txt (e.g., C:/...): ")

    # Read the sequence file
    try:
        with open(path, 'r') as file:
            sequence = file.read().strip().upper()  # Read, strip, and convert to uppercase
        print("Sequence successfully read and converted to uppercase!")
        print(f"Length of the sequence: {len(sequence)}")
    except FileNotFoundError:
        print("Error: File not found. Please check the file path.")
        continue  # Restart the loop if the file is not found

    # Check if a DNA or protein sequence was provided
    sequence_type = None
    while sequence_type is None:
        sequence_type = input("Is this a DNA or an Amino Acid (AA) Sequence? Enter 'DNA' or 'AA': ").strip().upper()

        if sequence_type == "DNA":
            is_dna = True
            print("DNA sequence confirmed.")
        elif sequence_type == "AA":
            is_dna = False
            print("Amino Acid sequence confirmed.")
        else:
            print("Invalid input! Please enter 'DNA' or 'AA'.")
            sequence_type = None  # Reset to loop again for valid input

    

    # Check for invalid characters
    if is_dna:
        invalid_chars = [base for base in sequence if base not in valid_bases]
        if invalid_chars:
            print(f"Invalid characters in DNA sequence: {set(invalid_chars)}")
            restart = input("Invalid sequence detected. Would you like to try again? (Y/N): ").strip().upper()
            if restart == 'Y':
                continue  # Restart the program
            else:
                print("Exiting program.")
                break  # Exit the program
        else:
            print("DNA sequence is valid.")
    else:
        invalid_chars = [aa for aa in sequence if aa not in valid_AA]
        if invalid_chars:
            print(f"Invalid characters in Amino Acid sequence: {set(invalid_chars)}")
            restart = input("Invalid sequence detected. Would you like to try again? (Y/N): ").strip().upper()
            if restart == 'Y':
                continue  # Restart the program
            else:
                print("Exiting program.")
                break  # Exit the program
        else:
            print("Amino Acid sequence is valid.")
    
    # If sequence is valid, break out of the main loop to proceed further in the code
    break


# Translation to protein sequence

# IF DNA 
if is_dna:
    
            
    while True:
        # Find Start Codon
        start_positions = []
        start_codons = ["ATG","CAT"]
        position = 0
        while position <= len(sequence)-3:
            codon = sequence[position:position+3]
            if codon in start_codons:
                start_positions.append(position)
            position += 1
        if not start_positions:
            append_start=input('Append a Start in 5 (ATG) or 3 (CAT) prime end? 5/3').strip().upper()
            if append_start == "5":
                sequence = "ATG" + sequence
                start_positions.append(0)     
            elif append_start == "3":
                sequence = sequence + "CAT"
                start_positions.append(len(sequence)-1)
        print("Those potential start codons (including reverse complement) were found:")
        print(start_positions)
        ORF_start = None
        while ORF_start is None:
            try:
                ORF_start = int(input("what's your start position?").strip())
            except ValueError:
                print("Invalid entry")
                continue
            if ORF_start not in start_positions:
                print("Invalid")
                ORF_start = None
            else:
                print("Accepted")
                print(sequence[ORF_start:ORF_start+3])
        
        #potentially reverse complement if the start is wrong directed
        if sequence[ORF_start:ORF_start+3] == "CAT":
            print("Reverse start detected. Trying to reverse complement the sequence!")
            sequence = reverse_complement(sequence)
            print("Sequence has been reverse complemented. New Sequence:")
            print(sequence)
            #restart the finding start
            continue
        else:
            print("No reverse complement needed. Continue with translation.")
            #break out of the finding start loop
            break    
# Translate to Amino Acids
    
    # Translate DNA sequence to protein sequence
    protein_sequence = []
    position = ORF_start

    while position <= len(sequence) - 3:
        codon = sequence[position:position+3]
        amino_acid = codon_table.get(codon, "")
        
        if amino_acid == "*":  # Stop codon
            print(f"Stop codon encountered at position {position}. Ending translation.")
            protein_sequence.append(amino_acid)
            break
        
        if amino_acid:  # Valid codon found
            protein_sequence.append(amino_acid)
        else:
            print(f"Invalid codon '{codon}' at position {position}. Skipping.")

        position += 3  # Move to the next codon

    # Join the amino acids into a protein sequence string
    protein_sequence = "".join(protein_sequence)     
           
#handle protein sequence
elif not is_dna:
    # Find Methionin
    start_positions = []
    start_codons = ["M"]
    position = 0
    while position <= len(sequence)-1:
        aa = sequence[position]
        if aa in start_codons:
            start_positions.append(position)
        position += 1
    if not start_positions:
        append_start=input('No potential starts found. Append a Methionin in the beginning? y/n').strip().upper()
        if append_start == "Y":
            sequence = "M" + sequence
        start_positions.append(0)  
                
    print("Those potential start codons (Methionins) were found:")
    print(start_positions)
    ORF_start = None
    
    while ORF_start is None:
        try:
            ORF_start = int(input("what's your start position?").strip())
        except ValueError:
            print("Invalid entry")
            continue
        if ORF_start not in start_positions:
            print("Invalid")
            ORF_start = None
        else:
            print("Accepted")
            
    protein_sequence = []
    position = ORF_start
    aa = sequence[position]
    while position <= len(sequence)-1 and aa != "*":
        aa = sequence[position]
        protein_sequence.append(aa)
        position += 1  # Move to the next aa

    # Join the amino acids into a protein sequence string
    protein_sequence = "".join(protein_sequence)
#check for stop codon
if protein_sequence[-1] != "*":
    append_stop = input("No Stop found, should a stop codon be appended? y/n").strip().upper()
    if append_stop == "Y":
        protein_sequence= protein_sequence+"*"
        
    
print("Generated protein sequence:")
print(protein_sequence)
# Optimization Parameters

# Prompt user to select a codon usage profile
codon_usage = None
while codon_usage is None:
    try:
        # User input
        codon_usage = int(input("Which codon usage do you want to apply?\n 1 - Human\n 2 - EColi\n 3 - Yeast\n 4 - Mouse\n 5 - Insect\n").strip())
        
    # Check if the input is valid
        if codon_usage in codon_usage_profiles:
            print(f"Selected codon usage profile: {codon_usage_profiles[codon_usage]}")
            
            # Load the appropriate codon usage file as a DataFrame
            file_path = codon_usage_files[codon_usage]
            print(file_path)
            codon_df = pd.read_excel(file_path)
            print(f"Loaded {codon_usage_profiles[codon_usage]} codon usage data.")
            break
        else:
            print("Invalid input. Please enter a number between 1 and 5.")
            codon_usage = None  # Reset to loop again
    except ValueError:
        print("Invalid input. Please enter a valid number.")


#Choose algorithm
algorithm = None
while algorithm is None:
    print("Available algorithms:\n")
    print("Optimization:\n 1 - Most Frequent\n 2 - Probability Frequency Distribution\n 3 - Enforced Frequency Distribution\n")
    print("Deoptimization:\n 4 - Least Frequent\n 5 - inverted Probability Frequency Distribution\n 6 - inverted Enforced Frequency Distribution\n")
    try:
        algorithm = int(input("Please choose an optimization algorithm:"))
    except ValueError:
        print("please provide a valid option")
        continue
    if algorithm != 1 and algorithm != 2 and algorithm != 3 and algorithm != 4 and algorithm != 5 and algorithm != 6:
        print("please provide a valid option")
        algorithm = None
    else:
        print(f"Algorithm {algorithm} selected.")
        break


# most frequent algorithm
if algorithm == 1:
    # Step 1: Find the most frequent codon for each amino acid
    most_frequent_codons = (
        codon_df.loc[codon_df.groupby('Amino')['Frequency'].idxmax()]
        .set_index('Amino')['Codon']
        .to_dict()
    )

    # Step 2: Translate the protein sequence using the most frequent codons
    optimized_dna_sequence = ''.join(most_frequent_codons[aa] for aa in protein_sequence)
# probability algorithm
elif algorithm == 2:
    #step 1: Generate a Probabalistic pool to draw the codons for every amino acid.
    # Step 1: Create codon pools based on frequency
    codon_pools = {}

    for amino_acid, group in codon_df.groupby('Amino'):
        pool = []
        for _, row in group.iterrows():
            # Scale frequency to a count for a 100-size pool
            count = int(row['Frequency'] * 100)
            pool.extend([row['Codon']] * count)
        codon_pools[amino_acid] = pool

    # Step 2: Translate protein sequence to DNA using probabilistic codon selection
    optimized_dna_sequence = ''.join(np.random.choice(codon_pools[aa]) for aa in protein_sequence)
elif algorithm == 3:
    # Step 1: Count amino acids in the protein sequence
    amino_acid_counts = Counter(protein_sequence)
    print(amino_acid_counts)
    # Step 2: Create a fixed codon pool for each amino acid based on frequency, slightly overfilled to ensure enough codons
    def create_codon_pool(amino_acid, count):
        pool = []
        amino_df = codon_df[codon_df['Amino'] == amino_acid]

        # Use math.ceil to ensure the required number of codons
        required_codons = math.ceil(count)
        for _, row in amino_df.iterrows():
            codon_count = math.ceil(row['Frequency'] * required_codons)
            pool.extend([row['Codon']] * codon_count)
            
        np.random.shuffle(pool)
        return pool

    # Initialize codon pools
    codon_pools = {amino_acid: create_codon_pool(amino_acid, count) for amino_acid, count in amino_acid_counts.items()}

    # Step 3: Translate protein sequence by drawing unique codons from pools
    optimized_dna_sequence = []

    for aa in protein_sequence:
        # Refill the pool if it’s empty by recreating it with the original count
        if not codon_pools[aa]:  
            codon_pools[aa] = create_codon_pool(aa, amino_acid_counts[aa])  

        # Draw codon from pool
        if codon_pools[aa]:
            codon = codon_pools[aa].pop()
        else:
            # Fallback to the most frequent codon if pool is unexpectedly empty
            codon = codon_df[codon_df['Amino'] == aa].sort_values(by='Frequency', ascending=False).iloc[0]['Codon']
        
        optimized_dna_sequence.append(codon)

    # Join list into final DNA sequence
    optimized_dna_sequence = ''.join(optimized_dna_sequence)
if algorithm == 4:
    # Algorithm 4: Least Frequent Codon
    optimized_dna_sequence = []
    for aa in protein_sequence:
        # Select the codon with the lowest frequency for the amino acid
        least_frequent_codon = codon_df[codon_df['Amino'] == aa].sort_values(by='Frequency').iloc[0]['Codon']
        optimized_dna_sequence.append(least_frequent_codon)
    # Join list into final DNA sequence
    optimized_dna_sequence = ''.join(optimized_dna_sequence)

elif algorithm == 5:
    # Algorithm 5: Inverted Probability Frequency Distribution
    # Invert frequencies in codon_df
    codon_df['Inverted_Frequency'] = 1 - codon_df['Frequency']
    
    # Replace any NaN values with 0
    codon_df['Inverted_Frequency'].fillna(0, inplace=True)
    # For codons with a frequency of 1, set Inverted_Frequency to 1 (highest preference)
    codon_df.loc[codon_df['Frequency'] == 1, 'Inverted_Frequency'] = 1
    codon_df['Inverted_Frequency'] /= codon_df.groupby('Amino')['Inverted_Frequency'].transform('sum')  # Normalize within each amino acid        
    optimized_dna_sequence = []
    for aa in protein_sequence:
        aa_df = codon_df[codon_df['Amino'] == aa]
        codons = aa_df['Codon'].values
        probabilities = aa_df['Inverted_Frequency'].values
        
        # Randomly choose codon based on inverted frequency distribution
        chosen_codon = np.random.choice(codons, p=probabilities)
        optimized_dna_sequence.append(chosen_codon)
    # Join list into final DNA sequence
    optimized_dna_sequence = ''.join(optimized_dna_sequence)

elif algorithm == 6:
    # Algorithm 6: Inverted Enforced Frequency Distribution
    # Step 1: Count amino acids in the protein sequence
    amino_acid_counts = Counter(protein_sequence)
    # Invert frequencies in codon_df
    codon_df['Inverted_Frequency'] = 1 - codon_df['Frequency']
           
    # Replace any NaN values with 0
    codon_df['Inverted_Frequency'].fillna(0, inplace=True)
    # For codons with a frequency of 1, set Inverted_Frequency to 1 (highest preference)
    codon_df.loc[codon_df['Frequency'] == 1, 'Inverted_Frequency'] = 1      
    codon_df['Inverted_Frequency'] /= codon_df.groupby('Amino')['Inverted_Frequency'].transform('sum')  # Normalize within each amino acid 
    # Step 2: Create a codon pool for each amino acid based on inverted frequency, using math.ceil to ensure enough codons
    def create_inverted_codon_pool(amino_acid, count):
        pool = []
        amino_df = codon_df[codon_df['Amino'] == amino_acid]
        
        # Use math.ceil to ensure the required number of codons
        required_codons = math.ceil(count)
        for _, row in amino_df.iterrows():
            codon_count = math.ceil(row['Inverted_Frequency'] * required_codons)
            pool.extend([row['Codon']] * codon_count)
            
        np.random.shuffle(pool)
        return pool
    
    # Initialize codon pools
    codon_pools = {amino_acid: create_inverted_codon_pool(amino_acid, count) for amino_acid, count in amino_acid_counts.items()}
    
    # Step 3: Translate protein sequence by drawing unique codons from pools
    optimized_dna_sequence = []
    
    for aa in protein_sequence:
        # Refill the pool if it’s empty by recreating it with the original count
        if not codon_pools[aa]:  
            codon_pools[aa] = create_inverted_codon_pool(aa, amino_acid_counts[aa])  
        
        # Draw codon from pool
        if codon_pools[aa]:
            codon = codon_pools[aa].pop()
        else:
            # Fallback to the least frequent codon if pool is unexpectedly empty
            codon = codon_df[codon_df['Amino'] == aa].sort_values(by='Frequency').iloc[0]['Codon']
        
        optimized_dna_sequence.append(codon)   
    # Join list into final DNA sequence
    optimized_dna_sequence = ''.join(optimized_dna_sequence) 
    
print("Initial optimized/deoptimized DNA Sequence:", optimized_dna_sequence)
# Calculate GC content for the optimized DNA sequence
gc_content = calculate_gc_content(optimized_dna_sequence)
print(f"GC Content of optimized DNA sequence: {gc_content:.2f}%")


def avoid_codon_duplets(dna_sequence, codon_df):
    """
    Checks for codon duplets within a sliding window and mutates one of the duplicates if found.

    Parameters:
    dna_sequence (str): The DNA sequence to modify.
    codon_df (DataFrame): Dataframe containing codon frequency data.

    Returns:
    str: Modified DNA sequence with reduced codon duplets.
    """
    # Convert sequence to a list of codons
    codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]
    
    # Sliding window approach to check for duplicates within a 6-base window (2 codons)
    for i in range(len(codons) - 1):
        if codons[i] == codons[i + 1]:  # Found duplicate codons in the window
            print(f"Duplicate codon found: {codons[i]} at positions {i} and {i + 1}")
            
            # Decide randomly which of the duplicates to mutate (flip a coin)
            to_mutate = i if np.random.rand() < 0.5 else i + 1
            
            # Identify the amino acid for this codon
            amino_acid = codon_df.loc[codon_df['Codon'] == codons[to_mutate], 'Amino'].values[0]
            
            # Get alternative codons for this amino acid sorted by frequency
            alternatives = codon_df[codon_df['Amino'] == amino_acid].sort_values(by='Frequency', ascending=False)
            alternative_codons = alternatives['Codon'].values
            
            # Find the next most favorable codon (skipping the current duplicate codon)
            for alt_codon in alternative_codons:
                if alt_codon != codons[to_mutate]:
                    codons[to_mutate] = alt_codon
                    print(f"Mutated codon at position {to_mutate} to {alt_codon}")
                    break
    
    # Join modified list back into a single DNA sequence
    return ''.join(codons)

# Ask if user wants to avoid codon duplets
avoid_duplets = input("Do you want to avoid codon duplets? (y/n): ").strip().lower() == 'y'
if avoid_duplets:
    optimized_dna_sequence = avoid_codon_duplets(optimized_dna_sequence, codon_df)
    print("Codon duplets minimized in the sequence.")
    print("New optimized sequence:")
    print(optimized_dna_sequence)
    gc_content = calculate_gc_content(optimized_dna_sequence)
    print(f"GC Content of optimized DNA sequence: {gc_content:.2f}%")

else:
    print("Codon duplets were not avoided.")
# slippery site avoidance
def avoid_slippery_sites(dna_sequence, codon_df):
    """
    Detects and avoids slippery sites (repeated bases) in the DNA sequence by mutating one of the affected codons.

    Parameters:
    dna_sequence (str): The DNA sequence to modify.
    codon_df (DataFrame): Dataframe containing codon frequency data.

    Returns:
    str: Modified DNA sequence with reduced slippery sites.
    """
    # Convert sequence to a list of codons
    codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]
    sequence_length = len(dna_sequence)
    
    # Slide a window of 4 bases through the sequence, one base at a time
    for i in range(sequence_length - 3):
        # Get the 4-base window
        window = dna_sequence[i:i+4]
        
        # Check if all four bases in the window are the same, indicating a slippery site
        if len(set(window)) == 1:  # All characters in the window are the same
            slippery_base = window[0]
            print(f"Slippery site detected: {window} at position {i}")
            
            # Identify the affected codons (two codons containing the slippery site)
            codon_start_1 = i // 3  # First codon in the frame
            codon_start_2 = (i + 3) // 3  # Next codon in the frame
            
            # Choose one of the affected codons to mutate, prioritizing higher frequency alternatives
            for codon_start in [codon_start_1, codon_start_2]:
                if codon_start >= len(codons):
                    continue  # Skip if out of range
                
                current_codon = codons[codon_start]
                amino_acid = codon_df.loc[codon_df['Codon'] == current_codon, 'Amino'].values[0]
                
                # Get alternatives for the amino acid, sorted by frequency
                alternatives = codon_df[(codon_df['Amino'] == amino_acid) & (codon_df['Codon'] != current_codon)]
                alternatives = alternatives.sort_values(by='Frequency', ascending=False)
                
                # Find the best alternative codon that does not introduce a new slippery site
                for alt_codon in alternatives['Codon']:
                    if slippery_base * 4 not in (alt_codon * 3):  # Ensure it doesn't reintroduce a slippery site
                        codons[codon_start] = alt_codon
                        print(f"Mutated codon at position {codon_start} to {alt_codon} to avoid slippery site")
                        break
            # Update dna_sequence to reflect the modified codons
            dna_sequence = ''.join(codons)
    
    return dna_sequence

# Ask if user wants to avoid slippery sites
avoid_slippery_sites_option = input("Do you want to avoid slippery sites? (y/n): ").strip().lower() == 'y'
if avoid_slippery_sites_option:
    optimized_dna_sequence = avoid_slippery_sites(optimized_dna_sequence, codon_df)
    print("Slippery sites minimized in the sequence.")
    print("New optimized sequence:")
    print(optimized_dna_sequence)
    gc_content = calculate_gc_content(optimized_dna_sequence)
    print(f"GC Content of optimized DNA sequence: {gc_content:.2f}%")
else:
    print("Slippery site avoidance not applied.")

# Generate a mapping of amino acids to their potential one-to-stop codons
amino_to_one_to_stop = {}
for aa in codon_df['Amino'].unique():
    aa_codons = codon_df[codon_df['Amino'] == aa]['Codon'].tolist()
    one_to_stop_codons = []
    for codon in aa_codons:
        for pos in range(3):
            for base in "ATGC":
                if base != codon[pos]:  # Ensure mutation
                    mutated_codon = codon[:pos] + base + codon[pos+1:]
                    if mutated_codon in stop_codons:
                        one_to_stop_codons.append(codon)
                        break
    amino_to_one_to_stop[aa] = list(set(one_to_stop_codons))  # Remove duplicates

#print("Amino Acid to One-to-Stop Codon Mapping:")
#for aa, codons in amino_to_one_to_stop.items():
#    print(f"{aa}: {codons}")

def one_to_stop_functionality(dna_sequence,mode="s", implementation_degree=50):
    """
    Modify codons based on one-to-stop functionality.
    
    :param mode: 'avoid' or 'implement'
    :param dna_sequence: DNA sequence as a string
    :param implementation_degree: Percentage of codons to mutate when implementing one-to-stop
    :return: Modified DNA sequence as a string
    """
    # Split DNA sequence into codons (triplets)
    codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]
    modified_codons = codons.copy()

    if mode == "a":
        for index, codon in enumerate(modified_codons):
            amino_acid = codon_df[codon_df['Codon'] == codon].iloc[0]['Amino']
            if codon in amino_to_one_to_stop.get(amino_acid, []):  # Codon is one-to-stop
                # Find alternative codons that are not one-to-stop
                aa_codons = codon_df[codon_df['Amino'] == amino_acid].sort_values(by='Frequency', ascending=False)['Codon'].tolist()
                alternative_codons = [alt for alt in aa_codons if alt not in amino_to_one_to_stop.get(amino_acid, [])]
                
                # Replace with the most frequent alternative, or keep the original if no alternative exists
                if alternative_codons:
                    print(f"Avoiding one-to-stop codon {codon} at position {index}. Replacing with {alternative_codons[0]}")
                    modified_codons[index] = alternative_codons[0]
                else:
                    print(f"No alternative codon available for {codon} at position {index}. Keeping original.")
    elif mode == "i":
        for index, codon in enumerate(modified_codons):
            amino_acid = codon_df[codon_df['Codon'] == codon].iloc[0]['Amino']
            
            # Skip if codon is already a one-to-stop
            if codon in amino_to_one_to_stop.get(amino_acid, []):
                continue

            # Check if a one-to-stop codon exists for this amino acid
            possible_one_to_stops = amino_to_one_to_stop.get(amino_acid, [])
            if possible_one_to_stops and random.randint(0, 100) < implementation_degree:
                # Use the most frequent one-to-stop codon
                chosen_one_to_stop = max(possible_one_to_stops, key=lambda c: codon_df[codon_df['Codon'] == c]['Frequency'].iloc[0])
                print(f"Implementing one-to-stop for {codon} at position {index}. Replacing with {chosen_one_to_stop}")
                modified_codons[index] = chosen_one_to_stop
            else:
                print(f"No one-to-stop codon available for {codon} at position {index}. Keeping original.")
    else:
        print("one-to-stop functionality skipped...")
    # Rejoin modified codons into a DNA sequence string
    return ''.join(modified_codons)

# Ask user for input
one_to_stop_mode = input("Choose one-to-stop handling mode (s - skip/ a - avoid/ i - implement): ").strip().lower()
if one_to_stop_mode == "i":
    degree = int(input("Enter degree of implementation (1-100%): ").strip())
else:
    degree = 0  # No degree needed for "skip" or "avoid"
# Apply one-to-stop handling
optimized_dna_sequence = one_to_stop_functionality(dna_sequence=optimized_dna_sequence, mode=one_to_stop_mode, implementation_degree=degree)
print(optimized_dna_sequence)
gc_content = calculate_gc_content(optimized_dna_sequence)
print(f"GC Content of optimized DNA sequence: {gc_content:.2f}%")
print("One-to-stop handling applied.")


# Cryptic splice avoidance
    # ask user if cryptic splice sites should be avoided
    # read in DBASS 3 & 5 databases .csv files
    # scan for splice sites
def preprocess_dbass(file_path):
    """
    Preprocess the DBASS database to extract 24-base motifs by sliding across each splice site motif.
    :param file_path: str, path to the DBASS file.
    :return: list of 24-base valid splice site motifs.
    """
    try:
        dbass = pd.read_csv(file_path)
    except FileNotFoundError:
        print(f"Database file {file_path} not found.")
        return []

    # Preprocess the NucleotideSequence column
    def clean_sequence(seq):
        # Convert to uppercase and remove invalid characters (only GATC are valid)
        return ''.join([char for char in seq.upper() if char in "GATC"])

    # Clean sequences and extract the cleaned sequence
    dbass['CleanedSequence'] = dbass['NucleotideSequence'].apply(lambda seq: clean_sequence(seq))

    # Create a list to store 24-base substrings
    motifs_24_base = []

    for seq in dbass['CleanedSequence'].dropna():
        # Slide a 24-base window across the splice site motif
        for i in range(len(seq) - 23):  # Ensure we only get 24-base substrings
            motifs_24_base.append(seq[i:i+24])

    # Remove duplicates and return the list of unique 24-base motifs
    return list(set(motifs_24_base))


def cryptic_splice_site_avoidance(dna_sequence, codon_df, avoid_splice_sites=True):
    """
    Avoid cryptic splice sites in a DNA sequence.
    :param dna_sequence: str, the optimized DNA sequence.
    :param codon_df: pandas DataFrame containing codon frequency data.
    :param avoid_splice_sites: bool, whether to avoid splice sites.
    :return: str, the modified DNA sequence.
    """
    if not avoid_splice_sites:
        return dna_sequence

    # Load and preprocess DBASS databases
    dbass3_motifs = preprocess_dbass(os.path.join(script_dir, "DBASS3.csv"))
    dbass5_motifs = preprocess_dbass(os.path.join(script_dir, "DBASS5.csv"))

    # Combine known splice site motifs
    splice_sites = dbass3_motifs + dbass5_motifs
    splice_sites = [site.upper() for site in splice_sites]  # Ensure uppercase for consistency

    print(f"Loaded {len(splice_sites)} cryptic splice site motifs.")

    # Sliding window approach
    codon_length = 3
    window_size = 24  # Adjust window size for 24 bases
    sequence_length = len(dna_sequence)

    for i in range(sequence_length - window_size + 1):  # Check 24-base windows
        window = dna_sequence[i:i + window_size]

        if window in splice_sites:
            print(f"Cryptic splice site identified: {window} at position {i}")

            # Determine the contributing codons
            codon_start_positions = [
                (i + offset) // codon_length * codon_length for offset in range(0, window_size, codon_length)
            ]

            affected_codons = [
                dna_sequence[pos:pos + codon_length] for pos in codon_start_positions if pos + codon_length <= sequence_length
            ]

            print(f"Affected codons: {affected_codons}")

            # Attempt to mutate one of the codons to avoid the splice site
            for codon in affected_codons:
                aa = codon_df[codon_df['Codon'] == codon]['Amino'].values[0]
                alt_codons = codon_df[(codon_df['Amino'] == aa) & (codon_df['Codon'] != codon)]

                if not alt_codons.empty:
                    new_codon = alt_codons.sort_values(by="Frequency", ascending=False)['Codon'].iloc[0]
                    codon_start = dna_sequence.find(codon)
                    dna_sequence = dna_sequence[:codon_start] + new_codon + dna_sequence[codon_start + codon_length:]
                    print(f"Mutated {codon} to {new_codon} to avoid cryptic splice site.")
                    break
            else:
                print(f"Could not mutate cryptic splice site {window} due to lack of alternative codons.")

    return dna_sequence

# Ask user whether to avoid cryptic splice sites
cryptic_avoidance = input("Should cryptic splice sites be avoided? (y/n): ").strip().upper()
if cryptic_avoidance == "Y":
    optimized_dna_sequence = cryptic_splice_site_avoidance(optimized_dna_sequence, codon_df)
    print(f"Optimized DNA sequence: {optimized_dna_sequence}")
    gc_content = calculate_gc_content(optimized_dna_sequence)
    print(f"GC Content of optimized DNA sequence: {gc_content:.2f}%")
else:
    print("Skipped cryptic splice site functionality...")

# avoid consensus splice.
def avoid_consensus_splice_sites(dna_sequence, codon_df):
    consensus_splice = "AGGTAAGT"
     # Sliding window approach
    codon_length = 3
    window_size = len(consensus_splice)  # Adjust window size to liking
    sequence_length = len(dna_sequence)
    for i in range(sequence_length - window_size + 1):  # Check 8-base windows
        window = dna_sequence[i:i + window_size]
        if window == consensus_splice:
            print(f"Splice site identified: {window} at position {i}")

            # Determine the contributing codons
            codon_start_positions = [
                (i + offset) // codon_length * codon_length for offset in range(0, window_size, codon_length)
            ]

            affected_codons = [
                dna_sequence[pos:pos + codon_length] for pos in codon_start_positions if pos + codon_length <= sequence_length
            ]

            print(f"Affected codons: {affected_codons}")

            # Attempt to mutate one of the codons to avoid the splice site
            for codon in affected_codons:
                aa = codon_df[codon_df['Codon'] == codon]['Amino'].values[0]
                alt_codons = codon_df[(codon_df['Amino'] == aa) & (codon_df['Codon'] != codon)]

                if not alt_codons.empty:
                    new_codon = alt_codons.sort_values(by="Frequency", ascending=False)['Codon'].iloc[0]
                    codon_start = dna_sequence.find(codon)
                    dna_sequence = dna_sequence[:codon_start] + new_codon + dna_sequence[codon_start + codon_length:]
                    print(f"Mutated {codon} to {new_codon} to avoid splice site.")
                    break
            else:
                print(f"Could not mutate splice site {window} due to lack of alternative codons.")
    return dna_sequence

# Ask user whether to avoid cryptic splice sites
consensus_splice_avoidance = input("Should consensus splice sites be avoided? (y/n): ").strip().upper()
if consensus_splice_avoidance == "Y":
    optimized_dna_sequence = avoid_consensus_splice_sites(optimized_dna_sequence, codon_df)
    print(f"Optimized DNA sequence: {optimized_dna_sequence}")
    gc_content = calculate_gc_content(optimized_dna_sequence)
    print(f"GC Content of optimized DNA sequence: {gc_content:.2f}%")
else:
    print("Skipped consensus splice site functionality...")
    
# choose motifs to avoid
motifs_to_avoid = []
enter_motifs = input("Do you want to enter motifs to avoid e.g. restriction sites? y/n\n").strip().upper()
if enter_motifs == "Y":
    enter_motifs = True
while enter_motifs is True:
    motif = input("Enter motif or 0 to exit\n").strip().upper()
    invalid = False
    if motif == "0":
       enter_motifs = False
    else:
        for letter in motif:
            if letter not in valid_bases:
                print("Invalid Motif")
                invalid = True
                continue      
        if not invalid:
            motifs_to_avoid.append(motif)
print(motifs_to_avoid) 
# scan for motifs to avoid ( e.g. restriction sites)
    # scan for the motif
    # define the three codons that make up the motif
    # vote which codon is favorable to mutate
