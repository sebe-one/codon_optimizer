
import pandas as pd
import openpyxl
import numpy as np
import os
from collections import Counter
import math
def reverse_complement(sequence):
    # Define complement mapping
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    
    # Reverse the sequence and map each base to its complement
    reversed_complemented_sequence = "".join(complement[base] for base in reversed(sequence))
    
    return reversed_complemented_sequence

# Define valid characters
valid_bases = ["G", "A", "T", "C", "*"]
valid_AA = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
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
    print("Optimization:\n  - Most Frequent\n 2 - Probability Frequency Distribution\n 3 - Enforced Frequency Distribution\n")
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
# inverted frequency for deoptimization
if algorithm > 3:
    # Invert frequencies in codon_df
    codon_df['Inverted_Frequency'] = 1 - codon_df['Frequency']
    codon_df['Inverted_Frequency'] /= codon_df.groupby('Amino')['Inverted_Frequency'].transform('sum')  # Normalize within each amino acid        
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
# avoid codon duplication 
    # ask if it is wanted
    # sliding window approach of 12 bases
    # if the same codon is used twice in the window, flip a coin and mutate either the first or the second one to the next favorable one in the codon frequence table

# one to stop avoidance or implementation
    # ask if it it should be skipped, avoided or implemented
    # if avoid: 
        # scan for one to stop codons in frame or out of frame
        # define the (two) codons that make up the one to stop
        # if out of frame:
            # vote which codon is more favorable to mutate by checking the next best frequency of those two
        # if in frame: mutate to the next favorable option
    # if implement: 
        # ask for the degree of implementation (int between 1 and 100 %)
        # ask if only in frame or total
        # identify possible one to stop options
        # for each option choose a random number between 0 and 101, if number > than cutoff mutate the codon to a one to stop.

# Cryptic splice avoidance
    # read in DBASS 3 & 5 databases
    # scan for splice sites
# scan for motifs to avoid ( e.g. restriction sites)
    # scan for the motif
    # define the three codons that make up the motif
    # vote which codon is favorable to mutate
