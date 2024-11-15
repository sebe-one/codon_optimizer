def handle_one_to_stop(dna_sequence, codon_df, mode="skip", implementation_degree=0):
    """
    Handles one-to-stop codons by either avoiding or implementing them.

    Parameters:
    dna_sequence (str): The DNA sequence to modify.
    codon_df (DataFrame): Dataframe containing codon frequency data.
    mode (str): "skip", "avoid", or "implement".
    implementation_degree (int): Percentage chance to implement one-to-stop codons.

    Returns:
    str: Modified DNA sequence.
    """
    # Define stop codons
    stop_codons = {"TAA", "TAG", "TGA"}

    # Convert sequence to a list of codons
    codons = [dna_sequence[i:i+3] for i in range(0, len(dna_sequence), 3)]

    if mode == "a":
        # Avoid one-to-stop codons
        for index, codon in enumerate(codons):
            # Generate all possible single-base mutations for the codon
            mutations = []
            for pos in range(3):
                for base in "ATGC":
                    if base != codon[pos]:  # Ensure it's a mutation
                        mutated_codon = codon[:pos] + base + codon[pos+1:]
                        mutations.append(mutated_codon)

            # Check if any mutation results in a stop codon
            if any(mutated in stop_codons for mutated in mutations):
                print(f"One-to-stop codon detected: {codon} at position {index}")

                # Find an alternative codon for the amino acid
                amino_acid = codon_df.loc[codon_df['Codon'] == codon, 'Amino'].values[0]
                alternatives = codon_df[(codon_df['Amino'] == amino_acid) & (codon_df['Codon'] != codon)]
                alternatives = alternatives.sort_values(by='Frequency', ascending=False)

                # Replace with the first valid alternative that is not a one-to-stop
                for alt_codon in alternatives['Codon']:
                    alt_mutations = [alt_codon[:pos] + base + alt_codon[pos+1:] for pos in range(3) for base in "ATGC" if base != alt_codon[pos]]
                    if not any(mutated in stop_codons for mutated in alt_mutations):
                        codons[index] = alt_codon
                        print(f"Mutated codon at position {index} to {alt_codon} to avoid one-to-stop potential")
                        break

    elif mode == "i":
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
        
        # Print the mapping for debugging
        print("One-to-Stop Codon Mapping:")
        for aa, codons in amino_to_one_to_stop.items():
            print(f"{aa}: {codons}")

        # Implement one-to-stop codons in the DNA sequence
        for index, codon in enumerate(codons):
            # Skip if the current codon is already a one-to-stop
            if any(mutated_codon in stop_codons for mutated_codon in [
                    codon[:pos] + base + codon[pos+1:]
                    for pos in range(3) for base in "ATGC" if base != codon[pos]]):
                continue

            # Determine the amino acid encoded by the codon
            aa_row = codon_df[codon_df['Codon'] == codon]
            if not aa_row.empty:
                amino_acid = aa_row.iloc[0]['Amino']

                # Get the one-to-stop options for this amino acid
                possible_one_to_stops = amino_to_one_to_stop.get(amino_acid, [])

                # If there are valid one-to-stop options, mutate based on the threshold
                if possible_one_to_stops and random.randint(0, 100) < implementation_degree:
                    chosen_one_to_stop = random.choice(possible_one_to_stops)
                    print(f"Mutating codon {codon} at position {index} to one-to-stop: {chosen_one_to_stop}")
                    codons[index] = chosen_one_to_stop
    else:
        print("One to stop functionality skipped.")
    # Reconstruct the DNA sequence
    return ''.join(codons)
    
# Ask user for input
one_to_stop_mode = input("Choose one-to-stop handling mode (s - skip/ a - avoid/ i - implement): ").strip().lower()
if one_to_stop_mode == "i":
    degree = int(input("Enter degree of implementation (1-100%): ").strip())
else:
    degree = 0  # No degree needed for "skip" or "avoid"

# Apply one-to-stop handling
optimized_dna_sequence = handle_one_to_stop(optimized_dna_sequence, codon_df, mode=one_to_stop_mode, implementation_degree=degree)
