#imprts
import json
def extract_dna_sequences(textFile):
    """
    Extract DNA sequences from a text file.

    Input:
    - file: text file with DNA sequences.

    Returns:
    - DNA sequences.
    """

    try:
        with open(textFile, 'r') as file:
            sequences = [line.strip() for line in file.readlines()]
        return sequences
    except IOError as e:
        print(f"Error reading file {textFile}: {e}")
        return []


def compute_gc_content(sequence):
    """
    Compute the GC content of a DNA sequence.

    Input:
    - sequence: DNA sequence.

    Returns:
    - GC content: GC content of the DNA sequence.
    """
    amount_of_gc = 0
    for i in range(len(sequence)):
        if sequence[i] == 'G' or sequence[i] == 'C':
            amount_of_gc += 1
    return (amount_of_gc / float(len(sequence))) * 100


def find_dna_motifs(sequence, motif):
    """
    Find all occurrences of a given motif in a DNA sequence.
    NOTE: This function includes overlapping positions of the motif if they occur.

    Input:
    - sequence: DNA sequence.
    - motif: motif to be searched.

    Returns:
    - occurrences: list of occurrences of motif (location by index).
    """
    positions = []
    pos = sequence.find(motif)
    while pos != -1:
       positions.append(pos + 1)  # Convert to start from 1 instead of 0
       pos = sequence.find(motif, pos + 1)

    return positions


def compute_codon_frequency(sequence):
    """
    Calculate the frequency of each codon.

    Input:
    - sequence: DNA sequence.

    Returns:
    - dictionary: A dictionary with codons as keys and their frequencies as values.
    """
    codon_dictionary = {}
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i + 3]
        if codon in codon_dictionary:
            codon_dictionary[codon] += 1
        else:
            codon_dictionary[codon] = 1
    return codon_dictionary


def longest_common_subsequence(seq1, seq2):
    """
    Find the longest common subsequence between two DNA sequences.

    Input:
    - seq1: DNA sequences list.
    - seq2: DNA sequences list.

    Returns:
    - subsequence: The longest common subsequence.
    """
    # Initialize the matrix with zeros
    lengths = []
    for i in range(len(seq1) + 1):
        lengths.append([0] * (len(seq2) + 1))
    longest_length = 0
    longest_end_pos = 0

    # Populate the lengths matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            if seq1[i - 1] == seq2[j - 1]:
                lengths[i][j] = lengths[i - 1][j - 1] + 1
                if lengths[i][j] > longest_length:
                    longest_length = lengths[i][j]
                    longest_end_pos = i

    # Extract the longest common substring
    longest_common_substring = seq1[longest_end_pos - longest_length:longest_end_pos]
    return longest_common_substring


def main(textfile, motif):
    sequences = extract_dna_sequences(textfile)

    results = {
        "motif": motif,
        "Sequences": [],
        "most_common_codon": None,
        "longest_common_substring": {
            "value": "",
            "sequences": [],
            "length": 0
        }
    }


    #initialize parameters
    all_codons = {}
    longest_common_val = ""
    longest_common_length = 0

    for i, sequence in enumerate(sequences):
        gc_content = compute_gc_content(sequence)
        motif_positions = find_dna_motifs(sequence, motif)
        codon_freq = compute_codon_frequency(sequence)

        for codon, freq in codon_freq.items():
            if codon in all_codons:
                all_codons[codon] += freq
            else:
                all_codons[codon] = freq

        results["Sequences"].append({
            "gc_content": round(gc_content, 2),
            "motif_positions": motif_positions if motif_positions else None,
            "codons": codon_freq
        })


        for j in range(i+1, len(sequences)):
            lcs_value = longest_common_subsequence(sequence, sequences[j])
            if len(lcs_value) > longest_common_length:
                longest_common_val = lcs_value
                longest_common_length = len(lcs_value)
                results["longest_common_substring"] = {
                    "value": lcs_value,
                    "sequences": [i + 1, j + 1],
                    "length": len(lcs_value)
                }

    results["most_common_codon"] = max(all_codons, key=all_codons.get)

    print(json.dumps(results, indent=4))

# The following lines are for testing purposes
#if __name__ == "__main__":
#    main("tester.txt", "GCTAGCTAGC")
