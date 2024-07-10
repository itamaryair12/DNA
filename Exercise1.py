#imprts

def extract_DNA_sequences(textFile):
    """
    Extract DNA sequences from a text file.

    Input:
    - file: text file with DNA sequences.

    Returns:
    - DNA sequences.
    """

    with open(textFile, 'r') as file:
        sequences = [line.strip() for line in file.readlines()]
    return sequences


def compute_GC_content(sequence):
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


def find_DNA_motifs(sequence, motif):
    """
    Find all occurrences of a given motif in a DNA sequence.

    Input:
    - sequence: DNA sequence.
    - motif: motif to be searched.

    Returns:
    - occurrences: list of occurrences of motif (location by index).
    """
    positions = []
    pos = 0  # in case we want to accept overlap
    while pos < len(sequence):
        pos = sequence.find(motif, pos)
        if pos == -1:
            break
        positions.append(pos)
        pos += 1  # Move to the next character to allow overlapping matches

    # in case we want to ignore overlap
    #pos = sequence.find(motif)
    #while pos != -1:
    #    positions.append(pos)
    #    pos = sequence.find(motif, pos + 1)

    return positions


def compute_codon_frequency(sequence):
    """
    Find all occurrences of a given motif in a DNA sequence.

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
