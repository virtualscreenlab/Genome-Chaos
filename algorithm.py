gene_sequence = list("ATGGCCCACGCCGGCCGCACCGGCTACGACAACCGCGAGATCGTGATGAAGTACATCCACTACAAGCTGAGCCAGCGCGGCTACGAGTGGGACGCCGGCGACGTGGGCGCCGCCCCCCCCGGCGCCGCCCCCGCCCCCGGCATCTTCAGCAGCCAGCCCGGCCACACCCCCCACCCCGCCGCCAGCCGCGACCCCGTGGCCCGCACCAGCCCCCTGCAGACCCCCGCCGCCCCCGGCGCCGCCGCCGGCCCCGCCCTGAGCCCCGTGCCCCCCGTGGTGCACCTGACCCTGCGCCAGGCCGGCGACGACTTCAGCCGCCGCTACCGCCGCGACTTCGCCGAGATGAGCAGCCAGCTGCACCTGACCCCCTTCACCGCCCGCGGCCGCTTCGCCACCGTGGTGGAGGAGCTGTTCCGCGACGGCGTGAACTGGGGCCGCATCGTGGCCTTCTTCGAGTTCGGCGGCGTGATGTGCGTGGAGAGCGTGAACCGCGAGATGAGCCCCCTGGTGGACAACATCGCCCTGTGGATGACCGAGTACCTGAACCGCCACCTGCACACCTGGATCCAGGACAACGGCGGCTGGGACGCCTTCGTGGAGCTGTACGGCCCCAGCATGCGCCCCCTGTTCGACTTCAGCTGGCTGAGCCTGAAGACCCTGCTGAGCCTGGCCCTGGTGGGCGCCTGCATCACCCTGGGCGCCTACCTGGGCCACAAG")

# Calculate the number of nucleotides
num_nucleotides = len(gene_sequence)

print("Number of Nucleotides:", num_nucleotides)


import random

def mutate_gene(gene_sequence, mutation_rate):
    mutated_sequence = list(gene_sequence)
    
    # Iterate through each nucleotide
    for i in range(len(mutated_sequence)):
        # Check if mutation occurs based on the mutation rate
        if random.random() < mutation_rate:
            # Mutate the nucleotide (replace it with a random nucleotide)
            mutated_sequence[i] = random.choice("ACGT")
            return mutated_sequence, i, mutated_sequence[i], True  # Return mutated sequence, position, mutated nucleotide, and True if mutation occurs
    
    return mutated_sequence, None, None, False  # Return mutated sequence, None, None, and False if no mutation occurs

# Original gene sequence
gene_sequence = list("ATGGCCCACGCCGGCCGCACCGGCTACGACAACCGCGAGATCGTGATGAAGTACATCCACTACAAGCTGAGCCAGCGCGGCTACGAGTGGGACGCCGGCGACGTGGGCGCCGCCCCCCCCGGCGCCGCCCCCGCCCCCGGCATCTTCAGCAGCCAGCCCGGCCACACCCCCCACCCCGCCGCCAGCCGCGACCCCGTGGCCCGCACCAGCCCCCTGCAGACCCCCGCCGCCCCCGGCGCCGCCGCCGGCCCCGCCCTGAGCCCCGTGCCCCCCGTGGTGCACCTGACCCTGCGCCAGGCCGGCGACGACTTCAGCCGCCGCTACCGCCGCGACTTCGCCGAGATGAGCAGCCAGCTGCACCTGACCCCCTTCACCGCCCGCGGCCGCTTCGCCACCGTGGTGGAGGAGCTGTTCCGCGACGGCGTGAACTGGGGCCGCATCGTGGCCTTCTTCGAGTTCGGCGGCGTGATGTGCGTGGAGAGCGTGAACCGCGAGATGAGCCCCCTGGTGGACAACATCGCCCTGTGGATGACCGAGTACCTGAACCGCCACCTGCACACCTGGATCCAGGACAACGGCGGCTGGGACGCCTTCGTGGAGCTGTACGGCCCCAGCATGCGCCCCCTGTTCGACTTCAGCTGGCTGAGCCTGAAGACCCTGCTGAGCCTGGCCCTGGTGGGCGCCTGCATCACCCTGGGCGCCTACCTGGGCCACAAG")

# Reference sequence (example, replace it with your reference sequence)
gene_sequence = reference_sequence 

# Mutation rate
mutation_rate = 1.0e-5

# Number of iterations
num_iterations = 1000

# Store the number of divisions for each iteration
division_counts = []

for _ in range(num_iterations):
    # Perform Monte Carlo simulation to estimate the number of divisions
    num_divisions = 0
    mutated_sequence = gene_sequence
    mutation_occurred = False

    while not mutation_occurred:
        num_divisions += 1
        mutated_sequence, _, _, mutation_occurred = mutate_gene(mutated_sequence, mutation_rate)

    division_counts.append(num_divisions)

# Calculate and print the average number of divisions
average_divisions = sum(division_counts) / num_iterations
print(f"Average number of divisions after {num_iterations} iterations: {average_divisions}")

import random

def mutate_gene(gene_sequence, mutation_rate):
    mutated_sequence = list(gene_sequence)
    
    # Iterate through each nucleotide
    for i in range(len(mutated_sequence)):
        # Check if mutation occurs based on the mutation rate
        if random.random() < mutation_rate:
            # Mutate the nucleotide (replace it with a random nucleotide)
            mutated_sequence[i] = random.choice("ACGT")
            return mutated_sequence, i, mutated_sequence[i], True  # Return mutated sequence, position, mutated nucleotide, and True if mutation occurs
    
    return mutated_sequence, None, None, False  # Return mutated sequence, None, None, and False if no mutation occurs

# Mutation rate
mutation_rate = 1.0e-5

# Perform Monte Carlo simulation to estimate the number of divisions
num_divisions = 0
mutated_sequence = gene_sequence
mutation_occurred = False

while not mutation_occurred:
    num_divisions += 1
    mutated_sequence, position, mutated_nucleotide, mutation_occurred = mutate_gene(mutated_sequence, mutation_rate)

# Compare with the reference sequence and print point mutations
if mutation_occurred:
    print("\nMutated Sequence:")
    print("".join(mutated_sequence))
    print("Reference Sequence:")
    print("".join(reference_sequence))
    print("Number of divisions needed:", num_divisions)
    print(f"Point Mutation occurred at position {position + 1}: {gene_sequence[position]} -> {mutated_nucleotide}")
else:
    print("No mutation occurred.")
    

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length for comparison.")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

# Reference sequence
reference_sequence = "ATGGCCCACGCCGGCCGCACCGGCTACGACAACCGCGAGATCGTGATGAAGTACATCCACTACAAGCTGAGCCAGCGCGGCTACGAGTGGGACGCCGGCGACGTGGGCGCCGCCCCCCCCGGCGCCGCCCCCGCCCCCGGCATCTTCAGCAGCCAGCCCGGCCACACCCCCCACCCCGCCGCCAGCCGCGACCCCGTGGCCCGCACCAGCCCCCTGCAGACCCCCGCCGCCCCCGGCGCCGCCGCCGGCCCCGCCCTGAGCCCCGTGCCCCCCGTGGTGCACCTGACCCTGCGCCAGGCCGGCGACGACTTCAGCCGCCGCTACCGCCGCGACTTCGCCGAGATGAGCAGCCAGCTGCACCTGACCCCCTTCACCGCCCGCGGCCGCTTCGCCACCGTGGTGGAGGAGCTGTTCCGCGACGGCGTGAACTGGGGCCGCATCGTGGCCTTCTTCGAGTTCGGCGGCGTGATGTGCGTGGAGAGCGTGAACCGCGAGATGAGCCCCCTGGTGGACAACATCGCCCTGTGGATGACCGAGTACCTGAACCGCCACCTGCACACCTGGATCCAGGACAACGGCGGCTGGGACGCCTTCGTGGAGCTGTACGGCCCCAGCATGCGCCCCCTGTTCGACTTCAGCTGGCTGAGCCTGAAGACCCTGCTGAGCCTGGCCCTGGTGGGCGCCTGCATCACCCTGGGCGCCTACCTGGGCCACAAG"

# Calculate Hamming distance between reference and mutated sequences
distance = hamming_distance(reference_sequence, mutated_gene)

print("Reference Sequence:", reference_sequence)
print("Mutated Gene Sequence:", mutated_gene)
print("Hamming Distance:", distance)

# Genetic code dictionary
genetic_code = {
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
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Function to translate a DNA sequence to a protein sequence
def translate_dna_to_protein(dna_sequence):
    protein_sequence = []
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i + 3]
        amino_acid = genetic_code.get(codon, 'X')  # 'X' for unknown codons
        protein_sequence.append(amino_acid)
    return ''.join(protein_sequence)

# Translate the reference and mutated sequences
translated_reference_sequence = translate_dna_to_protein(reference_sequence)
translated_mutated_sequence = translate_dna_to_protein(mutated_gene)

print("Translated Ref Gene Sequence:", translated_reference_sequence)
print("Translated Mut Gene Sequence:", translated_mutated_sequence)

# Function to calculate the percentage of sequence identity
def sequence_identity(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length for comparison.")
    identical_positions = sum(c1 == c2 for c1, c2 in zip(seq1, seq2))
    identity_percentage = (identical_positions / len(seq1)) * 100
    return identity_percentage

# Example protein sequences
protein_sequence1 = translated_reference_sequence
protein_sequence2 = translated_mutated_sequence

# Calculate the percentage of sequence identity
identity_percentage = sequence_identity(protein_sequence1, protein_sequence2)

print("Protein Sequence 1:", protein_sequence1)
print("Protein Sequence 2:", protein_sequence2)
print(f"Sequence Identity: {identity_percentage:.2f}%")

from Bio.SubsMat import MatrixInfo as matlist

# Example protein sequences
protein_sequence1 = "MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK"
protein_sequence2 = "MAHAWRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK"

# Load the BLOSUM62 similarity matrix
matrix = matlist.blosum62

# Initialize the similarity score
similarity_score = 0

# Ensure both sequences have the same length
if len(protein_sequence1) != len(protein_sequence2):
    raise ValueError("Sequences must have the same length for comparison.")

# Calculate the similarity score
for a, b in zip(protein_sequence1, protein_sequence2):
    if (a, b) in matrix:
        similarity_score += matrix[(a, b)]
    else:
        similarity_score += matrix[(b, a)]  # Check for reverse lookup

print("Protein Sequence 1:", protein_sequence1)
print("Protein Sequence 2:", protein_sequence2)
print("Similarity Score:", similarity_score)

# Provided protein sequences
translated_ref_sequence = "MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK"
translated_mut_gene_sequence = "MAHAWRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK"

# Initialize a counter for SNPs
snp_count = 0

# Ensure both sequences have the same length
if len(translated_ref_sequence) != len(translated_mut_gene_sequence):
    raise ValueError("Sequences must have the same length for comparison.")

# Compare sequences position by position
for ref_aa, mut_aa in zip(translated_ref_sequence, translated_mut_gene_sequence):
    if ref_aa != mut_aa:
        snp_count += 1

print("Number of SNPs:", snp_count)

# Provided protein sequences
translated_ref_sequence = "MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK"
translated_mut_gene_sequence = "MAHAWRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMRPLFDFSWLSLKTLLSLALVGACITLGAYLGHK"

# Ensure both sequences have the same length
if len(translated_ref_sequence) != len(translated_mut_gene_sequence):
    raise ValueError("Sequences must have the same length for comparison.")

# Iterate through sequences to find SNPs
for i, (ref_aa, mut_aa) in enumerate(zip(translated_ref_sequence, translated_mut_gene_sequence), start=1):
    if ref_aa != mut_aa:
        print(f"SNP at amino acid position {i}: {ref_aa} -> {mut_aa}")

# Function to calculate the percentage of sequence identity
def sequence_identity(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length for comparison.")
    identical_positions = sum(c1 == c2 for c1, c2 in zip(seq1, seq2))
    identity_percentage = (identical_positions / len(seq1)) * 100
    return identity_percentage

# Example protein sequences
protein_sequence1 = "MAHAGRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMR"
protein_sequence2 = "MAHAWRTGYDNREIVMKYIHYKLSQRGYEWDAGDVGAAPPGAAPAPGIFSSQPGHTPHPAASRDPVARTSPLQTPAAPGAAAGPALSPVPPVVHLTLRQAGDDFSRRYRRDFAEMSSQLHLTPFTARGRFATVVEELFRDGVNWGRIVAFFEFGGVMCVESVNREMSPLVDNIALWMTEYLNRHLHTWIQDNGGWDAFVELYGPSMR"

# Calculate the percentage of sequence identity
identity_percentage = sequence_identity(protein_sequence1, protein_sequence2)

print("Protein Sequence 1:", protein_sequence1)
print("Protein Sequence 2:", protein_sequence2)
print(f"Sequence Identity: {identity_percentage:.2f}%")

