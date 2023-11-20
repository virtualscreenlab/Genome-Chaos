gene_sequence = list("ATGGCGCACGCTGGGAGAACAGGGTACGATAACCGGGAGATAGTGATGAAGTACATCCATTATAAGCTGTCGCAGAGGGGCTACGAGTGGGATGCGGGAGATGTGGGCGCCGCGCCCCCGGGGGCCGCCCCCGCACCGGGCATCTTCTCCTCCCAGCCCGGGCACACGCCCCATCCAGCCGCATCCCGGGACCCGGTCGCCAGGACCTCGCCGCTGCAGACCCCGGCTGCCCCCGGCGCCGCCGCGGGGCCTGCGCTCAGCCCGGTGCCACCTGTGGTCCACCTGACCCTCCGCCAGGCCGGCGACGACTTCTCCCGCCGCTACCGCCGCGACTTCGCCGAGATGTCCAGCCAGCTGCACCTGACGCCCTTCACCGCGCGGGGACGCTTTGCCACGGTGGTGGAGGAGCTCTTCAGGGACGGGGTGAACTGGGGGAGGATTGTGGCCTTCTTTGAGTTCGGTGGGGTCATGTGTGTGGAGAGCGTCAACCGGGAGATGTCGCCCCTGGTGGACAACATCGCCCTGTGGATGACTGAGTACCTGAACCGGCACCTGCACACCTGGATCCAGGATAACGGAGGCTGGGATGCCTTTGTGGAACTGTACGGCCCCAGCATGCGGCCTCTGTTTGATTTCTCCTGGCTGTCTCTGAAGACTCTGCTCAGTTTGGCCCTGGTGGGAGCTTGCATCACCCTGGGTGCCTATCTGGGCCACAAGTGA")

# Calculate the number of nucleotides
num_nucleotides = len(gene_sequence)

print("Number of Nucleotides:", num_nucleotides)

import random

# Mutation rate
mutation_rate = 1.0e-4

# Number of cell divisions
num_divisions = 60

# Initialize the gene sequence (as a list of characters, for example)
gene_sequence = list("ATGGCGCACGCTGGGAGAACAGGGTACGATAACCGGGAGATAGTGATGAAGTACATCCATTATAAGCTGTCGCAGAGGGGCTACGAGTGGGATGCGGGAGATGTGGGCGCCGCGCCCCCGGGGGCCGCCCCCGCACCGGGCATCTTCTCCTCCCAGCCCGGGCACACGCCCCATCCAGCCGCATCCCGGGACCCGGTCGCCAGGACCTCGCCGCTGCAGACCCCGGCTGCCCCCGGCGCCGCCGCGGGGCCTGCGCTCAGCCCGGTGCCACCTGTGGTCCACCTGACCCTCCGCCAGGCCGGCGACGACTTCTCCCGCCGCTACCGCCGCGACTTCGCCGAGATGTCCAGCCAGCTGCACCTGACGCCCTTCACCGCGCGGGGACGCTTTGCCACGGTGGTGGAGGAGCTCTTCAGGGACGGGGTGAACTGGGGGAGGATTGTGGCCTTCTTTGAGTTCGGTGGGGTCATGTGTGTGGAGAGCGTCAACCGGGAGATGTCGCCCCTGGTGGACAACATCGCCCTGTGGATGACTGAGTACCTGAACCGGCACCTGCACACCTGGATCCAGGATAACGGAGGCTGGGATGCCTTTGTGGAACTGTACGGCCCCAGCATGCGGCCTCTGTTTGATTTCTCCTGGCTGTCTCTGAAGACTCTGCTCAGTTTGGCCCTGGTGGGAGCTTGCATCACCCTGGGTGCCTATCTGGGCCACAAGTGA")

# Function to apply mutations
def apply_mutation(gene_sequence, mutation_rate):
    mutated_sequence = []
    for base in gene_sequence:
        if random.random() < mutation_rate:
            # Mutate the base (for simplicity, just change to a random base)
            mutated_base = random.choice("ATCG")
            mutated_sequence.append(mutated_base)
        else:
            mutated_sequence.append(base)
    return mutated_sequence

# Simulate cell divisions and mutations
for _ in range(num_divisions):
    gene_sequence = apply_mutation(gene_sequence, mutation_rate)

# Print the mutated gene sequence
mutated_gene = ''.join(gene_sequence)
print("Mutated Gene Sequence:", mutated_gene)

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length for comparison.")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

# Reference sequence
reference_sequence = "ATGGCGCACGCTGGGAGAACAGGGTACGATAACCGGGAGATAGTGATGAAGTACATCCATTATAAGCTGTCGCAGAGGGGCTACGAGTGGGATGCGGGAGATGTGGGCGCCGCGCCCCCGGGGGCCGCCCCCGCACCGGGCATCTTCTCCTCCCAGCCCGGGCACACGCCCCATCCAGCCGCATCCCGGGACCCGGTCGCCAGGACCTCGCCGCTGCAGACCCCGGCTGCCCCCGGCGCCGCCGCGGGGCCTGCGCTCAGCCCGGTGCCACCTGTGGTCCACCTGACCCTCCGCCAGGCCGGCGACGACTTCTCCCGCCGCTACCGCCGCGACTTCGCCGAGATGTCCAGCCAGCTGCACCTGACGCCCTTCACCGCGCGGGGACGCTTTGCCACGGTGGTGGAGGAGCTCTTCAGGGACGGGGTGAACTGGGGGAGGATTGTGGCCTTCTTTGAGTTCGGTGGGGTCATGTGTGTGGAGAGCGTCAACCGGGAGATGTCGCCCCTGGTGGACAACATCGCCCTGTGGATGACTGAGTACCTGAACCGGCACCTGCACACCTGGATCCAGGATAACGGAGGCTGGGATGCCTTTGTGGAACTGTACGGCCCCAGCATGCGGCCTCTGTTTGATTTCTCCTGGCTGTCTCTGAAGACTCTGCTCAGTTTGGCCCTGGTGGGAGCTTGCATCACCCTGGGTGCCTATCTGGGCCACAAGTGA"

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

