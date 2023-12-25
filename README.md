# Genome-Chaos
 Gene Mutation and Translation

The wild-type BCL2 gene sequence (ENST00000398117.1) was retrieved from the COSMIC database. A Python script was utilized to simulate genetic mutations with the following parameters: a mutation rate of 10-5 and about 86 cell divisions. 
The Python script simulates genetic mutations in a gene sequence using a Monte Carlo approach. The mutate_gene function introduces mutations based on a given muta-tion rate, altering nucleotides randomly. The main part of the script conducts a Monte Carlo simulation, repeatedly applying the mutation function to estimate the average number of iterations (divisions) required for a mutation to occur. After 1000 iterations, the average number of divisions was found to be 86.36.
A genetic code dictionary was employed to translate both the mutated gene sequence and a reference gene sequence into protein sequences. A translation function, trans-late_dna_to_protein, utilized this dictionary to convert both the reference and mutated gene sequences into protein sequences.

Protein Stability Prediction

The wild-type (wt) structure, a component of the Bcl2-BINDI complex, was sourced from the Protein Data Bank as a crystal structure (PDB ID: 5JSN) with a resolution of 2.1 Å. Predictions for mutant forms of the Bcl-2 protein were generated using the SWISS-MODEL algorithm according to the protocol published elsewhere. To assess protein stability, the standard Rosetta protocol was implemented, involving the calculation of the energy score (ES). This stability prediction method was executed by aligning root-mean-square-deviation (RMSD) values with Rosetta energy parameters. The differences in Rosetta energy scores for mutated forms (ΔES) were computed using the following equation:

𝛥𝐸𝑆=𝐸𝑆𝑚𝑢𝑡−𝐸𝑆𝑤𝑡 

where ESmut and ESwt are the Rosetta energy score values for the wild type and mutated forms, respectively.

Protein Function Prediction

The crystal structure of Bcl-2 in association with a Bax BH3 peptide was retrieved from the Protein Data Bank (PDB ID: 2XA0) with a resolution of 2.7 Å to determine the Bcl2-Bax binding site. The molecular interaction between the protein and the peptide was evaluated utilizing the ZDOCK molecular docking server for the prediction of protein-protein complexes and symmetric multimers. The binding energy (Ebind) was calculated by using the PPI-affinity tool designed to predict and optimize the binding affinity of protein-peptide and protein-protein complexes.
