# Genome-Chaos
 Gene Mutation and Translation

The BCL2 cDNA sequence (ENST00000398117.1) was retrieved from the COSMIC database. A Python script was utilized to simulate genetic mutations with the following parameters: a mutation rate of 10-4 and 60 cell divisions (Hayflick limit). The gene sequence was initialized as a list of characters, each representing a specific DNA sequence for the simulation. The apply_mutation function was used to randomly alter bases in the gene sequence according to the mutation rate. For each base, there was a probability of mutation, where the base could change to any nucleotide (A, T, C, or G). A genetic code dictionary was employed to translate both the mutated gene sequence and a reference gene sequence into protein sequences. A translation function, translate_dna_to_protein, uti-lized this dictionary to convert both the reference and mutated gene sequences into protein sequences.

Protein Stability Prediction

The wild-type (wt) structure, a component of the Bcl2-BINDI complex, was sourced from the Protein Data Bank as a crystal structure (PDB ID: 5JSN) with a resolution of 2.1 Å. Predictions for mutant forms of the Bcl-2 protein were generated using the SWISS-MODEL algorithm according to the protocol published elsewhere. To assess protein stability, the standard Rosetta protocol was implemented, involving the calculation of the energy score (ES). This stability prediction method was executed by aligning root-mean-square-deviation (RMSD) values with Rosetta energy parameters. The differences in Rosetta energy scores for mutated forms (ΔES) were computed using the following equation:

𝛥𝐸𝑆=𝐸𝑆𝑚𝑢𝑡−𝐸𝑆𝑤𝑡 

where ESmut and ESwt are the Rosetta energy score values for the wild type and mutated forms, respectively.

Protein Function Prediction

The crystal structure of Bcl-2 in association with a Bax BH3 peptide was retrieved from the Protein Data Bank (PDB ID: 2XA0) with a resolution of 2.7 Å to determine the Bcl2-Bax binding site. The molecular interaction between the protein and the peptide was evaluated utilizing the ZDOCK molecular docking server for the prediction of protein-protein complexes and symmetric multimers. The binding energy (Ebind) was calculated by using the PPI-affinity tool designed to predict and optimize the binding affinity of protein-peptide and protein-protein complexes.
