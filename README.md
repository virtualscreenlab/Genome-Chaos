# Genome-Chaos: Gene Mutation and Protein Stability Analysis

## Overview  
This repository contains tools for simulating gene mutations, predicting protein stability, and analyzing protein-protein interactions, with a focus on cancer-related genome instability. The project integrates bioinformatics and computational biology approaches to study the relationship between genetic mutations and protein function.

## Features  

- **Gene Mutation Simulation**: Monte Carlo simulation of genetic mutations with configurable mutation rates  
- **Protein Translation**: Conversion of nucleotide sequences to protein sequences using standard genetic code  
- **Protein Stability Prediction**: Rosetta protocol implementation for energy score calculations  
- **Protein-Protein Interaction Analysis**: Binding energy calculations for protein complexes  
- **Genome Instability Analysis**: Tools for studying chaos in cancer cell populations  

## Installation  

### Prerequisites  
- Python 3.7+  
- Rosetta (for protein stability predictions)  
- BioPython  
- NumPy  

### Installation Steps  
```bash  
git clone https://github.com/yourusername/Genome-Chaos.git  
cd Genome-Chaos  
pip install -r requirements.txt  
```  

## Usage  

### Gene Mutation Simulation  
```python  
from algorithm import simulate_mutations  

# Simulate mutations with default parameters  
results = simulate_mutations(gene_sequence, mutation_rate=1e-5, iterations=1000)  
print(f"Average divisions until mutation: {results['average_divisions']}")  
```  

### Protein Translation  
```python  
from algorithm import translate_sequence  

protein_sequence = translate_sequence(gene_sequence)  
```  

### Protein Stability Prediction  
```bash  
python analysis.py --pdb wildtype.pdb --mutations mut.dat  
```  

## File Structure  

```
Genome-Chaos/  
├── algorithm.py          # Core mutation and translation algorithms  
├── algorithm.ipynb       # Jupyter notebook with examples  
├── analysis.py           # Protein stability analysis tools  
├── models/               # Protein structure models  
│   ├── model_01_bxE2_nn.tpdb  
│   └── model_01_bxE2_wt.pdb  
├── data/                 # Genetic and mutation data  
│   ├── mut.dat  
│   └── mt.dat  
├── docs/                 # Documentation and research  
│   └── RCL2_gene.html  
└── README.md             # This file  
```  

## Research Background  

This work investigates the relationship between genome instability in cancer cells and protein function. The methodology includes:  

1. Retrieving wild-type gene sequences from COSMIC database  
2. Simulating mutations with Monte Carlo methods  
3. Predicting protein stability changes using Rosetta protocols  
4. Analyzing protein-protein interactions with molecular clustering  

Key equations:  
- ΔES = ESmut - ESwt (Energy score difference)  
- Binding energy calculations using PH-affinity tools  

## Contributing  

Contributions are welcome. Please fork the repository and submit pull requests. For major changes, please open an issue first to discuss proposed changes.  

## License  

This project is licensed under the MIT License - see the LICENSE file for details.  

## References  

1. COSMIC database (https://cancer.sanger.ac.uk/cosmic)  
2. Protein Data Bank (https://www.rcsb.org/)  
3. Rosetta Commons (https://www.rosettacommons.org/)  
4. SWISS-MODEL algorithm for protein structure prediction
