# AFM-ProteinMPNN

We describe a computational approach to perform AlphaFold-Multimer (AFM) predictions using paired MSAs generated with ProteinMPNN interface design calculations. This enables improved prediction of synthetic protein-protein interactions based on computational models from which synthetic MSAs are derived. This computational approach, its rationale and benchmarks across multiple datasets of de novo designed miniprotein inhibitors, nanobodies and monobodies are presented in the article: "Improved prediction of synthetic protein-protein interactions with AlphaFold guided by ProteinMPNN" by Lourdes Carcelén, Alexandre Casadesús and Enrique Marcos.

The provided scripts are intented for running AFM predictions with ColabFold, either as individual predictions of binary complexes or in competition of two binders against the same receptor.

Dependencies:
- ColabFold
- ProteinMPNN
- PyRosetta
- TM-align
