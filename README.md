# Phylogeny Reconstruction with Simulated Bifurcation

This repository contains the code and released data used for the phylogenetic reconstruction experiments reported in the associated paper. Amino-acid (`AA`) and nucleotide (`NT`) analyses are organized separately, while shared neighbor-joining (NJ) and final comparison procedures are provided under `common/`.

Methodological definitions, experimental rationale, and interpretation of the results are described in the paper. This README documents the released files and how they correspond to the computational analyses.

## Associated paper

Citation information will be added upon publication.

## Repository structure

```text
.
├── AA/
│   ├── 1_Simulate_Dataset_AA.ipynb
│   ├── 2_Build_Similarity_Matrix_AA.ipynb
│   ├── 2A_Build_Evolutionary_Distance_Matrix_AA.ipynb
│   ├── 3_Reconstruct_Tree_AA.ipynb
│   ├── 4_Evaluate_Tree_AA_with_NJ.ipynb
│   ├── data/
│   │   ├── manifests/
│   │   │   └── simulation_manifest.csv
│   │   ├── trees/
│   │   │   └── *.nwk
│   │   ├── sim/
│   │   │   └── *.fa
│   │   └── nj_distances/
│   │       └── <variant>/*.csv
│   ├── models/
│   │   └── wag.dat
│   └── results_evaluation/
│       └── evaluation_aa_20260805_143349_222176.csv
├── NT/
│   ├── 1_Simulate_Dataset_NT.ipynb
│   ├── 2_Build_Similarity_Matrix_NT.ipynb
│   ├── 2A_Build_Evolutionary_Distance_Matrix_NT.ipynb
│   ├── 3_Reconstruct_Tree_NT.ipynb
│   ├── 4_Evaluate_Tree_NT_with_NJ.ipynb
│   ├── data/
│   │   ├── manifests/
│   │   │   └── simulation_manifest.csv
│   │   ├── trees/
│   │   │   └── *.nwk
│   │   ├── sim/
│   │   │   └── *.fa
│   │   └── nj_distances/
│   │       └── <variant>/*.csv
│   └── results_evaluation/
│       └── evaluation_nt_20260805_143621_077355.csv
└── common/
    ├── 5_Analyze_NJ_Comparison.ipynb
    └── Run_Neighbor_Joining.R
```

## Released data

### Simulation data

For both AA and NT:

- `data/trees/*.nwk`: simulated reference trees in Newick format.
- `data/sim/*.fa`: simulated sequence datasets in FASTA format.
- `data/manifests/simulation_manifest.csv`: simulation conditions, random seeds, and tree–sequence associations.

The AA analysis additionally uses:

- `AA/models/wag.dat`: WAG model parameters used in the AA evolutionary-distance calculation.

### NJ input matrices

For both AA and NT:

- `data/nj_distances/<variant>/*.csv`: distance matrices read directly by `common/Run_Neighbor_Joining.R`.

These files are needed when regenerating the NJ trees. They are not needed when reproducing the final aggregate analysis directly from the released evaluation CSV files.

### Reference evaluation results

The following files are the reference evaluation datasets used by the final analysis notebook:

```text
AA/results_evaluation/evaluation_aa_20260805_143349_222176.csv
NT/results_evaluation/evaluation_nt_20260805_143621_077355.csv
```

They contain the per-dataset results from which the reported aggregate comparisons are calculated.

## Reproducibility routes

### Final aggregate analysis

Run:

```text
common/5_Analyze_NJ_Comparison.ipynb
```

using the two released evaluation CSV files. This reproduces the aggregate tables and figures without rerunning matrix construction, SB reconstruction, or NJ.

### Tree reconstruction and evaluation

Using the released reference trees, FASTA files, manifests, and NJ distance matrices, run the following for AA:

```text
AA/2_Build_Similarity_Matrix_AA.ipynb
AA/2A_Build_Evolutionary_Distance_Matrix_AA.ipynb
AA/3_Reconstruct_Tree_AA.ipynb
common/Run_Neighbor_Joining.R
AA/4_Evaluate_Tree_AA_with_NJ.ipynb
```

Run the corresponding NT notebooks for the nucleotide analysis.

This route regenerates the derived matrices, reconstructed trees, NJ trees, and evaluation tables.

### Simulation regeneration

Run:

```text
AA/1_Simulate_Dataset_AA.ipynb
NT/1_Simulate_Dataset_NT.ipynb
```

The released manifests retain the simulation seeds and file associations used for the released datasets. Exact byte-for-byte regeneration may depend on the versions and behavior of R, IQ-TREE/AliSim, and their dependencies.

## File map

| File | Role |
|---|---|
| `1_Simulate_Dataset_*.ipynb` | Generate reference trees, simulated sequences, and manifests |
| `2_Build_Similarity_Matrix_*.ipynb` | Build BLAST-derived similarity matrices |
| `2A_Build_Evolutionary_Distance_Matrix_*.ipynb` | Build model-based evolutionary-distance matrices |
| `3_Reconstruct_Tree_*.ipynb` | Reconstruct trees using recursive normalized-cut SB optimization |
| `Run_Neighbor_Joining.R` | Generate NJ trees from released distance-matrix CSV files |
| `4_Evaluate_Tree_*_with_NJ.ipynb` | Compare reconstructed and NJ trees with the reference trees |
| `5_Analyze_NJ_Comparison.ipynb` | Produce the final AA/NT aggregate comparisons |

The notebooks are intended to be run with the repository structure shown above. Locations of external executables can be configured in the relevant notebook cells.

## Public release manifests

The released manifests are portability-oriented versions of the original execution manifests.

The following changes were made for release:

- environment-specific absolute paths were converted to paths relative to the corresponding `AA/` or `NT/` directory;
- execution-only fields, including log-file paths and run messages, were removed;
- the AA and NT manifests were normalized to a common schema;
- `seqtype` was recorded explicitly as `AA` or `DNA`;
- simulation parameters, tags, random seeds, status values, and tree–sequence associations were retained unchanged.

The common schema is:

```text
tag
generator
n_taxa
target_mean_branch_length
seqtype
model
seq_length
tree_seed
seq_seed
tree_path
fasta_path
status
```

For example:

```text
AA/data/sim/rtree_n30_bl0.125_rep001.fa
```

is stored in the AA manifest as:

```text
data/sim/rtree_n30_bl0.125_rep001.fa
```

and is resolved relative to `AA/`.

## Expected data counts

The released simulation design contains:

- 2 tree generators;
- 5 target mean branch-length conditions;
- 100 replicates per condition;
- 1,000 AA datasets;
- 1,000 NT datasets.

Each sequence type therefore contains:

```text
1,000 manifest rows
1,000 reference-tree files
1,000 FASTA files
```

The released evaluation datasets contain:

```text
AA: 13,000 rows
NT: 10,000 rows
```

These values can be used as basic completeness checks.

## Software

The workflow uses:

### Python

- Jupyter
- NumPy
- pandas
- SciPy
- Biopython
- matplotlib
- PyTorch
- `simulated-bifurcation`

### R

- R
- `ape`
- `TreeSim`

### External programs

- IQ-TREE 3 / AliSim (`iqtree3`)
- NCBI BLAST+ (`makeblastdb`, `blastp`, and `blastn`)

The software versions used for the reported analyses are documented in the associated paper.

## Randomness and numerical reproducibility

Simulation seeds are stored in the released manifests.

The SB reconstruction settings are recorded in the reconstruction notebooks. Rerunning SB may not produce byte-identical intermediate outputs or identical individual reconstructed trees across different software versions, hardware, numerical precision, or execution backends.

The released evaluation CSV files are therefore provided as the reference outputs used for the reported final analysis.

## Generated files not included

The following directories contain intermediate or derived outputs and are regenerated by the workflow:

```text
AA/data/blastdb/
AA/data/matrices/
AA/data/evo_distances/
AA/results_reconstruct/
AA/results_nj/

NT/data/blastdb/
NT/data/matrices/
NT/data/evo_distances/
NT/results_reconstruct/
NT/results_nj/

paper_analysis_nj_rep100/
```

They are not required in the initial repository state.

## Relationship to the paper

The associated paper provides:

- definitions of the reconstruction objectives and distance variants;
- the rationale for the experimental design;
- methodological details of normalized-cut and SB reconstruction;
- evaluation metrics and statistical interpretation;
- discussion of the biological and computational implications.

This repository provides the corresponding executable code, released inputs, random seeds, and reference evaluation outputs.

## Citation

Citation information will be added upon publication.

The WAG model parameters are not original to this repository. The original WAG model and the source of the distributed parameter file should be cited as described in the associated paper.

## License

See `LICENSE` for the terms covering the code and released data.
