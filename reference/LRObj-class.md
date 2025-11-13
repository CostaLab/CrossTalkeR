# Run and Generate all LR Downstream analysis

This function loads the single conditions LR outputs and return the LR
network It assumes that the table present the following columns Ligand,
measure

## Slots

- `graphs`:

  All Cell Cell Interaction Networks

- `tables`:

  All tables from single condition

- `max_iter`:

  Max meanLR from all

- `max_nodes`:

  All Celltype in the experiment

- `coords`:

  Cell Cell Interaction Plots

- `colors`:

  Cell type colors

- `rankings`:

  Ranking of cells and Genes

- `loadings`:

  CCI values to remove multiple times genes

- `pca`:

  PCA results

- `annot`:

  Annotation Results

- `stats`:

  Statistics of the analysis
