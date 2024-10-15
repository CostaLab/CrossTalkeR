## New Features v1.4.0

- Splitted generate_report function in two parts:
  - analise_LR() to only run the analysis without generating the CrossTalkeR report
  - make_report() to only generate a new CrossTalkeR report for existing CrossTalkeR results
- Added node types to the network:
  - we now consider the annotation of a gene as ligand (L) or receptor (R) to consider the biological function 
- Less constrains on the cell cluster name annotation (only '$' must be avoided in the cluster naming)
- Integration with liana-py for ligand-receptor interaction predictions
  

## Old Features v1.3.0

- Single and Comparative Reports
   - Cell Cell Interaction visualization
   - Sending and Receiving Cells Ranking
   - Gene Target based Sankey Plots
   - CCI and GCI PCA ranking
      - All measures and PC table
      - PC1 and PC2 based barplot
   - Leimkühler et. al. [2] data were added to the package
   - Fisher Test were implemented to highlight the CCI edges significance (new) 🔥**NEW**🔥
   - **Change input format: Please see the Documentation** 🔥**NEW**🔥
      - A python3 notebook are available to cast the old input to the new input.
   - Liana (Dimitrov et. al. [3]) Output can be used as CrossTalkeR input. 🔥**NEW**🔥
   - LR pair visualization plot can be done using a Seurat Object 🔥**NEW**🔥
