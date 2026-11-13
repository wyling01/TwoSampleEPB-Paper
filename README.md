# TwoSampleEPB_Paper
This repo contains R scripts needed to reproduce the TwoSample EPB paper.

Dir 'EPB main' contains core implementations of VREPB and DVEPB (both for Microarrays Data and scRNA-seq Pseudo Bulk Data).

Command To run simulation with K_A = [k_A], K_B = [k_A] in the terminal: Rscript simulation.R [k_A] [k_B]

## Paper figures

From the repository root, choose an output directory and run:

```sh
Rscript figures/render_simulations.R /tmp/twosample-figures
Rscript figures/render_applications.R /tmp/twosample-figures
```

Both scripts use base R and read the existing, tracked results; they do not rerun simulations or fit models. The simulation script requires R with Cairo support. The application script also requires `pdflatex` on `PATH`, with the `standalone`, TikZ/PGF, `amsmath`, and Latin Modern packages (available in TeX Live/MacTeX).

The four vector PDFs are:

- `unequal_ka3.pdf`: unequal-variance simulation with $K_A=3$.
- `overview.pdf`: the other four main-text simulation settings.
- `macrophage_pair.pdf` and `chipseq_pair.pdf`: VREPB variance-ratio CDF beside the DVEPB joint variance prior.

Simulation inputs are the six FDR/Power CSVs under `Simulation/Simulation Result Summary/*/Full Data/`. Application inputs are the `1D_*.csv` and `2D_*.csv` fits under `Real Dataset/GSE9101/Result/` and `Real Dataset/chip-seq/Result/`.

Application axes use natural logs. CDFs accumulate all fitted mass, including tails outside the window. The joint-prior windows show 98.365954% of the macrophage mass and 99.834451% of the ChIP-seq mass, without renormalization; stems below 0.01% mass fade proportionally. Gray dashed lines indicate equal variances. The copper ChIP-seq reference is the MAnorm2 variance ratio 1.2756994798302 from `fit.info$ratio.var` in `chipseq.R` (a point-mass CDF on the left and the corresponding line on the right).
