# TwoSampleEPB

Reproduction code for “*Empirical partially Bayes two sample testing*”. The repository contains the core
VREPB and DVEPB implementations, three data applications, the simulation study,
and scripts that render the paper figures.

## Code provenance

The core methods are implemented in `EPB main/EPB.R`. Weighted versions used
with voom and VoombyGroup are implemented in `EPB main/voom_EPB.R`.
`helper functions/voomByGroup.R` contains the VoombyGroup implementation used by
the PBMC1 application.

External application data are not stored in the repository. Each application
script downloads its public input before preprocessing. The ChIP-seq and PBMC1
downloads are pinned to specific Git commits; the GSE9101 application downloads
the public GEO supplementary archive.

## Setup

The tested R, Bioconductor, package, and external-software versions are recorded
in [`reproducibility.yml`](reproducibility.yml).

NPMLE fitting uses MOSEK through `Rmosek` and requires a working MOSEK
installation and license. If the license is not in MOSEK's default location,
set:

```bash
export MOSEKLM_LICENSE_FILE=/absolute/path/to/mosek.lic
```

The application and simulation scripts currently use paths relative to their
own directories. Run each command from the directory shown below.

## Data applications

### GSE9101 microarray application

The input consists of the 12 raw CEL files for GEO accession
[GSE9101](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9101). The script
downloads the data directly from NCBI GEO.

```bash
cd "Real Dataset/GSE9101"
Rscript GSE9101.R
cd ../..
```

### ChIP-seq application

The normalized H3K4me3 inputs for GM12891 and GM12892 come from the
[MAnorm2 repository](https://github.com/tushiqi/MAnorm2). The script downloads
the files from the commit recorded in `reproducibility.yml`.

```bash
cd "Real Dataset/chip-seq"
Rscript chipseq.R
cd ../..
```

### PBMC1 pseudobulk application

The PBMC pseudobulk input comes from the
[voomByGroup repository](https://github.com/YOU-k/voomByGroup). The script
downloads the RDS file from the commit recorded in `reproducibility.yml`.

```bash
cd "Real Dataset/PBMC1"
Rscript PBMC1.R
cd ../..
```

## Simulations

Run simulations from the `Simulation/` directory. The required arguments are:

- `--var_struct`: `0` for unequal variances, `1` for equal variances, or `2`
  for diffuse variance ratios;
- `--KA` and `--KB`: group sample sizes;
- `--rounds`: number of simulation replicates;
- `--num_workers`: number of parallel workers.

A one-replicate smoke test is:

```bash
cd Simulation
Rscript simulation.R --var_struct 0 --KA 5 --KB 5 --rounds 1 --num_workers 1
```

The paper results use 50 replicates for each point. The required grids can be
generated with:

```bash
run_grid () {
  scenario="$1"
  ka="$2"
  for kb in $(seq "$ka" "$((ka + 6))"); do
    Rscript simulation.R \
      --var_struct "$scenario" --KA "$ka" --KB "$kb" \
      --rounds 50 --num_workers 1
  done
}

run_grid 0 3
run_grid 0 5
run_grid 1 3
run_grid 1 5
run_grid 2 5
cd ..
```

Increase `--num_workers` to run replicates in parallel. Results are written to
`Simulation/Simulation_result/{unequal,equal,diffuse}/`. Each replicate uses a
deterministic seed derived from the base seed and replicate number.

## Paper plots

Paper plots are rendered separately from model fitting. First generate or place
the required application and simulation result CSVs in their `Result/` folders,
then run from the repository root:

```bash
Rscript paper_plots/render_applications.R
Rscript paper_plots/render_pbmc1.R
Rscript paper_plots/render_simulations.R
```

The generated vector PDFs are written to `paper_plots/figures/`:

- `macrophage_pair.pdf` and `chipseq_pair.pdf`;
- `pbmc1_pair.pdf`;
- `overview.pdf` and `unequal_ka3.pdf`.

## Repository layout

```text
EPB main/                 Core and weighted TwoSampleEPB implementations
helper functions/         VoombyGroup helper
Real Dataset/             GSE9101, ChIP-seq, and PBMC1 applications
Simulation/               Simulation driver and replicate-level results
paper_plots/              Standalone paper-figure rendering scripts
reproducibility.yml       Tested versions, data provenance, and checksums
```

## Reproducibility notes

`reproducibility.yml` records the tested R, Bioconductor, R-package, MOSEK, and
TeX versions; application data sources, pinned revisions, and checksums;
expected application results; and the simulation seed and settings.
