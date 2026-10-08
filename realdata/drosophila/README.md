# DGRP analysis (reviewer comment 7, minors 3 and 7)

| file | what it does |
|---|---|
| `dgrp_ascend.R` | the original analysis of the paper (network, hubs, Manhattan plots); unchanged apart from paths |
| `dgrp_data.R` | the same filtering as `dgrp_ascend.R` (same 500 SNPs, same 250 genes), cached once as `dgrp_input.rds`, plus the counts of every filtering step |
| `dgrp_controls.R` | one ASCEND run of the controls: the observed network, 100 genotype permutations, 20 random tier assignments, 100 subsamples of 80% of lines, 20 random gene sets (241 runs) |
| `dgrp_controls_summarise.R` | every number of the paper's DGRP controls paragraph, in `DGRP_CONTROLS.md`, with CSV tables and two figures |
| `slurm/submit_dgrp.sh` | submits everything on CREATE |

## Data (not in git)

Put these two files in one folder, for example `/scratch/users/$USER/dgrp`:

- `GSE117850_DGRP_GEO_Table_4_Gene_Male_Line_Means.txt.gz` from GEO series GSE117850
- `dgrp2.tgeno.txt` from `data.tar.gz` at https://zenodo.org/records/14871341

## Run on CREATE

```bash
cd Ascend                                    # the unpacked repository
export DGRP_DATA=/scratch/users/$USER/dgrp   # the folder with the two files
bash realdata/drosophila/slurm/submit_dgrp.sh
```

This runs the observed network first (it also caches the filtered data), then the 240 control runs as an array (4 runs per task), then the summary. Results are in `$DGRP_DATA/results`: open `DGRP_CONTROLS.md` first.

For gene symbols and the GO enrichment, install `org.Dm.eg.db` once:
`Rscript -e 'BiocManager::install("org.Dm.eg.db")'`. Without it the enrichment is skipped (or give your own list of FlyBase IDs with `DGRP_GENESET=file`).

## Check before using the numbers

- `DGRP_CONTROLS.md`, "Observed network": the directed and undirected counts and density must match the paper (880, 61, 3.02%). If they differ, the input files or the filtering differ from the original run.
- The flow table gives the SNP and gene counts for Supplementary Fig. S8.

## Test without the real data

```bash
Rscript tests/make_fake_dgrp.R /tmp/fake_dgrp
DGRP_DATA=/tmp/fake_dgrp DGRP_OUT=/tmp/fake_dgrp/results Rscript realdata/drosophila/dgrp_controls.R 1
DGRP_DATA=/tmp/fake_dgrp DGRP_OUT=/tmp/fake_dgrp/results Rscript realdata/drosophila/dgrp_controls.R 2
DGRP_OUT=/tmp/fake_dgrp/results DGRP_GENESET=/tmp/fake_dgrp/immune_genes.txt Rscript realdata/drosophila/dgrp_controls_summarise.R
```
