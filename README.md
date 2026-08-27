# Trowbridgii ddRAD population genomics

Code and final analysis datasets for the reference-mapped ddRAD study of
*Sorex trowbridgii*. Reads were mapped to the *Sorex ornatus* draft reference
genome (NCBI assembly GCA_041430635.1). The final downstream dataset excludes
`MVZ216210`, which was identified as an extreme missing-data outlier.

## Repository contents

- `code/`: publication-release scripts and analysis configuration files,
  grouped by analytical method.
- `data/metadata/`: population assignments and normalized geographic
  coordinates.
- `data/stacks_exports/`: final STACKS SNP exports in VCF, STRUCTURE, PLINK,
  PHYLIP, Treemix, and summary-statistic formats.
- `data/unlinked_snps/`: PLINK dataset containing one SNP per locus.
- `data/phylogeny/`: IQ-TREE alignment, maximum-likelihood trees, run report,
  and branch-support table.
- `data/population_structure/`: best-K ancestry/membership matrices and
  model-selection tables.
- `data/species_delimitation/`: reported SPEEDEMON, delimitR, and BPP inputs
  and compact outputs.
- `data/eems/`: self-contained EEMS inputs, fixed grid, final surface objects,
  and map scale ranges.
- `results/`: compact final quality-control, population-statistic, SFS, IBD,
  EEMS, and model-selection tables.
- `docs/`: data dictionary, reproducibility notes, and release checklist.

## Core dataset

| Dataset | Individuals | Sites/markers |
|---|---:|---:|
| Final filtered STACKS VCF | 99 | 23,482 variants |
| One-SNP-per-locus PLINK data | 99 | 2,557 SNPs |
| Biallelic IQ-TREE alignment | 99 | 11,584 sites |
| Coordinate-based analyses | 98 | 14,578 filtered EEMS SNPs |

The six a priori population labels are North, North_Coast, Sierra_1,
Sierra_2, Sierra_3, and South_Coast.

## Reproducing analyses

The code preserves the original parameter settings, random seeds, and input
filenames. Many scripts also retain the absolute project paths used during the
original analysis. Before execution on another computer, replace the project
root near the beginning of each script or recreate the expected directory
layout. See `code/README.md` for the recommended workflow order and
`code/SOFTWARE_REQUIREMENTS.md` for dependencies.

Compressed text files can be read without permanent decompression, for
example:

```bash
gzip -dc data/stacks_exports/populations.snps.vcf.gz | less
```

Run the release audit after cloning:

```bash
bash validate_repository.sh
```


## Citation and licensing

The final article citation will be provided upon publication.

