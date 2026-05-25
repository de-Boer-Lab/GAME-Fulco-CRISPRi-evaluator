## Evaluator Data

This folder contains the raw input data and processed files used by the Fulco CRISPRi evaluator. 

---

## Files

### Raw inputs (from the Fulco et al. 2019 supplementary data)

| File | Description |
|---|---|
| `41588_2019_538_MOESM3_ESM.xlsx` | Original supplementary data file from the publication. |
| `Supplementary Table 3a. Data for E-G pairs.csv` | 4,914 enhancer–gene (E–G) pairs with coordinates, measured fractional change in gene expression, significance, and effect-size statistics. |
| `Supplementary Table 5b. Human gene annotations from RefSeq collapsed to a single consensus annotation per gene.csv` | Consensus RefSeq annotations (chr, start, end, gene, strand) for the 30 tested genes. |

### Processed files (produced by `Parse_Fulco_inputData.py`)

| File | Description |
|---|---|
| `fulco_evaluator_coordinates.parquet` | 4,882 E-G pairs with element coordinates (`chr`, `start`, `end`), measured fractional expression change, target gene coordinates (`gene_chr`, `gene_start`, `gene_end`, `strand`), distance from element to gene, and the 4 Mb sequence window (`sequence_start`, `sequence_end`) centered on the gene's start position. The drop from 4,914 to 4,882 reflects E-G pairs whose elements fall outside the ±2 Mb window. |
| `fulco_reference_gene_sequence.parquet` | One reference sequence per tested gene (30 rows). Each sequence is 4,000,000 bp, pulled from `gene_start − 2,000,000` to `gene_start + 2,000,000` on hg19. Stored as `Gene`, `sequence`. |

### Scripts and plots

| File | Description |
|---|---|
| `Parse_Fulco_inputData.py` | Generates the two parquet files above from the supplementary tables. Requires `pysam` and a local hg19 FASTA (`hg19.fa`). |
| `eg_scatter.png` | Diagnostic scatter plot of element vs. fractional expression change, colored by target gene. |

---

## Reference Genome

All coordinates and reference sequences use **hg19 / GRCh37**. The parse script reads sequences from a local FASTA at `hg19.fa` via `pysam`; update the `fasta_path` variable at the top of `Parse_Fulco_inputData.py` to point to your copy.

---

## Regenerating the Processed Files

```bash
python Parse_Fulco_inputData.py
```

This will rebuild `fulco_evaluator_coordinates.parquet` and `fulco_reference_gene_sequence.parquet` from the two supplementary CSVs. 
