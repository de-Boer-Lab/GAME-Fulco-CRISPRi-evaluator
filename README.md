# Fulco CRISPRi Evaluator

To evaluate models’ abilities to predict long range interactions for synthetic modifications we included the Fulco et al. CRISPRi dataset (https://www.nature.com/articles/s41588-019-0538-0#Sec19). Endogenous enhancers are synthetically repressed and their effects on 30 unique genes are measured in K562 cells.

---

## How It Works

For all 30 genes in the dataset, the Evaluator:

1. Builds a request containing the reference gene sequence and one shuffled-enhancer variant per enhancer–gene pair
2. Negotiates serialization format (JSON or MsgPack) with the predictor
3. POSTs sequences to the Predictor's `/predict` endpoint
4. Computes the predicted fractional change in expression: `(predicted - reference) / reference`
5. Calculates Pearson *r* against measured CRISPRi fractional changes
6. Saves per-gene metrics and a final averaged Pearson *r* across all genes


### Run Evaluator using Apptainer container

Download container and evaluator data from Hugging Face: https://huggingface.co/datasets/deBoerLab/FulcoCRISPRi

```bash
apptainer run --containall \
  -B /absolute/path/to/evaluator_data:/evaluator_data \
  -B /absolute/path/to/predictions:/predictions \
  fulco_evaluator.sif <predictor_ip> <predictor_port> /predictions
```

---

## Arguments

| Argument | Description |
|---|---|
| `predictor_ip` | IP address or hostname of the predictor server |
| `predictor_port` | Port the predictor is listening on |
| `output_dir` | Directory where prediction JSONs and metric CSVs are written |

---

## Request Structure

One POST request is sent per gene (30 total), due to large sequence sizes. Each request contains:

- **One reference sequence** — a 4 Mb window centered on the gene's TSS (±2 Mb) `fulco_reference_gene_sequence.parquet`
- **One shuffled-enhancer sequence per enhancer–gene pair** — The file `fulco_evaluator_coordinates.parquet` contains the coordinates for each enhancer. The Evaluator `data_loader.py` code will shuffle the bases of the target enhancer in-place to represent "inactivation" (the rest of the sequence remains unchanged). A fixed random seed (`seed=42`) is used for reproducibility.

**Strand handling:** If the gene is on the negative strand, the reverse complement of the full 4 Mb window is taken — both for the reference sequence and for each shuffled-enhancer sequence — after the enhancer shuffling is applied. This ensures sequences are always presented in the 5'→3' direction of the gene.

`prediction_ranges` specifies the gene body coordinates relative to the sequence window for each sequence (2Mbp -> 2Mbp + gene length). This requests predictions at the gene in its native context or with the effect of the mutated enhancer (shuffled).


Each request payload follows this structure:

```json
{
  "readout": "point",
  "prediction_tasks": [
    {
      "name": "GENE",
      "type": "expression",
      "cell_type": "K562",
      "scale": "linear",
      "species": "homo_sapiens"
    }
  ],
  "sequences": {
    "GENE_reference_seq": "<4Mb sequence>",
    "GENE_element_seq": "<4Mb sequence with shuffled enhancer>",
    ...
  },
  "prediction_ranges": {
    "GENE_reference_seq": [start,end],
    "GENE_element_seq": [start,end],
    ...
  }
}
```

---

## Outputs

All outputs are written to `output_dir`:

| File | Description |
|---|---|
| `Fulco_CRISPRi_<Gene>.json` | Raw predictions returned by the predictor for each gene |
| `correlation_summary_Fulco_CRISPRi_<Gene>.csv` | Per-gene Pearson *r* |
| `correlation_summary_Fulco_CRISPRi_final.csv` | Average Pearson *r* across all genes |

---

## Directory Structure

```
Fulco_CRISPRi/
├── Fulco_evaluator.py                  # Main entry point
├── config.py                           # Settings (name, paths, formats, retries)
├── data_loader.py                      # Loads parquet files, builds request payload
├── evaluator_content_handler.py        # Format negotiation, HTTP POST, deserialization
├── evaluator_metrics_calculator.py     # Pearson r calculation and CSV output
├── fulco_evaluator.def                 # Apptainer container definition
└── evaluator_data/
    ├── fulco_reference_gene_sequence.parquet
    └── fulco_evaluator_coordinates.parquet
```

---


