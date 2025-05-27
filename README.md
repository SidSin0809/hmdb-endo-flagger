# hmdb-endo-flagger
A Python toolkit to identify and score endogenous human metabolites from HMDB XML metadata

![License](https://img.shields.io/github/license/SidSin0809/hmdb-endo-flagger)

HMDB Endogenous Metabolite Flagger is a Python-based, memory-efficient engine that streams HMDB’s massive XML metabolite dump via xml.etree.iterparse, parsing one <metabolite> node at a time to keep peak RAM usage in the low tens of megabytes. For each compound, it harvests text from fields such as origin, description, biospecimen annotations, tissue locations, and ontological tags, then computes a raw score by summing weighted keyword occurrences drawn from curated positive (e.g. “endogenous”, “human”) and negative (e.g. “xenobiotic”, “plant”) dictionaries. That raw score is passed through a logistic sigmoid to yield an intuitive 0.0–1.0 confidence metric, and any metabolite exceeding a user-configurable threshold (default 0.90) is flagged as human-endogenous in the output TSV.

## Features

- Parses HMDB XML (or zipped XML) records in a streaming fashion  
- Extracts biochemical metadata (formula, weights, InChI, ontology terms)  
- Applies keyword-weighted scoring with sigmoid transform for endogenous confidence  
- Flags metabolites as human-endogenous based on customizable threshold  
- Outputs a TSV with full metadata, confidence scores, and binary flags  

## Installation

```bash
# Clone the repo
git clone https://github.com/SidSin0809/hmdb-endo-flagger.git

# Install dependencies
  HMDB Endogenous Flagger
  No external dependencies; uses only Python 3 standard library (Requires Python ≥ 3.7).
``` 

## Usage
python hmdb_endo_flagger.py INPUT.XML[.zip] \
    --out endogenous_report.tsv \
    --threshold 0.90

1. INPUT.XML[.zip]: HMDB metabolite XML file (or ZIP containing it)
2. --out: path to output TSV (default: hmdb_output.tsv)
3. --threshold: confidence cutoff (0.0–1.0) for flagging endogenous metabolites

## Output Format
| Column                  | Description                              |
| ----------------------- | ---------------------------------------- |
| `HMDB_ID`               | Unique HMDB accession (e.g. HMDB0000123) |
| `Name`                  | Metabolite name                          |
| `Chemical_Formula`      | Molecular formula                        |
| …                       | …                                        |
| `Endogenous_Confidence` | Floating-point score (0.000–1.000)       |
| `Human_Endogenous_Flag` | `1` if confidence ≥ threshold, else `0`  |

## How It Works
1. Streaming parse of XML to minimize memory use
2. Keyword dictionaries of positive/negative indicators with weights
3. Aggregate text from origin, description, biospecimens, tissues, classification, ontology
4. Compute raw score by summing weights for keyword occurrences
5. Sigmoid transform to map raw score to [0,1] confidence
6. Threshold to produce binary flag
