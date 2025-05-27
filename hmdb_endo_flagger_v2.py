#!/usr/bin/env python3
"""

Detection of endogenous human metabolites using multi-field text mining
and confidence scoring based solely on HMDB XML metadata.

Usage:
    python hmdb_endo_flagger.py INPUT.XML[.zip] --out OUTPUT.tsv [--threshold 0.85]

Outputs TSV with columns:
  HMDB_ID, Name, Chemical_Formula, Average_Molecular_Weight,
  Monoisotopic_Molecular_Weight, SMILES, InChI, InChIKey,
  Status, Origin, Ontology_Term_Count, All_Ontology_Terms,
  Description, Biospecimens, Tissues, Classification,
  Endogenous_Confidence (0.0-1.0), Human_Endogenous_Flag (0/1)
"""
import argparse
import csv
import math
import sys
import zipfile
from pathlib import Path
from xml.etree.ElementTree import iterparse

# Define indicator keywords and weights
POSITIVE_KEYWORDS = {
    "endogenous": 2,
    "biosynthesized": 2,
    "de novo": 1,
    "present in human": 3,
    "human": 1,
    "homo sapiens": 3,
    "mammalian": 1,
    "blood": 2,
    "plasma": 2,
    "serum": 2,
    "urine": 2,
    "csf": 2,
    "tissue": 1,
    "liver": 1,
    "kidney": 1,
    "brain": 1,
    "metabolized by": 2,
    "synthesized": 2,
    "biosynthesis": 2
}
NEGATIVE_KEYWORDS = {
    "exogenous": 2,
    "xenobiotic": 2,
    "pharmaceutical": 2,
    "drug": 2,
    "toxin": 1,
    "toxic": 2,
    "pesticide": 2,
    "bacteria": 1,
    "fungus": 2,
    "dietary": 1,
    "plant": 1,
    "microbial": 1,
    "environmental": 1
}

THRESHOLD_DEFAULT = 0.85


def local_name(tag: str) -> str:
    return tag.rsplit('}', 1)[-1] if '}' in tag else tag


def extract_text(elem, tag_name):
    node = elem.find(f'.//{{*}}{tag_name}')
    return (node.text or '').strip().lower() if node is not None and node.text else ''


def collect_terms(elem):
    return {node.text.strip().lower()
            for node in elem.iter()
            if local_name(node.tag).lower() == 'term' and node.text}


def compute_confidence(text: str) -> float:
    score = 0.0
    for kw, w in POSITIVE_KEYWORDS.items():
        if kw in text:
            score += w
    for kw, w in NEGATIVE_KEYWORDS.items():
        if kw in text:
            score -= w
    # sigmoid transform
    try:
        confidence = 1 / (1 + math.exp(-score))
    except OverflowError:
        confidence = 1.0 if score > 0 else 0.0
    return confidence


def extract_metabolites(fileobj):
    for _, elem in iterparse(fileobj, events=('end',)):
        if local_name(elem.tag).lower() != 'metabolite':
            continue

        # Basic fields
        mid = extract_text(elem, 'accession').upper()
        if not mid:
            elem.clear()
            continue
        name = extract_text(elem, 'name')
        chem_formula = extract_text(elem, 'chemical_formula')
        avg_wt = extract_text(elem, 'average_molecular_weight')
        mono_wt = extract_text(elem, 'monisotopic_molecular_weight')
        smiles = extract_text(elem, 'smiles')
        inchi = extract_text(elem, 'inchi')
        inchikey = extract_text(elem, 'inchikey')
        status = extract_text(elem, 'status')
        origin = extract_text(elem, 'origin')

        # New metadata fields
        description = extract_text(elem, 'description')
        biospecimens = ' '.join(
            (node.text or '').strip().lower()
            for node in elem.findall('.//{*}biospecimen') if node.text)
        tissues = ' '.join(
            (node.text or '').strip().lower()
            for node in elem.findall('.//{*}tissue') if node.text)
        classification = ' '.join([
            extract_text(elem, 'super_class'),
            extract_text(elem, 'sub_class')
        ])

        # Ontology terms
        terms = collect_terms(elem)

        # Build combined text
        combined_text = ' '.join([
            origin, description, biospecimens, tissues, classification,
            ' '.join(terms)
        ])

        # Compute confidence & flag
        conf = compute_confidence(combined_text)
        flag = int(conf >= THRESHOLD_DEFAULT)

        yield {
            'HMDB_ID': mid,
            'Name': name,
            'Chemical_Formula': chem_formula,
            'Average_Molecular_Weight': avg_wt,
            'Monoisotopic_Molecular_Weight': mono_wt,
            'SMILES': smiles,
            'InChI': inchi,
            'InChIKey': inchikey,
            'Status': status,
            'Origin': origin,
            'Ontology_Term_Count': len(terms),
            'All_Ontology_Terms': ';'.join(sorted(terms)),
            'Description': description,
            'Biospecimens': biospecimens,
            'Tissues': tissues,
            'Classification': classification,
            'Endogenous_Confidence': f"{conf:.3f}",
            'Human_Endogenous_Flag': flag
        }

        elem.clear()


def open_input(path: Path):
    if path.suffix.lower() == '.zip':
        zf = zipfile.ZipFile(path)
        xml_files = [n for n in zf.namelist() if n.lower().endswith('.xml')]
        if not xml_files:
            raise RuntimeError(f"No .xml inside {path}")
        return zf.open(xml_files[0])
    return path.open('rb')


def main():
    parser = argparse.ArgumentParser(
        description="HMDB endogenous human metabolite detector v10"
    )
    parser.add_argument('input', type=Path, help='HMDB XML or .zip')
    parser.add_argument('--out', type=Path,
                        default=Path('hmdb_output.tsv'),
                        help='Output TSV (default: hmdb_output.tsv)')
    parser.add_argument('--threshold', type=float,
                        default=THRESHOLD_DEFAULT,
                        help='Confidence threshold for flag (0-1)')
    args = parser.parse_args()

    fields = [
        'HMDB_ID', 'Name', 'Chemical_Formula',
        'Average_Molecular_Weight', 'Monoisotopic_Molecular_Weight',
        'SMILES', 'InChI', 'InChIKey', 'Status', 'Origin',
        'Ontology_Term_Count', 'All_Ontology_Terms',
        'Description', 'Biospecimens', 'Tissues', 'Classification',
        'Endogenous_Confidence', 'Human_Endogenous_Flag'
    ]

    # Open output with UTF-8 encoding to support all characters
    with args.out.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        count = 0
        try:
            with open_input(args.input) as inp:
                for rec in extract_metabolites(inp):
                    writer.writerow(rec)
                    count += 1
                    if count % 100 == 0:
                        print(f"Processed {count} metabolites...",
                              file=sys.stderr)
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    print(f"Done: {count} entries written to {args.out}")


if __name__ == '__main__':
    main()
