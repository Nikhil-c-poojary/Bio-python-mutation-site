# Gene Mutation Possibility Analyzer

A small Flask website that uses Biopython to parse FASTA DNA sequences and rank possible point mutations.

## Run

Install Python 3.10+ first if `python --version` does not work.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python app.py
```

Open `http://127.0.0.1:5000`.

## Run With Docker

```powershell
docker build -t mutation-analyzer .
docker run --rm -p 5000:5000 mutation-analyzer
```

Open `http://127.0.0.1:5000`.

## Test

```powershell
pip install -r requirements.txt pytest
pytest
```

## Features

- Paste FASTA text or upload `.fa`, `.fasta`, `.fna`, or `.txt`.
- Parses FASTA records with `Bio.SeqIO`.
- Tests every base against the three possible point mutations.
- Shows mutation point, reference/alternate base, codon change, amino-acid change, and effect.
- Ranks mutation possibility using transition likelihood, GC/CpG context, local entropy, and homopolymer runs.

The scoring is educational and heuristic. It is not a clinical or research-grade pathogenicity model.
