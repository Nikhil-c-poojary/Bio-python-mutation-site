from io import StringIO
from math import log2

from Bio import SeqIO
from Bio.Seq import Seq
from flask import Flask, render_template, request


app = Flask(__name__)

BASES = ("A", "C", "G", "T")
STOP_CODONS = {"TAA", "TAG", "TGA"}
TRANSITION = {
    "A": "G",
    "G": "A",
    "C": "T",
    "T": "C",
}
HOTSPOTS = {
    "CG": "CpG dinucleotide hotspot",
    "GC": "GC-rich local context",
    "GG": "G-run slippage context",
    "CC": "C-run slippage context",
    "AAA": "A-run slippage context",
    "TTT": "T-run slippage context",
}


def parse_fasta(raw_text):
    handle = StringIO(raw_text.strip())
    records = list(SeqIO.parse(handle, "fasta"))
    if not records:
        raise ValueError("Please provide at least one FASTA record.")
    return records


def clean_dna(sequence):
    seq = str(sequence).upper().replace("U", "T")
    allowed = set(BASES)
    invalid = sorted({base for base in seq if base not in allowed})
    if invalid:
        joined = ", ".join(invalid)
        raise ValueError(f"Only DNA bases A, C, G, T are supported. Invalid characters: {joined}")
    if len(seq) < 3:
        raise ValueError("Sequence must contain at least 3 bases.")
    return seq


def shannon_entropy(window):
    counts = {base: window.count(base) for base in BASES}
    total = len(window)
    entropy = 0
    for count in counts.values():
        if count:
            probability = count / total
            entropy -= probability * log2(probability)
    return entropy


def detect_context(seq, index):
    reasons = []
    left = max(0, index - 3)
    right = min(len(seq), index + 4)
    context = seq[left:right]
    for pattern, label in HOTSPOTS.items():
        if pattern in context:
            reasons.append(label)
    if seq[index] in {"G", "C"}:
        reasons.append("GC base")
    if 0 < index < len(seq) - 1 and seq[index - 1] == seq[index] == seq[index + 1]:
        reasons.append("homopolymer run")
    return sorted(set(reasons))


def codon_effect(seq, index, alt_base):
    codon_start = (index // 3) * 3
    if codon_start + 3 > len(seq):
        return {
            "codon": "partial",
            "mutated_codon": "partial",
            "amino_acid": "-",
            "mutated_amino_acid": "-",
            "effect": "non-coding/partial codon",
        }

    codon = seq[codon_start : codon_start + 3]
    mutated = codon[: index - codon_start] + alt_base + codon[index - codon_start + 1 :]
    aa = str(Seq(codon).translate())
    mutated_aa = str(Seq(mutated).translate())

    if aa == mutated_aa:
        effect = "silent"
    elif mutated in STOP_CODONS:
        effect = "nonsense"
    elif codon in STOP_CODONS and mutated not in STOP_CODONS:
        effect = "stop-loss"
    else:
        effect = "missense"

    return {
        "codon": codon,
        "mutated_codon": mutated,
        "amino_acid": aa,
        "mutated_amino_acid": mutated_aa,
        "effect": effect,
    }


def estimate_probability(seq, index, alt_base):
    base = seq[index]
    context = detect_context(seq, index)
    window = seq[max(0, index - 5) : min(len(seq), index + 6)]
    entropy = shannon_entropy(window)

    score = 0.08
    if TRANSITION.get(base) == alt_base:
        score += 0.18
    if "CpG dinucleotide hotspot" in context:
        score += 0.24
    if "homopolymer run" in context:
        score += 0.14
    if base in {"G", "C"}:
        score += 0.08
    if entropy < 1.35:
        score += 0.08

    return min(score, 0.92), context


def probability_label(score):
    if score >= 0.55:
        return "High"
    if score >= 0.32:
        return "Medium"
    return "Low"


def analyze_record(record):
    seq = clean_dna(record.seq)
    mutations = []
    for index, base in enumerate(seq):
        for alt_base in BASES:
            if alt_base == base:
                continue
            probability, context = estimate_probability(seq, index, alt_base)
            effect = codon_effect(seq, index, alt_base)
            mutations.append(
                {
                    "position": index + 1,
                    "reference": base,
                    "alternate": alt_base,
                    "change": f"{base}{index + 1}{alt_base}",
                    "probability": round(probability * 100, 1),
                    "label": probability_label(probability),
                    "context": ", ".join(context) or "ordinary local context",
                    **effect,
                }
            )

    mutations.sort(key=lambda item: item["probability"], reverse=True)
    gc = round(((seq.count("G") + seq.count("C")) / len(seq)) * 100, 1)
    return {
        "id": record.id,
        "description": record.description,
        "length": len(seq),
        "gc": gc,
        "sequence": seq,
        "top_mutations": mutations[:50],
        "mutation_count": len(mutations),
    }


@app.route("/", methods=["GET", "POST"])
def index():
    default_fasta = """>BRCA1_demo_fragment
ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAATGTCCTTTGTGTAAGAATG"""
    result = None
    error = None
    fasta_text = request.form.get("fasta_text", default_fasta)

    if request.method == "POST":
        uploaded = request.files.get("fasta_file")
        if uploaded and uploaded.filename:
            fasta_text = uploaded.read().decode("utf-8", errors="replace")
        try:
            records = parse_fasta(fasta_text)
            result = [analyze_record(record) for record in records]
        except ValueError as exc:
            error = str(exc)

    return render_template(
        "index.html",
        default_fasta=default_fasta,
        fasta_text=fasta_text,
        result=result,
        error=error,
    )


if __name__ == "__main__":
    app.run(debug=True)
