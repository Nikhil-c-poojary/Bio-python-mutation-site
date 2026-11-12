
from flask import Flask, render_template, request
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from io import StringIO
from collections import Counter
import warnings

warnings.filterwarnings("ignore")

# -----------------------------------------
# FLASK APPLICATION
# -----------------------------------------

app = Flask(__name__)


# -----------------------------------------
# AMINO ACID INFORMATION
# -----------------------------------------

AMINO_ACIDS = {
    "A": "Alanine",
    "R": "Arginine",
    "N": "Asparagine",
    "D": "Aspartic acid",
    "C": "Cysteine",
    "E": "Glutamic acid",
    "Q": "Glutamine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "L": "Leucine",
    "K": "Lysine",
    "M": "Methionine",
    "F": "Phenylalanine",
    "P": "Proline",
    "S": "Serine",
    "T": "Threonine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
    "V": "Valine",
    "*": "Stop"
}


# -----------------------------------------
# FASTA PARSING AND VALIDATION
# -----------------------------------------

def parse_fasta(sequence_text):

    if not sequence_text or not sequence_text.strip():
        raise ValueError("Please enter a DNA FASTA sequence.")

    sequence_text = sequence_text.strip()

    if not sequence_text.startswith(">"):
        sequence_text = ">Sequence\n" + sequence_text

    try:
        records = list(
            SeqIO.parse(StringIO(sequence_text), "fasta")
        )

        if not records:
            raise ValueError("Invalid FASTA format.")

        record = records[0]

        sequence = str(record.seq).upper()
        sequence = "".join(sequence.split())

        if not sequence:
            raise ValueError("The DNA sequence is empty.")

        invalid_bases = set(sequence) - set("ATGCN")

        if invalid_bases:
            raise ValueError(
                "Invalid DNA characters found: "
                + ", ".join(sorted(invalid_bases))
            )

        return {
            "id": record.id,
            "description": record.description,
            "sequence": sequence
        }

    except ValueError:
        raise

    except Exception:
        raise ValueError("Unable to read the FASTA sequence.")


# -----------------------------------------
# GC AND AT CONTENT
# -----------------------------------------

def calculate_sequence_statistics(sequence):

    counts = Counter(sequence)

    length = len(sequence)

    gc_count = counts["G"] + counts["C"]
    at_count = counts["A"] + counts["T"]

    gc_content = (
        gc_count / length * 100 if length else 0
    )

    at_content = (
        at_count / length * 100 if length else 0
    )

    return {
        "length": length,
        "gc": round(gc_content, 2),
        "at": round(at_content, 2),
        "gc_count": gc_count,
        "at_count": at_count,
        "a_count": counts["A"],
        "t_count": counts["T"],
        "g_count": counts["G"],
        "c_count": counts["C"],
        "n_count": counts["N"]
    }


# -----------------------------------------
# CODON AND AMINO ACID ANALYSIS
# -----------------------------------------

def translate_codon(sequence, position):

    start = ((position - 1) // 3) * 3

    codon = sequence[start:start + 3]

    if len(codon) != 3 or "N" in codon:
        return codon, "Unknown", "X"

    amino_acid_code = str(Seq(codon).translate())

    amino_acid_name = AMINO_ACIDS.get(
        amino_acid_code,
        "Unknown"
    )

    return codon, amino_acid_name, amino_acid_code


def analyse_codon(reference, sample, ref_position, sample_position):

    ref_codon, ref_aa, ref_code = translate_codon(
        reference, ref_position
    )

    sample_codon, sample_aa, sample_code = translate_codon(
        sample, sample_position
    )

    if len(ref_codon) != 3 or len(sample_codon) != 3:
        effect = "Incomplete codon"

    elif "N" in ref_codon or "N" in sample_codon:
        effect = "Unknown"

    elif ref_code == sample_code:
        effect = "Synonymous"

    elif sample_code == "*":
        effect = "Nonsense"

    elif ref_code == "*":
        effect = "Stop-loss"

    else:
        effect = "Nonsynonymous"

    return {
        "codon_change": f"{ref_codon} → {sample_codon}",
        "ref_aa": ref_aa,
        "sample_aa": sample_aa,
        "protein_change": f"{ref_aa} → {sample_aa}",
        "effect": effect
    }


# -----------------------------------------
# DNA SEQUENCE ALIGNMENT
# -----------------------------------------

def align_sequences(reference, sample):

    alignments = pairwise2.align.globalms(
        reference,
        sample,
        2,
        -1,
        -5,
        -0.5,
        one_alignment_only=True
    )

    if not alignments:
        raise ValueError(
            "Unable to align the DNA sequences."
        )

    alignment = alignments[0]

    return (
        alignment.seqA,
        alignment.seqB,
        alignment.score
    )


# -----------------------------------------
# MUTATION DETECTION
# -----------------------------------------

def detect_mutations(
    reference,
    sample,
    reference_aligned,
    sample_aligned
):

    mutations = []
    aligned_bases = []

    reference_position = 0
    sample_position = 0

    substitution_count = 0
    insertion_count = 0
    deletion_count = 0

    for ref_base, sample_base in zip(
        reference_aligned,
        sample_aligned
    ):

        if ref_base != "-":
            reference_position += 1

        if sample_base != "-":
            sample_position += 1

        is_match = ref_base == sample_base

        aligned_bases.append({
            "reference": ref_base,
            "sample": sample_base,
            "match": is_match
        })

        if is_match:
            continue

        if ref_base == "-":

            mutation_type = "Insertion"
            insertion_count += 1

            position = max(1, reference_position)

            codon_info = {
                "codon_change": "Not applicable",
                "ref_aa": "Unknown",
                "sample_aa": "Unknown",
                "protein_change": "Not determined",
                "effect": "Requires frameshift analysis"
            }

        elif sample_base == "-":

            mutation_type = "Deletion"
            deletion_count += 1

            position = reference_position

            codon_info = {
                "codon_change": "Not applicable",
                "ref_aa": "Unknown",
                "sample_aa": "Unknown",
                "protein_change": "Not determined",
                "effect": "Requires frameshift analysis"
            }

        else:

            mutation_type = "Substitution"
            substitution_count += 1

            position = reference_position

            codon_info = analyse_codon(
                reference,
                sample,
                reference_position,
                sample_position
            )

        position_percent = round(
            position / len(reference) * 100,
            2
        ) if reference else 0

        mutation = {
            "position": position,
            "position_percent": position_percent,
            "change": f"{ref_base} → {sample_base}",
            "type": mutation_type,
            **codon_info
        }

        mutations.append(mutation)

    alignment_text = (
        f"Reference: {reference_aligned}\n"
        f"           "
        + "".join(
            "|" if a == b else " "
            for a, b in zip(
                reference_aligned,
                sample_aligned
            )
        )
        + f"\nSample:    {sample_aligned}"
    )

    return {
        "mutations": mutations,
        "aligned_bases": aligned_bases,
        "reference_aligned": reference_aligned,
        "sample_aligned": sample_aligned,
        "alignment": alignment_text,
        "substitution_count": substitution_count,
        "insertion_count": insertion_count,
        "deletion_count": deletion_count
    }


# -----------------------------------------
# COMPLETE SEQUENCE ANALYSIS
# -----------------------------------------

def analyse_sequences(reference_text, sample_text):

    reference_record = parse_fasta(reference_text)
    sample_record = parse_fasta(sample_text)

    reference = reference_record["sequence"]
    sample = sample_record["sequence"]

    reference_stats = calculate_sequence_statistics(
        reference
    )

    sample_stats = calculate_sequence_statistics(
        sample
    )

    (
        reference_aligned,
        sample_aligned,
        alignment_score
    ) = align_sequences(reference, sample)

    mutation_results = detect_mutations(
        reference,
        sample,
        reference_aligned,
        sample_aligned
    )

    mutations = mutation_results["mutations"]

    effect_counts = Counter(
        mutation["effect"]
        for mutation in mutations
    )

    mutation_type_counts = Counter(
        mutation["type"]
        for mutation in mutations
    )

    result = {

        "id": sample_record["id"],
        "description": sample_record["description"],
        "reference_id": reference_record["id"],

        "reference_length": len(reference),
        "length": len(sample),

        "sequence": sample,

        "gc": sample_stats["gc"],
        "at": sample_stats["at"],

        "gc_count": sample_stats["gc_count"],
        "at_count": sample_stats["at_count"],

        "a_count": sample_stats["a_count"],
        "t_count": sample_stats["t_count"],
        "g_count": sample_stats["g_count"],
        "c_count": sample_stats["c_count"],
        "n_count": sample_stats["n_count"],

        "reference_gc": reference_stats["gc"],
        "reference_at": reference_stats["at"],

        "reference_aligned": reference_aligned,
        "sample_aligned": sample_aligned,

        "aligned_bases": mutation_results["aligned_bases"],
        "aligned_length": len(reference_aligned),

        "alignment": mutation_results["alignment"],
        "alignment_score": round(alignment_score, 2),

        "mutation_count": len(mutations),

        "substitution_count": mutation_results[
            "substitution_count"
        ],

        "insertion_count": mutation_results[
            "insertion_count"
        ],

        "deletion_count": mutation_results[
            "deletion_count"
        ],

        "mutation_type_counts": dict(mutation_type_counts),
        "effect_counts": dict(effect_counts),

        "top_mutations": mutations[:25],

        "codon_mutations": [
            mutation
            for mutation in mutations
            if mutation["type"] == "Substitution"
        ][:25]
    }

    return [result]


# -----------------------------------------
# FLASK ROUTE
# -----------------------------------------

@app.route("/", methods=["GET", "POST"])
def index():

    result = None
    error = None

    if request.method == "POST":

        reference_sequence = request.form.get(
            "reference_sequence", ""
        )

        sample_sequence = request.form.get(
            "sample_sequence", ""
        )

        try:

            result = analyse_sequences(
                reference_sequence,
                sample_sequence
            )

        except ValueError as e:

            error = str(e)

        except Exception:

            error = (
                "An unexpected error occurred during analysis. "
                "Please check your FASTA sequences and try again."
            )

    return render_template(
        "index.html",
        result=result,
        error=error
    )


# -----------------------------------------
# LOCAL EXECUTION
# -----------------------------------------

if __name__ == "__main__":

    app.run(
        host="0.0.0.0",
        port=5000,
        debug=True
    )
