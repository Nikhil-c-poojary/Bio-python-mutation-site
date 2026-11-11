from flask import Flask, render_template, request
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from io import StringIO


# ==========================================
# FLASK APPLICATION
# ==========================================

app = Flask(__name__)


# ==========================================
# FASTA PARSING
# ==========================================

def parse_fasta(sequence_text):

    if not sequence_text or not sequence_text.strip():
        raise ValueError("Please enter a FASTA sequence.")

    try:
        records = list(SeqIO.parse(
            StringIO(sequence_text.strip()), "fasta"
        ))

        if not records:
            raise ValueError(
                "Invalid FASTA format. Please include a header starting with >."
            )

        record = records[0]

        sequence = str(record.seq).upper().replace(" ", "").replace("\n", "")

        if not sequence:
            raise ValueError("The FASTA sequence is empty.")

        invalid = set(sequence) - set("ATGCN")

        if invalid:
            raise ValueError(
                "Invalid DNA characters detected: "
                + ", ".join(sorted(invalid))
            )

        record.seq = Seq(sequence)

        return record

    except ValueError:
        raise

    except Exception:
        raise ValueError("Unable to read the FASTA sequence.")


# ==========================================
# GC CONTENT
# ==========================================

def calculate_gc(sequence):

    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")

    return round((gc_count / len(sequence)) * 100, 2)


# ==========================================
# AMINO ACID ANALYSIS
# ==========================================

def amino_acid_name(codon):

    if len(codon) != 3 or "-" in codon or "N" in codon:
        return "Unknown"

    try:
        amino = str(Seq(codon).translate())

        names = {
            "A": "Alanine (Ala)",
            "R": "Arginine (Arg)",
            "N": "Asparagine (Asn)",
            "D": "Aspartic acid (Asp)",
            "C": "Cysteine (Cys)",
            "Q": "Glutamine (Gln)",
            "E": "Glutamic acid (Glu)",
            "G": "Glycine (Gly)",
            "H": "Histidine (His)",
            "I": "Isoleucine (Ile)",
            "L": "Leucine (Leu)",
            "K": "Lysine (Lys)",
            "M": "Methionine (Met)",
            "F": "Phenylalanine (Phe)",
            "P": "Proline (Pro)",
            "S": "Serine (Ser)",
            "T": "Threonine (Thr)",
            "W": "Tryptophan (Trp)",
            "Y": "Tyrosine (Tyr)",
            "V": "Valine (Val)",
            "*": "Stop codon"
        }

        return names.get(amino, "Unknown")

    except Exception:
        return "Unknown"


# ==========================================
# MUTATION CLASSIFICATION
# ==========================================

def analyse_codon(reference, sample, position):

    result = {
        "codon_change": "Not applicable",
        "ref_aa": "Not analysed",
        "sample_aa": "Not analysed",
        "protein_change": "Not analysed",
        "effect": "Not classified"
    }

    # Amino acid analysis is performed for substitutions
    # using the original reference reading frame.

    codon_start = ((position - 1) // 3) * 3

    ref_codon = reference[codon_start:codon_start + 3]
    sample_codon = sample[codon_start:codon_start + 3]

    if (
        len(ref_codon) != 3
        or len(sample_codon) != 3
        or "N" in ref_codon
        or "N" in sample_codon
    ):
        return result

    ref_aa = str(Seq(ref_codon).translate())
    sample_aa = str(Seq(sample_codon).translate())

    result["codon_change"] = f"{ref_codon} → {sample_codon}"

    result["ref_aa"] = amino_acid_name(ref_codon)

    result["sample_aa"] = amino_acid_name(sample_codon)

    result["protein_change"] = f"{ref_aa} → {sample_aa}"

    if ref_aa == sample_aa:
        result["effect"] = "Synonymous"

    elif sample_aa == "*":
        result["effect"] = "Stop gained"

    elif ref_aa == "*":
        result["effect"] = "Stop lost"

    else:
        result["effect"] = "Nonsynonymous"

    return result


# ==========================================
# MUTATION DETECTION
# ==========================================

def detect_mutations(reference, sample):

    alignments = pairwise2.align.globalms(
        reference,
        sample,
        2,
        -1,
        -5,
        -0.5
    )

    if not alignments:
        raise ValueError("Unable to align the sequences.")

    aligned = alignments[0]

    ref_aligned = aligned.seqA
    sample_aligned = aligned.seqB

    mutations = []

    ref_position = 0

    for ref_base, sample_base in zip(
        ref_aligned,
        sample_aligned
    ):

        if ref_base != "-":
            ref_position += 1

        if ref_base == sample_base:
            continue

        if ref_base != "-" and sample_base != "-":

            mutation_type = "Substitution"

            change = f"{ref_base} → {sample_base}"

            codon_info = analyse_codon(
                reference,
                sample,
                ref_position
            )

        elif ref_base == "-":

            mutation_type = "Insertion"

            change = f"- → {sample_base}"

            codon_info = {
                "codon_change": "Not analysed",
                "ref_aa": "Not analysed",
                "sample_aa": "Not analysed",
                "protein_change": "Not analysed",
                "effect": "Not classified"
            }

        else:

            mutation_type = "Deletion"

            change = f"{ref_base} → -"

            codon_info = {
                "codon_change": "Not analysed",
                "ref_aa": "Not analysed",
                "sample_aa": "Not analysed",
                "protein_change": "Not analysed",
                "effect": "Not classified"
            }

        mutation = {
            "change": change,
            "position": ref_position,
            "type": mutation_type,
            "label": "Detected",

            **codon_info
        }

        mutation["position_percent"] = round(
            (ref_position / max(len(reference), 1)) * 100,
            2
        )

        mutations.append(mutation)

    alignment_text = (
        ref_aligned
        + "\n"
        + "".join(
            "|" if a == b else " "
            for a, b in zip(ref_aligned, sample_aligned)
        )
        + "\n"
        + sample_aligned
    )

    return (
        mutations,
        alignment_text,
        aligned.score
    )


# ==========================================
# MAIN ANALYSIS
# ==========================================

def analyse_sequences(sample_text, reference_text):

    sample_record = parse_fasta(sample_text)

    reference_record = parse_fasta(reference_text)

    sample = str(sample_record.seq)

    reference = str(reference_record.seq)

    mutations, alignment_text, score = detect_mutations(
        reference,
        sample
    )

    record = {

        "id": sample_record.id,

        "description": sample_record.description,

        "reference_id": reference_record.id,

        "reference_length": len(reference),

        "length": len(sample),

        "gc": calculate_gc(sample),

        "sequence": sample,

        "alignment": alignment_text,

        "alignment_score": round(score, 2),

        "mutation_count": len(mutations),

        "top_mutations": mutations[:25]
    }

    return [record]


# ==========================================
# HOME PAGE
# ==========================================

@app.route("/", methods=["GET", "POST"])
def index():

    result = None

    error = None

    if request.method == "POST":

        sample_sequence = request.form.get(
            "sample_sequence", ""
        )

        reference_sequence = request.form.get(
            "reference_sequence", ""
        )

        try:

            result = analyse_sequences(
                sample_sequence,
                reference_sequence
            )

        except ValueError as e:

            error = str(e)

        except Exception as e:

            app.logger.exception("Sequence analysis failed")

            error = (
                "An unexpected error occurred during analysis. "
                "Please check your FASTA sequences and try again."
            )

    return render_template(
        "index.html",
        result=result,
        error=error
    )


# ==========================================
# VERCEL ENTRY POINT
# ==========================================

if __name__ == "__main__":

    app.run(
        host="0.0.0.0",
        port=5000,
        debug=False
    )
