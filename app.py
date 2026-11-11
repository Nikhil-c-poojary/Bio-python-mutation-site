from flask import Flask, render_template, request
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from io import StringIO

app = Flask(__name__)

VALID_BASES = set("ATGCN")

# Amino acid names and abbreviations
AMINO_ACIDS = {
    "A": ("Alanine", "Ala"),
    "R": ("Arginine", "Arg"),
    "N": ("Asparagine", "Asn"),
    "D": ("Aspartic acid", "Asp"),
    "C": ("Cysteine", "Cys"),
    "Q": ("Glutamine", "Gln"),
    "E": ("Glutamic acid", "Glu"),
    "G": ("Glycine", "Gly"),
    "H": ("Histidine", "His"),
    "I": ("Isoleucine", "Ile"),
    "L": ("Leucine", "Leu"),
    "K": ("Lysine", "Lys"),
    "M": ("Methionine", "Met"),
    "F": ("Phenylalanine", "Phe"),
    "P": ("Proline", "Pro"),
    "S": ("Serine", "Ser"),
    "T": ("Threonine", "Thr"),
    "W": ("Tryptophan", "Trp"),
    "Y": ("Tyrosine", "Tyr"),
    "V": ("Valine", "Val"),
    "*": ("Stop codon", "Ter")
}


# Read and validate FASTA sequence
def read_fasta(content):

    if not content or not content.strip():
        raise ValueError("Please enter a FASTA sequence.")

    records = list(SeqIO.parse(StringIO(content), "fasta"))

    if not records:
        raise ValueError("Invalid FASTA format.")

    record = records[0]

    sequence = str(record.seq).upper()

    if not sequence:
        raise ValueError("FASTA sequence is empty.")

    invalid = set(sequence) - VALID_BASES

    if invalid:
        raise ValueError(
            f"Invalid DNA characters: {', '.join(sorted(invalid))}"
        )

    record.seq = Seq(sequence)

    return record


# Translate codon into amino acid
def translate_codon(codon):

    if len(codon) != 3 or "N" in codon:
        return "Unknown", "-"

    amino_acid = str(Seq(codon).translate())

    return AMINO_ACIDS.get(
        amino_acid, ("Unknown", "-")
    )


# Analyze codon and amino acid changes
def classify_effect(reference_seq, sample_seq, ref_pos, sample_pos):

    ref_start = ((ref_pos - 1) // 3) * 3
    sample_start = ((sample_pos - 1) // 3) * 3

    ref_codon = reference_seq[ref_start:ref_start + 3]
    sample_codon = sample_seq[sample_start:sample_start + 3]

    if len(ref_codon) != 3 or len(sample_codon) != 3:

        return {
            "effect": "Unknown (incomplete codon)",
            "ref_codon": ref_codon,
            "sample_codon": sample_codon,
            "ref_aa": "-",
            "sample_aa": "-",
            "protein_change": "Not determined"
        }

    ref_name, ref_abbr = translate_codon(ref_codon)
    sample_name, sample_abbr = translate_codon(sample_codon)

    if ref_abbr == "-" or sample_abbr == "-":

        effect = "Unknown (ambiguous base)"
        protein_change = "Not determined"

    elif ref_abbr == sample_abbr:

        effect = "Synonymous"
        protein_change = f"{ref_abbr} → {sample_abbr} (No change)"

    elif sample_abbr == "Ter":

        effect = "Stop gained"
        protein_change = f"{ref_abbr} → Ter"

    elif ref_abbr == "Ter":

        effect = "Stop lost"
        protein_change = f"Ter → {sample_abbr}"

    else:

        effect = "Nonsynonymous"
        protein_change = f"{ref_abbr} → {sample_abbr}"

    return {
        "effect": effect,
        "ref_codon": ref_codon,
        "sample_codon": sample_codon,
        "ref_aa": ref_name,
        "sample_aa": sample_name,
        "protein_change": protein_change
    }


@app.route("/", methods=["GET", "POST"])
def index():

    result = None
    error = None

    if request.method == "POST":

        try:

            sample_text = request.form.get(
                "sample_sequence", ""
            )

            reference_text = request.form.get(
                "reference_sequence", ""
            )

            sample_record = read_fasta(sample_text)
            reference_record = read_fasta(reference_text)

            sample_seq = str(sample_record.seq).upper()
            reference_seq = str(reference_record.seq).upper()

            # Global sequence alignment
            alignments = pairwise2.align.globalms(
                reference_seq,
                sample_seq,
                2,
                -1,
                -5,
                -0.5
            )

            if not alignments:
                raise ValueError("Unable to align sequences.")

            alignment = alignments[0]

            aligned_ref = alignment.seqA
            aligned_sample = alignment.seqB

            alignment_visual = []

            for a, b in zip(aligned_ref, aligned_sample):

                if a == b:
                    alignment_visual.append("|")
                else:
                    alignment_visual.append(" ")

            alignment_string = (
                aligned_ref + "\n"
                + "".join(alignment_visual) + "\n"
                + aligned_sample
            )

            mutations = []

            ref_pos = 0
            sample_pos = 0

            has_indels = (
                "-" in aligned_ref or "-" in aligned_sample
            )

            # Detect sequence differences
            for ref_base, sample_base in zip(
                aligned_ref, aligned_sample
            ):

                if ref_base != "-":
                    ref_pos += 1

                if sample_base != "-":
                    sample_pos += 1

                if ref_base == sample_base:
                    continue

                mutation_type = ""
                position = 0

                mutation_info = {
                    "effect": "Not determined",
                    "ref_codon": "-",
                    "sample_codon": "-",
                    "ref_aa": "-",
                    "sample_aa": "-",
                    "protein_change": "Not determined"
                }

                if ref_base == "-":

                    mutation_type = "Insertion"
                    change = f"-→{sample_base}"
                    position = ref_pos + 1

                    mutation_info["effect"] = "Insertion"

                elif sample_base == "-":

                    mutation_type = "Deletion"
                    change = f"{ref_base}→-"
                    position = ref_pos

                    mutation_info["effect"] = "Deletion"

                else:

                    mutation_type = "Substitution"
                    change = f"{ref_base}→{sample_base}"
                    position = ref_pos

                    if has_indels:

                        mutation_info["effect"] = (
                            "Not assigned (alignment contains indels)"
                        )

                    else:

                        mutation_info = classify_effect(
                            reference_seq,
                            sample_seq,
                            ref_pos,
                            sample_pos
                        )

                codon_change = (
                    f'{mutation_info["ref_codon"]}'
                    f'→{mutation_info["sample_codon"]}'
                )

                mutations.append({

                    "change": change,
                    "position": position,
                    "type": mutation_type,
                    "label": "DETECTED",

                    "codon_change": codon_change,

                    "ref_aa": mutation_info["ref_aa"],
                    "sample_aa": mutation_info["sample_aa"],

                    "protein_change": mutation_info["protein_change"],

                    "effect": mutation_info["effect"],

                    "position_percent": (
                        position / max(len(reference_seq), 1)
                    ) * 100

                })

            # Calculate GC percentage
            gc_count = (
                sample_seq.count("G")
                + sample_seq.count("C")
            )

            gc_percent = round(
                (gc_count / len(sample_seq)) * 100,
                2
            )

            result = [{

                "id": sample_record.id,
                "description": sample_record.description,

                "length": len(sample_seq),
                "gc": gc_percent,

                "mutation_count": len(mutations),

                "sequence": sample_seq,
                "alignment": alignment_string,
                "alignment_score": alignment.score,

                "top_mutations": mutations[:25],

                "reference_id": reference_record.id,
                "reference_length": len(reference_seq)

            }]

        except ValueError as e:

            error = str(e)

        except Exception as e:

            error = f"Analysis error: {str(e)}"

    return render_template(
        "index.html",
        result=result,
        error=error
    )


if __name__ == "__main__":
    app.run(debug=True)
