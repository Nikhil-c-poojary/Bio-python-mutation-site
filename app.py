from flask import Flask, render_template, request
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from io import StringIO

app = Flask(__name__)

VALID_BASES = set("ATGCN")


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
            f"Invalid DNA characters found: {', '.join(sorted(invalid))}"
        )

    return record


def classify_effect(reference_seq, sample_seq, ref_pos, sample_pos):

    ref_start = ((ref_pos - 1) // 3) * 3
    sample_start = ((sample_pos - 1) // 3) * 3

    ref_codon = reference_seq[ref_start:ref_start + 3]
    sample_codon = sample_seq[sample_start:sample_start + 3]

    if len(ref_codon) != 3 or len(sample_codon) != 3:
        return "Unknown (incomplete codon)", "-", "-"

    if "N" in ref_codon or "N" in sample_codon:
        return "Unknown (ambiguous base)", ref_codon, sample_codon

    ref_aa = str(Seq(ref_codon).translate())
    sample_aa = str(Seq(sample_codon).translate())

    if ref_aa == sample_aa:
        effect = "Synonymous"

    elif sample_aa == "*":
        effect = "Stop gained"

    elif ref_aa == "*":
        effect = "Stop lost"

    else:
        effect = "Nonsynonymous"

    return effect, ref_codon, sample_codon


@app.route("/", methods=["GET", "POST"])
def index():

    result = None
    error = None

    if request.method == "POST":

        try:

            sample_text = request.form.get("sample_sequence", "")
            reference_text = request.form.get("reference_sequence", "")

            sample_record = read_fasta(sample_text)
            reference_record = read_fasta(reference_text)

            sample_seq = str(sample_record.seq).upper()
            reference_seq = str(reference_record.seq).upper()

            # Perform sequence alignment

            alignments = pairwise2.align.globalms(
                reference_seq,
                sample_seq,
                2,
                -1,
                -5,
                -0.5
            )

            if not alignments:
                raise ValueError("Unable to align the sequences.")

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
                aligned_ref
                + "\n"
                + "".join(alignment_visual)
                + "\n"
                + aligned_sample
            )

            mutations = []

            ref_pos = 0
            sample_pos = 0

            for ref_base, sample_base in zip(
                aligned_ref, aligned_sample
            ):

                if ref_base != "-":
                    ref_pos += 1

                if sample_base != "-":
                    sample_pos += 1

                if ref_base == sample_base:
                    continue

                if ref_base == "-":

                    mutation_type = "Insertion"
                    change = f"-→{sample_base}"
                    effect = "Insertion"
                    codon_change = "-"

                    position = ref_pos + 1

                elif sample_base == "-":

                    mutation_type = "Deletion"
                    change = f"{ref_base}→-"
                    effect = "Deletion"
                    codon_change = "-"

                    position = ref_pos

                else:

                    mutation_type = "Substitution"

                    change = f"{ref_base}→{sample_base}"

                    position = ref_pos

                    effect, ref_codon, sample_codon = (
                        classify_effect(
                            reference_seq,
                            sample_seq,
                            ref_pos,
                            sample_pos
                        )
                    )

                    codon_change = (
                        f"{ref_codon}→{sample_codon}"
                    )

                mutations.append({

                    "change": change,

                    "position": position,

                    "label": "DETECTED",

                    "probability": "Not estimated",

                    "effect": effect,

                    "type": mutation_type,

                    "codon_change": codon_change,

                    "context": "Observed sequence difference",

                    "position_percent": (
                        position / max(len(reference_seq), 1)
                    ) * 100

                })

            # GC percentage

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
