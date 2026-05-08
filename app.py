from flask import Flask, render_template, request
from Bio import SeqIO, pairwise2
from io import StringIO

app = Flask(__name__)


def read_fasta(content):
    handle = StringIO(content)

    records = list(SeqIO.parse(handle, "fasta"))

    if not records:
        raise ValueError("No valid FASTA sequence found")

    return records[0]


@app.route("/", methods=["GET", "POST"])
def index():

    result = None
    error = None
    fasta_text = ""

    if request.method == "POST":

        try:

            sample_text = request.form.get("sample_sequence")
            reference_text = request.form.get("reference_sequence")

            if not sample_text or not reference_text:
                raise ValueError(
                    "Both sample and reference sequences are required"
                )

            fasta_text = sample_text

            sample_record = read_fasta(sample_text)
            reference_record = read_fasta(reference_text)

            sample_seq = str(sample_record.seq).upper()
            reference_seq = str(reference_record.seq).upper()

            alignment = pairwise2.align.globalxx(
                reference_seq,
                sample_seq
            )[0]

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

            pos = 0

            for i in range(len(aligned_ref)):

                ref_base = aligned_ref[i]
                sample_base = aligned_sample[i]

                if ref_base != "-":
                    pos += 1

                if ref_base != sample_base:

                    mutation = {
                        "change": f"{ref_base}→{sample_base}",
                        "position": pos,
                        "probability": 95,
                        "label": "HIGH",
                        "effect": "Substitution",
                        "context": "Detected through sequence alignment",
                        "position_percent": (
                            pos / max(len(reference_seq), 1)
                        ) * 100,
                    }

                    mutations.append(mutation)

            gc_count = (
                sample_seq.count("G")
                + sample_seq.count("C")
            )

            gc_percent = round(
                (gc_count / len(sample_seq)) * 100,
                2
            )

            result = [
                {
                    "id": sample_record.id,
                    "description": sample_record.description,
                    "length": len(sample_seq),
                    "gc": gc_percent,
                    "mutation_count": len(mutations),
                    "sequence": sample_seq,
                    "alignment": alignment_string,
                    "top_mutations": mutations[:25],
                }
            ]

        except Exception as e:
            error = str(e)

    return render_template(
        "index.html",
        result=result,
        error=error,
        fasta_text=fasta_text,
    )


if __name__ == "__main__":
    app.run(debug=True)
