from flask import Flask, render_template, request
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
                    probability = 95
                    label = "HIGH"

                    mutation = {
                        "change": f"{ref_base}→{sample_base}",
                        "position": pos,
                        "probability": probability,
                        "label": label,
                        "codon": "---",
                        "mutated_codon": "---",
                        "amino_acid": "Unknown",
                        "mutated_amino_acid": "Unknown",
                        "effect": "Substitution",
                        "context": "Detected through sequence alignment",
                        "position_percent": (
                            pos / max(len(reference_seq), 1)
                        ) * 100,
                    }

                    mutations.append(mutation)

            gc_count = sample_seq.count("G") + sample_seq.count("C")
            gc_percent = round((gc_count / len(sample_seq)) * 100, 2)

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
