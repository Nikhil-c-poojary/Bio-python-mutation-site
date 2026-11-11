<!DOCTYPE html>
<html lang="en">

<head>

    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <title>Gene Mutation Possibility Analyzer</title>

    <link rel="stylesheet"
          href="{{ url_for('static', filename='styles.css') }}">

</head>

<body>

<main class="shell">

    <!-- HERO SECTION -->

    <section class="hero">

        <div>

            <p class="eyebrow">BIOPYTHON • GENOMIC ANALYSIS</p>

            <h1>Gene Mutation Possibility Analyzer</h1>

            <p class="lede">
                Analyze DNA sequences, identify mutations,
                compare codon changes and predict possible
                amino acid effects using Biopython.
            </p>

        </div>

        <div class="metrics">

            <span>POWERED BY</span>

            <strong>Flask + Biopython</strong>

            <small>Sequence Analysis Platform</small>

        </div>

    </section>


    <!-- INPUT WORKSPACE -->

    <section class="workspace">

        <form method="POST" class="input-panel" id="sequenceForm">

            <h2>Sequence Input</h2>

            <p class="small">
                Enter reference and sample DNA sequences in FASTA format.
            </p>

            <label for="sample_sequence">
                Sample FASTA Sequence
            </label>

            <textarea
                id="sample_sequence"
                name="sample_sequence"
                placeholder=">Sample_Sequence&#10;ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
                required>{{ request.form.get('sample_sequence', '') }}</textarea>


            <label for="reference_sequence">
                Reference FASTA Sequence
            </label>

            <textarea
                id="reference_sequence"
                name="reference_sequence"
                placeholder=">Reference_Sequence&#10;ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
                required>{{ request.form.get('reference_sequence', '') }}</textarea>


            <div class="controls">

                <button type="submit">
                    Analyze DNA Sequence →
                </button>

                <button
                    type="button"
                    class="secondary-button"
                    onclick="clearSequences()">

                    Clear

                </button>

            </div>

        </form>


        <!-- FEATURES PANEL -->

        <aside class="note-panel">

            <p class="eyebrow">ANALYSIS MODULES</p>

            <h2>Platform Features</h2>

            <div class="feature-item">
                <span>01</span>
                <p>FASTA Sequence Parsing</p>
            </div>

            <div class="feature-item">
                <span>02</span>
                <p>Reference Sequence Alignment</p>
            </div>

            <div class="feature-item">
                <span>03</span>
                <p>Point Mutation Detection</p>
            </div>

            <div class="feature-item">
                <span>04</span>
                <p>Codon and Amino Acid Analysis</p>
            </div>

            <div class="feature-item">
                <span>05</span>
                <p>GC Content Calculation</p>
            </div>

            <div class="feature-item">
                <span>06</span>
                <p>Mutation Distribution Mapping</p>
            </div>

            <p class="small">
                Designed for educational bioinformatics
                and DNA sequence analysis.
            </p>

        </aside>

    </section>


    <!-- ERROR MESSAGE -->

    {% if error %}

    <section class="alert">

        <h3>Analysis Error</h3>

        <p>{{ error }}</p>

    </section>

    {% endif %}


    <!-- RESULTS -->

    {% if result %}

    {% for record in result %}

    <section class="results">

        <!-- RESULT HEADER -->

        <div class="result-head">

            <div>

                <p class="eyebrow">SEQUENCE ANALYSIS REPORT</p>

                <h2>{{ record.id }}</h2>

                <p>{{ record.description }}</p>

                <p class="small">
                    Reference: {{ record.reference_id }}
                    | Reference Length: {{ record.reference_length }} bp
                </p>

            </div>


            <!-- STATISTICS -->

            <div class="stat-grid">

                <div>
                    <span>Sequence Length</span>
                    <strong>{{ record.length }} bp</strong>
                </div>

                <div>
                    <span>GC Content</span>
                    <strong>{{ record.gc }}%</strong>
                </div>

                <div>
                    <span>Detected Differences</span>
                    <strong>{{ record.mutation_count }}</strong>
                </div>

            </div>

        </div>


        <!-- DNA SEQUENCE -->

        <section class="sequence-panel">

            <h2>Sample DNA Sequence</h2>

            <div class="sequence-strip">
                {{ record.sequence }}
            </div>

        </section>


        <!-- ALIGNMENT -->

        <section class="alignment-panel">

            <h2>Sequence Alignment</h2>

            <p class="small">
                Identical bases are represented by vertical bars.
            </p>

            <div class="alignment-box">

                <pre>{{ record.alignment }}</pre>

            </div>

            <p class="small">
                Alignment Score: {{ record.alignment_score }}
            </p>

        </section>


        <!-- MUTATION DISTRIBUTION -->

        <section class="graph-panel">

            <h2>Mutation Distribution</h2>

            <p class="small">
                Distribution of detected sequence differences
                across the reference sequence.
            </p>

            <div class="graph-bar">

                {% for mutation in record.top_mutations %}

                <div
                    class="mutation-marker"
                    style="left: {{ mutation.position_percent }}%;"
                    title="{{ mutation.change }} at position {{ mutation.position }}">
                </div>

                {% endfor %}

            </div>

            <div class="graph-labels">

                <span>Start</span>

                <span>Reference Sequence</span>

                <span>End</span>

            </div>

        </section>


        <!-- MUTATION TABLE -->

        <section class="table-section">

            <h2>Detected Mutations</h2>

            <p class="small">
                Sequence differences identified through
                reference and sample alignment.
            </p>

            {% if record.top_mutations %}

            <div class="table-wrap">

                <table>

                    <thead>

                        <tr>

                            <th>No.</th>
                            <th>Mutation</th>
                            <th>Position</th>
                            <th>Type</th>
                            <th>Codon Change</th>
                            <th>Reference Amino Acid</th>
                            <th>Sample Amino Acid</th>
                            <th>Protein Change</th>
                            <th>Effect</th>
                            <th>Status</th>

                        </tr>

                    </thead>

                    <tbody>

                        {% for mutation in record.top_mutations %}

                        <tr>

                            <td>{{ loop.index }}</td>

                            <td>
                                <strong>{{ mutation.change }}</strong>
                            </td>

                            <td>{{ mutation.position }}</td>

                            <td>{{ mutation.type }}</td>

                            <td>{{ mutation.codon_change }}</td>

                            <td>{{ mutation.ref_aa }}</td>

                            <td>{{ mutation.sample_aa }}</td>

                            <td>
                                <strong>
                                    {{ mutation.protein_change }}
                                </strong>
                            </td>

                            <td>{{ mutation.effect }}</td>

                            <td>
                                <span class="badge detected">
                                    {{ mutation.label }}
                                </span>
                            </td>

                        </tr>

                        {% endfor %}

                    </tbody>

                </table>

            </div>

            {% if record.mutation_count > 25 %}

            <p class="small">
                Showing the first 25 detected differences
                out of {{ record.mutation_count }}.
            </p>

            {% endif %}

            {% else %}

            <div class="success-message">

                <h3>No Sequence Differences Detected</h3>

                <p>
                    The aligned reference and sample sequences
                    show no differences.
                </p>

            </div>

            {% endif %}

        </section>


        <!-- AMINO ACID ANALYSIS -->

        <section class="amino-acid-panel">

            <h2>Amino Acid Analysis</h2>

            <p class="small">
                Predicted protein-level effects of detected DNA mutations.
            </p>

            {% for mutation in record.top_mutations %}

            <div class="feature-item">

                <span>{{ loop.index }}</span>

                <p>

                    <strong>{{ mutation.codon_change }}</strong>

                    |

                    {{ mutation.ref_aa }}

                    →

                    {{ mutation.sample_aa }}

                    |

                    {{ mutation.protein_change }}

                </p>

            </div>

            {% endfor %}

        </section>


        <!-- SCIENTIFIC NOTE -->

        <section class="scientific-note">

            <h3>Scientific Interpretation</h3>

            <p>
                Detected sequence differences represent variations
                between the supplied reference and sample DNA sequences.
                Synonymous and nonsynonymous classifications are based
                on codon translation.
            </p>

            <p>
                Amino acid predictions use the standard genetic code
                and assume that the supplied sequence begins at the
                first nucleotide of a coding region.
            </p>

            <p>
                Mutation detection does not establish pathogenicity
                or clinical significance. Results depend on sequence
                orientation, coding frame and alignment accuracy.
            </p>

        </section>

    </section>

    {% endfor %}

    {% endif %}


    <!-- FOOTER -->

    <footer class="footer">

        <p>
            Gene Mutation Possibility Analyzer
        </p>

        <span>
            Developed using Python, Flask and Biopython
        </span>

        <p class="small">
            Educational Bioinformatics Project by
            NIKHIL.C [Final year B.Sc Alva's]
        </p>

    </footer>

</main>


<!-- CLEAR BUTTON SCRIPT -->

<script>

function clearSequences() {

    document.getElementById("sample_sequence").value = "";

    document.getElementById("reference_sequence").value = "";

    // Remove previous results and errors
    const results = document.querySelector(".results");
    const alert = document.querySelector(".alert");

    if (results) {
        results.remove();
    }

    if (alert) {
        alert.remove();
    }

}

</script>

</body>
</html>
