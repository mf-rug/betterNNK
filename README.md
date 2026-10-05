# betterNNK

A standalone, browser-based multi-site saturation mutagenesis primer designer. No cloud services, accounts, dependencies, sequence uploads or build step.

## Run locally

From this directory:

```sh
npm start
```

Open http://127.0.0.1:49173. Requires Python 3; npm is only a convenience wrapper. Alternatively:

```sh
python3 -m http.server 49173 --bind 127.0.0.1 --directory web
```

Stop with Ctrl+C. If the port is occupied, use another port in the Python command. Serve the `web` directory rather than opening index.html directly (JavaScript modules require HTTP).

## Workflow

1. Paste an unmodified A/C/G/T primer, written 5′ to 3′. Whitespace and lowercase are normalized. Invalid characters are rejected, never silently removed.
2. Confirm one of six reading frames. The first stop-free frame is suggested; multiple candidates are explicitly flagged. Complete codons only; the standard genetic code is used.
3. Hover and click triplets to toggle multiple selections; Clear all removes all targets. Changing the frame clears selections. Keyboard Tab and Enter/Space work too.
4. Choose 22c, 20c (Tang), NNK, or Other. Other offers NNS, NDT, VMA, VRK, NNN, NNB, VHG, NYC, KST, RRK and a custom three-symbol DNA IUPAC codon. The bar chart includes all 20 amino acids and stops, with absolute theoretical percentages.
5. Copy/download full design/mixing TSV including pool IDs (default), or two-column ordering TSV. Header rows are optional.

A visible Primer name field supplies the common name for all oligos, without an `_input` suffix. Order settings provide optional full reverse-complement oligos. Each strand is a separate mixing pool. This is not a general PCR primer-pair optimizer. Input positions are 1-based and refer to the pasted primer; amino acid indices are local to the selected frame, not absolute protein numbering.

## Scientific rules

- 22c: NDT + VHG + TGG, molar parts 12:9:1. 22 distinct codons, 20 amino acids; Leu and Val each have two codons.
- Tang: NDT + VMA + ATG + TGG, molar parts 12:6:1:1. 20 distinct codons, one per canonical amino acid.
- Reverse-frame targets use AHN/CDB/CCA or AHN/TKB/CAT/CCA in the input-oriented sequence.
- Every non-target base is preserved. Incomplete terminal codons are not selectable.
- The app does not infer a biological reading frame with certainty or validate annealing, specificity, hairpins, flanking length, or experimental library balance.

See [research/METHODS.md](research/METHODS.md) for paper review and existing implementation findings.

## Validation

```sh
npm test
npm run check
```

Tests use Node's built-in runner (Node 18+), with no installed dependencies. They verify genetic-code coverage, exact degeneracy and mixing weights, both published primer examples, IUPAC reverse complements, all six frames across different sequence lengths, unchanged flanks, input validation, and TSV structure.

## Self-hosting later

Upload only `web/` to any static HTTP server (nginx, Caddy, Apache, or a static hosting service). The app needs no backend, database, API keys, package installation or OpenAI services. Serve `.js` as JavaScript. Use HTTPS for clipboard access outside localhost; download and manual-copy fallback remain available. All asset paths are relative, allowing deployment under a subdirectory. Python's development server is for local use; use a production web server for public hosting. Streamlit could embed this UI, but is not required.

## Structure

- `web/primer.js`: framework-independent scientific logic and TSV serialization.
- `web/app.js`: browser state, selection and export interactions.
- `web/index.html`, `web/style.css`: responsive interface.
- `tests/primer.test.js`: executable scientific and edge-case checks.
- `literature/`: original user-provided PDFs, kept outside the served directory.

Future useful controls: absolute protein-residue numbering, a concentration-aware mixing calculator.

## Multi-site libraries

The shared method is the default. Each selected position can override it with 22c, Tang, NNK, a preset or a custom IUPAC codon. Overrides survive shared-default changes and can be reset together; deselecting a site or changing the frame clears its override. Separate mode generates an independent library for each target (all other bases unchanged). Combinatorial mode generates the Cartesian product of method components, replacing every selected site in every oligo. Each site contributes its own component choices. Component molar weights multiply, and fractions normalize within each strand pool. Two 22c sites yield 6 separate-library oligos or 9 combinatorial oligos; complementary oligos double these counts.

The preview reports ordered oligo count, library count, possible sense amino-acid combinations per library, and stop probability at randomized sites. The distribution chart has a position selector to inspect each site or the shared default. Combinatorial stop probability is 1 minus the product of each site’s non-stop probability, assuming independent incorporation. It excludes pre-existing stops at unselected positions. Possible diversity does not guarantee experimental coverage.

Exports identify pools and list target coordinates and codon components in ascending input-position order even on reverse frames. Keep separate site pools separate. Order-only TSV omits pool/mixing columns; use the default full TSV as the mixing record. Generation is capped at 2,048 total oligos to bound browser work; the count is checked before creating combinations.
