# Method review and implementation decisions

## Supplied papers

1. Tang et al. (2012), *Construction of “small-intelligent” focused mutagenesis libraries using well-designed combinatorial degenerate primers*, DOI https://doi.org/10.2144/000113820. Methods and Results (p.150) give NDT/VMA/ATG/TGG and 12:6:1:1. The supplied forward example is TCTTCAACXXXCCGACATA. The reverse example uses AHN/TKB/CAT/CCA. These facts are direct fixtures in the tests. DC-Analyzer's arbitrary-amino-acid optimization is unnecessary for this fixed complete alphabet.
2. Kille et al. (2013), *Reducing Codon Redundancy and Screening Effort of Combinatorial Protein Libraries Created by Saturation Mutagenesis*, DOI https://doi.org/10.1021/sb300037w. Results pp.84–85 specify NDT/VHG/TGG, yielding 22 codons, with only Leu and Val represented twice. The PCR strategy determines whether one strand or both strands are ordered.
3. Kille supplementary information `sb300037w_si_001.pdf`: Table S1 explicitly gives 12:9:1 for one site, product weights for multiple sites; Table S4 supplies forward CHMO-426 CATGTTCATGGTGxxxGGACCGAATGGC and reverse GCCATTCGGTCCyyyCACCATGAACATG with xxx=NDT/VHG/TGG and yyy=AHN/CDB/CCA. The tests reproduce this pair. Some multi-site strings and numbering in the extracted SI appear inconsistent; the app derives complements algebraically rather than copying multi-site examples. Table S3 also appears to list ACT under NDT/Ser in extracted text; correct NDT expansion gives AGT. No table transcription is used as the genetic code.
4. MDC-Analyzer (2014), DOI https://doi.org/10.2144/000114177: an optimization strategy for desired combinations of amino acids at contiguous sites, allowing extra sequences to trade library size against primer number. Uses greedy best-first search and pruning. This is a different problem from the initial single-site full-alphabet recipes and is not implemented. Its table confirms primer counts grow as 4^n for Tang and 3^n for 22c, before adding complementary strands.
5. Qu et al. (2019), *The Crucial Role of Methodology Development in Directed Evolution of Selective Enzymes*: section 2.3.3 discusses 22c and Tang, confirms 22:20 vs 20:20 and the increased primer count for Tang. The review provides context for future reduced alphabets and alternatives, not a different single-site recipe.

## Existing implementations and online cross-checks

Searches covered GitHub and publications for 22c, Tang, MDC-Analyzer, SwiftLib, CodonGenie and Mutation Maker.

- SwiftLib authors' JavaScript repository: https://github.com/aleaverfay/swiftlib_javascript ; accompanying paper https://pmc.ncbi.nlm.nih.gov/articles/PMC4357694/. A dynamic-programming optimizer for degenerate libraries with broader constraints. It is a useful future reference but does not need to become a dependency for two fixed single-site recipes.
- Mutation Maker: https://github.com/ra100/Mutation_Maker (mirror of the Merck project), associated publication https://doi.org/10.1021/acssynbio.0c00542. A much broader mutagenic oligo design platform with backend services; not adopted here.
- DegenPrimer: https://github.com/allista/DegenPrimer . A degenerate-primer quality/PCR simulation tool rather than the requested fixed recipe generator; not used as an oracle.
- Independent experimental/economic comparison: https://www.nature.com/articles/srep10654 confirms both mixtures and theoretical amino-acid probabilities. Its methods paragraph appears to say VHA in one location while its earlier recipe specifies VMA; VHA would not give Tang's six-codon component. The original Tang paper and exhaustive codon enumeration resolve this in favor of VMA.

No third-party algorithm code is copied or assumed correct. Verification uses the original papers' explicit single-site fixtures plus exhaustive expansion of the implemented IUPAC codons with the standard genetic code. The external projects have not undergone a full correctness audit and are not claimed as validation oracles.

## Coordinate contract

Input is always the ordered oligo direction, 5′→3′. Forward frame offsets 0/1/2 count from its left edge. Reverse frame offsets 0/1/2 count from the 5′ edge of the reverse complement. A reverse-complement codon starting at zero-based position p maps to input interval [length-p-3, length-p). The interface keeps the input strand on screen, translates from the chosen coding strand, and labels reverse coding direction. A reverse selection replaces the input triplet with the reverse complement of the method's coding triplet. A complementary oligo, if requested, is the reverse complement of the entire generated oligo.

## Deliberate bounds

Shared-method multi-site selection in separate or combinatorial mode; standard genetic code; no existing ambiguity symbols in input; no silent removal of invalid symbols. No host-specific codon optimization, PCR temperature calculation, primer specificity prediction or wild-type exclusion. Mixing percentages are per strand pool and assume equal incorporation within degenerate positions. Biological/synthesis biases can still occur. Frame auto-selection is a deterministic heuristic, not proof.


## Multi-site extension

Implemented Cartesian-product component generation using Kille SI Table S1 product weights. Two-site 22c weights are 144, 108, 12, 108, 81, 9, 12, 9, 1 (sum 484). Separate-library generation retains all unselected positions, including other selected sites belonging to other libraries. Explicit pool IDs distinguish libraries and complementary strands. Tests verify these weights, within-pool normalization, reverse-strand placement, counts, combined stop probability and the 2,048-oligo generation limit.
