# NblockMatcher
Find sequence motif combinations according to custom specification.

The motif search itself is currently made by `FIMO` from `MEME suite`.

### Basic description:
Given motif combination, the script will find all motifs individually with `FIMO` and then filter them, so their reported
 combinations follows the given criteria (order, strand, distance).
The score of the combined motif satisfying given criteria is taken as a sum of the scores. The p-value for the combined 
 motif can be optionally computed with `pytfmpval` package (separate installation is needed at this moment).

### Refs:
- Fimo: Charles E. Grant, Timothy L. Bailey and William Stafford Noble, "FIMO: Scanning for occurrences of a given motif", Bioinformatics 27(7):1017-1018, 2011.
- pytfmpval: Touzet, H., Varré, JS. Efficient and accurate P-value computation for Position Weight Matrices. Algorithms Mol Biol 2, 15 (2007). https://doi.org/10.1186/1748-7188-2-15