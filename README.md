# LOSAT🚀
LOSAT (LOcal Sequence Alignment Tool) aims to provide a miniaturized reimplementation of NCBI BLAST algorithm with webassembly compatibility.

## Motivation
[gbdraw](https://github.com/satoshikawato/gbdraw/) visualizes pairwise alignments of neighbboring genomes, but the problem was that the user has to provide BLAST alignment results generated outside of gbdraw. This is particularly problematic for the [webapp version](https://gbdraw.app/).
As there seems to be no WebAssembly-comaptible builds for BLAST, I reimplemented BLAST algorithms in Rust.


## Currently supported modes
- BLASTN
- TBLASTX

## Performance
### Hit distribution
![overall_trend_comparison.png](./LOSAT/tests/plots/overall_trend_comparison.png)
### Execution time
![execution_time_comparison_all.png](./LOSAT/tests/plots/execution_time_comparison_all.png)

## References
[NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

[Wikipedia](https://en.wikipedia.org/wiki/MGM-166_LOSAT)
