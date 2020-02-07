Winnowmap
========================================================================

Winnowmap is a new long-read mapping algorithm, and a result of our exploration into superior minimizer sampling techniques. Minimizer sampling was originally introduced by [Roberts et al.](http://www.cs.toronto.edu/~wayne/research/papers/minimizers.pdf). This technique yields reduced representation of reference genome, enabling fast mapping. A key property of the minimizer technique is that if two sequences share a substring of a specified length, then they can be guaranteed to have a matching minimizer. However, because the k-mer distribution in eukaryotic genomes is highly uneven, minimizer-based long read mappers (e.g., [minimap2](https://github.com/lh3/minimap2/), [mashmap](https://github.com/marbl/MashMap)) opt to discard the most frequently occurring minimizers from the genome in order to avoid excessive false positive seed hits. By doing so, the underlying guarantee is lost and accuracy is reduced in repetitive genomic regions (e.g., long tandem repeats).

To address the above problem, Winnowmap implements a novel **weighted minimizer** sampling algorithm. A unique feature of the proposed algorithm is that it performs minimizer sampling while taking into account a weight for each k-mer; i.e, the higher the weight of a k-mer, the more likely it is to be selected. Rather than masking k-mers, Winnowmap opts to down-weight frequently occurring k-mers, thus reducing their chance of getting selected as minimizers. Winnowmap implements the new minimizer sampling and indexing algorithm, and borrows [minimap2’s](https://github.com/lh3/minimap2/) highly efficient anchor chaining and gapped alignment routines. The user-interface of Winnowmap is maintained similar to minimap2.

Comparing Winnowmap to minimap2, we observe a reduction in the mapping error-rate from 0.14% to 0.06% in the recently finished [human X chromosome](https://github.com/nanopore-wgs-consortium/CHM13), and from 3.6% to 0% within the highly repetitive X centromere (3.1 Mbp). Winnowmap improves mapping accuracy within repeats and achieves these results with sparser sampling, leading to better index compression and competitive runtimes.

## Compile

Winnowmap requires c++11 to build, which is available in GCC >= 4.8. To compile Winnowmap, run the `make` command. Expect two executables `computeHighFreqKmers` and `winnowmap`.

## Usage

* Step 1: compute a set of highly frequent k-mers:
  ```sh
  computeHighFreqKmers 19 1 1024 ref.fa bad_kmers.txt
  ```
  The above executable `computeHighFreqKmers` expects five arguments in the following order: k-mer length, binary integer 0/1 indicating whether to enable homopolymer compression, minimum k-mer frequency, reference genome and output file.  

* Step 2: Map long reads to reference:
  ```sh
  winnowmap -W bad_kmers.txt -cx map-pb ref.fa pacbio.fq.gz > output.paf
  ```
  Except the `-W` parameter above needed by Winnowmap, the remaining options are consistent with [minimap2 usage](https://github.com/lh3/minimap2/blob/master/README.md).
  
  Users should keep k-mer length and homopolymer compression parameters consistent in Steps 1 and 2. For example, `map-pb` preset in Step 2 uses k-mer length 19 with homopolymer compression enabled. The other popular preset `map-ont` uses k-mer length 15 without homopolymer compression. Later, the two steps will be merged into one in a subsequent release of Winnowmap for user convenience.

