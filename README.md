Winnowmap
========================================================================

Winnowmap is a long-read mapping algorithm, and a result of our exploration into superior minimizer sampling techniques. Minimizer sampling was originally introduced by [Roberts et al.](http://www.cs.toronto.edu/~wayne/research/papers/minimizers.pdf) This technique yields reduced representation of reference genome, enabling fast mapping analyses. A key property of the minimizer technique is that if two sequences share a substring of a specified length, then they can be guaranteed to have a matching minimizer. However, because the k-mer distribution in eukaryotic genomes is highly uneven, minimizer-based long read mappers (e.g., [minimap2](https://github.com/lh3/minimap2/), [mashmap](https://github.com/marbl/MashMap)) opt to discard the most frequently occurring minimizers from the genome in order to avoid excessive false positive seed hits. By doing so, the underlying guarantee is lost and accuracy is reduced in repetitive genomic regions (e.g., long tandem repeats).

To address the above problem, Winnowmap implements a novel **weighted minimizer** sampling algorithm. A unique feature of Winnowmap is that it performs minimizer sampling while taking into account a weight for each k-mer; i.e, the higher the weight of a k-mer, the more likely it is to be selected. Rather than masking k-mers, Winnowmap down-weights frequently occurring k-mers, thus reducing their chance of getting selected as minimizers. Winnowmap implements the new minimizer sampling and indexing algorithm, and borrows [minimap2’s](https://github.com/lh3/minimap2/) highly efficient anchor chaining and gapped alignment routines. The user-interface of Winnowmap is maintained similar to minimap2.

## Compile

Winnowmap requires c++11 to build, which is available in GCC >= 4.8. 
  ```sh
	git clone --recursive https://github.com/marbl/Winnowmap.git
	cd Winnomap
	make
  ```
Expect two executables `computeHighFreqKmers` and `winnowmap`.

## Usage

* Step 1: compute a set of highly frequent k-mers:
  ```sh
  computeHighFreqKmers 19 1 1024 ref.fa bad_Hk19_mers.txt (OR)
  computeHighFreqKmers 15 0 1024 ref.fa bad_k15_mers.txt
  ```
  The above executable `computeHighFreqKmers` expects five arguments in the following order: k-mer length, binary number 0/1 indicating whether to enable homopolymer compression, minimum k-mer frequency, reference genome and output file.  

* Step 2: Map long reads to reference:
  ```sh
  winnowmap -W bad_Hk19_mers.txt -cx map-pb ref.fa pacbio.fq.gz > output.paf (OR)
  winnowmap -W bad_k15_mers.txt -cx map-ont ref.fa ont.fq.gz > output.paf
  ```
  Except the `-W` parameter above needed by Winnowmap, the remaining options are consistent with [minimap2 usage](https://github.com/lh3/minimap2/blob/master/README.md).
  
  Users should keep k-mer length and homopolymer compression parameters consistent in Steps 1 and 2. For example, `map-pb` preset in Step 2 uses k-mer length 19 with homopolymer compression enabled. The other popular preset `map-ont` uses k-mer length 15 without homopolymer compression. In near future, these two steps will be merged into one  for user convenience.

## Benchmarking

Comparing Winnowmap to minimap2, we observed a reduction in the mapping error-rate from 0.14% to 0.06% in the recently finished [human X chromosome](https://github.com/nanopore-wgs-consortium/CHM13), and from 3.6% to 0% within the highly repetitive X centromere (3.1 Mbp). Winnowmap improves mapping accuracy within repeats and achieves these results with sparser sampling, leading to better index compression and competitive runtimes. By avoiding masking, we show that Winnowmap maintains uniform minimizer density.

<p align="center">
<img src="https://1aaaa1f6-a-62cb3a1a-s-sites.googlegroups.com/site/chirgjain/readme-winnowmap-density.jpg?attachauth=ANoY7cost_TsHo3yjf_COK13C-JBDQIio-GCb_hNSAdMQ92aRqISg21pJsg5dMKD5yMalcAugwI5vkqf9Cdu3sVk-xBz-SkRMkuyWAk3vK06_LEF2ay1pNSzCxU6nUNywhTYb5li8moC-YzRMmJZt7r3KFvcI34IbD7rktjXAPn_5Jba86E19uXq2o6zjAEDmsfjrKxqAdbsnPL3bU8L4wHwsH9gyv6170wD7WFJ_8pfFjeWam0v2uY%3D&attredirects=0" width=400px"> <br>
Minimizer sampling density using a human X chromosome as the reference, with the centromere positioned between 58 Mbp and 61 Mbp. ‘Standard’ method refers to the classic minimizer sampling algorithm from <a href="http://www.cs.toronto.edu/~wayne/research/papers/minimizers.pdf">Roberts et al.</a>, without any masking or modification.
</p>
