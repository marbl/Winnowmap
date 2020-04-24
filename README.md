Winnowmap
========================================================================

Winnowmap is a long-read mapping algorithm, and a result of our exploration into superior minimizer sampling techniques. Minimizer sampling was originally introduced by [Roberts et al.](http://www.cs.toronto.edu/~wayne/research/papers/minimizers.pdf) This technique yields reduced representation of reference genome, enabling fast mapping analyses. A key property of the minimizer technique is that if two sequences share a substring of a specified length, then they can be guaranteed to have a matching minimizer. However, because the k-mer distribution in eukaryotic genomes is highly uneven, minimizer-based long read mappers (e.g., [minimap2](https://github.com/lh3/minimap2/), [mashmap](https://github.com/marbl/MashMap)) opt to discard the most frequently occurring minimizers from the genome in order to avoid excessive false positive seed hits. By doing so, the underlying guarantee is lost and accuracy is reduced in repetitive genomic regions (e.g., long tandem repeats).

To address the above problem, Winnowmap implements a novel **weighted minimizer** sampling algorithm. A unique feature of Winnowmap is that it performs minimizer sampling while taking into account a weight for each k-mer; i.e, the higher the weight of a k-mer, the more likely it is to be selected. Rather than masking k-mers, Winnowmap down-weights frequently occurring k-mers, thus reducing their chance of getting selected as minimizers. Users can refer to our [preprint](https://doi.org/10.1101/2020.02.11.943241) for algorithmic details. Winnowmap implements the new minimizer sampling and indexing algorithm, and borrows [minimap2’s](https://github.com/lh3/minimap2/) highly efficient anchor chaining and gapped alignment routines. The user-interface of Winnowmap is maintained similar to minimap2.

## Compile

Winnowmap requires c++11 to build, which is available in GCC >= 4.8. 
  ```sh
	git clone --recursive https://github.com/marbl/Winnowmap.git
	cd Winnowmap
	make -j8
  ```
Expect `winnowmap` and `meryl` executables in `bin` folder. The `recursive` option used above is necessary to download all submodules .

## Usage

For either mapping long reads or computing whole-genome alignments, Winnowmap requires pre-computing high frequency k-mers (e.g., top 0.02% most frequent) in a reference. Winnowmap uses [meryl](https://github.com/marbl/meryl) k-mer counting tool for this purpose.  

*  Mapping ONT reads
  ```sh
	meryl count k=15 output merylDB_k15 ref.fa
	meryl print greater-than distinct=0.9998 merylDB_k15 > bad_k15_mers.txt

	winnowmap -W bad_k15_mers.txt -cx map-ont ref.fa ont.fq.gz > output.paf
  ```
  Except the `-W` parameter above needed by Winnowmap, the remaining options are consistent with [minimap2 usage](https://github.com/lh3/minimap2/blob/a79cb3e9912a89d823a90b3502f57798afe299f8/README.md).

*  Mapping PacBio (CLR / HiFi) reads

  ```sh
	meryl count compress k=19 output merylDB_Hk19 ref.fa
	meryl print greater-than distinct=0.9998 merylDB_Hk19 > bad_Hk19_mers.txt

	winnowmap -W bad_Hk19_mers.txt -cx map-pb ref.fa pacbio.fq.gz > output.paf
  ```

*  Mapping genome assemblies

  ```sh
	meryl count k=19 output merylDB_k19 asm1.fa
	meryl print greater-than distinct=0.9998 merylDB_k19 > bad_k19_mers.txt

	winnowmap -W bad_k19_mers.txt -cx asm5 asm1.fa asm2.fa > output.paf
  ```

  The k-mer counting step using meryl must be consistent with the final mapping step. For example, ONT reads are mapped with k-mer length 15, where as PacBio reads are mapped with homopolymer-compressed k-mers of length 19. This is kept consistent with `map-ont` and `map-pb` presets of minimap2. In all use-cases, pre-computing repetitive k-mers using [meryl](https://github.com/marbl/meryl) is quite fast, it typically takes 2-3 minutes for the human genome reference. 

## Benchmarking

Comparing Winnowmap to minimap2, we observed a reduction in the mapping error-rate from 0.14% to 0.06% in the recently finished [human X chromosome](https://github.com/nanopore-wgs-consortium/CHM13), and from 3.6% to 0% within the highly repetitive X centromere (3.1 Mbp). Winnowmap improves mapping accuracy within repeats and achieves these results with sparser sampling, leading to better index compression and competitive runtimes. By avoiding masking, we show that Winnowmap maintains uniform minimizer density.

<p align="center">
<img src="https://1aaaa1f6-a-62cb3a1a-s-sites.googlegroups.com/site/chirgjain/readme-winnowmap-density.jpg?attachauth=ANoY7cost_TsHo3yjf_COK13C-JBDQIio-GCb_hNSAdMQ92aRqISg21pJsg5dMKD5yMalcAugwI5vkqf9Cdu3sVk-xBz-SkRMkuyWAk3vK06_LEF2ay1pNSzCxU6nUNywhTYb5li8moC-YzRMmJZt7r3KFvcI34IbD7rktjXAPn_5Jba86E19uXq2o6zjAEDmsfjrKxqAdbsnPL3bU8L4wHwsH9gyv6170wD7WFJ_8pfFjeWam0v2uY%3D&attredirects=0" width=400px"> <br>
Minimizer sampling density using a human X chromosome as the reference, with the centromere positioned between 58 Mbp and 61 Mbp. ‘Standard’ method refers to the classic minimizer sampling algorithm from <a href="http://www.cs.toronto.edu/~wayne/research/papers/minimizers.pdf">Roberts et al.</a>, without any masking or modification.
</p>

## Publication

- **Chirag Jain, Arang Rhie, Haowen Zhang, Chaudia Chu, Brian Walenz, Sergey Koren and Adam Phillippy**. "Weighted minimizer sampling improves long read mapping". *Bioinformatics (ISMB proceedings)*, 2020.
