Winnowmap
========================================================================

Winnowmap is a long-read mapping algorithm optimized for mapping ONT and PacBio reads to repetitive reference sequences. Winnowmap development began on top of [minimap2](https://github.com/lh3/minimap2/) codebase, and since then we have incorporated the following two ideas to improve mapping accuracy within repeats. 

- Winnowmap implements a novel **weighted minimizer** sampling algorithm (>=v1.0). This optimization was motivated by the need to avoid masking of frequently occurring k-mers during the seeding stage in an efficient manner, and achieve better mapping accuracy in complex repeats (e.g., long tandem repeats) of the human genome. Using weighted minimizers, Winnowmap down-weights frequently occurring k-mers, thus reducing their chance of getting selected as minimizers. Users can refer to [this paper](https://doi.org/10.1093/bioinformatics/btaa435) for more details. This idea is helpful to preserve the [theoretical guarantee](http://www.cs.toronto.edu/~wayne/research/papers/minimizers.pdf) of minimizer sampling technique, i.e., if two sequences share a substring of a specified length, then they must be guaranteed to have a matching minimizer.   

- We noticed that the highest scoring alignment doesn't necessarily correspond to correct placement of reads in repetitive regions of [T2T human chromosomes](https://github.com/nanopore-wgs-consortium/CHM13). In the presence of a non-reference allele within a repeat, a read sampled from that region could be mapped to an incorrect repeat copy because the standard pairwise sequence alignment scoring system penalizes true variants. This is also sometimes referred to as allelic bias. To address this bias, we introduced and implemented an idea of using **minimal confidently alignable substrings** (>=v2.0). These are minimal-length substrings in a read that align end-to-end to a reference with [mapping quality score](https://genome.sph.umich.edu/wiki/Mapping_Quality_Scores) above a user-specified threshold. This approach treats each read mapping as a collection of confident sub-alignments, which is more tolerant of structural variation and more sensitive to paralog-specific variants (PSVs). Our [most recent paper](https://doi.org/10.1101/2020.11.01.363887) desribes this concept and benchmarking results.    

## Compile

Clone source code from master branch or download the [latest release](https://github.com/marbl/Winnowmap/releases/latest).
  ```sh
	git clone https://github.com/marbl/Winnowmap.git
  ```
Winnowmap compilation requires C++ compiler with c++11 and openmp, which are available by default in GCC >= 4.8.
  ```sh
	cd Winnowmap
	make -j8
  ```
Expect `winnowmap` and `meryl` executables in `bin` folder.

## Usage

For either mapping long reads or computing whole-genome alignments, Winnowmap requires pre-computing high frequency k-mers (e.g., top 0.02% most frequent) in a reference. Winnowmap uses [meryl](https://github.com/marbl/meryl) k-mer counting tool for this purpose.  

*  Mapping ONT or PacBio-hifi WGS reads
  ```sh
	meryl count k=15 output merylDB ref.fa
	meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

	winnowmap -W repetitive_k15.txt -ax map-ont ref.fa ont.fq.gz > output.sam  [OR]
	winnowmap -W repetitive_k15.txt -ax map-pb ref.fa hifi.fq.gz > output.sam
  ```

*  Mapping genome assemblies

  ```sh
	meryl count k=19 output merylDB asm1.fa
	meryl print greater-than distinct=0.9998 merylDB > repetitive_k19.txt

	winnowmap -W repetitive_k19.txt -ax asm20 asm1.fa asm2.fa > output.sam
  ```
  For the genome-to-genome use case, it may be useful to visualize the dot plot. This [perl script](https://github.com/marbl/MashMap/blob/master/scripts) can be used to generate a dot plot from [paf](https://github.com/lh3/miniasm/blob/master/PAF.md)-formatted output. In both usage cases, pre-computing repetitive k-mers using [meryl](https://github.com/marbl/meryl) is quite fast, e.g., it typically takes 2-3 minutes for the human genome reference.

## Benchmarking

When comparing Winnowmap (v1.0) to minimap2 (v2.17-r954), we observed a reduction in the mapping error-rate from 0.14% to 0.06% in the recently finished [human X chromosome](https://github.com/nanopore-wgs-consortium/CHM13), and from 3.6% to 0% within the highly repetitive X centromere (3.1 Mbp). Winnowmap improves mapping accuracy within repeats and achieves these results with sparser sampling, leading to better index compression and competitive runtimes. By avoiding masking, we show that Winnowmap maintains uniform minimizer density.

<p align="center">
<img src="https://i.postimg.cc/MKtqBYPn/readme-winnowmap-density.jpg" width=400px"> <br>
Minimizer sampling density using a human X chromosome as the reference, with the centromere positioned between 58 Mbp and 61 Mbp. ‘Standard’ method refers to the classic minimizer sampling algorithm from <a href="http://www.cs.toronto.edu/~wayne/research/papers/minimizers.pdf">Roberts et al.</a>, without any masking or modification.
</p>

## Publications

- **Chirag Jain, Arang Rhie, Nancy Hansen, Sergey Koren and Adam Phillippy**. "[Long-read mapping to repetitive reference sequences using Winnowmap2](https://doi.org/10.1038/s41592-022-01457-8)". *Nature Methods*, 2022.
- **Chirag Jain, Arang Rhie, Haowen Zhang, Chaudia Chu, Brian Walenz, Sergey Koren and Adam Phillippy**. "[Weighted minimizer sampling improves long read mapping](https://doi.org/10.1093/bioinformatics/btaa435)". *Bioinformatics (ISMB proceedings)*, 2020.
