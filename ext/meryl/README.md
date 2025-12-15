# Meryl

This is 'meryl', a near total rewrite of 'meryl' that appeared in both
[project kmer](http://kmer.sourceforge.net/) and
[Celera Assembler](http://wgs-assembler.sourceforge.net/).

*IMPORTANT*: Get the latest meryl code from this repo. This is not compatible with old meryl dbs built from canu 1.8 or earlier. The new meryl is significantly faster than the previous version.

meryl dbs are no longer in `.mcdat` and `.mcidx` file format. Meryl db is now designed as a DIRECTORY, containing 64 binaries + 64 indexes (128 files).

### Dependency
* gcc 10.2.0 or higher

### Installation

* Release version: download a stable [release](https://github.com/marbl/meryl/releases/tag/v1.4.1) version
```shell
# Example for Linux-amd64
wget https://github.com/marbl/meryl/releases/download/v1.4.1/meryl-1.4.1.Linux-amd64.tar.xz
tar -xJf meryl-1.4.1.Linux-amd64.tar.xz
export PATH=/path/to/meryl-1.4.1/bin:$PATH
``` 

* Experimental tip (use git 2.25.1 or higher):
```shell
git clone https://github.com/marbl/meryl.git

# build
cd meryl/src
make -j 24
export PATH=/path/to/meryl/*/bin:$PATH
```

## Evaluate assemblies with k-mers and more

See [Merqury](https://github.com/marbl/merqury).

## Citing Meryl

We didn't want to suffer the world with yet another k-mer counting paper, so we stuffed meryl into the merqury methods:
>Rhie, A., Walenz, B.P., Koren, S. et al. Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies. Genome Biol 21, 245 (2020). https://doi.org/10.1186/s13059-020-02134-9
