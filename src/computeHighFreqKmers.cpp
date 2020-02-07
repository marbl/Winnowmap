//g++ -std=c++11 computeHighFreqKmers.cpp -lz -o computeHighFreqKmers

//Input parameters
// <arg1> : kmer length 
// <arg2> : frequency threshold x  (any kmer with frequency >= x will be printed)
// <arg3> : 0: no compression, 1: homopolymer compression
// <arg4> : input fasta/fastq file
// <arg5> : output file name
// Goal   : output list of highly repetitive kmers
//          each ouptut kmer is encoded with 2k bits, using similar encoding as minimap2
//          each output kmer is represented by its canonical representation in the output


#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

//mapping characters to 2-bit (copied from minimap2)
	unsigned char seq_nt4_table[256] = {
		0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
		4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
	};

/**
 * @brief                     helper function to add kmer counts into the table
 * @param[in]     len         length of sequence
 * @param[in]     str         sequence
 * @param]in]     k           kmer length
 * @param[in]     hpc         compression mode
 * @param[in/out] table       kmer frequence table
 */
void doCounting (int len, const char *str, int k, int hpc, std::unordered_map <uint64_t, uint32_t> &table)
{
	assert(len > 0 && (k > 0 && k <= 28)); 
	uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, kmer[2] = {0,0};

	for (int i = 0, l = 0; i < len; ++i) 
	{
		int c = seq_nt4_table[(uint8_t)str[i]];
		if (c < 4) 
		{ 
			// not an ambiguous base
			int z;
			if (hpc) 
			{
				int skip_len = 1;
				if (i + 1 < len && seq_nt4_table[(uint8_t)str[i + 1]] == c) {
					for (skip_len = 2; i + skip_len < len; ++skip_len)
						if (seq_nt4_table[(uint8_t)str[i + skip_len]] != c)
							break;
					i += skip_len - 1; // put $i at the end of the current homopolymer run
				}
			} 

			kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
			kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer
			if (kmer[0] == kmer[1]) continue;  // skip "symmetric k-mers" as we don't know their strand
			z = kmer[0] < kmer[1]? 0 : 1; // strand
			++l;

			if (l >= k) 
			{
				if (table.find(kmer[z]) == table.end())
				{
					table.emplace(kmer[z], 1);  //frequency = 1
				}
				else
				{
					table[kmer[z]] += 1;
				}
			}
		}
		else
			l=0;
	}
}

int main(int argc, char *argv[])
{

	if (argc != 6)
	{
		fprintf(stderr, "%s expects five argument parameters\n", argv[0]);
		fprintf(stderr, "Usage: %s [k] [homopolymer compr: 0/1] [freq cutoff] ref.fa output.txt\n", argv[0]);
		exit(1);
	}

	int k = atoi(argv[1]);
	int hpc = atoi(argv[2]);
	int cutoff = atoi(argv[3]);
	const char* filename = argv[4];
	const char* outfilename = argv[5];

	assert (hpc == 0 || hpc == 1);


	std::unordered_map <uint64_t, uint32_t> kmerFreqTable;

	std::cout << "kmer size = " << k << std::endl;
	std::cout << "hpc on/off? = " << hpc << std::endl;
	std::cout << "starting frequency = " << cutoff << std::endl;
	std::cout << "input filename = " << filename << std::endl;
	std::cout << "output filename = " << outfilename << std::endl;

	//READ
	{
		// open query file for reading; you may use your favorite FASTA/Q parser
		gzFile f = gzopen(filename, "r");
		assert(f);
		kseq_t *ks = kseq_init(f);

		while (kseq_read(ks) >= 0) 
		{ 
			doCounting(ks->seq.l, ks->seq.s, k, hpc, kmerFreqTable);
		}

		kseq_destroy(ks); // close the query file
		gzclose(f);
	}

	//STATISTICS
	{
		std::cout << "#entries in the table = " << kmerFreqTable.size() << std::endl;
	}

	//WRITE
	{
		std::ofstream outFile (outfilename);
		std::size_t kmers_printed=0;

		for (const auto &e : kmerFreqTable) 
		{
			if ((int)(e.second) >= cutoff)  
			{
				outFile << e.first << "\n";
				kmers_printed++;
			}
		}
		std::cout << "finished writing " << kmers_printed << " kmers into the file" << std::endl;
	}

	return 0;
}

