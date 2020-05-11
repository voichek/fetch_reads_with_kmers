///
///      @file  fetch_reads.cpp
///     @brief  Giving a list of k-mers and reads file, filter the reads that have one of the
//		k-mers in them
///
///    @author  Yoav Voichek (YV), yoav.voichek@tuebingen.mpg.de
///
///  @internal
///    Created  01/20/19
///   Compiler  gcc/g++
///    Company  Max Planck Institute for Developmental Biology Dep 6
///  Copyright  Copyright (c) 2019, Yoav Voichek
///
/// This source code is released for free distribution under the terms of the
/// GNU General Public License as published by the Free Software Foundation.
///=====================================================================================
///


#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>

#include "./include/klib/khash.h"
#include "./include/klib/kseq.h"

// https://bioinformatics.stackexchange.com/questions/5359/what-is-the-most-compact-data-structure-for-canonical-k-mers-with-the-fastest-lo?noredirect=1&lq=1
static inline uint64_t hash_64(uint64_t key)
{ // more sophisticated hash function to reduce collisions
	key = (~key + (key << 21)); // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)); // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)); // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31));
	return key;
}

KHASH_INIT(64, khint64_t, int, 1, hash_64, kh_int64_hash_equal)

unsigned char seq_nt4_table[128] = { // Table to change "ACGTN" to 01234
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


static void insert_seq(khash_t(64) *h, int k, int len, char *seq)
{
	int i, l;
	uint64_t x[2], mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2, tot = 0;
	for (i = l = 0, x[0] = x[1] = 0; i < len; ++i) {
		int c = (uint8_t)seq[i] < 128? seq_nt4_table[(uint8_t)seq[i]] : 4;
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
			if (++l >= k) { // we find a k-mer
				khint_t itr;
				int absent;
				itr = kh_put(64, h, x[0], &absent);
				if (absent) kh_val(h, itr) = 0;
				tot += ++kh_val(h, itr);
				
				itr = kh_put(64, h, x[1], &absent);
				if (absent) kh_val(h, itr) = 0;
				tot += ++kh_val(h, itr);
			}
		} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
	}
}
static bool look_up(khash_t(64) *h, int k, int len, char *seq)
{
	int i, l;
	uint64_t x, mask = (1ULL<<k*2) - 1, shift = (k - 1) * 2;
	for (i = l = 0, x = 0; i < len; ++i) {
		int c = (uint8_t)seq[i] < 128? seq_nt4_table[(uint8_t)seq[i]] : 4;
		if (c < 4) { // not an "N" base
			x = (x << 2 | c) & mask;                  // forward strand
			if (++l >= k) { // we find a k-mer
				khint_t itr;
				itr = kh_get(64, h, x);
				if(itr != kh_end(h))
					return true;
			}
		} else l = 0, x = 0; // if there is an "N", restart
	}
	return false;
}

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	if(argc != 6) {
		fprintf(stderr, "Usage: %s <R1.fastq.gz> <R2.fastq.gz> <kmers.fasta> <kmer length> <output base name>\n", argv[0]);
		return 1;
	}
	khash_t(64) *h; // hash table
	int k = atoi(argv[4]); // len of k-mer
	fprintf(stderr, "kmer length = %d\n", k);
	if(k<1 || k>32) {
		fprintf(stderr, "Error: kmers length should be in the range 1-31\n");
		return 0;
	}
	int i;

	h = kh_init(64); // 
	// Insert k-mers to hash table
	uint64_t cnt_insertions = 0;
	gzFile fp_kmers = gzopen(argv[3], "r");
	kseq_t *ks_kmers;
	ks_kmers = kseq_init(fp_kmers);
	kh_resize(64,h,(1<<16)); // I  don't get any results if I put this before inserting
	while(kseq_read(ks_kmers) >= 0) {
		insert_seq(h, k, ks_kmers->seq.l, ks_kmers->seq.s); 
		cnt_insertions++;
	}
	fprintf(stderr, "[Insert kmers]: %d kmers\n", int(cnt_insertions));
	kseq_destroy(ks_kmers);
	gzclose(fp_kmers);

	// Look for k-mers
	uint64_t cnt_found = 0;
	gzFile fp_R1 = gzopen(argv[1], "r");
	gzFile fp_R2 = gzopen(argv[2], "r");
	kseq_t *ks_R1;	     // read from R1
	kseq_t *ks_R2;	     // read from R2
	ks_R1 = kseq_init(fp_R1);
	ks_R2 = kseq_init(fp_R2);

	// Open output files
	FILE *out_R1, *out_R2;
	char R1_out_fn[256];
	char R2_out_fn[256];
	sprintf(R1_out_fn, "%s_R1.fastq", argv[5]);
	sprintf(R2_out_fn, "%s_R2.fastq", argv[5]);
	out_R1 = fopen(R1_out_fn, "w");
	out_R2 = fopen(R2_out_fn, "w");

	// We assume R1 & R2 are of the same size
	int counter = 0;
	while((kseq_read(ks_R1) >=0 ) && (kseq_read(ks_R2))>=0) {
		if((counter % 1000000) == 0)
			fprintf(stderr, "counter = %d\n", counter);
		counter++;
		bool hit_R1 = look_up(h, k, ks_R1->seq.l, ks_R1->seq.s);
		bool hit_R2 = look_up(h, k, ks_R2->seq.l, ks_R2->seq.s);
		if (hit_R1 || hit_R2) {
			cnt_found++;
			fprintf(out_R1, "@%s\n%s\n+\n%s\n", ks_R1->name.s, ks_R1->seq.s, ks_R1->qual.s);
			fprintf(out_R2, "@%s\n%s\n+\n%s\n", ks_R2->name.s, ks_R2->seq.s, ks_R2->qual.s);
		}
	}
	kseq_destroy(ks_R1);
	kseq_destroy(ks_R2);
	gzclose(fp_R1);
	gzclose(fp_R2);
	fclose(out_R1);
	fclose(out_R2);
	kh_destroy(64, h);
	fprintf(stderr, "[found reads]: %d reads\n", int(cnt_found));
	return 0;
}

