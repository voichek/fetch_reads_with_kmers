# fetch_reads_with_kmers
A small program to fetch reads from paired-end sequences which contain predefined set of k-mers

usage:

`./fetch_reads <R1.fastq.gz> <R2.fastq.gz> <kmers.fasta> <kmer length> <output base name>`

* k-mers should have length <= 31
* read files can be compressed (.gz) or not
* kmers should be in fasta format, there is an example in `example_kmer_file.fasta`

