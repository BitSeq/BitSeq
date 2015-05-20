Binary .prob file format specification
======================================
version 0.1, 2015-05-20
-----------------------

HEADER
======

char[6] "bprob\0"
int32   file format version = 0
int64   Nmap      // Number of mapped reads
int64   Nall      // Total number of reads
int64   Naligns   // Total number of alignments
int64   M         // Number of transcripts
int64   G         // Number of genes

int64   pointer to transcript_names[0]
int64   pointer to gene_names[0]
int64   pointer to transcript_genes[0]
int64   pointer to transcript_len[0]
int64   pointer to transcript_efflen[0]
int64   pointer to read_names[0]
int64   pointer to read_first_alignment[0]
int64   pointer to alignment_tid[0]
int64   pointer to alignment_tr[0]


DATA
====

(char *)[M]      transcript_names (null-terminated strings)
(char *)[G]      gene_names (null-terminated strings)
int64[M]         transcript_genes   // pointers to gene_names
float[M]         transcript_len
float[M]         transcript_efflen
(char *)[Nmap]   read_names (null-terminated strings)
int64[Nmap]      read_first_alignment
int64[Naligns]   alignment_tid
double[Naligns]  alignment_pr
