CC=		gcc
# CFLAGS=		-Wall -O3 #-m64 #-arch ppc
# Edit: we don't need warnings and prefere -O2 optimization, with architecture specific tuning.
CFLAGS=		-Wno-unused-result -O2 -mtune=generic
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=1
INCLUDES=	-I.

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all: sam.o bam.o sam_header.o faidx.o bgzf.o knetfile.o kstring.o razf.o bam_aux.o bam_import.o bam_pileup.o

sam.o:sam.h bam.h
bam.o:bam.h razf.h bam_endian.h kstring.h sam_header.h
sam_header.o:sam_header.h khash.h
faidx.o:faidx.h razf.h khash.h
bam_import.o:bam.h kseq.h khash.h razf.h
bam_pileup.o:bam.h razf.h ksort.h

clean: 
	rm *.o
