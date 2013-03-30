CC=		gcc
# CFLAGS=		-Wall -O3 #-m64 #-arch ppc
# Edit: we don't need warnings and prefere -O2 optimization, with architecture specific tuning.
CFLAGS=		-Wno-unused-result -O2 -mtune=generic
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_CURSES_LIB=1
INCLUDES=	-I.

.c.o:
	$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all: \
   sam.o \
   bam.o \
   sam_header.o \
   kstring.o \
   faidx.o \
   razf.o \
   bgzf.o \
   bam_import.o \
   bam_aux.o \
   bam_pileup.o

clean: 
	rm *.o
