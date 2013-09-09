CC=		gcc
# CFLAGS=		-Wall -O3 #-m64 #-arch ppc
# Edit: we don't need warnings and prefer -O2 optimization, with architecture specific tuning.
CFLAGS=		-O2 -mtune=generic -ffunction-sections -fdata-sections 
# removed -Wno-unused-result as it is not recognized prior to gcc 4.5.0
INCLUDES=	-I.

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

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
