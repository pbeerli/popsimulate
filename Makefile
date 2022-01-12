# simulator suite
# Peter Beerli 2004, Tallahassee
# beerli@csit.fsu.edu
#-----------------------------------------------------
.SUFFIXES:
.SUFFIXES: .o .c .h

prefix = /usr/local
exec_prefix = ${prefix}
INSTALLDIR = ${exec_prefix}/bin #/usr/local/bin
MANPAGEDIR = ${prefix}/man/man1 #/usr/local/man/man1
#----------------------------------------------------
CC   =  cc
INSTALL = /usr/bin/installbsd -c
CFLAGS = -O3 -ffast-math -funroll-loops -finline-functions -fstrength-reduce -frerun-cse-after-loop -fomit-frame-pointer -Wall  -DTESTING -DHIGHBITS  -DHAVE_LIBM=1 -DSTDC_HEADERS=1 -DHAVE_LIMITS_H=1 -DHAVE_STRINGS_H=1 -DRETSIGTYPE=void -DHAVE_STRFTIME=1 -DHAVE_STRCSPN=1 -DHAVE_LGAMMA=1  -DHAVE_MALLOCWRAP
DEBUG_CFLAGS = -g -Wall  -DTESTING  -DHIGHBITS  -DHAVE_LIBM=1 -DSTDC_HEADERS=1 -DHAVE_LIMITS_H=1 -DHAVE_STRINGS_H=1 -DRETSIGTYPE=void -DHAVE_STRFTIME=1 -DHAVE_STRCSPN=1 -DHAVE_LGAMMA=1  -DHAVE_MALLOCWRAP
OTHERLIBS =   -lm
#----------------------------------------------------
#
CFILES= migtree.c migdata.c
HFILES= definitions.h  migdata.h  migration.h  random.h  ranlib.h  sighandler.h  sort.h  tools.h
RANLIBFILES = randomlib/src/com.c  randomlib/src/ranlib.c randomlib/linpack/linpack.c
RANLIBOFILES = com.o  ranlib.o linpack.o
SOURCEFILES = $(CFILES)
ALLSOURCES = $(CFILES) $(HFILES) Makefile README
OFILES = $(CFILES:.c=.o)
INSTALL_DEPENDS = $(INSTALLDIR) all 
PRODUCT_DEPENDS = $(OFILES)

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o
#
all::
#	$(MAKE) ranlib
	$(MAKE) migtree
	$(MAKE) migdata

#$(MAKE) solexa
#$(MAKE) solexadate

migtree::   
	$(CC) -g migtree.c $(RANLIBFILES) -o migtree  -lm
migtree_old::   
	$(CC) -g migtree_old.c $(RANLIBFILES) -o migtree_old  -lm
bottleneck::   
	$(CC) -g -arch i386 -m32 bottleneck.c $(RANLIBFILES) -o bottleneck  -lm
#	dsymutil bottleneck
migdata::   
	$(CC) -g migdata.c -o migdata  -lm

migdatanew::   
	$(CC) -g migdatanew.c -o migdatanew  -lm

solexa::   
	$(CC) -g solexa.c -o solexa  -lm

solexadate::   
	$(CC) -g solexa.c -o solexadate  -lm

migdata_old::   
	$(CC) -g migdata-old.c -o migdata_old  -lm

ranlib::
	cd randomlib;
	$(CC) -O -c $(RANLIBFILES)
	ar -r librandom.a $(RANLIBOFILES)
	cp librandom.a .. ;
	cd .. ;
	ranlib librandom.a 
install:: $(INSTALL_DEPENDS)
	$(INSTALL) $(IFLAGS) $(NAME) $(INSTALLDIR)
	$(INSTALL) $(IFLAGS) $(MANPAGE) $(MANPAGEDIR)

clean::
	-/bin/rm -rf  $(GARBAGE) *.o 

zip::
	-/bin/rm -i beerlisim.tar.gz 
	mkdir /tmp/beerlisim
	tar cf -  $(ALLSOURCES) randomlib | (cd /tmp/beerlisim; tar xf - )
	(cd /tmp; tar cf - beerlisim) | gzip > beerlisim.tar.gz
	rm -rf /tmp/beerlisim

help::
	@echo '  all (the default)'
	@echo '  migtree'
	@echo '  migdata'
	@echo '  ranlib '
	@echo '  clean'
	@echo '  install'

#----------------------------------------------------------
#
$(OFILE_DIR) $(INSTALLDIR):
	mkdir $@;

SRCROOT:
	@if [ -n "${$@}" ]; then exit 0; \
	else echo Must define $@; exit 1; fi








