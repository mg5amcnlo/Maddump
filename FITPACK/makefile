# ----------------------------------------------------------------------------
#
# Makefile for fitpack library
#
# ----------------------------------------------------------------------------

# Check for ../make_opts
ifeq ($(wildcard ../make_opts), ../make_opts)
  include ../make_opts
else
  FFLAGS+= -ffixed-line-length-132
  FC=gfortran
endif

FITPACK = bispev.o curfit.o dblint.o fpback.o fpbisp.o fpbspl.o fpchec.o \
	  fpcurf.o fpdisc.o fpgivs.o fpintb.o fpknot.o fporde.o fprank.o \
	  fprati.o fprota.o fpsurf.o profil.o splev.o splint.o surfit.o

LIBDIR=../../lib/
LIBRARY=libfitpack.$(libext)

all: $(LIBDIR)$(LIBRARY)


$(LIBRARY): $(FITPACK)
	ar cru $(LIBRARY) $(FITPACK)
	ranlib $(LIBRARY)

$(LIBDIR)$(LIBRARY): $(FITPACK)
	$(call CREATELIB, $@, $^)

clean: 
	$(RM) *.o $(LIBDIR)$(LIBRARY)
