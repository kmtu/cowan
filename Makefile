# Defining variables ===============================
FC = gfortran
MODFLAG = -J
FCFLAGS = $(MODFLAG)$(MODDIR)

SRCDIR = ./src/
#OBJDIR = ./obj/
MODDIR = ./mod/
#LIBDIR = $(SRCDIR)lib/
OUTDIR = ./out/

PROGRAM = cowan

vpath %.f90 $(SRCDIR)
vpath %.o $(SRCDIR)
vpath % $(OUTDIR)

# Program ==================================
all : $(PROGRAM)
.PHONY : all 
# ---------------------------
cowan : cowan.o
	$(FC) -o $(OUTDIR)$(@F) $(SRCDIR)$(<F)

cowan.o : cowan.f90
	$(FC) -o $(SRCDIR)$(@F) -c $< $(FCFLAGS)

# Library ===================================
#nrtype.o : nrtype.f90
#	$(FC) -o $(LIBDIR)$(@F) -c $< $(FCFLAGS)
#==================================================
.PHONY: clean
clean :
	rm -f $(SRCDIR)*.o $(MODDIR)*.mod
	rm -f $(addprefix $(OUTDIR), $(PROGRAM) )
	rm -f *~ $(SRCDIR)*~
