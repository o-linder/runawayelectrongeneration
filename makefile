.KEEP_STATE:

# -----	Set up compiler and flags -------------------------------------|
FC = /afs/@cell/common/soft/intel/ics2013/14.0/bin/intel64/bin/ifort
FFLAGS  = -u -m64 -O -openmp 

# ----- Declare build directories -------------------------------------|
MODDIR  = ./mod
OBJDIR	= ./obj
SRCDIR  = ./src

# ----- Main modules and demonstration program ------------------------|
REG		= runawayelectrongeneration
REGC	= $(REG)_complete
DEMO	= hot_tail_demo

# ----- Required modules in order of dependence -----------------------|
OBJS 	= 	double \
			physical_constants \
			file_io \
			Coulomb_logarithms \
			collision_frequencies \
			electric_fields \
			calculate_hot_tail_population \
			calculate_Dreicer_growthrate \
			calculate_avalanche_growthrate

OBJSPTH = $(patsubst %, $(OBJDIR)/%.o, $(OBJS))

# ----- Compile fortran files -----------------------------------------|
$(OBJDIR)/%.o: $(SRCDIR)/%.f
	@# Compiles a single .f-file into an object
	@mkdir -p $$(dirname $@)
	@mkdir -p $(MODDIR)
	$(FC) -c -module $(MODDIR) -o $@ $^

$(REGC).o: $(OBJSPTH) $(OBJDIR)/$(REG).o
	@# Creates a single object from all required objects and stores it
	@# and the main module of this project in the main directory
	@cp $(MODDIR)/$(REG).mod ./.
	ld -relocatable $^ -o $@
	cp $@ $(OBJDIR)/$@

$(DEMO): $(REGC).o $(OBJDIR)/$(DEMO).o
	@# Builds the demonstration program
	$(FC) $(FFLAGS) -o demo/$@ $^

.PHONY: project
project:
	make $(REGC).o

.PHONY: all
all:
	make project
	make $(DEMO)

# -----	Clean up ------------------------------------------------------|
.PHONY: clean
clean:
	@rm -f $(REGC).o $(REG).mod
	@rm -rf $(MODDIR)
	@rm -rf $(OBJDIR)
	   
# ----- Done ----------------------------------------------------------|
