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
DEMO	= hot_tail_demo

# ----- Required modules in order of dependence -----------------------|
MODS 	= 	double \
			physical_constants \
			math \
			var_args \
			file_io \
			Coulomb_logarithms \
			collision_frequencies \
			electric_fields \
			calculate_Dreicer_growthrate \
			calculate_hot_tail_population \
			calculate_avalanche_growthrate

MODOBJ 	= $(patsubst %, $(OBJDIR)/%.o, $(MODS))

# ----- Compile fortran files -----------------------------------------|
$(OBJDIR)/%.o: $(SRCDIR)/%.f
	@# Compiles a single .f-file into an object
	@mkdir -p $$(dirname $@)
	@mkdir -p $(MODDIR)
	$(FC) -c -module $(MODDIR) -o $@ $^

$(REG):	$(MODOBJ) $(OBJDIR)/$(REG).o
	@# Creates the main module of this project
	@cp $(MODDIR)/$(REG).mod ./.

$(DEMO): $(MODOBJ) $(OBJDIR)/$(REG).o $(OBJDIR)/$(DEMO).o
	@# Builds the demonstration program
	$(FC) $(FFLAGS) -o demo/$@ $^

.PHONY: all
all:
	make $(DEMO)

# -----	Clean up ------------------------------------------------------|
.PHONY: clean
clean:
	@rm -rf $(MODDIR)
	@rm -rf $(OBJDIR)
	   
# ----- DONE ----------------------------------------------------------|
