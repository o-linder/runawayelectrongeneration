# ---------------------------------------------------------------------|
#                    Written for GNU make 4.2.1                        |
# ---------------------------------------------------------------------|

# ----- Set up compiler and flags -------------------------------------|
FC = /afs/@cell/common/soft/intel/ics2013/14.0/bin/intel64/bin/ifort
FFLAGS  = -u -m64 -module $(MODDIR) -O -openmp

# ----- Declare build directories -------------------------------------|
MODDIR	= mod/
OBJDIR	= obj/
SRCDIR	= src/

# ----- Main modules and demonstration program ------------------------|
REG		= runawayelectrongeneration
REGC	= $(REG)_complete
DEMO	= hot_tail_demo
MODA	= mods.tar.gz

# ----- Phony rules ---------------------------------------------------|
.PHONY : all clean module

all :
	make module
	make $(DEMO)

clean :
	rm -f $(REGC).o $(MODA)
	rm -rf $(MODDIR)
	rm -rf $(OBJDIR)

module :
	make $(REGC).o

# ----- Main rule to compile fortran files ----------------------------|
$(OBJDIR)%.o : \
		$(SRCDIR)%.f
	@# Compiles a singe .f-file into an object
	@mkdir -p $(dir $@)
	@mkdir -p $(MODDIR)
	$(FC) $(FFLAGS) -c -o $@ $<

# ----- Combine all objects and modules into single files -------------|
$(REGC).o : $(OBJDIR)$(REG).o
	@# Creates a single object from all required objects and stores it
	@# and an archive of the required modules in the main directory
	@tar czf $(MODA) $(MODDIR)*.mod
	ld -relocatable $(OBJDIR)*.o -o $@
	cp $@ $(OBJDIR)$@

# ----- Hot-tail demonstration program --------------------------------|
$(DEMO) : $(REGC).o $(OBJDIR)$(DEMO).o
	@# Builds the demonstration program
	$(FC) $(FFLAGS) -o demo/$@ $^

# ----- List prerequisites --------------------------------------------|
#	# Runawayelectrongeneration module
$(OBJDIR)$(REG).o : \
		$(OBJDIR)calculate_hot_tail_population.o \
		$(OBJDIR)calculate_Dreicer_growthrate.o \
		$(OBJDIR)calculate_avalanche_growthrate.o \
		$(OBJDIR)electric_fields.o \
		$(OBJDIR)collision_frequencies.o \
		$(OBJDIR)Coulomb_logarithms.o \

#	# Module for hot-tail calculations
$(OBJDIR)calculate_hot_tail_population.o : \
		$(OBJDIR)double.o \
        $(OBJDIR)physical_constants.o \
        $(OBJDIR)file_io.o \
        $(OBJDIR)Coulomb_logarithms.o \
        $(OBJDIR)collision_frequencies.o \

#	# Module for Dreicer generation
$(OBJDIR)calculate_Dreicer_growthrate.o : \
		$(OBJDIR)double.o \
		$(OBJDIR)collision_frequencies.o \
		$(OBJDIR)electric_fields.o \
		$(SRCDIR)inc/parameters_CODE_neural_network.inc

#	# Module for avalanche generation
$(OBJDIR)calculate_avalanche_growthrate.o : \
		$(OBJDIR)double.o \
		$(OBJDIR)physical_constants.o \
		$(OBJDIR)Coulomb_logarithms.o \
		$(OBJDIR)collision_frequencies.o \
		$(OBJDIR)electric_fields.o \

#	# Module for electric fields
$(OBJDIR)electric_fields.o : \
		$(OBJDIR)double.o \
		$(OBJDIR)physical_constants.o \
		$(OBJDIR)Coulomb_logarithms.o \
		$(OBJDIR)collision_frequencies.o

#	# Module for collision frequencies
$(OBJDIR)collision_frequencies.o : \
		$(OBJDIR)double.o \
		$(OBJDIR)physical_constants.o \
		$(OBJDIR)Coulomb_logarithms.o \
		$(SRCDIR)inc/parameters_get_length_scale.inc \
		$(SRCDIR)inc/parameters_get_Ihat.inc

#	# Module for Coulomb logarithms
$(OBJDIR)Coulomb_logarithms.o : \
		$(OBJDIR)double.o \
		$(OBJDIR)physical_constants.o

#	# Module for file I/O
$(OBJDIR)file_io.o : \
		$(OBJDIR)double.o

#	# Module for physical constants
$(OBJDIR)physical_constants.o : \
		$(OBJDIR)double.o

# ----- Done ----------------------------------------------------------|
