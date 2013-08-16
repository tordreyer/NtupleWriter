# Package information
LIBRARY = Ntuple
OBJDIR  = obj
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = include

INCLUDES += -I$(FASTJETDIR)/../include

# execute configure.sh if necessary:
dummy := $(shell [ -d ../UHHAnalysis ] || ./configure.sh)

# Include the generic compilation rules
include $(SFRAME_DIR)/Makefile.common
