# Definizione delle variabili
ifndef FF_ROOT 
FF_ROOT = ./fastflow
endif
CXX = g++ -std=c++2a
MPI = mpicxx -std=c++2a
OPTFLAGS = -Wall -O3 -ffast-math
INCLUDES = -I. -I./include -I $(FF_ROOT)
LIBS = -pthread -fopenmp

SRCDIR = src
OBJDIR = obj
LOGDIR = log

# every source inside SRCDIR
SOURCES = $(filter-out $(SRCDIR)/mpi%.cpp, $(filter-out $(SRCDIR)/no_%.cpp, $(wildcard $(SRCDIR)/*.cpp)))
SOURCES_MPI = $(filter-out $(SRCDIR)/no_%.cpp, $(wildcard $(SRCDIR)/mpi*.cpp))

# every possible makefile target
TARGET = $(notdir $(basename $(SOURCES)))
TARGET_MPI = $(notdir $(basename $(SOURCES_MPI)))

# every possible obj inside OBJDIR
OBJECTS = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES:.cpp=.o)))
OBJECTS_MPI = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES_MPI:.cpp=.o)))

# debug option if selected
ifdef DEBUG
    OPTFLAGS += -DDEBUG=$(DEBUG)
endif

# default target, compile all executable
all: $(TARGET) $(TARGET_MPI)

allmpi: $(TARGET_MPI)

allnmpi: $(TARGET)

# create the OBJDIR folder if not present
$(OBJDIR):
	mkdir -p $(OBJDIR)

# general rule to compile each source file in the corresponding obj file
$(OBJDIR)/mpi%.o: $(SRCDIR)/mpi%.cpp | $(OBJDIR)
	$(MPI) $(OPTFLAGS) -c $< -o $@ $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CXX) $(OPTFLAGS) $(INCLUDES) -c $< -o $@ $(LIBS)


# general rule to link each executable (target)
$(TARGET): %: $(OBJDIR)/%.o
	$(CXX) $(OPTFLAGS) $< -o $@ $(LIBS)

# Regola per collegare gli eseguibili MPI
$(TARGET_MPI): %: $(OBJDIR)/%.o
	$(MPI) $(OPTFLAGS) $< -o $@ $(LIBS)


# remove all obj and executable files
clean:
	-rm -fr $(OBJECTS) $(OBJECTS_MPI) $(TARGET) $(TARGET_MPI)
# remove the logs
cleanlog: 
	-rm -fr $(LOGDIR)/*.csv
# remove also the obj folder
cleanall: clean cleanlog
	-rm -fr $(OBJDIR) $(LOGDIR) 

.PHONY: all allmpi allnmpi clean cleanlog cleanall run
