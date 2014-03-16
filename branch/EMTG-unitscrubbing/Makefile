##########################################################
# Really New Makefile
#
# I deserve only 25% credit.  This is based on a project I work on that originates back to 1991
# 
#

PROJECT = EMTG

# The names of the sub-projects that the  executable depends on.  These names are used later on to define targets, location, libraries, etc.
DEPS =

# The names of the sub-projects that the  test executable depends on.
EXE_DEPS =

SRC = $(wildcard src/*.cpp)

LOCAL_O_FILES       = $(SRC:.cpp=.o)

LOCAL_TEST_O_FILES  =
LOCAL_LIB_NAME      = EMTGlib
LOCAL_LIB_DIR       = .
LOCAL_TEST_LIB_NAME =
LOCAL_TEST_LIB_DIR  =
LOCAL_INCLUDE_ROOT  = .
EXE_NAME            = emtg

#point these at your local places correctly.
# _DIR is the include dir, _LIB_DIR is the library path
GCC_INSTALL_DIR = /usr/bin
GCC_INCLUDE_DIR = $(GCC_INSTALL_DIR)/include
GCC_LIB_DIR = $(GCC_INSTALL_DIR)/lib
BOOST_DIR = /usr/local/include/boost
BOOST_LIB_DIR = /usr/local/lib 
SNOPT_DIR = /home/jaengla2/libraries/snopt7.2-11/snopt7/cppsrc
SNOPT_LIB_DIR = $(SNOPT_DIR)/../lib
CSPICE_DIR = cspice/include
CSPICE_LIB_DIR = cspice/lib




#These are the compilers of choice.  Change these to be MPICC or NVCC as appropriate for the target
CC          = /usr/lib64/openmpi/bin/mpicc
CCC         = /usr/lib64/openmpi/bin/mpic++
#F77         = $(MPIC_DIR)/bin/mpif77
CLINKER     = $(CC)
CCLINKER    = $(CCC)
#FLINKER     = $(MPIC_DIR)/bin/mpif90
#F90         = $(MPIC_DIR)/bin/mpif90
#F90LINKER   = $(MPIC_DIR)/bin/mpif90
MAKE        = make --no-print-directory
SHELL       = /bin/sh


#PROF = -ggdb -O2 
PROF = -O3 -DEMTG_MPI -D_EMTG_proprietary -DQUIET_SNOPT -DBACKGROUND_MODE -D_STONEAGECplusplus
ADDITIONAL_LIBS = 
ADDITIONAL_LIBS = -lcspice -lsnopt_cpp -lsnopt_c -lsnprint_c  -lf2c -lboost_filesystem -lboost_system -lboost_serialization -lboost_mpi
SED  = $(shell which sed)
GREP = $(shell which grep)




# Expand SRCDIR so that it turns into "-I <srcdir>" for each listed directory
INCLUDE_DIRS = $(foreach dep, $(UNIQUE_DEPS), -I$($(dep)_INCLUDE_ROOT)) -I$(SNOPT_DIR) -I$(LOCAL_INCLUDE_ROOT) -I$(BOOST_DIR) -I$(CSPICE_DIR) -I$(GCC_INCLUDE_DIR)

#C_FLAGS  = -DROM $(PROF) $(WARNINGS) $(INCLUDE_DIRS)
C_FLAGS = -DROM $(PROF) $(INCLUDE_DIRS) 
L_FLAGS  = $(foreach dir, $(OBJC_ROOT) $(DEP_LIB_DIRS), -L$(dir)) -L$(BOOST_LIB_DIR) -L$(GCC_LIB_DIR) -L$(SNOPT_LIB_DIR) -L$(CSPICE_LIB_DIR) $(PROF) $(ADDITIONAL_LIBS)

TEST_L_FLAGS  = $(foreach dir, $(OBJC_ROOT) $(TEST_DEP_LIB_DIRS), -L$(dir)) -L/usr/lib  -Xlinker --whole-archive $(TEST_REQUIRED_LIBS) -Xlinker --no-whole-archive

REQUIRED_LIBS = $(foreach dep, $(DEPS) LOCAL, $(if $($(dep)_LIB_NAME), -l$($(dep)_LIB_NAME)))
TEST_REQUIRED_LIBS = $(foreach dep, $(TEST_DEPS) LOCAL LOCAL_TEST, $(if $($(dep)_LIB_NAME), -l$($(dep)_LIB_NAME)))

.SUFFIXES: 
.SUFFIXES: .d .cc .cpp .c .o .cu

ASSEMBLE_TARGETS = $(foreach dep,$(DEPS),$(dep)_LIB_NAME)

ASSEMBLE_TEST_TARGETS = $(foreach dep,$(TEST_DEPS),$(dep)_LIB_NAME)

CLEAN_TARGETS = $(foreach dep, $(UNIQUE_DEPS), $(dep)_CLEAN)

LOCAL_D_FILES = $(LOCAL_O_FILES:.o=.d)

#LOCAL_TEST_D_FILES = $(LOCAL_TEST_O_FILES:.o=.d) $(TEST_RUNNER:.o=.d)

DEP_LIB_DIRS = $(foreach dep, $(DEPS) LOCAL, $($(dep)_LIB_DIR))

#TEST_DEP_LIB_DIRS = $(foreach dep, $(TEST_DEPS) LOCAL LOCAL_TEST, $($(dep)_LIB_DIR))

UNIQUE_ASSEMBLE_TARGETS = $(sort $(ASSEMBLE_TARGETS) $(ASSEMBLE_TEST_TARGETS))

UNIQUE_DEPS = $(sort $(DEPS) $(TEST_DEPS))

LOCAL_LIB = $(if $(LOCAL_LIB_NAME), $(LOCAL_LIB_DIR)/lib$(LOCAL_LIB_NAME).a)

#LOCAL_TEST_LIB = $(if $(LOCAL_TEST_LIB_NAME), $(LOCAL_TEST_LIB_DIR)/lib$(LOCAL_TEST_LIB_NAME).a)

EXE_TARGETS = $(foreach dep,$(EXE_DEPS),$(dep)_EXE_NAME)

INSTALL_TARGETS = $(foreach dep,$(EXE_DEPS),$(dep)_INSTALL)

export PROJECTS += $(PROJECT)

# PHONY targets get remade even if there exists an up-to-date file with the same name.
.PHONY: default clean compile assemble link submit 

# Targets for mortals
default: link

clean: $(CLEAN_TARGETS) LOCAL_CLEAN

compile: $(DEPS) $(LOCAL_D_FILES) $(LOCAL_O_FILES)

assemble: compile $(ASSEMBLE_TARGETS) $(LOCAL_LIB)

link: assemble $(EXE_TARGETS) $(EXE_NAME)

install: $(INSTALL_TARGETS) $(if $(EXE_NAME), LOCAL_INSTALL)

submit: link
	qsub $(PROJECT).qsub
	

# Targets for make

# Invoking sub projects
$(UNIQUE_DEPS): MAKE_TARGET = $(if $(filter %_TEST, $@), compile-tests, compile)
$(UNIQUE_DEPS):
	$(if $(findstring $(DEP_NAME), $(PROJECTS)),,cd $($@_DIR) && $(MAKE) $(MAKE_TARGET))

# First, remove the _LIB_NAME attached by ASSEMBLE_TARGETS, a
$(UNIQUE_ASSEMBLE_TARGETS): DEP_NAME = $(@:_LIB_NAME=)
$(UNIQUE_ASSEMBLE_TARGETS): MAKE_TARGET = $(if $(filter %_TEST, $(DEP_NAME)), assemble-tests, assemble)
$(UNIQUE_ASSEMBLE_TARGETS):
	$(if $(findstring $(DEP_NAME), $(PROJECTS)),,cd $($(DEP_NAME)_DIR) && $(MAKE) $(MAKE_TARGET))

# First, remove the _EXE_NAME attached by EXE_TARGETS, a
$(EXE_TARGETS): DEP_NAME = $(@:_EXE_NAME=)
$(EXE_TARGETS):
	$(if $(findstring $(DEP_NAME), $(PROJECTS)),,cd $($(DEP_NAME)_DIR) && $(MAKE) link)

$(CLEAN_TARGETS): DEP_NAME = $(@:_CLEAN=)
$(CLEAN_TARGETS):
	$(if $(findstring $(DEP_NAME), $(PROJECTS)),,cd $($(DEP_NAME)_DIR) && $(MAKE) clean)

$(INSTALL_TARGETS): DEP_NAME = $(@:_INSTALL=)
$(INSTALL_TARGETS):
	$(if $(findstring $(DEP_NAME), $(PROJECTS)),,cd $($(DEP_NAME)_DIR) && $(MAKE) install)

#Local stuff

# The rule to change either a '.c' or a '.m' to a '.o'
%.o : %.c %.d %.cc %.cpp
.cc.o .cpp.o .c.o .cu.o:
	@echo "Compiling \"$<\""
	@$(CCC) -c $(C_FLAGS) $< -o $@

.cc.d .cpp.d .c.d .cu.d:
	@echo "Creating dependency for \"$<\""
#	@$(CCC) $(WARNINGS) $(INCLUDE_DIRS) -MM $< -o $@
# This foogly hack because gcc seems to have issues with emitting the correct target/dependency names of
# files in sub-dirs of the current dir (specifically, it doesn't add the sub-dir to the target
# name, and insists on adding the directory to the dependencies) which ends up breaking dependencies...)
	@dependLine=$$( $(CCC) $(C_FLAGS) $(INCLUDE_DIRS) -MM $< ); \
	dirName=$$( dirname $< | $(SED) "s/\//\\\\\//g" ); \
	dependLine=$$( echo $${dependLine} | $(SED) "s/ $${dirName}\// /g" ); \
	oFile=$$( echo $${dependLine} | $(SED) "s/:.*//" ); \
	dependencies=$$( echo $${dependLine} | $(SED) "s/.*://" ); \
	echo $${oFile} $${oFile%.o}.d: $${dependencies} | $(SED) "s/ \\\//g" > $@

$(LOCAL_LIB): compile $(ASSEMBLE_TARGETS)
	@echo Assembling $(LOCAL_LIB_NAME)...
	@ar rs $(LOCAL_LIB) $(LOCAL_O_FILES)

# Create the executable
$(EXE_NAME): assemble
	@echo Linking $(EXE_NAME)...
	@$(CCC) $(LOCAL_LIB) $(L_FLAGS) -o $(EXE_NAME)

# Erase all object files, the dependencies file, the core file and the executable, then rebuild everything
LOCAL_CLEAN:
	@echo Cleaning $(PROJECT)...
	@rm -f $(LOCAL_O_FILES) $(LOCAL_TEST_O_FILES) $(LOCAL_LIB) $(LOCAL_TEST_LIB) $(EXE_NAME) $(TEST_EXE_NAME) $(LOCAL_D_FILES) $(LOCAL_TEST_D_FILES) $(TEST_RUNNER) $(TEST_RUNNER:.o=.d) core*


ifeq (,$(findstring clean,$(MAKECMDGOALS)))
-include $(LOCAL_O_FILES:.o=.d) $(LOCAL_TEST_O_FILES:.o=.d)
endif

