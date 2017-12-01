.PHONY: ALL CLEAN
ALL: TARGETS_2D TARGETS_3D

# Build separate executables for 2D and 3D codes
SRCDIR = src
CFLAGS = ${PETSC_CC_INCLUDES}
CXXFLAGS = -std=c++11

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# Use specified source files
COMMON_SOURCES = $(SRCDIR)/IO_tools.cpp $(SRCDIR)/linear_sys.cpp $(SRCDIR)/field.cpp \
                 $(SRCDIR)/nonlocal_field.cpp

# 2D version
EXECUTABLE_2D=bin/solve_poisson_2D
SOURCES_2D = $(COMMON_SOURCES) $(SRCDIR)/main_2D.cpp $(SRCDIR)/poisson_solver_2D.cpp 
OBJECTS_2D = $(SOURCES_2D:%.cpp=%.o)
# Note: chkopts is important for PETSc!
TARGETS_2D: $(EXECUTABLE_2D) chkopts

$(EXECUTABLE_2D) : $(OBJECTS_2D)
	# Make bin directory if necessary
	mkdir -p bin/
	${CLINKER} ${CXXFLAGS} -o $@ $^ ${PETSC_LIB}
	
# 3D version
EXECUTABLE_3D=bin/solve_poisson_3D
SOURCES_3D = $(COMMON_SOURCES) $(SRCDIR)/main_3D.cpp $(SRCDIR)/poisson_solver_3D.cpp 
OBJECTS_3D = $(SOURCES_3D:%.cpp=%.o)
# Note: chkopts is important for PETSc!
TARGETS_3D: $(EXECUTABLE_3D) chkopts

$(EXECUTABLE_3D) : $(OBJECTS_3D)
	# Make bin directory if necessary
	mkdir -p bin/
	${CLINKER} ${CXXFLAGS} -o $@ $^ ${PETSC_LIB}
	
# Use uppercase to avoid conflicting with PETSc's built-in make clean...
CLEAN:
	$(RM) bin/*
	$(RM) $(OBJECTS_2D) $(OBJECTS_3D)