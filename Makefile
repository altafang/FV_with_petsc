EXECUTABLE=bin/solve_poisson

# Automatically detect source files
# Adapted from: http://hiltmon.com/blog/2013/07/03/a-simple-c-plus-plus-project-structure/
SRCDIR = src
SRCEXT = cpp
SOURCES = $(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
OBJECTS = $(patsubst $(SRCDIR)/%,$(SRCDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

# Note: chkopts is important for PETSc!
ALL: $(EXECUTABLE) chkopts
CFLAGS = ${PETSC_CC_INCLUDES}
CXXFLAGS = -std=c++11

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

$(EXECUTABLE) : $(OBJECTS)
	# Make bin directory if necessary
	mkdir -p bin/
	${CLINKER} ${CXXFLAGS} -o $@ $^ ${PETSC_LIB}
	${RM} $(SRCDIR)/*.o

