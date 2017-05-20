ALL: main
CFLAGS = ${PETSC_CC_INCLUDES}
CLEANFILES = main
CXXFLAGS = -std=c++11

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: main.o chkopts
	${CLINKER} ${CXXFLAGS} -o $@ main.o ${PETSC_LIB}
	${RM} main.o
