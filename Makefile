# ******************** MLCP MAKEFILE **************************
$(if $(wildcard arch.mk),,$(error Error: Please create arch.mk from config/ directory for machine-dependent configuration))
include arch.mk

#-----------------------------------------------------------------------
#                        STRIP ALL SPACES
#-----------------------------------------------------------------------

# Strip leading and trailing spaces from all variables.
FC := $(strip ${FC})
FOPTS := $(strip ${FOPTS})
OPTFLG := $(strip ${OPTFLG})
DEBUG := $(strip ${DEBUG})
DEBUGFLG := $(strip ${DEBUGFLG})
MODULEFLG := $(strip ${MODULEFLG})
MPIFLG := $(strip ${MPIFLG})
USEOMP := $(strip ${USEOMP})
OMPFLG := $(strip ${OMPFLG})
USEACC := $(strip ${USEACC})
ACCFLG := $(strip ${ACCFLG})
PREPROCFLG := $(strip ${PREPROCFLG})

LAPACKLIB := $(strip ${LAPACKLIB})
CULIB := $(strip ${CULIB})

#-----------------------------------------------------------------------
#              Setup linking and compilation flags
#-----------------------------------------------------------------------

# Compiler flags and libraries
COMPILEFLG =
COMPILEFLG += ${FOPTS} 

# If debugging set the appropriate flags
ifeq (${DEBUG}, yes)
    COMPILEFLG += ${DEBUGFLG}
else
    COMPILEFLG += ${OPTFLG}
endif

COMPILEFLG += ${MPIFLG}

ifeq (${USEOMP}, yes)
    COMPILEFLG += ${OMPFLG}
    PREPROCFLG += -DOMP_ENABLED=1
else
    PREPROCFLG += -DOMP_ENABLED=0
endif

# Build OpenACC code for GPUs
ifeq (${USEACC},yes)
    COMPILEFLG += ${ACCFLG}
    PREPROCFLG += -DACC_ENABLED=1
    LIBFLG += ${CULIB}
else
    PREPROCFLG += -DACC_ENABLED=0
endif

COMPILEFLG += ${PREPROCFLG}

LIBFLG += ${LAPACKLIB}

#-----------------------------------------------------------------------
#                         DIRECTORIES
#-----------------------------------------------------------------------

SRCDIR = src
OBJDIR = obj

#-----------------------------------------------------------------------
#                      List of object files
#-----------------------------------------------------------------------

# Define list of object from the list of all fortran files in the directory

# Common objects
COBJS = \
	${OBJDIR}/ErrorTrap.o \
	${OBJDIR}/Utils.o \
	${OBJDIR}/MyMPI.o \
	${OBJDIR}/MyACC.o \
	${OBJDIR}/Random.o \
	${OBJDIR}/DSORTPLUSDEP.o \
	${OBJDIR}/ChebLib.o \
	${OBJDIR}/LinAlg.o \
	${OBJDIR}/LinAlg8.o \
	${OBJDIR}/Munkres.o \
	${OBJDIR}/TargetedStates.o \
	${OBJDIR}/InputCP.o \
	${OBJDIR}/InputCS.o \
	${OBJDIR}/NodeTree.o \
	${OBJDIR}/ModeComb.o \
	${OBJDIR}/SepdRepn.o \
	${OBJDIR}/CPr8.o \
	${OBJDIR}/CPVV8.o \
	${OBJDIR}/CPMM8.o \
	${OBJDIR}/CPLS8.o \
	${OBJDIR}/TestCPr8.o \
	${OBJDIR}/CPConfig.o \
	${OBJDIR}/FFPES.o \
	${OBJDIR}/MODVECVECML.o \
	${OBJDIR}/CPMM.o \
	${OBJDIR}/REDORTHO.o \
	${OBJDIR}/ALSOO.o \
	${OBJDIR}/ALSOO8.o \
	${OBJDIR}/ALS.o \
	${OBJDIR}/ALS8.o \
	${OBJDIR}/Reduction.o

# MLCP objects
MOBJS = \
	${OBJDIR}/OpFuncs.o \
	${OBJDIR}/HamilSetup.o \
	${OBJDIR}/ALSPow.o \
	${OBJDIR}/ALSUtils.o \
	${OBJDIR}/BlockUtils.o \
	${OBJDIR}/FEAST8.o \
	${OBJDIR}/Restart.o \
	${OBJDIR}/Guess.o \
	${OBJDIR}/Updater.o \
	${OBJDIR}/Analyzer.o \
	${OBJDIR}/ModeH.o \
	${OBJDIR}/BlockPower.o \
	${OBJDIR}/LinSolver.o \
	${OBJDIR}/Solver.o \
	${OBJDIR}/Solver_CP8.o \
	${OBJDIR}/Timings.o \
	${OBJDIR}/MLmain.o

#-----------------------------------------------------------------------
#       Construct the compile and link variables
#-----------------------------------------------------------------------

# Compile command: ${COMPILE} <source>
COMPILE                 = ${FC} ${COMPILEFLG} ${MODULEFLG} ${OBJDIR}

#-----------------------------------------------------------------------
#                         MAKE RULES
#-----------------------------------------------------------------------

.SUFFIXES: .f90 .o .x

MLEXEFILE = mlcp.x

# Make target to build all the object files and assemble them
all : ${MLEXEFILE} 

mlcp : ${MLEXEFILE}

${MLEXEFILE}: ${COBJS} ${MOBJS}
	${COMPILE} -o ${MLEXEFILE} ${COBJS} ${MOBJS} ${LIBFLG}
	mv *.mod ${OBJDIR}

# Make a target object file by compiling the fortran code
${OBJDIR}/%.o : ${SRCDIR}/%.f90
	${COMPILE} -c ${SRCDIR}/$*.f90
	mv *.o ${OBJDIR}
${OBJDIR}/%.o : ${SRCDIR}/%.f
	${COMPILE} -c ${SRCDIR}/$*.f
	mv *.o ${OBJDIR}

# Make target to build required directories
directories : ${OBJDIR}
	mkdir -p ${OBJDIR}

# Remove compiled objects and related stuff
clean :
	rm -rf ${OBJDIR}/*

# ----------------------------------------------------------------------
# ----------------------- DEPENDENCIES----------------------------------
# ----------------------------------------------------------------------

COMMONDEP1 = ${OBJDIR}/DSORTPLUSDEP.o ${OBJDIR}/ErrorTrap.o \
             ${OBJDIR}/Utils.o ${OBJDIR}/MyMPI.o ${OBJDIR}/MyACC.o \
             ${OBJDIR}/ChebLib.o Makefile

COMMONDEP2 = ${OBJDIR}/Random.o ${OBJDIR}/LinAlg.o ${OBJDIR}/Munkres.o \
	     ${OBJDIR}/InputCP.o ${OBJDIR}/InputCS.o ${OBJDIR}/ModeComb.o \
	     ${OBJDIR}/SepdRepn.o ${OBJDIR}/CPConfig.o ${OBJDIR}/LinAlg8.o \
	     ${OBJDIR}/CPr8.o ${OBJDIR}/MODVECVECML.o ${OBJDIR}/CPMM.o \
	     ${OBJDIR}/ALSOO.o ${OBJDIR}/FFPES.o ${OBJDIR}/NodeTree.o \
	     ${OBJDIR}/REDORTHO.o ${OBJDIR}/Reduction.o \
	     ${COMMONDEP1}

# Sort vectors
${OBJDIR}/DSORTPLUSDEP.o : ${SRCDIR}/DSORTPLUSDEP.f Makefile

# Trap errors (a common dep.)
${OBJDIR}/ErrorTrap.o    : ${SRCDIR}/ErrorTrap.f90 Makefile

# Utilities
${OBJDIR}/Utils.o        : ${SRCDIR}/Utils.f90 ${OBJDIR}/ErrorTrap.o Makefile

# MPI wrapper functions
${OBJDIR}/MyMPI.o        : ${SRCDIR}/MyMPI.f90 ${OBJDIR}/ErrorTrap.o Makefile

# MPI wrapper functions
${OBJDIR}/MyACC.o        : ${SRCDIR}/MyACC.f90 ${OBJDIR}/ErrorTrap.o \
	                   ${OBJDIR}/MyMPI.o Makefile

# Chebyshev library
${OBJDIR}/ChebLib.o      : ${SRCDIR}/ChebLib.f90 ${OBJDIR}/ErrorTrap.o Makefile

# Random initialization
${OBJDIR}/Random.o       : ${SRCDIR}/Random.f90 ${COMMONDEP1}

# Linear algebra wrappers
${OBJDIR}/LinAlg.o       : ${SRCDIR}/LinAlg.f90 ${COMMONDEP1}

# Linear algebra wrappers
${OBJDIR}/LinAlg8.o       : ${SRCDIR}/LinAlg8.f90 ${COMMONDEP1}

# Targeted states
${OBJDIR}/TargetedStates.o : ${SRCDIR}/TargetedStates.f90 ${COMMONDEP1}

# Hungarian algorithm matrix assignment
${OBJDIR}/Munkres.o      : ${SRCDIR}/Munkres.f90 ${COMMONDEP1}

# CP.inp input reading
${OBJDIR}/InputCP.o  : ${SRCDIR}/InputCP.f90 ${COMMONDEP1}

# CS.inp input reading
${OBJDIR}/InputCS.o  : ${SRCDIR}/InputCS.f90 ${COMMONDEP1}

# Node tree module
${OBJDIR}/NodeTree.o     : ${SRCDIR}/NodeTree.f90 ${COMMONDEP1}

# Mode combination module
${OBJDIR}/ModeComb.o     : ${SRCDIR}/ModeComb.f90 ${OBJDIR}/NodeTree.o ${COMMONDEP1}

# CP-format types
${OBJDIR}/SepdRepn.o     : ${SRCDIR}/SepdRepn.f90 ${COMMONDEP1}

# CP-format types
${OBJDIR}/CPr8.o         : ${SRCDIR}/CPr8.f90 ${OBJDIR}/SepdRepn.o ${COMMONDEP1}

# CP-format vector inner products
${OBJDIR}/CPVV8.o        : ${SRCDIR}/CPVV8.f90 ${OBJDIR}/CPr8.o ${COMMONDEP1}

# CP-format matrix multiply
${OBJDIR}/CPMM8.o        : ${SRCDIR}/CPMM8.f90 ${OBJDIR}/CPr8.o ${COMMONDEP1}

# CP-format linear systems
${OBJDIR}/CPLS8.o        : ${SRCDIR}/CPLS8.f90 ${OBJDIR}/CPr8.o ${COMMONDEP1}

# CP configuration module
${OBJDIR}/CPConfig.o     : ${SRCDIR}/CPConfig.f90 ${OBJDIR}/SepdRepn.o ${COMMONDEP1}

# Force Field PES
${OBJDIR}/FFPES.o        : ${SRCDIR}/FFPES.f90 ${OBJDIR}/SepdRepn.o \
                           ${OBJDIR}/CPConfig.o ${COMMONDEP1}

# Separated representation linear algebra
${OBJDIR}/MODVECVECML.o  : ${SRCDIR}/MODVECVECML.f90 ${COMMONDEP1}

# Hamiltonian matrix-vector product
${OBJDIR}/CPMM.o         : ${SRCDIR}/CPMM.f90 ${OBJDIR}/SepdRepn.o \
	                   ${OBJDIR}/MODVECVECML.o ${COMMONDEP1}

# Orthogonal basis reduction
${OBJDIR}/REDORTHO.o     : ${SRCDIR}/REDORTHO.f90 ${OBJDIR}/CPConfig.o \
                           ${OBJDIR}/MODVECVECML.o ${COMMONDEP1}

# ALS reduction of psi in separated representation
${OBJDIR}/Reduction.o    : ${SRCDIR}/Reduction.f90 ${OBJDIR}/LinAlg.o \
                           ${OBJDIR}/MODVECVECML.o ${OBJDIR}/REDORTHO.o ${COMMONDEP1}

# Object-oriented ALS code
${OBJDIR}/ALSOO.o        : ${SRCDIR}/ALSOO.f90 ${OBJDIR}/LinAlg.o \
                           ${OBJDIR}/MODVECVECML.o ${COMMONDEP1}

# Object-oriented ALS code
${OBJDIR}/ALSOO8.o       : ${SRCDIR}/ALSOO8.f90 ${OBJDIR}/CPr8.o ${OBJDIR}/LinAlg8.o \
                           ${OBJDIR}/CPVV8.o ${OBJDIR}/CPMM8.o ${OBJDIR}/CPLS8.o \
			   ${COMMONDEP1}

# Driver for object-oriented ALS code 
${OBJDIR}/ALS.o          : ${SRCDIR}/ALS.f90 ${OBJDIR}/ALSOO.o \
	                   ${OBJDIR}/LinAlg.o ${OBJDIR}/MODVECVECML.o \
                           ${COMMONDEP1}


# Driver for object-oriented ALS code
${OBJDIR}/ALS8.o         : ${SRCDIR}/ALS8.f90 ${OBJDIR}/ALSOO8.o ${OBJDIR}/CPr8.o \
                           ${OBJDIR}/CPVV8.o ${OBJDIR}/CPMM8.o ${OBJDIR}/CPLS8.o \
                           ${COMMONDEP1}

# Primitive operator functions
${OBJDIR}/OpFuncs.o      : ${SRCDIR}/OpFuncs.f90 ${COMMONDEP1}

# Hamiltonian setup
${OBJDIR}/HamilSetup.o   : ${SRCDIR}/HamilSetup.f90 ${OBJDIR}/OpFuncs.o \
                           ${COMMONDEP2}

# Hamiltonian matrix-vector product + ALS
${OBJDIR}/ALSPow.o       : ${SRCDIR}/ALSPow.f90 ${OBJDIR}/CPMM.o \
                           ${COMMONDEP2}

# Hamiltonian matrix-vector product + ALS
${OBJDIR}/ALSUtils.o     : ${SRCDIR}/ALSUtils.f90 ${COMMONDEP2}

# Gram-Schmidt orthogonalization of separated representation vectors
${OBJDIR}/BlockUtils.o   : ${SRCDIR}/BlockUtils.f90 ${OBJDIR}/CPMM.o \
                           ${OBJDIR}/ALSUtils.o ${OBJDIR}/ALSPow.o \
                           ${COMMONDEP2}

${OBJDIR}/FEAST8.o       : ${SRCDIR}/FEAST8.f90 ${COMMONDEP2}

# Restart a crashed calculation
${OBJDIR}/Restart.o      : ${SRCDIR}/Restart.f90 ${COMMONDEP2}

# Wavefunction initial guess
${OBJDIR}/Guess.o        : ${SRCDIR}/Guess.f90 ${OBJDIR}/HamilSetup.o \
                           ${OBJDIR}/TargetedStates.o ${COMMONDEP2}
                                                       
# Hamiltonian updates
${OBJDIR}/Updater.o      : ${SRCDIR}/Updater.f90 ${OBJDIR}/HamilSetup.o \
                           ${COMMONDEP2}

# Analysis
${OBJDIR}/Analyzer.o     : ${SRCDIR}/Analyzer.f90 ${COMMONDEP2}

# Block and Hamiltonian initialization
${OBJDIR}/ModeH.o        : ${SRCDIR}/ModeH.f90 ${OBJDIR}/HamilSetup.o \
                           ${COMMONDEP2}

# Block Power code
${OBJDIR}/BlockPower.o   : ${SRCDIR}/BlockPower.f90 ${OBJDIR}/CPMM.o \
                           ${OBJDIR}/LinSolver.o ${OBJDIR}/BlockUtils.o ${COMMONDEP2}

# Linear equation solver
${OBJDIR}/LinSolver.o    : ${SRCDIR}/LinSolver.f90 ${OBJDIR}/ALSPow.o \
	                   ${OBJDIR}/CPMM.o ${COMMONDEP2}

# Eigensolver CP8
${OBJDIR}/Solver_CP8.o   : ${SRCDIR}/Solver_CP8.f90 ${OBJDIR}/CPr8.o ${OBJDIR}/ALS8.o \
                           ${OBJDIR}/Restart.o ${OBJDIR}/BlockUtils.o \
                           ${OBJDIR}/FEAST8.o ${COMMONDEP2}

# Eigensolver
${OBJDIR}/Solver.o       : ${SRCDIR}/Solver.f90 ${OBJDIR}/BlockPower.o \
                           ${OBJDIR}/ALSPow.o ${OBJDIR}/Restart.o \
                           ${OBJDIR}/ALSUtils.o ${OBJDIR}/LinSolver.o \
			   ${OBJDIR}/Solver_CP8.o ${COMMONDEP2}

# Test CP-format types
${OBJDIR}/TestCPr8.o     : ${SRCDIR}/TestCPr8.f90 ${OBJDIR}/CPr8.o ${OBJDIR}/CPMM.o \
                           ${OBJDIR}/CPVV8.o ${OBJDIR}/CPMM8.o ${OBJDIR}/CPLS8.o \
                           ${OBJDIR}/ALSOO.o ${OBJDIR}/LinAlg8.o ${OBJDIR}/SepdRepn.o \
                           ${OBJDIR}/ALS.o ${OBJDIR}/ALSOO8.o ${OBJDIR}/ALS8.o \
			   ${OBJDIR}/LinSolver.o ${COMMONDEP1}

# Timings
${OBJDIR}/Timings.o      : ${SRCDIR}/Timings.f90 ${OBJDIR}/HamilSetup.o \
                           ${OBJDIR}/Restart.o ${OBJDIR}/ModeH.o \
                           ${OBJDIR}/Guess.o ${OBJDIR}/Solver.o \
                           ${OBJDIR}/Updater.o ${OBJDIR}/Analyzer.o \
                           ${OBJDIR}/ALSPow.o ${OBJDIR}/LinSolver.o \
                           ${COMMONDEP2}

# Main MLCP program
${OBJDIR}/MLmain.o       : ${SRCDIR}/MLmain.f90 ${OBJDIR}/HamilSetup.o \
                           ${OBJDIR}/Restart.o ${OBJDIR}/ModeH.o \
                           ${OBJDIR}/Guess.o ${OBJDIR}/Solver.o \
                           ${OBJDIR}/Updater.o ${OBJDIR}/Analyzer.o \
                           ${OBJDIR}/ALSPow.o ${OBJDIR}/LinSolver.o \
			   ${OBJDIR}/TestCPr8.o ${OBJDIR}/Timings.o ${COMMONDEP2}

