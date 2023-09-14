# ******************** MLCP MAKEFILE **************************

#-----------------------------------------------------------------------
#                      USER DEFINABLE OPTIONS
#-----------------------------------------------------------------------

# Compiler ( gfortran, ifort )
FC = gfortran

# Debugging options ( yes or no )
DEBUG = no

# Optimization level
OPTLEVEL = 3

# linking LAPACK and BLAS 
LAPACK = yes

# linking matrix ID package
MATID = yes

# Compile with standard real 8 (see details about the flags for each compiler...)
REAL8 = yes

# Compile for parallelization using OpenMP
OPENMP = yes

#-----------------------------------------------------------------------
#                        STRIP ALL SPACES
#-----------------------------------------------------------------------

# Strip leading and trailing spaces from all variables.
FC := $(strip ${FC})
DEBUG := $(strip ${DEBUG})
OPTLEVEL := $(strip ${OPTLEVEL})
LAPACK := $(strip ${LAPACK})
OPENMP := $(strip ${OPENMP})

#-----------------------------------------------------------------------
#                       Compiler specific statements.
#-----------------------------------------------------------------------

ifeq (${FC},gfortran)
 # Debug flags
   DEBUGFLG = -g -fcheck=all -fbacktrace
 # GNU LAPACK and BLAS flags
   LAPACKFLG = lib/liblapack.a lib/librefblas.a
 # matrix-ID flag
   MIDFLG = lib/id_lib.a
 # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
       DATAFLG = -fdefault-real-8 -fdefault-double-8
   endif
 # Flag to specify the position of mod files
   MODULEFLG = -I
 # Flag to specify OpenMP parallelization
   OMPFLG = -fopenmp
endif

ifeq (${FC},ifort)
 # Debug flags
   DEBUGFLG = -g -Wall -Wextra -pedantic -fimplicit-none -fcheck=all -fbacktrace
 # LAPACK and BLAS flags
   LAPACKFLG = lib/liblapack.a lib/librefblas.a
 # matrix-ID flag
   MIDFLG = lib/id_lib.a
 # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
       DATAFLG = -r8
   endif
 # Flag to specify the position of mod files
   MODULEFLG = -I
 # Flag to specify OpenMP parallelization
   OMPFLG = -openmp
endif

#-----------------------------------------------------------------------
#              Setup linking and compilation flags
#-----------------------------------------------------------------------

# initialize flags
COMPILEFLG =
LIBFLG     =

# if debugging set the appropriate flags
ifeq (${DEBUG}, yes)
    COMPILEFLG += ${DEBUGFLG}
endif

# Set flags for defining standard variable kinds
COMPILEFLG +=  ${DATAFLG}

# If matrixID, add the linking options
ifeq (${MATID}, yes)
    LIBFLG += ${MIDFLG}
endif

# If lapack, add the linking options
ifeq (${LAPACK}, yes)
    LIBFLG += ${LAPACKFLG}
endif

# If lapack, add the linking options
ifeq (${OPENMP}, yes)
    COMPILEFLG += ${OMPFLG}
endif

#-----------------------------------------------------------------------
#              Determine the optimization level to be used.
#-----------------------------------------------------------------------

ifeq (${OPTLEVEL},0)
    OPTFLAGS = -O0
endif
ifeq (${OPTLEVEL},1)
    OPTFLAGS = -O1
endif
ifeq (${OPTLEVEL},2)
    OPTFLAGS = -O2
endif
ifeq (${OPTLEVEL},3)
    OPTFLAGS = -O3
endif

# if debugging override input optimization level
ifeq (${DEBUG}, yes)
    OPTFLAGS = -O0
endif

COMPILEFLG += ${OPTFLAGS}

#-----------------------------------------------------------------------
#                         DIRECTORIES
#-----------------------------------------------------------------------

SRCDIR = Source
OBJDIR = Objects

#-----------------------------------------------------------------------
#                      List of object files
#-----------------------------------------------------------------------

# Define list of object from the list of all fortran files in the directory

# Common objects
COBJS = \
	${OBJDIR}/ErrorTrap.o \
	${OBJDIR}/Utils.o \
	${OBJDIR}/DSORTPLUSDEP.o \
	${OBJDIR}/ChebLib.o \
	${OBJDIR}/LinAlg.o \
	${OBJDIR}/Munkres.o \
	${OBJDIR}/TargetedStates.o \
	${OBJDIR}/InputFields.o \
	${OBJDIR}/ModeComb.o \
	${OBJDIR}/SepdRepn.o \
	${OBJDIR}/CPConfig.o \
	${OBJDIR}/FFPES.o \
	${OBJDIR}/MODVECVECML.o \
	${OBJDIR}/CPMM.o \
	${OBJDIR}/REDORTHO.o \
	${OBJDIR}/ALSOO.o \
	${OBJDIR}/ALS.o \
	${OBJDIR}/Reduction.o

# MLCP objects
MOBJS = \
	${OBJDIR}/OpFuncs.o \
	${OBJDIR}/HamilOpt.o \
	${OBJDIR}/HamilSetup.o \
	${OBJDIR}/ALSPow.o \
	${OBJDIR}/ALSUtils.o \
	${OBJDIR}/BlockUtils.o \
	${OBJDIR}/Restart.o \
	${OBJDIR}/Guess.o \
	${OBJDIR}/Updater.o \
	${OBJDIR}/Analyzer.o \
	${OBJDIR}/ModeH.o \
	${OBJDIR}/BlockPower.o \
	${OBJDIR}/LinSolver.o \
	${OBJDIR}/CPMath.o \
	${OBJDIR}/HGOrtho.o \
	${OBJDIR}/HG.o \
	${OBJDIR}/CPR.o \
	${OBJDIR}/TestCPR.o \
	${OBJDIR}/toy.o \
	${OBJDIR}/MSBII.o \
	${OBJDIR}/InvItn.o \
	${OBJDIR}/Solver.o \
	${OBJDIR}/MLmain.o

# CS-PES objects
CSOBJ = \
	${OBJDIR}/GenConfig.o \
	${OBJDIR}/CStools.o \
	${OBJDIR}/ModeGrid.o \
	${OBJDIR}/SAMPLEPES.o \
	${OBJDIR}/RestartCSPES.o \
	${OBJDIR}/BPsimpleCF.o \
	${OBJDIR}/BPsimpleCP.o \
	${OBJDIR}/CSPES.o

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
CSEXEFILE = cspes.x

# Make target to build all the object files and assemble them
all : ${MLEXEFILE} ${CSEXEFILE}

mlcp : ${MLEXEFILE}

cspes : ${CSEXEFILE}

${MLEXEFILE}: ${COBJS} ${MOBJS}
	${COMPILE} -o ${MLEXEFILE} ${COBJS} ${MOBJS} ${LIBFLG}
	mv *.mod ${OBJDIR}

${CSEXEFILE}: ${COBJS} ${CSOBJ}
	${COMPILE} -o ${CSEXEFILE} ${COBJS} ${CSOBJ} ${LIBFLG}
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
             ${OBJDIR}/Utils.o ${OBJDIR}/ChebLib.o Makefile

COMMONDEP2 = ${OBJDIR}/LinAlg.o ${OBJDIR}/Munkres.o \
	     ${OBJDIR}/InputFields.o ${OBJDIR}/ModeComb.o \
	     ${OBJDIR}/SepdRepn.o ${OBJDIR}/CPConfig.o \
	     ${OBJDIR}/MODVECVECML.o ${OBJDIR}/CPMM.o \
	     ${OBJDIR}/ALSOO.o ${OBJDIR}/FFPES.o \
	     ${OBJDIR}/REDORTHO.o ${OBJDIR}/Reduction.o \
	     ${COMMONDEP1}

# Sort vectors
${OBJDIR}/DSORTPLUSDEP.o : ${SRCDIR}/DSORTPLUSDEP.f Makefile

# Trap errors (a common dep.)
${OBJDIR}/ErrorTrap.o    : ${SRCDIR}/ErrorTrap.f90 Makefile

# Utilities
${OBJDIR}/Utils.o        : ${SRCDIR}/Utils.f90 ${OBJDIR}/ErrorTrap.o Makefile

# Chebyshev library
${OBJDIR}/ChebLib.o      : ${SRCDIR}/ChebLib.f90 ${OBJDIR}/ErrorTrap.o Makefile

# Linear algebra wrappers
${OBJDIR}/LinAlg.o       : ${SRCDIR}/LinAlg.f90 ${COMMONDEP1}

# Targeted states
${OBJDIR}/TargetedStates.o : ${SRCDIR}/TargetedStates.f90 ${COMMONDEP1}

# Hungarian algorithm matrix assignment
${OBJDIR}/Munkres.o      : ${SRCDIR}/Munkres.f90 ${COMMONDEP1}

# Input file reading
${OBJDIR}/InputFields.o  : ${SRCDIR}/InputFields.f90 ${COMMONDEP1}

# Mode combination module
${OBJDIR}/ModeComb.o     : ${SRCDIR}/ModeComb.f90 ${COMMONDEP1}

# CP-format types
${OBJDIR}/SepdRepn.o     : ${SRCDIR}/SepdRepn.f90 ${COMMONDEP1}

# CP configuration module
${OBJDIR}/CPConfig.o     : ${SRCDIR}/CPConfig.f90 ${OBJDIR}/SepdRepn.o ${COMMONDEP1}

# Force Field PES
${OBJDIR}/FFPES.o        : ${SRCDIR}/FFPES.f90 ${OBJDIR}/SepdRepn.o \
                           ${OBJDIR}/CPConfig.o ${COMMONDEP1}

# Separated representation linear algebra
${OBJDIR}/MODVECVECML.o  : ${SRCDIR}/MODVECVECML.f90 ${COMMONDEP1}

# Hamiltonian matrix-vector product
${OBJDIR}/CPMM.o         : ${SRCDIR}/CPMM.f90 ${OBJDIR}/SepdRepn.o ${COMMONDEP1}

# Orthogonal basis reduction
${OBJDIR}/REDORTHO.o     : ${SRCDIR}/REDORTHO.f90 ${OBJDIR}/CPConfig.o \
                           ${OBJDIR}/MODVECVECML.o ${COMMONDEP1}

# ALS reduction of psi in separated representation
${OBJDIR}/Reduction.o    : ${SRCDIR}/Reduction.f90 ${OBJDIR}/LinAlg.o \
                           ${OBJDIR}/MODVECVECML.o ${OBJDIR}/REDORTHO.o ${COMMONDEP1}

# Object-oriented ALS code
${OBJDIR}/ALSOO.o        : ${SRCDIR}/ALSOO.f90 ${OBJDIR}/LinAlg.o \
                           ${OBJDIR}/MODVECVECML.o ${COMMONDEP1}

# Driver for object-oriented ALS code 
${OBJDIR}/ALS.o          : ${SRCDIR}/ALS.f90 ${OBJDIR}/ALSOO.o \
	                   ${OBJDIR}/LinAlg.o ${OBJDIR}/MODVECVECML.o \
                           ${COMMONDEP1}

# Primitive operator functions
${OBJDIR}/OpFuncs.o      : ${SRCDIR}/OpFuncs.f90 ${COMMONDEP1}

# Hamiltonian optimization
${OBJDIR}/HamilOpt.o     : ${SRCDIR}/HamilOpt.f90 ${COMMONDEP2}

# Hamiltonian setup
${OBJDIR}/HamilSetup.o   : ${SRCDIR}/HamilSetup.f90 ${OBJDIR}/OpFuncs.o \
                           ${OBJDIR}/HamilOpt.o ${COMMONDEP2}

# Hamiltonian matrix-vector product + ALS
${OBJDIR}/ALSPow.o       : ${SRCDIR}/ALSPow.f90 ${OBJDIR}/CPMM.o \
                           ${COMMONDEP2}

# Hamiltonian matrix-vector product + ALS
${OBJDIR}/ALSUtils.o     : ${SRCDIR}/ALSUtils.f90 ${COMMONDEP2}

# Gram-Schmidt orthogonalization of separated representation vectors
${OBJDIR}/BlockUtils.o   : ${SRCDIR}/BlockUtils.f90 ${OBJDIR}/CPMM.o \
                           ${OBJDIR}/ALSUtils.o ${OBJDIR}/ALSPow.o \
                           ${COMMONDEP2}

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

# CP math operations
${OBJDIR}/CPMath.o       : ${SRCDIR}/CPMath.f90 ${OBJDIR}/LinSolver.o \
	                   ${OBJDIR}/CPMM.o ${OBJDIR}/ALS.o ${COMMONDEP2}

# CP orthogonalization
${OBJDIR}/HGOrtho.o      : ${SRCDIR}/HGOrtho.f90 ${OBJDIR}/CPMM.o \
                           ${OBJDIR}/CPMath.o ${OBJDIR}/ALS.o ${COMMONDEP2}

# CP diagonalization
${OBJDIR}/HG.o           : ${SRCDIR}/HG.f90 ${OBJDIR}/CPMath.o \
	                   ${OBJDIR}/LinSolver.o ${OBJDIR}/HGOrtho.o \
	                   ${OBJDIR}/CPMM.o \
	                   ${OBJDIR}/ALS.o ${COMMONDEP2}

# CP-in-rank format
${OBJDIR}/CPR.o          : ${SRCDIR}/CPR.f90 ${OBJDIR}/CPMath.o \
                           ${OBJDIR}/LinSolver.o ${OBJDIR}/HGOrtho.o \
                           ${OBJDIR}/CPMM.o \
                           ${OBJDIR}/ALS.o ${COMMONDEP2}

# CP-in-rank format
${OBJDIR}/TestCPR.o      : ${SRCDIR}/TestCPR.f90 ${OBJDIR}/CPMath.o \
                           ${OBJDIR}/LinSolver.o ${OBJDIR}/HGOrtho.o \
                           ${OBJDIR}/CPMM.o ${OBJDIR}/CPR.o \
                           ${OBJDIR}/ALS.o ${COMMONDEP2}

# CP diagonalization toy problem
${OBJDIR}/toy.o          : ${SRCDIR}/toy.f90 ${OBJDIR}/CPMath.o \
                           ${OBJDIR}/LinSolver.o ${OBJDIR}/HGOrtho.o \
                           ${OBJDIR}/CPMM.o ${OBJDIR}/ALS.o ${COMMONDEP2}

# MSBII solver
${OBJDIR}/MSBII.o        : ${SRCDIR}/MSBII.f90 ${OBJDIR}/LinSolver.o \
                           ${OBJDIR}/ALS.o ${OBJDIR}/BlockUtils.o ${COMMONDEP2}
${OBJDIR}/InvItn.o       : ${SRCDIR}/InvItn.f90 ${OBJDIR}/LinSolver.o \
                           ${OBJDIR}/ALS.o ${COMMONDEP2}

# Eigensolver
${OBJDIR}/Solver.o       : ${SRCDIR}/Solver.f90 ${OBJDIR}/BlockPower.o \
                           ${OBJDIR}/ALSPow.o ${OBJDIR}/Restart.o \
                           ${OBJDIR}/ALSUtils.o ${OBJDIR}/LinSolver.o \
			   ${OBJDIR}/CPMath.o ${OBJDIR}/HG.o \
			   ${OBJDIR}/toy.o ${OBJDIR}/MSBII.o \
			   ${OBJDIR}/InvItn.o ${COMMONDEP2} 

# Main MLCP program
${OBJDIR}/MLmain.o       : ${SRCDIR}/MLmain.f90 ${OBJDIR}/HamilSetup.o \
                           ${OBJDIR}/Restart.o ${OBJDIR}/ModeH.o \
                           ${OBJDIR}/Guess.o ${OBJDIR}/Solver.o \
                           ${OBJDIR}/Updater.o ${OBJDIR}/Analyzer.o \
                           ${OBJDIR}/ALSPow.o ${OBJDIR}/LinSolver.o \
			   ${OBJDIR}/TestCPR.o ${COMMONDEP2}

# Generating config lists
${OBJDIR}/GenConfig.o    : ${SRCDIR}/GenConfig.f90 ${COMMONDEP1}

# Tools for CS
${OBJDIR}/CStools.o      : ${SRCDIR}/CStools.f90 ${OBJDIR}/GenConfig.o \
                           ${COMMONDEP2}

# Reading grid input file
${OBJDIR}/ModeGrid.o     : ${SRCDIR}/ModeGrid.f90 ${COMMONDEP2}

# PES sampling
${OBJDIR}/SAMPLEPES.o    : ${SRCDIR}/SAMPLEPES.f90 ${COMMONDEP2}

# CSPES restarts
${OBJDIR}/RestartCSPES.o : ${SRCDIR}/RestartCSPES.f90 ${OBJDIR}/ModeGrid.o \
                           ${OBJDIR}/SAMPLEPES.o ${COMMONDEP2}

# Basis pursuit (CP)
${OBJDIR}/BPsimpleCP.o   : ${SRCDIR}/BPsimpleCP.f90 ${OBJDIR}/ModeGrid.o \
                           ${OBJDIR}/CStools.o ${OBJDIR}/RestartCSPES.o \
                           ${COMMONDEP2}

# Basis pursuit (configs)
${OBJDIR}/BPsimpleCF.o   : ${SRCDIR}/BPsimpleCF.f90 ${OBJDIR}/ModeGrid.o \
                           ${OBJDIR}/CStools.o ${OBJDIR}/RestartCSPES.o \
                           ${COMMONDEP2}

# Main CSPES program
${OBJDIR}/CSPES.o        : ${SRCDIR}/CSPES.f90 ${OBJDIR}/SAMPLEPES.o \
                           ${OBJDIR}/BPsimpleCP.o ${OBJDIR}/BPsimpleCF.o \
                           ${OBJDIR}/RestartCSPES.o ${COMMONDEP2}

