# ******************** MLCP MAKEFILE **************************

#-----------------------------------------------------------------------
#                      USER DEFINABLE OPTIONS
#-----------------------------------------------------------------------

# Compilers ( gfortran, ifort )
FC = gfortran

# CUDA compiler
CUC = nvcc

# Using CUDA
USECUDA = yes

# Double precision in CUDA
USEDOUBLES = yes

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
CUC := $(strip ${CUC})
USECUDA := $(strip ${USECUDA})
USEDOUBLES := $(strip ${USEDOUBLES})
DEBUG := $(strip ${DEBUG})
OPTLEVEL := $(strip ${OPTLEVEL})
LAPACK := $(strip ${LAPACK})
MINPACK := $(strip ${MINPACK})
OPENMP := $(strip ${OPENMP})

#-----------------------------------------------------------------------
#                       Compiler specific statements.
#-----------------------------------------------------------------------

ifeq (${FC},gfortran)
 # Debug flags
   DEBUGFLG = -g -fbounds-check
 # GNU Lapack and Blas flags
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
   DEBUGFLG = -g -check bounds
 # Lapack and Blas flags
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


ifeq (${CUC},nvcc)
 # Debug flags
   CUDEBUGFLG = -g -G -lineinfo
 # Flag to specify the position of mod files
   CUMODULEFLG = -I
 # Flag for CUDA libraries
#   CULIBFLG = -lstdc++ -lgcc -L/usr/local/cuda/lib64/ -lcuda -lcudart
   CULIBFLG = -lstdc++ -lgcc -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/cuda/lib64 -lcudart -lcuda
 # Flag to specify OpenMP parallelization
   CUOMPFLG = -Xcompiler -fopenmp
 # CUDA precision flag
   ifeq (${USEDOUBLES},yes)
      CUPRECFLG = -DUSE_DOUBLES=1
   else
      CUPRECFLG = -DUSE_DOUBLES=0
   endif
endif

#-----------------------------------------------------------------------
#              Setup linking and compilation flags
#-----------------------------------------------------------------------

# initialize flags
COMPILEFLG   = -cpp
CUCOMPILEFLG = -Wno-deprecated-gpu-targets -default-stream per-thread --ptxas-options=-v
LIBFLG       = 

# if debugging set the appropriate flags
ifeq (${DEBUG}, yes)
    COMPILEFLG += ${DEBUGFLG}
    CUCOMPILEFLG += ${CUDEBUGFLG}
endif

# Set flags for defining standard variable kinds
COMPILEFLG +=  ${DATAFLG}

# If using CUDA, add the linking options
ifeq (${USECUDA}, yes)
    LIBFLG += ${CULIBFLG}
endif

# If matrixID, add the linking options
ifeq (${MATID}, yes)
    LIBFLG += ${MIDFLG}
endif

# If lapack, add the linking options
ifeq (${LAPACK}, yes)
    LIBFLG += ${LAPACKFLG}
   ifeq (${USECUDA}, yes)
#       LIBFLG += -lcublas -lcusolver
       LIBFLG += -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/math_libs/lib64 -lcublas -lcusolver
   endif
endif

# If openmp, add the linking options
ifeq (${OPENMP}, yes)
    COMPILEFLG += ${OMPFLG}
    CUCOMPILEFLG += ${CUOMPFLG}
endif

ifeq (${USEDOUBLES},yes)
   COMPILEFLG += ${CUPRECFLG}
   CUCOMPILEFLG += ${CUPRECFLG}
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
CUCOMPILEFLG += ${OPTFLAGS}

#-----------------------------------------------------------------------
#                         DIRECTORIES
#-----------------------------------------------------------------------

SRCDIR = Source
OBJDIR = Objects

#-----------------------------------------------------------------------
#                      List of object files
#-----------------------------------------------------------------------

# Define list of object from the list of all fortran files in the directory

# CUDA objects
CUOBJS = \
	${OBJDIR}/utils.o \
	${OBJDIR}/matrix.o \
	${OBJDIR}/linalg.o \
	${OBJDIR}/Munkres.o \
	${OBJDIR}/cpvector.o \
	${OBJDIR}/hoplist.o \
	${OBJDIR}/modvecvec.o \
	${OBJDIR}/als.o \
	${OBJDIR}/intertw_mvp.o \
	${OBJDIR}/intertw_lcv.o \
	${OBJDIR}/intertwining.o

# CUDA shared objects
CUSOBJS = \
	${OBJDIR}/culink.obj 

# Common objects
COBJS = \
	${OBJDIR}/ErrorTrap.o \
	${OBJDIR}/Utils.o \
	${OBJDIR}/DSORTPLUSDEP.o \
	${OBJDIR}/ChebLib.o \
	${OBJDIR}/LinAlg.o \
	${OBJDIR}/GenConfig.o \
	${OBJDIR}/InputFields.o \
	${OBJDIR}/ModeComb.o \
	${OBJDIR}/SepdRepn.o \
	${OBJDIR}/CPConfig.o \
	${OBJDIR}/FFPES.o \
	${OBJDIR}/MODVECVECML.o \
	${OBJDIR}/REDORTHO.o \
	${OBJDIR}/REDUCTIONML.o

# MLCP objects
MOBJS = \
	${OBJDIR}/OpFuncs.o \
	${OBJDIR}/HamilOpt.o \
	${OBJDIR}/HamilSetup.o \
	${OBJDIR}/MODHVECML.o \
	${OBJDIR}/ALSPow.o \
	${OBJDIR}/ALSUtils.o \
	${OBJDIR}/GPUIterate.o \
	${OBJDIR}/BlockUtils.o \
	${OBJDIR}/Restart.o \
	${OBJDIR}/Guess.o \
	${OBJDIR}/Updater.o \
	${OBJDIR}/Analyzer.o \
	${OBJDIR}/ModeH.o \
	${OBJDIR}/BlockPower.o \
	${OBJDIR}/Cheb.o \
	${OBJDIR}/Solver.o \
	${OBJDIR}/MLmain.o

# CS-PES objects
CSOBJ = \
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

CUCOMPILE               = ${CUC} ${CUCOMPILEFLG} ${CUMODULEFLG} ${OBJDIR}

#-----------------------------------------------------------------------
#                         MAKE RULES
#-----------------------------------------------------------------------

.SUFFIXES: .f90 .o .x .obj .cu

MLEXEFILE = mlcp.x
CSEXEFILE = cspes.x

# Make target to build all the object files and assemble them
all : ${MLEXEFILE} ${CSEXEFILE}

mlcp : ${MLEXEFILE}

cspes : ${CSEXEFILE}

${MLEXEFILE}: ${CUOBJS} ${CUSOBJS} ${COBJS} ${MOBJS}
	${COMPILE} -o ${MLEXEFILE} ${CUOBJS} ${CUSOBJS} ${COBJS} ${MOBJS} ${LIBFLG}
	mv *.mod ${OBJDIR}

${CSEXEFILE}: ${COBJS} ${CSOBJ}
	${COMPILE} -o ${CSEXEFILE} ${COBJS} ${CSOBJ} ${LIBFLG}
	mv *.mod ${OBJDIR}

# Make a target object file by compiling the source code
${OBJDIR}/%.o : ${SRCDIR}/%.f90
	${COMPILE} -c ${SRCDIR}/$*.f90
	mv *.o ${OBJDIR}
${OBJDIR}/%.o : ${SRCDIR}/%.f
	${COMPILE} -c ${SRCDIR}/$*.f
	mv *.o ${OBJDIR}
${OBJDIR}/%.o : ${SRCDIR}/%.cu
	${CUCOMPILE} -dc ${SRCDIR}/$*.cu
	mv *.o ${OBJDIR}
${OBJDIR}/%.obj :  
	${CUCOMPILE} -dlink ${CUOBJS} -o $*.obj
	mv *.obj ${OBJDIR}
 		
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
             ${OBJDIR}/MODVECVECML.o ${OBJDIR}/FFPES.o \
             ${OBJDIR}/REDORTHO.o ${OBJDIR}/REDUCTIONML.o \
             ${COMMONDEP1}

# CUDA utiltiies
${OBJDIR}/utils.o        : ${SRCDIR}/utils.cu Makefile

# CUDA matrix functions
${OBJDIR}/matrix.o       : ${SRCDIR}/matrix.cu Makefile

# CUDA linear algebra wrappers
${OBJDIR}/linalg.o       : ${SRCDIR}/linalg.cu ${OBJDIR}/utils.o \
                           ${OBJDIR}/matrix.o Makefile

# Hungarian algorithm matrix assignment
${OBJDIR}/Munkres.o      : ${SRCDIR}/Munkres.f90 ${COMMONDEP1}

# CUDA CP functions
${OBJDIR}/cpvector.o     : ${SRCDIR}/cpvector.cu Makefile 

# CUDA Hamiltonian storage 
${OBJDIR}/hoplist.o      : ${SRCDIR}/hoplist.cu Makefile

# CUDA vector inner products and normalization
${OBJDIR}/modvecvec.o    : ${SRCDIR}/modvecvec.cu ${OBJDIR}/utils.o \
                           ${OBJDIR}/matrix.o ${OBJDIR}/cpvector.o \
                           Makefile

# CUDA ALS helper kernel
${OBJDIR}/als.o          : ${SRCDIR}/als.cu ${OBJDIR}/utils.o \
                           ${OBJDIR}/matrix.o ${OBJDIR}/cpvector.o \
                           ${OBJDIR}/modvecvec.o Makefile

# CUDA MVP+ALS intertwining
${OBJDIR}/intertw_mvp.o  : ${SRCDIR}/intertw_mvp.cu ${OBJDIR}/hoplist.o \
                           ${OBJDIR}/matrix.o ${OBJDIR}/cpvector.o \
                           ${OBJDIR}/modvecvec.o ${OBJDIR}/als.o \
                           ${OBJDIR}/utils.o Makefile

# CUDA Gram-Schmidt/Vector-Updates+ALS intertwining
${OBJDIR}/intertw_lcv.o  : ${SRCDIR}/intertw_lcv.cu ${OBJDIR}/matrix.o \
                           ${OBJDIR}/cpvector.o ${OBJDIR}/modvecvec.o \
                           ${OBJDIR}/als.o ${OBJDIR}/utils.o Makefile

# CUDA intertwining driver
${OBJDIR}/intertwining.o : ${SRCDIR}/intertwining.cu ${OBJDIR}/hoplist.o \
                           ${OBJDIR}/matrix.o ${OBJDIR}/cpvector.o \
                           ${OBJDIR}/intertw_mvp.o ${OBJDIR}/intertw_lcv.o \
                           ${OBJDIR}/utils.o ${OBJDIR}/als.o Makefile

# CUDA linking
${OBJDIR}/culink.obj     : ${OBJDIR}/utils.o ${OBJDIR}/matrix.o \
                           ${OBJDIR}/cpvector.o ${OBJDIR}/hoplist.o \
                           ${OBJDIR}/modvecvec.o ${OBJDIR}/als.o \
                           ${OBJDIR}/intertw_mvp.o ${OBJDIR}/intertw_lcv.o \
                           ${OBJDIR}/intertwining.o Makefile

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

# Orthogonal basis reduction
${OBJDIR}/REDORTHO.o     : ${SRCDIR}/REDORTHO.f90 ${OBJDIR}/CPConfig.o \
                           ${OBJDIR}/MODVECVECML.o ${COMMONDEP1}

# ALS reduction of psi in separated representation
${OBJDIR}/REDUCTIONML.o  : ${SRCDIR}/REDUCTIONML.f90 ${OBJDIR}/LinAlg.o \
                           ${OBJDIR}/MODVECVECML.o ${OBJDIR}/REDORTHO.o ${COMMONDEP1}

# Primitive operator functions
${OBJDIR}/OpFuncs.o      : ${SRCDIR}/OpFuncs.f90 ${COMMONDEP1}

# Hamiltonian optimization
${OBJDIR}/HamilOpt.o     : ${SRCDIR}/HamilOpt.f90 ${COMMONDEP2}

# Hamiltonian setup
${OBJDIR}/HamilSetup.o   : ${SRCDIR}/HamilSetup.f90 ${OBJDIR}/OpFuncs.o \
                           ${OBJDIR}/HamilOpt.o ${COMMONDEP2}

# Hamiltonian matrix-vector product
${OBJDIR}/MODHVECML.o    : ${SRCDIR}/MODHVECML.f90 ${COMMONDEP2}

# Hamiltonian matrix-vector product + ALS
${OBJDIR}/ALSPow.o       : ${SRCDIR}/ALSPow.f90 ${OBJDIR}/MODHVECML.o \
                           ${COMMONDEP2}

# Hamiltonian matrix-vector product + ALS
${OBJDIR}/ALSUtils.o     : ${SRCDIR}/ALSUtils.f90 ${COMMONDEP2}

# CUDA GPU code wrapper
${OBJDIR}/GPUIterate.o   : ${SRCDIR}/GPUIterate.f90 ${OBJDIR}/HamilSetup.o \
                           ${OBJDIR}/intertwining.o ${COMMONDEP2}

# Gram-Schmidt orthogonalization of separated representation vectors
${OBJDIR}/BlockUtils.o   : ${SRCDIR}/BlockUtils.f90 ${OBJDIR}/MODHVECML.o \
                           ${OBJDIR}/ALSUtils.o ${OBJDIR}/ALSPow.o \
                           ${COMMONDEP2}

# Restart a crashed calculation
${OBJDIR}/Restart.o      : ${SRCDIR}/Restart.f90 ${COMMONDEP2}

# Wavefunction initial guess
${OBJDIR}/Guess.o        : ${SRCDIR}/Guess.f90 ${OBJDIR}/HamilSetup.o \
                           ${COMMONDEP2}
                                                       
# Hamiltonian updates
${OBJDIR}/Updater.o      : ${SRCDIR}/Updater.f90 ${OBJDIR}/HamilSetup.o \
                           ${COMMONDEP2}

# Analysis
${OBJDIR}/Analyzer.o     : ${SRCDIR}/Analyzer.f90 ${COMMONDEP2}

# Block and Hamiltonian initialization
${OBJDIR}/ModeH.o        : ${SRCDIR}/ModeH.f90 ${OBJDIR}/HamilSetup.o \
                           ${COMMONDEP2}

# Block Power code
${OBJDIR}/BlockPower.o   : ${SRCDIR}/BlockPower.f90 ${OBJDIR}/MODHVECML.o \
                           ${OBJDIR}/BlockUtils.o ${COMMONDEP2}

# Chebyshev solver
${OBJDIR}/Cheb.o         : ${SRCDIR}/Cheb.f90 ${OBJDIR}/MODHVECML.o \
                           ${OBJDIR}/BlockUtils.o ${COMMONDEP2}

# Eigensolver
${OBJDIR}/Solver.o       : ${SRCDIR}/Solver.f90 ${OBJDIR}/BlockPower.o \
                           ${OBJDIR}/Cheb.o ${OBJDIR}/GPUIterate.o \
                           ${OBJDIR}/ALSPow.o ${OBJDIR}/Restart.o \
                           ${OBJDIR}/ALSUtils.o ${COMMONDEP2} 

# Main MLCP program
${OBJDIR}/MLmain.o       : ${SRCDIR}/MLmain.f90 ${OBJDIR}/HamilSetup.o \
                           ${OBJDIR}/Restart.o ${OBJDIR}/ModeH.o \
                           ${OBJDIR}/Guess.o ${OBJDIR}/Solver.o \
                           ${OBJDIR}/Updater.o ${OBJDIR}/Analyzer.o \
                           ${OBJDIR}/ALSPow.o ${OBJDIR}/GPUIterate.o \
                           ${COMMONDEP2} \

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

