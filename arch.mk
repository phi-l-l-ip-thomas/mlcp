# Fortran compiler
FC = mpifort

# Flags to always include
FOPTS = -Mr8

OPTFLG = -O3

# Debug flags
DEBUG = yes
DEBUGFLG = -O0 -traceback

# Flag to specify the position of mod files
MODULEFLG = -I

# Flag to specify Message Passing Interface parallelization
MPIFLG = 

# Flag to specify OpenMP parallelization
OMPFLG = -mp

# Preprocessor flag
PREPROCFLG = -Mpreprocess

# LAPACK and BLAS flags
#LAPACKLIB = lib/liblapack.a lib/librefblas.a
LAPACKLIB = -L/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/compilers/lib -lblas -llapack


