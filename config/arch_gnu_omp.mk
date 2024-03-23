# Fortran compiler
FC = gfortran

# Flags to always include
FOPTS = -fdefault-real-8 -fdefault-double-8

OPTFLG = -O3

# Debug flags
DEBUG = no
DEBUGFLG = -O0 -g -fcheck=all -fbacktrace

# Flag to specify the position of mod files
MODULEFLG = -I

# Flag to specify Message Passing Interface parallelization
MPIFLG = 

# Flag to specify OpenMP parallelization
OMPFLG = -fopenmp

# Preprocessor flag
PREPROCFLG = -cpp

# LAPACK and BLAS flags
LAPACKLIB = lib/liblapack.a lib/librefblas.a

# matrix-ID flag
MIDLIB = lib/id_lib.a

