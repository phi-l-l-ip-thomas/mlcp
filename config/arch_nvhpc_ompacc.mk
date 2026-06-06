# Fortran compiler
FC = nvfortran

# Flags to always include
FOPTS = -Mr8

OPTFLG = -O3

# Debug flags
DEBUG = no
DEBUGFLG = -O0 -traceback

# Flag to specify the position of mod files
MODULEFLG = -I

# Flag to specify Message Passing Interface parallelization
USEMPI = no
MPIFLG =

# Flag to specify OpenMP parallelization
USEOMP = yes
OMPFLG = -mp

# Flag to specify OpenACC parallelization
USEACC = yes
ACCFLG = -cuda -acc -fast -Minfo=accel -cudalib=cublas,cusolver,cutensor,nvtx3 -gpu=cc60

# Preprocessor flag
PREPROCFLG = -Mpreprocess

# LAPACK and BLAS flags
#LAPACKLIB = lib/liblapack.a lib/librefblas.a
LAPACKLIB = -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/compilers/lib -lblas -llapack

# CUDA flags
CULIB = #-L/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/math_libs/lib64 -lcublas -lcusolver -lcutensor

