# arch.mk for TACC Vista; to set up environment:
# 
# module use /home1/apps/nvidia/modulefiles
# module purge
# module load nvhpc/24.7
# module load hpctoolkit

# Fortran compiler
FC = mpifort

# Flags to always include
FOPTS = -Mr8

OPTFLG = -O3

# Debug flags
DEBUG = no
DEBUGFLG = -O0 -traceback

# Flag to specify the position of mod files
MODULEFLG = -I

# Flag to specify Message Passing Interface parallelization
USEMPI = yes
MPIFLG = 

# Flag to specify OpenMP parallelization
USEOMP = yes
OMPFLG = -mp

# Flag to specify OpenACC parallelization
USEACC = yes
ACCFLG = -cuda -acc -fast -Minfo=accel -cudalib=cublas,cusolver,cutensor,nvtx3

# Preprocessor flag
PREPROCFLG = -Mpreprocess

# LAPACK and BLAS flags
LAPACKLIB = -L$(NVHPC_ROOT)/compilers/lib -llapack -lblas

# CUDA flags
CULIB = 

