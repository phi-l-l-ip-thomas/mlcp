# arch.mk for NERSC Perlmutter; to set up environment:
#
# module load PrgEnv-nvidia

# Fortran compiler
FC = ftn

# Flags to always include
FOPTS = -Mr8

OPTFLG = -O3

# Debug flags
DEBUG = no
DEBUGFLG = -O0 -traceback

# Flag to specify the position of mod files
MODULEFLG = -I

# Flag to specify Message Passing Interface parallelization
MPIFLG = 

# Flag to specify OpenMP parallelization
USEOMP = yes
OMPFLG = -mp

# Flag to specify OpenACC parallelization
USEACC = yes
ACCFLG = -cuda -acc -fast -Minfo=accel -cudalib=cublas,cusolver,cutensor,nvtx3 -gpu=cc80

# Preprocessor flag
PREPROCFLG = -Mpreprocess

# LAPACK and BLAS flags
LAPACKLIB = 

# CUDA flags
CULIB = 

