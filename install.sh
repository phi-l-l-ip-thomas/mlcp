### MLCP installation script ###

# First remove libraries, if present
cd lib
rm -f *.a

# Recompile the libraries
cd INSTALL
make clean
make
cd ../BLAS
make clean
make
cd ../LAPACK
make clean
make

# Compile the code
cd ../..
#cp config/arch_gnu_mpiomp.mk arch.mk # Parallel MPI + OpenMP
cp config/arch_gnu_omp.mk arch.mk     # Serial + OpenMP
make clean
make
