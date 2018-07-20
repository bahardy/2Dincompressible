EXEC=2Dincomp
CC=mpicc

IDIR1=${PETSC_DIR}/${PETSC_ARCH}/include/
IDIR2=${PETSC_DIR}/include/
IDIR_LOCAL=${PWD}/
CFLAGS=-std=c99 -Wall   
CFLAGS+= -g -fopenmp


#LIBS=-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc

LIBS=/home/ucl/mema/jonathan/easybuild/software/PETSc/3.8.2-foss-2017a/lib/libpetsc.a -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/Hypre/2.11.1-foss-2017a/lib -L/home/ucl/mema/jonathan/easybuild/software/Hypre/2.11.1-foss-2017a/lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/ScaLAPACK/2.0.2-gompi-2017a-OpenBLAS-0.2.19-LAPACK-3.7.0/lib -L/home/ucl/mema/jonathan/easybuild/software/ScaLAPACK/2.0.2-gompi-2017a-OpenBLAS-0.2.19-LAPACK-3.7.0/lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/FFTW/3.3.6-gompi-2017a/lib -L/home/ucl/mema/jonathan/easybuild/software/FFTW/3.3.6-gompi-2017a/lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/UMFPACK/Lib -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/UMFPACK/Lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/KLU/Lib -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/KLU/Lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/CHOLMOD/Lib -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/CHOLMOD/Lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/BTF/Lib -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/BTF/Lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/CCOLAMD/Lib -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/CCOLAMD/Lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/COLAMD/Lib -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/COLAMD/Lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/CAMD/Lib -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/CAMD/Lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/AMD/Lib -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/AMD/Lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/SuiteSparse_config -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/SuiteSparse_config -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/OpenBLAS/0.2.19-GCC-6.3.0-2.27-LAPACK-3.7.0/lib -L/home/ucl/mema/jonathan/easybuild/software/OpenBLAS/0.2.19-GCC-6.3.0-2.27-LAPACK-3.7.0/lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/ParMETIS/4.0.3-foss-2017a/lib -L/home/ucl/mema/jonathan/easybuild/software/ParMETIS/4.0.3-foss-2017a/lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/METIS/5.1.0-foss-2017a/lib -L/home/ucl/mema/jonathan/easybuild/software/METIS/5.1.0-foss-2017a/lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/SCOTCH/6.0.4-foss-2017a/lib -L/home/ucl/mema/jonathan/easybuild/software/SCOTCH/6.0.4-foss-2017a/lib -L/home/ucl/mema/jonathan/easybuild/software/hwloc/1.11.5-GCC-6.3.0-2.27/lib -L/home/ucl/mema/jonathan/easybuild/software/OpenMPI/2.0.2-GCC-6.3.0-2.27/lib -L/home/users/j/o/jonathan/easybuild/software/GCCcore/6.3.0/lib/gcc/x86_64-pc-linux-gnu/6.3.0 -L/home/users/j/o/jonathan/easybuild/software/GCCcore/6.3.0/lib/gcc -L/home/ucl/mema/jonathan/easybuild/software/GCCcore/6.3.0/lib64 -L/home/users/j/o/jonathan/easybuild/software/GCCcore/6.3.0/lib64 -L/home/ucl/mema/jonathan/easybuild/software/ncurses/6.0-GCCcore-6.3.0/lib -L/home/ucl/mema/jonathan/easybuild/software/SuiteSparse/4.5.5-foss-2017a-ParMETIS-4.0.3/lib -L/home/ucl/mema/jonathan/easybuild/software/Boost/1.65.1-foss-2017a/lib -L/home/ucl/mema/jonathan/easybuild/software/zlib/1.2.11-GCCcore-6.3.0/lib -L/home/ucl/mema/jonathan/easybuild/software/bzip2/1.0.6-GCCcore-6.3.0/lib -L/home/ucl/mema/jonathan/easybuild/software/numactl/2.0.11-GCC-6.3.0-2.27/lib -L/home/ucl/mema/jonathan/easybuild/software/binutils/2.27-GCCcore-6.3.0/lib -L/home/ucl/mema/jonathan/easybuild/software/GCCcore/6.3.0/lib -L/home/users/j/o/jonathan/easybuild/software/GCCcore/6.3.0/lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/hwloc/1.11.5-GCC-6.3.0-2.27/lib -Wl,-rpath,/home/ucl/mema/jonathan/easybuild/software/OpenMPI/2.0.2-GCC-6.3.0-2.27/lib -lHYPRE -lscalapack -lopenblas -lgfortran -lfftw3_mpi -lfftw3 -lumfpack -lklu -lcholmod -lbtf -lccolamd -lcolamd -lcamd -lamd -lsuitesparseconfig -lopenblas -lgfortran -lparmetis -lmetis -lptesmumps -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lgfortran -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lgfortran -lm -lgfortran -lquadmath -lstdc++ -ldl -lm -lpthread -lrt -lmpi -lgcc_s -ldl -lm -lpthread

all: $(EXEC)

2Dincomp: main.o collision.o fields_memory.o flow_solver.o forces.o particle_motion.o poisson.o set_up.o write.o
	 $(CC) -I$(IDIR1) -I$(IDIR2) -I$(IDIR_LOCAL) -o $@ $^ $(CFLAGS) $(LIBS)

%.o: %.c %.h
	$(CC) -I$(IDIR_LOCAL) -c $< -o $@ $(CFLAGS)

#perio2: perio.o write.o fields.o poisson.o forces.o
#	$(CC) -I$(IDIR1) -I$(IDIR2) -I$(IDIR_LOCAL) -o $@ $^ $(CFLAGS) $(LIBS)

#perio.o: periodic.c main.h
#	$(CC) -I$(IDIR_LOCAL) -c $< -o $@ $(CFLAGS) $(LIBS)


clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
