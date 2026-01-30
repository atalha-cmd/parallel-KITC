

#module load gcc/11.2.0 openmpi
#mpirun -np 4 ./ScreenedCoulomb.exe

MPICXX = mpicxx
CXXFLAGS = -O2
CFLAGS = -std=c++11

ScreenedCoulomb.exe : MPI_driver.cpp kernel_treecode.cpp
#cyclicCoulomb.exe : cyclic_MPI.cpp kernel_treecode.cpp
	$(MPICXX) $(CXXFLAGS) $(CFLAGS) $^ -lm -o $@


clean:
	rm -f *.o *.out *.exe *~

