
CXX = g++
CXXFLAGS = -fopenmp -O2
CFLAGS = -std=c++11

ScreenedCoulomb.exe : OpenMP_driver.cpp kernel_treecode.cpp
	$(CXX) $(CXXFLAGS) $(CFLAGS) $^ -lm -o $@


clean:
	rm -f *.o *.out *.exe *~

