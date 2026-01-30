
CXX = g++
CXXFLAGS = -O2
CFLAGS = -std=c++11

ScreenedCoulomb.exe : ScreenedCoulomb.cpp kernel_treecode.cpp
	$(CXX) $(CXXFLAGS) $(CFLAGS) $^ -lm -o $@


clean:
	rm -f *.o *.out *.exe *~

