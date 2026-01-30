
#ifndef KERNEL_H
#define KERNEL_H

#include <vector>
using namespace std;

constexpr int P = 5;  // Order of Taylor approximation
constexpr int Pflat = (P + 1) * (P + 1) * (P + 1);  // Computed from P
constexpr int N0 = 500;
constexpr double sq_theta = 0.25; // theta = 0.5
constexpr double DEL = 0.02;
constexpr double kappa = 1.0;
constexpr bool UseSleep = false; 

struct xyz // particle coordinates (physical)
{
	double* x;
	double* y;
	double* z;
	size_t* index;
	size_t* old_index;
	size_t size;
	xyz(size_t N_cube_in)
	{
		size = N_cube_in;
		x = new double[size];
		y = new double[size];
		z = new double[size];
		index = new size_t[size];
		old_index = new size_t[size];
	}
	~xyz()
	{
		delete[] x;
		delete[] y;
		delete[] z;
		delete[] index;
		delete[] old_index;
	}
};

struct panel {
    size_t members[2];
    vector<size_t> children;
    double xinterval[2];
    double yinterval[2];
    double zinterval[2];
    double xc; // panel center x coordinate
    double yc; // panel center y coordinate
    double zc; // panel center z coordinate
    double MAC; // r^2 / theta^2
    double moments[Pflat];
    int moment_flag;
    double t1[P + 1]; // interpolation points in x direction
    double t2[P + 1];
    double t3[P + 1];

    panel() // initialization
    {
        moment_flag = 0;
        members[0] = 0;
        members[1] = -1;

        // for (size_t kk = 0; kk < Pflat + 1; kk++)
        for (size_t kk = 0; kk < Pflat; kk++)
            moments[kk] = 0.0;

        for (int i = 0; i < P + 1; i++) {
            t1[i] = 0.0;
            t2[i] = 0.0;
            t3[i] = 0.0;
        }
    }
};

// Global variables (extern declarations)
extern int max_level;
extern vector<panel> tree;
extern vector<size_t> leaf;
extern double xyzminmax[6];
extern char PositionFile[];

#endif // KERNEL_H
