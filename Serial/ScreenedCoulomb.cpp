#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sys/times.h>
#include <chrono>
#include <string>
#include <algorithm>
#include <cfloat>
#include "kernel.h"

using namespace std;
using namespace chrono;

// Constants and global variables
char PositionFile[] = "./data/scaling/rand_600000.txt";
//char PositionFile[] = "./data/TestFile.txt";
vector<panel> tree;
vector<size_t> leaf;
double xyzminmax[6];
int max_level = 0;


// Function declarations
int Count_Lines(char file_name[]);
double minval(double* x, int N_cube_in);
double maxval(double* x, int N_cube_in);
void build_tree_init();
void Swap(size_t i, size_t j, struct xyz &s);
void split_tree_node(size_t panel_index, struct xyz &particles);
void build_tree_3D_Recursive(size_t panel_index, struct xyz &particles, int level);
double mypow(double x, int n);
void Panel_Moment_B(size_t panel_index, double *lambda, struct xyz &particles, double m[Pflat]);
double Call_Treecode(double x, double y, double z, int panel_index);
double Comput_RBF(double *lambda, struct xyz &particles, size_t particle_index, size_t panel_index);
void Cluster_Chev_Points(size_t size);
double Call_Ds(int limit_1, int limit_2, int particle_index, double p_x, double p_y, double p_z,
               struct xyz &particles, double *lambda);
static const int N_cube = Count_Lines(PositionFile);
// Main program
int main() {
    // Initialize data structures
    tree.reserve(5000);
    leaf.reserve(5000);
    struct xyz particles(N_cube);
    
    // Print configuration information
    cout << "===== No box shrink ===========" << endl;
    cout << "P is " << P << endl;
    cout << "N_cube is " << N_cube << endl;
    cout << "theta is " << sqrt(sq_theta) << endl;
    cout << "N0 is " << N0 << endl;

    // Read particle positions
    FILE *fp = fopen(PositionFile, "r");
    if (fp == NULL) {
        cerr << "Cannot open random points file" << endl;
        return EXIT_FAILURE;
    }

    double x1, x2, x3;
    for (int count = 0; count < N_cube; count++) {
        if (fscanf(fp, "%lf%lf%lf", &x1, &x2, &x3) != 3) {
            cerr << "Error reading values from file at line " << count + 1 << endl;
            fclose(fp);
            return EXIT_FAILURE;
        }

        particles.x[count] = x1;
        particles.y[count] = x2;
        particles.z[count] = x3;
        particles.index[count] = -1;
        particles.old_index[count] = count;
    }
    fclose(fp);

    // Find bounds of Cartesian box
    xyzminmax[0] = minval(particles.x, N_cube);
    xyzminmax[1] = maxval(particles.x, N_cube);
    xyzminmax[2] = minval(particles.y, N_cube);
    xyzminmax[3] = maxval(particles.y, N_cube);
    xyzminmax[4] = minval(particles.z, N_cube);
    xyzminmax[5] = maxval(particles.z, N_cube);

    // Read lambda values
    double *lambda = new double[N_cube];
    char lambda_Str_data_file[64];
    snprintf(lambda_Str_data_file, sizeof(lambda_Str_data_file), "./data/scaling/lambda_%d.txt", N_cube);
    
    fp = fopen(lambda_Str_data_file, "r");
    if (fp == NULL) {
        cerr << "Cannot open lambda file" << endl;
        delete[] lambda;
        return EXIT_FAILURE;
    }

    for (int count = 0; count < N_cube; count++) {
        if (fscanf(fp, "%lf", &x1) != 1) {
            cerr << "Error reading lambda values from file at line " << count + 1 << endl;
            fclose(fp);
            delete[] lambda;
            return EXIT_FAILURE;
        }
        lambda[count] = x1;
    }
    fclose(fp);

    // Build the tree structure
    auto Start_total = system_clock::now();
    auto Start_btree = system_clock::now();
    
    build_tree_init();
    build_tree_3D_Recursive(0, particles, 0);
    
    auto End_btree = system_clock::now();
    duration<double> build_tree_time = End_btree - Start_btree;

    // Compute moments for each panel
    size_t size = tree.size();
    Cluster_Chev_Points(size);
    
    for (size_t i = 1; i < size; i++) { // skip root
        Panel_Moment_B(i, lambda, particles, tree[i].moments);
    }

    // Compute potential using treecode
    double *ptntl = new double[N_cube];
    double *ptntl_old = new double[N_cube];
    
    for (int i = 0; i < N_cube; i++) {
        ptntl[i] = Comput_RBF(lambda, particles, i, 0);
    }

    auto End_total = system_clock::now();
    duration<double> treecode_cpu_time = End_total - Start_total;

    // Output timing information
    cout << "Treecode cpu time: " << treecode_cpu_time.count() << " seconds" << endl;
    cout << "Build tree time is: " << build_tree_time.count() << " seconds" << endl;

    // Direct summation for verification
    double *p_true = new double[N_cube];
    auto Start_ds = system_clock::now();
    
    for (size_t i = 0; i < N_cube; i++) {
        double temp_x = particles.x[i];
        double temp_y = particles.y[i];
        double temp_z = particles.z[i];
        double sum = 0.0;
        
        for (size_t j = 0; j < N_cube; j++) {
            if (i == j) continue;
            
            double xx = temp_x - particles.x[j];
            double yy = temp_y - particles.y[j];
            double zz = temp_z - particles.z[j];
            
            double R = sqrt(xx * xx + yy * yy + zz * zz);
            sum += lambda[j] * (1.0/R) * exp(-kappa * R);
        }
        p_true[i] = sum;
    }
    
    auto End_ds = system_clock::now();
    duration<double> ds_cpu_time = End_ds - Start_ds;
    cout << "Direct summation time is " << ds_cpu_time.count() << " seconds" << endl;

    // Calculate errors
    for (size_t i = 0; i < N_cube; i++) {
        ptntl_old[i] = ptntl[i];
    }

    // L_infinity error
    double max_n = 0.0;
    double max_d = 0.0;
    
    for (size_t i = 0; i < N_cube; i++) {
        double diff = fabs(p_true[i] - ptntl_old[i]);
        double val = fabs(ptntl_old[i]);
        
        if (diff > max_n) max_n = diff;
        if (val > max_d) max_d = val;
    }
    
    double E = max_n / max_d;
    cout << "L_inf error is " << E << endl;
//    cout << "Max numerator is " << max_n << endl;
//    cout << "Max denominator is " << max_d << endl;

    // L2 error
    double sum_n = 0.0;
    double sum_d = 0.0;
    
    for (int i = 0; i < N_cube; i++) {
        double diff = p_true[i] - ptntl_old[i];
        sum_n += diff * diff;
        sum_d += p_true[i] * p_true[i];
    }

    double err2 = sqrt(sum_n / sum_d);
    cout << "L2 error is " << err2 << endl;
    cout << "Tree depth is " << max_level << endl;

    // Clean up
    delete[] lambda;
    delete[] ptntl;
    delete[] ptntl_old;
    delete[] p_true;
    
    cout << "Done" << endl;
    return EXIT_SUCCESS;
}
