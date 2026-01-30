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
#include <mpi.h>

using namespace std;
using namespace chrono;

// Constants and global variables
char PositionFile[] = "./data/scaling/rand_100000.txt";
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

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    
    int myid, numprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    // Timing variables
    double Start_total, Start_btree, End_btree, End_total, Start_ds, End_ds;
    double build_tree_time, treecode_cpu_time, ds_cpu_time;    
    if (myid == 0) {
        cout << "===== No box shrink ===========" << endl;
        cout << "P is " << P << endl;
        cout << "N_cube is " << N_cube << endl;
        cout << "theta is " << sqrt(sq_theta) << endl;
        cout << "N0 is " << N0 << endl;
        cout << "Parallel Starting using " << numprocs << " MPI processes" << endl;
    }

    // Initialize data structures
    tree.reserve(5000);
    leaf.reserve(5000);
    struct xyz particles(N_cube);
    double* lambda = new double[N_cube];

    // Read and broadcast particle data
    if (myid == 0) {
        FILE* fp = fopen(PositionFile, "r");
        if (!fp) {
            cerr << "Cannot open points file" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        for (int i = 0; i < N_cube; i++) {
            if (fscanf(fp, "%lf %lf %lf", &particles.x[i], &particles.y[i], &particles.z[i]) != 3) {
                cerr << "Error reading line " << i+1 << endl;
                fclose(fp);
                MPI_Abort(MPI_COMM_WORLD, 1);
                return 1;
            }
            particles.index[i] = -1;
            particles.old_index[i] = i;
        }
        fclose(fp);

        // Read lambda values
        char lambda_file[64];
        snprintf(lambda_file, sizeof(lambda_file), "./data/scaling/lambda_%d.txt", N_cube);
        fp = fopen(lambda_file, "r");
        if (!fp) {
            cerr << "Cannot open lambda file" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        for (int i = 0; i < N_cube; i++) {
            if (fscanf(fp, "%lf", &lambda[i]) != 1) {
                cerr << "Error reading lambda at line " << i+1 << endl;
                fclose(fp);
                MPI_Abort(MPI_COMM_WORLD, 1);
                return 1;
            }
        }
        fclose(fp);
    }

    // Broadcast data to all processes
    MPI_Bcast(particles.x, N_cube, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(particles.y, N_cube, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(particles.z, N_cube, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(particles.index, N_cube, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(particles.old_index, N_cube, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(lambda, N_cube, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Calculate min/max values
    xyzminmax[0] = minval(particles.x, N_cube);
    xyzminmax[1] = maxval(particles.x, N_cube);
    xyzminmax[2] = minval(particles.y, N_cube);
    xyzminmax[3] = maxval(particles.y, N_cube);
    xyzminmax[4] = minval(particles.z, N_cube);
    xyzminmax[5] = maxval(particles.z, N_cube);

    // Build tree
    if (myid == 0){
        Start_total = MPI_Wtime();
        Start_btree = MPI_Wtime();
    }
    
    build_tree_init();
    build_tree_3D_Recursive(0, particles, 0);

    if (myid == 0){
        End_btree = MPI_Wtime();
        build_tree_time = End_btree - Start_btree;
    }

    // Prepare tree for computation
    size_t tree_size = tree.size();
    Cluster_Chev_Points(tree_size);
    for (size_t i = 1; i < tree_size; i++) {
        Panel_Moment_B(i, lambda, particles, tree[i].moments);
    }

    // ========== Treecode Computation (Cyclic) ==========
    // Calculate how many particles each process handles
    int local_count = (N_cube + numprocs - 1) / numprocs;
    int actual_local = 0;
    for (int i = myid; i < N_cube; i += numprocs) {
        actual_local++;
    }

    double* local_ptntl = new double[actual_local]();
    int idx = 0;
    for (int i = myid; i < N_cube; i += numprocs) {
        local_ptntl[idx++] = Comput_RBF(lambda, particles, i, 0);
    }

    // Prepare gathering
    int* recvcounts = nullptr;
    int* displs = nullptr;
    double* ptntl = nullptr;
    double* ptntl_old = nullptr;

    if (myid == 0) {
        recvcounts = new int[numprocs];
        displs = new int[numprocs];
        ptntl = new double[N_cube]();
        ptntl_old = new double[N_cube]();

        // Calculate receive counts and displacements
        for (int p = 0; p < numprocs; p++) {
            recvcounts[p] = 0;
            for (int i = p; i < N_cube; i += numprocs) {
                recvcounts[p]++;
            }
            displs[p] = (p == 0) ? 0 : displs[p-1] + recvcounts[p-1];
        }
    }

    // Gather all results to root
    MPI_Gatherv(local_ptntl, actual_local, MPI_DOUBLE,
               ptntl, recvcounts, displs, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    // Reorder the gathered data on root
    if (myid == 0) {
        int* indices = new int[numprocs]();
        for (int i = 0; i < N_cube; i++) {
            int p = i % numprocs;
            int pos = displs[p] + indices[p]++;
            ptntl_old[i] = ptntl[pos];
        }
        delete[] indices;
    }
    if (myid == 0) {
        End_total = MPI_Wtime();
        treecode_cpu_time = End_total - Start_total;
        cout << "Treecode time: " << treecode_cpu_time << " seconds" << endl;
        cout << "Build tree time: " << build_tree_time << " seconds" << endl;
    }

    // ========== Direct Summation (Cyclic) ==========
    if (myid == 0){
        Start_ds = MPI_Wtime();
    }
    
    double* local_ds = new double[N_cube]();
    for (int i = myid; i < N_cube; i += numprocs) {
        double sum = 0.0;
        for (int j = 0; j < N_cube; j++) {
            if (i == j) continue;
            double dx = particles.x[i] - particles.x[j];
            double dy = particles.y[i] - particles.y[j];
            double dz = particles.z[i] - particles.z[j];
            double R = sqrt(dx*dx + dy*dy + dz*dz);
            sum += lambda[j] * exp(-kappa * R) / R;
        }
        local_ds[i] = sum;
    }

    double* p_true = nullptr;
    if (myid == 0) {
        p_true = new double[N_cube]();
    }
    MPI_Reduce(local_ds, p_true, N_cube, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Calculate errors on root
    if (myid == 0) {
        End_ds = MPI_Wtime();
        ds_cpu_time = End_ds - Start_ds;
        cout << "Direct summation time: " << ds_cpu_time << " seconds" << endl;

        // L_infinity error
        double max_n = 0.0, max_d = 0.0;
        for (int i = 0; i < N_cube; i++) {
            double diff = fabs(p_true[i] - ptntl_old[i]);
            double val = fabs(ptntl_old[i]);
            max_n = max(max_n, diff);
            max_d = max(max_d, val);
        }
        double E_inf = max_n / max_d;

        // L2 error
        double sum_n = 0.0, sum_d = 0.0;
        for (int i = 0; i < N_cube; i++) {
            double diff = p_true[i] - ptntl_old[i];
            sum_n += diff * diff;
            sum_d += p_true[i] * p_true[i];
        }
        double E_2 = sqrt(sum_n / sum_d);

        cout << "L_inf error: " << E_inf << endl;
        cout << "L2 error: " << E_2 << endl;
        cout << "Tree depth: " << max_level << endl;
    }

    // Clean up
    delete[] lambda;
    delete[] local_ptntl;
    delete[] local_ds;
    if (myid == 0) {
        delete[] ptntl;
        delete[] ptntl_old;
        delete[] p_true;
        delete[] recvcounts;
        delete[] displs;
    }

    MPI_Finalize();
    return 0;
}