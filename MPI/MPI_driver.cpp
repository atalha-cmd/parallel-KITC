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
    int ierr, myid, numprocs;

    ierr = MPI_Init(&argc, &argv);
    if (ierr != MPI_SUCCESS) {
        std::cerr << "Error in calling MPI_Init\n";
        return 1;
    }

    ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    if (ierr != MPI_SUCCESS) {
    std::cerr << "Error in calling MPI_Comm_size\n";
    return 1;
    }

    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (ierr != MPI_SUCCESS) {
    std::cerr << "Error in calling MPI_Comm_rank\n";
    return 1;
    }
    
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
    
    tree.reserve(5000);
    leaf.reserve(5000);
    struct xyz particles(N_cube);
    
    // Read particles data (only rank 0 reads and broadcasts)
    if (myid == 0) {
        FILE *fp = fopen(PositionFile, "r");
        if (fp == NULL) {
            cerr << "Cannot open random points file" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        double x1, x2, x3;
        for (int count = 0; count < N_cube; count++) {
            if (fscanf(fp, "%lf%lf%lf", &x1, &x2, &x3) != 3) {
                cerr << "Error reading values from file at line " << count + 1 << endl;
                fclose(fp);
                MPI_Abort(MPI_COMM_WORLD, 1);
                return 1;
            }

            particles.x[count] = x1;
            particles.y[count] = x2;
            particles.z[count] = x3;
            particles.index[count] = -1;
            particles.old_index[count] = count;
        }
        fclose(fp);
    }
    
    // Broadcast particles data to all processes
    ierr = MPI_Bcast(particles.x, N_cube, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (ierr != 0) {
        std::cerr << " error in MPI_Bcast = " << ierr << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    ierr = MPI_Bcast(particles.y, N_cube, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (ierr != 0) {
        std::cerr << " error in MPI_Bcast = " << ierr << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    ierr = MPI_Bcast(particles.z, N_cube, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (ierr != 0) {
        std::cerr << " error in MPI_Bcast = " << ierr << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    ierr = MPI_Bcast(particles.index, N_cube, MPI_INT, 0, MPI_COMM_WORLD);
    if (ierr != 0) {
        std::cerr << " error in MPI_Bcast = " << ierr << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    ierr = MPI_Bcast(particles.old_index, N_cube, MPI_INT, 0, MPI_COMM_WORLD);
    if (ierr != 0) {
        std::cerr << " error in MPI_Bcast = " << ierr << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    // Calculate min/max values (done by all ranks)
    xyzminmax[0] = minval(particles.x, N_cube);
    xyzminmax[1] = maxval(particles.x, N_cube);
    xyzminmax[2] = minval(particles.y, N_cube);
    xyzminmax[3] = maxval(particles.y, N_cube);
    xyzminmax[4] = minval(particles.z, N_cube);
    xyzminmax[5] = maxval(particles.z, N_cube);
    
    // Read lambda values (only rank 0 reads and broadcasts)
    double *lambda = new double[N_cube];
    if (myid == 0) {
        char lambda_Str_data_file[64];
        snprintf(lambda_Str_data_file, sizeof(lambda_Str_data_file), "./data/scaling/lambda_%d.txt", N_cube);

        FILE *fp = fopen(lambda_Str_data_file, "r");
        if (fp == NULL) {
            cerr << "Cannot open lambda file" << endl;
            delete[] lambda;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        double x1;
        for (int count = 0; count < N_cube; count++) {
            if (fscanf(fp, "%lf", &x1) != 1) {
                cerr << "Error reading lambda values from file at line " << count + 1 << endl;
                fclose(fp);
                delete[] lambda;
                MPI_Abort(MPI_COMM_WORLD, 1);
                return 1;
            }
            lambda[count] = x1;
        }
        fclose(fp);
    }
    ierr = MPI_Bcast(lambda, N_cube, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (ierr != 0) {
        std::cerr << " error in MPI_Bcast = " << ierr << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // Timing for tree building
    if (myid == 0) {
        Start_total = MPI_Wtime();
        Start_btree = MPI_Wtime();
    }
    
    // All ranks build the tree together 
    build_tree_init();
    build_tree_3D_Recursive(0, particles, 0);
    
    if (myid == 0) {
        End_btree = MPI_Wtime();
        build_tree_time = End_btree - Start_btree;
    }
    
    size_t size = tree.size();
    Cluster_Chev_Points(size);

    for (size_t i = 1; i < size; i++) {
        Panel_Moment_B(i, lambda, particles, tree[i].moments);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Compute this processor's interval (is to ie)
    int is = (int)(1.0 * N_cube / numprocs) * myid;
    int ie = (int)(1.0 * N_cube / numprocs) * (myid + 1) - 1;
    if (myid == numprocs - 1) ie = N_cube - 1;

    int local_count = ie - is + 1;

    // Allocate local potential array
    double *local_ptntl = new double[local_count]();

    // Compute local potentials
    for (int i = 0; i < local_count; i++) {
        int global_index = is + i;
        local_ptntl[i] = Comput_RBF(lambda, particles, global_index, 0);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Prepare MPI_Gatherv buffers
    double *ptntl = nullptr;
    double *ptntl_old = nullptr;
    int *recvcounts = nullptr;
    int *displs = nullptr;

    // Root process allocates and prepares buffers
    if (myid == 0) {
        ptntl = new double[N_cube]();        // Full potential array
        ptntl_old = new double[N_cube]();
        recvcounts = new int[numprocs];      // How many elements each process sends
        displs = new int[numprocs];          // Displacements for gathering

        for (int i = 0; i < numprocs; i++) {
            int js = (int)(1.0 * N_cube / numprocs) * i;
            int je = (int)(1.0 * N_cube / numprocs) * (i + 1) - 1;
            if (i == numprocs - 1) je = N_cube - 1;
            recvcounts[i] = je - js + 1;
            displs[i] = js;
        }
    }

    // Gather local potentials to root
    ierr = MPI_Gatherv(local_ptntl, local_count, MPI_DOUBLE, ptntl, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (ierr != 0) {
        std::cerr << "Error in MPI_Gatherv = " << ierr << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    if (myid == 0) {
        End_total = MPI_Wtime();
        treecode_cpu_time = End_total - Start_total;
        cout << "Treecode cpu time: " << treecode_cpu_time << " seconds" << endl;
        cout << "Build tree time is: " << build_tree_time << " seconds" << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Direct summation computation
    if (myid == 0){
        Start_ds = MPI_Wtime();
    }
    int start_particle = floor(1.0 * N_cube / numprocs) * myid ;
    int end_particle = floor(1.0 * N_cube / numprocs) * (myid + 1) ;
    if(myid == numprocs - 1) end_particle = N_cube ;
    
    double *local_ds = new double[N_cube]();
    for (int i = start_particle; i < end_particle; i++) {
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
    MPI_Barrier(MPI_COMM_WORLD);
    double *p_true = nullptr;
    if (myid == 0) {
        p_true = new double[N_cube]();
    }
    ierr = MPI_Reduce(local_ds, p_true, N_cube, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (ierr != 0) {
        std::cerr << "Error in MPI_Reduce = " << ierr << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // Calculate errors (only on rank 0)
    if (myid == 0) {
        End_ds = MPI_Wtime();
        ds_cpu_time = End_ds - Start_ds;

        cout << "Direct summation time is : " << ds_cpu_time << " seconds" << endl;
        
        // Copy potentials
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
        
        // L2 error
        double sum_n = 0.0;
        double sum_d = 0.0;
        
        for (int i = 0; i < N_cube; i++) {
            double diff = p_true[i] - ptntl_old[i];
            sum_n += diff * diff;
            sum_d += p_true[i] * p_true[i];
        }
        
        double err2 = sqrt(sum_n / sum_d);
        
        cout << "L_inf error is " << E << endl;
        cout << "L2 error is " << err2 << endl;
        cout << "Tree depth is " << max_level << endl;
        cout << "Done" << endl;
    }
    
    // Clean up
    delete[] lambda;
    delete[] local_ptntl;
    delete[] local_ds;
    if (myid == 0) {
        delete[] ptntl;
        delete[] p_true;
        delete[] ptntl_old;
        delete[] recvcounts;
        delete[] displs;
    }
    MPI_Finalize();
    return 0;
}
