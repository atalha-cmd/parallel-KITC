#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sys/times.h>
#include <string>
#include <algorithm> //*minmax
#include <cfloat> // DBL_MIN
#include "kernel.h"

using namespace std;

//*****************************************************************************//
int Count_Lines(char file_name[])
{
	int numLines = 0;
	ifstream input_file;
	string str;
	input_file.open(file_name, ios::in);
	if(input_file.fail())
	{
		cout << "The file does not exist under this directory." << endl;
		return 0;
	}
	else
	{
		while (getline(input_file,str))
		{
			if (str.empty()==false)
			{
				numLines++;
			}
		}
		input_file.close();
		return numLines;
	}
}

// Compute Interval and Size of 
double minval(double* x, int N_cube_in)
{
  double MinVal = x[0];
  for (int i = 1; i < N_cube_in; i++) {
    if (MinVal > x[i])
      MinVal = x[i];
  }
  MinVal = MinVal;
  
  return MinVal;
}

double maxval(double* x, int N_cube_in)
{
  double MaxVal = x[0];
  for (int i = 1; i < N_cube_in; i++) {
    if (MaxVal < x[i])
      MaxVal = x[i];
  }

  MaxVal = MaxVal;

  return MaxVal;
}

static size_t node_count = 0;
//*****************************************************************************//
void build_tree_init() 
{
	panel temp_panel;

	// indices of particles belonging to panel
	temp_panel.members[0] = 0;
	temp_panel.members[1] = Count_Lines(PositionFile) - 1;

	// Interval defining the panel
  	temp_panel.xinterval[0] = xyzminmax[0];
  	temp_panel.xinterval[1] = xyzminmax[1];
  	temp_panel.yinterval[0] = xyzminmax[2];
  	temp_panel.yinterval[1] = xyzminmax[3];
  	temp_panel.zinterval[0] = xyzminmax[4];
  	temp_panel.zinterval[1] = xyzminmax[5];
  	//cout << "y min is " << xyzminmax[2] << endl;
  	//cout << "y max is " << xyzminmax[3] << endl;

    temp_panel.xc = 0.5 * (temp_panel.xinterval[0] + temp_panel.xinterval[1]);
    temp_panel.yc = 0.5 * (temp_panel.yinterval[0] + temp_panel.yinterval[1]);
    temp_panel.zc = 0.5 * (temp_panel.zinterval[0] + temp_panel.zinterval[1]);

    double xL = temp_panel.xinterval[1] - temp_panel.xinterval[0];
  	double yL = temp_panel.yinterval[1] - temp_panel.yinterval[0];
 	double zL = temp_panel.zinterval[1] - temp_panel.zinterval[0];

 	//cout<< "yL is"<< yL <<endl;

 	double sq_r = xL * xL + yL * yL + zL * zL; // r^2
  	temp_panel.MAC = sq_r / sq_theta; // MAC = r^2 / theta^2

	tree.push_back(temp_panel);
	node_count = 1;
}

//*****************************************************************************//
void Swap(size_t i, size_t j, struct xyz &s)
{
	if (i == j)
		return;

	double x = s.x[i];
	double y = s.y[i];
	double z = s.z[i];
	size_t index = s.index[i];
	size_t old_index = s.old_index[i];

	s.x[i] = s.x[j];
	s.y[i] = s.y[j];
	s.z[i] = s.z[j];
	s.index[i] = s.index[j];
	s.old_index[i] = s.old_index[j];

	s.x[j] = x;
	s.y[j] = y;
	s.z[j] = z;
	s.index[j] = index;
	s.old_index[j] = old_index;
}
//*****************************************************************************//
void split_tree_node(size_t panel_index, struct xyz &particles)
{
	panel child[8];

	double tp_x0 = tree[panel_index].xinterval[0];
	double tp_x1 = tree[panel_index].xinterval[1];
	double tp_y0 = tree[panel_index].yinterval[0];
	double tp_y1 = tree[panel_index].yinterval[1];
	double tp_z0 = tree[panel_index].zinterval[0];
	double tp_z1 = tree[panel_index].zinterval[1];

	double midpointx = (tp_x0 + tp_x1) / 2.0;
	double midpointy = (tp_y0 + tp_y1) / 2.0;
    double midpointz = (tp_z0 + tp_z1) / 2.0;

	double xc0 = (tp_x0 + midpointx) / 2.0;
	double xc1 = (tp_x1 + midpointx) / 2.0;
	double yc0 = (tp_y0 + midpointy) / 2.0;
	double yc1 = (tp_y1 + midpointy) / 2.0;
	double zc0 = (tp_z0 + midpointz) / 2.0;
	double zc1 = (tp_z1 + midpointz) / 2.0;

	child[0].xinterval[0] = tp_x0;
	child[0].xinterval[1] = midpointx;
	child[0].yinterval[0] = tp_y0;
	child[0].yinterval[1] = midpointy;
	child[0].zinterval[0] = tp_z0;
	child[0].zinterval[1] = midpointz;
	child[0].xc = xc0;
	child[0].yc = yc0;
	child[0].zc = zc0;
	child[0].MAC = ((midpointx - xc0) * (midpointx - xc0) + (midpointy - yc0) * (midpointy - yc0) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;

	child[1].xinterval[0] = midpointx;
	child[1].xinterval[1] = tp_x1;
	child[1].yinterval[0] = tp_y0;
    child[1].yinterval[1] = midpointy;
	child[1].zinterval[0] = tp_z0;
	child[1].zinterval[1] = midpointz;
	child[1].xc = xc1;
	child[1].yc = yc0;
	child[1].zc = zc0;
	child[1].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (midpointy - yc0) * (midpointy - yc0) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;

	child[2].xinterval[0] = tp_x0;
	child[2].xinterval[1] = midpointx;
	child[2].yinterval[0] = midpointy;
	child[2].yinterval[1] = tp_y1;
	child[2].zinterval[0] = tp_z0;
	child[2].zinterval[1] = midpointz;
	child[2].xc = xc0;
	child[2].yc = yc1;
	child[2].zc = zc0;
    child[2].MAC = ((midpointx - xc0) * (midpointx - xc0) + (tp_y1 - yc1) * (tp_y1 - yc1) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;

	child[3].xinterval[0] = midpointx;
	child[3].xinterval[1] = tp_x1;
	child[3].yinterval[0] = midpointy;
	child[3].yinterval[1] = tp_y1;
	child[3].zinterval[0] = tp_z0;
	child[3].zinterval[1] = midpointz;
	child[3].xc = xc1;
	child[3].yc = yc1;
	child[3].zc = zc0;
    child[3].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (tp_y1 - yc1) * (tp_y1 - yc1) + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;

	child[4].xinterval[0] = tp_x0;
	child[4].xinterval[1] = midpointx;
	child[4].yinterval[0] = tp_y0;
	child[4].yinterval[1] = midpointy;
    child[4].zinterval[0] = midpointz;
	child[4].zinterval[1] = tp_z1;
	child[4].xc = xc0;
	child[4].yc = yc0;
	child[4].zc = zc1;
	child[4].MAC = ((midpointx - xc0) * (midpointx - xc0) + (midpointy - yc0) * (midpointy - yc0) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;

	child[5].xinterval[0] = midpointx;
	child[5].xinterval[1] = tp_x1;
	child[5].yinterval[0] = tp_y0;
    child[5].yinterval[1] = midpointy;
	child[5].zinterval[0] = midpointz;
	child[5].zinterval[1] = tp_z1;
	child[5].xc = xc1;
	child[5].yc = yc0;
	child[5].zc = zc1;
	child[5].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (midpointy - yc0) * (midpointy - yc0) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;

	child[6].xinterval[0] = tp_x0;
	child[6].xinterval[1] = midpointx;
	child[6].yinterval[0] = midpointy;
	child[6].yinterval[1] = tp_y1;
	child[6].zinterval[0] = midpointz;
	child[6].zinterval[1] = tp_z1;
	child[6].xc = xc0;
	child[6].yc = yc1;
	child[6].zc = zc1;
	child[6].MAC = ((midpointx - xc0) * (midpointx - xc0) + (tp_y1 - yc1) * (tp_y1 - yc1) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;

	child[7].xinterval[0] = midpointx;
	child[7].xinterval[1] = tp_x1;
	child[7].yinterval[0] = midpointy;
	child[7].yinterval[1] = tp_y1;
	child[7].zinterval[0] = midpointz;
	child[7].zinterval[1] = tp_z1;
	child[7].xc = xc1;
	child[7].yc = yc1;
	child[7].zc = zc1;
	child[7].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1) + (tp_y1 - yc1) * (tp_y1 - yc1) + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;

	vector<size_t> v[8];
	size_t start = tree[panel_index].members[0];
	size_t end = tree[panel_index].members[1];
	size_t* addr_table = new size_t[end - start + 1];

	size_t index;
	for (index = start; index <= end; index++)
	{
		particles.index[index] = index;
		addr_table[index - start] = index;

		if (particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
			particles.z[index] <= midpointz)
			v[0].push_back(index);
		else if (particles.x[index] > midpointx && particles.y[index] <= midpointy &&
				particles.z[index] <= midpointz )
			v[1].push_back(index);
		else if (particles.x[index] <= midpointx && particles.y[index] > midpointy &&
				particles.z[index]<= midpointz)
			v[2].push_back(index);
		else if (particles.x[index] > midpointx && particles.y[index] > midpointy &&
				particles.z[index] <= midpointz)
			v[3].push_back(index);
		else if(particles.x[index] <= midpointx && particles.y[index] <= midpointy &&
				particles.z[index] > midpointz )
			v[4].push_back(index);
		else if (particles.x[index] > midpointx && particles.y[index] <= midpointy &&
				particles.z[index] > midpointz)
			v[5].push_back(index);
		else if (particles.x[index] <= midpointx && particles.y[index] > midpointy &&
				particles.z[index] > midpointz)
			v[6].push_back(index);
		else if (particles.x[index] > midpointx && particles.y[index] > midpointy &&
				particles.z[index] > midpointz)
			v[7].push_back(index);
	}

	size_t seq = start;
	for (size_t j = 0; j < 8; j++)
	{
		size_t size = v[j].size();

		if (size >= 1)
		{
			for (size_t k = 0; k < size; k++)
			{
				if (k == 0)
					child[j].members[0] = seq;
				if (k == size - 1)
					child[j].members[1] = seq;

				index = v[j][k];
				
				// This uses an address table
				size_t pos = addr_table[index - start];
				size_t out = particles.index[seq];
				Swap(pos, seq, particles);
				addr_table[index - start] = seq;
				addr_table[out - start] = pos;

				seq++;
			}

			node_count++;
			tree[panel_index].children.push_back(node_count - 1);
			tree.push_back(child[j]);
			v[j].clear();
		}
	}

	delete[] addr_table;
    
}

//*****************************************************************************//
void build_tree_3D_Recursive(size_t panel_index, struct xyz &particles, int level)
{
	if (level > max_level)
		max_level = level;
	
	size_t n = tree[panel_index].members[1] - tree[panel_index].members[0] + 1; //-1-0+1=0?

	if (n >= (size_t)N0){
		split_tree_node(panel_index, particles);

		for (size_t i = 0; i < tree[panel_index].children.size(); i++)
		{
			size_t panel_index_new = tree[panel_index].children[i];
			build_tree_3D_Recursive(panel_index_new, particles, level + 1);
		}
	}
	else
		leaf.push_back(panel_index);
}

//*****************************************************************************//

double mypow(double x, int n)
{
    double result = 1.0;
    while (n > 0) {
        if (n & 1) 
            result = result * x;
        n = n >> 1; 
        x = x * x;
    }
    
    return result;
}

//*****************************************************************************//
void Panel_Moment_B(size_t panel_index, double *lambda, struct xyz &particles,
                    double m[Pflat])
{
  // Intput : panel_index
  // Output : m: moments for panel_index^th panel
	double t1[P + 1];
	double t2[P + 1];
	double t3[P + 1];

  	int i, j, k, kk;
  	int a1exactind, a2exactind, a3exactind;

  	for (i = 0; i < P + 1; i++) {
    	t1[i] = tree[panel_index].t1[i];
    	t2[i] = tree[panel_index].t2[i];
    	t3[i] = tree[panel_index].t3[i];
  	}

  	double w1i[P + 1];
  	double dj[P + 1];
  	dj[0] = 0.5;
  	dj[P] = 0.5;
  	for (j = 1; j < P; j++)
    	dj[j] = 1.0;
  	for (j = 0; j < P + 1; j++)
    	w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];

  	double a1i[P + 1];
  	double a2j[P + 1];
  	double a3k[P + 1];

  	double x, y, z;
  	double dx, dy, dz;
  	double SumA1;
  	double SumA2;
  	double SumA3;
  	double D;
  	double s;

  	size_t tp0 = tree[panel_index].members[0];
  	size_t tp1 = tree[panel_index].members[1];
  	size_t tp_j;

  	for (tp_j = tp0; tp_j <= tp1; tp_j++) {
    	x = particles.x[tp_j];
    	y = particles.y[tp_j];
    	z = particles.z[tp_j];

    	a1exactind = -1;
    	a2exactind = -1;
    	a3exactind = -1;

    	SumA1 = 0.0;
    	SumA2 = 0.0;
    	SumA3 = 0.0;

    	for (j = 0; j < P + 1; j++) {
      		dx = x - t1[j];
      		dy = y - t2[j];
      		dz = z - t3[j];

      		if (fabs(dx) <= DBL_MIN)
        		a1exactind = j;
      		else {
        		a1i[j] = w1i[j] / dx;
        		SumA1 += a1i[j];
      		}

      		if (fabs(dy) <= DBL_MIN)
        		a2exactind = j;
      		else {
        		a2j[j] = w1i[j] / dy;
        		SumA2 += a2j[j];
      		}

      		if (fabs(dz) <= DBL_MIN)
        		a3exactind = j;
      		else {
        		a3k[j] = w1i[j] / dz;
        		SumA3 += a3k[j];
      		}
    	}

    	if (a1exactind > -1) {
      		SumA1 = 1.0;
      		for (j = 0; j < P + 1; j++) a1i[j] = 0.0;
        	a1i[a1exactind] = 1.0;
    	}

    	if (a2exactind > -1) {
      		SumA2 = 1.0;
      		for (j = 0; j < P + 1; j++) a2j[j] = 0.0;
        	a2j[a2exactind] = 1.0;
    	}

    	if (a3exactind > -1) {
      		SumA3 = 1.0;
      		for (j = 0; j < P + 1; j++) a3k[j] = 0.0;
        	a3k[a3exactind] = 1.0;
    	}

    	D = 1.0 / (SumA1 * SumA2 * SumA3);

    	kk = -1;
    	for (i = 0; i < P + 1; i++) {
      		for (j = 0; j < P + 1; j++) {
        		for (k = 0; k < P + 1; k++) {
          			kk++;
          			s = a1i[i] * a2j[j] * a3k[k] * D;
          			m[kk] += s * lambda[tp_j];
        		}
      		}
  	  	}
  	}
}

//*****************************************************************************//
double Call_Treecode(double x, double y, double z, int panel_index)
{
    double potential = 0.0;
  
    double temp_i;
    double temp_j;
    double temp_k;
    double R, R2, invR; 
    
    double dx[P + 1];
    double dy[P + 1];
    double dz[P + 1];
    
    double s;
    double temp_moments;
    double G;
    double DD;
    
    for (int i = 0; i < P + 1; i++)
    {
        dx[i] = x - tree[panel_index].t1[i];
        dy[i] = y - tree[panel_index].t2[i];
        dz[i] = z - tree[panel_index].t3[i];
    }
    
    
    int kk = -1;
    for (int i = 0; i < P + 1; i++)
    {
        temp_i = dx[i] * dx[i];
        for (int j = 0; j < P + 1; j++)
        {
            temp_j =  dy[j] * dy[j];
            for (int k = 0; k < P + 1; k++)
            {
                kk = kk + 1;
                
                temp_moments = tree[panel_index].moments[kk];
                //if (panel_index==9)
                //	printf("moments: %d %d %d %12.6e %12.6e %12.6e \n",i,j,k,temp_moments[0],temp_moments[1],temp_moments[2]);
                
                temp_k = dz[k] * dz[k];
                R2 = temp_i + temp_j + temp_k;
                R = sqrt(R2);
                invR = 1/R;
                
                potential += temp_moments * invR * exp(-kappa * R);
                
            }
        }
    }
    //printf("call Treecode: %12.6e %12.6e %12.6e \n",ptntlcity [0],ptntlcity [1],ptntlcity [2]);
    
    return potential;
}

//*****************************************************************************//
 double Call_Ds(int limit_1, int limit_2, int particle_index, double p_x, double p_y, double p_z, struct xyz &particles, double *lambda)
{
    
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    
    double ff;
    double x1,y1,z1;
    double R,R2,Rinv;
    
    
    double potential = 0.0;

    double G;
    double DD;
    for (size_t jj = limit_1; jj <= limit_2; jj++)
    {
        if (jj == particle_index)
            continue;
        
        x = p_x - particles.x[jj];
        y = p_y - particles.y[jj];
        z = p_z - particles.z[jj];
        ff = lambda[jj];
        x1 = x * x;
        y1 = y * y;
        z1 = z * z;
        
        R2 = x1 + y1 + z1;
        R = sqrt(R2);
        Rinv = 1.0/R;

        potential += ff * Rinv * exp(-kappa * R);
    }
     return potential;

}

//*****************************************************************************//

double Comput_RBF(double *lambda, struct xyz &particles,size_t particle_index, size_t panel_index)
{
    // input :
	//         lambda : RBF coefficients
	//         particles : all particles coordinates
    //         particle_index
    //         panel_index
	// output : 
	//          ptntlcity in 3D
	
    double potential = 0.0;
	
    double p_x = 0.0;
    double p_y = 0.0;
    double p_z = 0.0;
	double xc = 0.0;
	double yc = 0.0;
    double zc = 0.0;
    size_t limit_1;
	size_t limit_2;
    double R_sq;
    
    limit_1 = tree[panel_index].members[0];
    limit_2 = tree[panel_index].members[1];
    
	
	p_x = particles.x[particle_index];
	p_y = particles.y[particle_index];
    p_z = particles.z[particle_index];
	
	xc = tree[panel_index].xc;
	yc = tree[panel_index].yc;
    zc = tree[panel_index].zc;
    
    double tpx = p_x - xc;
    double tpy = p_y - yc;
    double tpz = p_z - zc;
	
    R_sq = tpx * tpx + tpy * tpy + tpz * tpz;
    //printf("%zd %zd %zd \n",panel_index,limit_1, limit_2);

    if (tree[panel_index].MAC < R_sq)
	{
        double tree_result = Call_Treecode(p_x, p_y, p_z, panel_index);
        //printf("partile-cluster-MAC: %12.6e %12.6e %12.6e \n",tree_result [0],tree_result [1],tree_result [2]);
        //printf("partile-cluster-MAC: %zd %zd %zd \n",panel_index,limit_1, limit_2);
        //printf("partile-cluster-MAC: %12.6e %12.6e %12.6e %12.6e\n",xc,yc,zc,R_sq);
        potential = potential + tree_result;
    }

    else
	{
		if (limit_2 - limit_1 < N0) //  otherwise, if cluster is a leaf, use direct sum
		{
            double DS_result = Call_Ds(limit_1, limit_2, particle_index,p_x,  p_y, p_z, particles, lambda);
            //printf("dir_sum: %12.6e %12.6e %12.6e \n",DS_result [0],DS_result [1],DS_result [2]);
            potential = potential + DS_result;
		}//
		else // othervise, if cluster is not a leaf, look at children
		{
			potential = 0.0;
			size_t length = tree[panel_index].children.size();
			for (size_t i = 0; i < length; i++)
			{
				size_t index = tree[panel_index].children[i];
                double temp_result = Comput_RBF(lambda, particles, particle_index, index);
                //printf("partile-cluster-child: %12.6e %12.6e %12.6e \n",temp_result [0],temp_result [1],temp_result [2]);
                potential = potential + temp_result;
			}
            
		}
        
	}
    return potential;
}

//*****************************************************************************//
void Cluster_Chev_Points(size_t tree_size)
{
    double h;
    h = 3.14159265358979323846/P;
    double t[P + 1] = {0.0};
    for (int i = 0; i < P + 1; i++)
        t[i] = cos(i * h);  //Chebyshev interpolation points [-1,1]
    
    double x1,x2,y1,y2,z1,z2;
    size_t tree_index;
    
    for (tree_index = 0; tree_index < tree_size ; tree_index++)
    {
        x1 = tree[tree_index].xinterval[0];
        x2 = tree[tree_index].xinterval[1];
        y1 = tree[tree_index].yinterval[0];
        y2 = tree[tree_index].yinterval[1];
        z1 = tree[tree_index].zinterval[0];
        z2 = tree[tree_index].zinterval[1];
 
        for (int i = 0; i < P + 1; i++) // map to the cluster
        {
            tree[tree_index].t1[i] =  x1 + (t[i] + 1)/2 * (x2 - x1);
            tree[tree_index].t2[i] =  y1 + (t[i] + 1)/2 * (y2 - y1);
            tree[tree_index].t3[i] =  z1 + (t[i] + 1)/2 * (z2 - z1);
        }
    }
}
