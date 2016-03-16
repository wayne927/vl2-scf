#ifndef __SCF_DECOMPOSE_H__
#define __SCF_DECOMPOSE_H__

#define PI 3.141592653589793
#define ROOT_FOUR_PI 3.54490770181

#define HOST_SCALE 60.424659

#define HOST_FILENAME "/data/ngan/VLII/data/halo_0_full_host.std"

#define HOST_CENTER_X 34448.12179
#define HOST_CENTER_Y 28426.30833
#define HOST_CENTER_Z 20275.10167

#define HOST_CENTER_VX -1.24
#define HOST_CENTER_VY 34.43
#define HOST_CENTER_VZ 256.660004

#define LBOX 40000.0

#define MASS_UNIT 9.465553186e15

long long read_snapshot(char*, float*);
void read_ascii();

void scf_integral(double* mat_cos, double* mat_sin);
void scf_integral_spherical(float* mat_cos);

float* xx;
float* yy;
float* zz;
float* mm;

// int decompose_nmax;
// int decompose_lmax;

unsigned long long NumPart;
unsigned long long TotalNumPart;

int ThisTask;
int NTask;

#endif

