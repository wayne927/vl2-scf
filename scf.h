#ifndef __SCF_H__
#define __SCF_H__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* For subhalos: they are low enough order that the hard coded single precison
 * Cna that I wrote are faster than GSL, so don't use GSL
 */

/* For host: Host needs higher orders, and it's hard to hard code the high
 * orders (only up to n=7 for now). So in general, GSL is easier to use and
 * faster. But if my own Cna can be hard coded to sufficiently high order, (instead
 * of relying on recursion), then it's *probably* faster than GSL, but needs
 * to be tested. So, use GSL for host.
 */

// #define USE_GSL_GEGENBAUER_SUBHALO
#define USE_GSL_GEGENBAUER_HOST

#define PI 3.141592653589793
#define ROOT_FOUR_PI 3.54490770181

#define SUBHALO_DECOMPOSE_ORDER 10
#define SUBHALO_ORBIT_ORDER 4

#define G 43016.0e-10

int decompose_nmax;
int decompose_lmax;

int ind(int n, int l, int m);

float Cna(int n, float alpha, float xi);

/* basis functions */
#ifdef USE_GSL_GEGENBAUER_HOST
float Phi_nl(double Cna, int l, float r);
float dPhi_nl_dr(double Cna, float phinl, int n, int l, float r);
#else
float Phi_nl(int n, int l, float r);
float dPhi_nl_dr(float phinl, int n, int l, float r);
#endif

/* spherical routines */
#ifdef USE_GSL_GEGENBAUER_SUBHALO
float Phi_n0(double Cna, float r);
float dPhi_n0_dr(double Cna, float phin0, int n, float r);
#else
float Phi_n0(int n, float r);
float dPhi_n0_dr(float phin0, int n, float r);
#endif

/* These actually do the sum to compute the potential and acceleration */
float scf_potential(float* a, double* mat_cos, double* mat_sin, int orbit_nmax, int orbit_lmax, float x, float y, float z);
float scf_potential_spherical(float* a, float* mat_cos, int orbit_nmax, float x, float y, float z);

#endif

