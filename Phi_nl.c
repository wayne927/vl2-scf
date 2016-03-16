#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include "scf.h"

/* Computing Phi_nl(r), which is the radial part of the basis functions
 *
 * The Gegenbauer polynomials for the host halo and subhalos can either be computed
 * using my own hard coded implementation (see gegenbauer.c) or GSL's implementation.
 *
 */

double int_power(float x, int n)
{
    return gsl_pow_int(x, n);
}

#ifdef USE_GSL_GEGENBAUER_HOST
float Phi_nl(double Cna, int l, float r)
{
    return -1.0f * int_power(r,l)/int_power(1.0f+r, 2.0*l+1) * (float)Cna * ROOT_FOUR_PI;
}

float dPhi_nl_dr(double Cna, float phinl, int n, int l, float r)
{
    int two_l = 2*l;
	float r1 = r+1.f;
	
	float val = phinl * ( l/r - (two_l+1.0f)/r1 );
	
	if(n>0)
	{
		val -= ROOT_FOUR_PI * int_power(r,l)*4.0f*(two_l+1.5f) * (float)Cna / int_power(r1,two_l+3);
	}
	
	return val;
}

#else
float Phi_nl(int n, int l, float r)
{
    float xi = (r-1.0f)/(r+1.0f);
    
    return -1.0f * int_power(r,l)/int_power(1.0f+r, 2.0*l+1) * Cna(n,2.0f*l+1.5f,xi) * ROOT_FOUR_PI;
}

float dPhi_nl_dr(float phinl, int n, int l, float r)
{
    int two_l = 2*l;
	float r1 = r+1.f;
	
	float val = phinl * ( l/r - (two_l+1.0f)/r1 );
	
	if(n>0)
	{
		float xi = (r-1.f)/r1;
		
		val -= ROOT_FOUR_PI * int_power(r,l)*4.0f*(two_l+1.5f) * Cna(n-1,two_l+2.5f,xi) / int_power(r1,two_l+3);
	}
	
	return val;
}
#endif


#ifdef USE_GSL_GEGENBAUER_SUBHALO
float Phi_n0(double Cna, float r) 
{
    return -1.0f / (1.0f+r) * (float)Cna * ROOT_FOUR_PI;
}

float dPhi_n0_dr(double Cna, float phin0, int n, float r)
{
    float r1 = r+1.f;
    
    float val = phin0 * -1.f/r1;
    
    if(n>0)
    {
        val -= ROOT_FOUR_PI * 6.f/r1/r1/r1 * (float)Cna;
    }
    
    return val;
}

#else
float Phi_n0(int n, float r)
{
    float xi = (r-1.0f)/(r+1.0f);
    
    return -1.0f / (1.0f+r) * Cna(n,1.5f,xi) * ROOT_FOUR_PI;
}

float dPhi_n0_dr(float phin0, int n, float r)
{
    float r1 = r+1.f;
    
    float val = phin0 * -1.f/r1;
    
    if(n>0)
    {
		float xi = (r-1.f)/r1;
        val -= ROOT_FOUR_PI * 6.f/r1/r1/r1 * Cna(n-1,2.5f,xi);
    }
    
    return val;
}

#endif


