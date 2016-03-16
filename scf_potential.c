#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include "scf.h"


/* Given the coefficients stored mat_cos and mat_sin, this code does
 * the summation to compute the potential and accelerations
 *
 * scf_potential() is the full 3D treatment for the host halo
 *
 * scf_potential_spherical() is the spherically symmetric treatment
 * for the subhalos. This spherical code skips the loops over l,m
 * so it's simpler.
 */

float scf_potential(float* a, double* mat_cos, double* mat_sin, int orbit_nmax, int orbit_lmax, float x, float y, float z)
{
    double x2 = (double)x * (double)x;
    double y2 = (double)y * (double)y;
    double r = sqrt(x2+y2+(double)z*(double)z);
    double r_xy = sqrt(x2+y2);
    double phi = atan2((double)y, (double)x);

    double sin_theta = r_xy / r;
    double sin_phi = (double)y / r_xy;
    double cos_theta = (double)z / r;
    double cos_phi = (double)x / r_xy;
    
    double cos_theta_for_Plm = cos_theta;
    
#ifdef USE_GSL_GEGENBAUER_HOST
    double xi = (r-1.0)/(r+1.0);
#endif


    double a_r = 0.0;
    double a_t = 0.0;
    double a_p = 0.0;

    int n,l,m;

    float C_lm = 0.0f;
    float D_lm = 0.0f;
    float E_lm = 0.0f;
    float F_lm = 0.0f;

    float this_Phi_nl = 0.f;
    float this_dPhi_nl_dr = 0.f;

#ifdef USE_GSL_GEGENBAUER_HOST
    double* Cna_2l_32_array = (double*)malloc((orbit_nmax+1)*sizeof(double));
    double* Cna_2l_52_array = (double*)malloc(orbit_nmax*sizeof(double));
#endif
    
    double this_P_lm = 0.0;
    double this_dP_lm = 0.0;
    double* P_l;
    double* dP_l;
    int P_size;
    
    float cos_mphi;
    float sin_mphi;

    float pot = 0.f;
    a[0] = 0.f;
    a[1] = 0.f;
    a[2] = 0.f;
	
    for(m=0; m<=orbit_lmax; m++)
    {
        cos_mphi = cos(m*phi);
        sin_mphi = sin(m*phi);
        
        P_size = orbit_lmax - m + 1;
        P_l = (double*)malloc(P_size*sizeof(double));
        dP_l = (double*)malloc(P_size*sizeof(double));
        
        gsl_sf_legendre_Plm_deriv_array(orbit_lmax, m, cos_theta_for_Plm, P_l, dP_l);
                
        for(l=m; l<=orbit_lmax; l++)
        {
            C_lm = 0.f;
            D_lm = 0.f;
            E_lm = 0.f;
            F_lm = 0.f;

#ifdef USE_GSL_GEGENBAUER_HOST
            gsl_sf_gegenpoly_array(orbit_nmax, 2*l+1.5, xi, Cna_2l_32_array);
            if(orbit_nmax > 0)
                gsl_sf_gegenpoly_array(orbit_nmax-1, 2*l+2.5, xi, Cna_2l_52_array);
#endif
   
            for(n=0; n<=orbit_nmax; n++)
            {
#ifdef USE_GSL_GEGENBAUER_HOST
                this_Phi_nl = Phi_nl(Cna_2l_32_array[n], l, r);
                if(n>0)
                    this_dPhi_nl_dr = dPhi_nl_dr(Cna_2l_52_array[n-1], this_Phi_nl, n, l, r);
                else
                    this_dPhi_nl_dr = dPhi_nl_dr(1e99, this_Phi_nl, n, l, r);
#else           
                this_Phi_nl = Phi_nl(n,l,r);
                this_dPhi_nl_dr = dPhi_nl_dr(this_Phi_nl, n, l, r);
#endif
                
                C_lm += mat_cos[ind(n,l,m)] * this_Phi_nl;
                D_lm += mat_sin[ind(n,l,m)] * this_Phi_nl;
       
                E_lm += mat_cos[ind(n,l,m)] * this_dPhi_nl_dr;
                F_lm += mat_sin[ind(n,l,m)] * this_dPhi_nl_dr;
            }
   
            this_P_lm = P_l[l-m];
            this_dP_lm = dP_l[l-m] * -sin_theta;
               
            pot += this_P_lm * (C_lm*cos_mphi + D_lm*sin_mphi);

            a_r += this_P_lm * (E_lm*cos_mphi + F_lm*sin_mphi);
            a_t += this_dP_lm * (C_lm*cos_mphi + D_lm*sin_mphi);
            a_p += m * this_P_lm * (D_lm*cos_mphi - C_lm*sin_mphi);

        } // l loop
        
        free(P_l);
        free(dP_l);

    } // m loop

    a_r = -a_r;
    a_t = -a_t/r;
    a_p = -a_p/r_xy;
    
    a[0] = a_r*sin_theta*cos_phi + a_t*cos_theta*cos_phi - a_p*sin_phi;
    a[1] = a_r*sin_theta*sin_phi + a_t*cos_theta*sin_phi + a_p*cos_phi;
    a[2] = a_r*cos_theta - a_t*sin_theta;

    a[0] *= G;
    a[1] *= G;
    a[2] *= G;
    pot *= G;
    
#ifdef USE_GSL_GEGENBAUER_HOST
    free(Cna_2l_32_array);
    free(Cna_2l_52_array);
#endif

    return pot;

}


float scf_potential_spherical(float* a, float* mat_cos, int orbit_nmax, float x, float y, float z)
{
    float x2 = x*x;
    float y2 = y*y;
    
    float r = sqrt(x2+y2+z*z);
    float r_xy = sqrt(x2+y2);
    
#ifdef USE_GSL_GEGENBAUER_SUBHALO
    float xi = (r-1.0f)/(r+1.0f);
#endif

    float sin_theta = r_xy / r;
    float sin_phi = y / r_xy;
    float cos_theta = z / r;
    float cos_phi = x / r_xy;


    float this_Phi_n0 = 0.f;
    float this_dPhi_n0_dr = 0.f;
    
#ifdef USE_GSL_GEGENBAUER_SUBHALO
    double* Cna_32_array = (double*)malloc((orbit_nmax+1)*sizeof(double));
    double* Cna_52_array = (double*)malloc(orbit_nmax*sizeof(double));
    gsl_sf_gegenpoly_array(orbit_nmax, 1.5, xi, Cna_32_array);
    gsl_sf_gegenpoly_array(orbit_nmax-1, 2.5, xi, Cna_52_array);
#endif

    float C_lm = 0.0f;
    float E_lm = 0.0f;

    int n;
    for(n=0; n<=orbit_nmax; n++)
    {
#ifdef USE_GSL_GEGENBAUER_SUBHALO
        this_Phi_n0 = Phi_n0(Cna_32_array[n], r);
        if(n>0)
            this_dPhi_n0_dr = dPhi_n0_dr(Cna_52_array[n-1], this_Phi_n0, n, r);
        else
            this_dPhi_n0_dr = dPhi_n0_dr(1e99, this_Phi_n0, n, r);
#else
        this_Phi_n0 = Phi_n0(n,r);
        this_dPhi_n0_dr = dPhi_n0_dr(this_Phi_n0, n, r);
#endif

        C_lm += mat_cos[n] * this_Phi_n0;

        E_lm += mat_cos[n] * this_dPhi_n0_dr;
    }

    float pot = C_lm;

    float a_r = -E_lm;

    a[0] = G * a_r*sin_theta*cos_phi;
    a[1] = G * a_r*sin_theta*sin_phi;
    a[2] = G * a_r*cos_theta;
    
    pot = G * pot;

#ifdef USE_GSL_GEGENBAUER_SUBHALO
    free(Cna_32_array);
    free(Cna_52_array);
#endif

    return pot;

}
