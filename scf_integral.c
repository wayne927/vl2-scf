#include "scf.h"
#include "scf_decompose.h"
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gegenbauer.h>

/* computes a!/b! =  a(a-1)(a-2)...(b+1) */
double factorial_frac(int a, int b) 
{
    if(a<0 || b<0)
    {
        perror("factorial_frac: either a or b is negative!!\n");
        exit(1);
    }

    int recip = 0;
    if(a==b)
        return 1.0;
    else if(a<b)
    {
        int temp = b;
        b = a;
        a = temp;
        recip = 1;
    }

    double val = 1;
    
    while(a > b)
    {
        val = val * a;
        a--;
    }
    
    if(recip)
        val = 1.0/val;
    
    return val;
    
}

void scf_integral(double* mat_cos, double* mat_sin)
{	
    long long k;
    double* rr = (double*)malloc(NumPart*sizeof(double));
    double* cos_tt = (double*)malloc(NumPart*sizeof(double));
    double* pp = (double*)malloc(NumPart*sizeof(double));

#ifdef USE_GSL_GEGENBAUER_HOST
    double* xxi = (double*)malloc(NumPart*sizeof(double));
    double this_Cna;
#endif

    /* cos_tt is calculated directly by cos(theta)=z/r because
    we never need theta by itself, but only cos(theta) in the 
    legendre polynomials. Also, notice that cos_tt must be
    calculated in double precision because with single precision
    zz[k] can be greater than rr[k] by a tiny bit, which 
    the plgndr routine will complain later on */

    for(k=0; k<NumPart; k++)
    {
        rr[k] = sqrt((double)xx[k]*(double)xx[k] +
                      (double)yy[k]*(double)yy[k] +
                      (double)zz[k]*(double)zz[k]);

        cos_tt[k] = (double)zz[k]/rr[k];

        pp[k] = atan2(yy[k], xx[k]);
        
#ifdef USE_GSL_GEGENBAUER_HOST
        xxi[k] = (rr[k]-1.0) / (rr[k]+1.0);
#endif
    }

    double N_lm;
    double K_nl, I_nl, A_nl;
    double gamm1, gamm2;
   
    float this_Phi_nl;
    double* P_lm_k = (double*)malloc(NumPart*sizeof(double));
   
    double m_Phinlk_Plm, Nlm_Anl;
    double sum_cos = 0.0;
    double sum_sin = 0.0;
   
    int n, l, m;
    for(l=0; l<=decompose_lmax; l++)
    {
        for(m=0; m<=l; m++)
        {
            N_lm = (2*l+1)/4.0/PI * 2.0 * factorial_frac(l-m,l+m);
            if(m==0)
                N_lm = 0.5*N_lm;

            for (k=0; k<NumPart; k++)
            {
                P_lm_k[k] = gsl_sf_legendre_Plm(l,m, cos_tt[k]);
                
                /* Note: This can be computed more efficiently using gsl_sf_legendre_Plm_array
                which computes a list of legendre polynomials at constant m with varying l.
                but this changes the code and I don't have time to deal with that for now.
                Numerical Recipes's plgndr is actually faster, but single precision. */
            }
               
            for(n=0; n<=decompose_nmax; n++)
            {
                K_nl = 0.5*n*(n+4*l+3)+(l+1)*(2*l+1);
                gamm1 = factorial_frac(n+4*l+2,n);
                gamm2 = factorial_frac(2*l+1, 2*(2*l+1));
                I_nl = -K_nl/(n+2*l+1.5) * gamm1 * gamm2 * gamm2;
                A_nl = 1.0/I_nl;
               
                sum_cos = 0.0;
                sum_sin = 0.0;
               
                for(k=0; k<NumPart; k++)
                {
#ifdef USE_GSL_GEGENBAUER_HOST
                    this_Cna = gsl_sf_gegenpoly_n(n, 2*l+1.5, xxi[k]);
                    this_Phi_nl = Phi_nl(this_Cna, l, rr[k]); 
#else
                    this_Phi_nl = Phi_nl(n, l, rr[k]);
#endif
                   
                    m_Phinlk_Plm = mm[k] * this_Phi_nl * P_lm_k[k];
                   
                    sum_cos += m_Phinlk_Plm * cos(m*pp[k]);
                    sum_sin += m_Phinlk_Plm * sin(m*pp[k]);
                }
               
                Nlm_Anl = N_lm * A_nl;
               
                mat_cos[ind(n,l,m)] = sum_cos * Nlm_Anl;
                mat_sin[ind(n,l,m)] = sum_sin * Nlm_Anl;


            } // n loop
       
        } // m loop
       
    } // l loop
   
    free(rr);
    free(pp);
    free(cos_tt);
    free(P_lm_k);

}

void scf_integral_spherical(float* mat_cos)
{
    long long k;
    float* rr = (float*)malloc(NumPart*sizeof(float));
    double temp_r;

    for(k=0; k<NumPart; k++)
    {
        temp_r = sqrt(xx[k]*xx[k] + yy[k]*yy[k] + zz[k]*zz[k]);
        rr[k] = (float)temp_r;
    }

    double N00 = 0.07957747154; /* = 1/(4*pi) */
    double K_n0, I_n0, A_n0;
    
    float this_Phi_n0;

    double sum_cos = 0.0;
   
    int n;
               
    for(n=0; n<=decompose_nmax; n++)
    {
        K_n0 = 0.5*n*(n+3)+1;
        I_n0 = -K_n0/(n+1.5)/4.0 * (n+2) * (n+1);
        A_n0 = 1.0/I_n0;
       
        sum_cos = 0.0;
       
        for(k=0; k<NumPart; k++)
        {
            this_Phi_n0 = Phi_n0(n, rr[k]);
            sum_cos += mm[k] * this_Phi_n0;
        }
       
        mat_cos[n] = (float) (sum_cos * A_n0 * N00);

    } // n loop
       
   
    free(rr);

}

