#include <stdlib.h>

/* Compute Gegenbauer polynomials C_na
 * When n<=7, the polynomials are hard coded.
 * When n>7, it uses the recursion relation (implemented iteratively here).
 *
 * Note that in the scf.h header file you can choose whether to use this
 * code to compute the polynomials, or GSL's implementation.
 * #define USE_GSL_GEGENBAUER_SUBHALO
 * and
 * #define USE_GSL_GEGENBAUER_HOST
 * 
 * My experience is that this code is faster, but GSL is more accurate.
 * I use this code to compute forces of subhalos, and GSL to compute forces
 * of the host. This is why by default USE_GSL_GEGENBAUER_SUBHALO is off, but 
 * USE_GSL_GEGENBAUER_HOST is on in scf.h
 *
 */


float Cna(int n, float a, float x)
{
    switch(n)
    {
        case 0:
        return 1;
        break;
        
        case 1:
        return 2*a*x;
        break;
        
        case 2:
        return -a + 2*a*(1+a)*x*x;
        break;
        
        case 3:
        {
            float x3 = x*x*x;
            return -2*a*(1+a)*x + 1.333333333333f*a*(1+a)*(2+a)*x3;
            break;
        }
        
        case 4:
        {
            float x2 = x*x;
            float x4 = x2*x2;
            return 0.1666666666667f*a*(1+a) * ( 3 - 12*(2+a)*x2 + 4*(6+5*a+a*a)*x4 );
            break;
        }
        
        case 5:
        {
            float x2 = x*x;
            float x4 = x2*x2;;
            return 0.06666666667f*a*(1+a)*(2+a) * x * ( 15 - 20*(3+a)*x2 + 4*(12+7*a+a*a)*x4 );
            // 1/15 = 0.066666667
            
            break;
        }
        
        case 6:
        {
            float x2 = x*x;
            float x4 = x2*x2;
            float x6 = x4*x2;
            float a2 = a*a;
            float a3 = a2*a;
            return 0.01111111111f*a*(1+a)*(2+a) * ( -15 + 90*(3+a)*x2 - 60*(12+7*a+a2)*x4 + 8*(60+47*a+12*a2+a3)*x6 );
            // 1/90 = 0.01111111
            
            break;
        }
        case 7:
        {
            float x2 = x*x;
            float x4 = x2*x2;
            float x6 = x4*x2;
            float a2 = a*a;
            float a3 = a2*a;
            return 0.00317460317f*a*(1+a)*(2+a) * x * 
                ( -3*(5+2*a)*( 15-20*(3+a)*x2 + 4*(12+7*a+a2)*x4 ) + 
                (6+a)*( -15 + 90*(3+a)*x2 - 60*(12+7*a+a2)*x4 + 8*(60+47*a+12*a2+a3)*x6 ) );
            // 1/315 = 0.00317460317
            
            break;
        }
        
    }
    
    int highest_n = 7;

    float* C_array = (float*)malloc((n+1)*sizeof(float));
    
    C_array[highest_n-1] = Cna(highest_n-1,a,x);
    C_array[highest_n]   = Cna(highest_n,a,x);
    
    int i;
    for(i=highest_n+1; i<=n; i++)
    {
        C_array[i] = 2.f*(i+a-1.f)/i * x * C_array[i-1] - (i+2.f*a-2.f)/i * C_array[i-2];
    }

    float val = C_array[n];
    free(C_array);

    return val;

}

