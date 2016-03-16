#include "scf.h"

/* This function converts the triple index (n,l,m) into a single index in a C array.
 * Note that the way the coefficients were stored in the file and then read by the 
 * drive code made the entire matrix (nmax+1)*(lmax+1)^2 in size, but the triple sum
 * in the SCF method only requires (nmax+1)*(lmax+1)*(lmax+2)/2 elements. The
 * difference is wasted memory, but it makes things simpler and it's only a few
 * thousand numbers which is negligible compared to the entire SCF method.
 */

int ind(int n, int l, int m)
{
	return (m*(decompose_lmax+1) + l)*(decompose_nmax+1) + n;
}
