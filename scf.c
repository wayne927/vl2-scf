#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scf.h"
#include "scf_decompose.h"
#include "mpi.h"

/* Used to debug the results by printing the coefficient matrices
 */
void printmat(double* mat1, double* mat2)
{
    int n,l,m;
    for(n=0; n<=decompose_nmax; n++)
    {
        for(l=0; l<=decompose_lmax; l++)
        {
            for(m=0; m<=l; m++)
            {
                printf("%d %d %d:  %e", n,l,m, mat1[ind(n,l,m)]);
                if(mat2 != NULL)
                    printf("  %e", mat2[ind(n,l,m)]);
                printf("\n");
            }
        }
    }
}

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        printf("scf n_max l_max\n");
        exit(0);
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);

    /* global variables */
    decompose_nmax = atoi(argv[1]);
    decompose_lmax = atoi(argv[2]);
    
    /* This n_elements is NOT the number of terms in the triple sum. See documentation. */
    int n_elements = (decompose_nmax+1)*(decompose_lmax+1)*(decompose_lmax+1);
	
    /* Reads the N-body snapshot. Right now this read_snapshot() routine is written for
     * VL-2 snapshots. You will need to write your own routine for any other formats.
     * Returns how many particles this computing node is handling.
     */    
//     read_snapshot(NULL, NULL);
    read_ascii();

    printf("Task %d: NumPart = %lld\n", ThisTask, NumPart);

    MPI_Barrier(MPI_COMM_WORLD);

    double* mat_cos = (double*)malloc(n_elements*sizeof(double));
    double* mat_sin = (double*)malloc(n_elements*sizeof(double));

    /* This does the integral to compute the coefficients, given the positions and masses
     * of the particles. This routine is only summing the particles that the current node
     * responsible for. The results are stored in the place holder arrays mat_cos and mat_sin,
     * which will be combined with the results from all other nodes below.
     */
    scf_integral(mat_cos, mat_sin);

    MPI_Barrier(MPI_COMM_WORLD);

    double* mat_cos_reduced = (double*)malloc(n_elements*sizeof(double));
    double* mat_sin_reduced = (double*)malloc(n_elements*sizeof(double));

    /* Sum up the integrals from other nodes for the cosine terms */
    MPI_Reduce(mat_cos, mat_cos_reduced, n_elements, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
    MPI_Barrier(MPI_COMM_WORLD);
   
    /* Sum up the integrals from other nodes for the sine terms */
    MPI_Reduce(mat_sin, mat_sin_reduced, n_elements, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  
    char out_filename[200];
    FILE* fout;
  
    /* Now write the coefficients to file */
    if(ThisTask == 0)
    {
//         printmat(mat_cos_reduced, mat_sin_reduced);

        sprintf(out_filename, "scfmatrix_%d_%d.scf", decompose_nmax, decompose_lmax);

        fout = fopen(out_filename, "w");
        
        if(fout == NULL)
        {
            perror("Can't write to matrix file!\n");
            MPI_Finalize();
            exit(1);
        }

        fwrite(&decompose_nmax, sizeof(decompose_nmax), 1, fout);
        fwrite(&decompose_lmax, sizeof(decompose_lmax), 1, fout);
        fwrite(mat_cos, sizeof(double), n_elements, fout);
        fwrite(mat_sin, sizeof(double), n_elements, fout);
        fwrite(mat_cos_reduced, sizeof(double), n_elements, fout);
        fwrite(mat_sin_reduced, sizeof(double), n_elements, fout);

        fclose(fout);

    }


    MPI_Barrier(MPI_COMM_WORLD);

	free(xx);
	free(yy);
	free(zz);
	free(mm);
	free(mat_cos);
	free(mat_sin);
	free(mat_cos_reduced);
	free(mat_sin_reduced);
	
	MPI_Finalize();
	
	return 0;
}

