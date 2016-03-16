#include <stdio.h>
#include <string.h>
#include "scf.h"
#include "scf_decompose.h"
#include "mpi.h"

int top_sub;


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);

        
    decompose_nmax = 10;
    decompose_lmax = 0;
    int n_elements = decompose_nmax + 1;
    float* mat_cos = (float*)malloc(n_elements*sizeof(float));
    float* mat_cos_reduced = (float*)malloc(n_elements*sizeof(float));
    
    float center[] = {  2.223352e+01,  -2.108813e+02,  1.285330e+02 };
    float scale_radius = 3.424658e+00;
    
    float given[] = {center[0]+HOST_CENTER_X, center[1]+HOST_CENTER_Y, center[2]+HOST_CENTER_Z, scale_radius};
    
    char* in_filename = "/data/ngan/VLII/data/subhalo_np150_dir/subhalo_full_23.std";
    
    /* read_snapshot() just needs to allocate memory for and fill up xx,yy,zz (in current-halocentric
     * corrdinates), mm (mass of particles), and assign NumPart and TotalNumPart.
     * The "given" vector is laid out this way because of the format of the positions in the snapshot
     * is stored in cosmological box coordinates, not halocentric coordinates.
     */
    read_snapshot(in_filename, given);
    
    /* Do the integral */
    scf_integral_spherical(mat_cos);
    
    MPI_Barrier(MPI_COMM_WORLD);

    /* Sum up the integrals from other nodes for the cosine terms */
    MPI_Reduce(mat_cos, mat_cos_reduced, n_elements, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    char* out_filename = "./subhalo23.scf";

    FILE* fout;
    
    if(ThisTask == 0)
    {
        
        fout = fopen(out_filename, "w");

        if(fout == NULL)
        {
            perror("Can't write to matrix file!\n");
            MPI_Finalize();
            exit(1);
        }
    
        fwrite(&decompose_nmax, sizeof(decompose_nmax), 1, fout);
        fwrite(&decompose_lmax, sizeof(decompose_lmax), 1, fout);
        fwrite(mat_cos_reduced, sizeof(float), n_elements, fout);
        
        /* Note that sine terms are not stored because the coefficients for the sine terms are zero
         * in the spherical case (l = m = 0)
         */
    
        fclose(fout);
        
        printf("Wrote SCF matrix to %s\n", out_filename);
    }
    
    free(xx);
    free(yy);
    free(zz);
    free(mm);
    free(mat_cos);
    free(mat_cos_reduced);

    MPI_Finalize(); 
   
    return 0;
    
}
