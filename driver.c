#include <stdio.h>
#include <string.h>
#include <math.h>
#include "scf.h"

#define HOST_SCALE 60.424659
#define SCF_FILENAME "./data/host_smooth_scf30_30.scf"
#define SUBHALO_CATALOG "./data/subhalo_catalog.dat"

int Nsub;

typedef struct {
    
    /* SUBHALO_DECOMPOSE_ORDER is n_max for subhalos. l_max is always 0. 
     * If subhalos were not spherical, the number of elements in this array would be
     * (n_max+1)*(l_max+1)^2
     */
    float scf_mat[SUBHALO_DECOMPOSE_ORDER+1];
    
    int id;
    float pos_x;
    float pos_y;
    float pos_z;
    float vel_x;
    float vel_y;
    float vel_z;
    float scale;
    float rvir;
    float mass;
} Sub;

Sub* subhalos;

float get_accel_host(float* a, double* mat_cos, double* mat_sin, int orbit_nmax, int orbit_lmax, float x, float y, float z)
{
    /*       a[3] = pointer to store the resulting accelerations
     *    mat_cos = the array of coefficients for the cosine terms, read from the SCF matrix
     *    mat_sin = the array of coefficients for the sine terms, read from the SCF matrix
     * orbit_nmax = how many n orders we are using for our orbit
     * orbit_lmax = how many l orders we are using for our orbit
     *      x,y,z = physical position where we want to compute the accelerations
     */
    
    a[0] = 0.0f;
    a[1] = 0.0f;
    a[2] = 0.0f;
    
    float pot = scf_potential(a, mat_cos, mat_sin, orbit_nmax, orbit_lmax, x/HOST_SCALE, y/HOST_SCALE, z/HOST_SCALE);
    
    float s2 = HOST_SCALE * HOST_SCALE;
    a[0] /= s2;
    a[1] /= s2;
    a[2] /= s2;

    return pot/HOST_SCALE;
    
}

float get_accel_subhalos(float* a, int orbit_nmax, float x, float y, float z)
{
    /* The arguments are similar to get_accel_host, except orbit_lmax is always zero, and the
     * SCF coefficients of each subhalo are already stored inside the global Sub data structures.
     */
     
    if(Nsub == 0)
        return 0;
    
    /* total potential and accelation of the given particle due to ALL subhalos */
    float pot = 0.0f;
    a[0] = 0.0f;
    a[1] = 0.0f;
    a[2] = 0.0f;
    
    /* potential and acceleration of the given particle due to each subhalo */
    float subhalo_pot;
    float subhalo_a[3];
    
    float this_scale, this_scale2;
    
    float dx, dy, dz;
    
    int i;
    for(i=0; i<Nsub; i++)
    {
        this_scale = subhalos[i].scale;
        this_scale2 = this_scale*this_scale;
        dx = x - subhalos[i].pos_x;
        dy = y - subhalos[i].pos_y;
        dz = z - subhalos[i].pos_z;
        
        /* Note that for subhalos, we use the spherical version of scf_potential() which takes less arguments and
         * skips all the loops over the angular parts
         */
        subhalo_pot = scf_potential_spherical(subhalo_a, subhalos[i].scf_mat, orbit_nmax, dx/this_scale, dy/this_scale, dz/this_scale);
        subhalo_pot /= this_scale; 
        subhalo_a[0] /= this_scale2;
        subhalo_a[1] /= this_scale2;
        subhalo_a[2] /= this_scale2;
        
        /* To speed things up, you can approximate the subhalo accelerations like this. For simplicity
         * I'm not doing it here.
         * 
         *  if(distance from the subhalo is far enough)
         *  {
         *      treat the subhalo as a point source;
         *      subhalo_a = GM/r^2;
         *  }
         *  else
         *  {
         *      use scf_potential_spherical()
         *  }
         */         
        
        pot += subhalo_pot;
        a[0] += subhalo_a[0];
        a[1] += subhalo_a[1];
        a[2] += subhalo_a[2];
    }
    
    return pot;
}

void load_subhalo_scf()
{
    if(Nsub == 0)
        return;
    
    subhalos = (Sub*)malloc(Nsub*sizeof(Sub));

    FILE* fin_info = fopen(SUBHALO_CATALOG, "r");
    
    int dummy;
    
    FILE* fin_subhalo;
    char filename_subhalo[200];
    
    int i;
    for(i=0; i<Nsub; i++)
    {
        /* Read the subhalo catalogue */
        fscanf(fin_info, "%d%e%e%e%e%e%e%e%e%e", &subhalos[i].id, &subhalos[i].mass,
            &subhalos[i].pos_x, &subhalos[i].pos_y, &subhalos[i].pos_z,
            &subhalos[i].vel_x, &subhalos[i].vel_y, &subhalos[i].vel_z,
            &subhalos[i].rvir, &subhalos[i].scale);
        
        /* Take the subhalo's ID, now read the subhalo's SCF file */
        memset(filename_subhalo, '\0', 200);
        sprintf(filename_subhalo, "./data/subhalo_scf/subhalo_%d_scf10_0.scf", subhalos[i].id);
       
        fin_subhalo = fopen(filename_subhalo, "r");
        if(fin_subhalo == NULL)
        {
            perror("Error opening subhalo SCF file!!\n");
            exit(1);
        }
        
        /* Similar to the host SCF file, the first two integers of each subhalo SCF file are also
         * the decompose_nmax and decompose_lmax. In our papers though, we didn't actually care that much
         * about subhalo's decomposition orders, so here I just hard code them. Each subhalo's SCF file has
         * decompose_nmax = 10
         * decompose_lmax = 0
         */
        fread(&dummy, sizeof(dummy), 1, fin_subhalo);
        fread(&dummy, sizeof(dummy), 1, fin_subhalo);
        
        /* Without the angular parts, the radial decomposition has decompose_nmax + 1 = 11 terms */
        fread(subhalos[i].scf_mat, sizeof(float), SUBHALO_DECOMPOSE_ORDER+1, fin_subhalo);
        
        fclose(fin_subhalo);
        
    }
    
    fclose(fin_info);

}

int main(int argc, char* argv[])
{
    /* orbit_nmax and orbit_lmax are how many orders we want to use for our orbits.
     * In Ngan et al (2015), orbit_nmax = orbit_lmax = 10
     */
    int orbit_nmax = atoi(argv[1]);
    int orbit_lmax = atoi(argv[2]);
    
    Nsub = 0;
    if(argc > 3)
    {
        Nsub = atoi(argv[3]);
    }
    
    FILE* fin = fopen(SCF_FILENAME, "r");
    
    if(fin == NULL)
    {
        perror("Error opening the SCF matrix file!\n");
        exit(1);
    }
    
    /* decompose_nmax = how many orders are stored inside the matrix file (=30 in the included host_smooth_scf30_30.scf file)
     *     orbit_nmax = how many orders we want to use for our orbit, specified as a command line argument
     *
     * (orbit_nmax must be less than decompose_nmax)
     *
     * Same thing for lmax.
     */
    
    fread(&decompose_nmax, sizeof(int), 1, fin);
    fread(&decompose_lmax, sizeof(int), 1, fin);
    
    /* This n_elements is NOT the number of terms in the triple sum. See documentation. */
    int n_elements = (decompose_nmax+1)*(decompose_lmax+1)*(decompose_lmax+1);

    double* mat_cos = (double*)malloc(n_elements*sizeof(double));
    double* mat_sin = (double*)malloc(n_elements*sizeof(double));
    
    /* These are the complex coefficients, divided for cosine and sine terms as outlined in Hernquist & Ostriker (1992) */
    fread(mat_cos, sizeof(double), n_elements, fin);
    fread(mat_sin, sizeof(double), n_elements, fin);
    
    fclose(fin);
    
    load_subhalo_scf();
    
    /* NP = 1 for just one particle */
#define NP 1
    float xx[NP];
    float yy[NP];
    float zz[NP];
    float vvx[NP];
    float vvy[NP];
    float vvz[NP];
    float EE[NP];
    
    /* initial position in kpc */    
    xx[0] = 30.0;
    yy[0] = 0;
    zz[0] = 0;
 
    /* initial x = (30,0,0) kpc, v = (0,120,0) km/s is "orbit 1" in Ngan et al (2015) */
    vvx[0] = 0.0;
    vvy[0] = 120.0;
    vvz[0] = 0.0;
    
    /* initial x = (30,0,0) kpc, v = (0,41,113) km/s is "orbit 2" in Ngan et al (2015) */
//     vvx[0] = 0.0;
//     vvy[0] = 41.0;
//     vvz[0] = 113.0;
    
    
    float acc[3] = {0.0f, 0.0f, 0.0f};
    float acc_subhalos[3] = {0.0f, 0.0f, 0.0f};
    float potential;
    
    
    /* time goes from 0 to 10 Gyr, in steps of dt=0.005 Gyr */
    float t = 0.0;
    float dt = 0.001;
    
    int p, s;
    
    /* Note: This Euler integrator is meant for simplicity in this demo, NOT FOR ACCURACY!! */
    for(t=0.0; t<10.0; t+=dt)
    {
        for(p=0; p<NP; p++)
        {
            acc[0] = 0;
            acc[1] = 0;
            acc[2] = 0;
            potential = get_accel_host(acc, mat_cos, mat_sin, orbit_nmax, orbit_lmax, xx[p], yy[p], zz[p]);
            
            if(Nsub > 0)
            {
                potential = potential + get_accel_subhalos(acc_subhalos, SUBHALO_ORBIT_ORDER, xx[p], yy[p], zz[p]);
                acc[0] += acc_subhalos[0];
                acc[1] += acc_subhalos[1];
                acc[2] += acc_subhalos[2];
            }
            
            vvx[p] = vvx[p] + acc[0]*dt;
            vvy[p] = vvy[p] + acc[1]*dt;
            vvz[p] = vvz[p] + acc[2]*dt;
        
            xx[p] = xx[p] + vvx[p]*dt;
            yy[p] = yy[p] + vvy[p]*dt;
            zz[p] = zz[p] + vvz[p]*dt;

            EE[p] = potential + 0.5*(vvx[p]*vvx[p] + vvy[p]*vvy[p] + vvz[p]*vvz[p]);
        }
        
        /* print the particle's position and total energy */
        printf("%e  %e  %e  %e  %e\n", t, xx[0], yy[0], zz[0], EE[0]);

        
        /* move the subhalos too! reuse the acc array, now for subhalos */
        for(s=0; s<Nsub; s++)
        {
            acc[0] = 0;
            acc[1] = 0;
            acc[2] = 0;
            get_accel_host(acc, mat_cos, mat_sin, orbit_nmax, orbit_lmax, subhalos[s].pos_x, subhalos[s].pos_y, subhalos[s].pos_z);
                        
            subhalos[s].vel_x += acc[0]*dt;
            subhalos[s].vel_y += acc[1]*dt;
            subhalos[s].vel_z += acc[2]*dt;
            
            subhalos[s].pos_x += subhalos[s].vel_x*dt;
            subhalos[s].pos_y += subhalos[s].vel_y*dt;
            subhalos[s].pos_z += subhalos[s].vel_z*dt;
        }
        
        /* print the subhalo orbits too, if you like */
//         printf("%e  %e  %e  %e\n", t, subhalos[0].pos_x, subhalos[0].pos_y, subhalos[0].pos_z);
        
        

    }
    
    free(mat_cos);
    free(mat_sin);
    
    return 0;
}


