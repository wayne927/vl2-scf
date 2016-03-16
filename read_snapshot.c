#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "mpi.h"
#include "scf_decompose.h"

typedef struct {
    double dTime;
    unsigned nBodies;
    unsigned nDim;
    unsigned nSph;
    unsigned nDark;
    unsigned nStar;
    unsigned nPad;
} tipsyHdr;

typedef struct {
    FILE *fp;
    int bEncoding;
} XDR;


//#define XDR_DECODE 2

int xdr_double(XDR* xdr, double *d)
{
//     uint64_t v;
    unsigned char c[8];
    union {
        double   d;
        uint64_t v;
    } u;
    
    if (fread(c,sizeof(c),1,xdr->fp)!=1) return 0;
    u.v = 0;   u.v |= c[0];
    u.v <<= 8; u.v |= c[1];
    u.v <<= 8; u.v |= c[2];
    u.v <<= 8; u.v |= c[3];
    u.v <<= 8; u.v |= c[4];
    u.v <<= 8; u.v |= c[5];
    u.v <<= 8; u.v |= c[6];
    u.v <<= 8; u.v |= c[7];
    *d = u.d; 
    return 1;
}

int xdr_float(XDR *xdr, float *f) {
//     uint32_t v;
    unsigned char c[4];
    union {
        float    f;
        uint32_t v;
    } u;

    if (fread(c,sizeof(c),1,xdr->fp)!=1) return 0;
    u.v = (c[0]<<24) | (c[1]<<16) | (c[2]<<8) | c[3];
    *f = u.f;
    return 1;

}

int xdr_u_int(XDR *xdr, uint32_t *u) {
    uint32_t v;
    unsigned char c[4];

    if (fread(c,sizeof(c),1,xdr->fp)!=1) return 0;
    v = (c[0]<<24) | (c[1]<<16) | (c[2]<<8) | c[3];
    *u = v;
    return 1;

}


int xdrHeader(XDR *pxdr,tipsyHdr *ph) {
    if (!xdr_double(pxdr,&ph->dTime)) return 0;
    if (!xdr_u_int(pxdr,&ph->nBodies)) return 0;
    if (!xdr_u_int(pxdr,&ph->nDim)) return 0;
    if (!xdr_u_int(pxdr,&ph->nSph)) return 0;
    if (!xdr_u_int(pxdr,&ph->nDark)) return 0;
    if (!xdr_u_int(pxdr,&ph->nStar)) return 0;
    if (!xdr_u_int(pxdr,&ph->nPad)) return 0;
    return 1;
}

void tipsySussHeader(tipsyHdr *h,uint64_t *pN, uint64_t *pDark,uint64_t *pSph, uint64_t *pStar) {
    *pN = h->nPad & 0x000000ff;
    *pN <<= 32;
    *pN += h->nBodies;

    *pSph = h->nPad & 0x0000ff00;
    *pSph <<= 24;
    *pSph += h->nSph;

    *pDark = h->nPad & 0x00ff0000;
    *pDark <<= 16;
    *pDark += h->nDark;

    *pStar = h->nPad & 0xff000000;
    *pStar <<= 8;
    *pStar += h->nStar;

    /* There may be junk in the pad field, in which case we ignore it */
    if ( *pN != *pDark + *pSph + *pStar ) {
	*pN = h->nBodies;
	*pDark = h->nDark;
	*pSph  = h->nSph;
	*pStar = h->nStar;
	}
}

double norm2_double(double r[])
{
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}
float norm2_float(float r[])
{
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}


long long read_snapshot(char* filename_given, float* scales_given)
{
    float cen[3];
    float scale;
    char* filename;
    
    if(scales_given == NULL)
    {
        cen[0] = HOST_CENTER_X;
        cen[1] = HOST_CENTER_Y;
        cen[2] = HOST_CENTER_Z;
        scale = HOST_SCALE;
    }
    else
    {
        cen[0] = scales_given[0];
        cen[1] = scales_given[1];
        cen[2] = scales_given[2];
        scale = scales_given[3];
    }

    if(filename_given == NULL)
        filename = HOST_FILENAME;
    else
        filename = filename_given;

    FILE* fin = fopen(filename, "r");
    if(fin == NULL)
    {
        perror("read_snapshot: Snapshot file not found!!\n");
        MPI_Finalize();
        exit(1);
    }
    
    uint64_t N;
    uint64_t nDark;
    uint64_t nSph;
    uint64_t nStar;
   
    tipsyHdr h;
    fread(&h, sizeof(h), 1, fin);
    
    rewind(fin);
    XDR xdr;
    xdr.bEncoding = 0;
    xdr.fp = fin;
    xdrHeader(&xdr, &h);
    
    tipsySussHeader(&h, &N, &nDark, &nSph, &nStar);
    
    TotalNumPart = N;
    
    unsigned long long head_pos = ftell(fin);
    
    fseek(fin, 0L, SEEK_END);
    unsigned long long filesize = ftell(fin);

    int doublePos = 0;
    int doubleVel = 0;
    
    // file size is header + mass,soft,pot + positions + velocity
    if(head_pos + sizeof(float)*N*3 + sizeof(float)*3*N + sizeof(double)*3*N == filesize)
        doublePos = 1;
    else if(head_pos + sizeof(float)*N*3 + sizeof(double)*N*6 == filesize)
    {
        doublePos = 1;
        doubleVel = 1;
    }
    
    unsigned long long particle_size = 36; // mass+soft+pot+pos+vel
    if(doublePos) particle_size += 12;
    if(doubleVel) particle_size += 12;
    
    float mass; //, soft, pot;
    float r_float[3]; //, v_float[3];
    double r_double[3]; //, v_double[3];
    
    /* see where each node should start reading */
    unsigned long long n_per_task = N / NTask;

    unsigned long long* nn = (unsigned long long*)malloc(NTask*sizeof(unsigned long long));
    unsigned long long* id_start = (unsigned long long*)malloc(NTask*sizeof(unsigned long long));
    
    int task, taskj;
    for(task=0; task<NTask; task++)
    {
        nn[task] = n_per_task;
        if(task < (N%NTask))
            nn[task]++;
    }
   
    id_start[0] = 0;
    for(task=1; task<NTask; task++)
    {
    	id_start[task] = 0;
    	for(taskj=0; taskj<task; taskj++)
	    	id_start[task] += nn[taskj];
    }
    
    NumPart = nn[ThisTask];
    
    xx = (float*)malloc(NumPart*sizeof(float));
    yy = (float*)malloc(NumPart*sizeof(float));
    zz = (float*)malloc(NumPart*sizeof(float));
    mm = (float*)malloc(NumPart*sizeof(float));
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    fseek(fin, head_pos+id_start[ThisTask]*particle_size, SEEK_SET);
    
    unsigned long long skip = 20; // skip velocities + soft + pot
    if(doubleVel) skip += 12;
    
    
    unsigned long long i, p;
    for(p=0,i=0; i<NumPart; i++)
    {
        xdr_float(&xdr, &mass);
        mm[p] = mass * MASS_UNIT;
        
        if(doublePos)
        {
            xdr_double(&xdr, r_double);
            xdr_double(&xdr, r_double+1);
            xdr_double(&xdr, r_double+2);
            
            r_double[0] = (r_double[0] + 0.5)*LBOX - cen[0];
            r_double[1] = (r_double[1] + 0.5)*LBOX - cen[1];
            r_double[2] = (r_double[2] + 0.5)*LBOX - cen[2];
            
//             if(norm2_double(r_double) > 160000)
//             {
//                 fseek(xdr.fp, skip, SEEK_CUR);
//                 continue;
//             }
            
            xx[p] = (float)r_double[0]/scale;
            yy[p] = (float)r_double[1]/scale;
            zz[p] = (float)r_double[2]/scale;
        }
        else
        {
            xdr_float(&xdr, r_float);
            xdr_float(&xdr, r_float+1);
            xdr_float(&xdr, r_float+2);
            
            r_float[0] = (r_float[0] + 0.5)*LBOX - cen[0];
            r_float[1] = (r_float[1] + 0.5)*LBOX - cen[1];
            r_float[2] = (r_float[2] + 0.5)*LBOX - cen[2];
            
//             if(norm2_float(r_float) > r_lim2)
//             {
//                 fseek(xdr.fp, skip, SEEK_CUR);
//                 continue;
//             }
            
            xx[p] = (float)r_float[0]/scale;
            yy[p] = (float)r_float[1]/scale;
            zz[p] = (float)r_float[2]/scale;
            
        }
        
        /* --- swap y and z ---- */
//         temp = zz[p];
//         zz[p] = yy[p];
//         yy[p] = -temp;
        
        
        // skips velocities + soft + pot
        fseek(xdr.fp, skip, SEEK_CUR);
                
        p++;

    }
    
    free(nn);
    free(id_start);
    
    fclose(fin);
    
    return p;

}

/* If your snapshot is in ASCII it's a little simpler */
void read_ascii()
{
//     char* filename = HOST_FILENAME;
    char* filename = "sphere/hernquist_sphere.dat";
    TotalNumPart = NumPart = 1000000;
    xx = (float*)malloc(NumPart*sizeof(float));
    yy = (float*)malloc(NumPart*sizeof(float));
    zz = (float*)malloc(NumPart*sizeof(float));
    mm = (float*)malloc(NumPart*sizeof(float));
    
    float scale = HOST_SCALE;

	FILE* fin = fopen(filename, "r");
	
	float this_x, this_y, this_z;
	
	long long k;
	for(k=0; k<NumPart; k++)
	{
		fscanf(fin, "%e%e%e", &this_x, &this_y, &this_z);
		
// 		printf("%e  %e  %e\n", this_x, this_y, this_z);
		
		xx[k] = this_x/scale;
		yy[k] = this_y/scale;
		zz[k] = this_z/scale;
		mm[k] = 1.0e5;
// 		mm[k] = MASS_UNIT;
		
	}

}

