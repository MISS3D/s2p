#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "vvector.h"
#include "iio.h"
#include "rpc.h"
#include "read_matrix.c"


void applyHom(double outv[3], double M[3][3], double v[3]) {
    MAT_DOT_VEC_3X3(outv, M, v);
    outv[0] /= outv[2];
    outv[1] /= outv[2];
    outv[2] /= outv[2];
}

static float interpolate_bilinear(float a, float b, float c, float d,
					float dx, float dy)
{
    double dfx = c-a;
    double dfy = b-a;
    double dfxy = a+d-c-b;
    
    return dfx*dx+dfy*dy+dfxy*dx*dy+a;
}

int main_disp_to_h(int c, char *v[])
{
    if (c != 13) {
        fprintf(stderr, "usage:\n\t"
                "%s rpca rpcb Ha Hb dispAB mskAB out_heights  RPCerr  x  y  w   h"
              // 0   1   2    3   4   5      6       7           8    9 10  11  12
                "\n", *v);
        fprintf(stderr,"c = %d\n",c);
        return EXIT_FAILURE;
    }

    // read input data
    struct rpc *rpc_list;
    int N_rpc = 2;
    rpc_list = (struct rpc *) malloc(N_rpc*sizeof( struct rpc  ));
    
    double **q_list;
    q_list = (double **) malloc(N_rpc*sizeof( double * ));
    for(int i=0; i<N_rpc; i++)
        q_list[i] = (double *) malloc(3*sizeof( double ));
    
    read_rpc_file_xml(&rpc_list[0], v[1]);
    read_rpc_file_xml(&rpc_list[1], v[2]);
    double Ha[3][3], Hb[3][3];
    read_matrix(Ha, v[3]);
    read_matrix(Hb, v[4]);
    int nx, ny, nch;
    float *dispy;
    float *dispx = iio_read_image_float_split(v[5], &nx, &ny, &nch);
    if (nch > 1) dispy = dispx + nx*ny;
    else dispy = calloc(nx*ny, sizeof(*dispy));
    float *msk  = iio_read_image_float_split(v[6], &nx, &ny, &nch);
    char *fout_heights  = v[7];
    char *fout_err = v[8];
    int X = atoi(v[9]);
    int Y = atoi(v[10]);
    int width  = atoi(v[11]);
    int height  = atoi(v[12]);
    float *heightMap = calloc(width*height, sizeof(*heightMap));
    float *errMap = calloc(width*height, sizeof(*errMap));

    // invert homography Hb
    double det;
    double invHb[3][3];
    INVERT_3X3(invHb, det, Hb);

    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
        {
            // build height map directly
            // in ref image geometry
            int posH = x + width*y;
            
            // Replace point (x,y) in the height map (h.m.)
            // by its real position inside de ref image;
            // the top-left corner of the h.m. is placed at (X,Y)
            // so : (x,y) --> (x+X,y+Y)
            // TODO : take into account zoom
            double q0[3] = {x+X, y+Y, 1.};
            double p0[3];
            
            // Now, project q0 (from ref image)
            // to the disparity map : p0
            applyHom(p0, Ha, q0);
            
            // Rounded position inside the disparity map
            int xD=(int) (p0[0]+0.5);
            int yD=(int) (p0[1]+0.5);
            int posD = xD + nx*yD;
            
            // Bilinear interpolation uses a neighborhood of radius 1.
            if ( (xD >= 1) && (xD < nx-1) && (yD >= 1) && (yD < ny-1))
            {
            
                if (msk[posD] <= 0) 
                {
                    heightMap[posH] = NAN;
                    errMap[posH] = NAN;
                } 
                else 
                {
                    // extremal coord. around p0 point
                    int x1=(int) (p0[0]);
                    int y1=(int) (p0[1]);
                    int x2=x1+1;
                    int y2=y1+1;
                    
                    // positions of the four corners around p0 point
                    int posD11 = x1 + nx*y1;
                    int posD12 = x1 + nx*y2;
                    int posD21 = x2 + nx*y1;
                    int posD22 = x2 + nx*y2;
                    
                    // Diff between first corner with floating position
                    double Dx = (double) (p0[0]-x1);
                    double Dy = (double) (p0[1]-y1);
                    
                    // Interpolate disparity map
                    double dx = interpolate_bilinear(dispx[posD11],
                                                     dispx[posD12],
                                                     dispx[posD21],
                                                     dispx[posD22],
                                                     Dx,Dy);
                                                     
                    double dy = interpolate_bilinear(dispy[posD11],
                                                     dispy[posD12],
                                                     dispy[posD21],
                                                     dispy[posD22],
                                                     Dx,Dy);
                    
                    // Map point q0 (from ref image) 
                    // to its mate q1 (from slave image)
                    double p1[3] = {p0[0]+dx, p0[1]+dy, 1.};
                    double q1[3];
                    
                    applyHom(q1, invHb, p1);
                    
                    q_list[0] = q0;
                    q_list[1] = q1;
                    
                    // Compute the related height & rpc_err
                    double err, hg;
                    hg = rpc_height_geo(rpc_list, q_list, N_rpc,  &err);

                    // Output the result in height & rpc_err maps,
                    // both in original geometry
                    heightMap[posH] = hg;
                    errMap[posH] = err;
                    
                }
            }
            else
            {
                    heightMap[posH] = NAN;
                    errMap[posH] = NAN;
            } 
        }

    // save the height map and error map
    iio_save_image_float_vec(fout_heights, heightMap, width, height, 1);
    iio_save_image_float_vec(fout_err, errMap, width, height, 1);
    
    // free mem
    free(rpc_list);
    free(q_list);

    return 0;
}

int main(int c, char *v[])
{
    return main_disp_to_h(c, v);
}
