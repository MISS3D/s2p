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


int main_disp_to_heights(int c, char *v[])
{
    if (c < 12) {
        fprintf(stderr, "usage:\n\t"
                "%s global_out_dir, col, row, tw, th, z, trg_cons, thr_cons, full_outputs, rpc_ref.xml, rpc_slave_1.xml, [rpc_slave_2.xml ...]"
              // 0         1         2    3   4   5   6      7        8          9         10...
                "\n", *v);
        fprintf(stderr,"c = %d\n",c);
        return EXIT_FAILURE;
    }
    

    // read input data
    char tmp_path[1000];
    double det;
    
    char *global_out_dir=v[1];
    int X=atoi(v[2]);
    int Y=atoi(v[3]);
    int width=atoi(v[4]);
    int height=atoi(v[5]);
    double z = atoi(v[6]);
    bool trg_cons = atoi(v[7]);
    double thr_cons = atof(v[8]);
    bool full_outputs = atoi(v[9]);
    
    // build tile_dir path
    char tile_dir[1000];
    sprintf(tile_dir,"%s/tile_%d_%d_row_%d/col_%d",global_out_dir,width,height,Y,X);
    
    // take into account zoom
    width = (int) (width/z)+1;
    height = (int) (height/z)+1;
    
    // useful global variable
    int N_rpc = c-10;
    int nb_pairs = N_rpc-1;
    int nb_sights = N_rpc;


    // ################################
    // build the list of pairs
    // ################################
    
    list_of_pairs list_pairs;
    list_pairs.tot_size = nb_pairs;
    list_pairs.real_size = 0;
    list_pairs.data = (type_pair *) malloc(list_pairs.tot_size*sizeof(type_pair));
    // pid == pair id
    for(int pid=1;pid<=nb_pairs;pid++)
    {
        sprintf(tmp_path,"%s/pair_%d/H_ref.txt",tile_dir,pid);
        FILE *test = fopen(tmp_path,"r");
        if ( test != NULL)
            {
                int index = list_pairs.real_size;
                
                list_pairs.real_size++;
                list_pairs.data[index].sight_slave = pid+1;
                
                // read RPC
                read_rpc_file_xml(&list_pairs.data[index].rpc_master,v[10]);
                read_rpc_file_xml(&list_pairs.data[index].rpc_slave,v[10+pid]);
                
                // read homographies
                sprintf(tmp_path,"%s/pair_%d/H_ref.txt",tile_dir,pid);
                read_matrix(list_pairs.data[index].H_ref.val, tmp_path); 
                sprintf(tmp_path,"%s/pair_%d/H_sec.txt",tile_dir,pid);
                read_matrix(list_pairs.data[index].H_sec.val, tmp_path);
                
                // read global/local pointing corrections by pair
                sprintf(tmp_path,"%s/pair_%d/pointing.txt",tile_dir,pid);
                FILE *test2 = fopen(tmp_path,"r");
                if ( test2 == NULL)
                    sprintf(tmp_path,"%s/global_pointing_pair_%d.txt",global_out_dir,pid);
                fclose(test2);
                read_matrix(list_pairs.data[index].pointing_correc.val, tmp_path); 
                
                // read disparity maps by pair
                sprintf(tmp_path,"%s/pair_%d/rectified_disp.tif",tile_dir,pid);
                int *nx = &list_pairs.data[index].nx;
                int *ny = &list_pairs.data[index].ny;
                int *nch = &list_pairs.data[index].nch;
                list_pairs.data[index].dispx = iio_read_image_float_split(tmp_path,nx,ny,nch);
                if (*nch > 1) 
                    list_pairs.data[index].dispy = list_pairs.data[index].dispx + (*nx)*(*ny);
                else 
                    list_pairs.data[index].dispy = calloc((*nx)*(*ny), sizeof(float *));
                    
                // read mask
                sprintf(tmp_path,"%s/pair_%d/rectified_mask.png",tile_dir,pid);
                int nx_mask,ny_mask,nch_mask;
                list_pairs.data[index].mask = iio_read_image_float_split(tmp_path, &nx_mask, &ny_mask, &nch_mask);
                
                // take into account pointing corrections for slave homographies
                double invA[3][3];
                INVERT_3X3(invA, det, list_pairs.data[index].pointing_correc.val);
                MATRIX_PRODUCT_3X3(list_pairs.data[index].H_secB.val,
                                   list_pairs.data[index].H_sec.val,invA);
                 
                 // slave homography inversion                  
                INVERT_3X3(list_pairs.data[index].H_invSecB.val, det, 
                           list_pairs.data[index].H_secB.val);
                
                fclose(test);
            }
    }

    // need at least one pair to do something
    if (list_pairs.real_size == 0)
        return 1;

    // ##################################
    // build outputs paths and alloc mem
    // ##################################
    
    // * heights
    char fout_heights[1000];
    sprintf(fout_heights,"%s/height_map.tif",tile_dir);
    float *heightMap = (float *) calloc(width*height, sizeof(float));
    
    // * rpc errors
    typedef char tabchar[1000];
    int size_of_fout_err_tab;
    if (full_outputs)   size_of_fout_err_tab = nb_sights+1;
    else            size_of_fout_err_tab = 1;
    
    tabchar *fout_err_tab = (tabchar *) malloc( (size_of_fout_err_tab)*sizeof(tabchar) );
    for(int i=1; i<size_of_fout_err_tab; i++)
        sprintf(fout_err_tab[i],"%s/rpc_err_norm_sight_%d.tif",tile_dir,i);
    sprintf(fout_err_tab[0],"%s/rpc_err_rms_allsights.tif",tile_dir);
    float **errMap_tab = (float **) malloc((size_of_fout_err_tab)*sizeof(float *));
    for(int i=0; i<size_of_fout_err_tab; i++)
        errMap_tab[i] = (float *) calloc(width*height, sizeof(float));
        
    // Yet another interesting thing :
    // track the vector from an opt. point
    // to a given view
    tabchar *fout_vec_tab;
    float **img_vec_tab;
    if (full_outputs)
    {
        fout_vec_tab = (tabchar *) malloc( nb_sights*sizeof(tabchar) );
        for(int i=0; i<nb_sights; i++)
            sprintf(fout_vec_tab[i],"%s/rpc_err_vec_sight_%d.tif",tile_dir,i+1);
            
        img_vec_tab = (float **) malloc(nb_sights*sizeof( float * ));
        for(int i=0;i<nb_sights;i++)
            img_vec_tab[i] = (float *) malloc(3*width*height*sizeof( float ));
    }
    
    // Yet another interesting thing :
    // reproject the above vectors
    // into original geometry
    tabchar *fout_rpj_tab;
    float **rpj_vec_tab;
    if (full_outputs)
    {
        fout_rpj_tab = (tabchar *) malloc(nb_sights*sizeof(tabchar) );
        for(int i=0; i<nb_sights; i++)
            sprintf(fout_rpj_tab[i],"%s/rpc_err_rpjvec_sight_%d.tif",tile_dir,i+1);
            
        rpj_vec_tab = (float **) malloc(nb_sights*sizeof( float * ));
        for(int i=0;i<nb_sights;i++)
            rpj_vec_tab[i] = (float *) malloc(3*width*height*sizeof( float ));
    }
    
    // Yet another interesting thing :
    // output the 2D disparities
    // for a given image pair
    tabchar *fout_disp2D_tab;
    float ** disp2D_tab;
    if (full_outputs)
    {
        fout_disp2D_tab = (tabchar *) malloc(nb_pairs*sizeof(tabchar) );
        for(int i=0; i<nb_pairs; i++)
            sprintf(fout_disp2D_tab[i],"%s/pair_%d/disp2D.tif",tile_dir,i+1);
        
        disp2D_tab = (float **) malloc(nb_pairs*sizeof( float * ));
        for(int i=0;i<nb_pairs;i++)
            disp2D_tab[i] = (float *) malloc(3*width*height*sizeof( float ));
    }
    
    // Yet another interesting thing :
    // we wish we could know
    // the number of views used for each pixel
    // as well as which views have been used for each pixel
    int *nb_views; // view <=> sight
    int **img_selected_views;
    char fnb_views[1000];
    tabchar * fout_selected_views;
    if (full_outputs)
    {
        sprintf(fnb_views,"%s/nb_sights.tif",tile_dir);
        nb_views = (int *) calloc(width*height, sizeof(int));

        fout_selected_views = (tabchar *) malloc( nb_sights*sizeof(tabchar) );
        for(int i=0; i<nb_sights; i++)
            sprintf(fout_selected_views[i],"%s/selected_sight_%d.tif",tile_dir,i+1);
        img_selected_views = (int **) malloc(nb_sights*sizeof( int * ));
        for(int i=0;i<nb_sights;i++)
            img_selected_views[i] = (int *) calloc(width*height,sizeof( int ));
    }
    
    
    // ################################
    // time to build the height map 
    // directly in ref image geometry !
    // ################################
    for (int y = 0; y < height; y++) 
        for (int x = 0; x < width; x++) 
        {
            int local_nb_sights = list_pairs.real_size +1;
            
            // position inside height map
            int posH = x + width*y;
            
            // for each pair (pid == pair id)
            for(int pid=0;pid<list_pairs.real_size;pid++) 
            {
                // 2D disparity initialization
                if (full_outputs)
                  for(int t=0; t<3; t++)
                    disp2D_tab[pid][width*3*y+3*x+t] = NAN;
                
                // Replace point (x,y) in the height map (h.m.)
                // by its real position inside de ref image:
                // (x,y) --> (x.z+X,y.z+Y)
                list_pairs.data[pid].q0[0] = x*z+X;
                list_pairs.data[pid].q0[1] = y*z+Y;
                list_pairs.data[pid].q0[2] = 1.;
                
                // Now, project q0 (from ref image)
                // to the point p0 inside disparity map
                double p0[3];
                applyHom(p0, list_pairs.data[pid].H_ref.val, 
                            list_pairs.data[pid].q0);
                
                // Rounded position inside the disparity map
                int xD=(int) (p0[0]+0.5);
                int yD=(int) (p0[1]+0.5);
                int posD = xD + list_pairs.data[pid].nx*yD;
                
                // determine whether p0 is inside the disp map and not masked
                list_pairs.data[pid].process = false;
                if ( (xD >= 1) && (xD < list_pairs.data[pid].nx -1) 
                    && (yD >= 1) && (yD < list_pairs.data[pid].ny-1) 
                    && (list_pairs.data[pid].mask[posD]>0)
                   )
                {
                    list_pairs.data[pid].process = true;
                    
                    // Bilinear interpolation uses a neighborhood of radius 1
                    // *extremal coord. around p0 point
                    int x1=(int) (p0[0]);
                    int y1=(int) (p0[1]);
                    int x2=x1+1;
                    int y2=y1+1;
                    
                    // *positions of the four corners around p0 point
                    int posD11 = x1 + list_pairs.data[pid].nx*y1;
                    int posD12 = x1 + list_pairs.data[pid].nx*y2;
                    int posD21 = x2 + list_pairs.data[pid].nx*y1;
                    int posD22 = x2 + list_pairs.data[pid].nx*y2;
                    
                    // *diff between first corner with floating position
                    double Dx = (double) (p0[0]-x1);
                    double Dy = (double) (p0[1]-y1);
                    
                    // *interpolate disparity map
                    double dx = interpolate_bilinear(list_pairs.data[pid].dispx[posD11],
                                                     list_pairs.data[pid].dispx[posD12],
                                                     list_pairs.data[pid].dispx[posD21],
                                                     list_pairs.data[pid].dispx[posD22],
                                                     Dx,Dy);
                                                     
                    double dy = interpolate_bilinear(list_pairs.data[pid].dispy[posD11],
                                                     list_pairs.data[pid].dispy[posD12],
                                                     list_pairs.data[pid].dispy[posD21],
                                                     list_pairs.data[pid].dispy[posD22],
                                                     Dx,Dy);
                                                     
                    // Map point q0 (from ref image) 
                    // to its mate q1 (from slave image)
                    double p1[3] = {p0[0]+dx, p0[1]+dy, 1.};
                    applyHom(list_pairs.data[pid].q1, 
                             list_pairs.data[pid].H_invSecB.val, p1);
                    
                    // 2D disparities
                    if (full_outputs)
                      for(int t=0; t<2; t++)
                        disp2D_tab[pid][width*3*y+3*x+t] = list_pairs.data[pid].q1[t]
                                                         -list_pairs.data[pid].q0[t];
                }
                else
                    local_nb_sights--;
            }

            // some initializations
            heightMap[posH] = NAN;
            for(int i=0; i<size_of_fout_err_tab; i++)
                errMap_tab[i][posH] = NAN;

            if (full_outputs)
                for(int i=0; i<nb_sights; i++)
                    for(int t=0; t<3; t++)
                    {
                        img_vec_tab[i][width*3*y+3*x+t] = NAN;
                        rpj_vec_tab[i][width*3*y+3*x+t] = NAN;
                    }

            
            if (local_nb_sights>=2) // at least two sights are required
            {
                // Altitude
                double hg;

                // rpc_height_geo is able to filter out
                // aberrant views (including ref view !)
                int final_nb_sights;
                
                // list of sights
                type_sight *sights_list = (type_sight *) malloc(local_nb_sights*sizeof(type_sight));
                
                // Compute the related height & rpc_err
                hg = rpc_height_geo(&list_pairs, 
                                    local_nb_sights, &final_nb_sights,
                                    trg_cons, thr_cons, sights_list);
                
                // 1) errors
                double sum = 0.0;
                double nb_elt = 0.0;
                if (final_nb_sights>=2)
                {
                    for(int s=0; s<local_nb_sights; s++)
                    {
                        // The r.m.s error is computed
                        // only from the sights being part
                        // of the consensus !!!
                        if (sights_list[s].consensus)
                        {
                            sum += pow(sights_list[s].err,2.0);
                            nb_elt++;
                        }
                        if (full_outputs)
                        { 
                          // err by sight : dist from opt point to sight i 
                          // --> for all sights, even those not in the consensus !
                          // (index 0 of errMap_tab is dedicated to rms error)
                          errMap_tab[sights_list[s].ID][posH] = sights_list[s].err;
                          
                          // vec err by sight : vector from opt point to sight i 
                          // --> for all sights, even those not in the consensus !
                          // (sights numbered from 1 to N so index = ID-1)
                          for(int t=0; t<3; t++)
                            img_vec_tab[sights_list[s].ID-1][width*3*y+3*x+t] = 
                                 sights_list[s].err_vec_ptopt_to_sight[t+3]
                                 -sights_list[s].err_vec_ptopt_to_sight[t];
                                 
                          // same vectors as above, but reprojected into
                          // ref image (so two components : drow and dcol)
                          // (sights numbered from 1 to N so index = ID-1)
                          double lgt1,lat1,alt1,pos1[2];
                          double lgt2,lat2,alt2,pos2[2];
                          ECEF_to_lgt_lat_alt(sights_list[s].err_vec_ptopt_to_sight[0], 
                                            sights_list[s].err_vec_ptopt_to_sight[1], 
                                            sights_list[s].err_vec_ptopt_to_sight[2],
                                            &lgt1,&lat1,&alt1);                                     
                          ECEF_to_lgt_lat_alt(sights_list[s].err_vec_ptopt_to_sight[3], 
                                            sights_list[s].err_vec_ptopt_to_sight[4], 
                                            sights_list[s].err_vec_ptopt_to_sight[5],
                                            &lgt2,&lat2,&alt2);
                                            
                          eval_rpci(pos1, &list_pairs.data[0].rpc_master, lgt1, lat1, alt1);
                          eval_rpci(pos2, &list_pairs.data[0].rpc_master, lgt2, lat2, alt2);
                        
                          for(int t=0; t<2; t++)
                            rpj_vec_tab[sights_list[s].ID-1][width*3*y+3*x+t] = pos2[t]-pos1[t];
                          
                          // nb of sights used for current pixel
                          nb_views[posH] = final_nb_sights;
                          
                          // tells which sight has been used for current pixel
                          // (sights numbered from 1 to N so index = ID-1)
                          if (sights_list[s].consensus)
                            img_selected_views[sights_list[s].ID-1][posH]=1;
                        }
                    }
                    
                    // r.m.s. error
                    sum = sqrt(sum/nb_elt);
                  
                    // Output the results in original geometry 
                    // * height
                    heightMap[posH] = hg;
                    // * r.m.s. error (dedicated to index 0)
                    errMap_tab[0][posH] = sum;
                }
                
                free(sights_list);
            } 
        }
        
        
    // save the height map / error map / nb_views
    iio_save_image_float(fout_heights, heightMap, width, height);
    for(int i=0;i<size_of_fout_err_tab;i++)
        iio_save_image_float(fout_err_tab[i], errMap_tab[i], width, height);
    if (full_outputs)
    {
        iio_save_image_int(fnb_views, nb_views, width, height);
        for(int i=0;i<N_rpc;i++)
        {
            iio_save_image_int(fout_selected_views[i], img_selected_views[i], width, height);
            iio_save_image_float_vec(fout_vec_tab[i], img_vec_tab[i], width, height, 3);
            iio_save_image_float_vec(fout_rpj_tab[i], rpj_vec_tab[i], width, height, 3);
        }
        for(int i=0;i<nb_pairs;i++)
            iio_save_image_float_vec(fout_disp2D_tab[i], disp2D_tab[i], width, height, 3);
    }

    // clean mem
    free(heightMap);
    if (full_outputs)
        free(nb_views);
    for(int i=0;i<list_pairs.real_size;i++)
    {
        free(list_pairs.data[i].dispx);
        free(list_pairs.data[i].mask);
        if (list_pairs.data[i].nch == 1)
            free(list_pairs.data[i].dispy);
    }
    free(list_pairs.data);
    
    for(int i=0;i<size_of_fout_err_tab;i++)
        free(errMap_tab[i]);
    free(errMap_tab);
    free(fout_err_tab);
    if (full_outputs)
    {
        for(int i=0; i<N_rpc; i++)
        {
            free(img_selected_views[i]);
            free(img_vec_tab[i]);
            free(rpj_vec_tab[i]);
        }
        for(int i=0; i<nb_pairs; i++)
        {
            free(disp2D_tab[i]);
        }
        free(img_selected_views);
        
        free(img_vec_tab);
        free(rpj_vec_tab);
        free(fout_selected_views);
        free(fout_vec_tab);
        free(fout_rpj_tab);
        free(fout_disp2D_tab);
        free(disp2D_tab);
    }
    return 0;
}

int main(int c, char *v[])
{
    return main_disp_to_heights(c, v);
}
