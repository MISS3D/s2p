// take a series of ply files and produce a digital elevation map

#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <geotiff/xtiffio.h>
#include <geotiff/geotiffio.h>
#include <geotiff/geo_tiffp.h>

#include "iio.h"
#include "lists.c"
#include "fail.c"
#include "xmalloc.c"


// convert string like '28N' into a number like 32628, according to:
// WGS84 / UTM northern hemisphere: 326zz where zz is UTM zone number
// WGS84 / UTM southern hemisphere: 327zz where zz is UTM zone number
// http://www.remotesensing.org/geotiff/spec/geotiff6.html#6.3.3.1
static int get_utm_zone_index_for_geotiff(char *utm_zone)
{	
	int out = 32000;
	if (utm_zone[2] == 'N')
		out += 600;
	else if (utm_zone[2] == 'S')
		out += 700;
	else
		fprintf(stderr, "error: bad utm zone value: %s\n", utm_zone);
	utm_zone[2] = '\0';
	out += atoi(utm_zone);
	return out;
}

void set_geotif_header(char *tiff_fname, char *utm_zone, float xoff,
		float yoff, float scale)
{
	// open tiff file
	TIFF *tif = XTIFFOpen(tiff_fname, "r+");
	if (!tif)
		fail("failed in XTIFFOpen\n");

	GTIF *gtif = GTIFNew(tif);
	if (!gtif)
		fail("failed in GTIFNew\n");

	// write TIFF tags
	double pixsize[3] = {scale, scale, 0.0};
	TIFFSetField(tif, GTIFF_PIXELSCALE, 3, pixsize);

	double tiepoint[6] = {0.0, 0.0, 0.0, xoff, yoff, 0.0};
	TIFFSetField(tif, GTIFF_TIEPOINTS, 6, tiepoint);

	// write GEOTIFF keys
	int utm_ind = get_utm_zone_index_for_geotiff(utm_zone);
	GTIFKeySet(gtif, ProjectedCSTypeGeoKey, TYPE_SHORT, 1, utm_ind);
	GTIFWriteKeys(gtif);

	// free and close
	GTIFFree(gtif);
	XTIFFClose(tif);
}



struct ply_property {
	enum {UCHAR,FLOAT,DOUBLE,UNKNOWN} type;
	char name[0x100];
	size_t len;
};

static bool parse_property_line(struct ply_property *t, char *buf)
{
	char typename[0x100];
	bool r = 2 == sscanf(buf, "property %s %s\n", typename, t->name);
	t->type = UNKNOWN;
	if (0 == strcmp(typename, "uchar")) { t->type = UCHAR;  t->len = 1;}
	if (0 == strcmp(typename, "float")) { t->type = FLOAT;  t->len = 4;}
	if (0 == strcmp(typename, "double")){ t->type = DOUBLE; t->len = 8;}
	return r;
}



// fast forward "f" until "end_header" is found
// returns the number of 'properties' 
// the array of structures *t, contains the names and sizes 
// the properties in bytes, isbin is set if binary encoded
// and reads the utm zone
static size_t header_get_record_length_and_utm_zone(FILE *f_in, char *utm, 
		int *isbin, struct ply_property *t)
{
	size_t n = 0;
	*isbin = 0;

	char buf[FILENAME_MAX] = {0};
	while (fgets(buf, FILENAME_MAX, f_in)) {
		if (0 == strcmp(buf, "format binary_little_endian 1.0\n")) *isbin=1;
		else if (0 == strcmp(buf, "format ascii 1.0\n")) *isbin=0;
		else {
			if (parse_property_line(t+n, buf))
				n += 1;
			else if (0 == strncmp(buf, "comment projection:", 19)) {
				sscanf(buf, "comment projection: UTM %s", utm);
			}
		}
		if (0 == strcmp(buf, "end_header\n"))
			break;
	}
	return n;
}



// re-scale a float between 0 and w
static int rescale_float_to_int(double x, double min, double max, int w, bool *flag)
{
	*flag=1;
	int r = w * (x - min)/(max - min);
	if (r < 0) *flag=0;
	if (r >= w) *flag=0;
	return r;
}


typedef struct {
	float x;
	float y;
} Position;

struct images {
	float *cnt;
	float *pixel_value;
	float **heights;
	Position **pos;
	int w, h;
};

// Help to sort tabs of float
int compare (const void * a, const void * b)
{
  return ( *(float*)a - *(float*)b );
}

double weight(Position pos, Position center_pos, unsigned int flag,float pinterp)
{
    double eps=10e-3;
    
    double d = pow(center_pos.x-pos.x,2.0)+pow(center_pos.y-pos.y,2.0);
    d = sqrt(d);
    
    switch (flag) 
    {
	case 6: 
	    return 1.0/(d+eps);
	case 7:
	    return exp(-pinterp*d*d);
	case 8:
	    return 1/(1+exp(pinterp*d));
    }
}


// update the output images with a new height
static void add_height_to_images(struct images *x, int i, int j, float v, Position pos, int flag)
{
    uint64_t k = (uint64_t) x->w * j + i;
    
    switch (flag) 
    {
	case 0: // nominal case
	{
	    x->pixel_value[k] = (v + x->cnt[k] * x->pixel_value[k]) / (1 + x->cnt[k]);
	    x->cnt[k] += 1;
	}
	break;
	case -3: // just count the number of occurrences
	{
	    x->cnt[k] += 1;
	}
	break;
	case -2: // memory alloc and heights tab filling
	{
	    if (x->cnt[k])
	    {
		if (!x->heights[k])
		{
		    x->heights[k] = xmalloc(x->cnt[k]*sizeof(float));
		    x->cnt[k]=0;
		}
		
		x->heights[k][(int) x->cnt[k]] = v;
		x->cnt[k] += 1;
	    }
	}
	break;
	case -1: // memory alloc and tab filling for heights/pos
	{
	    if (x->cnt[k])
	    {
		if ( (!x->heights[k]) && (!x->pos[k]) ) 
		{
		    x->heights[k] = xmalloc(x->cnt[k]*sizeof(float));
		    x->pos[k] = xmalloc(x->cnt[k]*sizeof(Position));
		    x->cnt[k]=0;
		}
		
		x->heights[k][(int) x->cnt[k]] = v;

		x->pos[k][(int) x->cnt[k]] = pos;
		x->cnt[k] += 1;
	    }
	}
	break;
    }
}

static void synth_heights(struct images *x, int i, int j, Position center_pos, int flag, int radius, int minnonan, float pinterp)
{
    uint64_t k = (uint64_t) x->w * j + i;
    
    switch (flag) 
    {	   	    
        case 1: // average
        {
            if (x->cnt[k])
            {
            float sum=0.;
            for(int t=0;t<x->cnt[k];t++)
                sum += x->heights[k][t];
            x->pixel_value[k] = sum / ( (float) x->cnt[k]);
            }
        }
        break;
        case 2: // var
        {
            if (x->cnt[k])
            {
            double sum1=0.,sumC=0.;
            for(int t=0;t<x->cnt[k];t++)
            {
                sum1 += (double) x->heights[k][t];
                sumC += pow( (double) x->heights[k][t],2.0);
            }
            double m1 = sum1 / ( (double) x->cnt[k]);
            double mc = sumC / ( (double) x->cnt[k]);

            x->pixel_value[k] = mc-m1*m1;
            }
        }
        break;
        case 3: // min
        {
            if (x->cnt[k])
            {
            qsort (x->heights[k], (int) x->cnt[k], sizeof(float), compare);
            x->pixel_value[k] = x->heights[k][0];
            }
        }
        break;
        case 4: // max
        {
            if (x->cnt[k])
            {
            qsort (x->heights[k], (int) x->cnt[k], sizeof(float), compare);
            x->pixel_value[k] = x->heights[k][(int) x->cnt[k]-1];
            }
        }
        break;
        case 5: // median
        {
            if (x->cnt[k])
            {
            qsort (x->heights[k], (int) x->cnt[k], sizeof(float), compare);
            x->pixel_value[k] = x->heights[k][(int) x->cnt[k]/2];
            }
        }
        break;
    }
    
    if (flag>=6)// weighted by dist from center of cell
    {
        double w;
        double sum=0.0,weighted_moy=0.0;
	int found = 0;
        
        if (x->cnt[k]) // Do not interpolate
        {    
	    radius = 0; // But still use the distances from the cell center to weight heights
	    found = minnonan;
	}
        
        for(int ii=-radius; ii<=radius; ii++)
            for(int jj=-radius; jj<=radius; jj++)
            if ( (i+ii>=0) && (i+ii<x->w) && (j+jj>=0) && (j+jj<x->h) )
            {
                uint64_t kt = (uint64_t) x->w * (j+jj) + i+ii;
                if (x->cnt[kt])
                {
                    found++;
                    for(int t=0;t<x->cnt[kt];t++)
                    {
                        w = weight(x->pos[kt][t],center_pos,flag,pinterp);
                        sum += w;
                        weighted_moy += w * ( (double) x->heights[kt][t] );
                    }
                }
            }
        if (found>=minnonan)
            x->pixel_value[k] = weighted_moy/sum;
    }
}

int get_record(FILE *f_in, int isbin, struct ply_property *t, int n, double *data){
	int rec = 0;
	if(isbin) {
		for (int i = 0; i < n; i++) {
			switch(t[i].type) {
				case UCHAR: {
						    unsigned char X;
						    rec += fread(&X, 1, 1, f_in);
						    data[i] = X;
						    break; }
				case FLOAT: {
						    float X;
						    rec += fread(&X, sizeof(float), 1, f_in);
						    data[i] = X;
						    break; }
				case DOUBLE: {
						     double X;
						     rec += fread(&X, sizeof(double), 1, f_in);
						     data[i] = X;
						     break; }
			}
		}
	} else {
		int i=0;
		while (i < n && !feof(f_in)) {
			rec += fscanf(f_in,"%lf", &data[i]);  i++;
		}
	}
	return rec;
}


// open a ply file, and accumulate its points to the image
static void add_ply_points_to_images(struct images *x,
		float xmin, float xmax, float ymin, float ymax,
		char utm_zone[3], char *fname, int col_idx, unsigned int flag)
{
	FILE *f = fopen(fname, "r");
	if (!f) {
		fprintf(stderr, "WARNING (from add_ply_points_to_images) : can not open file \"%s\"\n", fname);
		return;
	}

	// check that the utm zone is the same as the provided one
	char utm[3];
	int isbin=1;
	struct ply_property t[100];
	size_t n = header_get_record_length_and_utm_zone(f, utm, &isbin, t);
	if (0 != strncmp(utm_zone, utm, 3))
		fprintf(stderr, "error: different UTM zones among ply files\n");

	if (col_idx < 2 || col_idx > 5)
		exit(fprintf(stderr, "error: bad col_idx %d\n", col_idx));


	double data[n];
	double center_x,center_y,d;
	while ( n == get_record(f, isbin, t, n, data) ) {
		bool flag1,flag2;
		int i = rescale_float_to_int(data[0], xmin, xmax, x->w, &flag1);
		int j = rescale_float_to_int(-data[1], -ymax, -ymin, x->h, &flag2);
		
		Position pos;
		pos.x=data[0];
		pos.y=-data[1];
				
		if ( (flag1) && (flag2))
		{
		    if (col_idx == 2) {
			    add_height_to_images(x, i, j, data[2], pos, flag);
			    assert(isfinite(data[2]));
		    }
		    else
		    {
			    unsigned int rgb = data[col_idx];
			    add_height_to_images(x, i, j, rgb, pos, flag);
		    }
		}
	}

	fclose(f);
}


void help(char *s)
{
	fprintf(stderr, "usage:\n\t"
			"%s [-c column] [-flag flag] [-radius radius] [-minnonan minnonan] [-param_inter param_inter] \
			resolution out_dsm xmin xmax ymin ymax \
			rowmin steprow rowmax colmin stepcol colmax tw th root_out_dir\n", s);
	fprintf(stderr, "\t the resolution is in meters per pixel\n");
}

#include "pickopt.c"

int main(int c, char *v[])
{
	int col_idx = atoi(pick_option(&c, &v, "c", "2"));
	int flag = atoi(pick_option(&c, &v, "flag", "0"));
	int radius = atoi(pick_option(&c, &v, "radius", "0"));
	int minnonan = atoi(pick_option(&c, &v, "minnonan", "0"));
	float param_inter = atof(pick_option(&c, &v, "pinterp", "1"));

	// process input arguments
	if (c != 16) {
		help(*v);
		return 1;
	}
	float resolution = atof(v[1]);
	char *out_dsm = v[2];
	
	float xmin = atof(v[3]);
	float xmax = atof(v[4]);
	float ymin = atof(v[5]);
	float ymax = atof(v[6]);
	float xmin_orig=xmin;
	float xmax_orig=xmax;
	float ymin_orig=ymin;
	float ymax_orig=ymax;
	fprintf(stderr, "xmin: %20f, xmax: %20f, ymin: %20f, ymax: %20f\n", xmin,xmax,ymin,ymax);
	int rowmin = atoi(v[7]);
	int steprow = atoi(v[8]);
	int rowmax = atoi(v[9]);
	int colmin = atoi(v[10]);
	int stepcol = atoi(v[11]);
	int colmax = atoi(v[12]);
	int tw = atoi(v[13]);
	int th = atoi(v[14]);
	char *root_out_dir=v[15];
	
	if (flag>=6) // interpolation
	{
	    xmin -= radius*resolution;
	    xmax += radius*resolution;
	    ymin -= radius*resolution;
	    ymax += radius*resolution;
	    fprintf(stderr, "interpolation --> xmin: %20f, xmax: %20f, ymin: %20f, ymax: %20f\n", xmin,xmax,ymin,ymax);
	}

	// process each filename to determine x, y extremas and store the
	// filenames in a list of strings, to be able to open the files again
	FILE* ply_extrema_file = NULL;
	
	char ply[1000];
	char ply_extrema[1000];
	char utm[3];
	
	float local_xmin,local_xmax,local_ymin,local_ymax;
	
	struct list *l = NULL;
	
	// From the list of tiles, find each ply file
	uint64_t nbply_pushed=0;
	int row,col;
	for(row=rowmin;row<=rowmax;row+=steprow)
	    for(col=colmin;col<=colmax;col+=stepcol)
	    {
	       //strtok(tile_dir, "\n");
	       sprintf(ply_extrema,"%s/tile_%d_%d_row_%d/col_%d/plyextrema.txt",root_out_dir,tw,th,row,col);
	       
	       // Now, find the extent of a given ply file, 
	       // specified by [local_xmin local_xmax local_ymin local_ymax]
	       ply_extrema_file = fopen(ply_extrema, "r");
	       if (ply_extrema_file)
	       {
		  fscanf(ply_extrema_file, "%f %f %f %f", &local_xmin, &local_xmax, &local_ymin, &local_ymax);
		  fclose(ply_extrema_file);

		  // Only add ply files that intersect the extent specified by [xmin xmax ymin ymax]
		  // The test below simply tells whether two rectancles overlap
		  if ( (local_xmin <= xmax) && (local_xmax >= xmin) && (local_ymin <= ymax) && (local_ymax >= ymin) )
		   {
		       sprintf(ply,"%s/tile_%d_%d_row_%d/col_%d/cloud.ply",root_out_dir,tw,th,row,col);
		       // Record UTM zone
		       FILE *ply_file = fopen(ply, "r");
			if (ply_file) 
			{
			    l = push(l, ply);
			    nbply_pushed++;
			    int isbin=0;
			    struct ply_property t[100];
			    size_t n = header_get_record_length_and_utm_zone(ply_file, utm, &isbin, t);
			    fclose(ply_file);
			}
			else
			    fprintf(stderr, "WARNING 2 : can not open file \"%s\"\n", ply);
		   }
	       }
	       else
		    fprintf(stderr,"WARNING 1 : can not open file %s\n",ply_extrema);
	    } //end for loops
		
	if (nbply_pushed == 0)
	{
		fprintf(stderr, "ERROR : no ply file pushed\n");
		return 1;
	}
		
	// compute output image dimensions
	int w = 1 + (xmax - xmin) / resolution;
	int h = 1 + (ymax - ymin) / resolution;

	// allocate and initialize output images
	struct images x;
	x.w = w;
	x.h = h;
	x.cnt = xmalloc((uint64_t) w*h*sizeof(float));
	x.pixel_value = xmalloc((uint64_t) w*h*sizeof(float));
	if (flag != 0)
	{
	    x.heights = xmalloc((uint64_t) w*h*sizeof(float *));
	    if (flag>=6) // interpolation
		x.pos = xmalloc((uint64_t) w*h*sizeof(Position *));
	}

	for (uint64_t i = 0; i < (uint64_t) w*h; i++)
	{
		x.cnt[i] = 0;
		x.pixel_value[i] = 0;
		if (flag != 0)
		{
		    x.heights[i] = NULL;
		    if (flag>=6) // interpolation
			x.pos[i] = NULL;
		}
	}


	// process each filename to accumulate points in the dem
	struct list *begin = l;
	
	if (flag==0)
	{
	    while (l != NULL)
	    {
		    // printf("FILENAME: \"%s\"\n", l->current);
		    add_ply_points_to_images(&x, xmin, xmax, ymin, ymax, utm, l->current, col_idx,flag);
		    l = l->next;
	    }

	    // set unknown values to NAN
	    for (uint64_t i = 0; i < (uint64_t) w*h; i++)
		    if (!x.cnt[i])
			    x.pixel_value[i] = NAN;
	}
	else
	{

	    l=begin;
	    while (l != NULL)
	    {
		// printf("FILENAME: \"%s\"\n", l->current);
		add_ply_points_to_images(&x, xmin, xmax, ymin, ymax, utm, l->current, col_idx,-3);
		l = l->next;
	    }
	    // set unknown values to NAN
	    for (uint64_t i = 0; i < (uint64_t) w*h; i++)
		    if (!x.cnt[i])
			    x.pixel_value[i] = NAN;
	    
	    l=begin;
	    while (l != NULL)
	    {
		// printf("FILENAME: \"%s\"\n", l->current);
		int n;
		if (flag <6) // no interpolation
		    n=-2;
		else
		    n=-1;
		add_ply_points_to_images(&x, xmin, xmax, ymin, ymax, utm, l->current, col_idx,n);
		l = l->next;
	    }
	    
	    // heights synthesis 
	    Position center_pos;
	    for (int i = 0; i < w; i++)
		for (int j = 0; j < h; j++)
		{
		    center_pos.x=xmin+resolution/2.0 +(xmax-xmin)/( (float) w)*i;
		    center_pos.y=-ymax+resolution/2.0 +(ymax-ymin)/( (float) h)*j;
		    synth_heights(&x,i,j,center_pos,flag,radius,minnonan,param_inter);
		}
	}
	
	// cleanup and exit (1)
	free(x.cnt);
	for (uint64_t i = 0; i < (uint64_t) w*h; i++)
	    if (flag != 0)
		{
		    free(x.heights[i]);
		    if (flag>=6) // interpolation
			free(x.pos[i]);
		}
	free(x.heights);
	free(x.pos);
	
	
	if ( (flag>=6) && (radius>0) )//interpolation, remove extra pixels due to neighboring
	{
	    struct images x_cut;
	    int w_orig = 1 + (xmax_orig - xmin_orig) / resolution;
	    int h_orig = 1 + (ymax_orig - ymin_orig) / resolution;
	    x_cut.pixel_value = xmalloc((uint64_t) w_orig*h_orig*sizeof(float));
	    
	    Position center_pos_orig;
	    
	    for (int ii = 0; ii < w_orig; ii++)
		for (int jj = 0; jj < h_orig; jj++)
		{
		    center_pos_orig.x=xmin_orig+resolution/2.0 +(xmax_orig-xmin_orig)/( (float) w_orig)*ii;
		    center_pos_orig.y=-ymax_orig+resolution/2.0 +(ymax_orig-ymin_orig)/( (float) h_orig)*jj;
		    
		    int iic = (int) ((center_pos_orig.x-xmin-resolution/2.0)*w/(xmax-xmin)+0.5);
		    int jjc = (int) ((center_pos_orig.y+ymax-resolution/2.0)*h/(ymax-ymin)+0.5);
		    
		    uint64_t k = (uint64_t) w * jjc + iic;
		    
		    if (k < (uint64_t) w*h)
		    {
			uint64_t k_orig = (uint64_t) w_orig * jj + ii;
			x_cut.pixel_value[k_orig] = x.pixel_value[k];
		    }
		}
	    
	    // save output image
	    iio_save_image_float(out_dsm, x_cut.pixel_value, w_orig, h_orig);
	    set_geotif_header(out_dsm, utm, xmin_orig, ymax_orig, resolution);
	    
	    // cleanup and exit (2)
	    free(x.pixel_value);
	    free(x_cut.pixel_value);
	    
	}
	else
	{
	    // save output image
	    iio_save_image_float(out_dsm, x.pixel_value, w, h);
	    set_geotif_header(out_dsm, utm, xmin, ymax, resolution);
	    
	    // cleanup and exit (2)
	    free(x.pixel_value);
	}

	return 0;
}
