// rational polynomial coefficient stuff

#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>


// rational polynomial coefficients (and those of the inverse model)
struct rpc {
	double numx[20];
	double denx[20];
	double numy[20];
	double deny[20];
	double scale[3], offset[3];

	double inumx[20];
	double idenx[20];
	double inumy[20];
	double ideny[20];
	double iscale[3], ioffset[3];

	double dmval[4];
	double imval[4];
};


// read an XML file specifying an RPC model
void read_rpc_file_xml(struct rpc *p, char *filename);

void print_rpc(FILE *f, struct rpc *p, char *n);

// evaluate the direct rpc model
void eval_rpc(double *result,
		struct rpc *p, double x, double y, double z);

// evaluate the inverse rpc model
void eval_rpci(double *result,
		struct rpc *p, double x, double y, double z);

// evaluate an epipolar correspondence
static void eval_rpc_pair(double xprime[2],
		struct rpc *a, struct rpc *b,
		double x, double y, double z);


// convert (long,lat,h) to ECEF coord sys. (X,Y,Z)
void geotedic_to_ECEF(double lg, double lt, double h, 
			    double *X, double *Y, double *Z);

// given (X,Y,Z) ECEF coord, computes the alt above the WGS 84 ellipsoid
// (faster than ECEF_to_lgt_lat_alt, but only gives the altitude)
double get_altitude_from_ECEF(double X, double Y, double Z);

// given (X,Y,Z) ECEF coord, computes the lat/long/alt 
// above the WGS 84 ellipsoid
void ECEF_to_lgt_lat_alt(double X, double Y, double Z,
						 double *lgt, double *lat, double *h);


// 3x3 matrix object 
struct mat3x3 {
    double val[3][3];
};

// struct that defines the required 
// information to process a pair of tiles 
typedef struct {
    
    int ID;
    int sight_slave;
    bool process;
    
    // RPC
    struct rpc rpc_master;
    struct rpc rpc_slave;
    
    // homographies
    struct mat3x3 H_ref;
    struct mat3x3 H_sec;
    
    // pointing corrections
    struct mat3x3 pointing_correc;
    
    // take into account
    // pointing corrections
    // H_secB = H_sec * A^(-1)
    struct mat3x3 H_secB;
    
    // inverse slave homographies
    struct mat3x3 H_invSecB;
    
    // disparity maps
    int nx,ny,nch;
    float *dispx;
    float *dispy;
    
    // mask
    float *mask;
    
    double q0[3]; // a point inside ref image
    double q1[3]; // a point inside a given slave image
    
} type_pair;

// A simple list of pairs of tiles.
typedef struct
{
    type_pair *data;
    int tot_size;
    int real_size;
} list_of_pairs;

// define sight object
typedef struct 
{
    // sight ID
    int ID;
    
    // boolean telling whether this sight
    // is part of a consensus
    bool consensus;
    
    // the sight passes through points "s" ans "p"
    // (in ECEF coord.)
    double s[3]; 
    double p[3]; 
    // its unit direction vector is "v"
    double v[3]; 
    
    // distance from the optimal point to this sight
    double err;
    
    // indices 0 to 2 : optimal point in  ECEF coord
    // indices 3 to 5 : closest point in this sight
    // to optimal point in  ECEF coord.
    // Finally, err-...[i+3] - err-...[i] gives
    // the ith component of the smallest vector starting
    // from the optimal point and ending to a point in this sight
    double err_vec_ptopt_to_sight[6];
    
} type_sight;

// distance between a 3D point P and a line (defined by its
// normalized direction vector V and passing through a 3D point S) 
double dist_line_point3D(double *V,double *S,double *P);

// compute the vector VEC between a 3D point P0 and a line (
// passing through two 3D points P and S)
double vec_pt_to_line3D(double *P,double *S,double *P0,double *VEC);

// find closest 3D point from from a set of 3D lines
void find_point_opt(type_sight *sights_list, int N, bool *take,
		double *point_opt,double *outerr);

// compute the height of a point given its location inside two images
// geometric solution
double rpc_height_geo(list_of_pairs *list_pairs, 
		int local_nb_sights, int *final_nb_sights, 
		bool findConsensus, double thres, type_sight *sights_list);
